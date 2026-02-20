#!/usr/bin/env python3
"""Targeted per-species protein searches for EF-Tu, SecY, and MuA.

Two-pass strategy to maximize coverage across all 173 GTDB representatives:

  Pass 1 — Genus search: For named genera (Staphylococcus, etc.), search
           NCBI protein by genus name + protein name.

  Pass 2 — Assembly/TaxID search: For MAGs with accession-style names,
           look up the NCBI assembly to get the TaxID, then search for
           proteins by TaxID. Even broad-level TaxIDs (phylum/class)
           yield representative EF-Tu/SecY sequences with correct
           lineage-level identity.

Reads:
    data/taxonomy.csv           (species list with genus/domain info)
    data/metadata.csv           (existing metadata with IS621 from BLAST cache)
    data/ecoli_proteins.fasta   (E. coli reference sequences)
    data/gtdb_cache/*.tsv.gz    (GTDB taxonomy for genome ID lookup)

Writes:
    data/metadata.csv           (updated with targeted protein identities)
    data/targeted_cache/        (JSON caches for re-runs)
"""

import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TAXONOMY_CSV = Path("data/taxonomy.csv")
METADATA_CSV = Path("data/metadata.csv")
ECOLI_FASTA = Path("data/ecoli_proteins.fasta")
GTDB_DIR = Path("data/gtdb_cache")
CACHE_DIR = Path("data/targeted_cache")
GENUS_CACHE = CACHE_DIR / "search_results.json"
ASSEMBLY_INFO_CACHE = CACHE_DIR / "assembly_info.json"
TAXID_CACHE = CACHE_DIR / "taxid_results.json"

Entrez.email = "is621tree@example.com"
REQUEST_DELAY = 0.35  # seconds between Entrez calls (3/sec limit)

PROTEINS = ["EF_Tu", "SecY", "MuA"]
PASS2_PROTEINS = ["EF_Tu", "SecY"]  # Skip MuA in assembly search (phage gene)

MIN_SEQ_LEN = {"EF_Tu": 200, "SecY": 200, "MuA": 300}

SEARCH_STRATEGIES = {
    "EF_Tu": {
        "Bacteria": [
            '"elongation factor Tu"[Protein Name]',
            'tufA[Gene Name]',
            '"elongation factor Tu"[Title]',
        ],
        "Archaea": [
            '"elongation factor 1-alpha"[Protein Name]',
            '"elongation factor Tu"[Protein Name]',
            '"elongation factor"[Title]',
        ],
    },
    "SecY": {
        "Bacteria": [
            '"preprotein translocase subunit SecY"[Protein Name]',
            'secY[Gene Name]',
            '"SecY"[Title] AND translocase[Title]',
        ],
        "Archaea": [
            '"preprotein translocase subunit SecY"[Protein Name]',
            'secY[Gene Name]',
            '"SecY"[Title]',
        ],
    },
    "MuA": {
        "Bacteria": [
            'MuA[Gene Name]',
            '"phage Mu transposase"[Protein Name]',
            '"Mu-like"[Title] AND transposase[Title]',
        ],
        "Archaea": [
            'MuA[Gene Name]',
        ],
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def sanitize_tip_label(name: str) -> str:
    """Make a tree-safe tip label (must match 03v2 logic exactly)."""
    label = re.sub(r"[^A-Za-z0-9_]", "_", name)
    label = re.sub(r"_+", "_", label)
    return label.strip("_")[:60]


def is_searchable_genus(genus: str) -> bool:
    """True if genus is a real Latin name, not an accession-style MAG ID."""
    clean = re.sub(r"_[A-Z]+$", "", genus)
    if not clean or len(clean) < 4:
        return False
    return bool(re.match(r"^[A-Z][a-z]{3,}$", clean))


def clean_genus(genus: str) -> str:
    """Remove GTDB suffix (_A, _B, etc.) from genus name."""
    return re.sub(r"_[A-Z]+$", "", genus)


def load_ecoli_references() -> dict[str, str]:
    """Load E. coli reference protein sequences."""
    refs = {}
    for record in SeqIO.parse(str(ECOLI_FASTA), "fasta"):
        for protein in PROTEINS:
            if protein in record.id:
                refs[protein] = str(record.seq)
    return refs


def parse_gtdb_genome_ids() -> dict[str, str]:
    """Parse GTDB taxonomy → {sanitized_tip_label: assembly_accession}."""
    mapping = {}
    for gz_file in sorted(GTDB_DIR.glob("*.tsv.gz")):
        with gzip.open(gz_file, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) != 2:
                    continue
                genome_id = parts[0]
                ranks = parts[1].split(";")
                if len(ranks) != 7:
                    continue
                species = ranks[6].strip()
                if species.startswith("s__"):
                    species = species[3:]
                if not species:
                    continue
                tip = sanitize_tip_label(species)
                # Extract assembly accession (strip GB_/RS_ prefix)
                accn = genome_id
                for pfx in ("GB_", "RS_"):
                    if accn.startswith(pfx):
                        accn = accn[len(pfx):]
                mapping[tip] = accn
    return mapping


def _entrez_search_and_fetch(query: str, min_len: int) -> str | None:
    """Run Entrez esearch + efetch, return best protein sequence or None."""
    try:
        time.sleep(REQUEST_DELAY)
        handle = Entrez.esearch(db="protein", term=query, retmax=5)
        results = Entrez.read(handle)
        handle.close()

        if not results["IdList"]:
            return None

        time.sleep(REQUEST_DELAY)
        handle = Entrez.efetch(
            db="protein",
            id=",".join(results["IdList"][:5]),
            rettype="fasta",
            retmode="text",
        )
        text = handle.read()
        handle.close()

        best_seq = None
        best_len = 0
        for record in SeqIO.parse(StringIO(text), "fasta"):
            seq = str(record.seq).rstrip("*")
            desc = record.description.lower()
            if len(seq) < min_len:
                continue
            is_partial = "partial" in desc or "fragment" in desc
            if best_seq is None or (not is_partial and len(seq) > best_len):
                best_seq = seq
                best_len = len(seq)

        return best_seq
    except Exception as e:
        print(f" [err:{e}]", end="", flush=True)
        time.sleep(2)
        return None


def search_by_organism(organism: str, domain: str, protein: str) -> str | None:
    """Search NCBI protein by organism name + protein name."""
    strategies = SEARCH_STRATEGIES[protein].get(
        domain, SEARCH_STRATEGIES[protein]["Bacteria"]
    )
    min_len = MIN_SEQ_LEN.get(protein, 100)

    for search_term in strategies:
        for db_filter in [" AND refseq[filter]", ""]:
            query = f'"{organism}"[Organism] AND ({search_term}){db_filter}'
            seq = _entrez_search_and_fetch(query, min_len)
            if seq:
                return seq
    return None


def search_by_taxid(taxid: str, domain: str, protein: str) -> str | None:
    """Search NCBI protein by TaxID + protein name."""
    strategies = SEARCH_STRATEGIES[protein].get(
        domain, SEARCH_STRATEGIES[protein]["Bacteria"]
    )
    min_len = MIN_SEQ_LEN.get(protein, 100)

    for search_term in strategies:
        for db_filter in [" AND refseq[filter]", ""]:
            query = f"txid{taxid}[Organism] AND ({search_term}){db_filter}"
            seq = _entrez_search_and_fetch(query, min_len)
            if seq:
                return seq
    return None


def get_assembly_info(accession: str) -> dict | None:
    """Look up NCBI assembly to get TaxID and organism name."""
    try:
        time.sleep(REQUEST_DELAY)
        handle = Entrez.esearch(db="assembly", term=accession)
        result = Entrez.read(handle)
        handle.close()

        if not result["IdList"]:
            return None

        time.sleep(REQUEST_DELAY)
        handle = Entrez.esummary(db="assembly", id=result["IdList"][0])
        summary = Entrez.read(handle)
        handle.close()

        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        return {
            "taxid": str(doc.get("Taxid", "")),
            "organism": str(doc.get("Organism", "")),
        }
    except Exception as e:
        print(f" [asm-err:{e}]", end="", flush=True)
        time.sleep(2)
        return None


def batch_compute_identities(
    refs: dict[str, str], all_sequences: dict[str, dict[str, str]]
) -> dict[str, dict[str, float]]:
    """Compute identities for all fetched proteins using local blastp."""
    identities: dict[str, dict[str, float]] = {}

    for protein in PROTEINS:
        ref_seq = refs.get(protein)
        if not ref_seq:
            continue

        sequences = {}
        for key, proteins_found in all_sequences.items():
            if protein in proteins_found:
                sequences[key] = proteins_found[protein]

        if not sequences:
            print(f"  {protein}: no sequences to align", flush=True)
            continue

        # Write reference FASTA and create BLAST db
        ref_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        )
        ref_tmp.write(f">ecoli_{protein}\n{ref_seq}\n")
        ref_tmp.close()

        db_path = ref_tmp.name + "_db"
        subprocess.run(
            ["makeblastdb", "-in", ref_tmp.name, "-dbtype", "prot", "-out", db_path],
            capture_output=True, check=True,
        )

        # Write all query sequences
        query_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        )
        for key, seq in sequences.items():
            clean_seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
            query_tmp.write(f">{key}\n{clean_seq}\n")
        query_tmp.close()

        # Run BLAST
        result = subprocess.run(
            [
                "blastp", "-query", query_tmp.name, "-db", db_path,
                "-outfmt", "6 qseqid pident", "-max_target_seqs", "1",
                "-evalue", "1e-3", "-num_threads", "4",
            ],
            capture_output=True, text=True, timeout=120,
        )

        n_hits = 0
        if result.returncode == 0:
            for line in result.stdout.strip().split("\n"):
                if not line.strip():
                    continue
                parts = line.split("\t")
                key = parts[0]
                identity = round(float(parts[1]), 1)
                if key not in identities:
                    identities[key] = {}
                identities[key][protein] = identity
                n_hits += 1

        print(f"  {protein}: {n_hits}/{len(sequences)} aligned", flush=True)

        # Cleanup
        Path(ref_tmp.name).unlink(missing_ok=True)
        Path(query_tmp.name).unlink(missing_ok=True)
        for ext in [".pdb", ".phr", ".pin", ".psq", ".pog",
                     ".pot", ".ptf", ".pto", ".pjs"]:
            Path(db_path + ext).unlink(missing_ok=True)

    return identities


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60, flush=True)
    print("Targeted protein search: EF-Tu, SecY, MuA", flush=True)
    print("=" * 60, flush=True)

    if not shutil.which("blastp"):
        print("ERROR: blastp not found.", flush=True)
        sys.exit(1)

    # ==================================================================
    # [1/6] Load data
    # ==================================================================
    print("\n[1/6] Loading data...", flush=True)
    refs = load_ecoli_references()
    for p, seq in refs.items():
        print(f"  {p}: {len(seq)} aa reference", flush=True)

    species_info: dict[str, dict] = {}
    with open(TAXONOMY_CSV) as f:
        for row in csv.DictReader(f):
            species_info[row["Tip_ID"]] = {
                "genus": row["Genus"],
                "domain": row["Domain"],
            }
    print(f"  {len(species_info)} species", flush=True)

    # Parse GTDB for assembly accessions
    print("  Parsing GTDB for genome IDs...", flush=True)
    tip_to_assembly = parse_gtdb_genome_ids()
    n_matched = sum(1 for t in species_info if t in tip_to_assembly)
    print(f"  {n_matched}/{len(species_info)} matched to GTDB assemblies", flush=True)

    # Identify searchable genera (pass 1)
    genera_to_search: dict[str, dict] = {}
    genus_tip_map: dict[str, list[str]] = {}
    for tip_id, info in species_info.items():
        if is_searchable_genus(info["genus"]):
            gc = clean_genus(info["genus"])
            if gc not in genera_to_search:
                genera_to_search[gc] = {"domain": info["domain"]}
                genus_tip_map[gc] = []
            genus_tip_map[gc].append(tip_id)
    print(f"  {len(genera_to_search)} searchable genera (pass 1)", flush=True)

    # Load caches
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    genus_cache: dict = json.loads(GENUS_CACHE.read_text()) if GENUS_CACHE.exists() else {}
    asm_info_cache: dict = json.loads(ASSEMBLY_INFO_CACHE.read_text()) if ASSEMBLY_INFO_CACHE.exists() else {}
    taxid_cache: dict = json.loads(TAXID_CACHE.read_text()) if TAXID_CACHE.exists() else {}

    # ==================================================================
    # [2/6] Pass 1: Genus-based search
    # ==================================================================
    n_uncached = sum(1 for g in genera_to_search if g not in genus_cache)
    print(
        f"\n[2/6] Pass 1: Genus search "
        f"({len(genera_to_search)} genera, {n_uncached} uncached)...",
        flush=True,
    )

    for i, (genus, info) in enumerate(sorted(genera_to_search.items()), 1):
        if genus in genus_cache:
            continue
        domain = info["domain"]
        print(f"  [{i}/{len(genera_to_search)}] {genus} ({domain})...",
              end="", flush=True)
        genus_results: dict[str, str] = {}
        for protein in PROTEINS:
            seq = search_by_organism(genus, domain, protein)
            if seq:
                genus_results[protein] = seq
                print(f" {protein}({len(seq)}aa)", end="", flush=True)
            else:
                print(f" {protein}✗", end="", flush=True)
        genus_cache[genus] = genus_results
        print(flush=True)
        if i % 10 == 0:
            GENUS_CACHE.write_text(json.dumps(genus_cache))

    GENUS_CACHE.write_text(json.dumps(genus_cache))

    # Determine which tip_ids are already covered by genus search
    genus_covered_tips = set()
    for gc, tips in genus_tip_map.items():
        if gc in genus_cache and ("EF_Tu" in genus_cache[gc] or "SecY" in genus_cache[gc]):
            genus_covered_tips.update(tips)

    # ==================================================================
    # [3/6] Pass 2: Assembly TaxID search
    # ==================================================================
    need_pass2 = [
        t for t in species_info
        if t not in genus_covered_tips and t in tip_to_assembly
    ]
    print(
        f"\n[3/6] Pass 2: Assembly/TaxID search ({len(need_pass2)} species)...",
        flush=True,
    )

    # Step A: Get assembly info (TaxID, organism) for species we need
    n_asm_uncached = sum(1 for t in need_pass2 if tip_to_assembly[t] not in asm_info_cache)
    print(f"  Looking up {n_asm_uncached} uncached assemblies...", flush=True)

    for i, tip_id in enumerate(sorted(need_pass2), 1):
        accn = tip_to_assembly[tip_id]
        if accn in asm_info_cache:
            continue
        info = get_assembly_info(accn)
        asm_info_cache[accn] = info or {"taxid": "", "organism": ""}
        if info:
            print(f"  [{i}] {accn} -> txid{info['taxid']} {info['organism'][:40]}",
                  flush=True)
        else:
            print(f"  [{i}] {accn} -> not found", flush=True)

        if i % 20 == 0:
            ASSEMBLY_INFO_CACHE.write_text(json.dumps(asm_info_cache))

    ASSEMBLY_INFO_CACHE.write_text(json.dumps(asm_info_cache))

    # Step B: Search proteins by TaxID (deduplicated — same TaxID = same search)
    taxid_to_tips: dict[str, list[str]] = {}
    for tip_id in need_pass2:
        accn = tip_to_assembly[tip_id]
        info = asm_info_cache.get(accn, {})
        txid = info.get("taxid", "")
        if txid:
            taxid_to_tips.setdefault(txid, []).append(tip_id)

    n_taxids_uncached = sum(1 for t in taxid_to_tips if t not in taxid_cache)
    print(
        f"  {len(taxid_to_tips)} unique TaxIDs "
        f"({n_taxids_uncached} uncached protein searches)...",
        flush=True,
    )

    for i, (txid, tips) in enumerate(sorted(taxid_to_tips.items()), 1):
        if txid in taxid_cache:
            continue

        # Get domain from first tip
        domain = species_info[tips[0]]["domain"]
        # Get organism name for display
        first_accn = tip_to_assembly[tips[0]]
        organism = asm_info_cache.get(first_accn, {}).get("organism", "?")[:35]

        print(
            f"  [{i}/{len(taxid_to_tips)}] txid{txid} ({organism}, "
            f"{len(tips)} species)...",
            end="", flush=True,
        )

        taxid_results: dict[str, str] = {}
        for protein in PASS2_PROTEINS:
            seq = search_by_taxid(txid, domain, protein)
            if seq:
                taxid_results[protein] = seq
                print(f" {protein}({len(seq)}aa)", end="", flush=True)
            else:
                print(f" {protein}✗", end="", flush=True)

        taxid_cache[txid] = taxid_results
        print(flush=True)

        if i % 10 == 0:
            TAXID_CACHE.write_text(json.dumps(taxid_cache))

    TAXID_CACHE.write_text(json.dumps(taxid_cache))

    # ==================================================================
    # [4/6] Merge all sequences into a unified collection
    # ==================================================================
    print("\n[4/6] Merging sequences...", flush=True)
    all_sequences: dict[str, dict[str, str]] = {}

    # Add genus results (keyed as "g_<genus>")
    for genus, data in genus_cache.items():
        if data:
            all_sequences[f"g_{genus}"] = data

    # Add taxid results (keyed as "t_<taxid>")
    for txid, data in taxid_cache.items():
        if data:
            all_sequences[f"t_{txid}"] = data

    print(
        f"  {len(all_sequences)} entries "
        f"({sum(1 for k in all_sequences if k.startswith('g_'))} genus + "
        f"{sum(1 for k in all_sequences if k.startswith('t_'))} taxid)",
        flush=True,
    )

    # ==================================================================
    # [5/6] Compute identities via local blastp
    # ==================================================================
    print("\n[5/6] Computing identities via local blastp...", flush=True)
    all_identities = batch_compute_identities(refs, all_sequences)

    # Build per-tip_id identity lookup
    tip_identities: dict[str, dict[str, float]] = {}

    # From genus results
    for gc, tips in genus_tip_map.items():
        key = f"g_{gc}"
        if key in all_identities:
            for tip_id in tips:
                tip_identities[tip_id] = dict(all_identities[key])

    # From taxid results (fill in what genus didn't cover)
    for txid, tips in taxid_to_tips.items():
        key = f"t_{txid}"
        if key in all_identities:
            for tip_id in tips:
                if tip_id not in tip_identities:
                    tip_identities[tip_id] = {}
                for protein, identity in all_identities[key].items():
                    if protein not in tip_identities[tip_id]:
                        tip_identities[tip_id][protein] = identity

    print(f"  {len(tip_identities)} species with at least one identity", flush=True)

    # ==================================================================
    # [6/6] Update metadata.csv
    # ==================================================================
    print("\n[6/6] Updating metadata.csv...", flush=True)

    rows = []
    updated = {p: 0 for p in PROTEINS}

    with open(METADATA_CSV) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            tip_id = row["Tip_ID"]

            if tip_id in tip_identities:
                for protein in PROTEINS:
                    col = f"{protein}_Identity"
                    if protein in tip_identities[tip_id]:
                        new_val = tip_identities[tip_id][protein]
                        if row[col] == "NA" or row[col] == "":
                            row[col] = new_val
                            updated[protein] += 1
            rows.append(row)

    with open(METADATA_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    n = {}
    for protein in PROTEINS:
        n[protein] = sum(1 for r in rows if r[f"{protein}_Identity"] != "NA")

    print(
        f"\n  Updated (NA → value): "
        f"EF_Tu +{updated['EF_Tu']}, SecY +{updated['SecY']}, "
        f"MuA +{updated['MuA']}",
        flush=True,
    )
    print(f"\n{'='*40}", flush=True)
    print("Final Coverage", flush=True)
    print(f"{'='*40}", flush=True)
    print(f"EF-Tu: {n['EF_Tu']}/{len(rows)} species", flush=True)
    print(f"SecY:  {n['SecY']}/{len(rows)} species", flush=True)
    print(f"MuA:   {n['MuA']}/{len(rows)} species", flush=True)
    print(f"IS621: {n.get('IS621', sum(1 for r in rows if r['IS621_Identity'] != 'NA'))}/{len(rows)} species (unchanged)", flush=True)

    print(f"\nDone! Regenerate the tree:", flush=True)
    print(f"  Rscript scripts/04_plot_tree.R", flush=True)


if __name__ == "__main__":
    main()
