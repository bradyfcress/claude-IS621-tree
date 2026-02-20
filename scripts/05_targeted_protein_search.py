#!/usr/bin/env python3
"""Targeted per-species protein searches for EF-Tu, SecY, and MuA.

For each of the ~173 GTDB representative species, searches NCBI protein
database for the organism's ortholog and computes % amino acid identity
vs the E. coli reference sequence using local blastp.

This replaces the sparse genus-matching from 03v2, which relied on a
5000-hit BLAST cache that couldn't cover most GTDB genera.

Reads:
    data/taxonomy.csv           (species list with genus/domain info)
    data/metadata.csv           (existing metadata with IS621 from BLAST cache)
    data/ecoli_proteins.fasta   (E. coli reference sequences)

Writes:
    data/metadata.csv           (updated with targeted EF_Tu, SecY, MuA)
    data/targeted_cache/        (JSON cache of fetched sequences)
"""

import csv
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
CACHE_DIR = Path("data/targeted_cache")
CACHE_FILE = CACHE_DIR / "search_results.json"

Entrez.email = "is621tree@example.com"
REQUEST_DELAY = 0.35  # seconds between Entrez requests (3/sec limit)

PROTEINS = ["EF_Tu", "SecY", "MuA"]

# Minimum sequence length to accept (filters fragments)
MIN_SEQ_LEN = {"EF_Tu": 200, "SecY": 200, "MuA": 300}

# Search terms per protein per domain, tried in order (first hit wins).
# Each tuple: (search_term, use_refseq_filter)
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
def is_searchable_genus(genus: str) -> bool:
    """True if genus is a real organism name, not an accession-style MAG ID."""
    clean = re.sub(r"_[A-Z]+$", "", genus)  # Strip GTDB suffix
    if not clean or len(clean) < 4:
        return False
    # Latin genus convention: uppercase first, then 3+ lowercase, all alphabetic
    return bool(re.match(r"^[A-Z][a-z]{3,}$", clean))


def clean_genus(genus: str) -> str:
    """Remove GTDB suffix from genus name."""
    return re.sub(r"_[A-Z]+$", "", genus)


def load_ecoli_references() -> dict[str, str]:
    """Load E. coli reference protein sequences."""
    refs = {}
    for record in SeqIO.parse(str(ECOLI_FASTA), "fasta"):
        for protein in PROTEINS:
            if protein in record.id:
                refs[protein] = str(record.seq)
    return refs


def entrez_search_protein(genus: str, domain: str, protein: str) -> str | None:
    """Search NCBI protein for a genus + protein. Returns sequence or None."""
    strategies = SEARCH_STRATEGIES[protein].get(
        domain, SEARCH_STRATEGIES[protein]["Bacteria"]
    )
    min_len = MIN_SEQ_LEN.get(protein, 100)

    for search_term in strategies:
        # Try RefSeq first (better quality), then all records
        for db_filter in [" AND refseq[filter]", ""]:
            query = f'"{genus}"[Organism] AND ({search_term}){db_filter}'
            try:
                time.sleep(REQUEST_DELAY)
                handle = Entrez.esearch(db="protein", term=query, retmax=5)
                results = Entrez.read(handle)
                handle.close()

                if not results["IdList"]:
                    continue

                time.sleep(REQUEST_DELAY)
                handle = Entrez.efetch(
                    db="protein",
                    id=",".join(results["IdList"][:5]),
                    rettype="fasta",
                    retmode="text",
                )
                text = handle.read()
                handle.close()

                # Pick best: longest non-partial sequence above min length
                best_seq = None
                best_len = 0
                for record in SeqIO.parse(StringIO(text), "fasta"):
                    seq = str(record.seq).rstrip("*")  # Strip stop codon
                    desc = record.description.lower()
                    if len(seq) < min_len:
                        continue
                    is_partial = "partial" in desc or "fragment" in desc
                    # Prefer non-partial; among equals, prefer longer
                    if best_seq is None or (not is_partial and len(seq) > best_len):
                        best_seq = seq
                        best_len = len(seq)

                if best_seq:
                    return best_seq

            except Exception as e:
                print(f"      Entrez error ({genus}/{protein}): {e}", flush=True)
                time.sleep(2)

    return None


def batch_compute_identities(
    refs: dict[str, str], cache: dict[str, dict[str, str]]
) -> dict[str, dict[str, float]]:
    """Compute identities for all fetched proteins using local blastp.

    For each protein, creates a single BLAST database from the E. coli
    reference and queries all fetched sequences at once. Much faster than
    individual blastp calls.
    """
    genus_identities: dict[str, dict[str, float]] = {}

    for protein in PROTEINS:
        ref_seq = refs.get(protein)
        if not ref_seq:
            continue

        # Collect all genera with this protein
        sequences = {}
        for genus, proteins_found in cache.items():
            if protein in proteins_found:
                sequences[genus] = proteins_found[protein]

        if not sequences:
            print(f"  {protein}: no sequences to align", flush=True)
            continue

        # Write reference FASTA
        ref_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        )
        ref_tmp.write(f">ecoli_{protein}\n{ref_seq}\n")
        ref_tmp.close()

        # Create BLAST database from reference
        db_path = ref_tmp.name + "_db"
        subprocess.run(
            [
                "makeblastdb", "-in", ref_tmp.name,
                "-dbtype", "prot", "-out", db_path,
            ],
            capture_output=True,
            check=True,
        )

        # Write all query sequences to a single FASTA
        query_tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        )
        for genus, seq in sequences.items():
            # Clean sequence: remove stop codons, non-standard chars
            clean_seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
            query_tmp.write(f">{genus}\n{clean_seq}\n")
        query_tmp.close()

        # Run BLAST (all queries at once against E. coli reference)
        result = subprocess.run(
            [
                "blastp", "-query", query_tmp.name, "-db", db_path,
                "-outfmt", "6 qseqid pident", "-max_target_seqs", "1",
                "-evalue", "1e-3", "-num_threads", "4",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )

        n_hits = 0
        if result.returncode == 0:
            for line in result.stdout.strip().split("\n"):
                if not line.strip():
                    continue
                parts = line.split("\t")
                genus = parts[0]
                identity = round(float(parts[1]), 1)
                if genus not in genus_identities:
                    genus_identities[genus] = {}
                genus_identities[genus][protein] = identity
                n_hits += 1

        print(
            f"  {protein}: {n_hits}/{len(sequences)} genera aligned",
            flush=True,
        )

        # Cleanup temp files
        Path(ref_tmp.name).unlink(missing_ok=True)
        Path(query_tmp.name).unlink(missing_ok=True)
        for ext in [
            ".pdb", ".phr", ".pin", ".psq", ".pog",
            ".pot", ".ptf", ".pto", ".pjs",
        ]:
            Path(db_path + ext).unlink(missing_ok=True)

    return genus_identities


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60, flush=True)
    print("Targeted protein search: EF-Tu, SecY, MuA", flush=True)
    print("=" * 60, flush=True)

    # Verify local blastp
    if not shutil.which("blastp"):
        print("ERROR: blastp not found in PATH. Install BLAST+ first.", flush=True)
        sys.exit(1)
    print("Local blastp: available", flush=True)

    # ------------------------------------------------------------------
    # [1/4] Load data
    # ------------------------------------------------------------------
    print("\n[1/4] Loading data...", flush=True)
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
    print(f"  {len(species_info)} species in taxonomy", flush=True)

    # Deduplicate searchable genera
    genera_to_search: dict[str, dict] = {}
    for tip_id, info in species_info.items():
        genus_raw = info["genus"]
        if is_searchable_genus(genus_raw):
            gc = clean_genus(genus_raw)
            if gc not in genera_to_search:
                genera_to_search[gc] = {
                    "domain": info["domain"],
                    "tip_ids": [],
                }
            genera_to_search[gc]["tip_ids"].append(tip_id)

    n_species_covered = sum(len(g["tip_ids"]) for g in genera_to_search.values())
    print(
        f"  {len(genera_to_search)} searchable genera "
        f"(covering {n_species_covered}/{len(species_info)} species)",
        flush=True,
    )

    # Load cache
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache: dict[str, dict[str, str]] = {}
    if CACHE_FILE.exists():
        with open(CACHE_FILE) as f:
            cache = json.load(f)
        print(f"  Loaded {len(cache)} cached genus results", flush=True)

    # ------------------------------------------------------------------
    # [2/4] Search NCBI
    # ------------------------------------------------------------------
    n_to_search = sum(1 for g in genera_to_search if g not in cache)
    print(
        f"\n[2/4] Searching NCBI protein for {len(genera_to_search)} genera "
        f"({n_to_search} uncached)...",
        flush=True,
    )

    for i, (genus, info) in enumerate(sorted(genera_to_search.items()), 1):
        domain = info["domain"]

        if genus in cache:
            hits = [p for p in PROTEINS if p in cache[genus]]
            print(
                f"  [{i}/{len(genera_to_search)}] {genus}: cached "
                f"({len(hits)}/3 proteins)",
                flush=True,
            )
            continue

        print(
            f"  [{i}/{len(genera_to_search)}] {genus} ({domain})...",
            end="",
            flush=True,
        )
        genus_results: dict[str, str] = {}

        for protein in PROTEINS:
            seq = entrez_search_protein(genus, domain, protein)
            if seq:
                genus_results[protein] = seq
                print(f" {protein}({len(seq)}aa)", end="", flush=True)
            else:
                print(f" {protein}✗", end="", flush=True)

        cache[genus] = genus_results
        print(flush=True)

        # Save cache every 10 genera
        if i % 10 == 0:
            with open(CACHE_FILE, "w") as f:
                json.dump(cache, f)

    # Final cache save
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f)
    print(f"\n  Cache saved: {len(cache)} genera", flush=True)

    # ------------------------------------------------------------------
    # [3/4] Compute identities via local blastp
    # ------------------------------------------------------------------
    print("\n[3/4] Computing identities via local blastp...", flush=True)
    genus_identities = batch_compute_identities(refs, cache)

    # ------------------------------------------------------------------
    # [4/4] Update metadata.csv
    # ------------------------------------------------------------------
    print("\n[4/4] Updating metadata.csv...", flush=True)

    rows = []
    updated = {p: 0 for p in PROTEINS}

    with open(METADATA_CSV) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            tip_id = row["Tip_ID"]
            info = species_info.get(tip_id, {})
            gc = clean_genus(info.get("genus", ""))

            if gc in genus_identities:
                for protein in PROTEINS:
                    col = f"{protein}_Identity"
                    if protein in genus_identities[gc]:
                        new_val = genus_identities[gc][protein]
                        if row[col] == "NA" or row[col] == "":
                            row[col] = new_val
                            updated[protein] += 1
            rows.append(row)

    with open(METADATA_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    n_eftu = sum(1 for r in rows if r["EF_Tu_Identity"] != "NA")
    n_secy = sum(1 for r in rows if r["SecY_Identity"] != "NA")
    n_mua = sum(1 for r in rows if r["MuA_Identity"] != "NA")
    n_is621 = sum(1 for r in rows if r["IS621_Identity"] != "NA")

    print(
        f"\n  Updated (NA → value): "
        f"EF_Tu +{updated['EF_Tu']}, "
        f"SecY +{updated['SecY']}, "
        f"MuA +{updated['MuA']}",
        flush=True,
    )
    print(f"\n{'='*40}", flush=True)
    print(f"Final Coverage", flush=True)
    print(f"{'='*40}", flush=True)
    print(f"EF-Tu: {n_eftu}/{len(rows)} species", flush=True)
    print(f"SecY:  {n_secy}/{len(rows)} species", flush=True)
    print(f"MuA:   {n_mua}/{len(rows)} species", flush=True)
    print(f"IS621: {n_is621}/{len(rows)} species (unchanged)", flush=True)

    print(f"\nDone! Regenerate the tree:", flush=True)
    print(f"  Rscript scripts/04_plot_tree.R", flush=True)


if __name__ == "__main__":
    main()
