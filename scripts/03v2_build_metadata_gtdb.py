#!/usr/bin/env python3
"""Build tree metadata using GTDB taxonomy as the backbone (V2/V3).

The V1 approach (03_build_metadata.py) started from BLAST hits, producing a
Pseudomonadota-biased tree with no Archaea. This V2 script flips the logic:

  1. Download GTDB taxonomy (bac120 + ar53)
  2. Select ~180 balanced representatives across all prokaryotic phyla
  3. Map existing BLAST results onto those representatives (genus-level)

V3 addition: --complete-only flag filters to Complete Genome / Chromosome
assemblies, ensuring high-quality genomes with reliable gene content.

Reads:  data/blast_cache/{protein}_hits.tsv (4 files, from V1)
        data/gtdb_cache/*_metadata.tsv.gz   (V3: assembly quality info)
Writes: data/taxonomy.csv   (tip taxonomy for R tree building)
        data/metadata.csv   (identity values for R heatmap rings)
"""

import argparse
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
GTDB_URLS = {
    "bacteria": "https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz",
    "archaea": "https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz",
}
GTDB_DIR = Path("data/gtdb_cache")
CACHE_DIR = Path("data/blast_cache")
OUTPUT_METADATA = Path("data/metadata.csv")
OUTPUT_TAXONOMY = Path("data/taxonomy.csv")

PROTEINS = ["IS621", "EF_Tu", "SecY", "MuA", "IS911"]

GTDB_METADATA_URLS = {
    "bacteria": "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz",
    "archaea": "https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz",
}
TARGET_TOTAL = 180

# Tier thresholds: (min_families, max_reps)
# Phyla are classified by how many distinct families they have in GTDB.
TIER_RULES = [
    (200, 15),  # Large phyla: up to 15 reps
    (50, 8),    # Medium phyla: up to 8 reps
    (10, 3),    # Small phyla: up to 3 reps
    (0, 2),     # Tiny phyla: up to 2 reps
]

RANK_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
RANK_NAMES = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


# ---------------------------------------------------------------------------
# GTDB download and parsing
# ---------------------------------------------------------------------------
def download_gtdb_taxonomy() -> dict[str, Path]:
    """Download GTDB taxonomy files if not already cached."""
    GTDB_DIR.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, url in GTDB_URLS.items():
        filename = url.rsplit("/", 1)[-1]
        dest = GTDB_DIR / filename
        if dest.exists() and dest.stat().st_size > 100_000:
            print(f"  Cache hit: {dest} ({dest.stat().st_size:,} bytes)", flush=True)
            paths[key] = dest
            continue
        print(f"  Downloading {key}: {url}", flush=True)
        try:
            urlretrieve(url, dest)
            print(f"  Saved: {dest} ({dest.stat().st_size:,} bytes)", flush=True)
            paths[key] = dest
        except (URLError, OSError) as e:
            print(f"  ERROR downloading {key}: {e}", flush=True)
            if dest.exists():
                print(f"  Using stale cache: {dest}", flush=True)
                paths[key] = dest
            else:
                print(f"  FATAL: No cached file and download failed. "
                      f"Manually download from {url}", flush=True)
                sys.exit(1)
    return paths


def parse_gtdb_taxonomy(paths: dict[str, Path]) -> list[dict]:
    """Parse GTDB taxonomy TSV.gz files into a list of species records.

    Each record: {genome_id, Domain, Phylum, Class, Order, Family, Genus, Species}
    """
    records = []
    for key, path in paths.items():
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                genome_id = parts[0]
                ranks = parts[1].split(";")
                if len(ranks) != 7:
                    continue

                rec = {"genome_id": genome_id}
                skip = False
                for prefix, name, value in zip(RANK_PREFIXES, RANK_NAMES, ranks):
                    cleaned = value.strip()
                    if cleaned.startswith(prefix):
                        cleaned = cleaned[len(prefix):]
                    # Skip entries with empty ranks (unclassified placeholders)
                    if not cleaned:
                        skip = True
                        break
                    rec[name] = cleaned
                if skip:
                    continue
                records.append(rec)

        print(f"  Parsed {key}: {len(records):,} total species so far", flush=True)
    return records


# ---------------------------------------------------------------------------
# GTDB metadata loading (assembly quality)
# ---------------------------------------------------------------------------
def load_gtdb_metadata() -> dict[str, dict]:
    """Load GTDB metadata files to get assembly quality info.

    Returns {accession: {ncbi_assembly_level, checkm2_completeness, checkm2_contamination}}.
    """
    metadata = {}
    for key, url in GTDB_METADATA_URLS.items():
        filename = url.rsplit("/", 1)[-1]
        path = GTDB_DIR / filename
        if not path.exists():
            print(f"  WARNING: {path} not found, download GTDB metadata first", flush=True)
            continue
        with gzip.open(path, "rt") as f:
            header = f.readline().strip().split("\t")
            accn_idx = header.index("accession")
            level_idx = header.index("ncbi_assembly_level")
            comp_idx = header.index("checkm2_completeness")
            cont_idx = header.index("checkm2_contamination")
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) <= max(accn_idx, level_idx, comp_idx, cont_idx):
                    continue
                accn = parts[accn_idx]
                try:
                    comp = float(parts[comp_idx]) if parts[comp_idx] else 0.0
                    cont = float(parts[cont_idx]) if parts[cont_idx] else 100.0
                except ValueError:
                    comp, cont = 0.0, 100.0
                metadata[accn] = {
                    "ncbi_assembly_level": parts[level_idx],
                    "completeness": comp,
                    "contamination": cont,
                }
        print(f"  Loaded {key} metadata: {len(metadata):,} genomes", flush=True)
    return metadata


def filter_complete_genomes(
    gtdb_records: list[dict],
    metadata: dict[str, dict],
) -> list[dict]:
    """Filter GTDB records to complete/high-quality genomes.

    Primary: assembly_level in ('Complete Genome', 'Chromosome')
    Fallback: If a phylum has no complete genomes, allow high-quality scaffolds
    (completeness >= 95% AND contamination <= 5%).
    """
    # Classify each record
    complete = []
    hq_scaffold = []
    rest = []

    for rec in gtdb_records:
        meta = metadata.get(rec["genome_id"])
        if not meta:
            rest.append(rec)
            continue
        level = meta["ncbi_assembly_level"]
        if level in ("Complete Genome", "Chromosome"):
            rec["_assembly_level"] = level
            complete.append(rec)
        elif meta["completeness"] >= 95 and meta["contamination"] <= 5:
            rec["_assembly_level"] = f"{level} (HQ)"
            hq_scaffold.append(rec)
        else:
            rest.append(rec)

    # Find phyla represented only in HQ scaffolds (not in complete)
    complete_phyla = {r["Phylum"] for r in complete}
    fallback = [r for r in hq_scaffold if r["Phylum"] not in complete_phyla]

    filtered = complete + fallback

    # Stats
    filtered_phyla = {r["Phylum"] for r in filtered}
    all_phyla = {r["Phylum"] for r in gtdb_records}
    lost_phyla = all_phyla - filtered_phyla

    print(f"  Complete/Chromosome genomes: {len(complete):,}", flush=True)
    print(f"  HQ scaffold fallback (phyla w/o complete): {len(fallback):,}", flush=True)
    print(f"  Filtered total: {len(filtered):,} (from {len(gtdb_records):,})", flush=True)
    print(f"  Phyla retained: {len(filtered_phyla)} (lost {len(lost_phyla)})", flush=True)

    return filtered


# ---------------------------------------------------------------------------
# BLAST data loading
# ---------------------------------------------------------------------------
def load_blast_hits() -> dict[str, dict[str, float]]:
    """Load BLAST hit TSVs. Returns {protein: {organism: best_identity_pct}}."""
    hits = {}
    for protein in PROTEINS:
        tsv = CACHE_DIR / f"{protein}_hits.tsv"
        if not tsv.exists():
            print(f"  WARNING: {tsv} not found, skipping {protein}", flush=True)
            hits[protein] = {}
            continue

        protein_hits = {}
        skip_orgs = {"Bacteria", "Enterobacteriaceae", "unclassified sequences",
                     "Multispecies"}
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                org = row["organism"].strip()
                if not org or org in skip_orgs:
                    continue
                identity = float(row["identity_pct"])
                if org not in protein_hits or identity > protein_hits[org]:
                    protein_hits[org] = identity
        hits[protein] = protein_hits
        print(f"  {protein}: {len(protein_hits):,} unique organisms", flush=True)
    return hits


# ---------------------------------------------------------------------------
# Genus normalization and indexing
# ---------------------------------------------------------------------------
def normalize_genus(name: str) -> str:
    """Normalize a genus name for cross-database matching.

    Handles:
      - GTDB alphabetic suffixes: "Prevotella_A" -> "Prevotella"
      - Candidatus prefix: kept as-is
      - Case: lowered for matching
    """
    # Strip GTDB-style alphabetic suffixes (e.g., _A, _B, _AC)
    name = re.sub(r"_[A-Z]+$", "", name)
    return name.lower().strip()


def extract_genus_from_organism(organism: str) -> str:
    """Extract genus from an NCBI organism name."""
    parts = organism.strip().split()
    if not parts:
        return ""
    if parts[0].lower() == "candidatus" and len(parts) > 1:
        return f"Candidatus {parts[1]}"
    return parts[0]


def build_blast_genus_index(
    blast_hits: dict[str, dict[str, float]]
) -> dict[str, dict[str, float]]:
    """Build {normalized_genus: {protein: best_identity}} from BLAST data."""
    genus_index: dict[str, dict[str, float]] = defaultdict(dict)

    for protein, org_hits in blast_hits.items():
        for org, identity in org_hits.items():
            genus = normalize_genus(extract_genus_from_organism(org))
            if not genus:
                continue
            if protein not in genus_index[genus] or identity > genus_index[genus][protein]:
                genus_index[genus][protein] = identity

    return dict(genus_index)


# ---------------------------------------------------------------------------
# Representative selection
# ---------------------------------------------------------------------------
def classify_tier(n_families: int) -> int:
    """Return max representatives for a phylum based on its family count."""
    for min_fam, max_reps in TIER_RULES:
        if n_families >= min_fam:
            return max_reps
    return 1


def select_representatives(
    gtdb_records: list[dict],
    blast_genus_index: dict[str, dict[str, float]],
) -> list[dict]:
    """Select ~TARGET_TOTAL balanced representatives from GTDB taxonomy.

    Algorithm:
      1. Group species by phylum, count families per phylum
      2. Allocate reps per phylum by tier (family count)
      3. Within each phylum, select one rep per family, prioritizing BLAST hits
      4. Force E. coli inclusion
      5. Trim to TARGET_TOTAL if over budget
    """
    # --- Group by phylum ---
    phylum_families: dict[str, dict[str, list[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for rec in gtdb_records:
        phylum_families[rec["Phylum"]][rec["Family"]].append(rec)

    phylum_family_counts = {
        p: len(families) for p, families in phylum_families.items()
    }

    # Pre-compute normalized genus for all GTDB genera (avoid repeated regex)
    gtdb_genus_to_norm: dict[str, str] = {}
    for rec in gtdb_records:
        g = rec["Genus"]
        if g not in gtdb_genus_to_norm:
            gtdb_genus_to_norm[g] = normalize_genus(g)

    blast_genera = set(blast_genus_index.keys())

    # Compute per-phylum BLAST relevance
    phylum_blast_score = {}
    for phylum, families in phylum_families.items():
        score = 0
        for family, species_list in families.items():
            for sp in species_list:
                if gtdb_genus_to_norm.get(sp["Genus"], "") in blast_genera:
                    score += 1
                    break
        phylum_blast_score[phylum] = score

    # --- Filter phyla ---
    # GTDB has 200+ phyla; many are tiny candidate lineages with accession names
    # (e.g., "JAFGOL01"). Keep only phyla with ≥5 families, or ≥2 for Archaea,
    # or any BLAST relevance.
    MIN_FAMILIES = 5
    MIN_FAMILIES_ARCHAEA = 2
    included_phyla = set()
    for phylum, n_families in phylum_family_counts.items():
        # Check if archaeal
        first_family = next(iter(phylum_families[phylum]))
        is_archaeal = phylum_families[phylum][first_family][0]["Domain"] == "Archaea"
        min_fam = MIN_FAMILIES_ARCHAEA if is_archaeal else MIN_FAMILIES
        has_blast = phylum_blast_score.get(phylum, 0) > 0

        if n_families >= min_fam or has_blast:
            included_phyla.add(phylum)

    print(f"\n  GTDB phyla (total): {len(phylum_family_counts)}", flush=True)
    print(f"  Included phyla: {len(included_phyla)}", flush=True)
    print(f"  Excluded (tiny candidate phyla): "
          f"{len(phylum_family_counts) - len(included_phyla)}", flush=True)

    # --- Allocate per included phylum ---
    allocations = {}
    for phylum in included_phyla:
        n_families = phylum_family_counts[phylum]
        max_reps = classify_tier(n_families)
        allocations[phylum] = min(max_reps, n_families)

    total_allocated = sum(allocations.values())
    print(f"  Initial allocation: {total_allocated} across {len(allocations)} phyla",
          flush=True)

    # Scale down if over budget
    if total_allocated > TARGET_TOTAL:
        scale = TARGET_TOTAL / total_allocated
        for phylum in allocations:
            allocations[phylum] = max(1, round(allocations[phylum] * scale))
        total_allocated = sum(allocations.values())

        # Fine-tune: trim largest phyla first
        while total_allocated > TARGET_TOTAL:
            trimmable = {p: n for p, n in allocations.items() if n > 1}
            if not trimmable:
                break
            largest = max(trimmable, key=trimmable.get)
            allocations[largest] -= 1
            total_allocated -= 1

    print(f"  Final allocation: {total_allocated} across {len(allocations)} phyla",
          flush=True)

    # --- Select within each phylum ---
    all_reps = []
    for phylum, n_reps in sorted(allocations.items(), key=lambda x: -x[1]):
        families = phylum_families[phylum]

        # Score each family by BLAST coverage (check one species per genus)
        family_scores = []
        for family, species_list in families.items():
            best_score = 0.0
            best_species = species_list[0]
            seen_genera = set()
            for sp in species_list:
                genus_norm = gtdb_genus_to_norm.get(sp["Genus"], "")
                if genus_norm in seen_genera:
                    continue
                seen_genera.add(genus_norm)
                genus_data = blast_genus_index.get(genus_norm, {})
                # Prioritize IS621, then EF_Tu, then any hit
                is621 = genus_data.get("IS621", 0)
                eftu = genus_data.get("EF_Tu", 0)
                score = (1000 if is621 > 0 else 0) + is621 + (100 if eftu > 0 else 0) + eftu
                if score > best_score:
                    best_score = score
                    best_species = sp
            family_scores.append((family, best_species, best_score))

        # Sort: families with BLAST hits first, then by name for reproducibility
        family_scores.sort(key=lambda x: (-x[2], x[0]))

        # Take top n_reps
        selected = family_scores[:n_reps]
        for family, species, score in selected:
            all_reps.append(species)

    # --- Force E. coli ---
    ecoli_present = any(
        r["Species"] == "Escherichia coli" or
        (r["Genus"] == "Escherichia" and "coli" in r.get("Species", ""))
        for r in all_reps
    )
    if not ecoli_present:
        # Find E. coli in GTDB
        ecoli_rec = None
        for rec in gtdb_records:
            if rec["Genus"] == "Escherichia" and "coli" in rec.get("Species", ""):
                ecoli_rec = rec
                break
        if ecoli_rec:
            # Replace the lowest-scoring Pseudomonadota rep
            for i in range(len(all_reps) - 1, -1, -1):
                if all_reps[i]["Phylum"] == ecoli_rec["Phylum"]:
                    all_reps[i] = ecoli_rec
                    break
            else:
                all_reps.append(ecoli_rec)
            print("  Forced E. coli into representative set", flush=True)

    print(f"\n  Selected {len(all_reps)} total representatives", flush=True)

    # Summarize
    phylum_counts = defaultdict(int)
    domain_counts = defaultdict(int)
    for r in all_reps:
        phylum_counts[r["Phylum"]] += 1
        domain_counts[r["Domain"]] += 1

    for domain in sorted(domain_counts):
        print(f"  {domain}: {domain_counts[domain]} species", flush=True)
    print(flush=True)
    for phylum in sorted(phylum_counts, key=lambda p: -phylum_counts[p]):
        print(f"    {phylum}: {phylum_counts[phylum]}", flush=True)

    return all_reps


# ---------------------------------------------------------------------------
# BLAST mapping
# ---------------------------------------------------------------------------
def map_blast_to_reps(
    representatives: list[dict],
    blast_genus_index: dict[str, dict[str, float]],
) -> list[dict]:
    """Map BLAST identity values onto GTDB representatives via genus matching.

    Returns list of metadata dicts ready for CSV output.
    """
    metadata = []
    n_matched = defaultdict(int)

    for rep in representatives:
        genus_norm = normalize_genus(rep["Genus"])
        genus_data = blast_genus_index.get(genus_norm, {})

        is_ecoli = (
            rep["Genus"] == "Escherichia"
            and "coli" in rep.get("Species", "")
        )

        tip_id = sanitize_tip_label(rep["Species"])

        row = {
            "Tip_ID": tip_id,
            "Phylum": rep["Phylum"],
            "Ecoli_Reference": "TRUE" if is_ecoli else "FALSE",
        }

        for protein in PROTEINS:
            col = f"{protein}_Identity"
            if is_ecoli:
                row[col] = 100.0
            elif protein in genus_data:
                row[col] = round(genus_data[protein], 1)
                n_matched[protein] += 1
            else:
                row[col] = "NA"

        metadata.append(row)

    print("\n  BLAST mapping results:", flush=True)
    for protein in PROTEINS:
        print(f"    {protein}: {n_matched[protein]}/{len(representatives)} species matched",
              flush=True)

    return metadata


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def sanitize_tip_label(name: str) -> str:
    """Make a tree-safe tip label from a species name."""
    label = re.sub(r"[^A-Za-z0-9_]", "_", name)
    label = re.sub(r"_+", "_", label)
    return label.strip("_")[:60]


def write_outputs(
    representatives: list[dict],
    metadata: list[dict],
) -> None:
    """Write taxonomy.csv and metadata.csv for the R script."""
    # Build taxonomy rows, ensuring unique Tip_IDs
    seen_tips = set()
    taxonomy_rows = []
    for rep, meta in zip(representatives, metadata):
        tip_id = meta["Tip_ID"]
        if tip_id in seen_tips:
            suffix = 2
            while f"{tip_id}_{suffix}" in seen_tips:
                suffix += 1
            tip_id = f"{tip_id}_{suffix}"
            meta["Tip_ID"] = tip_id
        seen_tips.add(tip_id)

        taxonomy_rows.append({
            "Tip_ID": tip_id,
            "Domain": rep["Domain"],
            "Phylum": rep["Phylum"],
            "Class": rep["Class"],
            "Order": rep["Order"],
            "Family": rep["Family"],
            "Genus": rep["Genus"],
            "Species": tip_id,
        })

    # Write taxonomy CSV
    OUTPUT_TAXONOMY.parent.mkdir(parents=True, exist_ok=True)
    tax_fields = ["Tip_ID", "Domain", "Phylum", "Class", "Order", "Family",
                  "Genus", "Species"]
    with open(OUTPUT_TAXONOMY, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tax_fields)
        writer.writeheader()
        writer.writerows(taxonomy_rows)
    print(f"\n  Wrote {len(taxonomy_rows)} rows to {OUTPUT_TAXONOMY}", flush=True)

    # Write metadata CSV
    meta_fields = ["Tip_ID", "Phylum", "Ecoli_Reference",
                   "EF_Tu_Identity", "SecY_Identity", "MuA_Identity",
                   "IS621_Identity", "IS911_Identity"]
    with open(OUTPUT_METADATA, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=meta_fields)
        writer.writeheader()
        writer.writerows(metadata)
    print(f"  Wrote {len(metadata)} rows to {OUTPUT_METADATA}", flush=True)

    # Summary statistics
    phyla = set(r["Phylum"] for r in taxonomy_rows)
    domains = set(r["Domain"] for r in taxonomy_rows)
    print(f"\n=== Summary ===", flush=True)
    print(f"Representatives: {len(taxonomy_rows)}", flush=True)
    print(f"Domains: {', '.join(sorted(domains))}", flush=True)
    print(f"Phyla: {len(phyla)}", flush=True)
    for protein in PROTEINS:
        col = f"{protein}_Identity"
        n_hits = sum(1 for m in metadata if m.get(col) not in ("NA", None))
        print(f"{protein} hits: {n_hits}/{len(metadata)}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="GTDB-first tree metadata builder")
    parser.add_argument("--complete-only", action="store_true",
                        help="Filter to Complete Genome / Chromosome assemblies (V3)")
    args = parser.parse_args()

    mode = "V3 (complete genomes)" if args.complete_only else "V2 (all assemblies)"
    print("=" * 60, flush=True)
    print(f"GTDB-first phylogenetic tree metadata builder — {mode}", flush=True)
    print("=" * 60, flush=True)

    # Step 1: Download GTDB taxonomy
    print("\n[1/6] Downloading GTDB taxonomy...", flush=True)
    gtdb_paths = download_gtdb_taxonomy()

    # Step 2: Parse GTDB taxonomy
    print("\n[2/6] Parsing GTDB taxonomy...", flush=True)
    gtdb_records = parse_gtdb_taxonomy(gtdb_paths)
    print(f"  Total GTDB species with complete taxonomy: {len(gtdb_records):,}", flush=True)

    # Step 3: Filter to complete genomes (V3)
    if args.complete_only:
        print("\n[3/6] Filtering to complete genomes...", flush=True)
        gtdb_metadata = load_gtdb_metadata()
        gtdb_records = filter_complete_genomes(gtdb_records, gtdb_metadata)
    else:
        print("\n[3/6] Skipping assembly filter (use --complete-only for V3)...", flush=True)

    # Step 4: Load BLAST data and build genus index
    print("\n[4/6] Loading BLAST data...", flush=True)
    blast_hits = load_blast_hits()
    blast_genus_index = build_blast_genus_index(blast_hits)
    print(f"  Unique genera in BLAST data: {len(blast_genus_index)}", flush=True)

    # Step 5: Select representatives
    print("\n[5/6] Selecting representatives...", flush=True)
    representatives = select_representatives(gtdb_records, blast_genus_index)

    # Step 6: Map BLAST and write output
    print("\n[6/6] Mapping BLAST data and writing output...", flush=True)
    metadata = map_blast_to_reps(representatives, blast_genus_index)
    write_outputs(representatives, metadata)

    print("\nDone! Run the R script to generate the tree:", flush=True)
    print("  Rscript scripts/04_plot_tree.R", flush=True)


if __name__ == "__main__":
    main()
