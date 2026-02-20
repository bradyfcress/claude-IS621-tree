#!/usr/bin/env python3
"""Build tree metadata from BLAST results and NCBI taxonomy.

Reads: data/blast_cache/{protein}_hits.tsv (4 files)
Outputs: data/metadata.csv   (tip metadata for R)
         data/taxonomy.csv   (taxonomy table for tree building in R)

Strategy:
  1. Collect all organisms from BLAST hits
  2. Fetch NCBI taxonomy lineages via Entrez
  3. Select ~100-150 representative species (one per family)
  4. Build metadata CSV with identity values
"""

import os
import sys
import csv
import json
import time
import re
from pathlib import Path
from collections import defaultdict
from Bio import Entrez

Entrez.email = "brady.cress@gmail.com"

CACHE_DIR = Path("data/blast_cache")
TAXONOMY_CACHE = CACHE_DIR / "taxonomy_cache.json"
OUTPUT_METADATA = Path("data/metadata.csv")
OUTPUT_TAXONOMY = Path("data/taxonomy.csv")

PROTEINS = ["IS621", "EF_Tu", "SecY", "MuA"]

# Desired taxonomy ranks
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_MAP = {"superkingdom": "Domain"}


def load_blast_hits() -> dict:
    """Load all BLAST hit TSVs. Returns {protein: {organism: identity_pct}}."""
    hits = {}
    for protein in PROTEINS:
        tsv = CACHE_DIR / f"{protein}_hits.tsv"
        if not tsv.exists():
            print(f"WARNING: {tsv} not found, skipping {protein}")
            hits[protein] = {}
            continue

        protein_hits = {}
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                org = row["organism"].strip()
                if not org:
                    continue
                identity = float(row["identity_pct"])
                # Keep best hit per organism
                if org not in protein_hits or identity > protein_hits[org]:
                    protein_hits[org] = identity
        hits[protein] = protein_hits
        print(f"  {protein}: {len(protein_hits)} unique organisms")
    return hits


def search_taxid(organism: str) -> int | None:
    """Search NCBI taxonomy for an organism name, return taxid."""
    try:
        handle = Entrez.esearch(db="taxonomy", term=f'"{organism}"[Scientific Name]')
        result = Entrez.read(handle)
        handle.close()
        if result["IdList"]:
            return int(result["IdList"][0])
    except Exception:
        pass
    return None


def fetch_taxonomy_batch(taxids: list[int]) -> dict:
    """Fetch taxonomy lineages for a batch of taxids."""
    lineages = {}
    batch_size = 200
    for i in range(0, len(taxids), batch_size):
        batch = taxids[i:i + batch_size]
        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(str(t) for t in batch))
            records = Entrez.read(handle)
            handle.close()

            for rec in records:
                taxid = int(rec["TaxId"])
                lineage = {}
                for item in rec.get("LineageEx", []):
                    rank = item["Rank"]
                    if rank in RANKS:
                        lineage[rank] = item["ScientificName"]
                # Add the organism's own rank
                if rec.get("Rank") in RANKS:
                    lineage[rec["Rank"]] = rec["ScientificName"]
                lineages[taxid] = lineage

            time.sleep(0.5)
        except Exception as e:
            print(f"  Error fetching taxonomy batch: {e}")
            time.sleep(2)

    return lineages


def load_taxonomy_cache() -> dict:
    """Load cached taxonomy data."""
    if TAXONOMY_CACHE.exists() and TAXONOMY_CACHE.stat().st_size > 0:
        with open(TAXONOMY_CACHE) as f:
            return json.load(f)
    return {}


def save_taxonomy_cache(cache: dict) -> None:
    """Save taxonomy cache."""
    with open(TAXONOMY_CACHE, "w") as f:
        json.dump(cache, f, indent=2)


def sanitize_tip_label(name: str) -> str:
    """Make a tree-safe tip label from an organism name."""
    # Replace spaces and special chars with underscores
    label = re.sub(r'[^A-Za-z0-9_]', '_', name)
    label = re.sub(r'_+', '_', label)
    return label.strip('_')[:60]


def main():
    # Check if outputs already exist
    if OUTPUT_METADATA.exists() and OUTPUT_TAXONOMY.exists():
        meta_size = OUTPUT_METADATA.stat().st_size
        tax_size = OUTPUT_TAXONOMY.stat().st_size
        if meta_size > 0 and tax_size > 0:
            print(f"Cache hit: {OUTPUT_METADATA} and {OUTPUT_TAXONOMY} already exist.")
            return

    print("Loading BLAST hits...")
    all_hits = load_blast_hits()

    # Collect all unique organisms
    all_organisms = set()
    for protein_hits in all_hits.values():
        all_organisms.update(protein_hits.keys())
    print(f"\nTotal unique organisms across all BLAST results: {len(all_organisms)}")

    # Load/build taxonomy cache
    print("\nFetching taxonomy...")
    tax_cache = load_taxonomy_cache()
    # Cache format: {"organism_name": {"taxid": X, "lineage": {...}}}

    # Find organisms not yet cached
    uncached = [org for org in all_organisms if org not in tax_cache]
    print(f"  Cached: {len(tax_cache)}, Need to fetch: {len(uncached)}")

    if uncached:
        # First, search for taxids
        print(f"  Searching for taxonomy IDs (this may take a while)...")
        org_taxids = {}
        for i, org in enumerate(uncached):
            if i % 50 == 0 and i > 0:
                print(f"    {i}/{len(uncached)} searched...")
                time.sleep(1)
            taxid = search_taxid(org)
            if taxid:
                org_taxids[org] = taxid
            time.sleep(0.4)  # Rate limit

        print(f"  Found taxids for {len(org_taxids)}/{len(uncached)} organisms")

        # Fetch lineages in batches
        if org_taxids:
            all_taxids = list(set(org_taxids.values()))
            print(f"  Fetching {len(all_taxids)} lineages...")
            lineages = fetch_taxonomy_batch(all_taxids)

            # Map back to organism names
            for org, taxid in org_taxids.items():
                if taxid in lineages:
                    tax_cache[org] = {
                        "taxid": taxid,
                        "lineage": lineages[taxid],
                    }

        save_taxonomy_cache(tax_cache)
        print(f"  Taxonomy cache updated: {len(tax_cache)} entries")

    # Select representatives: one per family
    print("\nSelecting representative species...")
    family_reps = defaultdict(list)  # family -> [(org, identity_scores)]

    for org in all_organisms:
        if org not in tax_cache:
            continue
        lineage = tax_cache[org].get("lineage", {})
        family = lineage.get("family", "Unknown")
        if family == "Unknown":
            continue

        # Score: prefer species with IS621 hits, then highest EF-Tu identity
        is621_id = all_hits.get("IS621", {}).get(org, 0)
        eftu_id = all_hits.get("EF_Tu", {}).get(org, 0)
        score = (is621_id > 0, is621_id, eftu_id)

        family_reps[family].append((org, score, lineage))

    # Pick the best representative per family
    representatives = []
    for family, candidates in family_reps.items():
        candidates.sort(key=lambda x: x[1], reverse=True)
        best_org, _, lineage = candidates[0]
        representatives.append((best_org, lineage))

    # Ensure E. coli is included
    ecoli_included = any("Escherichia coli" in org for org, _ in representatives)
    if not ecoli_included:
        for org in all_organisms:
            if "Escherichia coli" in org:
                if org in tax_cache:
                    representatives.append((org, tax_cache[org].get("lineage", {})))
                    break

    print(f"Selected {len(representatives)} representatives from {len(family_reps)} families")

    # Limit to ~150 if too many (keep diverse set)
    if len(representatives) > 200:
        # Keep all IS621-containing reps, sample the rest
        has_is621 = [(o, l) for o, l in representatives
                     if all_hits.get("IS621", {}).get(o, 0) > 0]
        no_is621 = [(o, l) for o, l in representatives
                    if all_hits.get("IS621", {}).get(o, 0) == 0]

        # Ensure phylum diversity in the non-IS621 set
        phylum_groups = defaultdict(list)
        for o, l in no_is621:
            phylum_groups[l.get("phylum", "Unknown")].append((o, l))

        sampled = []
        target = 150 - len(has_is621)
        per_phylum = max(1, target // max(1, len(phylum_groups)))
        for phylum, members in sorted(phylum_groups.items()):
            sampled.extend(members[:per_phylum])

        representatives = has_is621 + sampled[:target]
        print(f"  Trimmed to {len(representatives)} representatives")

    # Build output tables
    print("\nBuilding output files...")
    taxonomy_rows = []
    metadata_rows = []

    for org, lineage in representatives:
        tip_id = sanitize_tip_label(org)

        # Taxonomy
        tax_row = {
            "Tip_ID": tip_id,
            "Domain": lineage.get("superkingdom", "Bacteria"),
            "Phylum": lineage.get("phylum", "Unknown"),
            "Class": lineage.get("class", "Unknown"),
            "Order": lineage.get("order", "Unknown"),
            "Family": lineage.get("family", "Unknown"),
            "Genus": lineage.get("genus", "Unknown"),
            "Species": tip_id,  # Use sanitized label as species-level ID
        }
        taxonomy_rows.append(tax_row)

        # Metadata
        is_ecoli = "Escherichia coli" in org
        meta_row = {
            "Tip_ID": tip_id,
            "Phylum": lineage.get("phylum", "Unknown"),
            "Ecoli_Reference": "TRUE" if is_ecoli else "FALSE",
        }
        for protein in PROTEINS:
            val = all_hits.get(protein, {}).get(org, None)
            col = f"{protein}_Identity"
            if val is not None and val > 0:
                meta_row[col] = round(val, 1)
            else:
                meta_row[col] = "NA"

        metadata_rows.append(meta_row)

    # Write taxonomy CSV
    OUTPUT_TAXONOMY.parent.mkdir(parents=True, exist_ok=True)
    tax_fields = ["Tip_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    with open(OUTPUT_TAXONOMY, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tax_fields)
        writer.writeheader()
        writer.writerows(taxonomy_rows)
    print(f"  Wrote {len(taxonomy_rows)} rows to {OUTPUT_TAXONOMY}")

    # Write metadata CSV
    meta_fields = ["Tip_ID", "Phylum", "Ecoli_Reference",
                   "EF_Tu_Identity", "SecY_Identity", "MuA_Identity", "IS621_Identity"]
    with open(OUTPUT_METADATA, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=meta_fields)
        writer.writeheader()
        writer.writerows(metadata_rows)
    print(f"  Wrote {len(metadata_rows)} rows to {OUTPUT_METADATA}")

    # Summary
    phyla = set(r["Phylum"] for r in taxonomy_rows)
    print(f"\n=== Summary ===")
    print(f"Representatives: {len(representatives)}")
    print(f"Phyla covered: {len(phyla)}")
    for p in sorted(phyla):
        n = sum(1 for r in taxonomy_rows if r["Phylum"] == p)
        print(f"  {p}: {n}")


if __name__ == "__main__":
    main()
