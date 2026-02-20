#!/usr/bin/env python3
"""Build tree metadata from BLAST results and NCBI taxonomy.

Reads: data/blast_cache/{protein}_hits.tsv (4 files)
Outputs: data/metadata.csv   (tip metadata for R)
         data/taxonomy.csv   (taxonomy table for tree building in R)

Optimized strategy:
  1. Collect all organisms from BLAST hits
  2. Pre-filter to one organism per genus (reduces taxonomy lookups 10-20x)
  3. Fetch NCBI taxonomy lineages via Entrez (batched)
  4. Select ~100-200 representative species (one per family)
  5. Build metadata CSV with identity values
"""

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
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_blast_hits() -> dict:
    """Load all BLAST hit TSVs. Returns {protein: {organism: identity_pct}}."""
    hits = {}
    for protein in PROTEINS:
        tsv = CACHE_DIR / f"{protein}_hits.tsv"
        if not tsv.exists():
            print(f"WARNING: {tsv} not found, skipping {protein}", flush=True)
            hits[protein] = {}
            continue

        protein_hits = {}
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                org = row["organism"].strip()
                if not org:
                    continue
                # Skip multi-species entries
                if org in ("Bacteria", "Enterobacteriaceae", "unclassified sequences"):
                    continue
                identity = float(row["identity_pct"])
                if org not in protein_hits or identity > protein_hits[org]:
                    protein_hits[org] = identity
        hits[protein] = protein_hits
        print(f"  {protein}: {len(protein_hits)} unique organisms", flush=True)
    return hits


def extract_genus(organism: str) -> str:
    """Extract genus from an organism name (first word of binomial)."""
    parts = organism.strip().split()
    if parts:
        # Handle names like "Candidatus Xyz abc" → "Candidatus Xyz"
        if parts[0] == "Candidatus" and len(parts) > 1:
            return f"{parts[0]} {parts[1]}"
        return parts[0]
    return organism


def pre_filter_by_genus(all_hits: dict, all_organisms: set) -> list[str]:
    """Reduce organisms to one per genus, preferring IS621 hits then EF-Tu."""
    genus_best = {}  # genus → (organism, score)

    for org in all_organisms:
        genus = extract_genus(org)
        is621_id = all_hits.get("IS621", {}).get(org, 0)
        eftu_id = all_hits.get("EF_Tu", {}).get(org, 0)
        score = (is621_id > 0, is621_id, eftu_id)

        if genus not in genus_best or score > genus_best[genus][1]:
            genus_best[genus] = (org, score)

    selected = [org for org, _ in genus_best.values()]
    return selected


def search_taxid(organism: str) -> int | None:
    """Search NCBI taxonomy for an organism name, return taxid."""
    try:
        handle = Entrez.esearch(db="taxonomy", term=f'"{organism}"[Scientific Name]',
                                retmax=1)
        result = Entrez.read(handle)
        handle.close()
        if result["IdList"]:
            return int(result["IdList"][0])
    except Exception:
        pass
    # Fallback: try without quotes and field
    try:
        handle = Entrez.esearch(db="taxonomy", term=organism, retmax=1)
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
            handle = Entrez.efetch(db="taxonomy",
                                   id=",".join(str(t) for t in batch))
            records = Entrez.read(handle)
            handle.close()

            for rec in records:
                taxid = int(rec["TaxId"])
                lineage = {}
                for item in rec.get("LineageEx", []):
                    rank = item["Rank"]
                    if rank in RANKS:
                        lineage[rank] = item["ScientificName"]
                if rec.get("Rank") in RANKS:
                    lineage[rec["Rank"]] = rec["ScientificName"]
                lineages[taxid] = lineage

            time.sleep(0.5)
        except Exception as e:
            print(f"  Error fetching taxonomy batch: {e}", flush=True)
            time.sleep(2)

    return lineages


def load_taxonomy_cache() -> dict:
    if TAXONOMY_CACHE.exists() and TAXONOMY_CACHE.stat().st_size > 0:
        with open(TAXONOMY_CACHE) as f:
            return json.load(f)
    return {}


def save_taxonomy_cache(cache: dict) -> None:
    with open(TAXONOMY_CACHE, "w") as f:
        json.dump(cache, f, indent=2)


def sanitize_tip_label(name: str) -> str:
    """Make a tree-safe tip label from an organism name."""
    label = re.sub(r'[^A-Za-z0-9_]', '_', name)
    label = re.sub(r'_+', '_', label)
    return label.strip('_')[:60]


def main():
    if OUTPUT_METADATA.exists() and OUTPUT_TAXONOMY.exists():
        if OUTPUT_METADATA.stat().st_size > 0 and OUTPUT_TAXONOMY.stat().st_size > 0:
            print(f"Cache hit: {OUTPUT_METADATA} and {OUTPUT_TAXONOMY} already exist.",
                  flush=True)
            return

    print("Loading BLAST hits...", flush=True)
    all_hits = load_blast_hits()

    # Collect all unique organisms
    all_organisms = set()
    for protein_hits in all_hits.values():
        all_organisms.update(protein_hits.keys())
    print(f"\nTotal unique organisms across all BLAST results: {len(all_organisms)}",
          flush=True)

    # Pre-filter: one organism per genus (huge speedup for taxonomy lookups)
    genus_reps = pre_filter_by_genus(all_hits, all_organisms)
    print(f"After genus pre-filter: {len(genus_reps)} organisms to look up",
          flush=True)

    # Load/build taxonomy cache
    print("\nFetching taxonomy...", flush=True)
    tax_cache = load_taxonomy_cache()

    uncached = [org for org in genus_reps if org not in tax_cache]
    print(f"  Cached: {len(tax_cache)}, Need to fetch: {len(uncached)}", flush=True)

    if uncached:
        print(f"  Searching for taxonomy IDs ({len(uncached)} organisms)...", flush=True)
        org_taxids = {}
        failed = 0
        for i, org in enumerate(uncached):
            if i % 25 == 0:
                print(f"    [{i}/{len(uncached)}] Looking up: {org[:50]}...", flush=True)
            taxid = search_taxid(org)
            if taxid:
                org_taxids[org] = taxid
            else:
                failed += 1
            time.sleep(0.35)  # NCBI rate limit

            # Save cache periodically (every 100 lookups)
            if i > 0 and i % 100 == 0:
                # Partial save of what we have so far
                for o, t in list(org_taxids.items()):
                    if o not in tax_cache:
                        tax_cache[o] = {"taxid": t, "lineage": {}}
                save_taxonomy_cache(tax_cache)

        print(f"  Found taxids for {len(org_taxids)}/{len(uncached)} organisms "
              f"({failed} failed)", flush=True)

        # Fetch lineages in batches
        if org_taxids:
            all_taxids = list(set(org_taxids.values()))
            print(f"  Fetching {len(all_taxids)} lineages in batches...", flush=True)
            lineages = fetch_taxonomy_batch(all_taxids)

            for org, taxid in org_taxids.items():
                if taxid in lineages:
                    tax_cache[org] = {
                        "taxid": taxid,
                        "lineage": lineages[taxid],
                    }

        save_taxonomy_cache(tax_cache)
        print(f"  Taxonomy cache updated: {len(tax_cache)} entries", flush=True)

    # Now re-expand: for ALL organisms (not just genus reps), assign taxonomy
    # from the genus representative's cached lineage
    print("\nMapping all organisms to taxonomy via genus representatives...", flush=True)
    genus_to_lineage = {}
    for org in genus_reps:
        if org in tax_cache:
            genus = extract_genus(org)
            genus_to_lineage[genus] = tax_cache[org].get("lineage", {})

    # Build org_taxonomy for all organisms
    org_taxonomy = {}
    for org in all_organisms:
        if org in tax_cache:
            org_taxonomy[org] = tax_cache[org].get("lineage", {})
        else:
            genus = extract_genus(org)
            if genus in genus_to_lineage:
                org_taxonomy[org] = genus_to_lineage[genus]

    print(f"  Taxonomy available for {len(org_taxonomy)}/{len(all_organisms)} organisms",
          flush=True)

    # Select representatives: one per family
    print("\nSelecting representative species...", flush=True)
    family_reps = defaultdict(list)

    for org, lineage in org_taxonomy.items():
        family = lineage.get("family", "Unknown")
        if family == "Unknown":
            continue
        # Skip eukaryotes (mitochondrial/chloroplast BLAST hits)
        eukaryotic_phyla = {"Chordata", "Arthropoda", "Nematoda", "Mollusca",
                            "Cnidaria", "Streptophyta", "Chlorophyta", "Ascomycota",
                            "Basidiomycota", "Apicomplexa", "Euglenozoa"}
        phylum = lineage.get("phylum", "")
        if phylum in eukaryotic_phyla:
            continue

        is621_id = all_hits.get("IS621", {}).get(org, 0)
        eftu_id = all_hits.get("EF_Tu", {}).get(org, 0)
        score = (is621_id > 0, is621_id, eftu_id)
        family_reps[family].append((org, score, lineage))

    # Pick best per family
    representatives = []
    for family, candidates in family_reps.items():
        candidates.sort(key=lambda x: x[1], reverse=True)
        best_org, _, lineage = candidates[0]
        representatives.append((best_org, lineage))

    # Force E. coli as the Enterobacteriaceae representative (it's the reference)
    ecoli_org = None
    for org in all_organisms:
        if org == "Escherichia coli" and org in org_taxonomy:
            ecoli_org = org
            break
    if ecoli_org is None:
        # Try partial match
        for org in all_organisms:
            if "Escherichia coli" in org and org in org_taxonomy:
                ecoli_org = org
                break

    if ecoli_org:
        ecoli_lineage = org_taxonomy[ecoli_org]
        # Remove any existing Enterobacteriaceae representative
        representatives = [(o, l) for o, l in representatives
                          if l.get("family") != "Enterobacteriaceae"]
        representatives.append((ecoli_org, ecoli_lineage))
        print(f"  Forced E. coli ({ecoli_org}) as Enterobacteriaceae representative",
              flush=True)
    else:
        print("  WARNING: Could not find E. coli in BLAST results!", flush=True)

    print(f"Selected {len(representatives)} representatives from {len(family_reps)} families",
          flush=True)

    # Trim to ~200 max if needed
    if len(representatives) > 200:
        has_is621 = [(o, l) for o, l in representatives
                     if all_hits.get("IS621", {}).get(o, 0) > 0]
        no_is621 = [(o, l) for o, l in representatives
                    if all_hits.get("IS621", {}).get(o, 0) == 0]

        phylum_groups = defaultdict(list)
        for o, l in no_is621:
            phylum_groups[l.get("phylum", "Unknown")].append((o, l))

        sampled = []
        target = 180 - len(has_is621)
        per_phylum = max(2, target // max(1, len(phylum_groups)))
        for phylum, members in sorted(phylum_groups.items()):
            sampled.extend(members[:per_phylum])

        representatives = has_is621 + sampled[:target]
        print(f"  Trimmed to {len(representatives)} representatives", flush=True)

    # Build output tables
    print("\nBuilding output files...", flush=True)
    taxonomy_rows = []
    metadata_rows = []
    seen_tip_ids = set()

    for org, lineage in representatives:
        tip_id = sanitize_tip_label(org)
        # Ensure unique tip IDs
        if tip_id in seen_tip_ids:
            tip_id = tip_id + "_2"
        seen_tip_ids.add(tip_id)

        tax_row = {
            "Tip_ID": tip_id,
            "Domain": lineage.get("superkingdom", "Bacteria"),
            "Phylum": lineage.get("phylum", "Unknown"),
            "Class": lineage.get("class", "Unknown"),
            "Order": lineage.get("order", "Unknown"),
            "Family": lineage.get("family", "Unknown"),
            "Genus": lineage.get("genus", "Unknown"),
            "Species": tip_id,
        }
        taxonomy_rows.append(tax_row)

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
    tax_fields = ["Tip_ID", "Domain", "Phylum", "Class", "Order", "Family",
                  "Genus", "Species"]
    with open(OUTPUT_TAXONOMY, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tax_fields)
        writer.writeheader()
        writer.writerows(taxonomy_rows)
    print(f"  Wrote {len(taxonomy_rows)} rows to {OUTPUT_TAXONOMY}", flush=True)

    # Write metadata CSV
    meta_fields = ["Tip_ID", "Phylum", "Ecoli_Reference",
                   "EF_Tu_Identity", "SecY_Identity", "MuA_Identity", "IS621_Identity"]
    with open(OUTPUT_METADATA, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=meta_fields)
        writer.writeheader()
        writer.writerows(metadata_rows)
    print(f"  Wrote {len(metadata_rows)} rows to {OUTPUT_METADATA}", flush=True)

    # Summary
    phyla = set(r["Phylum"] for r in taxonomy_rows)
    print(f"\n=== Summary ===", flush=True)
    print(f"Representatives: {len(representatives)}", flush=True)
    print(f"Phyla covered: {len(phyla)}", flush=True)
    for p in sorted(phyla):
        n = sum(1 for r in taxonomy_rows if r["Phylum"] == p)
        print(f"  {p}: {n}", flush=True)

    n_is621 = sum(1 for r in metadata_rows if r["IS621_Identity"] != "NA")
    n_mua = sum(1 for r in metadata_rows if r["MuA_Identity"] != "NA")
    print(f"\nIS621 hits: {n_is621}/{len(metadata_rows)}", flush=True)
    print(f"MuA hits: {n_mua}/{len(metadata_rows)}", flush=True)


if __name__ == "__main__":
    main()
