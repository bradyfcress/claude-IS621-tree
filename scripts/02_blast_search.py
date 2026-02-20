#!/usr/bin/env python3
"""BLAST E. coli reference proteins against NCBI refseq_protein.

Reads: data/ecoli_proteins.fasta
Outputs: data/blast_cache/{protein}_blast.xml  (raw BLAST XML)
         data/blast_cache/{protein}_hits.tsv    (parsed: organism, identity, accession)

Caches results — skips any protein whose XML file already exists.
"""

import os
import sys
import time
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

FASTA_FILE = Path("data/ecoli_proteins.fasta")
CACHE_DIR = Path("data/blast_cache")


def blast_protein(name: str, sequence: str, xml_path: Path) -> None:
    """Run BLAST against refseq_protein and cache the XML result."""
    print(f"  BLASTing {name} ({len(sequence)} aa) against refseq_protein...")
    print(f"  This may take 5-15 minutes...")
    start = time.time()

    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="refseq_protein",
        sequence=sequence,
        hitlist_size=5000,
        expect=1e-3,
        format_type="XML",
    )

    xml_data = result_handle.read()
    result_handle.close()

    with open(xml_path, "w") as f:
        f.write(xml_data)

    elapsed = time.time() - start
    print(f"  Done in {elapsed:.0f}s. Saved to {xml_path}")


def parse_blast_xml(xml_path: Path, tsv_path: Path) -> None:
    """Parse BLAST XML into a TSV of (organism, identity, accession, evalue)."""
    with open(xml_path) as f:
        records = NCBIXML.parse(f)
        record = next(records)

    hits = []
    for alignment in record.alignments:
        # Extract organism from hit description
        desc = alignment.hit_def
        # Organism is typically in brackets at the end
        org = ""
        if "[" in desc and "]" in desc:
            org = desc[desc.rfind("[") + 1 : desc.rfind("]")]

        # Best HSP
        best_hsp = max(alignment.hsps, key=lambda h: h.identities)
        identity_pct = 100.0 * best_hsp.identities / best_hsp.align_length
        evalue = best_hsp.expect

        accession = alignment.accession

        hits.append({
            "organism": org,
            "identity_pct": round(identity_pct, 2),
            "accession": accession,
            "evalue": evalue,
            "hit_def": desc[:200],  # truncate long descriptions
        })

    # Write TSV
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["organism", "identity_pct", "accession", "evalue", "hit_def"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(hits)

    print(f"  Parsed {len(hits)} hits → {tsv_path}")


def main():
    if not FASTA_FILE.exists():
        print(f"ERROR: {FASTA_FILE} not found. Run 01_fetch_sequences.py first.")
        sys.exit(1)

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Read all reference sequences
    sequences = {}
    for record in SeqIO.parse(FASTA_FILE, "fasta"):
        # ID format: Ecoli_IS621, Ecoli_EF_Tu, Ecoli_SecY, Ecoli_MuA
        name = record.id.replace("Ecoli_", "")
        sequences[name] = str(record.seq)

    print(f"Loaded {len(sequences)} reference sequences: {list(sequences.keys())}")

    for name, seq in sequences.items():
        xml_path = CACHE_DIR / f"{name}_blast.xml"
        tsv_path = CACHE_DIR / f"{name}_hits.tsv"

        # Check cache
        if xml_path.exists() and xml_path.stat().st_size > 0:
            print(f"Cache hit: {xml_path} exists, skipping BLAST for {name}.")
        else:
            blast_protein(name, seq, xml_path)
            time.sleep(5)  # Be polite between searches

        # Always re-parse (fast) to ensure TSV is up to date
        if xml_path.exists():
            parse_blast_xml(xml_path, tsv_path)


if __name__ == "__main__":
    main()
