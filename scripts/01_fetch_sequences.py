#!/usr/bin/env python3
"""Fetch E. coli reference protein sequences from NCBI.

Outputs: data/ecoli_proteins.fasta with 4 proteins:
  - IS621 (user-provided)
  - EF-Tu / TufA (NP_417798.1)
  - SecY (NP_417759.1)
  - MuA (NP_050634.1, phage Mu)
"""

import os
import sys
import time
from pathlib import Path
from Bio import Entrez, SeqIO
from io import StringIO

# NCBI requires an email
Entrez.email = "brady.cress@gmail.com"

OUTFILE = Path("data/ecoli_proteins.fasta")

# User-provided IS621 sequence
IS621_SEQ = (
    "MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEPVAECLYD"
    "AGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRHQALTDMHTQELNR"
    "TETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKTSAVLLAYIGLKDRFAHARQ"
    "FAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGRAFRDRLAANGKKGKVILGAMMRKLAQ"
    "VAYGVLKSGVPFDASRHNPVAA"
)

# Proteins to fetch from NCBI
NCBI_PROTEINS = {
    "EF_Tu":  {"accession": "NP_417798.1", "desc": "elongation factor Tu [Escherichia coli K-12]"},
    "SecY":   {"accession": "NP_417759.1", "desc": "preprotein translocase subunit SecY [Escherichia coli K-12]"},
    "MuA":    {"accession": "NP_050634.1", "desc": "transposase A [Enterobacteria phage Mu]"},
}


def fetch_protein(accession: str, retries: int = 3) -> str:
    """Fetch a protein FASTA from NCBI Entrez."""
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db="protein", id=accession,
                                   rettype="fasta", retmode="text")
            result = handle.read()
            handle.close()
            return result.strip()
        except Exception as e:
            if attempt < retries - 1:
                print(f"  Retry {attempt+1} for {accession}: {e}")
                time.sleep(3)
            else:
                raise


def main():
    if OUTFILE.exists() and OUTFILE.stat().st_size > 0:
        print(f"Cache hit: {OUTFILE} already exists, skipping fetch.")
        return

    OUTFILE.parent.mkdir(parents=True, exist_ok=True)

    records = []

    # Add IS621 (user-provided)
    print("Adding IS621 (user-provided sequence)...")
    records.append(f">Ecoli_IS621 IS621 transposase [Escherichia coli]\n{IS621_SEQ}")

    # Fetch from NCBI
    for name, info in NCBI_PROTEINS.items():
        acc = info["accession"]
        print(f"Fetching {name} ({acc}) from NCBI...")
        fasta = fetch_protein(acc)

        # Replace the NCBI header with our standardized format
        lines = fasta.split("\n")
        seq = "".join(lines[1:])
        records.append(f">Ecoli_{name} {info['desc']}\n{seq}")
        time.sleep(1)  # Be polite to NCBI

    with open(OUTFILE, "w") as f:
        f.write("\n\n".join(records) + "\n")

    print(f"Wrote {len(records)} sequences to {OUTFILE}")


if __name__ == "__main__":
    main()
