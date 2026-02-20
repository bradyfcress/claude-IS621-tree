# IS621 Phylogenetic Distribution Across Prokaryotes

Publication-quality circular phylogenetic tree comparing the evolutionary distribution
of four proteins across the prokaryotic tree of life, with *E. coli* as the reference.

## What This Shows

Four concentric heatmap rings show % amino acid identity to *E. coli* orthologs:

| Ring | Protein | Distribution | Biological Role |
|------|---------|-------------|-----------------|
| 1 (inner) | **EF-Tu** | Universal | Translation elongation factor |
| 2 | **SecY** | Universal | Protein translocation channel |
| 3 | **MuA** | Patchy | Phage Mu transposase |
| 4 (outer) | **IS621** | Very sparse | IS element transposase |

The contrast between the fully-blue inner rings (universal proteins) and the
mostly-grey outer rings (mobile elements) illustrates how insertion sequences
have a fundamentally different phylogenetic distribution than core cellular machinery.

## Quick Start

```bash
# Generate figure with mock data (instant, for testing)
bash scripts/run_pipeline.sh mock

# Generate figure with real NCBI BLAST data (20-60 min)
bash scripts/run_pipeline.sh real

# Re-run visualization only (after data is generated)
bash scripts/run_pipeline.sh plot
```

Output: `output/IS621_phylo_tree.pdf` (vector) and `.png` (raster)

## Pipeline Steps

1. **01_fetch_sequences.py** — Fetch *E. coli* reference proteins from NCBI
2. **02_blast_search.py** — BLAST each protein against NCBI refseq_protein
3. **03_build_metadata.py** — Parse BLAST hits, fetch taxonomy, select representatives
4. **04_plot_tree.R** — Render circular tree with `ggtree` + `gheatmap()`

## Dependencies

**Python 3.11+**: `biopython`
**R 4.4+**: `ggtree`, `ape`, `ggplot2`, `ggnewscale`, `treeio`, `dplyr`, `tidyr`

## Reference Sequences

| Protein | Accession | Length |
|---------|-----------|--------|
| IS621 | User-provided | 326 aa |
| EF-Tu (TufA) | NP_417798.1 | 394 aa |
| SecY | NP_417759.1 | 443 aa |
| MuA | NP_050634.1 | 663 aa |
