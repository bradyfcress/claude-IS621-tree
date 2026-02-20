#!/bin/bash
# ============================================================================
# run_pipeline.sh — Orchestrate the IS621 phylogenetic tree pipeline
#
# Usage:
#   bash scripts/run_pipeline.sh mock    # Use mock data (fast, for testing)
#   bash scripts/run_pipeline.sh real    # Use real NCBI BLAST data (V1: BLAST-first)
#   bash scripts/run_pipeline.sh v2      # Use GTDB taxonomy backbone (V2: tree-first)
#   bash scripts/run_pipeline.sh plot    # Only re-run visualization
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

MODE="${1:-mock}"

echo "=== IS621 Phylogenetic Tree Pipeline ==="
echo "Mode: $MODE"
echo "Working directory: $PROJECT_DIR"
echo ""

if [ "$MODE" = "mock" ]; then
    echo "--- Generating mock data ---"
    python3 scripts/generate_mock_data.py
    echo ""

elif [ "$MODE" = "real" ]; then
    echo "--- Step 1: Fetching E. coli reference sequences ---"
    python3 scripts/01_fetch_sequences.py
    echo ""

    echo "--- Step 2: Running BLAST searches (this may take 20-60 minutes) ---"
    python3 scripts/02_blast_search.py
    echo ""

    echo "--- Step 3: Building tree metadata ---"
    python3 scripts/03_build_metadata.py
    echo ""

elif [ "$MODE" = "v2" ]; then
    echo "--- Step 1: Fetching E. coli reference sequences ---"
    python3 scripts/01_fetch_sequences.py
    echo ""

    echo "--- Step 2: Running BLAST searches (uses cache if available) ---"
    python3 scripts/02_blast_search.py
    echo ""

    echo "--- Step 3v2: Building GTDB-first tree metadata ---"
    python3 scripts/03v2_build_metadata_gtdb.py
    echo ""

elif [ "$MODE" = "plot" ]; then
    echo "--- Skipping data generation, re-running visualization only ---"
    if [ ! -f data/metadata.csv ] || [ ! -f data/taxonomy.csv ]; then
        echo "ERROR: data/metadata.csv or data/taxonomy.csv not found."
        echo "Run with 'mock' or 'real' first to generate data."
        exit 1
    fi

else
    echo "ERROR: Unknown mode '$MODE'. Use 'mock', 'real', 'v2', or 'plot'."
    exit 1
fi

echo "--- Step 4: Generating visualization ---"
Rscript scripts/04_plot_tree.R
echo ""

echo "=== Pipeline complete ==="
echo "Output: output/IS621_phylo_tree.pdf"
echo "Output: output/IS621_phylo_tree.png"
