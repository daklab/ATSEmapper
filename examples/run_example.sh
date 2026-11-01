#!/bin/bash
# run_example.sh — end-to-end ATSEmapper demo on 20 chr19 cells
#
# Run from the repo root:
#   bash examples/run_example.sh
#
# What it does:
#   1. Checks that genome reference files exist (runs setup if needed)
#   2. Runs the Python demo (junction parsing → QC → ATSE mapping)
#   3. Runs the unit test suite to confirm everything is working

set -euo pipefail
REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO_ROOT"

SAMPLE_DIR="examples/sample_data"
GTF="$SAMPLE_DIR/chr19.gtf"
FASTA="$SAMPLE_DIR/chr19.fa"

echo "================================================"
echo "ATSEmapper — example run"
echo "================================================"
echo "Repo: $REPO_ROOT"
echo ""

# ── Step 1: genome reference files ────────────────────────────────────────────
if [[ ! -f "$GTF" || ! -f "$FASTA" ]]; then
    echo "Genome reference files not found. Running setup..."
    bash examples/setup_sample_data.sh
    echo ""
fi

echo "Reference files:"
ls -lh "$GTF" "$FASTA"
echo ""

# ── Step 2: unit tests ────────────────────────────────────────────────────────
echo "Running unit tests..."
python -m pytest tests/test_atsemapper.py -v --tb=short
echo ""

# ── Step 3: demo pipeline ─────────────────────────────────────────────────────
echo "Running demo pipeline..."
python examples/atsemapper_demo.py
echo ""

echo "================================================"
echo "Done. Check examples/sample_data/ for outputs."
echo "================================================"
