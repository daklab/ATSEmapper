#!/bin/bash
# setup_sample_data.sh — generate chr19 genome reference files for the ATSEmapper demo.
#
# The junction BED files are already in examples/sample_data/ (chr19 only, 20 cells).
# This script generates the two large reference files that are NOT tracked in git:
#   - examples/sample_data/chr19.gtf   (~18 MB)
#   - examples/sample_data/chr19.fa    (~60 MB)
#
# Usage (with a local full-genome GENCODE vM19 GTF and FASTA):
#   GENOME_GTF=/path/genes.gtf GENOME_FASTA=/path/genome.fa bash examples/setup_sample_data.sh
#
# Usage (on any machine with internet access):
#   bash examples/setup_sample_data.sh --download
#
# After running this, the demo is fully functional:
#   python examples/atsemapper_demo.py
#   # or
#   bash examples/run_example.sh

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SAMPLE_DIR="$SCRIPT_DIR/sample_data"
mkdir -p "$SAMPLE_DIR"

DOWNLOAD=false
for arg in "$@"; do
    [[ "$arg" == "--download" ]] && DOWNLOAD=true
done

# ── Option A: a local full-genome GENCODE vM19 GTF and FASTA (set GENOME_GTF / GENOME_FASTA) ──
LOCAL_GTF=${GENOME_GTF:-}
LOCAL_FASTA=${GENOME_FASTA:-}

if [[ -n "$LOCAL_GTF" && -f "$LOCAL_GTF" && "$DOWNLOAD" == "false" ]]; then
    echo "=== Using local genome files ==="

    echo "Extracting chr19 from GTF..."
    grep "^chr19	" "$LOCAL_GTF" > "$SAMPLE_DIR/chr19.gtf"
    echo "  -> $SAMPLE_DIR/chr19.gtf ($(wc -l < "$SAMPLE_DIR/chr19.gtf") lines)"

    if command -v samtools &>/dev/null; then
        echo "Extracting chr19 FASTA with samtools..."
        samtools faidx "$LOCAL_FASTA" chr19 > "$SAMPLE_DIR/chr19.fa"
        samtools faidx "$SAMPLE_DIR/chr19.fa"
        echo "  -> $SAMPLE_DIR/chr19.fa ($(du -h "$SAMPLE_DIR/chr19.fa" | cut -f1))"
    else
        echo "WARNING: samtools not found — cannot extract chr19.fa"
        echo "  Install samtools or use --download to fetch from UCSC"
    fi

# ── Option B: Download from UCSC/Ensembl ──────────────────────────────────────
else
    echo "=== Downloading chr19 genome files from public sources ==="

    if ! command -v wget &>/dev/null && ! command -v curl &>/dev/null; then
        echo "ERROR: wget or curl required for download mode"
        exit 1
    fi

    DL() { command -v wget &>/dev/null && wget -q -O "$2" "$1" || curl -sL -o "$2" "$1"; }

    # GENCODE vM25 chr19 FASTA (mm10/GRCm38)
    echo "Downloading chr19 FASTA from UCSC (mm10)..."
    DL "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz" \
       "$SAMPLE_DIR/chr19.fa.gz"
    gunzip -f "$SAMPLE_DIR/chr19.fa.gz"
    echo "  -> $SAMPLE_DIR/chr19.fa"

    # GENCODE vM25 chr19 GTF
    echo "Downloading GENCODE vM25 GTF..."
    DL "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" \
       "$SAMPLE_DIR/gencode.vM25.full.gtf.gz"
    echo "Extracting chr19..."
    zcat "$SAMPLE_DIR/gencode.vM25.full.gtf.gz" | grep "^chr19	" > "$SAMPLE_DIR/chr19.gtf"
    rm "$SAMPLE_DIR/gencode.vM25.full.gtf.gz"
    echo "  -> $SAMPLE_DIR/chr19.gtf ($(wc -l < "$SAMPLE_DIR/chr19.gtf") lines)"

    # Index FASTA
    if command -v samtools &>/dev/null; then
        samtools faidx "$SAMPLE_DIR/chr19.fa"
    fi
fi

echo ""
echo "=== Setup complete ==="
echo "Files ready:"
ls -lh "$SAMPLE_DIR/chr19.gtf" "$SAMPLE_DIR/chr19.fa" 2>/dev/null || true
echo ""
echo "Run the demo:"
echo "  python examples/atsemapper_demo.py"
echo "  # or"
echo "  bash examples/run_example.sh"
