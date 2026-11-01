"""
ATSEmapper demo — chr19 subset of 20 Tabula Muris Senis Smart-seq2 cells.

Before running this script, generate the genome reference files:
    GENOME_GTF=... GENOME_FASTA=... bash examples/setup_sample_data.sh   # from local genome files
    bash examples/setup_sample_data.sh --download  # anywhere with internet

Then run from the repo root:
    python examples/atsemapper_demo.py
"""

import os
import sys
import pandas as pd
from atsemapper.atsemapper.main import run_atsemapper

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "sample_data")
GTF_FILE    = os.path.join(EXAMPLES_DIR, "chr19.gtf")
FASTA_FILE  = os.path.join(EXAMPLES_DIR, "chr19.fa")
JUNC_LIST   = os.path.join(EXAMPLES_DIR, "junction_files.txt")

# ── Pre-flight checks ──────────────────────────────────────────────────────────
missing = [p for p in [GTF_FILE, FASTA_FILE, JUNC_LIST] if not os.path.exists(p)]
if missing:
    print("ERROR: missing required files:")
    for p in missing:
        print(f"  {p}")
    print()
    if GTF_FILE in missing or FASTA_FILE in missing:
        print("Generate genome reference files first:")
        print("  bash examples/setup_sample_data.sh")
    sys.exit(1)

with open(JUNC_LIST) as f:
    n_cells = sum(1 for l in f if l.strip())
print(f"Demo: {n_cells} cells, chr19 only")
print(f"  Junction list: {JUNC_LIST}")
print(f"  GTF:   {GTF_FILE}")
print(f"  FASTA: {FASTA_FILE}")
print()


# ── Configure ATSEmapper ───────────────────────────────────────────────────────
class Args:
    def __init__(self):
        # Input: text file listing one junction BED path per line
        self.input      = JUNC_LIST
        self.output     = os.path.join(os.path.dirname(__file__), "output")   # examples/output (git-ignored)
        self.annotation = GTF_FILE
        self.genome     = FASTA_FILE
        self.db_path    = None          # created automatically

        # Junction filters — relaxed for a small chr19 demo
        self.min_intron = 50
        self.max_intron = 500000
        self.min_reads  = 5             # low because only 20 cells
        self.min_cells  = 2

        # Processing
        self.batch_size  = 5
        self.num_workers = 2

        # Single-cell Smart-seq2 data
        self.sequencing_type   = "single_cell"
        self.annotation_status = "both"
        self.only_canonical    = True
        self.tolerance         = 1
        self.min_splice_site_usage = 0.01

        self.sample_size = None
        self.verbose     = True
        self.log_file    = None


# ── Run ───────────────────────────────────────────────────────────────────────
print("Running ATSEmapper...")
atse_file = run_atsemapper(Args())

# ── Summarise results ─────────────────────────────────────────────────────────
atse_df = pd.read_csv(atse_file, sep="\t", compression="gzip")
n_events = atse_df["event_id"].nunique()
n_junctions = len(atse_df)

print()
print(f"Results saved to: {atse_file}")
print(f"  {n_events} ATSE events  ({n_junctions} junction-event rows)")
print()
print("Top events by junction count:")
top = (atse_df.groupby("event_id")
       .size()
       .sort_values(ascending=False)
       .head(10)
       .reset_index(name="n_junctions"))
print(top.to_string(index=False))
print()
print("Sample rows:")
print(atse_df.head(5).to_string(index=False))
