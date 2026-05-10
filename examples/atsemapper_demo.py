import os
import pandas as pd
from atsemapper.atsemapper.main import run_atsemapper

# Paths are relative to the repo root. Run this script from there:
#   python examples/atsemapper_demo.py
EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "sample_data")

gtf_file      = os.path.join(EXAMPLES_DIR, "chr19.gtf")
fasta_file    = os.path.join(EXAMPLES_DIR, "chr19.fa")
junctions_dir = os.path.join(EXAMPLES_DIR, "junctions")

# =================================================================
# Configure parameters
# =================================================================
class Args:
    def __init__(self):
        self.input        = junctions_dir
        self.output       = None            # auto-generated with timestamp
        self.annotation   = gtf_file
        self.genome       = fasta_file
        self.db_path      = None            # created automatically if missing

        self.min_intron   = 50
        self.max_intron   = 500000
        self.min_reads    = 5
        self.min_cells    = 2

        self.batch_size   = 10
        self.num_workers  = 4

        self.sequencing_type    = "bulk"
        self.annotation_status  = "both"
        self.only_canonical     = True

        self.min_splice_site_usage = 0.01
        self.sample_size  = None
        self.tolerance    = 1
        self.verbose      = False
        self.log_file     = None

# =================================================================
# Run
# =================================================================
print("Starting ATSEmapper analysis...")
atse_file = run_atsemapper(Args())
print(f"Results saved to: {atse_file}")

atse_df = pd.read_csv(atse_file, sep="\t", compression="gzip")
print(f"\n{len(atse_df)} alternative splicing events found.")
print(atse_df.head())
