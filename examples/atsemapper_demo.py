# Import necessary libraries
import pandas as pd
import os
import sys

# Add the parent directory to the path to enable imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import LeafletFA utils module for ATSEmapper
from leafletfa_utils.atsemapper.main import run_atsemapper

# =================================================================
# STEP 1: Set up file paths
# =================================================================
# Path to the gene annotation file (GTF format)
# This file contains information about gene structure, exons, etc.
gtf_file = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/genome_files/gencode.vM19/genes/genes.gtf"

# Path to the genome reference file (FASTA format)
# This contains the actual genomic sequence
fasta_file = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/genome_files/gencode.vM19/fasta/genome.fa"

# Provide existing database (annotation file) for ATSEmapper 
# ADD EXAMPLE CODE FOR HOW TO CREATE THIS FILE... 
db_path = "/gpfs/commons/home/kisaev/LeafletFA-utils/examples/LeafletFA_ATSE_mapper_output_20250401_150841/annotation.db"

# File containing paths to junction files (one per line)
# These junctions come from RNA-seq data and represent splicing events
junction_files_path = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/Leaflet/ATSEmap/output/junction_files.txt"

# =================================================================
# STEP 2: Preview the junction files to understand the data
# =================================================================
# Read and display the first few lines of the junction files list
# This helps to verify the input format
print("Preview of junction files list:")
with open(junction_files_path, "r") as f:
    lines = f.readlines()
    for line in lines[:5]:
        print(line.strip())

# =================================================================
# STEP 3: Configure the ATSEmapper analysis parameters
# =================================================================
class Args:
    def __init__(self):
        # Input/output parameters
        self.input = junction_files_path    # Path to junction files
        self.output = None           # If None, auto-generated output directory
        self.annotation = gtf_file          # Gene annotation file
        self.genome = fasta_file            # Genome sequence file
        self.db_path = db_path                   # Existing database for annotation (can always be left blank in which case it will be created)
        
        # Junction filtering parameters
        self.min_intron = 50                # Minimum intron length (bp)
        self.max_intron = 500000            # Maximum intron length (bp)
        self.min_reads = 200                # Min total read count for a junction
        self.min_cells = 5                  # Min number of cells with the junction
        
        # Processing parameters
        self.batch_size = 32                # Num of files to process at once
        self.num_workers = 10               # Parallel processing threads
        
        # Data type and annotation parameters
        self.sequencing_type = "single_cell"  # "single_cell" or "bulk"
        self.annotation_status = "unanno_also" # Include unannotated junctions
        self.only_canonical = True          # Keep only GT-AG, GC-AG, AT-AC sites
        
        # Alternative splicing event parameters
        self.min_splice_site_usage = 0.01   # Min proportion for a splice site
        
        # Testing and debugging parameters
        self.sample_size = 50               # Randomly sample N junction files
        self.tolerance = 100                # Allowance for matching splice sites (bp)
        self.verbose = True                 # Detailed logging output
        self.log_file = None                # Custom log file path (None = auto)

# =================================================================
# STEP 4: Run the ATSEmapper analysis
# =================================================================
print("Starting ATSEmapper analysis...")

# Parse arguements 
args = Args()  # Use your custom Args class

# run_atsemapper performs the following steps:
# 1. Loads and processes junction files
# 2. Filters junctions based on quality criteria
# 3. Checks splice site sequences against genome
# 4. Annotates junctions using the GTF file
# 5. Builds splice graphs for each gene
# 6. Identifies alternative splicing events
# 7. Saves results to the output directory
atse_file = run_atsemapper(args)
print(f"Analysis complete. Results saved to: {atse_file}")

# =================================================================
# STEP 5: Load and examine the results
# =================================================================
print("Loading ATSEmapper results...")
# The output is a tab-separated file with information about
# alternative splicing events found in the data
atse_df = pd.read_csv(atse_file, sep="\t")

# Display the first few rows of results
print("\nFirst 5 alternative splicing events found:")
print(atse_df.head())