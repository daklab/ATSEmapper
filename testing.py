# Import necessary libraries
import os
import pandas as pd
import sys
import datetime
import gffutils  # type: ignore

# Import your modules
from leafletfa_utils.atsemapper.main import run_atsemapper
from leafletfa_utils.atsemapper.junction_parser import JunctionReader
from leafletfa_utils.atsemapper.event_detection import ATSEAnalyzer
from leafletfa_utils.atsemapper.genome_utils import GenomeDB, JunctionAnalyzer
import leafletfa_utils.atseviz as ja

# Set up paths for genome files and junction files
gtf_file = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/genome_files/gencode.vM19/genes/genes.gtf"
fasta_file = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/genome_files/gencode.vM19/fasta/genome.fa"
output_path = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/Leaflet/ATSEmap/output/"
junction_files_path = "/gpfs/commons/groups/knowles_lab/Karin/Leaflet-analysis-WD/TabulaSenis/Leaflet/ATSEmap/output/junction_files.txt"

# Create args object
class Args:
    def __init__(self):
        self.input = junction_files_path
        self.output = output_path
        self.annotation = gtf_file
        self.genome = fasta_file
        self.min_intron = 50
        self.max_intron = 500000
        self.min_reads = 200
        self.min_cells = 20
        self.batch_size = 32
        self.num_workers = 10
        self.sequencing_type = "single_cell"
        self.annotation_status = "unanno_also"
        self.only_canonical = True
        self.min_splice_site_usage = 0.01
        self.sample_size = 75
        self.tolerance = 100
        self.verbose = True
        self.log_file = None

# Run the analysis
args = Args()
atse_file = run_atsemapper(args)

# Visualization part
p = pd.read_csv(atse_file, sep="\t")
db = gffutils.FeatureDB("gencodeVM19", keep_order=True)
atse_event = p.sample(1)["event_id"].values[0]
juncs = p[p["event_id"]==atse_event]
juncs["usage_ratio"] = 0
juncs["Cluster"] = juncs["event_id"]
splice_junctions = ja.convert_junction_ids(juncs)

# Continue with visualization...