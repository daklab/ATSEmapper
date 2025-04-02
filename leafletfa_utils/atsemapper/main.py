#!/usr/bin/env python3
import os
import argparse
import logging
from pathlib import Path
import time
import datetime
import os
import datetime

from .junction_parser import JunctionReader
from .event_detection import ATSEAnalyzer
from .genome_utils import GenomeDB, JunctionAnalyzer

def parse_arguments():
    parser = argparse.ArgumentParser(description="ATSEmapper: Map alternative splicing events from splice junction files")
    
    # Required parameters
    parser.add_argument("--input", required=True, 
                        help="Directory containing junction files or path to a file listing junction files")
    parser.add_argument("--output", required=False, default=None,
                        help="Directory where output files will be saved (default: LeafletFA_ATSE_mapper_output_DATE)")
    parser.add_argument("--annotation", required=True, 
                        help="GTF/GFF3 file with genome annotation")
    parser.add_argument("--genome", required=True, 
                        help="FASTA file with genome sequence")
    parser.add_argument("--db_path", required=False, default=None,
                        help="Path to SQLite database for genome annotation (default: {output}/annotation.db)")
    
    # Junction filtering options
    parser.add_argument("--min-intron", type=int, default=50, 
                        help="Minimum intron length (default: 50)")
    parser.add_argument("--max-intron", type=int, default=500000, 
                        help="Maximum intron length (default: 500000)")
    parser.add_argument("--min-cells", type=int, default=2, 
                        help="Minimum cells with a junction (default: 2)")
    parser.add_argument("--min-reads", type=int, default=100, 
                        help="Minimum total reads for a junction (default: 100)")
    
    # Processing options
    parser.add_argument("--batch-size", type=int, default=10, 
                        help="Number of files to process in a batch (default: 10)")
    parser.add_argument("--num-workers", type=int, default=4, 
                        help="Number of worker threads (default: 4)")
    parser.add_argument("--sequencing-type", choices=["single_cell", "bulk"], default="single_cell", 
                        help="Sequencing data type (default: single_cell)")
    
    # Annotation options
    parser.add_argument("--annotation-status", choices=["both", "either", "unanno_also"], default="both",
                       help="Junction annotation filter (default: both)")
    parser.add_argument("--only-canonical", action="store_true", 
                        help="Keep only canonical splice sites")
    parser.add_argument("--tolerance", type=int, default=1,
                        help="Tolerance for matching splice sites to annotated exons (default: 1)")
    
    # ATSE options
    parser.add_argument("--min-splice-site-usage", type=float, default=0.01, 
                       help="Minimum proportion of reads a junction must have at a splice site (default: 0.01)")
    
    # Other options
    parser.add_argument("--sample-size", type=int, 
                        help="Sample N junction files (for testing)")
    parser.add_argument("--log-file", 
                        help="Path to log file (default: {output}/atsemapper_{timestamp}.log)")
    parser.add_argument("--verbose", action="store_true", 
                        help="Print verbose output")
    
    args = parser.parse_args()
    
    # Handle output directory creation if not specified
    if not args.output:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        args.output = f"LeafletFA_ATSE_mapper_output_{timestamp}"
        print(f"No output directory specified. Using: {args.output}")

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    return args

def run_atsemapper(args=None):
    """Main entry point for ATSEmapper"""
    
    # Parse arguments if not provided
    if args is None:
        args = parse_arguments()
    else: 
        # Handle the case when args is an object but output might be None
        if not getattr(args, 'output', None):
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            args.output = f"LeafletFA_ATSE_mapper_output_{timestamp}"
            print(f"No output directory specified. Using: {args.output}")
            
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Set up logging
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = args.log_file if args.log_file else os.path.join(args.output, f"atsemapper_{timestamp}.log")
    
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logging.info(f"Starting ATSEmapper with parameters: {vars(args)}")
    
    # Initialize genome database
    # Check if db_path is provided, if not use default 
    if args.db_path:
        db_path = args.db_path
        # Ensure file actually exists
        if not os.path.exists(db_path):
            logging.error(f"Database path {db_path} does not exist.")
            raise FileNotFoundError(f"Database path {db_path} does not exist.")
    else:
        # Create a default path for the database in the output directory
        db_path = os.path.join(args.output, "annotation.db")
        # Check if the database already exists
        if os.path.exists(db_path):
            logging.warning(f"Database already exists at {db_path}. Using existing database.")
        else:
            # Create the database if it doesn't exist
            logging.info(f"Creating new database at {db_path}")

    genome_db = GenomeDB(db_path, args.annotation, args.genome)
    
    # Get list of junction files
    junction_files = []
    if os.path.isdir(args.input):
        logging.info(f"Looking for junction files in directory: {args.input}")
        junction_files = [os.path.join(args.input, f) for f in os.listdir(args.input) 
                         if f.endswith('.bed') or f.endswith('.junc')]
    elif os.path.isfile(args.input) and args.input.endswith('.txt'):
        logging.info(f"Reading junction file list from: {args.input}")
        with open(args.input, 'r') as f:
            junction_files = [line.strip() for line in f if line.strip()]
    else:
        junction_files = [args.input]
        
    logging.info(f"Found {len(junction_files)} junction files")
    
    # Sample files if requested
    if args.sample_size and args.sample_size < len(junction_files):
        import random
        random.seed(42)  # For reproducibility
        junction_files = random.sample(junction_files, args.sample_size)
        logging.info(f"Sampled {args.sample_size} junction files for analysis")
    
    # Initialize junction reader
    logging.info("Initializing junction reader")
    reader = JunctionReader(
        min_intron=args.min_intron,
        max_intron=args.max_intron,
        sequencing_type=args.sequencing_type,
        batch_size=args.batch_size,
        min_cells=args.min_cells,
        min_reads=args.min_reads,
        num_workers=args.num_workers
    )
    
    # Process junction files
    logging.info("Processing junction files")
    start_time = time.time()
    junctions = reader.process_files(junction_files)
    elapsed = time.time() - start_time
    logging.info(f"Processed {len(junction_files)} files in {elapsed:.2f} seconds")
    logging.info(f"Found {len(junctions)} unique junctions")
    
    # Filter junctions based on QC parameters
    logging.info("Performing basic junction QC")
    filtered_junctions = reader.SJ_QC(junctions)
    logging.info(f"After QC: {len(filtered_junctions)} junctions remain")
    
    # Initialize junction analyzer
    logging.info("Initializing junction analyzer")
    analyzer = JunctionAnalyzer(
        fasta_file=args.genome, 
        db=genome_db.get_db(), 
        tolerance=args.tolerance
    )
    
    # Check splice sites
    logging.info("Checking splice site sequences")
    junctions_with_motifs = analyzer.check_splice_sites(filtered_junctions)
    
    # Filter for canonical splice sites if requested
    if args.only_canonical:
        logging.info("Filtering for canonical splice sites")
        canonical_junctions = analyzer.filter_canonical(junctions_with_motifs)
    else:
        canonical_junctions = junctions_with_motifs
        logging.info("Keeping all splice sites (canonical and non-canonical)")
    
    logging.info(f"After splice site filtering: {len(canonical_junctions)} junctions remain")
    
    # Add junction annotation
    logging.info("Adding junction annotation")
    start_time = time.time()
    annotated_junctions = analyzer.check_junction_annotation(canonical_junctions)
    elapsed = time.time() - start_time
    logging.info(f"Annotation completed in {elapsed:.2f} seconds")
    
    # Filter based on annotation status
    logging.info(f"Filtering junctions based on annotation status: {args.annotation_status}")
    annotation_filtered = analyzer.filter_annotated(annotated_junctions, args.annotation_status)
    logging.info(f"After annotation filtering: {len(annotation_filtered)} junctions remain")
    
    # Build splice graphs and find ATSEs
    logging.info("Building splice graphs")
    atse_analyzer = ATSEAnalyzer()
    gene_graphs, graph_stats = atse_analyzer.build_splice_graph(annotation_filtered)
    
    logging.info("Finding ATSE groups")
    atse_groups, junction_counts = atse_analyzer.find_atse_groups(
        gene_graphs, min_splice_site_usage=args.min_splice_site_usage
    )
    logging.info(f"Found {len(atse_groups)} ATSE groups")
    
    # Classify events
    logging.info("Classifying ATSE events")
    classified_events, event_counts = atse_analyzer.classify_events(gene_graphs, atse_groups)
    for event_type, count in event_counts.items():
        logging.info(f"  {event_type}: {count}")
    
    # Save results
    today = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    output_file = os.path.join(args.output, f"atse_events_{today}.tsv.gz")
    logging.info(f"Saving ATSE events to {output_file}")
    atse_analyzer.save_atse_file(classified_events, annotation_filtered, output_file)    
    logging.info("ATSEmapper completed successfully")
    return output_file

def main():
    """Command-line entry point"""
    output_file = run_atsemapper()
    return output_file

if __name__ == "__main__":
    main()