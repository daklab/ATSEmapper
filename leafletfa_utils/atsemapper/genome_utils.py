#!/usr/bin/env python3
import os
import time
import sqlite3
from typing import Dict
from pyfaidx import Fasta 
import gffutils 
import pandas as pd 
from tqdm import tqdm 

class GenomeDB:
    def __init__(self, db_name: str, gtf_file: str = None, fasta_file: str = None, max_retries: int = 3):
        self.db_name = db_name
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        self.genome = None
        self.db = None
        
        if self.fasta_file:
            print(f"Loading genome from {self.fasta_file}")
            self.genome = Fasta(self.fasta_file)
            
        if self.gtf_file and self.db_name:
            if not os.path.exists(self.db_name):
                print(f"Creating database {self.db_name} from {self.gtf_file}")
                self.db = self.create_db()
            else:
                print(f"Loading existing database {self.db_name}")
                for attempt in range(max_retries):
                    try:
                        self.db = gffutils.FeatureDB(self.db_name)
                        break
                    except sqlite3.OperationalError as e:
                        if attempt == max_retries - 1:
                            raise e
                        print(f"Database locked, retrying... (attempt {attempt + 1}/{max_retries})")
                        time.sleep(1)  # Wait before retry
        
    def create_db(self) -> gffutils.FeatureDB:
        return gffutils.create_db(
            self.gtf_file,
            dbfn=self.db_name,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
        )

    def get_db(self) -> gffutils.FeatureDB:
        return self.db

class JunctionAnalyzer:
    def __init__(self, fasta_file: str, db: gffutils.FeatureDB, window_size: int = 2, tolerance: int = 1):
        self.genome = Fasta(fasta_file)
        self.db = db
        self.window_size = window_size
        self.tolerance = tolerance
        self.canonical_motifs = {
            "GT-AG": ("GT", "AG"),
            "AT-AC": ("AT", "AC"),
            "GC-AG": ("GC", "AG")
        }

    @staticmethod
    def reverse_complement(seq: str) -> str:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def get_splice_sequences(self, junction: Dict):
        chrom = self.check_chromosome(junction["chrom"])
        
        if junction["strand"] == "+":
            donor_seq = str(self.genome[chrom][junction["start"]:junction["start"]+self.window_size]).upper()
            acceptor_seq = str(self.genome[chrom][junction["end"]-self.window_size:junction["end"]]).upper()
        else:
            donor_seq = self.reverse_complement(str(self.genome[chrom][junction["end"]-self.window_size:junction["end"]]).upper())
            acceptor_seq = self.reverse_complement(str(self.genome[chrom][junction["start"]:junction["start"]+self.window_size]).upper())
            
        return donor_seq, acceptor_seq

    def check_chromosome(self, chrom: str) -> str:
        if chrom not in self.genome:
            return f"chr{chrom}" if not chrom.startswith("chr") else chrom.replace("chr", "")
        return chrom

    def get_motif_type(self, donor_seq: str, acceptor_seq: str) -> str:
        for motif_name, (donor, acceptor) in self.canonical_motifs.items():
            if donor_seq == donor and acceptor_seq == acceptor:
                return motif_name
        return "non-canonical"

    def check_splice_sites(self, junctions: Dict[str, Dict]) -> Dict[str, Dict]:
        """Add splice site annotations to junction dictionary"""
        for j_id, junction in junctions.items():
            donor_seq, acceptor_seq = self.get_splice_sequences(junction)
            junction["splice_motif"] = self.get_motif_type(donor_seq, acceptor_seq)
            junction["donor_seq"] = donor_seq
            junction["acceptor_seq"] = acceptor_seq
        return junctions

    def filter_canonical(self, junctions: Dict[str, Dict], only_canonical: bool = True) -> Dict[str, Dict]:
       """Filter junctions based on splice site motifs"""
       if not only_canonical:
           return junctions

       initial_count = len(junctions)
       motif_counts = {}
    
       for j_data in junctions.values():
           motif = j_data["splice_motif"]
           motif_counts[motif] = motif_counts.get(motif, 0) + 1
    
       canonical_junctions = {j_id: j_data for j_id, j_data in junctions.items() 
                            if j_data["splice_motif"] != "non-canonical"}
    
       filtered_count = len(canonical_junctions)
    
       print(f"Initial junctions before splice site filtering: {initial_count}")
       print("\nSplice site motif breakdown:")
       for motif, count in motif_counts.items():
           print(f"{motif}: {count} ({(count/initial_count)*100:.2f}%)")
       print(f"\nCanonical junctions kept: {filtered_count}")
       print(f"Non-canonical removed: {initial_count - filtered_count}")
       print(f"Percentage canonical: {(filtered_count/initial_count)*100:.2f}%")
    
       return canonical_junctions
    
    def check_junction_annotation(self, junctions: Dict[str, Dict]) -> Dict[str, Dict]:
        """Add annotation status for both ends of each junction with improved transcript mapping"""

        # Group junctions by chromosome to minimize database queries
        chrom_groups = {}
        for j_id, j in junctions.items():
            chrom_groups.setdefault(j["chrom"], []).append((j_id, j))

        for chrom, junc_group in tqdm(chrom_groups.items(), desc="Processing chromosomes"):
            # Get all relevant transcripts for this chromosome group
            min_pos = min(j["start"] for _, j in junc_group)
            max_pos = max(j["end"] for _, j in junc_group)

            # Single query for all transcripts in range
            transcripts = list(self.db.region(
                region=(chrom, min_pos - 1000, max_pos + 1000),
                featuretype="transcript"
            ))

            # Cache exon data and transcript data 
            exon_cache = {t.id: list(self.db.children(t, featuretype="exon", order_by="start"))
                         for t in transcripts}

            # Get gene info and transcripts
            gene_info = {}
            for t in transcripts:
                gene = list(self.db.parents(t, featuretype="gene"))[0]
                gene_info[t.id] = {
                    'gene_id': gene.id,
                    'gene_name': gene.attributes.get('gene_name', [None])[0],
                    'gene_type': gene.attributes.get('gene_biotype', 
                                  gene.attributes.get('gene_type', [None]))[0]
                }

            # Get trancsript types 
            transcript_types = {}
            for t in transcripts:
                transcript_types[t.id] = t.attributes.get('transcript_type', [None])[0]

            # Process each junction
            for j_id, junction in tqdm(junc_group, desc="Processing junctions"):
                start = junction["start"]
                end = junction["end"]
                strand = junction["strand"]

                # Initialize overall annotation labels
                label_5_prime, label_3_prime = "unannotated on 5'", "unannotated on 3'"
                position_off_5_prime, position_off_3_prime = None, None

                # Initialize transcript categorization
                both_ends_transcripts = set()
                only_5_prime_transcripts = set()
                only_3_prime_transcripts = set()

                # Track overall for all transcripts
                transcript_types_junc = set()
                genes_found = set()
                gene_types_found = set()  # New: Track gene types
                transcripts_junc_5 = set()
                transcripts_junc_3 = set()

                for transcript in transcripts:
                    exons = exon_cache[transcript.id]
                    found_5_prime = False
                    found_3_prime = False

                    if strand == "+":
                        # Check 5' end (start position)
                        for exon in exons:
                            if abs(exon.end - start) <= self.tolerance:
                                label_5_prime = "annotated on 5'"
                                position_off_5_prime = exon.end - start
                                transcripts_junc_5.add(transcript.id)
                                found_5_prime = True
                                break 
                        # Check 3' end (end position)
                        for exon in exons:
                            if abs(exon.start - end) <= self.tolerance:
                                label_3_prime = "annotated on 3'"
                                position_off_3_prime = exon.start - end
                                found_3_prime = True
                                transcripts_junc_3.add(transcript.id)
                                break
                    else: # strand == "-"
                        # Check 5' end (end position for negative strand)
                        for exon in exons:
                            if abs(exon.start - end) <= self.tolerance:
                                label_5_prime = "annotated on 5'"
                                position_off_5_prime = exon.start - end
                                transcripts_junc_5.add(transcript.id)
                                found_5_prime = True
                                break

                        # Check 3' end (start position for negative strand)
                        for exon in exons:
                            if abs(exon.end - start) <= self.tolerance:
                                label_3_prime = "annotated on 3'"
                                position_off_3_prime = exon.end - start
                                found_3_prime = True
                                transcripts_junc_3.add(transcript.id)
                                break

                    # Add gene and transcript info if any end matches
                    if found_5_prime or found_3_prime:
                        if transcript.id in transcript_types:
                            transcript_types_junc.add(transcript_types[transcript.id])
                        if transcript.id in gene_info:
                            gene_data = gene_info[transcript.id]
                            genes_found.add((gene_data['gene_id'], gene_data['gene_name']))
                            if gene_data['gene_type']:
                                gene_types_found.add(gene_data['gene_type'])

                    # Categorize transcripts based on which junction ends they match
                    if found_5_prime and found_3_prime:
                        both_ends_transcripts.add(transcript.id)
                    elif found_5_prime:
                        only_5_prime_transcripts.add(transcript.id)
                    elif found_3_prime:
                        only_3_prime_transcripts.add(transcript.id)

                # Get all transcripts that match at least one end
                all_matching_transcripts = both_ends_transcripts.union(only_5_prime_transcripts, only_3_prime_transcripts)

                junction.update({
                    "label_5_prime": label_5_prime,
                    "label_3_prime": label_3_prime,
                    "position_off_5_prime": position_off_5_prime,
                    "position_off_3_prime": position_off_3_prime,
                    "transcripts": list(all_matching_transcripts),
                    "both_ends_transcripts": list(both_ends_transcripts),
                    "only_5_prime_transcripts": list(only_5_prime_transcripts),
                    "only_3_prime_transcripts": list(only_3_prime_transcripts),
                    "transcript_types": list(transcript_types_junc),
                    "gene_ids": [g[0] for g in genes_found],
                    "gene_names": [g[1] for g in genes_found],
                    "transcripts_junc_5": list(transcripts_junc_5),
                    "transcripts_junc_3": list(transcripts_junc_3),
                    "gene_types": list(gene_types_found)  # Added gene types to output
                })

        return junctions
    
    def filter_annotated(self, junctions: Dict[str, Dict], annotation_status_include: str = 'both') -> Dict[str, Dict]:
        """
        Filter junctions based on annotation status.
    
        Args:
            junctions: Dictionary of junctions
            annotation_status_include: Filtering option for junctions
                - 'both': Keep only junctions where both ends are annotated
                - 'either': Keep junctions where at least one end is annotated
                - 'unanno_also': Keep all junctions regardless of annotation

        Returns:
            Filtered junction dictionary
        """
        if annotation_status_include not in ['both', 'either', 'unanno_also']:
            raise ValueError("annotation_status_include must be one of: 'both', 'either', 'unanno_also'")
        
        initial_count = len(junctions)

        annotation_stats = {
        "both_ends": 0,
        "five_prime_only": 0,
        "three_prime_only": 0,
        "unannotated": 0,
        "multi_gene": 0
        }

        filtered_junctions = {}
        multi_gene_junctions = {}

        # Process each junction
        for j_id, j_data in junctions.items():
            
            # Initialize flags
            five_prime = False
            three_prime = False

            # Check position off for 5' and 3' ends if not None
            if j_data["position_off_5_prime"] is not None:
                five_prime = (0 <= j_data["position_off_5_prime"] <= 1)

            if j_data["position_off_3_prime"] is not None:
                three_prime = (0 <= j_data["position_off_3_prime"] <= 1)
    
            # Handle multi-gene junctions
            if len(j_data["gene_ids"]) > 1:
                annotation_stats["multi_gene"] += 1
                multi_gene_junctions[j_id] = j_data
                continue

            # Categorize junction
            if five_prime and three_prime:
                annotation_stats["both_ends"] += 1
                category = "both"
            elif five_prime:
                annotation_stats["five_prime_only"] += 1
                category = "five_prime"
            elif three_prime:
                annotation_stats["three_prime_only"] += 1
                category = "three_prime"
            else:
                annotation_stats["unannotated"] += 1
                category = "unannotated"

            # Apply filtering based on annotation_status_include
            keep_junction = False
            if annotation_status_include == 'both':
                keep_junction = (category == "both")
            elif annotation_status_include == 'either':
                keep_junction = (category in ["both", "five_prime", "three_prime"])
            elif annotation_status_include == 'unanno_also':
                keep_junction = True

            if keep_junction:
                # Add junction to filtered dictionary and include its category
                j_data["annotation_status"] = category
                filtered_junctions[j_id] = j_data

        filtered_count = len(filtered_junctions)
        
        # Save multi-gene junctions to file
        if multi_gene_junctions:
            multi_gene_file = f'multi_gene_junctions_{time.strftime("%Y%m%d-%H%M%S")}.csv'
            pd.DataFrame.from_dict(multi_gene_junctions, orient='index').to_csv(multi_gene_file)
            print(f"\nMulti-gene junctions saved to {multi_gene_file}")

        # Print statistics
        print(f"\nInitial junctions: {initial_count}")
        print("\nAnnotation breakdown:")
        for status, count in annotation_stats.items():
            print(f"{status}: {count} ({(count/initial_count)*100:.2f}%)")

        print(f"\nFiltering summary:")
        print(f"Annotation status filter: {annotation_status_include}")
        print(f"Junctions kept: {filtered_count}")
        print(f"Junctions removed: {initial_count - filtered_count}")
        print(f"Percentage kept: {(filtered_count/initial_count)*100:.2f}%")

        return filtered_junctions