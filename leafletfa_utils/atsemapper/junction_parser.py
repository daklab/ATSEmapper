#!/usr/bin/env python3
import pandas as pd
import logging
from typing import Dict, List
from tqdm import tqdm # type: ignore
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

class JunctionReader:
    def __init__(self, 
                min_intron: int = 50, 
                max_intron: int = 500000, 
                sequencing_type: str = "single_cell", 
                batch_size: int = 10,
                min_cells: int = 2,
                min_reads: int = 100, 
                num_workers: int = 4):
        
        self.min_intron = min_intron
        self.max_intron = max_intron
        self.sequencing_type = sequencing_type
        self.dtypes = {0: str, 1: 'int32', 2: 'int32', 3: str, 4: 'int32', 5: str,
                      6: 'int32', 7: 'int32', 8: str, 9: 'int32', 10: str, 11: str}
        self.batch_size = batch_size
        self.min_cells = min_cells
        self.min_reads = min_reads
        self.num_workers = num_workers
        self.combined_junctions = {}

    def parse_file(self, file_path: str) -> Dict[str, Dict]:
        """
        Read and process a junction file, returning a dictionary of junctions with cell information.

        Args:
            file_path (str): Path to the junction file

        Returns:
            Dict[str, Dict]: Dictionary containing junction information with cell-specific details
        """
        try:
            juncs = pd.read_csv(file_path, sep="\t", header=None, dtype=self.dtypes)

            col_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]
            if self.sequencing_type == "single_cell":
                col_names += ["num_cells_wjunc", "cell_readcounts"]
            juncs.columns = col_names[:len(juncs.columns)]

            juncs[['block_add_start', 'block_subtract_end']] = (
                juncs["blockSizes"].str.extract(r'(\d+),(\d+)').astype(int)
            )
            juncs["chromStart"] += juncs['block_add_start']
            juncs["chromEnd"] -= juncs['block_subtract_end']
            juncs["intron_length"] = juncs["chromEnd"] - juncs["chromStart"]

            juncs = juncs[
                (juncs["intron_length"] >= self.min_intron) & 
                (juncs["intron_length"] <= self.max_intron)
            ]
            standard_chromosomes_pattern = r'^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)$'
            juncs = juncs[juncs['chrom'].str.match(standard_chromosomes_pattern)]

            juncs['junction_id'] = (
                juncs['chrom'] + '_' + 
                juncs['chromStart'].astype(str) + '_' +
                juncs['chromEnd'].astype(str) + '_' + 
                juncs['strand']
            )

            junc_dict = {}
            for _, row in juncs.iterrows():
                junction_id = row['junction_id']
                if junction_id not in junc_dict:
                    junc_dict[junction_id] = {
                        'cells': 1,
                        'total_score': row['score'],
                        'chrom': row['chrom'],
                        'start': row['chromStart'],
                        'end': row['chromEnd'],
                        'strand': row['strand']
                    }
                else:
                    junc_dict[junction_id]['cells'] += 1
                    junc_dict[junction_id]['total_score'] += row['score']

            return junc_dict
        
        except Exception as e:
            logging.error(f"Could not read in {file_path}: {e}")
            return {}
        
    def process_files(self, file_list: List[str]) -> Dict[str, Dict]:
        """
        Process multiple junction files in parallel batches.

        Args:
            file_list: List of file paths

        Returns:
            Combined junction dictionary
        """
        def process_batch(batch_files):
            batch_junctions = {}
            for file_path in batch_files:
                try:
                    file_junctions = self.parse_file(file_path)
                    for j_id, j_data in file_junctions.items():
                        if j_id not in batch_junctions:
                            batch_junctions[j_id] = j_data
                        else:
                            batch_junctions[j_id]['cells'] += j_data['cells']
                            batch_junctions[j_id]['total_score'] += j_data['total_score']
                except Exception as e:
                    logging.error(f"Error processing file {file_path}: {str(e)}")
            return batch_junctions

        # Create batches
        batches = [file_list[i:i + self.batch_size] for i in range(0, len(file_list), self.batch_size)]
        total_batches = len(batches)

        print(f"\nProcessing {len(file_list)} files in {total_batches} batches")
        print(f"Using {self.num_workers} workers, batch size: {self.batch_size}\n")

        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            futures = [executor.submit(process_batch, batch) for batch in batches]

            with tqdm(total=total_batches, desc="Processing batches") as pbar:
                for future in as_completed(futures):
                    try:
                        batch_result = future.result()

                        # Merge results
                        for j_id, j_data in batch_result.items():
                            if j_id not in self.combined_junctions:
                                self.combined_junctions[j_id] = j_data
                            else:
                                self.combined_junctions[j_id]['cells'] += j_data['cells']
                                self.combined_junctions[j_id]['total_score'] += j_data['total_score']

                        pbar.update(1)

                    except Exception as e:
                        logging.error(f"Batch processing error: {str(e)}")
                        pbar.update(1)

        return self.combined_junctions

    def SJ_QC(self, junctions: Dict[str, Dict]) -> Dict[str, Dict]:
        """
        Filter junctions based on minimum cell and read requirements.

        Args:
            junctions: Dictionary of junctions
            min_cells: Minimum number of cells required
            min_reads: Minimum number of reads required

        Returns:
            Filtered junction dictionary
        """
        initial_count = len(junctions)

        # Filter based on cells and reads
        filtered_junctions = {
            j_id: j_data for j_id, j_data in junctions.items()
            if j_data['cells'] >= self.min_cells and j_data['total_score'] >= self.min_reads
        }

        # Log filtering results
        filtered_count = len(filtered_junctions)
        # To log also add which filters were used
        print(f"Minimum cells per junction to pass set to: {self.min_cells}")
        print(f"Minimum total reads per junction to pass set to: {self.min_reads}")
        print(f"Initial junctions before basic QC: {initial_count}")
        print(f"Filtered junctions: {filtered_count}")
        print(f"Removed junctions: {initial_count - filtered_count}")
        print(f"Percentage kept: {(filtered_count/initial_count)*100:.2f}%")

        return filtered_junctions

    def clear_combined_junctions(self):
        """Reset the combined junctions dictionary"""
        self.combined_junctions = {}