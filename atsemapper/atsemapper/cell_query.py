#!/usr/bin/env python3
import os
import logging
import pandas as pd
from typing import Optional, Dict, Any, List


class CellJunctionQuery:
    """
    Query per-cell raw junction counts from regtools output files.

    Usage:
        query = CellJunctionQuery()
        query.build_index_from_manifest("junction_files.txt")
        result = query.query_junction("H3-D041894-3_9_M-1-1", "chr1_5135937_5143749_+")
    """

    def __init__(self):
        self.cell_index: Dict[str, str] = {}  # cell_id -> file_path

    # ------------------------------------------------------------------
    # Index building
    # ------------------------------------------------------------------

    def build_index_from_manifest(self, manifest_path: str) -> int:
        """
        Build cell_id -> file_path index from a manifest file (one path per line).
        The cell ID is inferred from the immediate parent directory of each file.

        Returns the number of cells indexed.
        """
        count = 0
        with open(manifest_path) as fh:
            for line in fh:
                file_path = line.strip()
                if file_path:
                    cell_id = os.path.basename(os.path.dirname(file_path))
                    self.cell_index[cell_id] = file_path
                    count += 1
        logging.info(f"Indexed {count} cells from {manifest_path}")
        return count

    def build_index_from_directory(
        self,
        junctions_dir: str,
        filename: str = "junctions_with_barcodes.bed",
    ) -> int:
        """
        Build index by walking a directory tree and registering every file
        named *filename*. The parent directory name becomes the cell ID.

        Returns the number of cells indexed.
        """
        count = 0
        for root, _dirs, files in os.walk(junctions_dir):
            if filename in files:
                cell_id = os.path.basename(root)
                self.cell_index[cell_id] = os.path.join(root, filename)
                count += 1
        logging.info(f"Indexed {count} cells from {junctions_dir}")
        return count

    # ------------------------------------------------------------------
    # Lookup helpers
    # ------------------------------------------------------------------

    def get_cell_path(self, cell_id: str) -> Optional[str]:
        """Return the regtools output file path for *cell_id*, or None."""
        return self.cell_index.get(cell_id)

    def list_cells(self) -> List[str]:
        """Return all indexed cell IDs."""
        return list(self.cell_index.keys())

    # ------------------------------------------------------------------
    # Junction queries
    # ------------------------------------------------------------------

    def query_junction(
        self, cell_id: str, junction_id: str
    ) -> Optional[Dict[str, Any]]:
        """
        Retrieve the raw read count for *junction_id* in *cell_id*.

        junction_id must be in the same format used by JunctionReader:
            <chrom>_<adjusted_start>_<adjusted_end>_<strand>
        e.g. "chr1_5135937_5143749_+"

        Returns a dict with keys:
            junction_id, raw_count, cell_id, file_path, [cell_readcounts]
        or None if the cell or junction is not found.
        """
        file_path = self.get_cell_path(cell_id)
        if file_path is None:
            logging.warning(f"Cell '{cell_id}' not in index")
            return None
        if not os.path.exists(file_path):
            logging.error(f"File missing: {file_path}")
            return None

        df = self._load_junction_file(file_path)
        if df is None:
            return None

        match = df[df["junction_id"] == junction_id]
        if match.empty:
            return None

        row = match.iloc[0]
        result: Dict[str, Any] = {
            "junction_id": junction_id,
            "raw_count": int(row["score"]),
            "cell_id": cell_id,
            "file_path": file_path,
        }
        if "cell_readcounts" in df.columns:
            result["cell_readcounts"] = row["cell_readcounts"]
        return result

    def query_junction_batch(
        self, cell_ids: List[str], junction_id: str
    ) -> pd.DataFrame:
        """
        Query *junction_id* across multiple cells.

        Returns a DataFrame with columns:
            cell_id, junction_id, raw_count, file_path
        Cells where the junction is absent get raw_count = 0.
        """
        rows = []
        for cell_id in cell_ids:
            r = self.query_junction(cell_id, junction_id)
            if r:
                rows.append(
                    {
                        "cell_id": cell_id,
                        "junction_id": junction_id,
                        "raw_count": r["raw_count"],
                        "file_path": r["file_path"],
                    }
                )
            else:
                rows.append(
                    {
                        "cell_id": cell_id,
                        "junction_id": junction_id,
                        "raw_count": 0,
                        "file_path": self.get_cell_path(cell_id),
                    }
                )
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _load_junction_file(file_path: str) -> Optional[pd.DataFrame]:
        """
        Read a regtools BED file and return a DataFrame with a 'junction_id'
        column using the same coordinate adjustment logic as JunctionReader.
        """
        dtypes = {
            0: str, 1: "int32", 2: "int32", 3: str, 4: "int32", 5: str,
            6: "int32", 7: "int32", 8: str, 9: "int32", 10: str, 11: str,
        }
        try:
            df = pd.read_csv(file_path, sep="\t", header=None, dtype=dtypes)
        except Exception as e:
            logging.error(f"Could not read {file_path}: {e}")
            return None

        col_names = [
            "chrom", "chromStart", "chromEnd", "name", "score", "strand",
            "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
        ]
        if len(df.columns) >= 13:
            col_names.append("num_cells_wjunc")
        if len(df.columns) >= 14:
            col_names.append("cell_readcounts")
        df.columns = col_names[: len(df.columns)]

        # Mirror the coordinate adjustment in JunctionReader.parse_file
        extracted = df["blockSizes"].str.extract(r"(\d+),(\d+)").astype(int)
        df["chromStart"] = df["chromStart"] + extracted[0]
        df["chromEnd"] = df["chromEnd"] - extracted[1]

        df["junction_id"] = (
            df["chrom"] + "_"
            + df["chromStart"].astype(str) + "_"
            + df["chromEnd"].astype(str) + "_"
            + df["strand"]
        )
        return df
