import os
import tempfile
import pytest
from atsemapper.atsemapper.junction_parser import JunctionReader


# Minimal valid regtools BED line: blockSizes "5,3," means
# actual start = chromStart + 5, actual end = chromEnd - 3
def _make_bed_line(chrom="chr1", start=1000, end=1200, score=50, strand="+",
                   block_sizes="5,3,"):
    return f"{chrom}\t{start}\t{end}\tJUNC\t{score}\t{strand}\t{start}\t{end}\t255,0,0\t2\t{block_sizes}\t0,197\n"


def _write_bed(lines, suffix=".bed"):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    f.writelines(lines)
    f.flush()
    f.close()
    return f.name


# ---------------------------------------------------------------------------
# parse_file
# ---------------------------------------------------------------------------

class TestParseFile:
    def setup_method(self):
        self.reader = JunctionReader(min_intron=50, max_intron=500000,
                                     sequencing_type="bulk", min_cells=1, min_reads=1)

    def test_returns_dict_with_expected_keys(self):
        path = _write_bed([_make_bed_line()])
        try:
            result = self.reader.parse_file(path)
            assert len(result) == 1
            jdata = next(iter(result.values()))
            for key in ("cells", "total_score", "chrom", "start", "end", "strand"):
                assert key in jdata
        finally:
            os.unlink(path)

    def test_junction_id_format(self):
        path = _write_bed([_make_bed_line(chrom="chr1", start=1000, end=1200, strand="+")])
        try:
            result = self.reader.parse_file(path)
            jid = next(iter(result))
            parts = jid.split("_")
            assert parts[0] == "chr1"
            assert parts[-1] == "+"
            assert parts[1].isdigit() and parts[2].isdigit()
        finally:
            os.unlink(path)

    def test_block_size_adjustment(self):
        # blockSizes "5,3," → start += 5, end -= 3
        path = _write_bed([_make_bed_line(start=1000, end=1200, block_sizes="5,3,")])
        try:
            result = self.reader.parse_file(path)
            jdata = next(iter(result.values()))
            assert jdata["start"] == 1005
            assert jdata["end"] == 1197
        finally:
            os.unlink(path)

    def test_short_intron_filtered_out(self):
        # adjusted length = (1030 - 3) - (1000 + 5) = 22 < min_intron 50
        path = _write_bed([_make_bed_line(start=1000, end=1030, block_sizes="5,3,")])
        try:
            result = self.reader.parse_file(path)
            assert len(result) == 0
        finally:
            os.unlink(path)

    def test_long_intron_filtered_out(self):
        reader = JunctionReader(max_intron=100, sequencing_type="bulk")
        path = _write_bed([_make_bed_line(start=1000, end=2000, block_sizes="5,3,")])
        try:
            result = reader.parse_file(path)
            assert len(result) == 0
        finally:
            os.unlink(path)

    def test_non_standard_chromosome_filtered(self):
        path = _write_bed([_make_bed_line(chrom="chrUn_gl000220")])
        try:
            result = self.reader.parse_file(path)
            assert len(result) == 0
        finally:
            os.unlink(path)

    def test_standard_chromosomes_kept(self):
        lines = [
            _make_bed_line(chrom="chr1"),
            _make_bed_line(chrom="chrX"),
            _make_bed_line(chrom="chrY"),
            _make_bed_line(chrom="chrMT"),
        ]
        path = _write_bed(lines)
        try:
            result = self.reader.parse_file(path)
            assert len(result) == 4
        finally:
            os.unlink(path)

    def test_duplicate_junctions_merged(self):
        # Same junction twice in one file → cells=2, score summed
        line = _make_bed_line(chrom="chr1", start=1000, end=1200, score=10)
        path = _write_bed([line, line])
        try:
            result = self.reader.parse_file(path)
            assert len(result) == 1
            jdata = next(iter(result.values()))
            assert jdata["cells"] == 2
            assert jdata["total_score"] == 20
        finally:
            os.unlink(path)

    def test_missing_file_returns_empty_dict(self):
        result = self.reader.parse_file("/nonexistent/path/file.bed")
        assert result == {}


# ---------------------------------------------------------------------------
# SJ_QC
# ---------------------------------------------------------------------------

class TestSJQC:
    def _make_junctions(self, entries):
        return {
            f"chr1_{i}_1200_+": {
                "cells": cells, "total_score": reads,
                "chrom": "chr1", "start": i, "end": 1200, "strand": "+"
            }
            for i, (cells, reads) in enumerate(entries)
        }

    def test_filters_below_min_cells(self):
        reader = JunctionReader(min_cells=3, min_reads=10)
        juncs = self._make_junctions([(1, 100), (3, 100)])
        result = reader.SJ_QC(juncs)
        assert len(result) == 1

    def test_filters_below_min_reads(self):
        reader = JunctionReader(min_cells=1, min_reads=50)
        juncs = self._make_junctions([(5, 10), (5, 100)])
        result = reader.SJ_QC(juncs)
        assert len(result) == 1

    def test_passes_junction_meeting_both_thresholds(self):
        reader = JunctionReader(min_cells=2, min_reads=20)
        juncs = self._make_junctions([(2, 20)])
        result = reader.SJ_QC(juncs)
        assert len(result) == 1

    def test_junction_below_both_thresholds_removed(self):
        reader = JunctionReader(min_cells=2, min_reads=50)
        juncs = self._make_junctions([(1, 10)])
        result = reader.SJ_QC(juncs)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# process_files
# ---------------------------------------------------------------------------

class TestProcessFiles:
    def test_two_files_same_junction_merged(self):
        line = _make_bed_line(chrom="chr2", start=5000, end=6000, score=30)
        paths = [_write_bed([line]), _write_bed([line])]
        try:
            reader = JunctionReader(sequencing_type="bulk")
            result = reader.process_files(paths)
            assert len(result) == 1
            jdata = next(iter(result.values()))
            assert jdata["cells"] == 2
            assert jdata["total_score"] == 60
        finally:
            for p in paths:
                os.unlink(p)

    def test_two_files_distinct_junctions(self):
        path1 = _write_bed([_make_bed_line(chrom="chr1", start=1000, end=1200)])
        path2 = _write_bed([_make_bed_line(chrom="chr2", start=2000, end=2400)])
        try:
            reader = JunctionReader(sequencing_type="bulk")
            result = reader.process_files([path1, path2])
            assert len(result) == 2
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_clear_resets_state(self):
        path = _write_bed([_make_bed_line()])
        try:
            reader = JunctionReader(sequencing_type="bulk")
            reader.process_files([path])
            assert len(reader.combined_junctions) > 0
            reader.clear_combined_junctions()
            assert reader.combined_junctions == {}
        finally:
            os.unlink(path)
