"""Tests for vfam_trees.concat — pure helpers (no MAFFT / trimAl dependency)."""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees.concat import (
    _safe_charset_name,
    concatenate_aligned_markers,
    remove_per_marker_length_outliers,
    write_partition_file_nexus,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rec(acc: str, length: int) -> SeqRecord:
    return SeqRecord(Seq("M" * length), id=acc, description="")


def _aln(acc: str, sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id=acc, description="")


# ---------------------------------------------------------------------------
# remove_per_marker_length_outliers
# ---------------------------------------------------------------------------

class TestPerMarkerLengthOutliers:
    def test_normal_lengths_kept(self):
        # Three genomes each with 2 markers, all lengths near each marker's median
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000), "MCP": _rec("YP_2", 600)},
            "NC_002": {"polB": _rec("YP_3", 1010), "MCP": _rec("YP_4", 590)},
            "NC_003": {"polB": _rec("YP_5",  990), "MCP": _rec("YP_6", 610)},
        }
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 0
        assert stats["n_short_dropped"] == 0
        assert all(len(m) == 2 for m in updated.values())

    def test_long_outlier_dropped(self):
        # NC_002 polB is 4× the median → dropped (default hi_mult = 3.0)
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000)},
            "NC_002": {"polB": _rec("YP_2", 4500)},  # huge
            "NC_003": {"polB": _rec("YP_3",  990)},
        }
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 1
        assert "polB" not in updated["NC_002"]
        assert "polB" in updated["NC_001"]
        assert "polB" in updated["NC_003"]

    def test_short_outlier_dropped(self):
        # NC_002 polB is 1/4 the median → dropped (default lo_mult = 1/3)
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000)},
            "NC_002": {"polB": _rec("YP_2", 250)},
            "NC_003": {"polB": _rec("YP_3", 990)},
        }
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_short_dropped"] == 1
        assert "polB" not in updated["NC_002"]

    def test_refseq_genome_protected(self, caplog):
        # NC_002 is RefSeq; its outlier marker is kept and a warning logged
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000)},
            "NC_002": {"polB": _rec("YP_2", 4500)},  # outlier
            "NC_003": {"polB": _rec("YP_3", 990)},
        }
        with caplog.at_level("WARNING"):
            updated, stats = remove_per_marker_length_outliers(
                genomes, refseq_genome_ids={"NC_002"},
            )
        assert stats["n_long_dropped"] == 0
        assert stats["n_refseq_protected"] == 1
        assert "polB" in updated["NC_002"]
        assert any("NC_002" in m and "protected" in m for m in caplog.messages)

    def test_per_marker_independent(self):
        # polB outlier in NC_002, MCP normal — only polB dropped for that genome
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000), "MCP": _rec("YP_2", 600)},
            "NC_002": {"polB": _rec("YP_3", 4500), "MCP": _rec("YP_4", 590)},
            "NC_003": {"polB": _rec("YP_5",  990), "MCP": _rec("YP_6", 610)},
        }
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 1
        assert "polB" not in updated["NC_002"]
        assert "MCP" in updated["NC_002"]

    def test_per_marker_median_recorded(self):
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000), "MCP": _rec("YP_2", 600)},
            "NC_002": {"polB": _rec("YP_3", 1010), "MCP": _rec("YP_4", 590)},
            "NC_003": {"polB": _rec("YP_5",  990), "MCP": _rec("YP_6", 610)},
        }
        _, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["per_marker_median"]["polB"] == 1000.0
        assert stats["per_marker_median"]["MCP"] == 600.0

    def test_disabled_hi_mult(self):
        # hi_mult=0 disables upper bound; the long sequence is kept.
        genomes = {
            "NC_001": {"polB": _rec("YP_1", 1000)},
            "NC_002": {"polB": _rec("YP_2", 9999)},
            "NC_003": {"polB": _rec("YP_3", 990)},
        }
        updated, stats = remove_per_marker_length_outliers(genomes, set(), hi_mult=0)
        assert stats["n_long_dropped"] == 0
        assert "polB" in updated["NC_002"]


# ---------------------------------------------------------------------------
# concatenate_aligned_markers
# ---------------------------------------------------------------------------

class TestConcatenateAlignedMarkers:
    def test_full_coverage_no_padding(self):
        aligned = {
            "polB": {
                "NC_001": _aln("NC_001", "MAAAA"),
                "NC_002": _aln("NC_002", "MBBBB"),
            },
            "MCP": {
                "NC_001": _aln("NC_001", "QQQ"),
                "NC_002": _aln("NC_002", "RRR"),
            },
        }
        concat, partition = concatenate_aligned_markers(
            aligned, ["NC_001", "NC_002"], ["polB", "MCP"],
        )
        assert str(concat["NC_001"].seq) == "MAAAAQQQ"
        assert str(concat["NC_002"].seq) == "MBBBBRRR"
        assert partition == {"polB": (1, 5), "MCP": (6, 8)}

    def test_missing_marker_gap_padded(self):
        # NC_002 lacks MCP — gets 3 dashes for that block.
        aligned = {
            "polB": {
                "NC_001": _aln("NC_001", "MAAAA"),
                "NC_002": _aln("NC_002", "MBBBB"),
            },
            "MCP": {
                "NC_001": _aln("NC_001", "QQQ"),
            },
        }
        concat, partition = concatenate_aligned_markers(
            aligned, ["NC_001", "NC_002"], ["polB", "MCP"],
        )
        assert str(concat["NC_001"].seq) == "MAAAAQQQ"
        assert str(concat["NC_002"].seq) == "MBBBB---"
        assert partition == {"polB": (1, 5), "MCP": (6, 8)}

    def test_empty_marker_skipped_in_partition(self):
        # No genome has 'MCP' at all → block_length 0 → marker omitted from partition.
        aligned = {
            "polB": {
                "NC_001": _aln("NC_001", "MAAAA"),
                "NC_002": _aln("NC_002", "MBBBB"),
            },
            "MCP": {},  # no aligned records
        }
        concat, partition = concatenate_aligned_markers(
            aligned, ["NC_001", "NC_002"], ["polB", "MCP"],
        )
        assert str(concat["NC_001"].seq) == "MAAAA"
        assert "MCP" not in partition
        assert partition == {"polB": (1, 5)}

    def test_marker_order_respected(self):
        aligned = {
            "polB": {"NC_001": _aln("NC_001", "AAAA")},
            "MCP":  {"NC_001": _aln("NC_001", "CC")},
            "hel":  {"NC_001": _aln("NC_001", "HHH")},
        }
        # MCP first, then hel, then polB
        concat, partition = concatenate_aligned_markers(
            aligned, ["NC_001"], ["MCP", "hel", "polB"],
        )
        assert str(concat["NC_001"].seq) == "CCHHHAAAA"
        assert partition == {"MCP": (1, 2), "hel": (3, 5), "polB": (6, 9)}

    def test_partition_coordinates_one_based_inclusive(self):
        aligned = {
            "polB": {"g1": _aln("g1", "AAAAA")},  # 5 cols → 1-5
            "MCP":  {"g1": _aln("g1", "CCC")},    # 3 cols → 6-8
        }
        _, partition = concatenate_aligned_markers(
            aligned, ["g1"], ["polB", "MCP"],
        )
        assert partition["polB"] == (1, 5)
        assert partition["MCP"] == (6, 8)

    def test_inconsistent_marker_alignment_raises(self):
        aligned = {
            "polB": {
                "g1": _aln("g1", "AAAA"),
                "g2": _aln("g2", "AAA"),  # different length — alignment broken
            },
        }
        with pytest.raises(ValueError, match="inconsistent column counts"):
            concatenate_aligned_markers(aligned, ["g1", "g2"], ["polB"])

    def test_genome_ids_order_preserved(self):
        aligned = {
            "polB": {
                "NC_001": _aln("NC_001", "AA"),
                "NC_002": _aln("NC_002", "BB"),
                "NC_003": _aln("NC_003", "CC"),
            },
        }
        concat, _ = concatenate_aligned_markers(
            aligned, ["NC_003", "NC_001", "NC_002"], ["polB"],
        )
        assert list(concat.keys()) == ["NC_003", "NC_001", "NC_002"]


# ---------------------------------------------------------------------------
# write_partition_file_nexus
# ---------------------------------------------------------------------------

class TestPartitionFile:
    def test_basic_nexus_format(self, tmp_path):
        partition = {"polB": (1, 1234), "MCP": (1235, 2456)}
        out = tmp_path / "partitions.nex"
        write_partition_file_nexus(partition, out)
        text = out.read_text()
        assert text.startswith("#nexus\n")
        assert "begin sets;" in text
        assert "charset polB = 1-1234;" in text
        assert "charset MCP = 1235-2456;" in text
        assert text.rstrip().endswith("end;")

    def test_marker_names_with_spaces_sanitized(self, tmp_path):
        partition = {"DNA polymerase": (1, 1234), "major capsid protein": (1235, 2456)}
        out = tmp_path / "partitions.nex"
        write_partition_file_nexus(partition, out)
        text = out.read_text()
        # NEXUS charset names may not contain whitespace
        assert "charset DNA_polymerase = 1-1234;" in text
        assert "charset major_capsid_protein = 1235-2456;" in text

    def test_marker_names_with_special_chars_sanitized(self, tmp_path):
        # The sanitizer strips anything outside [A-Za-z0-9_] for cross-tool
        # safety (some downstream tools choke on hyphens in charset names).
        partition = {"poly(A)-pol": (1, 100), "VLTF-3": (101, 200)}
        out = tmp_path / "partitions.nex"
        write_partition_file_nexus(partition, out)
        text = out.read_text()
        assert "charset poly_A_pol = 1-100;" in text
        assert "charset VLTF_3 = 101-200;" in text

    def test_creates_parent_directory(self, tmp_path):
        out = tmp_path / "deep" / "nested" / "partitions.nex"
        write_partition_file_nexus({"polB": (1, 100)}, out)
        assert out.exists()

    def test_safe_charset_name_helper(self):
        assert _safe_charset_name("DNA polymerase") == "DNA_polymerase"
        assert _safe_charset_name("RPB1/RPB2") == "RPB1_RPB2"
        assert _safe_charset_name("UL30") == "UL30"
        assert _safe_charset_name("") == "marker"
        assert _safe_charset_name("___") == "marker"
