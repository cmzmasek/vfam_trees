"""Tests for vfam_trees.concat — pure helpers (no MAFFT / trimAl dependency)."""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees import concat as concat_module
from vfam_trees.concat import (
    _safe_charset_name,
    build_raw_concat,
    cluster_and_merge_genomes,
    concatenate_aligned_markers,
    identify_refseq_genomes,
    is_refseq_genome,
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

# Realistic per-marker variation so MAD on log-lengths is non-zero.
def _normal_polB(idx: int) -> SeqRecord:
    lengths = (980, 990, 995, 1000, 1000, 1005, 1010, 1015, 1020, 1025)
    return _rec(f"YP_polB_{idx}", lengths[idx % len(lengths)])


def _normal_MCP(idx: int) -> SeqRecord:
    lengths = (580, 590, 595, 600, 600, 605, 610, 615, 620, 625)
    return _rec(f"YP_MCP_{idx}", lengths[idx % len(lengths)])


def _genomes_with_normal_markers(n: int) -> dict[str, dict[str, SeqRecord]]:
    return {
        f"NC_{i:03d}": {"polB": _normal_polB(i), "MCP": _normal_MCP(i)}
        for i in range(n)
    }


class TestPerMarkerLengthOutliers:
    def test_normal_lengths_kept(self):
        genomes = _genomes_with_normal_markers(10)
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 0
        assert stats["n_short_dropped"] == 0
        assert all(len(m) == 2 for m in updated.values())

    def test_long_outlier_dropped(self):
        # One polB at ~30× median is clearly out of bounds at any reasonable k.
        genomes = _genomes_with_normal_markers(10)
        genomes["NC_999"] = {"polB": _rec("YP_huge", 30000)}
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 1
        assert "polB" not in updated["NC_999"]

    def test_short_outlier_dropped(self):
        # One polB at ~1/30 the median is dropped.
        genomes = _genomes_with_normal_markers(10)
        genomes["NC_999"] = {"polB": _rec("YP_tiny", 30)}
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_short_dropped"] == 1
        assert "polB" not in updated["NC_999"]

    def test_refseq_genome_protected(self, caplog):
        genomes = _genomes_with_normal_markers(10)
        genomes["NC_999"] = {"polB": _rec("YP_huge", 30000)}
        with caplog.at_level("WARNING"):
            updated, stats = remove_per_marker_length_outliers(
                genomes, refseq_genome_ids={"NC_999"},
            )
        assert stats["n_long_dropped"] == 0
        assert stats["n_refseq_protected"] == 1
        assert "polB" in updated["NC_999"]
        assert any("NC_999" in m and "protected" in m for m in caplog.messages)

    def test_per_marker_independent(self):
        # polB outlier on NC_999, MCP normal → only polB dropped for that genome
        genomes = _genomes_with_normal_markers(10)
        genomes["NC_999"] = {"polB": _rec("YP_huge", 30000), "MCP": _normal_MCP(0)}
        updated, stats = remove_per_marker_length_outliers(genomes, set())
        assert stats["n_long_dropped"] == 1
        assert "polB" not in updated["NC_999"]
        assert "MCP" in updated["NC_999"]

    def test_per_marker_median_recorded(self):
        genomes = _genomes_with_normal_markers(10)
        _, stats = remove_per_marker_length_outliers(genomes, set())
        # Both marker medians fall on the central pair (1000 / 600).
        assert stats["per_marker_median"]["polB"] == 1002.5
        assert stats["per_marker_median"]["MCP"] == 602.5

    def test_disabled_k_and_floor(self):
        # k=0 disables MAD; min_lo_mult=0 + max_hi_mult=0 disables the floor.
        # All together → filter is fully off, the long sequence is kept.
        genomes = _genomes_with_normal_markers(10)
        genomes["NC_999"] = {"polB": _rec("YP_huge", 30000)}
        updated, stats = remove_per_marker_length_outliers(
            genomes, set(), k=0, min_lo_mult=0, max_hi_mult=0,
        )
        assert stats["n_long_dropped"] == 0
        assert "polB" in updated["NC_999"]


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


# ---------------------------------------------------------------------------
# is_refseq_genome / identify_refseq_genomes
# ---------------------------------------------------------------------------

class TestIsRefseqGenome:
    def test_classic_refseq_prefixes(self):
        assert is_refseq_genome("NC_001234.1") is True
        assert is_refseq_genome("NZ_CP012345.1") is True
        assert is_refseq_genome("AC_000003.2") is True

    def test_bare_genbank_accession_rejected(self):
        assert is_refseq_genome("MN908947.3") is False
        assert is_refseq_genome("OP123456") is False

    def test_lowercase_prefix_rejected(self):
        assert is_refseq_genome("nc_001234.1") is False

    def test_empty_or_short_input(self):
        assert is_refseq_genome("") is False
        assert is_refseq_genome("X") is False
        assert is_refseq_genome(None) is False  # type: ignore[arg-type]


class TestIdentifyRefseqGenomes:
    def test_picks_only_refseq_accessions(self):
        species_genomes = {
            "Vaccinia virus": {
                "NC_006998.1": {"DNA polymerase": SeqRecord(Seq("M" * 1000), id="YP_1")},
                "MW123456.1":  {"DNA polymerase": SeqRecord(Seq("M" * 1000), id="YP_2")},
            },
            "Variola virus": {
                "NC_001611.1": {"DNA polymerase": SeqRecord(Seq("M" * 1000), id="YP_3")},
            },
        }
        refseqs = identify_refseq_genomes(species_genomes)
        assert refseqs == {"NC_006998.1", "NC_001611.1"}

    def test_empty_input(self):
        assert identify_refseq_genomes({}) == set()


# ---------------------------------------------------------------------------
# build_raw_concat
# ---------------------------------------------------------------------------

class TestBuildRawConcat:
    def test_full_coverage_concatenates_in_marker_order(self):
        proteins = {
            "polB": SeqRecord(Seq("AAAAA"), id="YP_1"),
            "MCP":  SeqRecord(Seq("CCC"),  id="YP_2"),
            "hel":  SeqRecord(Seq("HH"),   id="YP_3"),
        }
        rec = build_raw_concat(proteins, ["polB", "MCP", "hel"], "NC_001")
        assert rec.id == "NC_001"
        assert str(rec.seq) == "AAAAACCCHH"

    def test_marker_order_respected(self):
        proteins = {
            "polB": SeqRecord(Seq("AAA"), id="YP_1"),
            "MCP":  SeqRecord(Seq("CCC"), id="YP_2"),
        }
        rec = build_raw_concat(proteins, ["MCP", "polB"], "NC_001")
        assert str(rec.seq) == "CCCAAA"

    def test_missing_marker_omitted_no_padding(self):
        # Unaligned concat skips missing markers entirely (no gap padding).
        proteins = {"polB": SeqRecord(Seq("AAAAA"), id="YP_1")}  # MCP missing
        rec = build_raw_concat(proteins, ["polB", "MCP"], "NC_001")
        assert str(rec.seq) == "AAAAA"

    def test_empty_genome_yields_empty_seq(self):
        rec = build_raw_concat({}, ["polB", "MCP"], "NC_001")
        assert str(rec.seq) == ""
        assert rec.id == "NC_001"


# ---------------------------------------------------------------------------
# cluster_and_merge_genomes — uses monkeypatching to avoid MMseqs2 dependency
# ---------------------------------------------------------------------------

class TestClusterAndMergeGenomes:
    def _genome_dict(self, species_count=2, per_species=3):
        """Build a synthetic species_genomes dict with N species × M genomes."""
        species_genomes = {}
        for i in range(species_count):
            sp = f"Species{i}"
            genomes = {}
            # First genome of each species is a RefSeq (NC_ prefix)
            for j in range(per_species):
                if j == 0:
                    genome_id = f"NC_00{i}{j}.1"
                else:
                    genome_id = f"MN0000{i}{j}.1"
                genomes[genome_id] = {
                    "polB": SeqRecord(Seq("M" * 100), id=f"YP_{i}{j}_polB"),
                    "MCP":  SeqRecord(Seq("M" * 50),  id=f"YP_{i}{j}_MCP"),
                }
            species_genomes[sp] = genomes
        return species_genomes

    def _patch_no_op_clustering(self, monkeypatch):
        """Make absorb_into_refseqs and adaptive_cluster_species pass-throughs.

        Skips the actual MMseqs2 invocation so the test runs offline and
        deterministically.
        """
        def fake_absorb(records, refseq_ids, threshold, seq_type, work_dir,
                        clustering_tool="mmseqs2", **kwargs):
            return records, 0
        def fake_cluster(records, max_reps, threshold_min, threshold_max,
                         seq_type, work_dir, clustering_tool="mmseqs2", **kwargs):
            # Just return up to max_reps records
            return records[:max_reps], threshold_max
        monkeypatch.setattr(concat_module, "absorb_into_refseqs", fake_absorb)
        monkeypatch.setattr(concat_module, "adaptive_cluster_species", fake_cluster)

    def test_basic_end_to_end(self, tmp_path, monkeypatch):
        self._patch_no_op_clustering(monkeypatch)
        species_genomes = self._genome_dict(species_count=2, per_species=3)
        selected, stats = cluster_and_merge_genomes(
            species_genomes=species_genomes,
            marker_order=["polB", "MCP"],
            target_n=4,
            max_reps_per_species=3,
            threshold_min=0.7,
            threshold_max=0.99,
            work_dir=tmp_path,
        )
        # All 6 genomes are unique reps; merge picks 4 → 2 per species.
        assert len(selected) == 4
        assert stats["n_species_with_reps"] == 2
        assert stats["n_total_reps"] == 6
        assert stats["n_selected"] == 4

    def test_refseqs_prioritized_when_oversubscribed(self, tmp_path, monkeypatch):
        self._patch_no_op_clustering(monkeypatch)
        # 4 species × 3 genomes = 12 reps; target = 4 → must pick 1 per species,
        # and within each species the RefSeq (NC_ prefix) must win the slot.
        species_genomes = self._genome_dict(species_count=4, per_species=3)
        selected, stats = cluster_and_merge_genomes(
            species_genomes=species_genomes,
            marker_order=["polB", "MCP"],
            target_n=4,
            max_reps_per_species=3,
            threshold_min=0.7,
            threshold_max=0.99,
            work_dir=tmp_path,
        )
        assert len(selected) == 4
        # Every selected genome should be a RefSeq (NC_ prefix)
        assert all(s.startswith("NC_") for s in selected), selected

    def test_target_caps_total_selection(self, tmp_path, monkeypatch):
        self._patch_no_op_clustering(monkeypatch)
        species_genomes = self._genome_dict(species_count=3, per_species=5)  # 15 reps
        selected, stats = cluster_and_merge_genomes(
            species_genomes=species_genomes,
            marker_order=["polB", "MCP"],
            target_n=6,
            max_reps_per_species=5,
            threshold_min=0.7,
            threshold_max=0.99,
            work_dir=tmp_path,
        )
        assert len(selected) == 6

    def test_empty_species_skipped(self, tmp_path, monkeypatch):
        self._patch_no_op_clustering(monkeypatch)
        species_genomes = self._genome_dict(species_count=2, per_species=2)
        species_genomes["EmptySpecies"] = {}
        selected, stats = cluster_and_merge_genomes(
            species_genomes=species_genomes,
            marker_order=["polB", "MCP"],
            target_n=4,
            max_reps_per_species=2,
            threshold_min=0.7,
            threshold_max=0.99,
            work_dir=tmp_path,
        )
        # Only the 2 non-empty species contribute.
        assert stats["n_species_with_reps"] == 2

    def test_absorption_disabled_skips_call(self, tmp_path, monkeypatch):
        # When absorption is disabled, absorb_into_refseqs must NOT be invoked.
        called = []
        def boom_absorb(*a, **k):
            called.append(True)
            return a[0], 0
        def fake_cluster(records, max_reps, threshold_min, threshold_max,
                         seq_type, work_dir, clustering_tool="mmseqs2", **kwargs):
            return records[:max_reps], threshold_max
        monkeypatch.setattr(concat_module, "absorb_into_refseqs", boom_absorb)
        monkeypatch.setattr(concat_module, "adaptive_cluster_species", fake_cluster)

        species_genomes = self._genome_dict(species_count=1, per_species=2)
        cluster_and_merge_genomes(
            species_genomes=species_genomes,
            marker_order=["polB", "MCP"],
            target_n=2,
            max_reps_per_species=2,
            threshold_min=0.7,
            threshold_max=0.99,
            work_dir=tmp_path,
            refseq_absorption_enabled=False,
        )
        assert called == []
