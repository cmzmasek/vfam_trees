"""Tests for vfam_trees.branch_outliers — the median + factor × MAD criterion
shared by single-protein and concat pipelines.
"""
from __future__ import annotations

import pytest

from vfam_trees.branch_outliers import (
    branch_length_stats,
    find_branch_length_outliers,
)


def _write_newick(tmp_path, name: str, newick: str):
    path = tmp_path / f"{name}.nwk"
    path.write_text(newick)
    return path


# ---------------------------------------------------------------------------
# branch_length_stats
# ---------------------------------------------------------------------------

class TestBranchLengthStats:
    def test_simple_tree_yields_correct_map(self, tmp_path):
        nwk = _write_newick(tmp_path, "ok",
                            "(A:0.1,B:0.1,(C:0.1,D:0.1):0.05);")
        stats = branch_length_stats(nwk)
        assert set(stats["bl_map"].keys()) == {"A", "B", "C", "D"}
        assert stats["median"] == pytest.approx(0.1)
        # MAD of identical values is 0
        assert stats["mad"] == pytest.approx(0.0)

    def test_outlier_inflates_neither_median_nor_mad_too_much(self, tmp_path):
        nwk = _write_newick(tmp_path, "outlier",
                            "(A:0.1,B:0.1,C:0.1,(D:0.1,Z:5.0):0.05);")
        stats = branch_length_stats(nwk)
        # Median is robust: should still be near 0.1 (not 5.0)
        assert stats["median"] == pytest.approx(0.1)
        # MAD is also robust: median absolute deviation from 0.1 is 0
        # for {0.1,0.1,0.1,0.1} but the outlier 5.0 contributes 4.9 to one;
        # of 5 values four are 0.0, one is 4.9 → MAD = 0.0
        assert stats["mad"] == pytest.approx(0.0)
        assert stats["bl_map"]["Z"] == pytest.approx(5.0)

    def test_too_few_branches_returns_zero_stats(self, tmp_path):
        nwk = _write_newick(tmp_path, "small", "(A:0.1,B:0.1);")
        stats = branch_length_stats(nwk)
        assert stats["median"] == 0.0
        assert stats["mad"] == 0.0

    def test_unparseable_file_returns_zero_stats(self, tmp_path):
        bad = tmp_path / "garbage.nwk"
        bad.write_text("not a newick string at all")
        stats = branch_length_stats(bad)
        assert stats["median"] == 0.0
        assert stats["mad"] == 0.0
        assert stats["bl_map"] == {}


# ---------------------------------------------------------------------------
# find_branch_length_outliers
# ---------------------------------------------------------------------------

class TestFindBranchLengthOutliers:
    def test_no_outlier_returns_empty(self, tmp_path):
        nwk = _write_newick(tmp_path, "uniform",
                            "(A:0.1,B:0.1,C:0.1,(D:0.1,E:0.1):0.05);")
        # All zero MAD → returns empty regardless of factor
        assert find_branch_length_outliers(nwk, factor=20.0) == set()

    def test_obvious_outlier_detected(self, tmp_path):
        # 5 leaves: four around 0.1, one at 5.0.  MAD with one outlier may be
        # zero, in which case the outlier is technically not detected by
        # threshold = median + factor*0 = median.  Use mixed branch lengths
        # to give a non-zero MAD so the threshold is meaningful.
        nwk = _write_newick(tmp_path, "spread",
                            "(A:0.05,B:0.10,C:0.15,(D:0.20,Z:5.0):0.05);")
        outliers = find_branch_length_outliers(nwk, factor=5.0)
        assert "Z" in outliers
        assert "A" not in outliers and "D" not in outliers

    def test_returns_strings(self, tmp_path):
        nwk = _write_newick(tmp_path, "spread",
                            "(A:0.05,B:0.10,C:0.15,(D:0.20,Z:5.0):0.05);")
        outliers = find_branch_length_outliers(nwk, factor=5.0)
        assert all(isinstance(o, str) for o in outliers)

    def test_factor_controls_strictness(self, tmp_path):
        nwk = _write_newick(tmp_path, "spread",
                            "(A:0.05,B:0.10,C:0.15,(D:0.20,Z:5.0):0.05);")
        # Very strict factor — Z still flagged
        assert "Z" in find_branch_length_outliers(nwk, factor=2.0)
        # Very lenient factor — Z still flagged because it's so far above
        assert "Z" in find_branch_length_outliers(nwk, factor=10.0)

    def test_nonexistent_file_returns_empty(self, tmp_path):
        assert find_branch_length_outliers(tmp_path / "missing.nwk", factor=20.0) == set()


# ---------------------------------------------------------------------------
# Backward compat: pipeline.py re-exports these as private names
# ---------------------------------------------------------------------------

def test_pipeline_private_aliases_exist():
    """The single-protein code uses underscore-prefixed names; verify the
    re-exports are still in place after the extraction.
    """
    from vfam_trees.pipeline import (
        _branch_length_stats,
        _find_branch_length_outliers,
    )
    assert _branch_length_stats is branch_length_stats
    assert _find_branch_length_outliers is find_branch_length_outliers
