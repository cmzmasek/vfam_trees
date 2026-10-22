"""Tests for vfam_trees.branch_outliers — the median + factor × MAD criterion
shared by single-protein and concat pipelines.
"""
from __future__ import annotations

import pytest

from vfam_trees.branch_outliers import (
    MIN_TREE_LEAVES,
    QUICK_STAGE,
    branch_length_stats,
    find_branch_length_outliers,
    outlier_round,
    select_outliers,
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
# Both pipelines make the outlier decision through the shared helpers
# ---------------------------------------------------------------------------

def test_pipelines_use_the_shared_round_helpers():
    from vfam_trees import pipeline, pipeline_concat
    for mod in (pipeline, pipeline_concat):
        assert mod.select_outliers is select_outliers
        assert mod.outlier_round is outlier_round


class _Log:
    def __init__(self):
        self.lines = []

    def info(self, msg, *a):
        self.lines.append(("info", msg % a))

    def warning(self, msg, *a):
        self.lines.append(("warning", msg % a))


def _star(tmp_path, name, long_tips=("L0",), n=60):
    """n ordinary tips (0.05–0.07) plus very long ones (5.0)."""
    tips = [f"T{i}:{0.05 + 0.0003 * i:.4f}" for i in range(n)]
    tips += [f"{t}:5.0" for t in long_tips]
    path = tmp_path / name
    path.write_text("(" + ",".join(tips) + ");")
    return path


def _select(path, protected=frozenset(), n_current=61, min_seqs=40, log=None, stage=""):
    return select_outliers(
        path, 20.0, set(protected), n_current, min_seqs,
        log=log or _Log(), label="500", iteration=1,
        display=lambda s: f"<{s}>", noun="sequence", stage=stage,
    )


class TestSelectOutliers:
    def test_removes_long_tips_and_logs_names(self, tmp_path):
        log = _Log()
        out = _select(_star(tmp_path, "a.nwk", ("L0", "L1")), n_current=62, log=log)
        assert out == {"L0", "L1"}
        removed = [m for lvl, m in log.lines if "removing outlier" in m]
        assert len(removed) == 2 and "<L0>" in removed[0]
        assert all(m.startswith("tree_500 iteration 1:") for _, m in log.lines)

    def test_no_outliers_means_stop(self, tmp_path):
        path = tmp_path / "b.nwk"
        path.write_text("(" + ",".join(f"T{i}:{0.05 + 0.001 * i}" for i in range(50)) + ");")
        log = _Log()
        assert _select(path, log=log) == set()
        assert "no branch-length outliers" in log.lines[-1][1]

    def test_protected_are_kept_and_warned(self, tmp_path):
        log = _Log()
        out = _select(_star(tmp_path, "c.nwk", ("L0",)), protected={"L0"}, log=log)
        assert out == set()
        assert any(lvl == "warning" and "<L0>" in m and "KEEPING" in m for lvl, m in log.lines)
        assert "only protected outliers remain" in log.lines[-1][1]

    def test_min_seqs_stops_removal(self, tmp_path):
        log = _Log()
        out = _select(_star(tmp_path, "d.nwk", ("L0", "L1")), n_current=41, min_seqs=40, log=log)
        assert out == set()
        assert "would leave only 39 sequences (min_seqs=40)" in log.lines[-1][1]

    def test_stage_tags_every_line(self, tmp_path):
        log = _Log()
        _select(_star(tmp_path, "e.nwk"), log=log, stage=QUICK_STAGE)
        assert log.lines and all(m.startswith(f"tree_500 iteration 1{QUICK_STAGE}:") for _, m in log.lines)


class TestOutlierRound:
    """The quick tree only ever saves the configured search; the final tree is
    always the configured one, checked the same way."""

    def _builders(self, quick_result, final_result):
        calls = []

        def build_tree():
            calls.append("configured")
            return "configured.nwk"

        def build_quick():
            calls.append("quick")
            return "quick.nwk"

        def decide(tf, stage):
            calls.append(("decide", tf, stage))
            return set(quick_result if tf == "quick.nwk" else final_result)

        return calls, build_tree, build_quick, decide

    def test_quick_tree_outliers_skip_the_configured_search(self):
        calls, bt, bq, dec = self._builders({"x"}, set())
        tree, remove = outlier_round(bt, dec, check=True, n_current=50, build_quick_tree=bq)
        assert (tree, remove) == (None, {"x"})
        assert "configured" not in calls

    def test_clean_quick_tree_runs_and_checks_the_configured_tree(self):
        calls, bt, bq, dec = self._builders(set(), {"y"})
        tree, remove = outlier_round(bt, dec, check=True, n_current=50, build_quick_tree=bq)
        assert (tree, remove) == ("configured.nwk", {"y"})
        assert calls == ["quick", ("decide", "quick.nwk", QUICK_STAGE),
                         "configured", ("decide", "configured.nwk", "")]

    def test_no_check_builds_only_the_configured_tree(self):
        calls, bt, bq, dec = self._builders({"x"}, {"y"})
        tree, remove = outlier_round(bt, dec, check=False, n_current=50, build_quick_tree=bq)
        assert (tree, remove) == ("configured.nwk", set())
        assert calls == ["configured"]

    def test_without_quick_builder_behaves_as_before(self):
        calls, bt, _bq, dec = self._builders({"x"}, {"y"})
        tree, remove = outlier_round(bt, dec, check=True, n_current=50)
        assert (tree, remove) == ("configured.nwk", {"y"})
        assert calls == ["configured", ("decide", "configured.nwk", "")]

    def test_quick_removal_that_would_leave_too_few_falls_through(self):
        n = MIN_TREE_LEAVES + 1
        calls, bt, bq, dec = self._builders({"a", "b"}, set())
        tree, remove = outlier_round(bt, dec, check=True, n_current=n, build_quick_tree=bq)
        assert tree == "configured.nwk" and remove == set()
        assert "configured" in calls

    def test_real_trees_quick_then_clean(self, tmp_path):
        quick = _star(tmp_path, "quick.nwk", ("L0",))
        final = tmp_path / "final.nwk"
        final.write_text("(" + ",".join(f"T{i}:{0.05 + 0.0003 * i:.4f}" for i in range(60)) + ");")
        decide = lambda tf, stage: _select(tf, n_current=61, stage=stage)
        tree, remove = outlier_round(lambda: final, decide, check=True, n_current=61,
                                     build_quick_tree=lambda: quick)
        assert (tree, remove) == (None, {"L0"})


# ---------------------------------------------------------------------------
# IQ-TREE writes zero-length tips as its minimum branch length (1e-6), not 0.0.
# Those must be left out of the median / MAD exactly like FastTree's zeros,
# or a tree with many near-identical tips gets a tiny threshold and ordinary
# tips are flagged.
# ---------------------------------------------------------------------------

def _outbreak_like(tmp_path, name, zero_value):
    """30 identical tips, 30 near-identical ones, 40 ordinary ones."""
    tips = [f"Z{i}:{zero_value}" for i in range(30)]
    tips += [f"N{i}:{0.0001 + 0.00001 * i:.6f}" for i in range(30)]
    tips += [f"T{i}:{0.05 + 0.004 * i:.4f}" for i in range(40)]
    path = tmp_path / name
    path.write_text("(" + ",".join(tips) + ");")
    return path


class TestIqtreeMinimumBranchLength:
    def test_iqtree_minimum_treated_like_fasttree_zero(self, tmp_path):
        fasttree = branch_length_stats(_outbreak_like(tmp_path, "ft.nwk", "0.0"))
        iqtree = branch_length_stats(_outbreak_like(tmp_path, "iq.nwk", "0.0000010000"))
        assert iqtree["median"] == pytest.approx(fasttree["median"])
        assert iqtree["mad"] == pytest.approx(fasttree["mad"])

    def test_identical_tips_do_not_flag_ordinary_tips(self, tmp_path):
        nwk = _outbreak_like(tmp_path, "iq.nwk", "0.0000010000")
        # Counting the 1e-6 tips puts the median among the near-identical
        # ones (~3e-4), so median + 20 × MAD fell below every ordinary tip
        # and all 40 were flagged.
        assert find_branch_length_outliers(nwk, factor=20.0) == set()
