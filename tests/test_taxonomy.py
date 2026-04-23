"""Tests for vfam_trees.taxonomy — lca_min_rank filtering."""
import pytest

from vfam_trees.taxonomy import _lineage_reaches_rank, _RANK_DEPTH


def _e(rank):
    return {"name": f"taxon_{rank}", "rank": rank}


# ---------------------------------------------------------------------------
# _lineage_reaches_rank
# ---------------------------------------------------------------------------

def test_none_mode_always_true():
    assert _lineage_reaches_rank([], "none") is True
    assert _lineage_reaches_rank([_e("family")], "none") is True


def test_empty_min_rank_always_true():
    assert _lineage_reaches_rank([_e("family")], "") is True


def test_genus_reached_when_genus_present():
    lin = [_e("family"), _e("genus"), _e("species")]
    assert _lineage_reaches_rank(lin, "genus") is True


def test_genus_not_reached_when_only_subfamily():
    lin = [_e("family"), _e("subfamily")]
    assert _lineage_reaches_rank(lin, "genus") is False


def test_species_reached_when_species_present():
    lin = [_e("family"), _e("genus"), _e("species")]
    assert _lineage_reaches_rank(lin, "species") is True


def test_species_not_reached_when_only_genus():
    lin = [_e("family"), _e("genus")]
    assert _lineage_reaches_rank(lin, "species") is False


def test_subfamily_reached_when_subfamily_present():
    lin = [_e("family"), _e("subfamily")]
    assert _lineage_reaches_rank(lin, "subfamily") is True


def test_no_rank_entries_do_not_satisfy_min():
    lin = [_e("no rank"), _e("clade")]
    assert _lineage_reaches_rank(lin, "genus") is False


def test_subspecies_satisfies_genus_min():
    # subspecies depth > genus depth
    lin = [_e("subspecies")]
    assert _lineage_reaches_rank(lin, "genus") is True


def test_unknown_min_rank_treated_as_zero_always_true():
    lin = [_e("family")]
    assert _lineage_reaches_rank(lin, "superkingdom_xyz") is True


def test_rank_depth_ordering():
    assert _RANK_DEPTH["genus"] < _RANK_DEPTH["species"]
    assert _RANK_DEPTH["family"] < _RANK_DEPTH["genus"]
    assert _RANK_DEPTH["subfamily"] < _RANK_DEPTH["genus"]
