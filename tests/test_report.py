"""Tests for report.py visualization helpers."""
import types

import pytest

from vfam_trees.report import _internal_label, _SHOW_INTERNAL_RANKS


def _make_clade(name, is_term, rank=None):
    """Return a minimal mock clade."""
    c = types.SimpleNamespace()
    c.name = name
    c.is_terminal = lambda: is_term
    if rank is not None:
        c._taxonomy_rank = rank
    return c


class TestInternalLabel:
    def test_terminal_always_shown(self):
        c = _make_clade("Species A|strain|ACC|host", True)
        assert _internal_label(c) == "Species A|strain|ACC|host"

    def test_terminal_empty_name(self):
        c = _make_clade("", True)
        assert _internal_label(c) == ""

    def test_internal_genus_shown(self):
        c = _make_clade("Alphainfluenzavirus", False, rank="genus")
        assert _internal_label(c) == "Alphainfluenzavirus"

    def test_internal_subgenus_shown(self):
        c = _make_clade("Sarbecovirus", False, rank="subgenus")
        assert _internal_label(c) == "Sarbecovirus"

    def test_internal_subfamily_shown(self):
        c = _make_clade("Orthocoronavirinae", False, rank="subfamily")
        assert _internal_label(c) == "Orthocoronavirinae"

    def test_internal_family_shown(self):
        c = _make_clade("Orthomyxoviridae", False, rank="family")
        assert _internal_label(c) == "Orthomyxoviridae"

    def test_internal_species_suppressed(self):
        c = _make_clade("Influenza A virus", False, rank="species")
        assert _internal_label(c) == ""

    def test_internal_unranked_suppressed(self):
        c = _make_clade("Influenza A virus", False, rank="")
        assert _internal_label(c) == ""

    def test_internal_no_rank_attr_suppressed(self):
        c = _make_clade("Some virus", False)  # no _taxonomy_rank set
        assert _internal_label(c) == ""

    def test_internal_order_suppressed(self):
        c = _make_clade("Nidovirales", False, rank="order")
        assert _internal_label(c) == ""

    def test_show_ranks_contents(self):
        assert _SHOW_INTERNAL_RANKS == {"genus", "subgenus", "subfamily", "family"}
