"""Tests for report.py visualization helpers."""
import types

import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees.report import (
    _draw_tree_fig,
    _draw_unrooted_tree_fig,
    _internal_label,
    _SHOW_INTERNAL_RANKS,
)


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


def _make_unladdered_tree():
    """Return a tree whose child ordering will change under ladderize(reverse=True)."""
    # Left child has 1 leaf, right child has 3 leaves → ladderize(reverse=True)
    # reorders so the larger subtree comes first.
    small = Clade(name="small", clades=[Clade(name="l0", branch_length=0.1)])
    big = Clade(name="big", clades=[
        Clade(name="l1", branch_length=0.1),
        Clade(name="l2", branch_length=0.1),
        Clade(name="l3", branch_length=0.1),
    ])
    root = Clade(clades=[small, big])
    return Tree(root=root)


class TestDrawingFunctionsDoNotMutateInput:
    def test_draw_tree_fig_does_not_ladderize_input(self):
        tree = _make_unladdered_tree()
        original_order = [c.name for c in tree.root.clades]
        assert original_order == ["small", "big"]

        fig = _draw_tree_fig(tree, family="F", label="100")
        try:
            # Caller's tree must retain its original child order
            assert [c.name for c in tree.root.clades] == original_order
        finally:
            if fig is not None:
                import matplotlib.pyplot as plt
                plt.close(fig)

    def test_draw_unrooted_tree_fig_does_not_ladderize_input(self):
        tree = _make_unladdered_tree()
        original_order = [c.name for c in tree.root.clades]

        fig = _draw_unrooted_tree_fig(tree, family="F", label="100")
        try:
            assert [c.name for c in tree.root.clades] == original_order
        finally:
            if fig is not None:
                import matplotlib.pyplot as plt
                plt.close(fig)
