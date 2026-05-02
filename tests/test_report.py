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


# ---------------------------------------------------------------------------
# generate_family_report — concat-mode marker-coverage page
# CONCAT_DESIGN.md §6.3
# ---------------------------------------------------------------------------

class TestMarkerCoveragePage:
    """Verifies the new optional per-marker coverage page."""

    @staticmethod
    def _read_pdf_page_count(path) -> int:
        # Best-effort page counter — counts /Type /Page objects in the PDF.
        # Avoids a hard dependency on PyPDF2 / pypdf in the test env.
        text = path.read_bytes()
        return text.count(b"/Type /Page") - text.count(b"/Type /Pages")

    def _basic_summary_row(self) -> dict:
        return {
            "ncbi_taxid": 10240,
            "lineage": "Viruses",
            "molecule_region": "protein, concatenated",
            "species_discovered": 5,
            "species_with_seqs": 5,
            "tree500_seq_type": "protein",
            "tree500_msa_tool": "mafft",
            "tree500_msa_options": "",
            "tree500_tree_tool": "fasttree",
            "tree500_tree_model": "LG+G",
            "tree500_tree_options": "",
            "tree500_leaves": 100,
            "tree500_msa_length": 4321,
            "tree500_msa_gap_pct": 5.0,
            "tree500_cluster_thresh_min": 0.85,
            "tree500_cluster_thresh_max": 0.99,
            "tree500_support_type": "SH_like",
            "tree500_support_median": 92,
            "tree500_support_iqr": 8,
            "tree100_seq_type": "protein",
            "tree100_msa_tool": "mafft",
            "tree100_msa_options": "",
            "tree100_tree_tool": "iqtree",
            "tree100_tree_model": "MFP",
            "tree100_tree_options": "-p partitions.nex -B 1000",
            "tree100_leaves": 30,
            "tree100_msa_length": 4321,
            "tree100_msa_gap_pct": 5.0,
            "tree100_cluster_thresh_min": 0.85,
            "tree100_cluster_thresh_max": 0.99,
            "tree100_support_type": "SH_aLRT",
            "tree100_support_median": 95,
            "tree100_support_iqr": 4,
        }

    def test_coverage_page_added_when_marker_coverage_present(self, tmp_path):
        from vfam_trees.report import generate_family_report

        out_with    = tmp_path / "with_coverage.pdf"
        out_without = tmp_path / "without_coverage.pdf"
        common = dict(
            family="Poxviridae",
            summary_row=self._basic_summary_row(),
            seq_lengths=[200, 250, 300, 350],
            tree_seq_lengths={"500": [4321, 4321, 4321], "100": [4321, 4321]},
        )
        generate_family_report(output_pdf=out_without, **common)
        generate_family_report(
            output_pdf=out_with,
            **common,
            marker_coverage={
                "500": {"polB": 100, "MCP": 98, "hel": 95},
                "100": {"polB": 30, "MCP": 30, "hel": 28},
            },
            concat_min_fraction=0.7,
        )
        assert out_with.exists() and out_without.exists()
        # The coverage version must have exactly one more page.
        diff = self._read_pdf_page_count(out_with) - self._read_pdf_page_count(out_without)
        assert diff == 1, f"expected 1 extra page, got {diff}"

    def test_coverage_page_skipped_when_marker_coverage_empty(self, tmp_path):
        from vfam_trees.report import generate_family_report

        out_a = tmp_path / "a.pdf"
        out_b = tmp_path / "b.pdf"
        common = dict(
            family="Poxviridae",
            summary_row=self._basic_summary_row(),
            seq_lengths=[200, 250, 300],
        )
        generate_family_report(output_pdf=out_a, **common)
        generate_family_report(output_pdf=out_b, **common, marker_coverage={})
        # Empty marker_coverage → no extra page
        assert (self._read_pdf_page_count(out_a) ==
                self._read_pdf_page_count(out_b))

    def test_coverage_page_handles_missing_label(self, tmp_path):
        # Only tree_500 has marker coverage data — should still render one panel.
        from vfam_trees.report import generate_family_report

        out = tmp_path / "single_label.pdf"
        generate_family_report(
            family="Poxviridae",
            output_pdf=out,
            summary_row=self._basic_summary_row(),
            seq_lengths=[200, 250],
            marker_coverage={"500": {"polB": 100, "MCP": 95}},
            concat_min_fraction=0.7,
        )
        assert out.exists()
        # PDF should be parseable and contain at least one /Type /Page
        assert self._read_pdf_page_count(out) >= 1
