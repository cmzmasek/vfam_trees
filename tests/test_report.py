"""Tests for report.py visualization helpers."""
import types

import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees import __version__ as VFAM_VERSION
from vfam_trees.report import (
    _build_tree_caption_info,
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


# ---------------------------------------------------------------------------
# _build_tree_caption_info — figure-legend caption embedded in tree PDFs/PNGs.
# Regression coverage for v1.2.10's caption helper, including the bug where
# pipeline_concat called save_tree_images without summary_row, leaving info=""
# and only the bold "(n external nodes)" suptitle visible.
# ---------------------------------------------------------------------------

class TestBuildTreeCaptionInfo:
    def _full_row(self) -> dict:
        # Mirrors a populated summary_row for label="100" (UFBoot range 0-100).
        return {
            "species_with_seqs": 89,
            "species_discovered": 120,
            "tree100_tree_tool": "iqtree",
            "tree100_tree_model": "GTR+G",
            "tree100_msa_tool": "mafft",
            "tree100_msa_options": "--auto",
            "tree100_msa_length": 10432,
            "tree100_msa_gap_pct": 14.0,
            "tree100_trim_tool": "trimal",
            "tree100_trim_options": "-automated1",
            "tree100_support_type": "UFBoot",
            "tree100_support_median": 95,
            "tree100_support_q1": 80,
            "tree100_support_q3": 100,
            "tree100_n_genera": 4,
            "tree100_n_subfamilies": 2,
            "tree100_n_outliers_removed": 3,
        }

    def test_none_summary_row_returns_empty(self):
        assert _build_tree_caption_info(None, "100") == ""

    def test_empty_summary_row_returns_empty(self):
        assert _build_tree_caption_info({}, "100") == ""

    def test_full_row_two_lines(self):
        info = _build_tree_caption_info(self._full_row(), "100")
        assert "\n" in info
        line1, line2 = info.split("\n", 1)
        # Line 1: method · alignment
        assert "ML phylogeny" in line1
        assert "IQ-TREE" in line1               # tree_tool display mapping
        assert "GTR+G" in line1
        assert "10,432 cols" in line1           # thousands separator
        assert "mafft --auto" in line1
        assert "trimal -automated1" in line1
        assert "14% gaps" in line1
        # Line 2: support · taxonomy · species · outliers · version
        assert "UFBoot median 95 (IQR 80–100)" in line2
        assert "4 genera / 2 subfamilies" in line2
        assert "89 / 120 species" in line2
        assert "3 outliers removed" in line2
        assert f"vfam_trees v{VFAM_VERSION}" in line2

    def test_fasttree_support_uses_fractional_format(self):
        # FastTree SH-like values are in [0,1] — must format with two decimals
        # rather than rounding to "median 1 (IQR 1–1)".
        row = {
            "tree500_tree_tool": "fasttree",
            "tree500_tree_model": "GTR+G",
            "tree500_support_type": "SH_like",
            "tree500_support_median": 0.92,
            "tree500_support_q1": 0.80,
            "tree500_support_q3": 0.98,
        }
        info = _build_tree_caption_info(row, "500")
        assert "FastTree" in info  # display mapping (not "Fasttree")
        assert "median 0.92" in info
        assert "IQR 0.80–0.98" in info

    def test_ufboot_support_uses_integer_format(self):
        row = {
            "tree100_support_type": "UFBoot",
            "tree100_support_median": 95,
            "tree100_support_q1": 80,
            "tree100_support_q3": 100,
        }
        info = _build_tree_caption_info(row, "100")
        # Integer formatting — no decimals when value > 1
        assert "median 95" in info
        assert "0.95" not in info

    def test_singular_pluralization(self):
        row = {
            "tree100_n_genera": 1,
            "tree100_n_subfamilies": 1,
            "tree100_n_outliers_removed": 1,
        }
        info = _build_tree_caption_info(row, "100")
        assert "1 genus" in info
        assert "1 subfamily" in info
        assert "1 outlier removed" in info
        assert "outliers" not in info  # singular only

    def test_plural_pluralization(self):
        row = {
            "tree100_n_genera": 5,
            "tree100_n_subfamilies": 2,
            "tree100_n_outliers_removed": 3,
        }
        info = _build_tree_caption_info(row, "100")
        assert "5 genera" in info
        assert "2 subfamilies" in info
        assert "3 outliers removed" in info

    def test_zero_outliers_omitted(self):
        row = {"tree100_n_outliers_removed": 0}
        info = _build_tree_caption_info(row, "100")
        assert "outlier" not in info

    def test_concat_markers_render(self):
        row = {
            "tree100_tree_tool": "iqtree",
            "tree100_msa_tool": "mafft",
            "tree100_concat_n_markers_used": 8,
            "tree100_concat_n_markers_target": 10,
            "tree100_msa_length": 4321,
        }
        info = _build_tree_caption_info(row, "100")
        assert "concat 8/10 markers" in info

    def test_label_500_uses_tree500_prefix(self):
        # The same key under tree100_* must NOT leak into a label="500" caption.
        row = {
            "tree100_tree_tool": "iqtree",
            "tree100_tree_model": "GTR+G",
            "tree500_tree_tool": "fasttree",
            "tree500_tree_model": "JTT",
        }
        info = _build_tree_caption_info(row, "500")
        assert "FastTree" in info
        assert "JTT" in info
        assert "IQ-TREE" not in info

    def test_version_always_included(self):
        # Even with an otherwise-empty row (one truthy key to bypass the
        # `if not summary_row` early-return), the version line should appear.
        info = _build_tree_caption_info({"_marker": 1}, "100")
        assert f"vfam_trees v{VFAM_VERSION}" in info


# ---------------------------------------------------------------------------
# save_tree_images — accepts summary_row, threaded through to _draw_tree_fig.
# Regression for the concat-pipeline omission that blanked the caption.
# ---------------------------------------------------------------------------

class TestSaveTreeImagesSummaryRow:
    def test_summary_row_passed_through_to_draw(self, tmp_path, monkeypatch):
        from vfam_trees import report as report_mod

        seen: dict[str, dict | None] = {}

        def _fake_draw(tree, family, label="100", **kwargs):
            seen[label] = kwargs.get("summary_row")
            return None  # skip image rendering

        monkeypatch.setattr(report_mod, "_draw_tree_fig", _fake_draw)
        monkeypatch.setattr(report_mod, "_draw_unrooted_tree_fig", _fake_draw)

        # A minimal stand-in for a Bio.Phylo Tree (only needs to be non-None
        # to enter the per-label loop in save_tree_images).
        fake_trees = {"100": object(), "500": object()}
        sentinel = {"tree100_tree_tool": "iqtree", "marker": "ok"}

        report_mod.save_tree_images(
            family="Testviridae",
            output_dir=tmp_path,
            bio_trees=fake_trees,
            summary_row=sentinel,
        )

        assert seen.get("100") is sentinel
        assert seen.get("500") is sentinel
