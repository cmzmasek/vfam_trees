"""Tests for vfam_trees.summary helpers."""
import csv

from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees.summary import (
    STATUS_COLUMNS,
    build_status_row,
    build_summary_row,
    compute_msa_stats,
    compute_seqlen_stats,
    compute_support_stats,
    format_molecule_region,
    load_family_annotations,
    write_status_row,
    write_summary_row,
)


# ---------------------------------------------------------------------------
# compute_seqlen_stats
# ---------------------------------------------------------------------------

class TestComputeSeqlenStats:
    def test_empty_returns_empty_strings(self):
        assert compute_seqlen_stats([]) == {
            "min": "", "q1": "", "median": "", "mean": "", "q3": "", "max": "", "iqr": "",
        }

    def test_single_value(self):
        stats = compute_seqlen_stats([100])
        assert stats["min"] == 100
        assert stats["max"] == 100
        assert stats["median"] == 100

    def test_three_values_quantiles_computed(self):
        stats = compute_seqlen_stats([100, 200, 300])
        assert stats["min"] == 100
        assert stats["max"] == 300
        assert stats["median"] == 200

    def test_larger_dataset(self):
        lengths = list(range(1, 101))
        stats = compute_seqlen_stats(lengths)
        assert stats["min"] == 1
        assert stats["max"] == 100
        assert stats["median"] == 50.5
        assert stats["iqr"] > 0

    def test_iqr_positive_when_data_spread(self):
        stats = compute_seqlen_stats([10, 20, 30, 40, 50])
        assert stats["iqr"] > 0


# ---------------------------------------------------------------------------
# compute_support_stats
# ---------------------------------------------------------------------------

class TestComputeSupportStats:
    def _tree_with_conf(self, values):
        leaves = [Clade(name=f"l{i}") for i in range(len(values) + 1)]
        clades = [Clade(clades=[leaves[i], leaves[i+1]], confidence=v)
                  for i, v in enumerate(values)]
        def chain(cs):
            return cs[0] if len(cs) == 1 else Clade(clades=[cs[0], chain(cs[1:])])
        return Tree(root=chain(clades))

    def test_empty_tree(self):
        tree = Tree(root=Clade(clades=[Clade(name="a"), Clade(name="b")]))
        stats = compute_support_stats(tree)
        assert stats == {k: "" for k in ("min", "q1", "median", "q3", "max", "iqr")}

    def test_range_captured(self):
        tree = self._tree_with_conf([10, 50, 90, 95, 99])
        stats = compute_support_stats(tree)
        assert stats["min"] == 10
        assert stats["max"] == 99

    def test_terminals_excluded_from_support(self):
        # Give a terminal a confidence value; it must not be included.
        leaves = [Clade(name=f"l{i}", confidence=42) for i in range(2)]
        root = Clade(clades=leaves, confidence=80)
        tree = Tree(root=root)
        stats = compute_support_stats(tree)
        assert stats["min"] == 80
        assert stats["max"] == 80


# ---------------------------------------------------------------------------
# compute_msa_stats
# ---------------------------------------------------------------------------

class TestComputeMsaStats:
    def test_empty_file(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        stats = compute_msa_stats(empty)
        assert stats == {"length": "", "gap_pct": ""}

    def test_simple_alignment(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">a\nATGC\n>b\nA-GC\n>c\nATGC\n")
        stats = compute_msa_stats(aln)
        assert stats["length"] == 4
        # 1 gap out of 12 chars = 8.3%
        assert stats["gap_pct"] == 8.3

    def test_no_gaps(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">a\nATGC\n>b\nATGC\n")
        stats = compute_msa_stats(aln)
        assert stats["gap_pct"] == 0.0

    def test_all_gaps(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">a\n----\n>b\n----\n")
        stats = compute_msa_stats(aln)
        assert stats["gap_pct"] == 100.0


# ---------------------------------------------------------------------------
# format_molecule_region
# ---------------------------------------------------------------------------

class TestFormatMoleculeRegion:
    def test_nuc_whole_genome(self):
        assert format_molecule_region("nucleotide", "whole_genome", None) == \
            "nucleotide, whole genome"

    def test_protein_gene(self):
        assert format_molecule_region("protein", "RdRp", None) == \
            "protein, gene: RdRp"

    def test_segment_overrides_region(self):
        assert format_molecule_region("nucleotide", "whole_genome", "segment S") == \
            "nucleotide, segment S"


# ---------------------------------------------------------------------------
# build_summary_row / write_summary_row
# ---------------------------------------------------------------------------

class TestBuildSummaryRow:
    def test_minimal_row_has_family_and_molecule(self):
        row = build_summary_row(
            family="Flaviviridae",
            family_taxid=11050,
            family_lineage=[{"name": "Viruses"}, {"name": "Riboviria"}],
            seq_type="nucleotide",
            region="whole_genome",
            segment=None,
            n_species_discovered=100,
            n_species_with_seqs=90,
            seqlen_stats={"min": 9000, "max": 11000},
            tree_stats={},
        )
        assert row["family"] == "Flaviviridae"
        assert row["ncbi_taxid"] == 11050
        assert "Viruses; Riboviria" in row["lineage"]
        assert row["molecule_region"] == "nucleotide, whole genome"

    def test_tree_stats_flattened_with_prefix(self):
        row = build_summary_row(
            family="F", family_taxid=1, family_lineage=[],
            seq_type="nucleotide", region="whole_genome", segment=None,
            n_species_discovered=0, n_species_with_seqs=0, seqlen_stats={},
            tree_stats={
                "500": {"leaves": 487, "tree_tool": "fasttree"},
                "100": {"leaves": 95, "tree_tool": "iqtree"},
            },
        )
        assert row["tree500_leaves"] == 487
        assert row["tree100_leaves"] == 95
        assert row["tree500_tree_tool"] == "fasttree"
        assert row["tree100_tree_tool"] == "iqtree"


class TestLoadFamilyAnnotations:
    def test_missing_file_returns_empty_dict(self, tmp_path):
        assert load_family_annotations(tmp_path / "nope.tsv") == {}

    def test_basic_load_case_insensitive_key(self, tmp_path):
        p = tmp_path / "ann.tsv"
        p.write_text(
            "family\tbaltimore_class\thost_range\n"
            "Flaviviridae\tIV\tVertebrates\n"
            "Poxviridae\tI\tVertebrates\n"
        )
        ann = load_family_annotations(p)
        # Keys are lowercased
        assert "flaviviridae" in ann
        assert ann["flaviviridae"]["baltimore_class"] == "IV"
        assert ann["poxviridae"]["baltimore_class"] == "I"

    def test_blank_family_rows_skipped(self, tmp_path):
        p = tmp_path / "ann.tsv"
        p.write_text(
            "family\tbaltimore_class\n"
            "\tX\n"
            "Flaviviridae\tIV\n"
        )
        ann = load_family_annotations(p)
        assert list(ann.keys()) == ["flaviviridae"]

    def test_missing_family_column_returns_empty(self, tmp_path):
        p = tmp_path / "ann.tsv"
        p.write_text("name\tbaltimore_class\nFlaviviridae\tIV\n")
        assert load_family_annotations(p) == {}


class TestBaltimoreInSummaryRow:
    def test_baltimore_populated_from_annotation(self):
        row = build_summary_row(
            family="Flaviviridae", family_taxid=11050, family_lineage=[],
            seq_type="nucleotide", region="whole_genome", segment=None,
            n_species_discovered=0, n_species_with_seqs=0, seqlen_stats={},
            tree_stats={},
            family_annotation={"baltimore_class": "IV", "host_range": "Vertebrates"},
        )
        assert row["baltimore_class"] == "IV"

    def test_baltimore_empty_when_no_annotation(self):
        row = build_summary_row(
            family="X", family_taxid=None, family_lineage=[],
            seq_type="nucleotide", region="whole_genome", segment=None,
            n_species_discovered=0, n_species_with_seqs=0, seqlen_stats={},
            tree_stats={},
        )
        assert row["baltimore_class"] == ""

    def test_baltimore_column_follows_lineage(self):
        from vfam_trees.summary import COLUMNS
        assert COLUMNS.index("baltimore_class") == COLUMNS.index("lineage") + 1


class TestStatusRow:
    def test_column_order_exact(self):
        assert STATUS_COLUMNS == [
            "family", "ncbi_taxid", "molecule_region", "status",
            "lineage", "baltimore_class",
        ]

    def test_ok_status_row(self):
        row = build_status_row(
            family="Flaviviridae", family_taxid=11050,
            family_lineage=[{"name": "Viruses"}, {"name": "Riboviria"}],
            seq_type="nucleotide", region="whole_genome", segment=None,
            status="OK",
            family_annotation={"baltimore_class": "IV"},
        )
        assert row["family"] == "Flaviviridae"
        assert row["ncbi_taxid"] == 11050
        assert row["molecule_region"] == "nucleotide, whole genome"
        assert row["status"] == "OK"
        assert "Riboviria" in row["lineage"]
        assert row["baltimore_class"] == "IV"

    def test_skip_reason_row(self):
        row = build_status_row(
            family="X", family_taxid=None, family_lineage=[],
            seq_type="nucleotide", region="whole_genome", segment=None,
            status="no species found in NCBI taxonomy",
        )
        assert row["ncbi_taxid"] == ""
        assert row["status"] == "no species found in NCBI taxonomy"
        assert row["lineage"] == ""
        assert row["baltimore_class"] == ""

    def test_empty_annotation_yields_empty_baltimore(self):
        row = build_status_row(
            family="Unknown", family_taxid=1, family_lineage=[],
            seq_type="protein", region="RdRp", segment=None,
            status="OK",
        )
        assert row["baltimore_class"] == ""
        assert row["molecule_region"] == "protein, gene: RdRp"


def test_write_status_row_creates_header_then_appends(tmp_path):
    path = tmp_path / "status.tsv"
    ok = build_status_row(
        family="A", family_taxid=1, family_lineage=[],
        seq_type="nucleotide", region="whole_genome", segment=None,
        status="OK",
        family_annotation={"baltimore_class": "IV"},
    )
    skip = build_status_row(
        family="B", family_taxid=None, family_lineage=[],
        seq_type="nucleotide", region="whole_genome", segment=None,
        status="no species found in NCBI taxonomy",
    )
    write_status_row(path, ok)
    write_status_row(path, skip)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames == STATUS_COLUMNS
        rows = list(reader)
    assert len(rows) == 2
    assert rows[0]["family"] == "A"
    assert rows[0]["status"] == "OK"
    assert rows[0]["baltimore_class"] == "IV"
    assert rows[1]["status"] == "no species found in NCBI taxonomy"
    assert rows[1]["baltimore_class"] == ""


def test_write_summary_row_creates_header_then_appends(tmp_path):
    path = tmp_path / "summary.tsv"
    row_a = build_summary_row("A", 1, [], "nucleotide", "whole_genome", None,
                              0, 0, {}, {})
    row_b = build_summary_row("B", 2, [], "protein", "polymerase", None,
                              0, 0, {}, {})
    write_summary_row(path, row_a)
    write_summary_row(path, row_b)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 2
    assert rows[0]["family"] == "A"
    assert rows[1]["family"] == "B"
