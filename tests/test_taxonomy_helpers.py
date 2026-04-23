"""Tests for taxonomy internal helpers: LCA, rank inference, label thinning, bootstrap normalization."""
from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees.taxonomy import (
    _find_lca_lineage,
    _infer_rank,
    _keep_deepest_labels,
    _normalize_bootstrap,
)


# ---------------------------------------------------------------------------
# _find_lca_lineage
# ---------------------------------------------------------------------------

class TestFindLcaLineage:
    def test_empty_input(self):
        assert _find_lca_lineage([]) == []

    def test_single_lineage_returned_unchanged(self):
        lin = [["Viruses", "Riboviria", "Flaviviridae", "Flavivirus"]]
        assert _find_lca_lineage(lin) == lin[0]

    def test_identical_lineages_full_match(self):
        l = ["Viruses", "Riboviria", "Flaviviridae"]
        assert _find_lca_lineage([l, l]) == l

    def test_common_prefix_to_family(self):
        lineages = [
            ["Viruses", "Riboviria", "Flaviviridae", "Flavivirus", "Dengue virus"],
            ["Viruses", "Riboviria", "Flaviviridae", "Pegivirus"],
        ]
        assert _find_lca_lineage(lineages) == ["Viruses", "Riboviria", "Flaviviridae"]

    def test_no_common_prefix(self):
        lineages = [["A", "B"], ["X", "Y"]]
        assert _find_lca_lineage(lineages) == []

    def test_common_prefix_stops_at_first_divergence(self):
        lineages = [
            ["A", "B", "C", "D"],
            ["A", "B", "Q"],
            ["A", "B", "C", "Z"],
        ]
        assert _find_lca_lineage(lineages) == ["A", "B"]


# ---------------------------------------------------------------------------
# _infer_rank (ICTV suffix-based rank inference)
# ---------------------------------------------------------------------------

class TestInferRank:
    def test_family_suffix_viridae(self):
        assert _infer_rank("Flaviviridae") == "family"

    def test_subfamily_suffix_virinae(self):
        assert _infer_rank("Orthocoronavirinae") == "subfamily"

    def test_order_suffix_virales(self):
        assert _infer_rank("Nidovirales") == "order"

    def test_class_suffix_viricetes(self):
        assert _infer_rank("Pisoniviricetes") == "class"

    def test_phylum_suffix_viricota(self):
        assert _infer_rank("Pisuviricota") == "phylum"

    def test_realm_suffix_viria(self):
        assert _infer_rank("Riboviria") == "realm"

    def test_kingdom_suffix_virae(self):
        assert _infer_rank("Orthornavirae") == "kingdom"

    def test_single_word_virus_is_genus(self):
        assert _infer_rank("Flavivirus") == "genus"

    def test_multi_word_virus_is_unranked(self):
        assert _infer_rank("Dengue virus") == ""

    def test_empty_name_returns_empty(self):
        assert _infer_rank("") == ""

    def test_longest_suffix_wins_viricetes_over_virae(self):
        # "Pisoniviricetes" ends with both "viricetes" and "virae"
        # The longest-first ordering must resolve to "class".
        assert _infer_rank("Pisoniviricetes") == "class"


# ---------------------------------------------------------------------------
# _keep_deepest_labels
# ---------------------------------------------------------------------------

class TestKeepDeepestLabels:
    def _tree_with_duplicate_labels(self):
        # Structure:
        #   root(Flaviviridae)
        #     ├── inner_a(Flavivirus)
        #     │     ├── leaf_1
        #     │     └── leaf_2
        #     └── inner_b(Flavivirus)  <-- duplicate of inner_a
        #           ├── leaf_3
        #           └── leaf_4
        leaves = [Clade(name=f"l{i}") for i in range(1, 5)]
        inner_a = Clade(clades=leaves[:2], name="Flavivirus")
        inner_b = Clade(clades=leaves[2:], name="Flavivirus")
        root = Clade(clades=[inner_a, inner_b], name="Flaviviridae")
        return Tree(root=root), inner_a, inner_b

    def test_one_crown_label_kept_duplicate_cleared(self):
        tree, a, b = self._tree_with_duplicate_labels()
        _keep_deepest_labels(tree)
        # Exactly one of the two "Flavivirus" nodes should retain its name.
        names = [a.name, b.name]
        assert names.count("Flavivirus") == 1
        assert names.count("") == 1

    def test_crown_is_larger_subtree(self):
        # Make inner_a the bigger subtree so it wins the tie.
        leaves = [Clade(name=f"l{i}") for i in range(1, 6)]
        inner_a = Clade(clades=leaves[:3], name="Dup")
        inner_b = Clade(clades=leaves[3:], name="Dup")
        root = Clade(clades=[inner_a, inner_b])
        tree = Tree(root=root)
        _keep_deepest_labels(tree)
        # Larger clade retains its label; smaller one is cleared
        assert inner_a.name == "Dup"
        assert inner_b.name == ""

    def test_distinct_labels_both_kept(self):
        leaves = [Clade(name=f"l{i}") for i in range(4)]
        inner_a = Clade(clades=leaves[:2], name="Alpha")
        inner_b = Clade(clades=leaves[2:], name="Beta")
        root = Clade(clades=[inner_a, inner_b])
        _keep_deepest_labels(Tree(root=root))
        assert inner_a.name == "Alpha"
        assert inner_b.name == "Beta"


# ---------------------------------------------------------------------------
# _normalize_bootstrap
# ---------------------------------------------------------------------------

class TestNormalizeBootstrap:
    def _tree_with_confidences(self, values):
        leaves = [Clade(name=f"l{i}") for i in range(len(values) + 1)]
        clades = [Clade(clades=[leaves[i], leaves[i+1]], confidence=v)
                  for i, v in enumerate(values)]
        # Chain them into a comb tree
        def chain(cs):
            if len(cs) == 1:
                return cs[0]
            return Clade(clades=[cs[0], chain(cs[1:])])
        return Tree(root=chain(clades))

    def test_zero_to_one_scaled_to_0_100(self):
        t = self._tree_with_confidences([0.85, 0.95, 0.1])
        _normalize_bootstrap(t)
        confs = sorted(c.confidence for c in t.find_clades() if c.confidence is not None)
        assert confs == [10, 85, 95]

    def test_already_0_100_rounded(self):
        t = self._tree_with_confidences([85.4, 95.6, 33.0])
        _normalize_bootstrap(t)
        confs = sorted(c.confidence for c in t.find_clades() if c.confidence is not None)
        assert confs == [33, 85, 96]

    def test_no_confidences_no_op(self):
        leaves = [Clade(name="a"), Clade(name="b")]
        t = Tree(root=Clade(clades=leaves))
        _normalize_bootstrap(t)   # must not raise
        assert all(c.confidence is None for c in t.find_clades())

    def test_max_exactly_1_treated_as_0_1_scale(self):
        # Edge case: the normalization threshold is max ≤ 1.0
        t = self._tree_with_confidences([0.5, 1.0])
        _normalize_bootstrap(t)
        confs = sorted(c.confidence for c in t.find_clades() if c.confidence is not None)
        assert confs == [50, 100]
