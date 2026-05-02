"""Tests for vfam_trees.markers — NameMatchIdentifier (v1)."""
from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from vfam_trees.markers import NameMatchIdentifier


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rec(acc: str, description: str, length: int = 500,
         locus_tag: str | None = None, product: str | None = None) -> SeqRecord:
    """Build a minimal protein-like SeqRecord for matcher input."""
    rec = SeqRecord(Seq("M" * length), id=acc, description=description)
    rec.annotations = {"organism": "Test virus"}
    if locus_tag:
        rec.annotations["locus_tag"] = locus_tag
    if product:
        feat = SeqFeature(FeatureLocation(0, length), type="Protein")
        feat.qualifiers = {"product": [product]}
        rec.features = [feat]
    else:
        rec.features = []
    return rec


# Common fixture: a tiny marker set
_POLB = {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase", "DNA pol"]}
_MCP  = {"name": "major capsid protein", "aliases": ["MCP"]}


# ---------------------------------------------------------------------------
# Basic name + alias matching
# ---------------------------------------------------------------------------

class TestBasicMatching:
    def setup_method(self):
        self.id = NameMatchIdentifier()

    def test_name_match(self):
        proteins = [_rec("YP_1", "DNA polymerase [Vaccinia virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert "DNA polymerase" in out
        assert out["DNA polymerase"].id == "YP_1"

    def test_alias_match(self):
        proteins = [_rec("YP_1", "DNA pol [Vaccinia virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_case_insensitive(self):
        proteins = [_rec("YP_1", "DNA POLYMERASE [some virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_missing_marker_omitted_from_result(self):
        proteins = [_rec("YP_1", "DNA polymerase")]
        out = self.id.identify(proteins, [_POLB, _MCP])
        assert "DNA polymerase" in out
        assert "major capsid protein" not in out

    def test_empty_proteins_returns_empty_dict(self):
        out = self.id.identify([], [_POLB, _MCP])
        assert out == {}

    def test_empty_marker_set_returns_empty_dict(self):
        out = self.id.identify([_rec("YP_1", "DNA polymerase")], [])
        assert out == {}

    def test_multiple_markers_all_found(self):
        proteins = [
            _rec("YP_1", "DNA polymerase"),
            _rec("YP_2", "major capsid protein"),
        ]
        out = self.id.identify(proteins, [_POLB, _MCP])
        assert out["DNA polymerase"].id == "YP_1"
        assert out["major capsid protein"].id == "YP_2"

    def test_match_against_protein_feature_product(self):
        # Description is generic; product qualifier carries the marker name.
        proteins = [_rec("YP_1", "hypothetical protein", product="DNA polymerase")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"


# ---------------------------------------------------------------------------
# Subfamily-aware aliases
# ---------------------------------------------------------------------------

class TestSubfamilyAwareAliases:
    def setup_method(self):
        self.id = NameMatchIdentifier()
        self.marker = {
            "name": "DNA polymerase",
            "aliases": ["DNA pol"],
            "aliases_Entomopoxvirinae": ["DNA polymerase B"],
        }

    def test_subfamily_alias_used_when_lineage_matches(self):
        proteins = [_rec("YP_1", "DNA polymerase B [Amsacta entomopoxvirus]")]
        lineage = [
            {"rank": "family", "name": "Poxviridae"},
            {"rank": "subfamily", "name": "Entomopoxvirinae"},
            {"rank": "genus", "name": "Alphaentomopoxvirus"},
        ]
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert out["DNA polymerase"].id == "YP_1"

    def test_subfamily_alias_ignored_when_lineage_is_other_subfamily(self):
        # Chordopoxvirinae lineage → entomopox-specific alias not added.
        # The protein only matches "DNA polymerase B", which is entomopox-specific.
        # Note: "DNA polymerase B" still contains "DNA pol" substring — so it
        # would still match via the base alias.  Use a description that does
        # NOT contain any base alias to test true exclusion.
        proteins = [_rec("YP_1", "polymerase variant Z")]  # no base alias hit
        lineage = [
            {"rank": "subfamily", "name": "Chordopoxvirinae"},
        ]
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert "DNA polymerase" not in out

    def test_subfamily_alias_ignored_when_no_lineage_provided(self):
        # No lineage → only base aliases used.
        proteins = [_rec("YP_1", "polymerase variant Z")]
        out = self.id.identify(proteins, [self.marker], species_lineage=None)
        assert "DNA polymerase" not in out

    def test_subfamily_alias_ignored_when_lineage_has_no_subfamily_rank(self):
        proteins = [_rec("YP_1", "polymerase variant Z")]
        lineage = [{"rank": "family", "name": "Poxviridae"}]  # no subfamily entry
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert "DNA polymerase" not in out


# ---------------------------------------------------------------------------
# Tiebreakers
# ---------------------------------------------------------------------------

class TestTiebreakers:
    def setup_method(self):
        self.id = NameMatchIdentifier()

    def test_locus_tag_hint_wins(self):
        # Two candidates both match the name; only one has the matching locus_tag.
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        proteins = [
            _rec("YP_1", "DNA polymerase", locus_tag="A24R"),
            _rec("YP_2", "DNA polymerase", locus_tag="E9L"),
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_locus_tag_hint_regex_alternation(self):
        marker = {**_POLB, "locus_tag_hint": r"E9L|polB"}
        proteins = [
            _rec("YP_1", "DNA polymerase", locus_tag="random"),
            _rec("YP_2", "DNA polymerase", locus_tag="some-polB-thing"),
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_locus_tag_hint_no_match_falls_through(self):
        # No candidate has the hint locus tag; fall through to next tiebreaker.
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        proteins = [
            _rec("YP_1", "DNA polymerase", length=900),
            _rec("YP_2", "DNA polymerase", length=1100),  # longest wins
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_length_range_midpoint(self):
        marker = {**_POLB, "length_range": [800, 1200]}  # midpoint = 1000
        proteins = [
            _rec("YP_1", "DNA polymerase", length=950),   # |950-1000| = 50
            _rec("YP_2", "DNA polymerase", length=1100),  # |1100-1000| = 100
            _rec("YP_3", "DNA polymerase", length=995),   # |995-1000| = 5  ← winner
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_3"

    def test_longest_wins_when_no_hints(self):
        proteins = [
            _rec("YP_1", "DNA polymerase", length=500),
            _rec("YP_2", "DNA polymerase", length=1000),
            _rec("YP_3", "DNA polymerase", length=800),
        ]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_2"

    def test_lowest_accession_breaks_length_ties(self):
        proteins = [
            _rec("YP_3", "DNA polymerase", length=900),
            _rec("YP_1", "DNA polymerase", length=900),
            _rec("YP_2", "DNA polymerase", length=900),
        ]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_tiebreaker_priority_locus_tag_then_length(self):
        # Locus-tag hint defined AND length range defined.  Locus-tag must win.
        marker = {
            **_POLB,
            "locus_tag_hint": r"E9L",
            "length_range": [900, 1100],   # midpoint = 1000
        }
        proteins = [
            # Hint match but length far from midpoint (700)
            _rec("YP_1", "DNA polymerase", length=700, locus_tag="E9L"),
            # No hint match but length right at midpoint (1000)
            _rec("YP_2", "DNA polymerase", length=1000, locus_tag="random"),
        ]
        out = self.id.identify(proteins, [marker])
        # Locus-tag hint is priority 1; YP_1 wins despite being further from midpoint.
        assert out["DNA polymerase"].id == "YP_1"

    def test_locus_tag_via_feature_qualifier(self):
        # When annotations don't carry locus_tag, fall back to feature qualifier.
        rec = _rec("YP_1", "DNA polymerase")
        feat = SeqFeature(FeatureLocation(0, 500), type="CDS")
        feat.qualifiers = {"locus_tag": ["E9L"]}
        rec.features = [feat]
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        out = self.id.identify([rec], [marker])
        assert out["DNA polymerase"].id == "YP_1"
