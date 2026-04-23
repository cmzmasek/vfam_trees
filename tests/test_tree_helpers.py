"""Tests for vfam_trees.tree helpers (no external tool calls)."""
import pytest

from vfam_trees.tree import (
    _fasttree_aa_flags,
    _fasttree_nuc_flags,
    is_model_finder_spec,
    parse_iqtree_best_model,
    validate_newick,
)


# ---------------------------------------------------------------------------
# is_model_finder_spec
# ---------------------------------------------------------------------------

class TestIsModelFinderSpec:
    def test_empty_is_model_finder(self):
        assert is_model_finder_spec("") is True

    def test_mf_is_model_finder(self):
        assert is_model_finder_spec("MF") is True

    def test_mfp_is_model_finder(self):
        assert is_model_finder_spec("MFP") is True

    def test_test_is_model_finder(self):
        assert is_model_finder_spec("TEST") is True
        assert is_model_finder_spec("TESTONLY") is True

    def test_concrete_model_not_model_finder(self):
        assert is_model_finder_spec("GTR+G") is False
        assert is_model_finder_spec("LG+G") is False
        assert is_model_finder_spec("JC69") is False

    def test_case_insensitive(self):
        assert is_model_finder_spec("mf") is True
        assert is_model_finder_spec("mfp") is True
        assert is_model_finder_spec("test") is True

    def test_model_with_rate_suffix(self):
        assert is_model_finder_spec("MFP+G") is True
        assert is_model_finder_spec("GTR+I+G") is False


# ---------------------------------------------------------------------------
# parse_iqtree_best_model
# ---------------------------------------------------------------------------

class TestParseIqtreeBestModel:
    def test_missing_file_returns_none(self, tmp_path):
        assert parse_iqtree_best_model(tmp_path / "nope.iqtree") is None

    def test_best_fit_line_preferred(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text(
            "Some header\n"
            "Best-fit model according to BIC: GTR+I+G4\n"
            "Other stuff\n"
            "Model of substitution: IGNORED\n"
        )
        assert parse_iqtree_best_model(log) == "GTR+I+G4"

    def test_fallback_to_model_of_substitution(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text("Model of substitution: LG+F+G\n")
        assert parse_iqtree_best_model(log) == "LG+F+G"

    def test_no_match_returns_none(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text("some unrelated log content\n")
        assert parse_iqtree_best_model(log) is None


# ---------------------------------------------------------------------------
# _fasttree_nuc_flags / _fasttree_aa_flags
# ---------------------------------------------------------------------------

class TestFasttreeNucFlags:
    def test_gtr_adds_gtr_flag(self):
        assert _fasttree_nuc_flags("GTR") == ["-gtr"]

    def test_gtr_plus_gamma_adds_both(self):
        flags = _fasttree_nuc_flags("GTR+G")
        assert "-gtr" in flags and "-gamma" in flags

    def test_jc_gives_default_no_flags(self):
        assert _fasttree_nuc_flags("JC") == []
        assert _fasttree_nuc_flags("JC69") == []

    def test_jc_plus_gamma(self):
        flags = _fasttree_nuc_flags("JC+G")
        assert "-gamma" in flags
        assert "-gtr" not in flags

    def test_unsupported_falls_back_to_gtr(self):
        # HKY is not supported → warn + use GTR
        flags = _fasttree_nuc_flags("HKY")
        assert "-gtr" in flags

    def test_empty_model_defaults_to_gtr(self):
        assert "-gtr" in _fasttree_nuc_flags("")


class TestFasttreeAaFlags:
    def test_lg_adds_lg(self):
        assert _fasttree_aa_flags("LG") == ["-lg"]

    def test_lg_plus_gamma(self):
        flags = _fasttree_aa_flags("LG+G")
        assert "-lg" in flags and "-gamma" in flags

    def test_wag_adds_wag(self):
        assert _fasttree_aa_flags("WAG") == ["-wag"]

    def test_jtt_is_default_no_flag(self):
        assert _fasttree_aa_flags("JTT") == []

    def test_empty_model_defaults_to_jtt(self):
        assert _fasttree_aa_flags("") == []

    def test_unsupported_falls_back_to_wag(self):
        flags = _fasttree_aa_flags("DAYHOFF")
        assert "-wag" in flags


# ---------------------------------------------------------------------------
# validate_newick
# ---------------------------------------------------------------------------

class TestValidateNewick:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            validate_newick(tmp_path / "missing.nwk")

    def test_empty_file_raises(self, tmp_path):
        empty = tmp_path / "e.nwk"
        empty.write_text("")
        with pytest.raises(FileNotFoundError):
            validate_newick(empty)

    def test_invalid_newick_raises(self, tmp_path):
        # Unbalanced parentheses make the Newick parser raise.
        bad = tmp_path / "bad.nwk"
        bad.write_text("(A:0.1,B:0.2,(C:0.3")
        with pytest.raises(ValueError):
            validate_newick(bad)

    def test_valid_newick_passes(self, tmp_path):
        good = tmp_path / "good.nwk"
        good.write_text("(A:0.1,B:0.2,(C:0.3,D:0.4):0.05);")
        validate_newick(good)  # should not raise
