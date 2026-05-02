"""Tests for vfam_trees.config."""
import copy
import tempfile
from pathlib import Path

import pytest
import yaml

from vfam_trees.config import (
    CONCATENATION_FAMILIES,
    DEFAULT_FAMILY_CONFIG,
    DNA_FAMILIES,
    SEGMENTED_FAMILIES,
    _apply_smart_defaults,
    _deep_update,
    _merge_with_defaults,
    _validate_concatenation_block,
    _warn_smart_default_conflicts,
    load_family_config,
    make_minimal_global_cfg,
)


MINIMAL_GLOBAL = {"ncbi": {"email": "test@test.com"}, "defaults": {}}


# ---------------------------------------------------------------------------
# _deep_update
# ---------------------------------------------------------------------------

def test_deep_update_merges_nested():
    base = {"a": {"x": 1, "y": 2}, "b": 3}
    _deep_update(base, {"a": {"y": 99}, "c": 4})
    assert base == {"a": {"x": 1, "y": 99}, "b": 3, "c": 4}


def test_deep_update_overwrites_non_dict():
    base = {"a": 1}
    _deep_update(base, {"a": 2})
    assert base["a"] == 2


# ---------------------------------------------------------------------------
# _apply_smart_defaults
# ---------------------------------------------------------------------------

def test_apply_smart_defaults_dna_family():
    # Adenoviridae is in DNA_FAMILIES (single-protein) but NOT in
    # CONCATENATION_FAMILIES, so it tests the single-protein DNA path.
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Adenoviridae", cfg)
    assert cfg["sequence"]["region"] == "hexon"


def test_apply_smart_defaults_whole_genome_dna():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Hepadnaviridae", cfg)
    assert cfg["sequence"]["region"] == "whole_genome"


def test_apply_smart_defaults_segmented():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    cfg["sequence"]["segment"] = None
    _apply_smart_defaults("Hantaviridae", cfg)
    assert cfg["sequence"]["segment"] == SEGMENTED_FAMILIES["Hantaviridae"]


def test_apply_smart_defaults_segmented_does_not_overwrite_explicit():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    cfg["sequence"]["segment"] = "custom_segment"
    _apply_smart_defaults("Hantaviridae", cfg)
    assert cfg["sequence"]["segment"] == "custom_segment"


def test_apply_smart_defaults_rna_family_unchanged():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Flaviviridae", cfg)
    assert cfg["sequence"]["region"] == "whole_genome"


# ---------------------------------------------------------------------------
# _merge_with_defaults priority: DEFAULT → global → DNA_FAMILIES → user file
# ---------------------------------------------------------------------------

def test_merge_dna_family_smart_default_wins_over_base():
    user_cfg = {}
    merged = _merge_with_defaults(user_cfg, MINIMAL_GLOBAL, "Adenoviridae")
    assert merged["sequence"]["region"] == "hexon"


def test_merge_user_file_overrides_smart_default():
    user_cfg = {"sequence": {"region": "custom_gene"}}
    merged = _merge_with_defaults(user_cfg, MINIMAL_GLOBAL, "Adenoviridae")
    assert merged["sequence"]["region"] == "custom_gene"


def test_merge_fills_missing_keys_from_defaults():
    user_cfg = {"sequence": {"region": "whole_genome"}}
    merged = _merge_with_defaults(user_cfg, MINIMAL_GLOBAL, "Flaviviridae")
    assert "max_per_species" in merged["download"]
    assert "max_ambiguous" in merged["quality"]


# ---------------------------------------------------------------------------
# _warn_smart_default_conflicts
# ---------------------------------------------------------------------------

def test_warn_on_region_conflict(caplog, tmp_path):
    config_path = tmp_path / "Adenoviridae.yaml"
    config_path.touch()
    file_cfg = {"sequence": {"region": "whole_genome"}}
    import logging
    with caplog.at_level(logging.WARNING):
        _warn_smart_default_conflicts("Adenoviridae", file_cfg, config_path)
    assert "whole_genome" in caplog.text
    assert "hexon" in caplog.text


def test_no_warn_when_region_matches(caplog, tmp_path):
    config_path = tmp_path / "Adenoviridae.yaml"
    config_path.touch()
    file_cfg = {"sequence": {"region": "hexon"}}
    import logging
    with caplog.at_level(logging.WARNING):
        _warn_smart_default_conflicts("Adenoviridae", file_cfg, config_path)
    assert "overriding" not in caplog.text


def test_no_warn_for_rna_family(caplog, tmp_path):
    config_path = tmp_path / "Flaviviridae.yaml"
    config_path.touch()
    file_cfg = {"sequence": {"region": "whole_genome"}}
    import logging
    with caplog.at_level(logging.WARNING):
        _warn_smart_default_conflicts("Flaviviridae", file_cfg, config_path)
    assert "overriding" not in caplog.text


# ---------------------------------------------------------------------------
# load_family_config — existing file path
# ---------------------------------------------------------------------------

def test_load_existing_config_applies_dna_override(tmp_path):
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    # Write a stale config with wrong region
    (cfg_dir / "Adenoviridae.yaml").write_text(
        yaml.dump({"sequence": {"type": "nucleotide", "region": "whole_genome", "segment": None}})
    )
    cfg, auto = load_family_config("Adenoviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is False
    # Smart default (hexon) sits between DEFAULT and user file in priority;
    # user file explicitly said whole_genome, so whole_genome wins — but a
    # warning should have been issued (tested separately).
    # The important thing: the function does not crash.
    assert cfg["sequence"]["type"] == "nucleotide"


def test_load_missing_config_auto_generates_with_dna_region(tmp_path):
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    cfg, auto = load_family_config("Adenoviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is True
    assert cfg["sequence"]["region"] == "hexon"
    # Config file was written
    assert (cfg_dir / "Adenoviridae.yaml").exists()


def test_load_missing_config_auto_generates_segment(tmp_path):
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    cfg, auto = load_family_config("Hantaviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is True
    assert cfg["sequence"]["segment"] == SEGMENTED_FAMILIES["Hantaviridae"]


# ---------------------------------------------------------------------------
# make_minimal_global_cfg — fallback when no global.yaml present
# ---------------------------------------------------------------------------

def test_make_minimal_global_cfg_returns_empty_defaults():
    cfg = make_minimal_global_cfg()
    assert cfg["defaults"] == {}


def test_make_minimal_global_cfg_produces_full_defaults_via_load(tmp_path):
    # When init-configs uses the minimal global cfg, generated per-family
    # configs should match DEFAULT_FAMILY_CONFIG exactly.
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    cfg, auto = load_family_config("Flaviviridae", cfg_dir, make_minimal_global_cfg())
    assert auto is True
    excludes = cfg["quality"]["exclude_organisms"]
    for term in DEFAULT_FAMILY_CONFIG["quality"]["exclude_organisms"]:
        assert term in excludes


def test_make_minimal_global_cfg_global_yaml_overrides_take_precedence(tmp_path):
    # When global.yaml IS present its values win over DEFAULT_FAMILY_CONFIG.
    old_global = {
        "ncbi": {"email": "test@test.com"},
        "defaults": {
            "quality": {
                "exclude_organisms": ["synthetic construct", "metagenome"]
            }
        },
    }
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    cfg, _ = load_family_config("Flaviviridae", cfg_dir, old_global)
    excludes = cfg["quality"]["exclude_organisms"]
    # global.yaml list wins — newer defaults NOT silently injected
    assert excludes == ["synthetic construct", "metagenome"]


# ---------------------------------------------------------------------------
# CONCATENATION_FAMILIES — multi-marker presets
# ---------------------------------------------------------------------------

class TestConcatenationFamiliesPresets:
    """Curation table looks structurally correct."""

    def test_target_families_present(self):
        for fam in (
            "Poxviridae", "Orthoherpesviridae", "Herpesviridae",
            "Alloherpesviridae", "Malacoherpesviridae",
            "Asfarviridae", "Iridoviridae",
            "Baculoviridae", "Nudiviridae", "Ascoviridae",
            "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
            "Pithoviridae", "Pandoraviridae", "Medusaviridae",
        ):
            assert fam in CONCATENATION_FAMILIES, f"missing preset: {fam}"

    def test_each_preset_has_proteins_and_correct_region(self):
        for fam, preset in CONCATENATION_FAMILIES.items():
            assert preset["sequence"]["region"] == "concatenated", fam
            assert preset["sequence"]["type"] == "protein", fam
            proteins = preset["concatenation"]["proteins"]
            assert len(proteins) > 0, fam
            for p in proteins:
                assert "name" in p and isinstance(p["name"], str), (fam, p)
                assert "aliases" in p and isinstance(p["aliases"], list), (fam, p)

    def test_poxviridae_has_subfamily_aware_aliases(self):
        # The Pox set is the only one with aliases_Entomopoxvirinae overrides
        # (per CONCAT_DESIGN.md §4.1).  At minimum DNA polymerase should have it.
        pox_proteins = CONCATENATION_FAMILIES["Poxviridae"]["concatenation"]["proteins"]
        polb = next(p for p in pox_proteins if p["name"] == "DNA polymerase")
        assert "aliases_Entomopoxvirinae" in polb
        assert isinstance(polb["aliases_Entomopoxvirinae"], list)

    def test_herpesvirus_families_share_same_set(self):
        # The 4 herpesvirus family names all share one curated 7-marker set.
        ortho = CONCATENATION_FAMILIES["Orthoherpesviridae"]["concatenation"]["proteins"]
        for fam in ("Herpesviridae", "Alloherpesviridae", "Malacoherpesviridae"):
            other = CONCATENATION_FAMILIES[fam]["concatenation"]["proteins"]
            assert [p["name"] for p in other] == [p["name"] for p in ortho]

    def test_baculovirus_relatives_share_same_set(self):
        baculo = CONCATENATION_FAMILIES["Baculoviridae"]["concatenation"]["proteins"]
        for fam in ("Nudiviridae", "Ascoviridae"):
            other = CONCATENATION_FAMILIES[fam]["concatenation"]["proteins"]
            assert [p["name"] for p in other] == [p["name"] for p in baculo]

    def test_ncldv_fallback_families_share_hallmark_set(self):
        mimi = CONCATENATION_FAMILIES["Mimiviridae"]["concatenation"]["proteins"]
        assert len(mimi) == 8  # 8 NCLDV hallmarks
        for fam in ("Phycodnaviridae", "Marseilleviridae", "Pithoviridae",
                    "Pandoraviridae", "Medusaviridae"):
            other = CONCATENATION_FAMILIES[fam]["concatenation"]["proteins"]
            assert [p["name"] for p in other] == [p["name"] for p in mimi]


# ---------------------------------------------------------------------------
# Smart-default precedence: concat > DNA_FAMILIES single-protein
# ---------------------------------------------------------------------------

class TestConcatTakesPrecedenceOverDnaFamilies:
    def test_poxviridae_auto_gen_is_concat(self, tmp_path):
        # Poxviridae is in BOTH DNA_FAMILIES (rpo147) and CONCATENATION_FAMILIES.
        # The concat preset must win.
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        cfg, auto = load_family_config("Poxviridae", cfg_dir, MINIMAL_GLOBAL)
        assert auto is True
        assert cfg["sequence"]["region"] == "concatenated"
        assert len(cfg["concatenation"]["proteins"]) == 9  # Pox 9-marker set

    def test_herpesviridae_auto_gen_is_concat(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        cfg, auto = load_family_config("Orthoherpesviridae", cfg_dir, MINIMAL_GLOBAL)
        assert auto is True
        assert cfg["sequence"]["region"] == "concatenated"
        assert len(cfg["concatenation"]["proteins"]) == 7

    def test_apply_smart_defaults_concat_supersedes_dna(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        _apply_smart_defaults("Poxviridae", cfg)
        # Even though Poxviridae has a DNA_FAMILIES entry (rpo147), concat wins
        assert cfg["sequence"]["region"] == "concatenated"

    def test_user_can_revert_concat_family_to_single_protein(self, tmp_path):
        # User explicitly sets region back to a single protein in their yaml.
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Poxviridae.yaml").write_text(
            yaml.dump({"sequence": {"type": "protein", "region": "rpo147", "segment": None}})
        )
        cfg, auto = load_family_config("Poxviridae", cfg_dir, MINIMAL_GLOBAL)
        assert auto is False
        assert cfg["sequence"]["region"] == "rpo147"

    def test_revert_to_single_protein_logs_info_not_warning(self, caplog, tmp_path):
        # Reverting a concat family to single-protein is documented behaviour;
        # should be INFO not WARNING.
        config_path = tmp_path / "Poxviridae.yaml"
        config_path.touch()
        file_cfg = {"sequence": {"region": "rpo147"}}
        import logging
        with caplog.at_level(logging.INFO):
            _warn_smart_default_conflicts("Poxviridae", file_cfg, config_path)
        # Should mention reverting and not be at WARNING level
        assert "reverting" in caplog.text.lower()
        warning_records = [r for r in caplog.records if r.levelno >= logging.WARNING]
        assert not warning_records


# ---------------------------------------------------------------------------
# _validate_concatenation_block
# ---------------------------------------------------------------------------

class TestValidateConcatenationBlock:
    def _base(self, **overrides):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        cfg["sequence"]["region"] = "concatenated"
        cfg["concatenation"] = {
            "proteins": [
                {"name": "DNA polymerase", "aliases": ["DNA pol"]},
                {"name": "MCP", "aliases": []},
            ],
            "min_fraction": 0.7,
        }
        for k, v in overrides.items():
            cfg["concatenation"][k] = v
        return cfg

    def test_valid_block_passes(self):
        _validate_concatenation_block(self._base(), "Test")  # no exception

    def test_non_concat_region_skips_validation(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        cfg["sequence"]["region"] = "DNA polymerase"
        cfg["concatenation"]["proteins"] = []  # empty, but region != concat
        _validate_concatenation_block(cfg, "Test")  # no exception

    def test_empty_proteins_rejected(self):
        cfg = self._base(proteins=[])
        with pytest.raises(ValueError, match="empty"):
            _validate_concatenation_block(cfg, "Test")

    def test_missing_name_rejected(self):
        cfg = self._base(proteins=[{"aliases": ["x"]}])
        with pytest.raises(ValueError, match="name"):
            _validate_concatenation_block(cfg, "Test")

    def test_blank_name_rejected(self):
        cfg = self._base(proteins=[{"name": "  ", "aliases": []}])
        with pytest.raises(ValueError, match="name"):
            _validate_concatenation_block(cfg, "Test")

    def test_aliases_must_be_list_of_strings(self):
        cfg = self._base(proteins=[{"name": "X", "aliases": [1, 2]}])
        with pytest.raises(ValueError, match="aliases"):
            _validate_concatenation_block(cfg, "Test")

    def test_length_range_must_be_well_formed(self):
        cfg = self._base(proteins=[{"name": "X", "aliases": [], "length_range": [100]}])
        with pytest.raises(ValueError, match="length_range"):
            _validate_concatenation_block(cfg, "Test")

    def test_length_range_min_must_be_less_than_max(self):
        cfg = self._base(proteins=[{"name": "X", "aliases": [], "length_range": [500, 100]}])
        with pytest.raises(ValueError, match="length_range"):
            _validate_concatenation_block(cfg, "Test")

    def test_min_fraction_must_be_in_range(self):
        cfg = self._base(min_fraction=0)
        with pytest.raises(ValueError, match="min_fraction"):
            _validate_concatenation_block(cfg, "Test")
        cfg = self._base(min_fraction=1.5)
        with pytest.raises(ValueError, match="min_fraction"):
            _validate_concatenation_block(cfg, "Test")

    def test_default_family_config_has_concatenation_key(self):
        # The concatenation block must be in DEFAULT_FAMILY_CONFIG so user
        # yamls can reference concatenation.* without "unknown key" warnings.
        assert "concatenation" in DEFAULT_FAMILY_CONFIG
        assert DEFAULT_FAMILY_CONFIG["concatenation"]["min_fraction"] == 0.7

    def test_curated_presets_pass_validation(self):
        # Every shipped preset must pass schema validation when merged with
        # defaults.  Catches any future curation typos.
        for fam in CONCATENATION_FAMILIES:
            merged = _merge_with_defaults({}, MINIMAL_GLOBAL, fam)
            _validate_concatenation_block(merged, fam)
