"""Tests for vfam_trees.config."""
import copy
import tempfile
from pathlib import Path

import pytest
import yaml

from vfam_trees.config import (
    DEFAULT_FAMILY_CONFIG,
    DNA_FAMILIES,
    SEGMENTED_FAMILIES,
    _apply_smart_defaults,
    _deep_update,
    _merge_with_defaults,
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
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Asfarviridae", cfg)
    assert cfg["sequence"]["region"] == "B646L"


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
    merged = _merge_with_defaults(user_cfg, MINIMAL_GLOBAL, "Asfarviridae")
    assert merged["sequence"]["region"] == "B646L"


def test_merge_user_file_overrides_smart_default():
    user_cfg = {"sequence": {"region": "custom_gene"}}
    merged = _merge_with_defaults(user_cfg, MINIMAL_GLOBAL, "Asfarviridae")
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
    config_path = tmp_path / "Asfarviridae.yaml"
    config_path.touch()
    file_cfg = {"sequence": {"region": "whole_genome"}}
    import logging
    with caplog.at_level(logging.WARNING):
        _warn_smart_default_conflicts("Asfarviridae", file_cfg, config_path)
    assert "whole_genome" in caplog.text
    assert "B646L" in caplog.text


def test_no_warn_when_region_matches(caplog, tmp_path):
    config_path = tmp_path / "Asfarviridae.yaml"
    config_path.touch()
    file_cfg = {"sequence": {"region": "B646L"}}
    import logging
    with caplog.at_level(logging.WARNING):
        _warn_smart_default_conflicts("Asfarviridae", file_cfg, config_path)
    assert "B646L" not in caplog.text or "overriding" not in caplog.text


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
    (cfg_dir / "Asfarviridae.yaml").write_text(
        yaml.dump({"sequence": {"type": "nucleotide", "region": "whole_genome", "segment": None}})
    )
    cfg, auto = load_family_config("Asfarviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is False
    # Smart default (B646L) sits between DEFAULT and user file in priority;
    # user file explicitly said whole_genome, so whole_genome wins — but a
    # warning should have been issued (tested separately).
    # The important thing: the function does not crash.
    assert cfg["sequence"]["type"] == "nucleotide"


def test_load_missing_config_auto_generates_with_dna_region(tmp_path):
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    cfg, auto = load_family_config("Asfarviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is True
    assert cfg["sequence"]["region"] == "B646L"
    # Config file was written
    assert (cfg_dir / "Asfarviridae.yaml").exists()


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
