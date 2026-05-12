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
    _validate_manual_block,
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
# _apply_smart_defaults — log suppression when file_cfg overrides a smart default
# ---------------------------------------------------------------------------

def test_auto_configured_segment_logged_when_file_is_silent(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Hantaviridae", cfg, file_cfg={})
    assert any("Auto-configured segment" in r.getMessage() for r in caplog.records)


def test_auto_configured_segment_suppressed_when_file_overrides_it(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    file_cfg = {"sequence": {"segment": "segment S"}}
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Hantaviridae", cfg, file_cfg=file_cfg)
    assert not any("Auto-configured segment" in r.getMessage() for r in caplog.records)


def test_auto_configured_dna_marker_logged_when_file_is_silent(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Adenoviridae", cfg, file_cfg={})
    assert any("Auto-configured DNA family" in r.getMessage() for r in caplog.records)


def test_auto_configured_dna_marker_suppressed_when_file_overrides_region(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    file_cfg = {"sequence": {"region": "penton"}}
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Adenoviridae", cfg, file_cfg=file_cfg)
    assert not any("Auto-configured DNA family" in r.getMessage() for r in caplog.records)


def test_auto_configured_concat_logged_when_file_is_silent(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Poxviridae", cfg, file_cfg={})
    assert any("Auto-configured concatenation mode" in r.getMessage() for r in caplog.records)


def test_auto_configured_concat_suppressed_when_file_reverts_to_single_protein(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    file_cfg = {"sequence": {"region": "rpo147"}}
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Poxviridae", cfg, file_cfg=file_cfg)
    assert not any("Auto-configured concatenation mode" in r.getMessage() for r in caplog.records)


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
        # Malacoherpesviridae was dropped from CONCATENATION_FAMILIES in the
        # 2026-05 audit (5 species, only MCP and DNA pol annotated — concat
        # is meaningless at that scale; falls back to single-marker MCP).
        for fam in (
            "Poxviridae", "Orthoherpesviridae", "Herpesviridae",
            "Alloherpesviridae",
            "Asfarviridae", "Iridoviridae",
            "Baculoviridae", "Nudiviridae", "Ascoviridae",
            "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
            "Pithoviridae", "Pandoraviridae", "Medusaviridae",
        ):
            assert fam in CONCATENATION_FAMILIES, f"missing preset: {fam}"
        assert "Malacoherpesviridae" not in CONCATENATION_FAMILIES

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

    def test_herpesviridae_alias_matches_ortho(self):
        # Herpesviridae (legacy ICTV family name) keeps the full 7-marker
        # Orthoherpesviridae core verbatim.  Alloherpesviridae uses a
        # 6-marker subset (drops ssDNA-binding) per the 2026-05 audit.
        ortho = CONCATENATION_FAMILIES["Orthoherpesviridae"]["concatenation"]["proteins"]
        legacy = CONCATENATION_FAMILIES["Herpesviridae"]["concatenation"]["proteins"]
        assert [p["name"] for p in legacy] == [p["name"] for p in ortho]
        assert len(ortho) == 7

    def test_alloherpesviridae_drops_ssdna_binding(self):
        ortho_names = {
            p["name"] for p in
            CONCATENATION_FAMILIES["Orthoherpesviridae"]["concatenation"]["proteins"]
        }
        allo_names = {
            p["name"] for p in
            CONCATENATION_FAMILIES["Alloherpesviridae"]["concatenation"]["proteins"]
        }
        # 6 markers: the Ortho 7 minus ssDNA-binding.
        assert len(allo_names) == 6
        assert ortho_names - allo_names == {"single-stranded DNA-binding protein"}

    def test_baculo_and_nudi_share_set_ascoviridae_does_not(self):
        baculo = CONCATENATION_FAMILIES["Baculoviridae"]["concatenation"]["proteins"]
        nudi = CONCATENATION_FAMILIES["Nudiviridae"]["concatenation"]["proteins"]
        assert [p["name"] for p in nudi] == [p["name"] for p in baculo]
        # Ascoviridae has its own reduced 3-marker preset (post 2026-05 audit).
        asco = CONCATENATION_FAMILIES["Ascoviridae"]["concatenation"]["proteins"]
        assert len(asco) == 3
        assert {p["name"] for p in asco} == {
            "DNA polymerase", "DNA helicase P143", "major capsid protein",
        }

    def test_ncldv_fallback_only_pandora_and_medusa(self):
        # After the 2026-05 audit only Pandoraviridae and Medusaviridae
        # remain on the 8-marker NCLDV-hallmark fallback; Mimi, Phyco,
        # Marseille, and Pitho got pruned per-family presets.
        pandora = CONCATENATION_FAMILIES["Pandoraviridae"]["concatenation"]["proteins"]
        medusa = CONCATENATION_FAMILIES["Medusaviridae"]["concatenation"]["proteins"]
        assert len(pandora) == 8
        assert [p["name"] for p in medusa] == [p["name"] for p in pandora]
        # Pruned-preset families have fewer markers.
        for fam, expected_n in (
            ("Phycodnaviridae", 4),
            ("Mimiviridae", 6),
            ("Marseilleviridae", 6),
            ("Pithoviridae", 6),
        ):
            n = len(CONCATENATION_FAMILIES[fam]["concatenation"]["proteins"])
            assert n == expected_n, f"{fam}: expected {expected_n} markers, got {n}"


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
        assert len(cfg["concatenation"]["proteins"]) == 8  # Pox 8-marker set (post 2026-05 audit, ssDNA-binding dropped)

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

    def test_source_nuc_min_length_frac_must_be_in_range(self):
        cfg = self._base(source_nuc_min_length_frac=-0.1)
        with pytest.raises(ValueError, match="source_nuc_min_length_frac"):
            _validate_concatenation_block(cfg, "Test")
        cfg = self._base(source_nuc_min_length_frac=1.0)
        with pytest.raises(ValueError, match="source_nuc_min_length_frac"):
            _validate_concatenation_block(cfg, "Test")

    def test_source_nuc_min_length_frac_zero_disables(self):
        cfg = self._base(source_nuc_min_length_frac=0)
        _validate_concatenation_block(cfg, "Test")  # no exception

    def test_default_family_config_has_source_nuc_min_length_frac(self):
        assert DEFAULT_FAMILY_CONFIG["concatenation"]["source_nuc_min_length_frac"] == 0.3

    def test_smart_default_concat_bumps_max_per_species(self):
        # Concat-mode families auto-get max_per_species=3000 because each
        # marker is fetched independently and partial single-gene submissions
        # would otherwise crowd out genome-scale RefSeq proteins under the
        # 300 default.
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        # Pick any family in CONCATENATION_FAMILIES — we only care about the
        # smart-default behaviour, not the marker preset.
        sample = next(iter(CONCATENATION_FAMILIES))
        _apply_smart_defaults(sample, cfg)
        assert cfg["download"]["max_per_species"] == 3000

    def test_smart_default_concat_respects_explicit_max_per_species(self):
        # Stale user file with max_per_species=300 → smart default still
        # bumps because 300 is the inherited DEFAULT_FAMILY_CONFIG value, not
        # an explicit user choice.  But _merge_with_defaults applies user file
        # AFTER smart defaults, so an explicit value wins overall.
        sample = next(iter(CONCATENATION_FAMILIES))
        merged = _merge_with_defaults(
            {"download": {"max_per_species": 1000}}, MINIMAL_GLOBAL, sample,
        )
        assert merged["download"]["max_per_species"] == 1000

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


# ---------------------------------------------------------------------------
# manual block — include/exclude lists
# ---------------------------------------------------------------------------

class TestManualBlock:
    def test_default_family_config_has_manual_key(self):
        assert "manual" in DEFAULT_FAMILY_CONFIG
        assert DEFAULT_FAMILY_CONFIG["manual"]["include"] == []
        assert DEFAULT_FAMILY_CONFIG["manual"]["exclude"] == []

    def test_empty_lists_pass(self):
        cfg = {"manual": {"include": [], "exclude": []}}
        _validate_manual_block(cfg, "Test")  # no exception

    def test_missing_block_treated_as_empty(self):
        cfg = {}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"] == {"include": [], "exclude": [], "include_fasta": [], "include_species": []}

    def test_include_and_exclude_overlap_rejected(self):
        cfg = {"manual": {"include": ["NC_001.1"], "exclude": ["NC_001.1"]}}
        with pytest.raises(ValueError, match="overlap"):
            _validate_manual_block(cfg, "Test")

    def test_non_string_entry_rejected(self):
        cfg = {"manual": {"include": [123], "exclude": []}}
        with pytest.raises(ValueError, match="string accession"):
            _validate_manual_block(cfg, "Test")

    def test_empty_string_entry_rejected(self):
        cfg = {"manual": {"include": ["  "], "exclude": []}}
        with pytest.raises(ValueError, match="empty"):
            _validate_manual_block(cfg, "Test")

    def test_non_list_block_rejected(self):
        cfg = {"manual": {"include": "NC_001.1", "exclude": []}}
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_whitespace_stripped(self):
        cfg = {"manual": {"include": ["  NC_001.1  "], "exclude": []}}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include"] == ["NC_001.1"]

    def test_duplicates_deduped(self):
        cfg = {"manual": {"include": ["NC_001.1", "NC_001.1"], "exclude": []}}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include"] == ["NC_001.1"]

    def test_loaded_config_has_empty_manual_block(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        # Stale config without a manual: block should still load fine.
        (cfg_dir / "Flaviviridae.yaml").write_text(
            yaml.dump({"sequence": {"type": "nucleotide", "region": "whole_genome"}})
        )
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["include"] == []
        assert cfg["manual"]["exclude"] == []

    def test_load_rejects_yaml_with_overlap(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": ["NC_X.1"], "exclude": ["NC_X.1"]},
        }))
        with pytest.raises(ValueError, match="overlap"):
            load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)


# ---------------------------------------------------------------------------
# manual.include_fasta — pasted-sequence entries
# ---------------------------------------------------------------------------

def _wrap_fasta(entries):
    """Helper: wrap include_fasta entries in a minimal cfg."""
    return {
        "sequence": {"type": "nucleotide", "region": "whole_genome"},
        "manual": {"include": [], "exclude": [], "include_fasta": entries},
    }


class TestManualIncludeFasta:
    def test_default_has_empty_include_fasta(self):
        assert DEFAULT_FAMILY_CONFIG["manual"]["include_fasta"] == []

    def test_empty_list_passes(self):
        cfg = _wrap_fasta([])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta"] == []

    def test_valid_entry_passes(self):
        cfg = _wrap_fasta([
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGTACGT"},
        ])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta"] == [
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGTACGT"},
        ]

    def test_sequence_whitespace_stripped_and_uppercased(self):
        cfg = _wrap_fasta([
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "acgt\n acgt\n"},
        ])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta"][0]["sequence"] == "ACGTACGT"

    def test_id_and_organism_stripped(self):
        cfg = _wrap_fasta([
            {"id": "  MySeq1  ", "organism": "  Foo virus  ", "sequence": "ACGT"},
        ])
        _validate_manual_block(cfg, "Test")
        entry = cfg["manual"]["include_fasta"][0]
        assert entry["id"] == "MySeq1"
        assert entry["organism"] == "Foo virus"

    def test_non_list_rejected(self):
        cfg = _wrap_fasta("not a list")
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_non_mapping_entry_rejected(self):
        cfg = _wrap_fasta(["not a mapping"])
        with pytest.raises(ValueError, match="must be a mapping"):
            _validate_manual_block(cfg, "Test")

    def test_missing_id_rejected(self):
        cfg = _wrap_fasta([{"organism": "Foo virus", "sequence": "ACGT"}])
        with pytest.raises(ValueError, match=r"include_fasta\[0\]\.id"):
            _validate_manual_block(cfg, "Test")

    def test_missing_organism_rejected(self):
        cfg = _wrap_fasta([{"id": "X", "sequence": "ACGT"}])
        with pytest.raises(ValueError, match=r"include_fasta\[0\]\.organism"):
            _validate_manual_block(cfg, "Test")

    def test_missing_sequence_rejected(self):
        cfg = _wrap_fasta([{"id": "X", "organism": "Foo"}])
        with pytest.raises(ValueError, match=r"include_fasta\[0\]\.sequence"):
            _validate_manual_block(cfg, "Test")

    def test_empty_sequence_rejected(self):
        cfg = _wrap_fasta([{"id": "X", "organism": "Foo", "sequence": "   \n  "}])
        with pytest.raises(ValueError, match="empty"):
            _validate_manual_block(cfg, "Test")

    def test_duplicate_id_within_fasta_rejected(self):
        cfg = _wrap_fasta([
            {"id": "X", "organism": "Foo", "sequence": "ACGT"},
            {"id": "X", "organism": "Bar", "sequence": "TGCA"},
        ])
        with pytest.raises(ValueError, match="duplicates an earlier entry"):
            _validate_manual_block(cfg, "Test")

    def test_id_collision_with_include_rejected(self):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {
                "include": ["NC_X.1"],
                "exclude": [],
                "include_fasta": [
                    {"id": "NC_X.1", "organism": "Foo", "sequence": "ACGT"},
                ],
            },
        }
        with pytest.raises(ValueError, match="overlap with manual.include"):
            _validate_manual_block(cfg, "Test")

    def test_id_collision_with_exclude_rejected(self):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {
                "include": [],
                "exclude": ["NC_X.1"],
                "include_fasta": [
                    {"id": "NC_X.1", "organism": "Foo", "sequence": "ACGT"},
                ],
            },
        }
        with pytest.raises(ValueError, match="overlap with manual.exclude"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_rejects_include_fasta(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {
                "include": [],
                "exclude": [],
                "include_fasta": [
                    {"id": "MySeq1", "organism": "Foo", "sequence": "MKL"},
                ],
            },
        }
        with pytest.raises(ValueError, match="not supported when sequence.region is 'concatenated'"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_with_empty_fasta_ok(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {"include": [], "exclude": [], "include_fasta": []},
        }
        _validate_manual_block(cfg, "Test")

    def test_loaded_stale_config_gets_empty_include_fasta(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": []},
        }))
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["include_fasta"] == []


# ---------------------------------------------------------------------------
# manual.include_species validation
# ---------------------------------------------------------------------------

class TestManualIncludeSpecies:
    def _base_cfg(self, species=None):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": [], "include_fasta": []},
        }
        if species is not None:
            cfg["manual"]["include_species"] = species
        return cfg

    def test_default_has_empty_include_species(self):
        cfg = self._base_cfg()
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == []

    def test_empty_list_passes(self):
        cfg = self._base_cfg(species=[])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == []

    def test_valid_name_string_accepted(self):
        cfg = self._base_cfg(species=["Dengue virus"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == ["Dengue virus"]

    def test_valid_taxid_integer_converted_to_string(self):
        cfg = self._base_cfg(species=[11103])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == ["11103"]

    def test_valid_taxid_digit_string_accepted(self):
        cfg = self._base_cfg(species=["11103"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == ["11103"]

    def test_non_list_rejected(self):
        cfg = self._base_cfg(species="Dengue virus")
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_empty_string_entry_rejected(self):
        cfg = self._base_cfg(species=[""])
        with pytest.raises(ValueError, match="is empty"):
            _validate_manual_block(cfg, "Test")

    def test_whitespace_only_string_rejected(self):
        cfg = self._base_cfg(species=["   "])
        with pytest.raises(ValueError, match="is empty"):
            _validate_manual_block(cfg, "Test")

    def test_wrong_type_float_rejected(self):
        cfg = self._base_cfg(species=[3.14])
        with pytest.raises(ValueError, match="must be a species name"):
            _validate_manual_block(cfg, "Test")

    def test_duplicate_entries_deduped_with_warning(self, caplog):
        import logging
        cfg = self._base_cfg(species=["Dengue virus", "Dengue virus"])
        with caplog.at_level(logging.WARNING):
            _validate_manual_block(cfg, "Flaviviridae")
        assert cfg["manual"]["include_species"] == ["Dengue virus"]
        assert "duplicate" in caplog.text.lower()

    def test_mixed_names_and_taxids_all_accepted(self):
        cfg = self._base_cfg(species=["Dengue virus", 11103, "1234567"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_species"] == ["Dengue virus", "11103", "1234567"]

    def test_loaded_stale_config_gets_empty_include_species(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": []},
        }))
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["include_species"] == []


# ---------------------------------------------------------------------------
# File config overrides all internal defaults
#
# Priority chain (lowest → highest):
#   DEFAULT_FAMILY_CONFIG
#   → global_cfg["defaults"]
#   → smart defaults (DNA_FAMILIES / CONCATENATION_FAMILIES / SEGMENTED_FAMILIES)
#   → family-specific config file
#
# A file present in configs/ must win at every layer for every key it sets.
# Keys absent from the file are filled from the defaults (deep merge, not replace).
# ---------------------------------------------------------------------------

class TestFileConfigOverridesAllDefaults:

    # --- DEFAULT_FAMILY_CONFIG layer ---

    def test_overrides_default_max_per_species(self):
        merged = _merge_with_defaults(
            {"download": {"max_per_species": 99}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["download"]["max_per_species"] == 99

    def test_overrides_default_clustering_tool(self):
        merged = _merge_with_defaults(
            {"clustering": {"tool": "cdhit"}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["clustering"]["tool"] == "cdhit"

    def test_overrides_default_clustering_thresholds(self):
        merged = _merge_with_defaults(
            {"clustering": {"threshold_min": 0.80, "threshold_max": 0.95}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["clustering"]["threshold_min"] == 0.80
        assert merged["clustering"]["threshold_max"] == 0.95

    def test_overrides_default_max_reps(self):
        merged = _merge_with_defaults(
            {"clustering": {"max_reps_500": 7, "max_reps_100": 2}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["clustering"]["max_reps_500"] == 7
        assert merged["clustering"]["max_reps_100"] == 2

    def test_overrides_default_quality_max_ambiguous(self):
        merged = _merge_with_defaults(
            {"quality": {"max_ambiguous": 0.05}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["quality"]["max_ambiguous"] == 0.05

    def test_overrides_default_quality_min_length(self):
        merged = _merge_with_defaults(
            {"quality": {"min_length": 500}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["quality"]["min_length"] == 500

    def test_overrides_default_targets(self):
        merged = _merge_with_defaults(
            {"targets": {"max_500": 300, "max_100": 50}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["targets"]["max_500"] == 300
        assert merged["targets"]["max_100"] == 50

    def test_overrides_default_msa_500_tool_and_options(self):
        merged = _merge_with_defaults(
            {"msa_500": {"tool": "muscle", "options_nuc": "--custom"}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["msa_500"]["tool"] == "muscle"
        assert merged["msa_500"]["options_nuc"] == "--custom"

    def test_overrides_default_msa_100_options(self):
        merged = _merge_with_defaults(
            {"msa_100": {"options_nuc": "--retree 1"}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["msa_100"]["options_nuc"] == "--retree 1"

    def test_overrides_default_tree_500_tool_and_model(self):
        merged = _merge_with_defaults(
            {"tree_500": {"tool": "iqtree", "model_nuc": "HKY+G"}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["tree_500"]["tool"] == "iqtree"
        assert merged["tree_500"]["model_nuc"] == "HKY+G"

    def test_overrides_default_tree_100_model(self):
        merged = _merge_with_defaults(
            {"tree_100": {"model_nuc": "K80+G"}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["tree_100"]["model_nuc"] == "K80+G"

    def test_overrides_default_msa_trim_disabled(self):
        merged = _merge_with_defaults(
            {"msa_trim": {"enabled": False}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["msa_trim"]["enabled"] is False

    def test_overrides_default_outlier_removal_factor(self):
        merged = _merge_with_defaults(
            {"outlier_removal": {"factor": 10.0, "max_iterations": 1}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["outlier_removal"]["factor"] == 10.0
        assert merged["outlier_removal"]["max_iterations"] == 1

    def test_overrides_default_length_outlier_k(self):
        merged = _merge_with_defaults(
            {"length_outlier": {"k": 3.0, "enabled": False}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["length_outlier"]["k"] == 3.0
        assert merged["length_outlier"]["enabled"] is False

    def test_overrides_default_coloring_genus_inference(self):
        merged = _merge_with_defaults(
            {"coloring": {"genus_inference": "none"}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["coloring"]["genus_inference"] == "none"

    def test_overrides_default_taxonomy_lca_min_rank(self):
        merged = _merge_with_defaults(
            {"taxonomy": {"lca_min_rank": "genus"}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["taxonomy"]["lca_min_rank"] == "genus"

    def test_overrides_default_refseq_absorption_disabled(self):
        merged = _merge_with_defaults(
            {"refseq_absorption": {"enabled": False}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        assert merged["refseq_absorption"]["enabled"] is False

    # --- global_cfg["defaults"] layer ---

    def test_overrides_global_defaults(self):
        global_with_defaults = {
            "ncbi": {"email": "test@test.com"},
            "defaults": {
                "clustering": {"tool": "cdhit"},
                "download": {"max_per_species": 500},
            },
        }
        merged = _merge_with_defaults(
            {"clustering": {"tool": "mmseqs2"}, "download": {"max_per_species": 77}},
            global_with_defaults, "Flaviviridae",
        )
        assert merged["clustering"]["tool"] == "mmseqs2"
        assert merged["download"]["max_per_species"] == 77

    def test_global_defaults_take_effect_when_file_is_silent(self):
        global_with_defaults = {
            "ncbi": {"email": "test@test.com"},
            "defaults": {"clustering": {"tool": "cdhit"}},
        }
        merged = _merge_with_defaults({}, global_with_defaults, "Flaviviridae")
        assert merged["clustering"]["tool"] == "cdhit"

    # --- smart defaults (DNA_FAMILIES) layer ---

    def test_overrides_dna_families_region(self):
        # Adenoviridae smart default is "hexon"; file says "penton"
        merged = _merge_with_defaults(
            {"sequence": {"region": "penton"}}, MINIMAL_GLOBAL, "Adenoviridae"
        )
        assert merged["sequence"]["region"] == "penton"

    def test_overrides_dna_families_type(self):
        # Adenoviridae smart default is protein; file switches to nucleotide
        merged = _merge_with_defaults(
            {"sequence": {"type": "nucleotide", "region": "hexon"}},
            MINIMAL_GLOBAL, "Adenoviridae",
        )
        assert merged["sequence"]["type"] == "nucleotide"

    # --- smart defaults (SEGMENTED_FAMILIES) layer ---

    def test_overrides_segmented_families_segment(self):
        # Hantaviridae smart default is "segment L"; file overrides it
        merged = _merge_with_defaults(
            {"sequence": {"segment": "segment S"}}, MINIMAL_GLOBAL, "Hantaviridae"
        )
        assert merged["sequence"]["segment"] == "segment S"

    # --- smart defaults (CONCATENATION_FAMILIES) layer ---

    def test_overrides_concat_families_region(self):
        # Poxviridae smart default is "concatenated"; file reverts to single gene
        merged = _merge_with_defaults(
            {"sequence": {"region": "rpo147", "type": "protein"}},
            MINIMAL_GLOBAL, "Poxviridae",
        )
        assert merged["sequence"]["region"] == "rpo147"

    def test_overrides_concat_families_min_fraction(self):
        merged = _merge_with_defaults(
            {"concatenation": {"min_fraction": 0.5}},
            MINIMAL_GLOBAL, "Poxviridae",
        )
        assert merged["concatenation"]["min_fraction"] == 0.5

    # --- deep-merge: absent keys are NOT wiped ---

    def test_absent_sibling_keys_retain_defaults(self):
        # File only sets clustering.tool; every other clustering key should
        # still come from DEFAULT_FAMILY_CONFIG, not disappear.
        merged = _merge_with_defaults(
            {"clustering": {"tool": "cdhit"}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        dc = DEFAULT_FAMILY_CONFIG["clustering"]
        assert merged["clustering"]["threshold_min"] == dc["threshold_min"]
        assert merged["clustering"]["threshold_max"] == dc["threshold_max"]
        assert merged["clustering"]["max_reps_500"] == dc["max_reps_500"]
        assert merged["clustering"]["max_reps_100"] == dc["max_reps_100"]

    def test_absent_top_level_keys_retain_defaults(self):
        # File only touches download; all other top-level sections should still
        # be present and match DEFAULT_FAMILY_CONFIG.
        merged = _merge_with_defaults(
            {"download": {"max_per_species": 1}}, MINIMAL_GLOBAL, "Flaviviridae"
        )
        for key in DEFAULT_FAMILY_CONFIG:
            assert key in merged, f"top-level key {key!r} missing from merged config"

    # --- end-to-end: load_family_config with file on disk ---

    def test_load_family_config_file_values_all_survive(self, tmp_path):
        """Full disk-to-merged-dict path: every non-default value in the file
        must appear unchanged in the final config dict."""
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        file_values = {
            "download":        {"max_per_species": 42},
            "clustering":      {"tool": "cdhit", "max_reps_500": 7, "max_reps_100": 2,
                                "threshold_min": 0.80, "threshold_max": 0.97},
            "quality":         {"max_ambiguous": 0.03, "min_length": 800},
            "targets":         {"max_500": 250, "max_100": 60},
            "msa_500":         {"tool": "muscle", "options_nuc": "--custom-nuc"},
            "msa_100":         {"options_nuc": "--retree 1"},
            "msa_trim":        {"enabled": False},
            "tree_500":        {"tool": "iqtree2", "model_nuc": "HKY+G"},
            "tree_100":        {"model_nuc": "K80+G", "options_nuc": "--fast"},
            "outlier_removal": {"factor": 15.0, "max_iterations": 1, "min_seqs": 20},
            "length_outlier":  {"enabled": False, "k": 3.0},
            "coloring":        {"genus_inference": "suffix"},
            "taxonomy":        {"lca_min_rank": "genus"},
            "refseq_absorption": {"enabled": False, "threshold": 0.95},
            "sequence":        {"type": "nucleotide", "region": "whole_genome"},
        }
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump(file_values))
        cfg, auto = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert auto is False
        assert cfg["download"]["max_per_species"] == 42
        assert cfg["clustering"]["tool"] == "cdhit"
        assert cfg["clustering"]["max_reps_500"] == 7
        assert cfg["clustering"]["max_reps_100"] == 2
        assert cfg["clustering"]["threshold_min"] == 0.80
        assert cfg["clustering"]["threshold_max"] == 0.97
        assert cfg["quality"]["max_ambiguous"] == 0.03
        assert cfg["quality"]["min_length"] == 800
        assert cfg["targets"]["max_500"] == 250
        assert cfg["targets"]["max_100"] == 60
        assert cfg["msa_500"]["tool"] == "muscle"
        assert cfg["msa_500"]["options_nuc"] == "--custom-nuc"
        assert cfg["msa_100"]["options_nuc"] == "--retree 1"
        assert cfg["msa_trim"]["enabled"] is False
        assert cfg["tree_500"]["tool"] == "iqtree2"
        assert cfg["tree_500"]["model_nuc"] == "HKY+G"
        assert cfg["tree_100"]["model_nuc"] == "K80+G"
        assert cfg["tree_100"]["options_nuc"] == "--fast"
        assert cfg["outlier_removal"]["factor"] == 15.0
        assert cfg["outlier_removal"]["max_iterations"] == 1
        assert cfg["outlier_removal"]["min_seqs"] == 20
        assert cfg["length_outlier"]["enabled"] is False
        assert cfg["length_outlier"]["k"] == 3.0
        assert cfg["coloring"]["genus_inference"] == "suffix"
        assert cfg["taxonomy"]["lca_min_rank"] == "genus"
        assert cfg["refseq_absorption"]["enabled"] is False
        assert cfg["refseq_absorption"]["threshold"] == 0.95
