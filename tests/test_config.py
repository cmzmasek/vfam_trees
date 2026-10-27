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
    RDRP_DOMAIN_PADDING,
    RDRP_FAMILIES,
    RDRP_FAMILY_PROFILES,
    RDRP_HMM_TOKENS,
    RDRP_PALM,
    SEGMENTED_FAMILIES,
    _apply_smart_defaults,
    _deep_update,
    _generate_default_family_config,
    _merge_with_defaults,
    _validate_additional_data_sources,
    _validate_concatenation_block,
    _validate_hmm_block,
    _validate_tree_blocks,
    _validate_manual_block,
    _validate_numeric_fields,
    _warn_smart_default_conflicts,
    load_family_config,
    make_minimal_global_cfg,
    sanitize_output_prefix,
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
        _apply_smart_defaults("Phasmaviridae", cfg, file_cfg={})
    assert any("Auto-configured segment" in r.getMessage() for r in caplog.records)


def test_auto_configured_segment_suppressed_when_file_overrides_it(caplog):
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    file_cfg = {"sequence": {"segment": "segment S"}}
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Phasmaviridae", cfg, file_cfg=file_cfg)
    assert not any("Auto-configured segment" in r.getMessage() for r in caplog.records)


def test_rdrp_family_logs_no_segment(caplog):
    # The RdRp preset clears the segment, so none is announced.
    import logging
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    with caplog.at_level(logging.INFO, logger="vfam_trees.config"):
        _apply_smart_defaults("Hantaviridae", cfg, file_cfg={})
    assert not any("Auto-configured segment" in r.getMessage() for r in caplog.records)
    assert cfg["sequence"]["segment"] is None


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


def test_apply_smart_defaults_family_without_preset_keeps_whole_genome():
    # Kolmioviridae (deltaviruses) is in no preset: template whole genome.
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Kolmioviridae", cfg)
    assert cfg["sequence"]["region"] == "whole_genome"
    assert cfg["sequence"]["type"] == "nucleotide"


# Phasmaviridae: segmented, and not on the RdRp preset (its RdRps are too
# divergent for the Pfam set), so it keeps nucleotide segment-L analysis.
def test_apply_smart_defaults_segmented():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    cfg["sequence"]["segment"] = None
    _apply_smart_defaults("Phasmaviridae", cfg)
    assert cfg["sequence"]["segment"] == SEGMENTED_FAMILIES["Phasmaviridae"]


def test_apply_smart_defaults_segmented_does_not_overwrite_explicit():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    cfg["sequence"]["segment"] = "custom_segment"
    _apply_smart_defaults("Phasmaviridae", cfg)
    assert cfg["sequence"]["segment"] == "custom_segment"


def test_apply_smart_defaults_rna_family_off_the_rdrp_preset_unchanged():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Narnaviridae", cfg)
    assert cfg["sequence"]["region"] == "whole_genome"
    assert cfg["hmm"]["enabled"] is False


def test_apply_smart_defaults_rdrp_family():
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Flaviviridae", cfg)
    assert cfg["sequence"]["type"] == "protein"
    assert cfg["sequence"]["region"] == "RdRp"
    assert cfg["sequence"]["hmms"] == list(RDRP_FAMILY_PROFILES["Flaviviridae"])
    assert cfg["sequence"]["domain_padding"] == RDRP_DOMAIN_PADDING
    assert cfg["hmm"]["enabled"] is True and cfg["hmm"]["database"] == "RdRp"
    assert cfg["hmm"]["best_hit_only"] is True


def test_rdrp_preset_clears_the_segment():
    # Segmented RNA families are on SEGMENTED_FAMILIES for nucleotide mode;
    # the HMM finds the RdRp on whichever segment carries it.
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _apply_smart_defaults("Hantaviridae", cfg)
    assert "Hantaviridae" in SEGMENTED_FAMILIES
    assert cfg["sequence"]["segment"] is None
    assert cfg["sequence"]["region"] == "RdRp"


class TestRdrpOptOut:
    """A family file that chooses its own analysis gets no RdRp preset at
    all — never a mix of nucleotide type with the preset's HMM protein fetch."""

    def test_type_nucleotide_alone_restores_nucleotide_with_segment(self):
        # e.g. the repo's pre-v1.3 configs/Hantaviridae.yaml: no hmm block.
        merged = _merge_with_defaults(
            {"sequence": {"type": "nucleotide"}}, MINIMAL_GLOBAL, "Hantaviridae")
        assert merged["sequence"]["type"] == "nucleotide"
        assert merged["sequence"]["region"] == "whole_genome"
        assert merged["sequence"]["segment"] == SEGMENTED_FAMILIES["Hantaviridae"]
        assert merged["hmm"]["enabled"] is False
        assert not merged["sequence"].get("hmms")

    def test_hmm_disabled_alone_opts_out(self):
        merged = _merge_with_defaults({"hmm": {"enabled": False}}, MINIMAL_GLOBAL, "Filoviridae")
        assert merged["sequence"]["type"] == "nucleotide"
        assert merged["sequence"]["region"] == "whole_genome"
        assert merged["hmm"]["enabled"] is False

    def test_v150_snapshot_runs_as_before(self):
        snapshot = {"sequence": {"type": "nucleotide", "region": "whole_genome",
                                 "segment": "segment L"},
                    "hmm": {"enabled": False, "database": None}}
        merged = _merge_with_defaults(snapshot, MINIMAL_GLOBAL, "Hantaviridae")
        assert merged["sequence"]["segment"] == "segment L"
        assert merged["hmm"]["enabled"] is False

    def test_v160_file_keeps_the_preset(self):
        file_cfg = {"sequence": {"type": "protein", "region": "RdRp"},
                    "hmm": {"enabled": True}}
        merged = _merge_with_defaults(file_cfg, MINIMAL_GLOBAL, "Hantaviridae")
        assert merged["hmm"]["database"] == "RdRp"
        assert merged["sequence"]["segment"] is None

    def test_unrelated_overrides_keep_the_preset(self):
        merged = _merge_with_defaults({"download": {"max_per_species": 50}},
                                      MINIMAL_GLOBAL, "Togaviridae")
        assert merged["sequence"]["type"] == "protein"
        assert merged["hmm"]["enabled"] is True

    def test_warning_names_what_runs(self, caplog, tmp_path):
        import logging
        path = tmp_path / "Hantaviridae.yaml"
        with caplog.at_level(logging.WARNING):
            _warn_smart_default_conflicts(
                "Hantaviridae", {"sequence": {"type": "nucleotide"}}, path)
        assert "sequence type nucleotide, region whole_genome" in caplog.text
        assert repr(SEGMENTED_FAMILIES["Hantaviridae"]) in caplog.text
        assert "keeps whole-sequence" not in caplog.text

    def test_no_warning_for_a_preset_file(self, caplog, tmp_path):
        import logging
        with caplog.at_level(logging.WARNING):
            _warn_smart_default_conflicts(
                "Hantaviridae", {"sequence": {"type": "protein", "region": "RdRp"}},
                tmp_path / "Hantaviridae.yaml")
        assert caplog.text == ""


class TestRdrpPresetExtras:
    def test_sibling_domains_count_for_best_hit(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        _apply_smart_defaults("Nairoviridae", cfg)
        assert cfg["sequence"]["hmms"] == ["Bunya_RdRp"]
        assert cfg["hmm"]["best_hit_also"] == ["L_PA-C-like", "L_thumb_ring_vir"]

    def test_no_siblings_listed_twice(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        _apply_smart_defaults("Flaviviridae", cfg)   # has both NS5 profiles already
        assert cfg["hmm"]["best_hit_also"] == []

    def test_min_protein_length(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        _apply_smart_defaults("Togaviridae", cfg)
        assert cfg["hmm"]["min_protein_length"] == 100
        assert DEFAULT_FAMILY_CONFIG["hmm"]["min_protein_length"] is None

    def test_template_documents_domain_padding(self):
        assert "domain_padding" in DEFAULT_FAMILY_CONFIG["sequence"]
        assert DEFAULT_FAMILY_CONFIG["sequence"]["domain_padding"] is None

    def test_new_hmm_keys_validated(self):
        for bad in ({"best_hit_also": "L_PA-C-like"}, {"best_hit_also": [""]},
                    {"min_protein_length": 0}, {"min_protein_length": True}):
            cfg = {"hmm": {"enabled": False, **bad}}
            with pytest.raises(ValueError):
                _validate_hmm_block(cfg, "X")
        cfg = {"hmm": {"enabled": False, "best_hit_also": [" L_PA-C-like "],
                       "min_protein_length": 150}}
        _validate_hmm_block(cfg, "X")
        assert cfg["hmm"]["best_hit_also"] == ["L_PA-C-like"]
        assert cfg["hmm"]["min_protein_length"] == 150


class TestInjectionTypeChecks:
    def test_external_sources_refused_for_protein_family(self):
        cfg = _merge_with_defaults(
            {"additional_data_sources": [{"source": "pathoplexus", "organism": "ebola-zaire"}]},
            MINIMAL_GLOBAL, "Filoviridae")
        with pytest.raises(ValueError, match="analysed as protein.*RdRp marker preset"):
            _validate_additional_data_sources(cfg, "Filoviridae")

    def test_external_sources_fine_for_nucleotide_opt_out(self):
        cfg = _merge_with_defaults(
            {"sequence": {"type": "nucleotide"}, "hmm": {"enabled": False},
             "additional_data_sources": [{"source": "pathoplexus", "organism": "ebola-zaire"}]},
            MINIMAL_GLOBAL, "Filoviridae")
        _validate_additional_data_sources(cfg, "Filoviridae")
        assert cfg["additional_data_sources"][0]["organism"] == "ebola-zaire"

    def test_pasted_nucleotide_refused_for_protein_family(self):
        cfg = _merge_with_defaults(
            {"manual": {"include_seq": [{"id": "x1", "organism": "Ebola virus",
                                         "sequence": "ACGT" * 100}]}},
            MINIMAL_GLOBAL, "Filoviridae")
        with pytest.raises(ValueError, match=r"include_seq\[0\].*nucleotide"):
            _validate_manual_block(cfg, "Filoviridae")

    def test_pasted_protein_fine_for_protein_family(self):
        cfg = _merge_with_defaults(
            {"manual": {"include_seq": [{"id": "x1", "organism": "Ebola virus",
                                         "sequence": "MSLEKTFQRIVDGAPWYCHN" * 10}]}},
            MINIMAL_GLOBAL, "Filoviridae")
        _validate_manual_block(cfg, "Filoviridae")


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
    # The preset turns on HMM (hexon); with the file's nucleotide type that
    # would fetch hexon proteins and reject them all in nucleotide QC, so
    # loading stops and says how to switch the family to nucleotide.
    with pytest.raises(ValueError, match="hmm.enabled: false"):
        load_family_config("Adenoviridae", cfg_dir, MINIMAL_GLOBAL)
    (cfg_dir / "Adenoviridae.yaml").write_text(yaml.dump({
        "sequence": {"type": "nucleotide", "region": "whole_genome", "segment": None},
        "hmm": {"enabled": False},
    }))
    cfg, auto = load_family_config("Adenoviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is False
    # Smart default (hexon) sits between DEFAULT and user file in priority;
    # user file explicitly said whole_genome, so whole_genome wins — but a
    # warning should have been issued (tested separately).
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
    cfg, auto = load_family_config("Phasmaviridae", cfg_dir, MINIMAL_GLOBAL)
    assert auto is True
    assert cfg["sequence"]["segment"] == SEGMENTED_FAMILIES["Phasmaviridae"]


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
            "Pithoviridae", "Pandoraviridae",
        ):
            assert fam in CONCATENATION_FAMILIES, f"missing preset: {fam}"
        assert "Malacoherpesviridae" not in CONCATENATION_FAMILIES
        # Mamonoviridae runs on one HMM marker (DNA_FAMILIES); "Medusaviridae"
        # was a name NCBI never used.
        assert "Mamonoviridae" not in CONCATENATION_FAMILIES
        assert "Medusaviridae" not in CONCATENATION_FAMILIES

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

    def test_ncldv_fallback_only_pandora(self):
        # After the 2026-05 audit only Pandoraviridae remains on the 8-marker
        # NCLDV-hallmark fallback; Mimi, Phyco, Marseille, and Pitho got
        # pruned per-family presets, Mamonoviridae a single HMM marker.
        pandora = CONCATENATION_FAMILIES["Pandoraviridae"]["concatenation"]["proteins"]
        assert len(pandora) == 8
        # Pruned-preset families have fewer markers.
        for fam, expected_n in (
            ("Phycodnaviridae", 4),
            ("Mimiviridae", 6),
            ("Marseilleviridae", 6),
            ("Pithoviridae", 7),     # HMM-identified since 1.6.2
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
        assert DEFAULT_FAMILY_CONFIG["manual"]["name"] == ""

    def test_empty_lists_pass(self):
        cfg = {"manual": {"include": [], "exclude": []}}
        _validate_manual_block(cfg, "Test")  # no exception

    def test_missing_block_treated_as_empty(self):
        cfg = {}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"] == {"include": [], "exclude": [], "include_seq": [], "include_fasta_files": [], "restrict_to_lineages": [], "name": ""}

    def test_name_string_is_accepted(self):
        cfg = {"manual": {"name": "  Hantaviridae 2026  "}}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["name"] == "Hantaviridae 2026"

    def test_name_empty_string_is_accepted(self):
        cfg = {"manual": {"name": ""}}
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["name"] == ""

    def test_name_non_string_rejected(self):
        cfg = {"manual": {"name": 42}}
        with pytest.raises(ValueError, match="manual.name must be a string"):
            _validate_manual_block(cfg, "Test")

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
# manual.include_seq — pasted-sequence entries
# ---------------------------------------------------------------------------

def _wrap_seq(entries):
    """Helper: wrap include_seq entries in a minimal cfg."""
    return {
        "sequence": {"type": "nucleotide", "region": "whole_genome"},
        "manual": {"include": [], "exclude": [], "include_seq": entries},
    }


class TestManualIncludeSeq:
    def test_default_has_empty_include_seq(self):
        assert DEFAULT_FAMILY_CONFIG["manual"]["include_seq"] == []

    def test_empty_list_passes(self):
        cfg = _wrap_seq([])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_seq"] == []

    def test_valid_entry_passes(self):
        cfg = _wrap_seq([
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGTACGT"},
        ])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_seq"] == [
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGTACGT"},
        ]

    def test_sequence_whitespace_stripped_and_uppercased(self):
        cfg = _wrap_seq([
            {"id": "MySeq1", "organism": "Foo virus", "sequence": "acgt\n acgt\n"},
        ])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_seq"][0]["sequence"] == "ACGTACGT"

    def test_id_and_organism_stripped(self):
        cfg = _wrap_seq([
            {"id": "  MySeq1  ", "organism": "  Foo virus  ", "sequence": "ACGT"},
        ])
        _validate_manual_block(cfg, "Test")
        entry = cfg["manual"]["include_seq"][0]
        assert entry["id"] == "MySeq1"
        assert entry["organism"] == "Foo virus"

    def test_non_list_rejected(self):
        cfg = _wrap_seq("not a list")
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_non_mapping_entry_rejected(self):
        cfg = _wrap_seq(["not a mapping"])
        with pytest.raises(ValueError, match="must be a mapping"):
            _validate_manual_block(cfg, "Test")

    def test_missing_id_rejected(self):
        cfg = _wrap_seq([{"organism": "Foo virus", "sequence": "ACGT"}])
        with pytest.raises(ValueError, match=r"include_seq\[0\]\.id"):
            _validate_manual_block(cfg, "Test")

    def test_missing_organism_rejected(self):
        cfg = _wrap_seq([{"id": "X", "sequence": "ACGT"}])
        with pytest.raises(ValueError, match=r"include_seq\[0\]\.organism"):
            _validate_manual_block(cfg, "Test")

    def test_missing_sequence_rejected(self):
        cfg = _wrap_seq([{"id": "X", "organism": "Foo"}])
        with pytest.raises(ValueError, match=r"include_seq\[0\]\.sequence"):
            _validate_manual_block(cfg, "Test")

    def test_empty_sequence_rejected(self):
        cfg = _wrap_seq([{"id": "X", "organism": "Foo", "sequence": "   \n  "}])
        with pytest.raises(ValueError, match="empty"):
            _validate_manual_block(cfg, "Test")

    def test_duplicate_id_within_seq_rejected(self):
        cfg = _wrap_seq([
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
                "include_seq": [
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
                "include_seq": [
                    {"id": "NC_X.1", "organism": "Foo", "sequence": "ACGT"},
                ],
            },
        }
        with pytest.raises(ValueError, match="overlap with manual.exclude"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_rejects_include_seq(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {
                "include": [],
                "exclude": [],
                "include_seq": [
                    {"id": "MySeq1", "organism": "Foo", "sequence": "MKL"},
                ],
            },
        }
        with pytest.raises(ValueError, match="not supported when sequence.region is 'concatenated'"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_with_empty_seq_ok(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {"include": [], "exclude": [], "include_seq": []},
        }
        _validate_manual_block(cfg, "Test")

    def test_loaded_stale_config_gets_empty_include_seq(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": []},
        }))
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["include_seq"] == []


# ---------------------------------------------------------------------------
# manual.include_fasta_files — external FASTA file paths
# ---------------------------------------------------------------------------

def _wrap_fasta_files(paths):
    """Helper: wrap include_fasta_files paths in a minimal cfg."""
    return {
        "sequence": {"type": "nucleotide", "region": "whole_genome"},
        "manual": {"include": [], "exclude": [], "include_fasta_files": paths},
    }


class TestManualIncludeFastaFiles:
    def test_default_has_empty_include_fasta_files(self):
        assert DEFAULT_FAMILY_CONFIG["manual"]["include_fasta_files"] == []

    def test_empty_list_passes(self):
        cfg = _wrap_fasta_files([])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta_files"] == []

    def test_non_empty_string_list_passes(self):
        cfg = _wrap_fasta_files(["/path/to/seqs.fasta", "/other/file.fa"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta_files"] == ["/path/to/seqs.fasta", "/other/file.fa"]

    def test_whitespace_stripped(self):
        cfg = _wrap_fasta_files(["  /path/to/seqs.fasta  "])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["include_fasta_files"] == ["/path/to/seqs.fasta"]

    def test_non_list_rejected(self):
        cfg = _wrap_fasta_files("/path/to/seqs.fasta")
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_non_string_entry_rejected(self):
        cfg = _wrap_fasta_files([123])
        with pytest.raises(ValueError, match="path string"):
            _validate_manual_block(cfg, "Test")

    def test_empty_string_entry_rejected(self):
        cfg = _wrap_fasta_files(["   "])
        with pytest.raises(ValueError, match="empty"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_rejects_nonempty_fasta_files(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {
                "include": [],
                "exclude": [],
                "include_fasta_files": ["/some/file.fa"],
            },
        }
        with pytest.raises(ValueError, match="not supported when sequence.region is 'concatenated'"):
            _validate_manual_block(cfg, "Test")

    def test_concatenated_region_with_empty_fasta_files_ok(self):
        cfg = {
            "sequence": {"type": "protein", "region": "concatenated"},
            "manual": {"include": [], "exclude": [], "include_fasta_files": []},
        }
        _validate_manual_block(cfg, "Test")

    def test_loaded_stale_config_gets_empty_include_fasta_files(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": []},
        }))
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["include_fasta_files"] == []


# ---------------------------------------------------------------------------
# manual.restrict_to_lineages validation
# ---------------------------------------------------------------------------

class TestManualRestrictToLineages:
    def _base_cfg(self, lineages=None):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": [], "include_seq": []},
        }
        if lineages is not None:
            cfg["manual"]["restrict_to_lineages"] = lineages
        return cfg

    def test_default_has_empty_restrict_to_lineages(self):
        cfg = self._base_cfg()
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == []

    def test_empty_list_passes(self):
        cfg = self._base_cfg(lineages=[])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == []

    def test_valid_name_string_accepted(self):
        cfg = self._base_cfg(lineages=["Orthohantavirus"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == ["Orthohantavirus"]

    def test_valid_taxid_integer_converted_to_string(self):
        cfg = self._base_cfg(lineages=[11103])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == ["11103"]

    def test_valid_taxid_digit_string_accepted(self):
        cfg = self._base_cfg(lineages=["11103"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == ["11103"]

    def test_non_list_rejected(self):
        cfg = self._base_cfg(lineages="Orthohantavirus")
        with pytest.raises(ValueError, match="must be a list"):
            _validate_manual_block(cfg, "Test")

    def test_empty_string_entry_rejected(self):
        cfg = self._base_cfg(lineages=[""])
        with pytest.raises(ValueError, match="is empty"):
            _validate_manual_block(cfg, "Test")

    def test_whitespace_only_string_rejected(self):
        cfg = self._base_cfg(lineages=["   "])
        with pytest.raises(ValueError, match="is empty"):
            _validate_manual_block(cfg, "Test")

    def test_wrong_type_float_rejected(self):
        cfg = self._base_cfg(lineages=[3.14])
        with pytest.raises(ValueError, match="must be a taxon name"):
            _validate_manual_block(cfg, "Test")

    def test_duplicate_entries_deduped_with_warning(self, caplog):
        import logging
        cfg = self._base_cfg(lineages=["Orthohantavirus", "Orthohantavirus"])
        with caplog.at_level(logging.WARNING):
            _validate_manual_block(cfg, "Hantaviridae")
        assert cfg["manual"]["restrict_to_lineages"] == ["Orthohantavirus"]
        assert "duplicate" in caplog.text.lower()

    def test_mixed_names_and_taxids_all_accepted(self):
        cfg = self._base_cfg(lineages=["Orthohantavirus", 11103, "1234567"])
        _validate_manual_block(cfg, "Test")
        assert cfg["manual"]["restrict_to_lineages"] == ["Orthohantavirus", "11103", "1234567"]

    def test_loaded_stale_config_gets_empty_restrict_to_lineages(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Flaviviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {"include": [], "exclude": []},
        }))
        cfg, _ = load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)
        assert cfg["manual"]["restrict_to_lineages"] == []

    def test_legacy_include_species_key_is_hard_error(self):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {
                "include": [], "exclude": [], "include_seq": [],
                "include_species": ["Dengue virus"],
            },
        }
        with pytest.raises(ValueError, match="include_species has been renamed"):
            _validate_manual_block(cfg, "Flaviviridae")

    def test_legacy_limit_lineages_key_is_hard_error(self):
        cfg = {
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "manual": {
                "include": [], "exclude": [], "include_seq": [],
                "limit_lineages": ["Orthohantavirus"],
            },
        }
        with pytest.raises(ValueError, match="limit_lineages has been renamed"):
            _validate_manual_block(cfg, "Hantaviridae")


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

    def test_overrides_default_tree_500_options_nuc_and_aa(self):
        merged = _merge_with_defaults(
            {"tree_500": {"options_nuc": "-fastest", "options_aa": "-lg"}},
            MINIMAL_GLOBAL, "Flaviviridae",
        )
        assert merged["tree_500"]["options_nuc"] == "-fastest"
        assert merged["tree_500"]["options_aa"] == "-lg"
        assert "options" not in merged["tree_500"]

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


# ---------------------------------------------------------------------------
# _validate_numeric_fields — catch YAML-typo'd numbers at load time
# ---------------------------------------------------------------------------

class TestValidateNumericFields:
    def test_default_config_passes(self):
        cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
        _validate_numeric_fields(cfg, "Test")

    def test_empty_cfg_passes(self):
        _validate_numeric_fields({}, "Test")

    def test_string_with_embedded_spaces_rejected(self):
        # The original bug: `max_100: 4  00` in YAML parses as the string "4  00".
        cfg = {"targets": {"max_100": "4  00"}}
        with pytest.raises(ValueError, match=r"targets\.max_100.*integer.*str.*'4  00'"):
            _validate_numeric_fields(cfg, "Filoviridae")

    def test_string_with_digits_includes_typo_hint(self):
        cfg = {"targets": {"max_500": "1 000"}}
        with pytest.raises(ValueError, match="YAML typo"):
            _validate_numeric_fields(cfg, "Test")

    def test_float_rejected_for_int_field(self):
        cfg = {"targets": {"max_500": 500.5}}
        with pytest.raises(ValueError, match=r"targets\.max_500.*integer.*float"):
            _validate_numeric_fields(cfg, "Test")

    def test_bool_rejected_for_int_field(self):
        cfg = {"targets": {"max_500": True}}
        with pytest.raises(ValueError, match="must be a int"):
            _validate_numeric_fields(cfg, "Test")

    def test_negative_int_rejected(self):
        cfg = {"targets": {"max_500": -1}}
        with pytest.raises(ValueError, match=r"targets\.max_500=-1.*must be >= 1"):
            _validate_numeric_fields(cfg, "Test")

    def test_zero_rejected_for_positive_field(self):
        cfg = {"targets": {"max_100": 0}}
        with pytest.raises(ValueError, match=r"max_100=0.*must be >= 1"):
            _validate_numeric_fields(cfg, "Test")

    def test_fraction_out_of_range_high(self):
        cfg = {"quality": {"max_ambiguous": 1.5}}
        with pytest.raises(ValueError, match=r"max_ambiguous=1.5.*must be <= 1"):
            _validate_numeric_fields(cfg, "Test")

    def test_clustering_threshold_strictly_positive(self):
        cfg = {"clustering": {"threshold_min": 0.0}}
        with pytest.raises(ValueError, match=r"threshold_min=0.*must be > 0"):
            _validate_numeric_fields(cfg, "Test")

    def test_float_accepted_for_number_field(self):
        cfg = {"length_outlier": {"k": 4.5}}
        _validate_numeric_fields(cfg, "Test")

    def test_int_accepted_for_number_field(self):
        cfg = {"length_outlier": {"k": 5}}
        _validate_numeric_fields(cfg, "Test")

    def test_quality_min_length_none_accepted(self):
        cfg = {"quality": {"min_length": None}}
        _validate_numeric_fields(cfg, "Test")

    def test_quality_min_length_string_rejected(self):
        cfg = {"quality": {"min_length": "8 00"}}
        with pytest.raises(ValueError, match=r"quality\.min_length.*integer.*str"):
            _validate_numeric_fields(cfg, "Test")

    def test_load_family_config_aborts_on_typo(self, tmp_path):
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        (cfg_dir / "Filoviridae.yaml").write_text(yaml.dump({
            "sequence": {"type": "nucleotide", "region": "whole_genome"},
            "targets": {"max_500": 2000, "max_100": "4  00"},
        }))
        with pytest.raises(ValueError, match=r"targets\.max_100"):
            load_family_config("Filoviridae", cfg_dir, MINIMAL_GLOBAL)

    def test_auto_generated_config_passes_validation(self, tmp_path):
        # The auto-generation path also runs validation — make sure the
        # defaults plus smart-default overlays never trip it.
        cfg_dir = tmp_path / "configs"
        cfg_dir.mkdir()
        load_family_config("Hantaviridae", cfg_dir, MINIMAL_GLOBAL)
        load_family_config("Adenoviridae", cfg_dir, MINIMAL_GLOBAL)
        load_family_config("Flaviviridae", cfg_dir, MINIMAL_GLOBAL)


class TestSanitizeOutputPrefix:
    def test_clean_family_name_unchanged(self):
        # Families without manual.name must keep their exact output names.
        assert sanitize_output_prefix("Filoviridae") == "Filoviridae"
        assert sanitize_output_prefix("Picornaviridae") == "Picornaviridae"

    def test_spaces_become_underscores(self):
        assert sanitize_output_prefix("Ebola virus disease") == "Ebola_virus_disease"

    def test_path_separators_stripped(self):
        assert "/" not in sanitize_output_prefix("a/b/c")
        assert sanitize_output_prefix("a/b") == "a_b"

    def test_pipe_and_brackets_collapsed(self):
        assert sanitize_output_prefix("Ebola [concatenated|3 markers]") == "Ebola_concatenated_3_markers"

    def test_runs_collapse_to_single_underscore(self):
        assert sanitize_output_prefix("a   b") == "a_b"
        assert sanitize_output_prefix("a___b") == "a_b"

    def test_leading_trailing_separators_stripped(self):
        assert sanitize_output_prefix("  Ebola  ") == "Ebola"
        assert sanitize_output_prefix("__Ebola__") == "Ebola"
        assert sanitize_output_prefix("...Ebola...") == "Ebola"

    def test_hyphen_dot_underscore_preserved(self):
        assert sanitize_output_prefix("SARS-CoV-2_v1.2") == "SARS-CoV-2_v1.2"

    def test_empty_or_garbage_falls_back(self):
        assert sanitize_output_prefix("") == "family"
        assert sanitize_output_prefix("   ") == "family"
        assert sanitize_output_prefix("///") == "family"


# ---------------------------------------------------------------------------
# _validate_additional_data_sources
# ---------------------------------------------------------------------------

class TestValidateAdditionalDataSources:
    def _cfg(self, entries, region="whole_genome"):
        return {"sequence": {"region": region}, "additional_data_sources": entries}

    def test_default_in_template(self):
        assert DEFAULT_FAMILY_CONFIG["additional_data_sources"] == []

    def test_none_normalises_to_empty_list(self):
        cfg = {"sequence": {}}
        _validate_additional_data_sources(cfg, "Fam")
        assert cfg["additional_data_sources"] == []

    def test_happy_path_fills_defaults_and_types(self):
        cfg = self._cfg([
            {"source": "pathoplexus", "organism": "ebola-zaire",
             "country": "Uganda", "date_from": "2024-01-01", "max_seqs": 10}
        ])
        _validate_additional_data_sources(cfg, "Fam")
        entry = cfg["additional_data_sources"][0]
        assert entry == {
            "source": "pathoplexus", "organism": "ebola-zaire",
            "country": "Uganda", "host": None,
            "date_from": "2024-01-01", "date_to": None,
            "max_seqs": 10, "dedup_vs_ncbi": True,
            "name_prefix": "", "outbreak_name": "",
        }

    def test_not_a_list_rejected(self):
        with pytest.raises(ValueError, match="must be a list"):
            _validate_additional_data_sources(
                {"sequence": {}, "additional_data_sources": {"source": "pathoplexus"}},
                "Fam",
            )

    def test_unknown_source_rejected(self):
        with pytest.raises(ValueError, match="source must be one of"):
            _validate_additional_data_sources(
                self._cfg([{"source": "genbank", "organism": "ebola-zaire"}]), "Fam"
            )

    def test_unknown_organism_rejected(self):
        with pytest.raises(ValueError, match="not a known Pathoplexus organism"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebolx"}]), "Fam"
            )

    def test_missing_organism_rejected(self):
        with pytest.raises(ValueError, match="organism must be a non-empty"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus"}]), "Fam"
            )

    def test_bad_date_rejected(self):
        with pytest.raises(ValueError, match="ISO date"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "date_to": "2024-13-40"}]), "Fam"
            )

    def test_non_positive_max_seqs_rejected(self):
        with pytest.raises(ValueError, match="positive integer"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "max_seqs": 0}]), "Fam"
            )

    def test_bool_max_seqs_rejected(self):
        with pytest.raises(ValueError, match="positive integer"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "max_seqs": True}]), "Fam"
            )

    def test_non_bool_dedup_rejected(self):
        with pytest.raises(ValueError, match="must be a boolean"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "dedup_vs_ncbi": "yes"}]), "Fam"
            )

    def test_name_prefix_and_outbreak_name_accepted(self):
        cfg = self._cfg([
            {"source": "pathoplexus", "organism": "ebola-bdbv",
             "name_prefix": "PATHOPLEXUS_", "outbreak_name": "  Bdbv-2026  "}
        ])
        _validate_additional_data_sources(cfg, "Fam")
        entry = cfg["additional_data_sources"][0]
        assert entry["name_prefix"] == "PATHOPLEXUS_"   # verbatim, not stripped
        assert entry["outbreak_name"] == "Bdbv-2026"     # stripped

    def test_presentation_fields_default_empty(self):
        cfg = self._cfg([{"source": "pathoplexus", "organism": "ebola-bdbv"}])
        _validate_additional_data_sources(cfg, "Fam")
        entry = cfg["additional_data_sources"][0]
        assert entry["name_prefix"] == ""
        assert entry["outbreak_name"] == ""

    def test_non_string_name_prefix_rejected(self):
        with pytest.raises(ValueError, match="name_prefix must be a string"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "name_prefix": 123}]), "Fam"
            )

    def test_non_string_outbreak_name_rejected(self):
        with pytest.raises(ValueError, match="outbreak_name must be a string"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire",
                            "outbreak_name": ["x"]}]), "Fam"
            )

    def test_concat_region_rejected(self):
        with pytest.raises(ValueError, match="not supported when"):
            _validate_additional_data_sources(
                self._cfg([{"source": "pathoplexus", "organism": "ebola-zaire"}],
                          region="concatenated"),
                "Fam",
            )

    def test_empty_list_allowed_in_concat(self):
        # the concat guard only triggers when there are entries
        cfg = self._cfg([], region="concatenated")
        _validate_additional_data_sources(cfg, "Fam")
        assert cfg["additional_data_sources"] == []
        assert sanitize_output_prefix(None) == "family"


class TestValidateHmmBlock:
    def _cfg(self, hmm=None, *, sequence=None, concatenation=None):
        cfg = {"sequence": sequence if sequence is not None else {"region": "x"}}
        if hmm is not None:
            cfg["hmm"] = hmm
        if concatenation is not None:
            cfg["concatenation"] = concatenation
        return cfg

    def test_missing_block_fills_disabled_default(self):
        cfg = self._cfg(None)
        _validate_hmm_block(cfg, "Fam")
        assert cfg["hmm"]["enabled"] is False
        assert cfg["hmm"]["database"] is None
        assert cfg["hmm"]["default_evalue"] == 1.0e-5
        assert cfg["hmm"]["relative_length_cutoff"] == 0.5
        assert cfg["hmm"]["multidomain_overlap_tolerance"] == 30

    def test_enabled_with_database_ok(self):
        cfg = self._cfg({"enabled": True, "database": " Poxviridae "})
        _validate_hmm_block(cfg, "Fam")
        assert cfg["hmm"]["enabled"] is True
        assert cfg["hmm"]["database"] == "Poxviridae"  # stripped

    def test_enabled_without_database_raises(self):
        cfg = self._cfg({"enabled": True})
        with pytest.raises(ValueError, match="hmm.database is not set"):
            _validate_hmm_block(cfg, "Fam")

    def test_non_bool_enabled_raises(self):
        with pytest.raises(ValueError, match="hmm.enabled must be a boolean"):
            _validate_hmm_block(self._cfg({"enabled": "yes"}), "Fam")

    def test_non_positive_evalue_raises(self):
        with pytest.raises(ValueError, match="default_evalue"):
            _validate_hmm_block(self._cfg({"default_evalue": 0}), "Fam")

    def test_relative_length_out_of_range_raises(self):
        with pytest.raises(ValueError, match="relative_length_cutoff"):
            _validate_hmm_block(self._cfg({"relative_length_cutoff": 1.5}), "Fam")

    def test_negative_overlap_tolerance_raises(self):
        with pytest.raises(ValueError, match="multidomain_overlap_tolerance"):
            _validate_hmm_block(
                self._cfg({"multidomain_overlap_tolerance": -1}), "Fam")

    def test_bad_threads_raises(self):
        with pytest.raises(ValueError, match="hmm.threads"):
            _validate_hmm_block(self._cfg({"threads": 0}), "Fam")

    def test_malformed_token_in_sequence_raises(self):
        cfg = self._cfg(
            {"enabled": False},
            sequence={"region": "hexon", "hmms": ["A----B"]},
        )
        with pytest.raises(ValueError, match="sequence.hmms token"):
            _validate_hmm_block(cfg, "Fam")

    def test_malformed_token_in_concatenation_raises(self):
        cfg = self._cfg(
            {"enabled": False},
            concatenation={"proteins": [{"name": "m", "hmms": [" --A"]}]},
        )
        with pytest.raises(ValueError, match=r"concatenation.proteins\[0\].hmms"):
            _validate_hmm_block(cfg, "Fam")

    def test_tokens_validated_even_when_disabled(self):
        # Disabled but a typo'd token still surfaces (master switch independent).
        cfg = self._cfg(
            {"enabled": False},
            sequence={"region": "x", "hmms": [123]},
        )
        with pytest.raises(ValueError, match="must be strings"):
            _validate_hmm_block(cfg, "Fam")

    def test_well_formed_multidomain_token_ok(self):
        cfg = self._cfg(
            {"enabled": True, "database": "Adenoviridae"},
            sequence={"region": "hexon", "hmms": ["Adeno_hexon--Adeno_hexon_C"]},
        )
        _validate_hmm_block(cfg, "Fam")  # no raise


class TestHmmPresets:
    def test_adenoviridae_preset_hmm_enabled(self):
        cfg = _generate_default_family_config("Adenoviridae", make_minimal_global_cfg())
        _validate_hmm_block(cfg, "Adenoviridae")
        assert cfg["hmm"]["enabled"] is True
        assert cfg["hmm"]["database"] == "Adenoviridae"
        # Architecture token first, single-domain fallback OR'd in.
        assert cfg["sequence"]["hmms"] == [
            "Adeno_hexon--Adeno_hexon_C", "Adeno_hexon",
        ]

    def test_poxviridae_preset_tokens(self):
        cfg = _generate_default_family_config("Poxviridae", make_minimal_global_cfg())
        _validate_hmm_block(cfg, "Poxviridae")
        assert cfg["hmm"]["enabled"] is True
        toks = [p.get("hmms") for p in cfg["concatenation"]["proteins"]]
        assert ["DNA_pol_B"] in toks and ["UDG"] in toks
        assert all(t for t in toks)  # every marker is HMM-gated

    def test_non_hmm_family_defaults_off(self):
        # Smacoviridae: name-based Rep (DNA_FAMILIES), no HMM.
        cfg = _generate_default_family_config("Smacoviridae", make_minimal_global_cfg())
        _validate_hmm_block(cfg, "Smacoviridae")
        assert cfg["hmm"]["enabled"] is False
        assert cfg["hmm"]["database"] is None


# ---------------------------------------------------------------------------
# v1.5.0 inference defaults (accuracy review, 2026-10-20)
# ---------------------------------------------------------------------------

class TestInferenceDefaults:
    def test_tree_500_is_iqtree_fast_nuc_capped_search_aa(self):
        # --fast for nucleotide; a 50-iteration search for protein, where
        # --fast lands in a different local optimum on every seed.
        t5 = DEFAULT_FAMILY_CONFIG["tree_500"]
        assert t5["tool"] == "iqtree"
        assert t5["options_nuc"].split() == ["--fast"]
        assert t5["options_aa"].split() == ["-n", "50"]

    def test_tree_100_full_search_with_both_supports(self):
        t1 = DEFAULT_FAMILY_CONFIG["tree_100"]
        assert t1["tool"] == "iqtree"
        for key in ("options_nuc", "options_aa"):
            parts = t1[key].split()
            assert "--fast" not in parts
            assert "-alrt" in parts and "-B" in parts and "--bnni" in parts

    def test_nucleotide_msa_adjusts_direction(self):
        for block in ("msa_500", "msa_100"):
            nuc = DEFAULT_FAMILY_CONFIG[block]["options_nuc"].split()
            aa = DEFAULT_FAMILY_CONFIG[block]["options_aa"].split()
            assert "--adjustdirection" in nuc
            # MAFFT accepts --adjustdirection for nucleotide input only.
            assert "--adjustdirection" not in aa
        assert "--retree" in DEFAULT_FAMILY_CONFIG["msa_500"]["options_nuc"]

    def test_global_template_mirrors_family_defaults(self):
        from vfam_trees.cli import GLOBAL_CONFIG_TEMPLATE
        tmpl = yaml.safe_load(GLOBAL_CONFIG_TEMPLATE)["defaults"]
        for block in ("msa_500", "msa_100", "tree_500", "tree_100"):
            assert tmpl[block] == DEFAULT_FAMILY_CONFIG[block], block



# ---------------------------------------------------------------------------
# _validate_tree_blocks — tool and options must match (load-time check)
# ---------------------------------------------------------------------------

class TestValidateTreeBlocks:
    def _cfg(self, override, family="Testviridae"):
        return _merge_with_defaults(copy.deepcopy(override), make_minimal_global_cfg(), family)

    def test_defaults_pass(self):
        for fam in ("Xinmoviridae", "Circoviridae", "Poxviridae"):
            _validate_tree_blocks(self._cfg({}, fam), fam)

    def test_fasttree_with_template_options_rejected(self):
        # The template's options_nuc ("--fast") wins over anything else, so
        # switching only the tool used to crash FastTree per family.
        with pytest.raises(ValueError, match=r"tree_500.*--fast.*FastTree rejects"):
            _validate_tree_blocks(self._cfg({"tree_500": {"tool": "fasttree"}}), "Testviridae")

    def test_old_style_fasttree_block_rejected(self):
        cfg = self._cfg({"tree_500": {"tool": "fasttree", "options": ""}})
        with pytest.raises(ValueError, match="options_nuc"):
            _validate_tree_blocks(cfg, "Testviridae")

    def test_fasttree_with_cleared_options_passes(self):
        cfg = self._cfg({"tree_500": {"tool": "fasttree", "options_nuc": "", "options_aa": ""}})
        _validate_tree_blocks(cfg, "Testviridae")

    def test_fasttree_own_options_pass(self):
        cfg = self._cfg({"tree_500": {"tool": "fasttree",
                                      "options_nuc": "-spr 4 -mlacc 2 -slownni",
                                      "options_aa": ""}})
        _validate_tree_blocks(cfg, "Testviridae")

    def test_iqtree_with_fasttree_options_rejected(self):
        cfg = self._cfg({"tree_100": {"options_nuc": "-gtr -gamma"}})
        with pytest.raises(ValueError, match=r"tree_100.*-gtr -gamma.*IQ-TREE rejects"):
            _validate_tree_blocks(cfg, "Testviridae")

    def test_only_the_options_key_in_use_is_checked(self):
        # Nucleotide family (no preset): options_aa is never used, so it is not checked.
        cfg = self._cfg({"tree_500": {"tool": "fasttree", "options_nuc": ""}})
        _validate_tree_blocks(cfg, "Testviridae")
        # Protein family: options_aa is the one in use.
        with pytest.raises(ValueError, match="options_aa"):
            _validate_tree_blocks(
                self._cfg({"tree_500": {"tool": "fasttree", "options_nuc": ""}}, "Circoviridae"),
                "Circoviridae",
            )

    def test_unknown_tool_rejected(self):
        with pytest.raises(ValueError, match="must be 'iqtree' or 'fasttree'"):
            _validate_tree_blocks(self._cfg({"tree_100": {"tool": "raxml"}}), "Testviridae")

    def test_load_family_config_runs_the_check(self, tmp_path):
        (tmp_path / "Xinmoviridae.yaml").write_text("tree_500:\n  tool: fasttree\n")
        with pytest.raises(ValueError, match="FastTree rejects"):
            load_family_config("Xinmoviridae", tmp_path, make_minimal_global_cfg())

    def test_fasttree_with_iqtree_iteration_count_rejected(self):
        # FastTree's -n is the number of alignments in the input, not iterations.
        cfg = self._cfg({"tree_500": {"tool": "fasttree", "options_nuc": ""}}, "Circoviridae")
        with pytest.raises(ValueError, match=r"tree_500\.options_aa.*-n.*FastTree rejects or reads differently"):
            _validate_tree_blocks(cfg, "Circoviridae")

    def test_quick_check_defaults_on(self):
        cfg = self._cfg({})
        assert cfg["outlier_removal"]["quick_check"] is True
        _validate_tree_blocks(self._cfg({"outlier_removal": {"quick_check": False}}), "Testviridae")

    @pytest.mark.parametrize("bad", ["yes", 1, None])
    def test_quick_check_must_be_boolean(self, bad):
        cfg = self._cfg({"outlier_removal": {"quick_check": bad}})
        with pytest.raises(ValueError, match="outlier_removal.quick_check must be true or false"):
            _validate_tree_blocks(cfg, "Testviridae")


# ---------------------------------------------------------------------------
# RDRP_FAMILIES — RNA-virus families on the HMM-identified RdRp marker
# ---------------------------------------------------------------------------

class TestRdrpPreset:
    def test_disjoint_from_other_presets(self):
        assert not RDRP_FAMILIES & set(DNA_FAMILIES)
        assert not RDRP_FAMILIES & set(CONCATENATION_FAMILIES)

    def test_tokens_are_the_bundled_profiles(self, tmp_path):
        from vfam_trees.hmm.database import parse_ga_cutoffs, resolve_database_path
        db = resolve_database_path("RdRp", cache_dir=tmp_path)
        names = set(parse_ga_cutoffs(db))
        assert set(RDRP_HMM_TOKENS) == names
        assert len(names) == 23

    def test_every_rdrp_family_config_validates(self):
        for fam in sorted(RDRP_FAMILIES):
            cfg = _merge_with_defaults({}, make_minimal_global_cfg(), fam)
            _validate_hmm_block(cfg, fam)
            _validate_tree_blocks(cfg, fam)
            assert cfg["sequence"]["type"] == "protein", fam
            assert cfg["hmm"]["use_ga_when_available"] is False
            assert cfg["hmm"]["default_evalue"] == pytest.approx(1e-3)
            assert cfg["hmm"]["relative_length_cutoff"] == pytest.approx(0.2)
            assert cfg["hmm"]["best_hit_only"] is True
            from vfam_trees.config import RDRP_NON_CORE_DOMAINS
            tokens = [p for p in RDRP_FAMILY_PROFILES[fam] if p not in RDRP_NON_CORE_DOMAINS]
            assert tokens, fam  # every family keeps a core RdRp profile
            assert cfg["sequence"]["hmms"] == tokens, fam
            # The non-core L domains still count as the family's own for best_hit_only.
            assert set(RDRP_FAMILY_PROFILES[fam]) - set(tokens) <= set(cfg["hmm"]["best_hit_also"])

    def test_non_core_l_domains_never_select(self):
        """L_PA-C-like / L_thumb_ring_vir select and crop nothing: a partial L
        covering only such a domain would enter the alignment cut to it."""
        from vfam_trees.config import RDRP_NON_CORE_DOMAINS, rdrp_preset
        assert RDRP_NON_CORE_DOMAINS == {"L_PA-C-like", "L_thumb_ring_vir"}
        for fam, tokens in (("Peribunyaviridae", ["Bunya_RdRp", "Mononeg_RNA_pol"]),
                            ("Phenuiviridae", ["Bunya_RdRp", "Flu_PB1"]),
                            ("Cruliviridae", ["Bunya_RdRp"])):
            preset = rdrp_preset(fam)
            assert preset["sequence"]["hmms"] == tokens, fam
            assert preset["hmm"]["best_hit_also"] == ["L_PA-C-like", "L_thumb_ring_vir"], fam

    def test_family_profiles_are_bundled_profiles(self):
        assert set(RDRP_FAMILY_PROFILES) == RDRP_FAMILIES
        used = set()
        for fam, profiles in RDRP_FAMILY_PROFILES.items():
            assert profiles, fam
            assert len(set(profiles)) == len(profiles), fam
            assert set(profiles) <= set(RDRP_HMM_TOKENS), fam
            used.update(profiles)
        # Every bundled profile serves some family, so a foreign RdRp is
        # recognised by its best hit whatever family it was filed under.
        assert used == set(RDRP_HMM_TOKENS)

    def test_palm_profiles_come_as_a_group(self):
        for fam, profiles in RDRP_FAMILY_PROFILES.items():
            palm = set(profiles) & set(RDRP_PALM)
            assert palm in (set(), set(RDRP_PALM)), fam

    def test_foreign_kinds_excluded(self):
        # The misfiled orbivirus / bunyavirus records seen under Alphavirus.
        assert "Orbi_VP1" not in RDRP_FAMILY_PROFILES["Togaviridae"]
        assert "Bunya_RdRp" not in RDRP_FAMILY_PROFILES["Togaviridae"]
        assert RDRP_FAMILY_PROFILES["Rhabdoviridae"] == ("Mononeg_RNA_pol",)

    def test_best_hit_only_off_by_default_and_validated(self):
        cfg = _merge_with_defaults({}, make_minimal_global_cfg(), "Testviridae")
        _validate_hmm_block(cfg, "Testviridae")
        assert cfg["hmm"]["best_hit_only"] is False
        cfg["hmm"]["best_hit_only"] = "yes"
        with pytest.raises(ValueError, match="best_hit_only"):
            _validate_hmm_block(cfg, "Testviridae")

    def test_user_file_wins_and_is_warned_about(self, tmp_path, caplog):
        import logging
        (tmp_path / "Hantaviridae.yaml").write_text(
            "sequence:\n  type: nucleotide\n  region: whole_genome\n  segment: segment L\n"
            "hmm:\n  enabled: false\n"
        )
        with caplog.at_level(logging.WARNING):
            cfg, auto = load_family_config("Hantaviridae", tmp_path, make_minimal_global_cfg())
        assert auto is False
        assert cfg["sequence"]["type"] == "nucleotide" and cfg["hmm"]["enabled"] is False
        assert cfg["sequence"]["segment"] == "segment L"
        assert any("RdRp" in r.getMessage() for r in caplog.records)

    def test_domain_padding_validated(self):
        cfg = _merge_with_defaults({}, make_minimal_global_cfg(), "Flaviviridae")
        cfg["sequence"]["domain_padding"] = -5
        with pytest.raises(ValueError, match="domain_padding"):
            _validate_hmm_block(cfg, "Flaviviridae")
        cfg["sequence"]["domain_padding"] = None
        _validate_hmm_block(cfg, "Flaviviridae")


class TestUnknownKeys:
    def test_nested(self):
        from vfam_trees.config import unknown_config_keys, DEFAULT_FAMILY_CONFIG
        cfg = {"length_outlier": {"k": 5.0, "hi_mult": 3.0}, "tree_100": {"tool": "iqtree"},
               "bogus": 1, "_family": "X",
               "concatenation": {"proteins": [{"name": "A", "anything": 1}]}}
        assert unknown_config_keys(cfg, DEFAULT_FAMILY_CONFIG) == [
            "length_outlier.hi_mult", "bogus"]

    def test_family_file_nested_key_warns_and_is_kept(self, tmp_path, caplog):
        import logging
        (tmp_path / "Geminiviridae.yaml").write_text(
            "length_outlier:\n  hi_mult: 3.0\n  lo_mult: 0.333\n")
        with caplog.at_level(logging.WARNING):
            cfg, _ = load_family_config("Geminiviridae", tmp_path, make_minimal_global_cfg())
        assert cfg["_unknown_keys"] == ["length_outlier.hi_mult", "length_outlier.lo_mult"]
        assert "length_outlier.hi_mult, length_outlier.lo_mult" in caplog.text
        assert cfg["length_outlier"]["k"] == 5.0  # the real settings still apply

    def test_generated_files_of_every_preset_family_load_cleanly(self, tmp_path, caplog):
        """Presets may only set keys the template knows (sequence.hmms was
        missing), or every generated family file would warn on reload."""
        import logging
        from vfam_trees import config as c
        fams = (set(c.DNA_FAMILIES) | set(c.CONCATENATION_FAMILIES)
                | set(c.SEGMENTED_FAMILIES) | set(c.HMM_MARKER_FAMILIES))
        g = make_minimal_global_cfg()
        for fam in sorted(fams):
            load_family_config(fam, tmp_path, g)  # generates the file
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            for fam in sorted(fams):
                cfg, auto = load_family_config(fam, tmp_path, g)
                assert not auto and cfg["_unknown_keys"] == [], fam
        assert "Unknown config key" not in caplog.text

    def test_global_retired_keys_explained(self, tmp_path, caplog):
        import logging
        from vfam_trees.config import load_global_config
        p = tmp_path / "global.yaml"
        p.write_text('ncbi:\n  email: "a@b.org"\noutput_dir: results\nlog_level: INFO\n'
                     "icon:\n  colour: red\ndefaults:\n  quality:\n    min_len: 5\n")
        with caplog.at_level(logging.WARNING):
            cfg = load_global_config(p)
        assert cfg["_unknown_keys"] == [
            "output_dir", "log_level", "icon.colour", "defaults.quality.min_len"]
        assert "output_dir (not read: the output folder is set with `vfam_trees run -o`)" in caplog.text
        assert "log_level (not read: the log level is set with" in caplog.text

    def test_clean_global_has_no_unknown_keys(self, tmp_path, caplog):
        import logging
        from vfam_trees.config import load_global_config
        p = tmp_path / "global.yaml"
        p.write_text('ncbi:\n  email: "a@b.org"\n  api_key: ""\nannotation_tsv: null\n'
                     "cache:\n  dir: ~/.vfam_cache\n  ttl_days: 90\n")
        with caplog.at_level(logging.WARNING):
            assert load_global_config(p)["_unknown_keys"] == []
        assert "Unknown config key" not in caplog.text


class TestGlobalSettings:
    def test_fills_missing_keys(self):
        from vfam_trees.config import global_settings
        assert global_settings({}, "icon") == {
            "size": 256, "bg_color": "#EAF3F2", "branch_color": "#000000"}
        assert global_settings({"icon": {"size": 64}}, "icon")["size"] == 64
        assert global_settings({"cache": None}, "cache") == {"dir": None, "ttl_days": None}

    def test_effective_global_config(self):
        from vfam_trees.config import effective_global_config, DEFAULT_FAMILY_CONFIG
        g = {"ncbi": {"email": "a@b.org", "api_key": "K"}, "_unknown_keys": [],
             "defaults": {"targets": {"max_500": 300}, "length_outlier": {"hi_mult": 3.0}}}
        eff = effective_global_config(g)
        assert eff["ncbi"] == {"email": "a@b.org", "api_key": "<set, not recorded>"}
        assert eff["visualization"] == {"branch_linewidth": 0.5}
        assert eff["defaults"]["targets"] == {"max_500": 300, "max_100": 100}
        assert eff["defaults"]["length_outlier"]["hi_mult"] == 3.0  # kept as given
        assert eff["defaults"]["hmm"] == DEFAULT_FAMILY_CONFIG["hmm"]
        assert "_unknown_keys" not in eff
        assert g["ncbi"]["api_key"] == "K"  # caller's dict untouched
        assert effective_global_config({"ncbi": {"api_key": ""}})["ncbi"]["api_key"] == ""


class TestHmmSequenceMode:
    """Code review 2026-10-21: partly edited RdRp files must not run as a mix
    of nucleotide settings and the HMM protein fetch."""

    def _generated(self, tmp_path, fam):
        load_family_config(fam, tmp_path, make_minimal_global_cfg())
        return yaml.safe_load((tmp_path / f"{fam}.yaml").read_text())

    def test_only_type_nucleotide_is_refused_with_instructions(self, tmp_path):
        cfg = self._generated(tmp_path, "Filoviridae")
        cfg["sequence"]["type"] = "nucleotide"
        (tmp_path / "Filoviridae.yaml").write_text(yaml.safe_dump(cfg))
        with pytest.raises(ValueError) as exc:
            load_family_config("Filoviridae", tmp_path, make_minimal_global_cfg())
        msg = str(exc.value)
        assert "sequence.region: whole_genome, hmm.enabled: false" in msg
        assert "Filoviridae.yaml" in msg

    def test_segmented_family_instructions_name_the_segment(self, tmp_path):
        (tmp_path / "Hantaviridae.yaml").write_text(
            "sequence:\n  type: nucleotide\nhmm:\n  enabled: true\n")
        with pytest.raises(ValueError, match="sequence.segment: 'segment L'"):
            load_family_config("Hantaviridae", tmp_path, make_minimal_global_cfg())

    def test_segment_with_hmm_marker_is_refused(self, tmp_path):
        (tmp_path / "Hantaviridae.yaml").write_text("sequence:\n  segment: segment M\n")
        with pytest.raises(ValueError, match="Protein titles name no segment"):
            load_family_config("Hantaviridae", tmp_path, make_minimal_global_cfg())

    def test_full_opt_out_warning_states_the_real_segment(self, tmp_path, caplog):
        import logging
        cfg = self._generated(tmp_path, "Hantaviridae")
        cfg["sequence"].update(type="nucleotide", region="whole_genome")
        cfg["hmm"]["enabled"] = False
        (tmp_path / "Hantaviridae.yaml").write_text(yaml.safe_dump(cfg))
        with caplog.at_level(logging.WARNING):
            merged, _ = load_family_config("Hantaviridae", tmp_path, make_minimal_global_cfg())
        assert merged["sequence"]["segment"] is None       # the file's null wins
        assert "no segment: all of its segments go into one alignment" in caplog.text
        assert "add sequence.segment: 'segment L'" in caplog.text
        assert "segment 'segment L'" not in caplog.text.replace("sequence.segment: 'segment L'", "")

    def test_full_opt_out_with_segment(self, tmp_path, caplog):
        import logging
        (tmp_path / "Hantaviridae.yaml").write_text(
            "sequence:\n  type: nucleotide\n  region: whole_genome\n  segment: segment L\n"
            "hmm:\n  enabled: false\n")
        with caplog.at_level(logging.WARNING):
            merged, _ = load_family_config("Hantaviridae", tmp_path, make_minimal_global_cfg())
        assert merged["sequence"]["segment"] == "segment L"
        assert "region whole_genome, segment 'segment L')" in caplog.text

    def test_generated_configs_never_trip_the_check(self, tmp_path):
        from vfam_trees import config as c
        fams = (set(c.DNA_FAMILIES) | set(c.CONCATENATION_FAMILIES)
                | set(c.SEGMENTED_FAMILIES) | set(c.HMM_MARKER_FAMILIES))
        for fam in sorted(fams):
            c._validate_sequence_mode(
                c._generate_default_family_config(fam, make_minimal_global_cfg()), fam)

    def test_preview_writes_nothing(self, tmp_path):
        from vfam_trees.config import preview_family_config
        cfg = preview_family_config("Togaviridae", tmp_path, make_minimal_global_cfg())
        assert cfg["hmm"]["enabled"] is True
        assert not list(tmp_path.iterdir())


# ---------------------------------------------------------------------------
# RT and Rep markers (v1.6.0): reverse-transcribing and CRESS-DNA families
# ---------------------------------------------------------------------------

def _bundled_names(marker):
    from vfam_trees.hmm import BUNDLED_HMMS_DIR
    names = set()
    for f in (BUNDLED_HMMS_DIR / marker).glob("*.hmm"):
        names |= {l.split()[1] for l in f.read_text().splitlines() if l.startswith("NAME")}
    return names


class TestRtRepMarkers:
    def test_bundled_sets_match_the_token_lists(self):
        from vfam_trees.config import RT_HMM_TOKENS, REP_HMM_TOKENS
        assert _bundled_names("RT") == set(RT_HMM_TOKENS)
        assert _bundled_names("Rep") == set(REP_HMM_TOKENS)

    def test_marker_families_are_disjoint_from_other_presets(self):
        from vfam_trees.config import (
            HMM_MARKER_FAMILIES, RT_FAMILY_PROFILES, REP_FAMILIES)
        assert not set(HMM_MARKER_FAMILIES) & set(DNA_FAMILIES)
        assert not set(HMM_MARKER_FAMILIES) & set(CONCATENATION_FAMILIES)
        assert not RDRP_FAMILIES & set(RT_FAMILY_PROFILES)
        assert not RDRP_FAMILIES & REP_FAMILIES
        assert not set(RT_FAMILY_PROFILES) & REP_FAMILIES
        assert len(HMM_MARKER_FAMILIES) == len(RDRP_FAMILIES) + 4 + 23

    def test_rt_preset(self):
        cfg = _merge_with_defaults({}, make_minimal_global_cfg(), "Retroviridae")
        seq = cfg["sequence"]
        assert (seq["type"], seq["region"], seq["segment"]) == ("protein", "RT", None)
        assert seq["hmms"] == ["RVT_1"] and seq["domain_padding"] == 100
        assert cfg["hmm"]["database"] == "RT" and cfg["hmm"]["enabled"] is True
        assert cfg["hmm"]["best_hit_only"] is True
        assert cfg["hmm"]["best_hit_also"] == ["RVT_thumb", "RVT_connect"]
        assert cfg["hmm"]["use_ga_when_available"] is False

    def test_rep_preset_keeps_whole_proteins(self):
        for fam in ("Geminiviridae", "Circoviridae", "Alphasatellitidae"):
            cfg = _merge_with_defaults({}, make_minimal_global_cfg(), fam)
            assert cfg["sequence"]["type"] == "protein" and cfg["sequence"]["region"] == "Rep", fam
            assert cfg["sequence"]["hmms"] == ["Viral_Rep", "RNA_helicase", "Gemini_AL1", "Gemini_AL1_M"]
            assert cfg["sequence"]["domain_padding"] is None
            assert cfg["hmm"]["database"] == "Rep" and cfg["hmm"]["best_hit_only"] is False

    def test_moved_from_dna_presets(self):
        # Circoviridae was a name-based Rep, Hepadnaviridae whole genome.
        assert _merge_with_defaults({}, MINIMAL_GLOBAL, "Hepadnaviridae")["sequence"]["region"] == "RT"
        assert _merge_with_defaults({}, MINIMAL_GLOBAL, "Circoviridae")["hmm"]["enabled"] is True
        smaco = _merge_with_defaults({}, MINIMAL_GLOBAL, "Smacoviridae")
        assert smaco["sequence"]["region"] == "Rep" and smaco["hmm"]["enabled"] is False

    def test_every_marker_family_config_validates(self):
        from vfam_trees.config import HMM_MARKER_FAMILIES, _validate_sequence_mode
        for fam, marker in sorted(HMM_MARKER_FAMILIES.items()):
            cfg = _merge_with_defaults({}, make_minimal_global_cfg(), fam)
            _validate_hmm_block(cfg, fam)
            _validate_sequence_mode(cfg, fam)
            assert cfg["sequence"]["region"] == marker == cfg["hmm"]["database"], fam
            used = set(cfg["sequence"]["hmms"]) | set(cfg["hmm"]["best_hit_also"])
            assert used <= _bundled_names(marker), fam

    def test_opt_out(self):
        from vfam_trees.config import hmm_marker_opt_out
        assert hmm_marker_opt_out({"sequence": {"region": "RdRp"}}, "RT")
        assert not hmm_marker_opt_out({"sequence": {"region": "RT"}}, "RT")
        cfg = _merge_with_defaults(
            {"sequence": {"type": "nucleotide", "region": "whole_genome"}, "hmm": {"enabled": False}},
            MINIMAL_GLOBAL, "Geminiviridae")
        assert cfg["sequence"]["type"] == "nucleotide" and cfg["sequence"]["hmms"] == []

    def test_bundled_databases_resolve(self, tmp_path):
        from vfam_trees.hmm import resolve_database_path, parse_ga_cutoffs
        for marker, n in (("RT", 4), ("Rep", 4)):
            db = resolve_database_path(marker, cache_dir=tmp_path)
            assert len(parse_ga_cutoffs(db)) == n


class TestMamonoviridaePreset:
    """Medusaviruses: family-B DNA polymerase by HMM, cut to the domain."""

    def test_auto_generated_config(self, tmp_path):
        cfg, auto = load_family_config("Mamonoviridae", tmp_path, MINIMAL_GLOBAL)
        assert auto is True
        seq = cfg["sequence"]
        assert (seq["type"], seq["region"], seq["hmms"], seq["domain_padding"]) == \
            ("protein", "DNA polymerase", ["DNA_pol_B"], 100)
        assert cfg["hmm"]["enabled"] is True and cfg["hmm"]["database"] == "Mamonoviridae"
        assert "concatenation" not in cfg or not cfg["concatenation"].get("proteins")

    def test_bundled_profile_present(self):
        import vfam_trees
        hmm = Path(vfam_trees.__file__).parent / "data" / "hmms" / "Mamonoviridae" / "DNA_pol_B.hmm"
        text = hmm.read_text()
        assert "NAME  DNA_pol_B" in text and "ACC   PF00136" in text


class TestNcldvHmmConcat:
    """Asfarviridae and Pithoviridae: concatenations of HMM-identified NCLDV
    core markers (bundled set "NCLDV"), not GenBank names."""

    def test_presets(self, tmp_path):
        for fam, n, has_mcp in (("Asfarviridae", 8, True), ("Pithoviridae", 7, False)):
            cfg, _ = load_family_config(fam, tmp_path, MINIMAL_GLOBAL)
            assert cfg["sequence"]["region"] == "concatenated", fam
            assert cfg["hmm"]["enabled"] is True and cfg["hmm"]["database"] == "NCLDV", fam
            assert cfg["hmm"]["relative_length_cutoff"] == 0.30, fam
            proteins = cfg["concatenation"]["proteins"]
            assert len(proteins) == n and all(p["hmms"] for p in proteins), fam
            assert any(p["hmms"] == ["Capsid_NCLDV"] for p in proteins) is has_mcp, fam
            _validate_hmm_block(cfg, fam)

    def test_every_token_is_a_bundled_profile(self):
        tokens = {tok for fam in ("Asfarviridae", "Pithoviridae")
                  for p in CONCATENATION_FAMILIES[fam]["concatenation"]["proteins"]
                  for tok in p["hmms"]}
        assert tokens == _bundled_names("NCLDV")

    def test_bundled_database_resolves(self, tmp_path):
        from vfam_trees.hmm import resolve_database_path, parse_ga_cutoffs
        db = resolve_database_path("NCLDV", cache_dir=tmp_path)
        assert len(parse_ga_cutoffs(db)) == 9
