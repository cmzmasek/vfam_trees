"""Config loading, validation, and auto-generation for vfam_trees."""
from __future__ import annotations

import copy
from datetime import date
from pathlib import Path

import yaml

from .logger import get_logger

log = get_logger(__name__)


# Segmented families and their default phylogenetic segment.
# The segment string is used as a title keyword in the NCBI query
# (e.g. "segment L" → AND "segment L"[Title]).
# Override per-family in the YAML config if needed.
SEGMENTED_FAMILIES: dict[str, str] = {
    # 2-segment negative-sense / ambisense
    "Arenaviridae": "segment L",
    # 2-segment dsRNA
    "Birnaviridae": "segment B",
    "Picobirnaviridae": "segment 1",
    "Partitiviridae": "RNA1",
    "Amalgaviridae": "RNA1",
    "Megabirnaviridae": "RNA1",
    # 3-segment negative-sense (Bunyavirales)
    "Hantaviridae": "segment L",
    "Nairoviridae": "segment L",
    "Peribunyaviridae": "segment L",
    "Phenuiviridae": "segment L",
    "Tospoviridae": "segment L",
    "Fimoviridae": "RNA1",
    "Leishbuviridae": "segment L",
    "Phasmaviridae": "segment L",
    "Wupedeviridae": "segment L",
    "Discoviridae": "segment L",
    "Cruliviridae": "segment L",
    "Blumeviridae": "segment L",
    "Konkoviridae": "segment L",
    "Mypoviridae": "segment L",
    "Tulasviridae": "segment L",
    "Steitzviridae": "segment L",
    # 4-segment dsRNA
    "Chrysoviridae": "RNA1",
    "Quadriviridae": "RNA1",
    # 6-8-segment negative-sense
    "Orthomyxoviridae": "PB1",
    # 9-12-segment dsRNA (Reovirales)
    "Reoviridae": "segment 1",
    "Sedoreoviridae": "segment 1",
    "Spinareoviridae": "segment 1",
}


# DNA virus families and their optimal search strategy.
# "region" is used with "<region>"[Protein Name] OR "<region>"[Gene] for
# protein markers, or a title-based "complete genome" search when set to
# whole_genome.  Nearly all families use a conserved marker protein — the
# only whole-genome default remaining is Hepadnaviridae (~3.2 kb dsDNA),
# whose overlapping-ORF architecture makes whole-genome nucleotide
# alignment workable at family scale.
DNA_FAMILIES: dict[str, dict] = {
    # ---- Small ssDNA / dsDNA families: marker-protein approach ----
    # CRESS-DNA families — ambisense organisation (Rep and Cap on opposite
    # strands around the origin of replication) makes whole-genome
    # alignment unreliable; Rep is the ICTV species-demarcation marker.
    "Circoviridae":       {"sequence": {"region": "Rep", "type": "protein"}},          # ~2 kb ssDNA circular
    "Smacoviridae":       {"sequence": {"region": "Rep", "type": "protein"}},          # ~2.3 kb ssDNA circular
    # Anelloviridae — ORF1 (capsid) is the standard phylogenetic marker;
    # ORF2/ORF3 are hypervariable and make whole-genome alignment poor.
    "Anelloviridae":      {"sequence": {"region": "ORF1", "type": "protein"}},         # ~2–4 kb ssDNA circular
    # Hepadnaviridae — small enough (3.2 kb) with heavily overlapping ORFs
    # that whole-genome nucleotide alignment still works at family scale.
    "Hepadnaviridae":     {"sequence": {"region": "whole_genome"}},                    # ~3.2 kb partial dsDNA
    # Parvoviridae — NS1 (replicase) is ICTV's recommended phylogenetic
    # marker; conserved enough to align across Parvovirinae/Densovirinae.
    "Parvoviridae":       {"sequence": {"region": "NS1", "type": "protein"}},          # ~5 kb ssDNA linear
    # Polyomaviridae — large T antigen is the standard cross-genus marker
    # (early/late ambisense organisation rules out whole-genome alignment).
    "Polyomaviridae":     {"sequence": {"region": "large T antigen", "type": "protein"}},  # ~5 kb dsDNA circular
    # Papillomaviridae — L1 (major capsid) is the ICTV gold-standard
    # marker: species demarcation = >10% L1 nucleotide divergence.
    "Papillomaviridae":   {"sequence": {"region": "L1", "type": "protein"}},           # ~8 kb dsDNA circular

    # ---- Medium–large dsDNA families: marker-protein approach ----
    # Adenoviridae (~35 kb) — hexon (major capsid protein)
    "Adenoviridae":       {"sequence": {"region": "hexon", "type": "protein"}},
    # Herpesviruses (~130–240 kb) — DNA polymerase (UL30/UL30-like)
    "Orthoherpesviridae": {"sequence": {"region": "DNA polymerase", "type": "protein"}},
    "Herpesviridae":      {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # legacy ICTV name
    "Alloherpesviridae":  {"sequence": {"region": "DNA polymerase", "type": "protein"}},
    "Malacoherpesviridae":{"sequence": {"region": "DNA polymerase", "type": "protein"}},
    # Poxviridae (~130–375 kb) — largest DNA-directed RNA polymerase subunit
    # (rpo147 / A24R homolog).  Single-marker cross-subfamily resolution is
    # inherently limited; concatenated core genes (~20–30) are the real
    # fix and sit outside this pipeline's single-gene model.  Chordopoxvirinae
    # are typically annotated gene=rpo147; some Entomopoxvirinae annotations
    # differ and may need an override in the per-family YAML.
    "Poxviridae":         {"sequence": {"region": "rpo147", "type": "protein"}},
    # Iridoviridae (~100–220 kb) — major capsid protein
    "Iridoviridae":       {"sequence": {"region": "major capsid protein", "type": "protein"}},
    # Asfarviridae (~190 kb) — B646L (p72, major capsid protein)
    "Asfarviridae":       {"sequence": {"region": "B646L", "type": "protein"}},
    # Nimaviridae (~300 kb, shrimp WSSV) — DNA polymerase
    "Nimaviridae":        {"sequence": {"region": "DNA polymerase", "type": "protein"}},
    # Hytrosaviridae (~120–190 kb, insect salivary gland hypertrophy viruses)
    "Hytrosaviridae":     {"sequence": {"region": "DNA polymerase", "type": "protein"}},
    # Baculoviridae / Nudiviridae / Ascoviridae (~80–230 kb) — lef-8 (late
    # expression factor 8) is a member of the ICTV 3-gene species-demarcation
    # set (polh, lef-8, lef-9) with the strongest single-marker phylogenetic
    # signal for these insect dsDNA families.
    "Baculoviridae":      {"sequence": {"region": "lef-8", "type": "protein"}},
    "Nudiviridae":        {"sequence": {"region": "lef-8", "type": "protein"}},
    "Ascoviridae":        {"sequence": {"region": "lef-8", "type": "protein"}},

    # ---- Nucleocytoplasmic large DNA viruses (NCLDVs) / giant dsDNA ----
    # DNA polymerase B is the universal NCLDV phylogenetic marker.
    "Phycodnaviridae":    {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~150–560 kb algal
    "Mimiviridae":        {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~1.2 Mb amoebal
    "Marseilleviridae":   {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~350 kb amoebal
    "Pandoraviridae":     {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~2 Mb amoebal
    "Pithoviridae":       {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~600 kb amoebal
    "Medusaviridae":      {"sequence": {"region": "DNA polymerase", "type": "protein"}},  # ~380 kb amoebal
}


DEFAULT_FAMILY_CONFIG: dict = {
    "download": {
        "max_per_species": 300,
    },
    "sequence": {
        "type": "nucleotide",
        "region": "whole_genome",
        "segment": None,   # set automatically for known segmented families
    },
    "quality": {
        "min_length": None,
        "max_ambiguous": 0.01,
        "exclude_organisms": [
            "synthetic construct",
            "metagenome",
            "MAG",
            "uncultured",
            "unverified",
            "vector",
            "recombinant",
            "patent",
        ],
    },
    "clustering": {
        "tool": "mmseqs2",
        "threshold_min": 0.70,
        "threshold_max": 0.99,
        "max_reps_500": 20,
        "max_reps_100": 5,
    },
    "targets": {
        "max_500": 500,
        "max_100": 100,
    },
    "msa_500": {
        "tool": "mafft",
        "options_nuc": "--6merpair --retree 1",
        "options_aa": "--6merpair --retree 1",
    },
    "msa_100": {
        "tool": "mafft",
        "options_nuc": "--retree 2",
        "options_aa": "--auto",
    },
    "msa_trim": {
        "enabled": True,
        "tool": "trimal",
        "options": "-automated1",
    },
    "tree_500": {
        "tool": "fasttree",
        "options": "",
        "model_nuc": "GTR+G",
        "model_aa": "LG+G",
    },
    "tree_100": {
        "tool": "iqtree",
        # Nucleotide: --fast (SH-aLRT branch support auto-added by the
        # IQ-TREE wrapper). Protein: UFBoot (-B 1000) produces stronger
        # support on divergent protein families; --fast is incompatible
        # with -B so a separate option is required.
        "options_nuc": "--fast",
        "options_aa": "-B 1000",
        "model_nuc": "GTR+G",
        "model_aa": "TEST",
    },
    "length_outlier": {
        "enabled": True,
        "hi_mult": 3.0,
        "lo_mult": 0.333,
    },
    "outlier_removal": {
        "enabled": True,
        "factor": 20.0,
        "max_iterations": 3,
        "min_seqs": 40,
    },
}


def load_global_config(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Global config not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    _validate_global_config(cfg)
    return cfg


def _validate_global_config(cfg: dict) -> None:
    import re
    email = cfg.get("ncbi", {}).get("email", "").strip()
    if not email:
        raise ValueError(
            "global.yaml: ncbi.email must be set. "
            "Run 'vfam_trees init' to generate a template config."
        )
    if email == "your.email@example.com":
        raise ValueError(
            "global.yaml: ncbi.email is still set to the placeholder value. "
            "Replace it with your actual email address."
        )
    if not re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", email):
        raise ValueError(
            f"global.yaml: ncbi.email does not look like a valid email address: '{email}'"
        )


_KNOWN_FAMILY_CONFIG_KEYS = frozenset(DEFAULT_FAMILY_CONFIG.keys())


def _warn_unknown_keys(cfg: dict, path: Path) -> None:
    """Warn about top-level keys in a user config that are not recognised."""
    unknown = [
        k for k in cfg
        if not k.startswith("_") and k not in _KNOWN_FAMILY_CONFIG_KEYS
    ]
    if unknown:
        log.warning(
            "Unknown config key(s) in %s: %s — these will be ignored. "
            "Known keys: %s",
            path,
            ", ".join(sorted(unknown)),
            ", ".join(sorted(_KNOWN_FAMILY_CONFIG_KEYS)),
        )


def load_family_config(family: str, configs_dir: Path, global_cfg: dict) -> tuple[dict, bool]:
    """Load per-family config, auto-generating it if missing.

    Returns (config_dict, was_auto_generated).
    """
    config_path = configs_dir / f"{family}.yaml"
    if config_path.exists():
        with open(config_path) as f:
            file_cfg = yaml.safe_load(f) or {}
        _warn_unknown_keys(file_cfg, config_path)
        _warn_smart_default_conflicts(family, file_cfg, config_path)
        cfg = _merge_with_defaults(file_cfg, global_cfg, family)
        return cfg, False
    else:
        cfg = _generate_default_family_config(family, global_cfg)
        _write_family_config(cfg, config_path)
        log.warning(
            "No config found for %s — auto-generated default at %s. "
            "Edit it to tune parameters for this family.",
            family,
            config_path,
        )
        return cfg, True


def _warn_smart_default_conflicts(family: str, file_cfg: dict, config_path: Path) -> None:
    """Warn when a config file overrides a DNA_FAMILIES or SEGMENTED_FAMILIES setting."""
    dna_overrides = DNA_FAMILIES.get(family)
    if dna_overrides:
        expected_region = (dna_overrides.get("sequence") or {}).get("region")
        file_region = (file_cfg.get("sequence") or {}).get("region")
        if expected_region and file_region and file_region != expected_region:
            log.warning(
                "%s: %s sets sequence.region=%r, overriding the recommended "
                "value %r for this family. If this is unintentional (e.g. a "
                "stale auto-generated file), delete %s and re-run to regenerate.",
                family, config_path.name, file_region, expected_region, config_path.name,
            )

    expected_segment = SEGMENTED_FAMILIES.get(family)
    if expected_segment:
        file_segment = (file_cfg.get("sequence") or {}).get("segment")
        if file_segment and file_segment != expected_segment:
            log.warning(
                "%s: %s sets sequence.segment=%r, overriding the recommended "
                "value %r for this family.",
                family, config_path.name, file_segment, expected_segment,
            )


def _merge_with_defaults(cfg: dict, global_cfg: dict, family: str = "") -> dict:
    """Fill missing keys from global defaults.

    Priority (lowest → highest):
      DEFAULT_FAMILY_CONFIG → global_defaults → DNA_FAMILIES/SEGMENTED_FAMILIES → cfg
    """
    global_defaults = global_cfg.get("defaults", {})
    merged = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _deep_update(merged, global_defaults)
    _apply_smart_defaults(family, merged)
    _deep_update(merged, cfg)
    return merged


def _apply_smart_defaults(family: str, cfg: dict) -> None:
    """Apply DNA_FAMILIES and SEGMENTED_FAMILIES as smart defaults in-place."""
    segment = SEGMENTED_FAMILIES.get(family)
    if segment and not cfg["sequence"].get("segment"):
        cfg["sequence"]["segment"] = segment
        log.info("Auto-configured segment '%s' for segmented family %s", segment, family)

    dna_overrides = DNA_FAMILIES.get(family)
    if dna_overrides:
        _deep_update(cfg, dna_overrides)
        region = cfg["sequence"].get("region", "whole_genome")
        if region == "whole_genome":
            log.info("Auto-configured DNA family %s: whole-genome search", family)
        else:
            log.info("Auto-configured DNA family %s: marker gene '%s'", family, region)


def _generate_default_family_config(family: str, global_cfg: dict) -> dict:
    global_defaults = global_cfg.get("defaults", {})
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _deep_update(cfg, global_defaults)
    _apply_smart_defaults(family, cfg)
    cfg["_family"] = family
    cfg["_generated"] = str(date.today())
    cfg["_auto_generated"] = True
    return cfg


def _write_family_config(cfg: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        f"# vfam_trees auto-generated config for {cfg.get('_family', 'unknown')}\n"
        f"# Generated: {cfg.get('_generated', '')}\n"
        "# Review and edit clustering thresholds and other parameters as needed.\n\n"
    )
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    with open(path, "w") as f:
        f.write(header)
        yaml.dump(clean, f, default_flow_style=False, sort_keys=False)


def _deep_update(base: dict, override: dict) -> None:
    for k, v in override.items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            _deep_update(base[k], v)
        else:
            base[k] = v
