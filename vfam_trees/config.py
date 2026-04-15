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
# "region" sets the NCBI [Gene] filter; whole_genome uses the title-based
# complete-genome search instead.  Families with large genomes (>~30 kb) use
# a phylogenetically informative marker gene; small-genome families use the
# full genome.
DNA_FAMILIES: dict[str, dict] = {
    # ---- Small genomes: whole-genome approach ----
    "Anelloviridae":      {"sequence": {"region": "whole_genome"}},  # ~2 kb ssDNA circular
    "Circoviridae":       {"sequence": {"region": "whole_genome"}},  # ~2 kb ssDNA circular
    "Smacoviridae":       {"sequence": {"region": "whole_genome"}},  # ~2 kb ssDNA circular
    "Hepadnaviridae":     {"sequence": {"region": "whole_genome"}},  # ~3.2 kb partial dsDNA
    "Parvoviridae":       {"sequence": {"region": "whole_genome"}},  # ~5 kb ssDNA linear
    "Polyomaviridae":     {"sequence": {"region": "whole_genome"}},  # ~5 kb dsDNA circular
    "Papillomaviridae":   {"sequence": {"region": "whole_genome"}},  # ~8 kb dsDNA circular
    # ---- Large genomes: marker gene approach ----
    # Adenoviridae (~35 kb) — hexon (major capsid protein, well annotated in NCBI)
    "Adenoviridae":       {"sequence": {"region": "hexon"}},
    # Herpesviruses (~130–240 kb) — DNA polymerase
    "Orthoherpesviridae": {"sequence": {"region": "DNA polymerase"}},
    "Herpesviridae":      {"sequence": {"region": "DNA polymerase"}},  # legacy ICTV name
    "Alloherpesviridae":  {"sequence": {"region": "DNA polymerase"}},
    "Malacoherpesviridae":{"sequence": {"region": "DNA polymerase"}},
    # Poxviridae (~130–375 kb) — DNA polymerase
    "Poxviridae":         {"sequence": {"region": "DNA polymerase"}},
    # Iridoviridae (~100–220 kb) — major capsid protein (conserved across genera)
    "Iridoviridae":       {"sequence": {"region": "major capsid protein"}},
    # Asfarviridae (~190 kb) — B646L (p72, major capsid protein of ASFV)
    "Asfarviridae":       {"sequence": {"region": "B646L"}},
    # Insect large dsDNA (~80–230 kb) — DNA polymerase
    "Baculoviridae":      {"sequence": {"region": "DNA polymerase"}},
    "Nudiviridae":        {"sequence": {"region": "DNA polymerase"}},
    "Ascoviridae":        {"sequence": {"region": "DNA polymerase"}},
}


DEFAULT_FAMILY_CONFIG: dict = {
    "download": {
        "max_per_species": 200,
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
        "options": "--6merpair --retree 1",
    },
    "msa_100": {
        "tool": "mafft",
        "options": "--retree 1",
    },
    "tree_500": {
        "tool": "fasttree",
        "options": "",
        "model_nuc": "GTR+G",
        "model_aa": "WAG+G",
    },
    "tree_100": {
        "tool": "iqtree",
        "options": "--fast",
        "model_nuc": "GTR+G",
        "model_aa": "WAG+G",
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
    if not cfg.get("ncbi", {}).get("email"):
        raise ValueError("global.yaml: ncbi.email must be set")


def load_family_config(family: str, configs_dir: Path, global_cfg: dict) -> tuple[dict, bool]:
    """Load per-family config, auto-generating it if missing.

    Returns (config_dict, was_auto_generated).
    """
    config_path = configs_dir / f"{family}.yaml"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
        cfg = _merge_with_defaults(cfg, global_cfg)
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


def _merge_with_defaults(cfg: dict, global_cfg: dict) -> dict:
    """Fill missing keys from global defaults."""
    global_defaults = global_cfg.get("defaults", {})
    merged = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _deep_update(merged, global_defaults)
    _deep_update(merged, cfg)
    return merged


def _generate_default_family_config(family: str, global_cfg: dict) -> dict:
    global_defaults = global_cfg.get("defaults", {})
    cfg = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _deep_update(cfg, global_defaults)

    # Auto-set segment for known segmented RNA families
    segment = SEGMENTED_FAMILIES.get(family)
    if segment:
        cfg["sequence"]["segment"] = segment
        log.info("Auto-configured segment '%s' for segmented family %s", segment, family)

    # Auto-configure DNA virus families
    dna_overrides = DNA_FAMILIES.get(family)
    if dna_overrides:
        _deep_update(cfg, dna_overrides)
        region = cfg["sequence"].get("region", "whole_genome")
        if region == "whole_genome":
            log.info("Auto-configured DNA family %s: whole-genome search", family)
        else:
            log.info("Auto-configured DNA family %s: marker gene '%s'", family, region)

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
