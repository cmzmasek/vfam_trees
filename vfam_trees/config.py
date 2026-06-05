"""Config loading, validation, and auto-generation for vfam_trees."""
from __future__ import annotations

import copy
import re
from datetime import date
from pathlib import Path

import yaml

from .datasources import (
    DEFAULT_MAX_SEQS as _ADS_DEFAULT_MAX_SEQS,
    KNOWN_PATHOPLEXUS_ORGANISMS,
    SOURCE_FETCHERS,
)
from .logger import get_logger

log = get_logger(__name__)


_PREFIX_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_output_prefix(name: str) -> str:
    """Return a filesystem-safe prefix for output files and directories.

    ``manual.name`` is free-text (it can contain spaces, ``|``, ``/``, ``[``,
    ``]`` — concat titles use ``[concatenated|N markers]`` style), so it must be
    sanitized before it is used in a path: runs of any character outside
    ``[A-Za-z0-9._-]`` collapse to a single underscore, and leading/trailing
    separators are stripped.  Clean inputs (e.g. ``Filoviridae``) are returned
    unchanged, so families without ``manual.name`` keep their current output
    names exactly.
    """
    s = _PREFIX_UNSAFE_RE.sub("_", (name or "").strip())
    s = re.sub(r"_+", "_", s).strip("_.")
    return s or "family"


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
    # Malacoherpesviridae — 5 species, only MCP and DNA pol annotated;
    # MCP is the better-covered single marker per 2026-05 cache audit.
    "Malacoherpesviridae":{"sequence": {"region": "major capsid protein", "type": "protein"}},
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


# Curated multi-marker presets for large DNA virus families where single-
# protein analysis carries insufficient phylogenetic signal.  Each marker
# spec: name (canonical), aliases (GenBank annotation variants), optional
# aliases_<subfamily> (subfamily-specific aliases applied at fetch time
# based on species lineage), optional length_range [min_aa, max_aa],
# optional locus_tag_hint (regex used as paralog tiebreaker).  When a
# family appears here it takes precedence over DNA_FAMILIES in auto-
# generated configs; the user can revert to single-protein mode by
# editing the per-family yaml (set sequence.region back to a single
# protein name and remove the concatenation block).
#
# Curation status: first pass from the literature; ASFV gene IDs are
# flagged TBV (verify against NC_001659).  See CONCAT_DESIGN.md §4 for
# references.
CONCATENATION_FAMILIES: dict[str, dict] = {

    # ---- Poxviridae (9 markers) -----------------------------------------
    # Refs: Upton et al. 2003; Hughes & Friedman 2005; ICTV Poxviridae.
    # Single 9-marker set covers both Chordopoxvirinae and Entomopoxvirinae;
    # subfamily-aware aliases handle entomopox annotation drift.
    "Poxviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {
                    "name": "DNA polymerase",
                    "aliases": ["DNA-directed DNA polymerase", "DNA pol"],
                    "aliases_Entomopoxvirinae": ["DNA polymerase B"],
                    "locus_tag_hint": r"E9L|polB",
                },
                {
                    "name": "DNA-directed RNA polymerase 147 kDa subunit",
                    "aliases": ["RPO147", "RNA polymerase subunit RPO147"],
                    "locus_tag_hint": r"A24R|RPO147|rpo147",
                },
                {
                    "name": "DNA-directed RNA polymerase 132 kDa subunit",
                    "aliases": ["RPO132", "RNA polymerase subunit RPO132"],
                    "locus_tag_hint": r"J6R|RPO132|rpo132",
                },
                {
                    "name": "mRNA capping enzyme large subunit",
                    "aliases": ["capping enzyme large subunit"],
                    "locus_tag_hint": r"D1R",
                },
                {
                    "name": "DNA helicase",
                    "aliases": ["NPH-II", "transcript release factor"],
                    "locus_tag_hint": r"A18R|NPH2",
                },
                {
                    "name": "poly(A) polymerase catalytic subunit",
                    "aliases": ["poly(A) polymerase large subunit"],
                    "locus_tag_hint": r"E1L|VP55",
                },
                {
                    "name": "late transcription factor VLTF-3",
                    "aliases": ["VLTF3", "late transcription factor 3"],
                    "locus_tag_hint": r"A1L|VLTF",
                },
                {
                    "name": "uracil-DNA glycosylase",
                    "aliases": ["UNG"],
                    "locus_tag_hint": r"D4R|UNG",
                },
                # ssDNA-binding protein (I3L) was previously included but a
                # 2026-05 cache audit (128 species) showed 1 % coverage —
                # the I3L name is essentially VACV-Cop-specific and not used
                # outside Chordopoxvirinae reference annotations.
            ],
        },
    },

    # ---- Herpesviridae and other herpesvirus families (7 markers) -------
    # Refs: McGeoch et al. 1995, 2006; Davison 2010.  Single family-wide
    # set; deliberately excludes glycoprotein B (UL27) — too divergent
    # across alpha/beta/gamma subfamilies for clean alignment.
    "Orthoherpesviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {
                    "name": "DNA polymerase catalytic subunit",
                    "aliases": ["DNA-directed DNA polymerase catalytic subunit", "DNA polymerase"],
                    "locus_tag_hint": r"UL30",
                },
                {
                    "name": "helicase-primase helicase subunit",
                    "aliases": ["DNA helicase"],
                    "locus_tag_hint": r"UL5\b",
                },
                {
                    "name": "helicase-primase primase subunit",
                    "aliases": ["primase"],
                    "locus_tag_hint": r"UL52",
                },
                {
                    "name": "major capsid protein",
                    "aliases": ["MCP", "capsid protein VP5"],
                    "locus_tag_hint": r"UL19|VP5\b",
                },
                {
                    "name": "capsid triplex subunit 2",
                    "aliases": ["VP23", "minor capsid protein"],
                    "locus_tag_hint": r"UL18|VP23",
                },
                {
                    "name": "DNA packaging terminase subunit 1",
                    "aliases": ["terminase ATPase subunit", "DNA packaging terminase ATPase"],
                    "locus_tag_hint": r"UL15",
                },
                {
                    "name": "single-stranded DNA-binding protein",
                    "aliases": ["major DNA-binding protein", "ICP8"],
                    "locus_tag_hint": r"UL29|ICP8",
                },
            ],
        },
    },

    # ---- Asfarviridae (6 markers — ASFV gene IDs TBV against NC_001659) -
    # Refs: Yutin & Koonin 2012; Iyer et al. 2006.
    "Asfarviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {
                    "name": "DNA polymerase B",
                    "aliases": ["DNA polymerase", "DNA-directed DNA polymerase"],
                    "locus_tag_hint": r"polB",
                },
                {
                    "name": "major capsid protein p72",
                    "aliases": ["major capsid protein", "p72", "capsid protein p72"],
                    "locus_tag_hint": r"B646L|p72",
                },
                {
                    "name": "packaging ATPase",
                    "aliases": ["A32-like ATPase", "FtsK-like ATPase"],
                    "locus_tag_hint": r"A32",
                },
                {
                    "name": "primase-helicase",
                    "aliases": ["D5-like helicase", "superfamily 3 helicase"],
                    "locus_tag_hint": r"D5|A18",
                },
                {
                    "name": "late transcription factor 3",
                    "aliases": ["VLTF-3"],
                    "locus_tag_hint": r"VLTF",
                },
                {
                    "name": "DNA-directed RNA polymerase subunit",
                    "aliases": ["RPB1-like subunit", "RNA polymerase largest subunit"],
                    "locus_tag_hint": r"NP1450L|RPB1",
                },
            ],
        },
    },

    # ---- Iridoviridae (6 markers) ---------------------------------------
    # Refs: Tidona & Darai 1997; Eaton et al. 2007; ICTV Iridoviridae.
    # Originally a 7-marker set; a first 2026-05 audit cut it to 4 because
    # packaging ATPase, D5-like helicase, and VLTF-3 each looked like 0 %
    # coverage.  After landing the [Title] field-fix in fetch.py
    # (NCBI's [Protein Name] index is sparsely populated for this family),
    # the re-audit showed packaging ATPase at 21 family-wide hits and D5
    # helicase at 12 — both worth keeping.  VLTF-3 stays dropped (still 0 %).
    "Iridoviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {"name": "major capsid protein", "aliases": ["MCP"], "locus_tag_hint": r"MCP"},
                {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase"], "locus_tag_hint": r"polB"},
                {"name": "packaging ATPase", "aliases": ["A32-like ATPase"], "locus_tag_hint": r"A32"},
                {"name": "ribonuclease III", "aliases": ["RNase III"], "locus_tag_hint": r"rnc"},
                {"name": "DNA helicase", "aliases": ["D5-like helicase"], "locus_tag_hint": r"D5|helicase"},
                {"name": "immediate-early protein ICP-46", "aliases": ["ICP46"], "locus_tag_hint": r"ICP46"},
            ],
        },
    },

    # ---- Baculoviridae and relatives (7 markers) ------------------------
    # Refs: Herniman et al. 2003; Jehle et al. 2006; Miele et al. 2011.
    # Same set used for Nudiviridae and Ascoviridae (closely related insect
    # dsDNA virus families that share the baculoviral core gene set).
    "Baculoviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase"], "locus_tag_hint": r"polB"},
                {"name": "late expression factor 8", "aliases": ["LEF-8", "RNA polymerase subunit LEF-8"], "locus_tag_hint": r"lef-?8"},
                {"name": "late expression factor 9", "aliases": ["LEF-9", "RNA polymerase subunit LEF-9"], "locus_tag_hint": r"lef-?9"},
                {"name": "DNA helicase P143", "aliases": ["p143", "helicase"], "locus_tag_hint": r"p143|helicase"},
                {"name": "per os infectivity factor 1", "aliases": ["PIF-1", "P74"], "locus_tag_hint": r"pif-?1|p74"},
                {"name": "per os infectivity factor 2", "aliases": ["PIF-2"], "locus_tag_hint": r"pif-?2"},
                {"name": "major capsid protein", "aliases": ["VP39", "capsid protein VP39"], "locus_tag_hint": r"vp39|MCP"},
            ],
        },
    },

    # ---- NCLDV hallmark fallback (8 markers) ----------------------------
    # Refs: Yutin & Koonin 2009, 2012; Koonin & Yutin 2019.  Applied to
    # large-DNA-virus families lacking a curated set (Mimi/Phyco/Marseille/
    # Pitho/Pandora/Medusaviridae and any future NCLDV-like family).
}

# Nudiviridae uses the Baculoviridae 7-gene set verbatim (same insect-dsDNA
# core, see Miele et al. 2011 §4.5 in CONCAT_DESIGN.md).  Ascoviridae was
# previously aliased to the same 7-gene set, but a cache audit (2026-05)
# showed that 0/27 NCBI-annotated ascovirus species carry recognizable
# lef-8/lef-9/pif-1/pif-2 protein names — so the concat fetch produced 0
# qualifying genomes at any usable min_fraction.  The 3-gene core below
# (DNA pol, MCP, P143 helicase) is the actual conserved annotated set in
# Ascoviridae and gives ~12 qualifying species at min_fraction 0.5.
CONCATENATION_FAMILIES["Nudiviridae"] = copy.deepcopy(CONCATENATION_FAMILIES["Baculoviridae"])

CONCATENATION_FAMILIES["Ascoviridae"] = {
    "sequence": {"region": "concatenated", "type": "protein"},
    "concatenation": {
        "proteins": [
            {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase"], "locus_tag_hint": r"polB"},
            {"name": "DNA helicase P143", "aliases": ["p143", "helicase"], "locus_tag_hint": r"p143|helicase"},
            {"name": "major capsid protein", "aliases": ["MCP", "capsid protein"], "locus_tag_hint": r"MCP"},
        ],
    },
}

# Herpesviridae (legacy ICTV family name) reuses the Orthoherpesviridae
# 7-gene core verbatim.
CONCATENATION_FAMILIES["Herpesviridae"] = copy.deepcopy(CONCATENATION_FAMILIES["Orthoherpesviridae"])

# Alloherpesviridae uses a 6-gene subset of the Orthoherpesviridae core:
# the single-stranded DNA-binding protein (UL29/ICP8) is dropped because a
# 2026-05 cache audit (24 species) showed 0 % coverage in fish-herpesvirus
# annotations.
CONCATENATION_FAMILIES["Alloherpesviridae"] = {
    "sequence": {"region": "concatenated", "type": "protein"},
    "concatenation": {
        "proteins": [
            {
                "name": "DNA polymerase catalytic subunit",
                "aliases": ["DNA-directed DNA polymerase catalytic subunit", "DNA polymerase"],
                "locus_tag_hint": r"UL30",
            },
            {
                "name": "helicase-primase helicase subunit",
                "aliases": ["DNA helicase"],
                "locus_tag_hint": r"UL5\b",
            },
            {
                "name": "helicase-primase primase subunit",
                "aliases": ["primase"],
                "locus_tag_hint": r"UL52",
            },
            {
                "name": "major capsid protein",
                "aliases": ["MCP", "capsid protein VP5"],
                "locus_tag_hint": r"UL19|VP5\b",
            },
            {
                "name": "capsid triplex subunit 2",
                "aliases": ["VP23", "minor capsid protein"],
                "locus_tag_hint": r"UL18|VP23",
            },
            {
                "name": "DNA packaging terminase subunit 1",
                "aliases": ["terminase ATPase subunit", "DNA packaging terminase ATPase"],
                "locus_tag_hint": r"UL15",
            },
        ],
    },
}

# Malacoherpesviridae has only 5 species and only MCP / DNA pol annotated
# (others all 0 %).  Concat is meaningless at that scale; the family
# instead falls back to single-marker mode via DNA_FAMILIES below
# (region: "major capsid protein").

# NCLDV-hallmark 8-marker fallback set.  Applied to large-DNA-virus
# families that lack a curated per-family set in CONCATENATION_FAMILIES
# but should still default to concatenated mode.
_NCLDV_HALLMARK_PROTEINS = [
    {
        "name": "DNA polymerase",
        "aliases": ["DNA-directed DNA polymerase", "polymerase B", "DNA polymerase B"],
        "locus_tag_hint": r"polB",
    },
    {
        "name": "major capsid protein",
        "aliases": ["MCP", "capsid protein"],
        "locus_tag_hint": r"MCP",
    },
    {
        "name": "packaging ATPase",
        "aliases": ["A32-like ATPase", "FtsK-like ATPase"],
        "locus_tag_hint": r"A32",
    },
    {
        "name": "primase-helicase",
        "aliases": ["D5-like helicase", "superfamily 3 helicase"],
        "locus_tag_hint": r"D5",
    },
    {
        "name": "late transcription factor 3",
        "aliases": ["VLTF-3"],
        "locus_tag_hint": r"VLTF",
    },
    {
        "name": "mRNA capping enzyme",
        "aliases": ["capping enzyme large subunit"],
    },
    {
        "name": "DNA-directed RNA polymerase subunit alpha",
        "aliases": ["RNA polymerase RPB1", "largest subunit RNA polymerase"],
        "locus_tag_hint": r"RPB1",
    },
    {
        "name": "DNA-directed RNA polymerase subunit beta",
        "aliases": ["RNA polymerase RPB2", "second-largest subunit RNA polymerase"],
        "locus_tag_hint": r"RPB2",
    },
]
for _ncldv_family in ("Pandoraviridae", "Medusaviridae"):
    CONCATENATION_FAMILIES[_ncldv_family] = {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {"proteins": copy.deepcopy(_NCLDV_HALLMARK_PROTEINS)},
    }

# NCLDV families with cache audits (2026-05) get pruned per-family presets
# instead of the 8-marker hallmark fallback.  In every case the markers
# removed had ≤ 4 % family-wide coverage in NCBI annotations and were just
# inflating the min_fraction denominator.

# Phycodnaviridae — first 2026-05 audit cut the NCLDV-8 set to 3 (DNA pol,
# MCP, capping enzyme).  After the [Title] fetch fix, packaging ATPase
# came back to 42 family-wide hits and is restored.  The other four
# hallmarks (primase-helicase, VLTF-3, RPB1, RPB2) stay dropped (still
# ≤ 1 % even with [Title]).
CONCATENATION_FAMILIES["Phycodnaviridae"] = {
    "sequence": {"region": "concatenated", "type": "protein"},
    "concatenation": {
        "proteins": [
            {
                "name": "DNA polymerase",
                "aliases": ["DNA-directed DNA polymerase", "polymerase B", "DNA polymerase B"],
                "locus_tag_hint": r"polB",
            },
            {
                "name": "major capsid protein",
                "aliases": ["MCP", "capsid protein"],
                "locus_tag_hint": r"MCP",
            },
            {
                "name": "packaging ATPase",
                "aliases": ["A32-like ATPase", "FtsK-like ATPase"],
                "locus_tag_hint": r"A32",
            },
            {
                "name": "mRNA capping enzyme",
                "aliases": ["capping enzyme large subunit"],
            },
        ],
    },
}

# Mimiviridae and Marseilleviridae share the same trim: drop primase-helicase
# (2–4 %) and VLTF-3 (0–4 %) from the NCLDV-8 set; the remaining six markers
# all have meaningful coverage in both families.  Note: Marseilleviridae DNA
# polymerase is at 0 % under the current aliases — flagged for follow-up
# alias investigation; the marker is kept in the preset on the assumption
# that the alias issue is fixable.
_NCLDV_TRIMMED_NO_HELICASE_NO_VLTF = [
    p for p in copy.deepcopy(_NCLDV_HALLMARK_PROTEINS)
    if p["name"] not in ("primase-helicase", "late transcription factor 3")
]
for _f in ("Mimiviridae", "Marseilleviridae"):
    CONCATENATION_FAMILIES[_f] = {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": copy.deepcopy(_NCLDV_TRIMMED_NO_HELICASE_NO_VLTF),
        },
    }

# Pithoviridae — different trim profile: packaging ATPase (0 %) and
# VLTF-3 (0 %) drop out, but primase-helicase is kept (it's at 4 % vs 0 %
# for the other dropouts and may improve with future alias work).  DNA
# polymerase is at 17 % — also flagged for follow-up alias investigation.
CONCATENATION_FAMILIES["Pithoviridae"] = {
    "sequence": {"region": "concatenated", "type": "protein"},
    "concatenation": {
        "proteins": [
            p for p in copy.deepcopy(_NCLDV_HALLMARK_PROTEINS)
            if p["name"] not in ("packaging ATPase", "late transcription factor 3")
        ],
    },
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
        "max_length": None,
        "max_ambiguous": 0.01,
        "exclude_organisms": [
            "synthetic construct",
            "metagenome",
            "MAG:",
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
        "options_nuc": "--6merpair --retree 2",
        "options_aa": "--auto",
    },
    "msa_100": {
        "tool": "mafft",
        "options_nuc": "--retree 3",
        "options_aa": "--maxiterate 1000 --localpair",
    },
    "msa_trim": {
        "enabled": True,
        "tool": "trimal",
        "options": "-automated1",
    },
    "tree_500": {
        "tool": "fasttree",
        "options_nuc": "",
        "options_aa": "",
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
    "concatenation": {
        # Multi-marker protein concatenation mode (sequence.region must be
        # set to "concatenated" for this block to take effect — see
        # CONCATENATION_FAMILIES for curated per-family marker presets,
        # and CONCAT_DESIGN.md for the full specification).
        #
        # proteins:        list of marker specs.  Each entry: {name (str),
        #                  aliases (list[str]), aliases_<subfamily> (optional
        #                  list[str]), length_range (optional [min_aa,
        #                  max_aa]), locus_tag_hint (optional regex)}.
        # min_fraction:    genomes with fewer than ceil(min_fraction × N)
        #                  markers are dropped.  Default 0.7.
        # partition_tree_100: partitioned IQ-TREE on tree_100 (one model per
        #                  marker).  Default true.
        # partition_tree_500: FastTree does not support partitioned analysis;
        #                  tree_500 in concat mode is single-model on the
        #                  full concatenation.  Default false (cannot be
        #                  enabled in MVP).
        "proteins": [],
        "min_fraction": 0.7,
        # Source-nuc length filter: drop source-nuc accessions whose parent
        # nucleotide record is shorter than this fraction of the longest
        # source-nuc parent in the per-species fetch.  Filters out single-
        # gene partial submissions (very common for ASFV p72, papillomavirus
        # L1, etc.) that otherwise crowd out genome-scale proteins under
        # max_per_species and prevent any genome from clearing min_fraction.
        # 0 disables the filter.
        "source_nuc_min_length_frac": 0.3,
        "partition_tree_100": True,
        "partition_tree_500": False,
    },
    "refseq_absorption": {
        # Pre-clustering step: drops non-RefSeq sequences that are near-
        # identical (≥ threshold) to a RefSeq within the same species.
        # Prevents the tree from showing redundant near-zero-branch
        # clusters of isolates around their RefSeq.  All RefSeqs are
        # always kept; only non-RefSeqs are absorbed.
        "enabled": True,
        "threshold": 0.99,
    },
    "length_outlier": {
        "enabled": True,
        "k": 5.0,
        "min_lo_mult": 0.20,
        "max_hi_mult": 5.0,
    },
    "outlier_removal": {
        "enabled": True,
        "factor": 20.0,
        "max_iterations": 3,
        "min_seqs": 40,
    },
    "labeling": {
        # Python-style format string controlling leaf labels in PhyloXML
        # <name> elements, PDF/PNG tree images, and FASTA display names.
        # Supported placeholders: {species}, {id} (accession), {host},
        # {strain}, {location}, {year}, {genus}.
        # Literal text (including separators like "|") is reproduced verbatim.
        "format": "{species}|{id}|{host}",
        # Replace spaces with underscores in field values (Newick-safe).
        "replace_whitespace": True,
        # When False (default): empty / "unknown" / "n/a" fields AND their
        # immediately preceding separator are dropped so no leading or
        # consecutive separators appear.  Set True to keep them.
        "keep_separator_on_empty": False,
    },
    "coloring": {
        # Genus inference strategy when no formal genus rank is present in the
        # NCBI lineage.  "none" keeps the current behaviour (grey).
        # "suffix"  — treat any single-word name ending in "virus" as a genus
        #             proxy (catches NCBI "no rank" nodes that are biologically
        #             genus-level).
        # "deepest" — suffix first, then fall back to the deepest lineage entry
        #             that is above species level (e.g. an unranked clade or a
        #             subfamily node), maximally aggressive.
        "genus_inference": "deepest",
    },
    "taxonomy": {
        # Minimum rank a leaf lineage must reach to participate in internal-node
        # LCA annotation.  Leaves whose lineage ends above this rank are excluded
        # from the LCA vote (they stay in the tree but do not blur internal
        # annotations).  "none" keeps the current behaviour (all leaves vote).
        # Typical values: "genus", "species", "subfamily".
        "lca_min_rank": "species",
    },
    "output": {
        # Emit an Auspice v2 JSON (Nextstrain interactive tree) alongside each
        # PhyloXML.  Divergence tree only — branch lengths are substitutions
        # per site, not time, so the tree is not time-resolved.  The "color by
        # genus" view reuses the same palette as the PDF/PNG/PhyloXML output.
        "auspice_json": True,
    },
    "manual": {
        # Curator overrides on per-family record selection.  include/exclude
        # hold exact accessions with version (e.g. "NC_002617.1").
        # include: force-keep — bypasses all QC (length, ambiguity, organism
        #          exclusion) and is protected through clustering, proportional
        #          merge, and length-outlier filtering (stronger than RefSeq at
        #          QC, equal to RefSeq downstream).
        # include_seq: paste sequences directly (e.g. records not yet in
        #          GenBank).  Each entry is a mapping with required keys
        #          'id', 'organism', and 'sequence'.  Injected after fetch,
        #          fully bypass QC, and receive the same downstream protection
        #          as manual.include records.  The 'id' must not collide with
        #          any accession returned by NCBI for this family or any id in
        #          manual.include / manual.exclude.  Not supported when
        #          sequence.region is 'concatenated'.
        # include_fasta_files: paths to one or more FASTA files whose sequences
        #          are injected identically to include_seq entries.  The FASTA
        #          id field (first whitespace-delimited token of the header) is
        #          used as the sequence id; the remainder of the header line
        #          becomes the organism/name.  A single-token header (no
        #          remainder) yields an empty organism, so the leaf label is
        #          just the header itself.  Subject to the same collision
        #          rules and concat-mode restriction as include_seq.
        # exclude: dropped immediately after fetch, before QC.
        # restrict_to_lineages: restrict the pipeline to species under one or more
        #          taxonomic lineages.  Each entry is a taxon at any rank
        #          (species, genus, subfamily, ...) given as a scientific name
        #          (matched against NCBI taxonomy) or a numeric taxid.  Only
        #          species whose taxid is a descendant of at least one listed
        #          taxon proceed to download; others are skipped entirely.
        #          When empty or absent the full discovered species list is used.
        # name: override the display name used in PDF/PNG titles and the PhyloXML
        #          <name> element.  When empty the biological family name is used.
        #          Also drives the output prefix: when set, every output file and
        #          the output directory are named from a filesystem-sanitized form
        #          of this value (e.g. results/Ebola_3052317/Ebola_tree_100.nwk)
        #          rather than the family name; the taxid suffix on the directory
        #          is retained.  NCBI queries, taxonomy, caching, and the
        #          summary.tsv 'family' column still use the biological family.
        "include": [],
        "include_seq": [],
        "include_fasta_files": [],
        "exclude": [],
        "restrict_to_lineages": [],
        "name": "",
    },
    # Supplement the tree with forced-include sequences fetched from an external
    # pathogen database (currently only Pathoplexus, via its GenSpectrum LAPIS
    # API).  Typical use: add the latest outbreak genomes to a family tree.
    # Like manual.include, these BYPASS QC + clustering/subsampling and are
    # protected from outlier removal; unlike pasted sequences they carry real
    # host/location/date/taxon metadata, so they colour by genus.
    # Nucleotide-only; not supported when sequence.region is 'concatenated'.
    # Each list entry is a mapping:
    #   source:        required; one of the supported sources ("pathoplexus").
    #   organism:      required; the source's organism slug (e.g. "ebola-zaire").
    #                  Pathoplexus is keyed by curated organism, not taxid, so a
    #                  genus tree lists several entries (ebola-zaire/sudan/bdbv).
    #   country:       optional, exact match (LAPIS geoLocCountry).
    #   host:          optional, exact match (LAPIS hostNameScientific).
    #   date_from:     optional ISO 'YYYY-MM-DD' (collection date lower bound).
    #   date_to:       optional ISO 'YYYY-MM-DD' (collection date upper bound).
    #   max_seqs:      optional cap (default 200); these skip subsampling.
    #   dedup_vs_ncbi: optional bool (default true); drop records whose INSDC
    #                  accession is already present in the NCBI download.
    #   name_prefix:   optional; prepended verbatim to each injected leaf label
    #                  (e.g. "PATHOPLEXUS_") so external tips stand out.
    #   outbreak_name: optional; tags each injected tip with this label as the
    #                  PhyloXML 'vipr:outbreak' property and an Auspice trait
    #                  (e.g. "Bdbv-2026") for colouring / filtering.
    "additional_data_sources": [],
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


def _validate_concatenation_block(cfg: dict, family: str) -> None:
    """Validate the concatenation block when sequence.region == "concatenated".

    Raises ValueError on malformed configuration so the run aborts before
    anything is fetched or aligned.  Validates only what's structurally
    required for downstream phases — semantic checks (e.g. "is this a real
    GenBank protein name") are left to the fetcher.
    """
    region = (cfg.get("sequence") or {}).get("region", "")
    if region != "concatenated":
        return

    block = cfg.get("concatenation") or {}
    proteins = block.get("proteins")
    if not isinstance(proteins, list) or len(proteins) == 0:
        raise ValueError(
            f"{family}: sequence.region is 'concatenated' but "
            "concatenation.proteins is empty.  Provide at least one marker "
            "spec (name + aliases) or revert sequence.region to a single "
            "protein name."
        )

    for i, protein in enumerate(proteins):
        if not isinstance(protein, dict):
            raise ValueError(
                f"{family}: concatenation.proteins[{i}] is not a dict — "
                "each marker must be a mapping with at least 'name' and "
                "'aliases' keys."
            )
        name = protein.get("name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError(
                f"{family}: concatenation.proteins[{i}].name is missing or "
                "empty — every marker spec needs a canonical protein name."
            )
        aliases = protein.get("aliases", [])
        if not isinstance(aliases, list) or not all(isinstance(a, str) for a in aliases):
            raise ValueError(
                f"{family}: concatenation.proteins[{i}].aliases must be a "
                "list of strings (use [] if no aliases)."
            )
        lr = protein.get("length_range")
        if lr is not None:
            if not (isinstance(lr, list) and len(lr) == 2 and
                    all(isinstance(x, (int, float)) for x in lr) and lr[0] < lr[1]):
                raise ValueError(
                    f"{family}: concatenation.proteins[{i}].length_range must "
                    "be a [min, max] list with min < max."
                )

    min_fraction = block.get("min_fraction", 0.7)
    if not (isinstance(min_fraction, (int, float)) and 0 < min_fraction <= 1):
        raise ValueError(
            f"{family}: concatenation.min_fraction must be a number in (0, 1] "
            f"(got {min_fraction!r})."
        )

    snf = block.get("source_nuc_min_length_frac", 0.3)
    if not (isinstance(snf, (int, float)) and 0 <= snf < 1):
        raise ValueError(
            f"{family}: concatenation.source_nuc_min_length_frac must be a "
            f"number in [0, 1) (got {snf!r})."
        )


# Declarative spec of numeric config fields and their permitted ranges.
# Each entry: (("section", "key"), kind, low, high, low_inclusive, high_inclusive)
# kind ∈ {"int", "number"}; low/high may be None for one-sided bounds.
# Catches typos like ``max_100: 4  00`` (YAML parses as the string "4  00")
# before they detonate inside a downstream ``%d`` log format or arithmetic op.
_NUMERIC_FIELD_SPECS: tuple = (
    (("download", "max_per_species"),         "int",    1,    None, True,  True),
    (("quality",  "max_ambiguous"),           "number", 0.0,  1.0,  True,  True),
    (("clustering", "threshold_min"),         "number", 0.0,  1.0,  False, True),
    (("clustering", "threshold_max"),         "number", 0.0,  1.0,  False, True),
    (("clustering", "max_reps_500"),          "int",    1,    None, True,  True),
    (("clustering", "max_reps_100"),          "int",    1,    None, True,  True),
    (("targets",    "max_500"),               "int",    1,    None, True,  True),
    (("targets",    "max_100"),               "int",    1,    None, True,  True),
    (("length_outlier", "k"),                 "number", 0.0,  None, True,  True),
    (("length_outlier", "min_lo_mult"),       "number", 0.0,  None, True,  True),
    (("length_outlier", "max_hi_mult"),       "number", 0.0,  None, True,  True),
    (("outlier_removal", "factor"),           "number", 0.0,  None, True,  True),
    (("outlier_removal", "max_iterations"),   "int",    0,    None, True,  True),
    (("outlier_removal", "min_seqs"),         "int",    1,    None, True,  True),
    (("refseq_absorption", "threshold"),      "number", 0.0,  1.0,  False, True),
)

# quality.min_length / quality.max_length are allowed to be None (= no bound).
_NULLABLE_NUMERIC_FIELDS: tuple = (
    (("quality", "min_length"), "int",    1,    None, True, True),
    (("quality", "max_length"), "int",    1,    None, True, True),
)


def _validate_numeric_fields(cfg: dict, family: str) -> None:
    """Validate that numeric config fields are actually numeric and in range.

    A common failure mode is a YAML typo like ``max_100: 4  00`` (two spaces),
    which PyYAML parses as the string ``"4  00"`` and then silently propagates
    until it crashes inside a downstream ``%d`` log format or arithmetic
    operation, often after a long fetch+MSA+tree-inference run.

    Raises ``ValueError`` with a precise location so the run aborts immediately
    at load time instead.
    """
    for (section, key), kind, lo, hi, lo_inc, hi_inc in _NUMERIC_FIELD_SPECS:
        if section not in cfg or key not in (cfg.get(section) or {}):
            continue
        val = cfg[section][key]
        _check_numeric(val, section, key, kind, lo, hi, lo_inc, hi_inc, family, allow_none=False)

    for (section, key), kind, lo, hi, lo_inc, hi_inc in _NULLABLE_NUMERIC_FIELDS:
        if section not in cfg or key not in (cfg.get(section) or {}):
            continue
        val = cfg[section][key]
        if val is None:
            continue
        _check_numeric(val, section, key, kind, lo, hi, lo_inc, hi_inc, family, allow_none=True)


def _check_numeric(
    val,
    section: str,
    key: str,
    kind: str,
    lo,
    hi,
    lo_inc: bool,
    hi_inc: bool,
    family: str,
    allow_none: bool,
) -> None:
    # bool is a subclass of int; reject it explicitly so True/False can't
    # masquerade as a numeric setting.
    if isinstance(val, bool):
        raise ValueError(
            f"{family}: {section}.{key} must be a {kind}, got bool {val!r}."
        )
    if kind == "int":
        if not isinstance(val, int):
            hint = _typo_hint(val)
            none_hint = " (use null or omit the key for no limit)" if allow_none else ""
            raise ValueError(
                f"{family}: {section}.{key} must be an integer, got "
                f"{type(val).__name__} {val!r}.{hint}{none_hint}"
            )
    else:  # "number"
        if not isinstance(val, (int, float)):
            hint = _typo_hint(val)
            raise ValueError(
                f"{family}: {section}.{key} must be a number, got "
                f"{type(val).__name__} {val!r}.{hint}"
            )
    if lo is not None:
        if (val < lo) if lo_inc else (val <= lo):
            op = ">=" if lo_inc else ">"
            raise ValueError(
                f"{family}: {section}.{key}={val!r} must be {op} {lo}."
            )
    if hi is not None:
        if (val > hi) if hi_inc else (val >= hi):
            op = "<=" if hi_inc else "<"
            raise ValueError(
                f"{family}: {section}.{key}={val!r} must be {op} {hi}."
            )


def _typo_hint(val) -> str:
    """Return a hint string when *val* looks like a YAML-typo'd number."""
    if isinstance(val, str) and any(c.isdigit() for c in val):
        return (
            f"  (this looks like a YAML typo — check for stray spaces or "
            f"non-digit characters in {val!r})"
        )
    return ""


def _validate_manual_block(cfg: dict, family: str) -> None:
    """Validate the optional manual.{include,include_seq,include_fasta_files,exclude} entries.

    include / exclude are lists of non-empty accession strings and must be
    disjoint.  include_seq is a list of {id, organism, sequence} mappings;
    include_fasta_files is a list of path strings; both share the same id
    collision rules and concat-mode restriction.
    Whitespace is stripped and accession-list duplicates are deduped in-place
    with a warning — typos are easier to spot at the line level than after a
    full pipeline run.
    """
    block = cfg.get("manual") or {}
    if not isinstance(block, dict):
        raise ValueError(
            f"{family}: 'manual' must be a mapping with 'include', "
            "'include_seq', 'include_fasta_files', and 'exclude' entries."
        )

    if block.get("include_fasta"):
        log.warning(
            "%s: manual.include_fasta is renamed to manual.include_seq — "
            "please rename the key in your config.",
            family,
        )

    if "include_species" in block:
        raise ValueError(
            f"{family}: manual.include_species has been renamed to "
            f"manual.restrict_to_lineages and now matches at any taxonomic "
            f"rank (species, genus, subfamily, ...).  Rename the key in your "
            f"per-family config."
        )

    if "limit_lineages" in block:
        raise ValueError(
            f"{family}: manual.limit_lineages has been renamed to "
            f"manual.restrict_to_lineages (the previous name implied the "
            f"wrong polarity).  Rename the key in your per-family config."
        )

    for key in ("include", "exclude"):
        raw = block.get(key) or []
        if not isinstance(raw, list):
            raise ValueError(
                f"{family}: manual.{key} must be a list of accession strings (got {type(raw).__name__})."
            )
        cleaned: list[str] = []
        seen: set[str] = set()
        for i, entry in enumerate(raw):
            if not isinstance(entry, str):
                raise ValueError(
                    f"{family}: manual.{key}[{i}] must be a string accession (got {type(entry).__name__})."
                )
            stripped = entry.strip()
            if not stripped:
                raise ValueError(
                    f"{family}: manual.{key}[{i}] is empty — remove the entry or fill it in."
                )
            if stripped in seen:
                log.warning("%s: duplicate accession %r in manual.%s — deduped.",
                            family, stripped, key)
                continue
            seen.add(stripped)
            cleaned.append(stripped)
        block[key] = cleaned

    overlap = set(block.get("include", [])) & set(block.get("exclude", []))
    if overlap:
        raise ValueError(
            f"{family}: manual.include and manual.exclude overlap on "
            f"{sorted(overlap)} — an accession cannot be both forced-in and dropped."
        )

    fasta_raw = block.get("include_seq") or []
    if not isinstance(fasta_raw, list):
        raise ValueError(
            f"{family}: manual.include_seq must be a list of mappings with "
            f"keys id/organism/sequence (got {type(fasta_raw).__name__})."
        )
    cleaned_fasta: list[dict] = []
    fasta_ids: set[str] = set()
    for i, entry in enumerate(fasta_raw):
        if not isinstance(entry, dict):
            raise ValueError(
                f"{family}: manual.include_seq[{i}] must be a mapping with "
                f"keys id/organism/sequence (got {type(entry).__name__})."
            )
        cleaned_entry: dict = {}
        for fkey in ("id", "organism", "sequence"):
            val = entry.get(fkey)
            if not isinstance(val, str):
                raise ValueError(
                    f"{family}: manual.include_seq[{i}].{fkey} must be a "
                    f"non-empty string (got {type(val).__name__})."
                )
            stripped = val.strip()
            if not stripped:
                raise ValueError(
                    f"{family}: manual.include_seq[{i}].{fkey} is empty — "
                    "remove the entry or fill it in."
                )
            cleaned_entry[fkey] = stripped
        seq_compact = "".join(cleaned_entry["sequence"].split()).upper()
        if not seq_compact:
            raise ValueError(
                f"{family}: manual.include_seq[{i}].sequence is empty after "
                "stripping whitespace."
            )
        cleaned_entry["sequence"] = seq_compact

        if cleaned_entry["id"] in fasta_ids:
            raise ValueError(
                f"{family}: manual.include_seq[{i}].id duplicates an earlier "
                f"entry ({cleaned_entry['id']!r})."
            )
        fasta_ids.add(cleaned_entry["id"])
        cleaned_fasta.append(cleaned_entry)

    clash_include = fasta_ids & set(block.get("include", []))
    if clash_include:
        raise ValueError(
            f"{family}: manual.include_seq ids overlap with manual.include "
            f"accessions: {sorted(clash_include)}."
        )
    clash_exclude = fasta_ids & set(block.get("exclude", []))
    if clash_exclude:
        raise ValueError(
            f"{family}: manual.include_seq ids overlap with manual.exclude "
            f"accessions: {sorted(clash_exclude)}."
        )

    files_raw = block.get("include_fasta_files") or []
    if not isinstance(files_raw, list):
        raise ValueError(
            f"{family}: manual.include_fasta_files must be a list of path "
            f"strings (got {type(files_raw).__name__})."
        )
    cleaned_files: list[str] = []
    for i, entry in enumerate(files_raw):
        if not isinstance(entry, str):
            raise ValueError(
                f"{family}: manual.include_fasta_files[{i}] must be a path "
                f"string (got {type(entry).__name__})."
            )
        stripped = entry.strip()
        if not stripped:
            raise ValueError(
                f"{family}: manual.include_fasta_files[{i}] is empty — "
                "remove the entry or fill it in."
            )
        cleaned_files.append(stripped)

    if cleaned_fasta or cleaned_files:
        region = (cfg.get("sequence") or {}).get("region", "")
        if region == "concatenated":
            raise ValueError(
                f"{family}: manual.include_seq / include_fasta_files are not "
                "supported when sequence.region is 'concatenated' — pasted "
                "sequences cannot be split into the configured marker proteins.  "
                "Remove the include_seq / include_fasta_files entries or switch "
                "the region."
            )

    block["include_seq"] = cleaned_fasta
    block["include_fasta_files"] = cleaned_files

    lineages_raw = block.get("restrict_to_lineages") or []
    if not isinstance(lineages_raw, list):
        raise ValueError(
            f"{family}: manual.restrict_to_lineages must be a list of taxon "
            f"names or taxids (got {type(lineages_raw).__name__})."
        )
    cleaned_lineages: list[str] = []
    seen_lineages: set[str] = set()
    for i, entry in enumerate(lineages_raw):
        if isinstance(entry, int):
            val = str(entry)
        elif isinstance(entry, str):
            val = entry.strip()
            if not val:
                raise ValueError(
                    f"{family}: manual.restrict_to_lineages[{i}] is empty — "
                    "remove the entry or fill it in."
                )
        else:
            raise ValueError(
                f"{family}: manual.restrict_to_lineages[{i}] must be a taxon "
                f"name (string) or a taxid (integer), got "
                f"{type(entry).__name__}."
            )
        if val in seen_lineages:
            log.warning(
                "%s: duplicate entry %r in manual.restrict_to_lineages — deduped.",
                family, val,
            )
            continue
        seen_lineages.add(val)
        cleaned_lineages.append(val)
    block["restrict_to_lineages"] = cleaned_lineages

    name_raw = block.get("name")
    if name_raw is None:
        block["name"] = ""
    elif not isinstance(name_raw, str):
        raise ValueError(
            f"{family}: manual.name must be a string (got {type(name_raw).__name__})."
        )
    else:
        block["name"] = name_raw.strip()

    cfg["manual"] = block


def _validate_additional_data_sources(cfg: dict, family: str) -> None:
    """Validate and normalise the optional ``additional_data_sources`` block.

    Each entry must name a supported ``source`` and a valid ``organism`` slug;
    optional ``country``/``host``/``name_prefix``/``outbreak_name`` are strings,
    ``date_from``/``date_to`` ISO dates, ``max_seqs`` a positive int,
    ``dedup_vs_ncbi`` a bool.  Defaults are
    filled in so the pipeline reads a fully-populated, typed entry.  No network
    I/O happens here.  Not supported when ``sequence.region`` is
    ``'concatenated'`` (mirrors the manual.include_seq restriction).
    """
    raw = cfg.get("additional_data_sources")
    if raw is None:
        cfg["additional_data_sources"] = []
        return
    if not isinstance(raw, list):
        raise ValueError(
            f"{family}: additional_data_sources must be a list of mappings "
            f"(got {type(raw).__name__})."
        )

    cleaned: list[dict] = []
    for i, entry in enumerate(raw):
        if not isinstance(entry, dict):
            raise ValueError(
                f"{family}: additional_data_sources[{i}] must be a mapping "
                f"(got {type(entry).__name__})."
            )
        source = entry.get("source")
        if not isinstance(source, str) or source not in SOURCE_FETCHERS:
            raise ValueError(
                f"{family}: additional_data_sources[{i}].source must be one of "
                f"{sorted(SOURCE_FETCHERS)} (got {source!r})."
            )
        organism = entry.get("organism")
        if not isinstance(organism, str) or not organism.strip():
            raise ValueError(
                f"{family}: additional_data_sources[{i}].organism must be a "
                "non-empty string (the source's organism slug)."
            )
        organism = organism.strip()
        if source == "pathoplexus" and organism not in KNOWN_PATHOPLEXUS_ORGANISMS:
            raise ValueError(
                f"{family}: additional_data_sources[{i}].organism {organism!r} "
                f"is not a known Pathoplexus organism. Valid organisms: "
                f"{sorted(KNOWN_PATHOPLEXUS_ORGANISMS)}."
            )

        for dkey in ("date_from", "date_to"):
            dval = entry.get(dkey)
            if dval is None:
                continue
            if not isinstance(dval, str):
                raise ValueError(
                    f"{family}: additional_data_sources[{i}].{dkey} must be an "
                    f"ISO date string 'YYYY-MM-DD' (got {type(dval).__name__})."
                )
            try:
                date.fromisoformat(dval.strip())
            except ValueError as exc:
                raise ValueError(
                    f"{family}: additional_data_sources[{i}].{dkey} must be an "
                    f"ISO date 'YYYY-MM-DD' (got {dval!r})."
                ) from exc

        for skey in ("country", "host", "name_prefix", "outbreak_name"):
            sval = entry.get(skey)
            if sval is not None and not isinstance(sval, str):
                raise ValueError(
                    f"{family}: additional_data_sources[{i}].{skey} must be a "
                    f"string (got {type(sval).__name__})."
                )

        max_seqs = entry.get("max_seqs", _ADS_DEFAULT_MAX_SEQS)
        if not isinstance(max_seqs, int) or isinstance(max_seqs, bool) or max_seqs <= 0:
            raise ValueError(
                f"{family}: additional_data_sources[{i}].max_seqs must be a "
                f"positive integer (got {max_seqs!r})."
            )
        dedup = entry.get("dedup_vs_ncbi", True)
        if not isinstance(dedup, bool):
            raise ValueError(
                f"{family}: additional_data_sources[{i}].dedup_vs_ncbi must be a "
                f"boolean (got {type(dedup).__name__})."
            )

        cleaned.append({
            "source": source,
            "organism": organism,
            "country": (entry.get("country") or "").strip() or None,
            "host": (entry.get("host") or "").strip() or None,
            "date_from": (entry.get("date_from") or "").strip() or None,
            "date_to": (entry.get("date_to") or "").strip() or None,
            "max_seqs": max_seqs,
            "dedup_vs_ncbi": dedup,
            # Presentation-only: prepended verbatim to each injected leaf label,
            # and emitted as the vipr:outbreak property / Auspice trait.
            "name_prefix": (entry.get("name_prefix") or ""),
            "outbreak_name": (entry.get("outbreak_name") or "").strip(),
        })

    if cleaned:
        region = (cfg.get("sequence") or {}).get("region", "")
        if region == "concatenated":
            raise ValueError(
                f"{family}: additional_data_sources are not supported when "
                "sequence.region is 'concatenated' — outbreak nucleotide "
                "sequences cannot be split into the configured marker proteins.  "
                "Remove the additional_data_sources entries or switch the region."
            )
    cfg["additional_data_sources"] = cleaned


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
        _validate_numeric_fields(cfg, family)
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
        _validate_additional_data_sources(cfg, family)
        return cfg, False
    else:
        cfg = _generate_default_family_config(family, global_cfg)
        _validate_numeric_fields(cfg, family)
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
        _validate_additional_data_sources(cfg, family)
        _write_family_config(cfg, config_path)
        log.warning(
            "No config found for %s — auto-generated default at %s. "
            "Edit it to tune parameters for this family.",
            family,
            config_path,
        )
        return cfg, True


def _warn_smart_default_conflicts(family: str, file_cfg: dict, config_path: Path) -> None:
    """Warn when a config file overrides a CONCATENATION_FAMILIES, DNA_FAMILIES,
    or SEGMENTED_FAMILIES setting.  Concat presets are the most specific so
    they're checked first; the DNA_FAMILIES check is suppressed when the
    family is in CONCATENATION_FAMILIES (the user is intentionally reverting
    concat → single-protein, which is documented behaviour, not a mistake).
    """
    file_region = (file_cfg.get("sequence") or {}).get("region")

    if family in CONCATENATION_FAMILIES:
        if file_region and file_region != "concatenated":
            log.info(
                "%s: %s sets sequence.region=%r, reverting from the recommended "
                "concatenated multi-marker mode to single-protein analysis. "
                "If this is unintentional, delete %s and re-run to regenerate.",
                family, config_path.name, file_region, config_path.name,
            )
    else:
        dna_overrides = DNA_FAMILIES.get(family)
        if dna_overrides:
            expected_region = (dna_overrides.get("sequence") or {}).get("region")
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
    _apply_smart_defaults(family, merged, file_cfg=cfg)
    _deep_update(merged, cfg)
    return merged


def _apply_smart_defaults(family: str, cfg: dict, file_cfg: dict | None = None) -> None:
    """Apply CONCATENATION_FAMILIES, DNA_FAMILIES, and SEGMENTED_FAMILIES as
    smart defaults in-place.  Precedence: concat > single-protein DNA
    overrides; segmentation is orthogonal and always applied when relevant.

    *file_cfg* is the user's family config dict (before it is merged).  When
    provided, "Auto-configured" log messages are suppressed for any value the
    file already overrides — those messages only make sense when the smart
    default is actually taking effect.
    """
    _file_seq = (file_cfg or {}).get("sequence") or {}

    segment = SEGMENTED_FAMILIES.get(family)
    if segment and not cfg["sequence"].get("segment"):
        cfg["sequence"]["segment"] = segment
        if not _file_seq.get("segment"):
            log.info("Auto-configured segment '%s' for segmented family %s", segment, family)

    concat_overrides = CONCATENATION_FAMILIES.get(family)
    if concat_overrides:
        _deep_update(cfg, concat_overrides)
        # Concat mode fetches one query per marker, so a single capped value
        # has to cover all markers' worth of partial submissions before the
        # complete-genome RefSeq proteins surface in the result set.  Bump
        # to 3000 so families with many partial single-gene submissions
        # (e.g. Asfarviridae's p72) don't crowd out genome-scale proteins.
        cfg.setdefault("download", {})
        if cfg["download"].get("max_per_species") in (None, 300):
            cfg["download"]["max_per_species"] = 3000
        n_markers = len(cfg.get("concatenation", {}).get("proteins", []))
        file_region = _file_seq.get("region")
        if not file_region or file_region == "concatenated":
            log.info(
                "Auto-configured concatenation mode for %s: %d markers",
                family, n_markers,
            )
        return  # concat presets supersede DNA_FAMILIES single-protein defaults

    dna_overrides = DNA_FAMILIES.get(family)
    if dna_overrides:
        _deep_update(cfg, dna_overrides)
        region = cfg["sequence"].get("region", "whole_genome")
        file_region = _file_seq.get("region")
        if not file_region or file_region == region:
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


def make_minimal_global_cfg() -> dict:
    """Return a minimal global config built entirely from hardcoded defaults.

    Used when no global.yaml is present (e.g. init-configs without a prior
    init run).  The empty ``defaults`` dict means _generate_default_family_config
    will use DEFAULT_FAMILY_CONFIG unchanged.
    """
    return {"ncbi": {}, "defaults": {}}


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
