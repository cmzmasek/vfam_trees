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
                {
                    "name": "single-stranded DNA-binding protein",
                    "aliases": ["ssDNA binding protein"],
                    "locus_tag_hint": r"I3L|ssb",
                },
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

    # ---- Iridoviridae (7 markers) ---------------------------------------
    # Refs: Tidona & Darai 1997; Eaton et al. 2007; ICTV Iridoviridae.
    "Iridoviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {
            "proteins": [
                {"name": "major capsid protein", "aliases": ["MCP"], "locus_tag_hint": r"MCP"},
                {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase"], "locus_tag_hint": r"polB"},
                {"name": "packaging ATPase", "aliases": ["A32-like ATPase"], "locus_tag_hint": r"A32"},
                {"name": "ribonuclease III", "aliases": ["RNase III"], "locus_tag_hint": r"rnc"},
                {"name": "DNA helicase", "aliases": ["D5-like helicase"], "locus_tag_hint": r"D5|helicase"},
                {"name": "late transcription factor 3", "aliases": ["VLTF-3"], "locus_tag_hint": r"VLTF"},
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

# Nudiviridae and Ascoviridae use the Baculoviridae 7-gene set verbatim
# (same insect-dsDNA core, see Miele et al. 2011 §4.5 in CONCAT_DESIGN.md).
CONCATENATION_FAMILIES["Nudiviridae"] = copy.deepcopy(CONCATENATION_FAMILIES["Baculoviridae"])
CONCATENATION_FAMILIES["Ascoviridae"] = copy.deepcopy(CONCATENATION_FAMILIES["Baculoviridae"])

# Herpesviridae (legacy ICTV family name) and Alloherpesviridae +
# Malacoherpesviridae use the same 7-gene herpesvirus core.
CONCATENATION_FAMILIES["Herpesviridae"]       = copy.deepcopy(CONCATENATION_FAMILIES["Orthoherpesviridae"])
CONCATENATION_FAMILIES["Alloherpesviridae"]   = copy.deepcopy(CONCATENATION_FAMILIES["Orthoherpesviridae"])
CONCATENATION_FAMILIES["Malacoherpesviridae"] = copy.deepcopy(CONCATENATION_FAMILIES["Orthoherpesviridae"])

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
for _ncldv_family in (
    "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
    "Pithoviridae", "Pandoraviridae", "Medusaviridae",
):
    CONCATENATION_FAMILIES[_ncldv_family] = {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {"proteins": copy.deepcopy(_NCLDV_HALLMARK_PROTEINS)},
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
    "manual": {
        # Curator overrides on per-family record selection.  Both lists hold
        # exact accessions with version (e.g. "NC_002617.1").
        # include: force-keep — bypasses all QC (length, ambiguity, organism
        #          exclusion) and is protected through clustering, proportional
        #          merge, and length-outlier filtering (stronger than RefSeq at
        #          QC, equal to RefSeq downstream).
        # exclude: dropped immediately after fetch, before QC.
        "include": [],
        "exclude": [],
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


def _validate_manual_block(cfg: dict, family: str) -> None:
    """Validate the optional manual.include / manual.exclude lists.

    Both must be lists of non-empty accession strings; the two lists must be
    disjoint.  Whitespace is stripped and duplicates within a list are deduped
    in-place with a warning — typos are easier to spot at the line level than
    after a full pipeline run.
    """
    block = cfg.get("manual") or {}
    if not isinstance(block, dict):
        raise ValueError(
            f"{family}: 'manual' must be a mapping with 'include' and 'exclude' lists."
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

    cfg["manual"] = block


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
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
        return cfg, False
    else:
        cfg = _generate_default_family_config(family, global_cfg)
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
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
    _apply_smart_defaults(family, merged)
    _deep_update(merged, cfg)
    return merged


def _apply_smart_defaults(family: str, cfg: dict) -> None:
    """Apply CONCATENATION_FAMILIES, DNA_FAMILIES, and SEGMENTED_FAMILIES as
    smart defaults in-place.  Precedence: concat > single-protein DNA
    overrides; segmentation is orthogonal and always applied when relevant.
    """
    segment = SEGMENTED_FAMILIES.get(family)
    if segment and not cfg["sequence"].get("segment"):
        cfg["sequence"]["segment"] = segment
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
        log.info(
            "Auto-configured concatenation mode for %s: %d markers",
            family, n_markers,
        )
        return  # concat presets supersede DNA_FAMILIES single-protein defaults

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
