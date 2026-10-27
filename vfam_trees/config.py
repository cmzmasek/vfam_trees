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
from .quality import MIN_LENGTH_AA
from .trim import DEFAULT_TRIM_OPTIONS

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
# whole_genome.  Nearly all families use a conserved marker protein.
# Circoviridae and Hepadnaviridae moved to the HMM-identified Rep and RT
# markers (HMM_MARKER_FAMILIES) in v1.6.0.
DNA_FAMILIES: dict[str, dict] = {
    # ---- Small ssDNA / dsDNA families: marker-protein approach ----
    # CRESS-DNA families — ambisense organisation (Rep and Cap on opposite
    # strands around the origin of replication) makes whole-genome
    # alignment unreliable; Rep is the ICTV species-demarcation marker.
    # The others are on the HMM-identified Rep (REP_FAMILIES); Smacoviridae
    # keeps the name-based Rep because Pfam finds its Rep in only 43 % of
    # RefSeq viruses (38/88, 2026-10-21).
    "Smacoviridae":       {"sequence": {"region": "Rep", "type": "protein"}},          # ~2.3 kb ssDNA circular
    # Anelloviridae — ORF1 (capsid) is the standard phylogenetic marker;
    # ORF2/ORF3 are hypervariable and make whole-genome alignment poor.
    "Anelloviridae":      {"sequence": {"region": "ORF1", "type": "protein"}},         # ~2–4 kb ssDNA circular
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
    # Adenoviridae (~35 kb) — hexon (major capsid protein).  HMM-identified by
    # the adenovirus hexon domain architecture Adeno_hexon (PF01065, N-term)
    # --> Adeno_hexon_C (PF03678, C-term), OR'd with the single Adeno_hexon
    # domain as a fallback.  The architecture token is tried first (exercises
    # N-to-C domain-order matching); the single-domain fallback recovers the
    # many real hexons it misses — a 2026-06 live run over 761 species found
    # the strict 2-domain gate is highly specific but under-sensitive, because
    # Adeno_hexon_C misses its GA cutoff on divergent adenoviruses and NCBI
    # carries fragmented/partial hexon deposits (the two domains land on
    # different records).  Tokens are OR'd: a protein satisfying EITHER wins.
    "Adenoviridae": {
        "sequence": {
            "region": "hexon", "type": "protein",
            "hmms": ["Adeno_hexon--Adeno_hexon_C", "Adeno_hexon"],
        },
        "hmm": {"enabled": True, "database": "Adenoviridae"},
    },
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
    # Mamonoviridae (medusaviruses, ~360–560 kb amoebal; "Medusaviridae" was
    # the proposed name, which NCBI never used, so that preset never applied
    # and the family ran as whole-genome nucleotide: MAFFT on 4 genomes ran
    # for hours).  Family-B DNA polymerase by HMM, cut to the domain ± 100
    # aa: two of the four genomes (2026-10) split the gene into two proteins
    # (Medusavirus 'euryale': exonuclease part / polymerase part;
    # Furtivovirus B10: within the polymerase domain, larger part 65 % of the
    # profile), and the spliced polymerases of M. medusae and M. sthenus
    # carry a ~680-aa N-terminal region the others lack.
    "Mamonoviridae": {
        "sequence": {
            "region": "DNA polymerase", "type": "protein",
            "hmms": ["DNA_pol_B"], "domain_padding": 100,
        },
        "hmm": {"enabled": True, "database": "Mamonoviridae"},
    },
}


# RNA-virus families (kingdom Orthornavirae) — RdRp protein marker.
#
# Whole-genome nucleotide alignment is saturated at family scale for most RNA
# viruses: a 2026-10-20 survey of earlier results found a median tip-to-tip
# distance above 1.5 subs/site in 84 of 117 nucleotide families (signal.py,
# CHANGELOG 1.6.0), i.e. their deep structure was noise.  The RNA-dependent
# RNA polymerase is the one protein every Orthornavirae member encodes and the
# standard family-level marker.  It is found by HMM homology against the
# bundled Pfam RdRp set ("RdRp", 23 profiles) whatever GenBank calls the
# protein (NIb, nsP4, 2a, lambda-3, L, PB1, polyprotein, …), and cut to the
# RdRp domain ± RDRP_DOMAIN_PADDING residues, so whole polyproteins do not
# dominate the alignment with variable regions.
#
# Hit rule (hmm block below): domain E-value <= 1e-3 and >= 20% profile
# coverage, no Pfam GA.  Every profile in the set is an RdRp, so the looser
# rule costs little: on 2026-10-20 RefSeq proteins it raised the families with
# >= 70% of viruses covered from 78 to 107 of 138, and an audit of the extra
# hits found them all RdRps under other names except ~4 in several thousand.
#
# Families listed here had >= 70% of their RefSeq viruses covered.  The other
# Orthornavirae families stay nucleotide (RefSeq viruses covered / total;
# 0/0 = no RefSeq proteins) — their RdRps are too divergent for Pfam:
# Alphaormycoviridae 0/0, Alphatetraviridae 4/6, Alvernaviridae 0/1,
# Amnoonviridae 0/1, Aspiviridae 0/6, Betaormycoviridae 0/0,
# Deltaormycoviridae 0/0, Discoviridae 3/5, Dumbiviridae 0/7, Euroniviridae
# 2/3, Gammaormycoviridae 0/0, Hypoviridae 5/30, Konkoviridae 0/0,
# Medioniviridae 1/2, Mesoniviridae 10/15, Mononiviridae 0/1, Mycoalphaviridae
# 0/0, Mypoviridae 0/1, Narnaviridae 5/15, Ourmiaviridae 3/12, Phasmaviridae
# 10/33, Qinviridae 3/8, Quambiviridae 0/4, Rhizouliviridae 0/5,
# Splipalmiviridae 0/0, Tosoviridae 0/1, Trimbiviridae 0/5, Tulasviridae 1/2,
# Unambiviridae 0/4, Wupedeviridae 0/1, Yueviridae 0/2
# Not here by design: reverse-transcribing families (Pararnavirae; no RdRp;
# they are on the RT marker, RT_FAMILY_PROFILES), Kolmioviridae (deltaviruses;
# host polymerase), satellites and viroids.
RDRP_HMM_TOKENS: list[str] = [
    "RdRP_1", "RdRP_2", "RdRP_3", "RdRP_4", "RdRP_5", "Mononeg_RNA_pol",
    "Bunya_RdRp", "Birna_RdRp_palm", "Arena_RNA_pol", "Flu_PB1", "Flavi_NS5",
    "Flavi_NS5_thumb", "RNA_replicase_B", "Mitovir_RNA_pol", "Orbi_VP1",
    "CoV_RPol_N", "Permu_RdRp_palm", "RNA-dir_pol_phytoreovirus",
    "Seadorna_RNAP", "CPV_RdRP_pol_dom", "Luteo_P1-P2", "L_thumb_ring_vir",
    "L_PA-C-like",
]
RDRP_DOMAIN_PADDING = 100

# The profiles each family's own RdRps hit (RDRP_FAMILY_PROFILES) are the
# family's marker tokens, and the preset turns on hmm.best_hit_only: a
# protein whose best-scoring profile in the whole set is not one of them is a
# different kind of RdRp and is left out.  This is what catches records filed
# under the wrong taxon; on 2026-10-20 NCBI had ~20 orbivirus and bunyavirus
# species (taxids 345xxxx) under Alphavirus, whose VP1 / L proteins hit
# Orbi_VP1 / Bunya_RdRp best and would otherwise have joined the Togaviridae
# tree (whole-genome nucleotide mode took them too).  The lists are every
# profile with a passing hit on the family's RefSeq RdRps in the coverage
# test (below), in order of how many viruses they hit.  RdRP_1–RdRP_4 are
# broad palm-domain profiles that overlap heavily (picornavirus RefSeqs hit
# all four), so a family with any of them gets all four (RDRP_PALM): telling
# them apart is not reliable from a family's few RefSeqs, and it keeps
# divergent metagenomic members of small families.  In the validation runs
# the rule left out the 29 misfiled Togaviridae genomes and nothing else
# (0 of ~8,600 genomes in Hantaviridae, Partitiviridae, Tombusviridae,
# Xinmoviridae).
RDRP_PALM: tuple[str, ...] = ("RdRP_1", "RdRP_2", "RdRP_3", "RdRP_4")

RDRP_FAMILY_PROFILES: dict[str, tuple[str, ...]] = {
    'Aliusviridae': ('Mononeg_RNA_pol',),
    'Alphaflexiviridae': RDRP_PALM,
    'Alternaviridae': RDRP_PALM,
    'Amalgaviridae': RDRP_PALM,
    'Ambiguiviridae': RDRP_PALM,
    'Arenaviridae': ('Arena_RNA_pol', 'Bunya_RdRp'),
    'Arteriviridae': RDRP_PALM,
    'Artiviridae': RDRP_PALM,
    'Artoviridae': ('Mononeg_RNA_pol',),
    'Astroviridae': RDRP_PALM,
    'Barnaviridae': RDRP_PALM,
    'Benyviridae': RDRP_PALM,
    'Betaflexiviridae': RDRP_PALM,
    'Birnaviridae': ('Birna_RdRp_palm', 'Permu_RdRp_palm', *RDRP_PALM),
    'Blumeviridae': ('RNA_replicase_B', 'Mitovir_RNA_pol'),
    'Bornaviridae': ('Mononeg_RNA_pol',),
    'Botourmiaviridae': ('Mitovir_RNA_pol',),
    'Botybirnaviridae': RDRP_PALM,
    'Bromoviridae': RDRP_PALM,
    'Caliciviridae': RDRP_PALM,
    'Carmotetraviridae': RDRP_PALM,
    'Chrysoviridae': RDRP_PALM,
    'Chuviridae': ('Mononeg_RNA_pol',),
    'Closteroviridae': RDRP_PALM,
    'Coronaviridae': ('CoV_RPol_N', *RDRP_PALM),
    'Cremegaviridae': RDRP_PALM,
    'Crepuscuviridae': ('Mononeg_RNA_pol',),
    'Cruliviridae': ('Bunya_RdRp', 'L_thumb_ring_vir'),
    'Curvulaviridae': RDRP_PALM,
    'Deltaflexiviridae': RDRP_PALM,
    'Dicistroviridae': RDRP_PALM,
    'Endornaviridae': RDRP_PALM,
    'Filoviridae': ('Mononeg_RNA_pol',),
    'Fimoviridae': ('Bunya_RdRp',),
    'Flaviviridae': ('Flavi_NS5', 'Flavi_NS5_thumb', *RDRP_PALM),
    'Fusagraviridae': RDRP_PALM,
    'Fusariviridae': RDRP_PALM,
    'Gammaflexiviridae': RDRP_PALM,
    'Giardiaviridae': RDRP_PALM,
    'Gresnaviridae': RDRP_PALM,
    'Hadakaviridae': RDRP_PALM,
    'Hantaviridae': ('Bunya_RdRp', 'Flu_PB1'),
    'Hepaciviridae': RDRP_PALM,
    'Hepeviridae': RDRP_PALM,
    'Iflaviridae': RDRP_PALM,
    'Inseviridae': RDRP_PALM,
    'Kitaviridae': RDRP_PALM,
    'Lebotiviridae': RDRP_PALM,
    'Leishbuviridae': ('Bunya_RdRp',),
    'Lispiviridae': ('Mononeg_RNA_pol',),
    'Marnaviridae': RDRP_PALM,
    'Matonaviridae': RDRP_PALM,
    'Mayoviridae': RDRP_PALM,
    'Megabirnaviridae': RDRP_PALM,
    'Megatotiviridae': RDRP_PALM,
    'Mitoviridae': ('Mitovir_RNA_pol', 'RNA_replicase_B'),
    'Monocitiviridae': RDRP_PALM,
    'Mymonaviridae': ('Mononeg_RNA_pol',),
    'Myriaviridae': ('Mononeg_RNA_pol',),
    'Nairoviridae': ('Bunya_RdRp',),
    'Nanghoshaviridae': RDRP_PALM,
    'Nanhypoviridae': RDRP_PALM,
    'Natareviridae': ('Mononeg_RNA_pol',),
    'Nodaviridae': RDRP_PALM,
    'Noraviridae': RDRP_PALM,
    'Nyamiviridae': ('Mononeg_RNA_pol',),
    'Olifoviridae': RDRP_PALM,
    'Ootiviridae': RDRP_PALM,
    'Orthomyxoviridae': ('Flu_PB1', 'Bunya_RdRp'),
    'Orthototiviridae': RDRP_PALM,
    'Parahypoviridae': RDRP_PALM,
    'Paramyxoviridae': ('Mononeg_RNA_pol',),
    'Partitiviridae': RDRP_PALM,
    'Peribunyaviridae': ('Bunya_RdRp', 'L_thumb_ring_vir', 'Mononeg_RNA_pol'),
    'Permutotetraviridae': ('Permu_RdRp_palm', 'Birna_RdRp_palm'),
    'Pestiviridae': ('Flavi_NS5_thumb', 'Flavi_NS5', *RDRP_PALM),
    'Phenuiviridae': ('Bunya_RdRp', 'L_PA-C-like', 'Flu_PB1'),
    'Phlegiviridae': RDRP_PALM,
    'Picobirnaviridae': ('Mitovir_RNA_pol', *RDRP_PALM),
    'Picornaviridae': RDRP_PALM,
    'Pistolviridae': RDRP_PALM,
    'Pneumoviridae': ('Mononeg_RNA_pol',),
    'Polycipiviridae': RDRP_PALM,
    'Polymycoviridae': RDRP_PALM,
    'Potyviridae': RDRP_PALM,
    'Pseudototiviridae': ('Permu_RdRp_palm', *RDRP_PALM),
    'Quadriviridae': RDRP_PALM,
    'Rhabdoviridae': ('Mononeg_RNA_pol',),
    'Roniviridae': RDRP_PALM,
    'Secoviridae': ('Permu_RdRp_palm', *RDRP_PALM),
    'Sedoreoviridae': ('Orbi_VP1', 'Seadorna_RNAP', 'RNA-dir_pol_phytoreovirus', *RDRP_PALM),
    'Sinhaliviridae': RDRP_PALM,
    'Solemoviridae': RDRP_PALM,
    'Solinviviridae': RDRP_PALM,
    'Spiciviridae': RDRP_PALM,
    'Spinareoviridae': ('CPV_RdRP_pol_dom', 'RdRP_5', *RDRP_PALM),
    'Steitzviridae': ('RNA_replicase_B', 'Mitovir_RNA_pol'),
    'Sunviridae': ('Mononeg_RNA_pol',),
    'Tobaliviridae': RDRP_PALM,
    'Tobaniviridae': ('CoV_RPol_N', *RDRP_PALM),
    'Togaviridae': RDRP_PALM,
    'Tombusviridae': ('Luteo_P1-P2', *RDRP_PALM),
    'Tospoviridae': ('Bunya_RdRp',),
    'Tymoviridae': RDRP_PALM,
    'Virgaviridae': RDRP_PALM,
    'Xinmoviridae': ('Mononeg_RNA_pol',),
    'Yadokariviridae': RDRP_PALM,
}

RDRP_FAMILIES: frozenset[str] = frozenset(RDRP_FAMILY_PROFILES)

# Profiles for other domains of the same polymerase as a family profile
# (hmm.best_hit_also).  A protein whose best hit is one of these is still the
# family's own kind of RdRp, so best_hit_only must not leave it out: 6 of 149
# Peribunyaviridae RefSeq L proteins fit L_thumb_ring_vir better than
# Bunya_RdRp, while Hantaviridae and Nairoviridae RefSeqs never hit it, so a
# divergent hanta- or nairovirus L could otherwise be dropped as "foreign".
RDRP_SIBLING_DOMAINS: dict[str, tuple[str, ...]] = {
    "Bunya_RdRp": ("L_PA-C-like", "L_thumb_ring_vir"),       # bunyavirus L
    "Arena_RNA_pol": ("L_PA-C-like", "L_thumb_ring_vir"),    # arenavirus L
    "Flavi_NS5": ("Flavi_NS5_thumb",),                       # flavivirus NS5
    "Flavi_NS5_thumb": ("Flavi_NS5",),
}

# Sibling domains outside the RdRp core (L_PA-C-like, L_thumb_ring_vir).
# They never become marker tokens, even where a family's RefSeq L proteins
# hit them (Peribunyaviridae, Phenuiviridae, Cruliviridae): a token selects
# and crops on its own, so a partial L covering only such a domain would
# enter the alignment cut to it, with no RdRp columns.  They only count for
# best_hit_only (best_hit_also).  Checked 2026-10-21 on the coverage scan:
# every RefSeq L protein of those families also hits Bunya_RdRp (151, 154
# and 3 proteins), and the median cut shrinks from ~1300 to ~900 aa (the
# Bunya_RdRp domain ± 100).
RDRP_NON_CORE_DOMAINS: frozenset[str] = frozenset(
    s for sibs in RDRP_SIBLING_DOMAINS.values() for s in sibs
) - frozenset(RDRP_SIBLING_DOMAINS)

# Proteins shorter than this are never marker candidates in RdRp mode
# (hmm.min_protein_length).  Equal to quality.MIN_LENGTH_AA: a shorter protein
# yields a shorter cropped domain, which quality control drops anyway, so the
# floor changes no tree.  It is applied to the downloaded records, not in the
# NCBI query: a [SLEN] clause made each protein search 4-6 s slower
# (2026-10-21), against 2-14 % fewer records to download.  A higher floor
# would also drop partial RdRp deposits (3-8 % of the selected proteins were
# < 200 aa).
RDRP_MIN_PROTEIN_LENGTH = MIN_LENGTH_AA


# ---------------------------------------------------------------------------
# Reverse-transcriptase and Rep markers (v1.6.0)
# ---------------------------------------------------------------------------
# The same approach for the two groups without an RdRp whose whole-genome
# trees were saturated or weak in the 2026-10-20 signal survey.  Bundled
# Pfam sets "RT" and "Rep" (names verified against InterPro 2026-10-21); hit
# rule, length floor and opt-out as for the RdRp.  Coverage of RefSeq
# viruses under that rule (2026-10-21, scratch coverage test):
#
# RT (reverse-transcribing viruses, Pararnavirae): RVT_1 in Caulimoviridae
# 108/108, Hepadnaviridae 27/27, Metaviridae 2/2, Retroviridae 64/89 — every
# uncovered retrovirus has no RT gene (transforming viruses with an onc gene
# in place of pol, endogenous-retrovirus fragments), so they leave the tree.
# The RT is cut to its domain ± 100 aa, as Gag-Pol, the hepadnavirus P
# protein and caulimovirus polyproteins are 800-1850 aa.  RVT_thumb and
# RVT_connect (other domains of the same RT) only count for best_hit_only;
# RVT_2 (Ty1/copia RT) is in the set so such a protein is recognised as
# foreign.  Pseudoviridae (0 RefSeq proteins) stays as it was.
RT_HMM_TOKENS: list[str] = ["RVT_1", "RVT_2", "RVT_thumb", "RVT_connect"]
RT_FAMILY_PROFILES: dict[str, tuple[str, ...]] = {
    "Caulimoviridae": ("RVT_1",),
    "Hepadnaviridae": ("RVT_1",),
    "Metaviridae": ("RVT_1",),
    "Retroviridae": ("RVT_1",),
}
RT_SIBLING_DOMAINS: tuple[str, ...] = ("RVT_thumb", "RVT_connect")

# Rep (CRESS-DNA viruses, Cressdnaviricota, and Alphasatellitidae, whose
# satellites encode a circovirus-like Rep): >= 94 % of RefSeq viruses in all
# 23 families.  Two kinds of Rep: geminivirus-like (Gemini_AL1 / AL1_M best;
# Geminiviridae, Genomoviridae, Geplanaviridae) and circovirus-like
# (Viral_Rep + SF3 helicase; the rest), but the profiles do not tell families
# apart, so every Rep profile is a token and best_hit_only is off.  Rep
# proteins are ~270-470 aa: kept whole (no cut).  Not here: Smacoviridae
# (38/88; its name-based Rep preset stays), Metaxyviridae (no RefSeq
# proteins), Tolecusatellitidae (betasatellites encode no Rep).
REP_HMM_TOKENS: list[str] = ["Viral_Rep", "RNA_helicase", "Gemini_AL1", "Gemini_AL1_M"]
REP_FAMILIES: frozenset[str] = frozenset({
    "Adamaviridae", "Alphasatellitidae", "Amesuviridae", "Anicreviridae",
    "Bacilladnaviridae", "Circoviridae", "Draupnirviridae", "Endolinaviridae",
    "Gandrviridae", "Geminiviridae", "Genomoviridae", "Geplanaviridae",
    "Kanorauviridae", "Kirkoviridae", "Mahapunaviridae", "Nanoviridae",
    "Naryaviridae", "Nenyaviridae", "Oomyviridae", "Ouroboviridae",
    "Pecoviridae", "Redondoviridae", "Vilyaviridae",
})

#: Every family on an HMM-identified marker preset → the marker (also the
#: bundled profile set's name and the ``sequence.region`` it runs under).
HMM_MARKER_FAMILIES: dict[str, str] = {
    **{f: "RdRp" for f in RDRP_FAMILIES},
    **{f: "RT" for f in RT_FAMILY_PROFILES},
    **{f: "Rep" for f in REP_FAMILIES},
}


def _hmm_marker_preset(
    marker: str, tokens: list[str], *, padding: int | None,
    best_hit_only: bool, also: list[str],
) -> dict:
    """Config overrides for one HMM-identified marker (shared hit rule)."""
    return {
        "sequence": {
            "type": "protein",
            "region": marker,
            "segment": None,   # the HMM finds the marker on whichever segment carries it
            "hmms": list(tokens),
            "domain_padding": padding,
        },
        "hmm": {
            "enabled": True,
            "database": marker,
            "use_ga_when_available": False,
            "default_evalue": 1.0e-3,
            "relative_length_cutoff": 0.2,
            "best_hit_only": best_hit_only,
            "best_hit_also": list(also),
            "min_protein_length": RDRP_MIN_PROTEIN_LENGTH,
        },
    }


def rdrp_preset(family: str) -> dict:
    """Config overrides that put *family* (in RDRP_FAMILIES) on the
    HMM-identified RdRp marker, restricted to the family's own profiles."""
    profiles = RDRP_FAMILY_PROFILES[family]
    tokens = [p for p in profiles if p not in RDRP_NON_CORE_DOMAINS] or list(profiles)
    siblings = sorted((
        {s for p in tokens for s in RDRP_SIBLING_DOMAINS.get(p, ())}
        | (set(profiles) - set(tokens))
    ) - set(tokens))
    return _hmm_marker_preset(
        "RdRp", tokens, padding=RDRP_DOMAIN_PADDING, best_hit_only=True, also=siblings,
    )


def hmm_marker_preset(family: str) -> dict:
    """The preset of a family in HMM_MARKER_FAMILIES (RdRp, RT or Rep)."""
    marker = HMM_MARKER_FAMILIES[family]
    if marker == "RdRp":
        return rdrp_preset(family)
    if marker == "RT":
        return _hmm_marker_preset(
            "RT", list(RT_FAMILY_PROFILES[family]), padding=RDRP_DOMAIN_PADDING,
            best_hit_only=True, also=list(RT_SIBLING_DOMAINS),
        )
    return _hmm_marker_preset(
        "Rep", REP_HMM_TOKENS, padding=None, best_hit_only=False, also=[],
    )


def hmm_marker_opt_out(file_cfg: dict | None, marker: str) -> bool:
    """True when a family file chooses its own analysis instead of its HMM
    marker preset (*marker*: RdRp, RT or Rep): a ``sequence.type`` other than
    protein, a ``sequence.region`` other than *marker*, or
    ``hmm.enabled: false``.

    The preset is then not applied at all, so the family gets the template
    (and SEGMENTED_FAMILIES segment) underneath the file's own values — never
    a mix such as nucleotide type with the preset's HMM protein fetch.  A
    v1.5.0 file (whole-genome nucleotide, hmm off) therefore runs as it did.
    """
    file_seq = (file_cfg or {}).get("sequence") or {}
    file_hmm = (file_cfg or {}).get("hmm") or {}
    seq_type = file_seq.get("type")
    region = file_seq.get("region")
    return (
        (seq_type is not None and seq_type != "protein")
        or (region is not None and region != marker)
        or file_hmm.get("enabled") is False
    )


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
# Curation status: first pass from the literature; Asfarviridae and
# Pithoviridae moved to HMM-identified markers in 1.6.2.  See
# CONCAT_DESIGN.md §4 for references.
# NCLDV core markers identified by HMM (bundled set "NCLDV"), shared by the
# Asfarviridae and Pithoviridae concatenations.  Tokens of one marker are
# OR'd: the capping enzyme carries both domains, and some genomes hit only one.
_NCLDV_HMM_MARKERS: list[dict] = [
    {"name": "DNA polymerase B", "aliases": ["DNA polymerase"], "hmms": ["DNA_pol_B"]},
    {"name": "major capsid protein", "aliases": ["p72"], "hmms": ["Capsid_NCLDV"]},
    {"name": "DNA-directed RNA polymerase largest subunit", "aliases": ["RPB1"],
     "hmms": ["RNA_pol_Rpb1_2"]},
    {"name": "DNA-directed RNA polymerase second subunit", "aliases": ["RPB2"],
     "hmms": ["RNA_pol_Rpb2_6"]},
    {"name": "primase-helicase", "aliases": ["D5-like helicase"], "hmms": ["D5_N"]},
    {"name": "late transcription factor 3", "aliases": ["VLTF-3"], "hmms": ["Pox_VLTF3"]},
    {"name": "mRNA capping enzyme", "aliases": ["capping enzyme"],
     "hmms": ["mRNA_G-N7_MeTrfase", "mRNA_cap_enzyme"]},
    {"name": "DNA topoisomerase II", "aliases": ["topoisomerase II"], "hmms": ["DNA_topoisoIV"]},
]

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
                    "hmms": ["DNA_pol_B"],                  # PF00136
                },
                {
                    "name": "DNA-directed RNA polymerase 147 kDa subunit",
                    "aliases": ["RPO147", "RNA polymerase subunit RPO147"],
                    "locus_tag_hint": r"A24R|RPO147|rpo147",
                    "hmms": ["RNA_pol_Rpb1_2"],             # PF00623
                },
                {
                    "name": "DNA-directed RNA polymerase 132 kDa subunit",
                    "aliases": ["RPO132", "RNA polymerase subunit RPO132"],
                    "locus_tag_hint": r"J6R|RPO132|rpo132",
                    "hmms": ["RNA_pol_Rpb2_6"],             # PF00562
                },
                {
                    "name": "mRNA capping enzyme large subunit",
                    "aliases": ["capping enzyme large subunit"],
                    "locus_tag_hint": r"D1R",
                    # PF03291 (the cap N7-methyltransferase domain carried by the
                    # D1 large subunit).  The generic capping families PF01331
                    # (catalytic) and PF03919 (C-term) do NOT match poxvirus D1
                    # at all — a 2026-06 live run scored 0/499 coverage with
                    # PF01331; PF03291 recovers it (score ~58 vs GA ~27).
                    "hmms": ["mRNA_G-N7_MeTrfase"],        # PF03291
                },
                {
                    "name": "DNA helicase",
                    "aliases": ["NPH-II", "transcript release factor"],
                    "locus_tag_hint": r"A18R|NPH2",
                    "hmms": ["NPH-II"],                     # PF12011
                },
                {
                    "name": "poly(A) polymerase catalytic subunit",
                    "aliases": ["poly(A) polymerase large subunit"],
                    "locus_tag_hint": r"E1L|VP55",
                    "hmms": ["Pox_polyA_pol"],              # PF03296
                },
                {
                    "name": "late transcription factor VLTF-3",
                    "aliases": ["VLTF3", "late transcription factor 3"],
                    "locus_tag_hint": r"A1L|VLTF",
                    "hmms": ["Pox_VLTF3"],                  # PF04947
                },
                {
                    "name": "uracil-DNA glycosylase",
                    "aliases": ["UNG"],
                    "locus_tag_hint": r"D4R|UNG",
                    "hmms": ["UDG"],                        # PF03167
                },
                # ssDNA-binding protein (I3L) was previously included but a
                # 2026-05 cache audit (128 species) showed 1 % coverage —
                # the I3L name is essentially VACV-Cop-specific and not used
                # outside Chordopoxvirinae reference annotations.
            ],
        },
        # relative_length_cutoff lowered to 0.40 (from the 0.50 default): the
        # divergent Avipoxvirus capping enzyme (D1) matches mRNA_G-N7_MeTrfase
        # with a strong score (~53 vs GA ~27) but covers only 47 % of the HMM
        # model, so the 0.50 gate dropped it; 0.40 admits it.  The GA bit-score
        # gate still filters non-homologs, and the 6 universal markers are
        # high-coverage so they are unaffected.
        "hmm": {
            "enabled": True, "database": "Poxviridae",
            "relative_length_cutoff": 0.40,
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

    # ---- Asfarviridae (8 markers, HMM-identified; 2026-10-27) -----------
    # Refs: Yutin & Koonin 2012; Iyer et al. 2006.  Until 1.6.2 the markers
    # were found by GenBank name, and ASFV names only its DNA polymerase and
    # p72 that way: every genome failed min_fraction and the family had no
    # tree.  Profiles from the shared NCLDV set, chosen from the Pfam domains
    # UniProt gives the ASFV BA71V proteins and checked on the pipeline's own
    # fetch (2026-10-27): all 35 genome-scale records of ASFV (6000-protein
    # ceiling) and of Abalone asfa-like virus keep ≥ 6 of 8.  The packaging
    # ATPase (Pox_A32) hits only the asfa-like virus and is left out.
    # relative_length_cutoff 0.30: HMMER splits the polymerase domain of
    # ASFV (hits of 0.47 / 0.45 of the profile) and of the asfa-like virus
    # (0.31 / 0.23 / 0.36) at an insertion, so no single hit reaches 0.50.
    "Asfarviridae": {
        "sequence": {"region": "concatenated", "type": "protein"},
        "concatenation": {"proteins": copy.deepcopy(_NCLDV_HMM_MARKERS)},
        "hmm": {"enabled": True, "database": "NCLDV", "relative_length_cutoff": 0.30},
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
    # Pitho/Pandoraviridae and any future NCLDV-like family).
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
CONCATENATION_FAMILIES["Pandoraviridae"] = {
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

# Pithoviridae — 7 HMM-identified markers (2026-10-27; name-based until
# 1.6.2, when 2 genomes passed min_fraction and the family had no tree).
# Every complete genome (18 records of 7 species) carries all 7; pithoviruses
# have no Pfam-detectable major capsid protein, so Capsid_NCLDV is left out.
# The metagenome-assembled genomes (LCPAC…) are split over contigs with at
# most 4 markers each and stay out, as before.  Cutoff as for Asfarviridae.
CONCATENATION_FAMILIES["Pithoviridae"] = {
    "sequence": {"region": "concatenated", "type": "protein"},
    "concatenation": {
        "proteins": [
            copy.deepcopy(p) for p in _NCLDV_HMM_MARKERS
            if p["name"] != "major capsid protein"
        ],
    },
    "hmm": {"enabled": True, "database": "NCLDV", "relative_length_cutoff": 0.30},
}


DEFAULT_FAMILY_CONFIG: dict = {
    "download": {
        "max_per_species": 300,
    },
    "sequence": {
        "type": "nucleotide",
        "region": "whole_genome",
        "segment": None,   # set automatically for known segmented families
        # HMM profile tokens of the single marker (hmm.enabled); set by the
        # RdRp, RT, Rep and Adenoviridae presets.  Empty: no profile gate.
        "hmms": [],
        # HMM single-marker mode (hmm.enabled with sequence.hmms tokens, set by
        # the RdRp, RT, Rep and Adenoviridae presets): keep only the span of the
        # marker's profile hits ± this many residues (e.g. the RdRp of a
        # polyprotein); null keeps whole proteins.
        "domain_padding": None,
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
        # --adjustdirection (nucleotide only): reverse-complements GenBank
        # records deposited on the opposite strand.  The MSA wrapper strips
        # MAFFT's "_R_" name prefix and the pipeline logs the flipped ids.
        "options_nuc": "--retree 2 --adjustdirection",
        "options_aa": "--auto",
    },
    "msa_100": {
        "tool": "mafft",
        "options_nuc": "--retree 3 --adjustdirection",
        "options_aa": "--maxiterate 1000 --localpair",
    },
    "msa_trim": {
        "enabled": True,
        "tool": "trimal",
        "options": DEFAULT_TRIM_OPTIONS,
    },
    "tree_500": {
        # Nucleotide: IQ-TREE --fast, same wall time as FastTree with a
        # clearly better likelihood (CHANGELOG 1.5.0); a full search gained
        # nothing there.  Protein: -n 50, IQ-TREE's normal search stopped
        # after 50 iterations — --fast lands in a different local optimum on
        # each seed for weak-signal protein sets, -n 50 comes as close to the
        # full search as a second full search does (CHANGELOG 1.6.0).  The
        # wrapper adds -alrt 1000, so the tree carries SH-aLRT support.
        # FastTree is still available: tool: fasttree with options_nuc /
        # options_aa set to "" (config validation rejects IQ-TREE options).
        "tool": "iqtree",
        "options_nuc": "--fast",
        "options_aa": "-n 50",
        "model_nuc": "GTR+G",
        "model_aa": "LG+G",
    },
    "tree_100": {
        "tool": "iqtree",
        # Full tree search with SH-aLRT and UFBoot; --bnni curbs UFBoot
        # overestimation under model violation.  Outputs carry the UFBoot
        # value; the SH-aLRT value stays in the IQ-TREE .treefile in _work.
        "options_nuc": "-alrt 1000 -B 1000 --bnni",
        "options_aa": "-alrt 1000 -B 1000 --bnni",
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
        # Check each round on a quick IQ-TREE --fast tree first; the
        # configured (slower) search runs once that finds nothing to remove,
        # and its tree is checked as before.  Off: every round runs the
        # configured search.
        "quick_check": True,
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
        # Write style:font_color (the genus palette) on every PhyloXML tip.
        # Off by default: Archaeopteryx and Archaeopteryx.js open a tree with
        # their own Color-by, which hides style colours, and these have no
        # legend.  PDF/PNG/Auspice keep the palette either way.
        "phyloxml_font_colors": False,
        # Write {prefix}_tree_{500,100}.html next to each PhyloXML: the tree
        # inlined in a page that opens it in Archaeopteryx.js (the viewer is
        # loaded from cdn.jsdelivr.net, so it needs an internet connection).
        "archaeopteryx_html": True,
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
    # HMM-based marker-protein identification (HMMER hmmscan + curated
    # family-specific Pfam-A profiles).  When enabled, markers are identified
    # by profile homology instead of GenBank protein names — far more robust
    # for large DNA viruses whose annotation is inconsistent.  Applies to both
    # the single-marker (DNA_FAMILIES) and concatenated-marker
    # (CONCATENATION_FAMILIES) paths via the MarkerIdentifier swap
    # (markers.make_identifier); a marker is HMM-gated when its spec carries an
    # 'hmms' token list (single "Name" or multidomain "A--B--C" in N-to-C order).
    #   enabled:        master switch (default false → name + alias matching).
    #   database:       a bundled per-family set NAME (e.g. "Poxviridae" →
    #                   vfam_trees/data/hmms/Poxviridae/), a directory of .hmm
    #                   files, or a single .hmm file path.  null when disabled.
    #   default_evalue: domain E-value gate for profiles lacking a curated GA.
    #   use_ga_when_available: prefer each profile's curated Pfam GA bit-score.
    #   relative_length_cutoff: minimum fraction of the HMM model a hit must
    #                   cover (in HMM-model coordinates), in [0, 1].
    #   multidomain_overlap_tolerance: residues two adjacent domains of a
    #                   multidomain "A--B" token may overlap (cleavage-site slack).
    #   best_hit_only:  a protein counts only when its best-scoring profile in
    #                   the database is one of the marker's tokens, so a protein
    #                   that merely cross-hits a token while another profile
    #                   fits it better is left out (on in the RdRp preset).
    #   best_hit_also:  profiles that also count as the marker's own for
    #                   best_hit_only without being tokens (other domains of
    #                   the same protein, e.g. the L thumb-ring domain).
    #   min_protein_length: proteins shorter than this (aa) are left out of
    #                   the downloaded candidates in HMM mode; null → no floor
    #                   (100 in the RdRp preset).
    #   threads:        hmmscan --cpu; null → use the run's thread count.
    "hmm": {
        "enabled": False,
        "database": None,
        "default_evalue": 1.0e-5,
        "use_ga_when_available": True,
        "relative_length_cutoff": 0.5,
        "multidomain_overlap_tolerance": 30,
        "best_hit_only": False,
        "best_hit_also": [],
        "min_protein_length": None,
        "threads": None,
    },
}


# Settings read from global.yaml outside its ``defaults:`` block, with the
# values the pipeline uses when a key is absent.  Single source of truth for
# those fallbacks (``global_settings``), the unknown-key check and the run
# record (``effective_global_config``).
GLOBAL_SETTINGS_DEFAULTS: dict = {
    "ncbi": {"email": "", "api_key": ""},
    # null → virus_families_annotation.tsv next to global.yaml
    "annotation_tsv": None,
    "cache": {"dir": None, "ttl_days": None},
    "visualization": {"branch_linewidth": 0.5},
    "icon": {"size": 256, "bg_color": "#EAF3F2", "branch_color": "#000000"},
}

# Keys older templates wrote that nothing reads, with where the setting lives.
RETIRED_GLOBAL_KEYS: dict[str, str] = {
    "output_dir": "the output folder is set with `vfam_trees run -o`",
    "log_level": "the log level is set with `vfam_trees run --log-level` "
                 "(or --verbose / --quiet)",
}


def global_settings(global_cfg: dict, block: str) -> dict:
    """One global.yaml block (``visualization``, ``icon``, ``cache``, ``ncbi``)
    with every key present: the file's values over GLOBAL_SETTINGS_DEFAULTS."""
    return {**GLOBAL_SETTINGS_DEFAULTS[block], **((global_cfg or {}).get(block) or {})}


def unknown_config_keys(cfg: dict, reference: dict, prefix: str = "") -> list[str]:
    """Dotted paths of keys in *cfg* that *reference* does not have.

    Recurses into blocks whose reference is a non-empty mapping; list values
    (marker specs, sources) are not checked.  Keys starting with "_" are
    internal and skipped."""
    found: list[str] = []
    for key, value in (cfg or {}).items():
        if not isinstance(key, str) or key.startswith("_"):
            continue
        if key not in reference:
            found.append(prefix + key)
        elif isinstance(value, dict) and isinstance(reference[key], dict) and reference[key]:
            found.extend(unknown_config_keys(value, reference[key], f"{prefix}{key}."))
    return found


def unknown_keys_message(path: Path | str, keys: list[str]) -> str:
    """WARNING text for keys a config file sets that this version ignores."""
    parts = []
    for key in keys:
        hint = RETIRED_GLOBAL_KEYS.get(key)
        parts.append(f"{key} (not read: {hint})" if hint else key)
    return (
        f"Unknown config key(s) in {path}: {', '.join(parts)}. This version "
        "does not read them, so they have no effect; remove them, or check the "
        "spelling against `vfam_trees init --full`."
    )


def _global_reference() -> dict:
    return {**GLOBAL_SETTINGS_DEFAULTS, "defaults": DEFAULT_FAMILY_CONFIG}


def load_global_config(path: Path) -> dict:
    """Read and validate global.yaml.  Keys this version does not know are
    logged (WARNING) and kept in ``cfg["_unknown_keys"]`` so the pipeline can
    repeat them in each family log."""
    if not path.exists():
        raise FileNotFoundError(f"Global config not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    _validate_global_config(cfg)
    unknown = unknown_config_keys(cfg, _global_reference())
    if unknown:
        log.warning(unknown_keys_message(path, unknown))
    cfg["_unknown_keys"] = unknown
    return cfg


def effective_global_config(global_cfg: dict, redact: bool = True) -> dict:
    """Every global setting a run uses: GLOBAL_SETTINGS_DEFAULTS and the
    template (``defaults:``) overlaid with the file's values.  Unknown keys
    are kept as given (they are flagged elsewhere); ``ncbi.api_key`` is
    replaced by a marker when *redact* so a results folder can be shared.
    Family presets and family files are applied later, per family."""
    cfg = {k: v for k, v in (global_cfg or {}).items() if not str(k).startswith("_")}
    out = copy.deepcopy(GLOBAL_SETTINGS_DEFAULTS)
    _deep_update(out, {k: v for k, v in cfg.items() if k != "defaults"})
    defaults = copy.deepcopy(DEFAULT_FAMILY_CONFIG)
    _deep_update(defaults, cfg.get("defaults") or {})
    out["defaults"] = defaults
    if redact and out.get("ncbi", {}).get("api_key"):
        out["ncbi"]["api_key"] = "<set, not recorded>"
    return out


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


# Options only one of the two tree tools understands.  Each list holds
# spellings the other tool does not have or reads differently (IQ-TREE's -n
# is the number of search iterations, FastTree's the number of alignments in
# the input), so a hit is always a mistake — typically a tool switched
# without clearing the options inherited from the template (whose
# options_nuc / options_aa keys always win over a bare legacy ``options`` key).
_IQTREE_ONLY_OPTIONS = frozenset({
    "--fast", "-alrt", "--alrt", "-abayes", "--abayes", "-lbp", "--lbp",
    "-B", "--ufboot", "-bb", "-b", "--boot", "--bnni", "-m", "-T", "-p",
    "-te", "--tree-fix", "--seqtype", "--prefix", "-s", "-n",
})
_FASTTREE_ONLY_OPTIONS = frozenset({
    "-gtr", "-gamma", "-wag", "-lg", "-cat", "-nocat", "-spr", "-mlacc",
    "-slownni", "-fastest", "-pseudo", "-nosupport", "-noml", "-nome",
    "-mllen", "-bionj", "-topm", "-close", "-refresh", "-gtrrates",
    "-gtrfreqs", "-intree", "-constraints",
})
_TREE_TOOLS = {"fasttree": "FastTree", "iqtree": "IQ-TREE", "iqtree2": "IQ-TREE"}


def _validate_tree_blocks(cfg: dict, family: str) -> None:
    """Check tree_500 / tree_100: a known tool, and options that tool accepts.

    Only the options key the run will use is checked (``options_aa`` for
    protein and concatenated runs, else ``options_nuc``; the legacy bare
    ``options`` key when neither is set), mirroring
    ``pipeline._resolve_tree_options``.  Raises ``ValueError`` at load time
    instead of failing after the fetch and alignment.  Also checks
    ``outlier_removal.quick_check`` (a boolean), which picks the tree the
    outlier rounds are checked on.
    """
    seq = cfg.get("sequence") or {}
    protein = seq.get("type") == "protein" or seq.get("region") == "concatenated"
    key = "options_aa" if protein else "options_nuc"
    for label in ("500", "100"):
        block_name = f"tree_{label}"
        block = cfg.get(block_name) or {}
        tool = str(block.get("tool", "")).lower().replace("-", "").replace("_", "")
        if tool not in _TREE_TOOLS:
            raise ValueError(
                f"{family}: {block_name}.tool must be 'iqtree' or 'fasttree' "
                f"(got {block.get('tool')!r})."
            )
        used_key = key if key in block else "options"
        options = block.get(used_key) or ""
        if not isinstance(options, str):
            raise ValueError(
                f"{family}: {block_name}.{used_key} must be a string of "
                f"command-line options (got {type(options).__name__})."
            )
        wrong = _IQTREE_ONLY_OPTIONS if tool == "fasttree" else _FASTTREE_ONLY_OPTIONS
        bad = [o for o in options.split() if o in wrong]
        if bad:
            other = "IQ-TREE" if tool == "fasttree" else "FastTree"
            raise ValueError(
                f"{family}: {block_name}.tool is {block.get('tool')!r} but "
                f"{block_name}.{used_key} has {other} option(s) {' '.join(bad)}, "
                f"which {_TREE_TOOLS[tool]} rejects or reads differently. When switching tools, set "
                f"{block_name}.options_nuc and {block_name}.options_aa as well "
                f"(\"\" for none)."
            )
    quick = (cfg.get("outlier_removal") or {}).get("quick_check", True)
    if not isinstance(quick, bool):
        raise ValueError(
            f"{family}: outlier_removal.quick_check must be true or false "
            f"(got {quick!r})."
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

    from .quality import sequence_type_mismatch  # local: keeps config import light
    seq_type = (cfg.get("sequence") or {}).get("type", "nucleotide")
    for i, entry in enumerate(cleaned_fasta):
        problem = sequence_type_mismatch(entry["sequence"], seq_type)
        if problem:
            raise ValueError(
                f"{family}: manual.include_seq[{i}] ({entry['id']!r}) {problem}."
            )

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
        seq_cfg = cfg.get("sequence") or {}
        region = seq_cfg.get("region", "")
        if region == "concatenated":
            raise ValueError(
                f"{family}: additional_data_sources are not supported when "
                "sequence.region is 'concatenated' — outbreak nucleotide "
                "sequences cannot be split into the configured marker proteins.  "
                "Remove the additional_data_sources entries or switch the region."
            )
        if seq_cfg.get("type") == "protein":
            marker = HMM_MARKER_FAMILIES.get(family)
            preset = f" (its HMM-identified {marker} marker preset)" if marker else ""
            raise ValueError(
                f"{family}: additional_data_sources inject nucleotide genomes, "
                f"but this family is analysed as protein{preset}.  For an "
                f"outbreak tree set sequence.type: nucleotide and hmm.enabled: "
                f"false in {family}.yaml (plus sequence.region / segment as "
                "needed), or remove the additional_data_sources entries."
            )
    cfg["additional_data_sources"] = cleaned


def _validate_hmm_block(cfg: dict, family: str) -> None:
    """Validate the optional ``hmm`` block and any marker ``hmms`` tokens.

    Fills defaults so the pipeline reads a fully-typed block.  No network or
    filesystem I/O — the database is resolved (and hmmscan invoked) only at
    runtime.  When ``hmm.enabled`` is true a ``database`` must be set, and every
    ``hmms`` token (on concatenation markers or the single-marker ``sequence``
    block) must be well-formed (``parse_hmm_token`` rejects empty components).
    """
    from .hmm import parse_hmm_token  # local import: cheap, avoids any cycle

    raw = cfg.get("hmm")
    if raw is None:
        raw = dict(DEFAULT_FAMILY_CONFIG["hmm"])
    if not isinstance(raw, dict):
        raise ValueError(
            f"{family}: hmm must be a mapping (got {type(raw).__name__})."
        )

    enabled = raw.get("enabled", False)
    if not isinstance(enabled, bool):
        raise ValueError(f"{family}: hmm.enabled must be a boolean.")

    database = raw.get("database")
    if database is not None and (
        not isinstance(database, str) or not database.strip()
    ):
        raise ValueError(
            f"{family}: hmm.database must be null or a non-empty string "
            "(a bundled set name, a directory of .hmm files, or a .hmm path)."
        )

    ev = raw.get("default_evalue", 1.0e-5)
    if isinstance(ev, bool) or not isinstance(ev, (int, float)) or ev <= 0:
        raise ValueError(f"{family}: hmm.default_evalue must be a positive number.")

    rel = raw.get("relative_length_cutoff", 0.5)
    if isinstance(rel, bool) or not isinstance(rel, (int, float)) or not (
        0.0 <= rel <= 1.0
    ):
        raise ValueError(
            f"{family}: hmm.relative_length_cutoff must be a number in [0, 1]."
        )

    use_ga = raw.get("use_ga_when_available", True)
    if not isinstance(use_ga, bool):
        raise ValueError(f"{family}: hmm.use_ga_when_available must be a boolean.")

    tol = raw.get("multidomain_overlap_tolerance", 30)
    if isinstance(tol, bool) or not isinstance(tol, int) or tol < 0:
        raise ValueError(
            f"{family}: hmm.multidomain_overlap_tolerance must be a "
            "non-negative integer."
        )

    best_only = raw.get("best_hit_only", False)
    if not isinstance(best_only, bool):
        raise ValueError(f"{family}: hmm.best_hit_only must be a boolean.")

    also = raw.get("best_hit_also") or []
    if not isinstance(also, list) or not all(
        isinstance(a, str) and a.strip() for a in also
    ):
        raise ValueError(
            f"{family}: hmm.best_hit_also must be a list of profile names."
        )

    min_len = raw.get("min_protein_length")
    if min_len is not None and (
        isinstance(min_len, bool) or not isinstance(min_len, int) or min_len <= 0
    ):
        raise ValueError(
            f"{family}: hmm.min_protein_length must be null or a positive integer."
        )

    threads = raw.get("threads")
    if threads is not None and (
        isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0
    ):
        raise ValueError(
            f"{family}: hmm.threads must be null or a positive integer."
        )

    if enabled and not database:
        raise ValueError(
            f"{family}: hmm.enabled is true but hmm.database is not set. Point "
            "it at a bundled family set name, a directory of .hmm files, or a "
            ".hmm file."
        )

    def _check_tokens(tokens, where: str) -> None:
        # Tokens are validated even when hmm is disabled, so a typo surfaces
        # regardless of the master switch.
        if tokens is None:
            return
        if not isinstance(tokens, list):
            raise ValueError(f"{family}: {where} must be a list of token strings.")
        for tok in tokens:
            if not isinstance(tok, str):
                raise ValueError(
                    f"{family}: {where} entries must be strings "
                    f"(got {type(tok).__name__})."
                )
            try:
                parse_hmm_token(tok)
            except ValueError as exc:
                raise ValueError(f"{family}: {where} token {tok!r}: {exc}") from exc

    _check_tokens((cfg.get("sequence") or {}).get("hmms"), "sequence.hmms")
    padding = (cfg.get("sequence") or {}).get("domain_padding")
    if padding is not None and (isinstance(padding, bool) or not isinstance(padding, int) or padding < 0):
        raise ValueError(
            f"{family}: sequence.domain_padding must be null (keep whole proteins) "
            f"or a non-negative integer (residues kept on each side of the marker domain)."
        )
    proteins = ((cfg.get("concatenation") or {}).get("proteins")) or []
    for i, spec in enumerate(proteins):
        if isinstance(spec, dict):
            _check_tokens(spec.get("hmms"), f"concatenation.proteins[{i}].hmms")

    cfg["hmm"] = {
        "enabled": enabled,
        "database": database.strip() if isinstance(database, str) else database,
        "default_evalue": float(ev),
        "use_ga_when_available": use_ga,
        "relative_length_cutoff": float(rel),
        "multidomain_overlap_tolerance": int(tol),
        "best_hit_only": best_only,
        "best_hit_also": [a.strip() for a in also],
        "min_protein_length": min_len,
        "threads": threads,
    }


def _warn_unknown_keys(cfg: dict, path: Path) -> list[str]:
    """Warn about keys in a family config file that this version does not
    read, at any depth (e.g. a retired ``length_outlier.hi_mult``).  Returns
    their dotted paths."""
    unknown = unknown_config_keys(cfg, DEFAULT_FAMILY_CONFIG)
    if unknown:
        log.warning(unknown_keys_message(path, unknown))
    return unknown


def preview_family_config(family: str, configs_dir: Path, global_cfg: dict) -> dict:
    """The merged config a run would use for *family*, without writing a
    file, warning or validating (for checks made before a run starts, such
    as which external tools it needs)."""
    path = configs_dir / f"{family}.yaml"
    if path.exists():
        with open(path) as f:
            file_cfg = yaml.safe_load(f) or {}
        return _merge_with_defaults(file_cfg, global_cfg, family)
    return _generate_default_family_config(family, global_cfg)


def load_family_config(family: str, configs_dir: Path, global_cfg: dict) -> tuple[dict, bool]:
    """Load per-family config, auto-generating it if missing.

    Returns (config_dict, was_auto_generated).
    """
    config_path = configs_dir / f"{family}.yaml"
    if config_path.exists():
        with open(config_path) as f:
            file_cfg = yaml.safe_load(f) or {}
        unknown = _warn_unknown_keys(file_cfg, config_path)
        cfg = _merge_with_defaults(file_cfg, global_cfg, family)
        _validate_sequence_mode(cfg, family, config_path)
        _warn_smart_default_conflicts(family, file_cfg, config_path, merged=cfg)
        # Repeated in the family log by the pipeline (logging is set up later).
        cfg["_unknown_keys"] = unknown
        _validate_numeric_fields(cfg, family)
        _validate_tree_blocks(cfg, family)
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
        _validate_additional_data_sources(cfg, family)
        _validate_hmm_block(cfg, family)
        return cfg, False
    else:
        cfg = _generate_default_family_config(family, global_cfg)
        _validate_sequence_mode(cfg, family, config_path)
        _validate_numeric_fields(cfg, family)
        _validate_tree_blocks(cfg, family)
        _validate_concatenation_block(cfg, family)
        _validate_manual_block(cfg, family)
        _validate_additional_data_sources(cfg, family)
        _validate_hmm_block(cfg, family)
        _write_family_config(cfg, config_path)
        log.warning(
            "No config found for %s — auto-generated default at %s. "
            "Edit it to tune parameters for this family.",
            family,
            config_path,
        )
        return cfg, True


def _validate_sequence_mode(cfg: dict, family: str, config_path: Path | None = None) -> None:
    """Reject HMM settings that cannot work with the rest of the merged config.

    HMM marker identification selects protein records, so it cannot run with
    ``sequence.type: nucleotide`` (the run would fetch proteins, then reject
    them all as too ambiguous for nucleotide QC), nor with a
    ``sequence.segment`` (the segment filter matches record titles, and
    protein titles name no segment, so every record would be dropped).  Both
    arise when a family file generated with the RdRp preset is edited only
    partly, e.g. just its ``sequence.type``.
    """
    if not (cfg.get("hmm") or {}).get("enabled"):
        return
    seq = cfg.get("sequence") or {}
    where = config_path.name if config_path else f"{family}.yaml"
    seg = SEGMENTED_FAMILIES.get(family)
    to_nucleotide = (
        "sequence.type: nucleotide, sequence.region: whole_genome, hmm.enabled: false"
        + (f" (and sequence.segment: {seg!r} for that segment)" if seg else "")
    )
    if seq.get("type") == "nucleotide":
        raise ValueError(
            f"{family}: hmm.enabled is true but sequence.type is nucleotide. HMM "
            "marker identification works on proteins only, so the run would fetch "
            f"{seq.get('region')!r} proteins and then reject them all in "
            f"nucleotide quality control. For whole genomes set {to_nucleotide} "
            f"in {where}; for the HMM protein marker set sequence.type: protein "
            f"(or delete {where} to regenerate it)."
        )
    if seq.get("segment"):
        raise ValueError(
            f"{family}: sequence.segment is {seq['segment']!r}, but {family} runs on "
            f"an HMM-identified protein marker ({seq.get('region')!r}), found on "
            "whichever segment carries it. Protein titles name no segment, so "
            "every record would be dropped. Remove sequence.segment from "
            f"{where}, or for that segment as nucleotide set sequence.type: "
            "nucleotide, sequence.region: whole_genome, hmm.enabled: false."
        )


def _warn_smart_default_conflicts(
    family: str, file_cfg: dict, config_path: Path, merged: dict | None = None,
) -> None:
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

    marker = HMM_MARKER_FAMILIES.get(family)
    if marker:
        if not hmm_marker_opt_out(file_cfg, marker):
            return  # the marker preset clears the segment; no segment check
        file_seq = file_cfg.get("sequence") or {}
        file_hmm = file_cfg.get("hmm") or {}
        # Report what will run (*merged*): a file's `segment: null` (as in a
        # file generated with the marker preset) wins over the usual segment.
        if merged is None:
            merged = _merge_with_defaults(file_cfg, make_minimal_global_cfg(), family)
        result = merged.get("sequence") or {}
        segment = result.get("segment")
        expected = SEGMENTED_FAMILIES.get(family)
        if segment:
            seg_text = f", segment {segment!r}"
        elif expected and result.get("type") == "nucleotide":
            seg_text = (
                f", no segment: all of its segments go into one alignment; add "
                f"sequence.segment: {expected!r} to {config_path.name} for that "
                "segment alone"
            )
        else:
            seg_text = ""
        chosen = [
            f"{key}={value!r}" for key, value in (
                ("sequence.type", file_seq.get("type")),
                ("sequence.region", file_seq.get("region")),
                ("hmm.enabled", file_hmm.get("enabled")),
            ) if value is not None
        ]
        log.warning(
            "%s: %s sets %s, so the recommended HMM-identified %s protein "
            "marker is not used: the file's own settings apply on top of the "
            "template defaults (sequence type %s, region %s%s). Whole-genome "
            "nucleotide alignment is saturated at family scale for most of the "
            "families on this marker. If this is a file from an earlier "
            "version, delete %s and re-run to regenerate.",
            family, config_path.name, ", ".join(chosen), marker,
            result.get("type"), result.get("region"), seg_text,
            config_path.name,
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
    """Apply CONCATENATION_FAMILIES, DNA_FAMILIES, HMM_MARKER_FAMILIES and
    SEGMENTED_FAMILIES as smart defaults in-place.  Precedence: concat >
    single-protein DNA overrides; HMM_MARKER_FAMILIES (disjoint from both)
    puts a family on its HMM-identified marker (RdRp, RT or Rep) with no
    segment, unless the family file opts out (`hmm_marker_opt_out`);
    otherwise segmentation is orthogonal and applied when relevant.

    *file_cfg* is the user's family config dict (before it is merged).  When
    provided, "Auto-configured" log messages are suppressed for any value the
    file already overrides — those messages only make sense when the smart
    default is actually taking effect.
    """
    _file_seq = (file_cfg or {}).get("sequence") or {}
    # The marker preset applies unless the family file chooses its own
    # analysis (hmm_marker_opt_out); it then clears the segment, so none is
    # auto-configured.
    marker = HMM_MARKER_FAMILIES.get(family)
    marker_on = marker is not None and not hmm_marker_opt_out(file_cfg, marker)

    segment = SEGMENTED_FAMILIES.get(family)
    if segment and not marker_on and not cfg["sequence"].get("segment"):
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

    if marker_on:
        preset = hmm_marker_preset(family)
        _deep_update(cfg, preset)
        if not _file_seq.get("type") and not _file_seq.get("region"):
            padding = preset["sequence"]["domain_padding"]
            log.info(
                "Auto-configured %s: %s protein marker (HMM-identified, %s)",
                family, marker,
                f"{marker} domain ± {padding} aa" if padding is not None else "whole protein",
            )


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
