"""Marker identification for single-marker and concatenated multi-protein mode.

A `MarkerIdentifier` maps a genome's protein records to its family's curated
marker set.  `NameMatchIdentifier` does name + alias substring matching with
locus-tag and length tiebreakers.  `HMMIdentifier` (the CONCAT_DESIGN §5.1/§8
"v2 swap") identifies markers by HMMER profile homology behind the same
Protocol, falling back to name matching for marker specs that declare no
`hmms`.  `make_identifier(cfg)` picks between them from `hmm.enabled`.

Marker specs are plain dicts (as delivered by the config layer); see
CONCAT_DESIGN.md §3.1 for the schema.  An HMM-gated spec additionally carries
`hmms: [token, ...]` — a single `"Name"` or multidomain `"A--B--C"` token
(HMMs in N-to-C order); tokens within one spec are alternatives (OR'd).
"""
from __future__ import annotations

import re
from typing import Protocol

from Bio.SeqRecord import SeqRecord

from .hmm import (
    cds_satisfies_token,
    get_ga_cutoffs,
    parse_hmm_token,
    passes_cutoffs,
    scan_proteins,
)
from .logger import get_logger

log = get_logger(__name__)

# Key under which passing-annotated HMM hits are cached on each SeqRecord's
# ``.annotations`` dict, so a genome's proteins are scanned at most once.
_HMM_HITS_KEY = "hmm_hits"


class MarkerIdentifier(Protocol):
    """Map proteins from one genome to the family's curated marker set.

    Implementations choose at most one protein per marker, applying any
    spec-level disambiguators (locus_tag_hint, length_range).  Markers
    with no candidate are simply omitted from the result — the caller
    decides whether the genome's coverage is acceptable.
    """

    def identify(
        self,
        proteins: list[SeqRecord],
        marker_set: list[dict],
        species_lineage: list[dict] | None = None,
    ) -> dict[str, SeqRecord]:
        ...


class NameMatchIdentifier:
    """v1 — name + alias substring matching, with tiebreaker priority:

    1. ``locus_tag_hint`` regex match (when the spec defines one).
    2. Sequence length closest to ``length_range`` midpoint (when defined).
    3. Longest sequence.
    4. Lowest accession (deterministic final tiebreaker).
    """

    def identify(
        self,
        proteins: list[SeqRecord],
        marker_set: list[dict],
        species_lineage: list[dict] | None = None,
    ) -> dict[str, SeqRecord]:
        subfamily = _subfamily_from_lineage(species_lineage)
        result: dict[str, SeqRecord] = {}
        for marker in marker_set:
            chosen = self._identify_one(proteins, marker, subfamily)
            if chosen is not None:
                result[marker["name"]] = chosen
        return result

    def _identify_one(
        self,
        proteins: list[SeqRecord],
        marker: dict,
        subfamily: str | None,
    ) -> SeqRecord | None:
        names = self._effective_aliases(marker, subfamily)
        candidates = [p for p in proteins if _record_matches_any(p, names)]
        if not candidates:
            return None
        return _apply_tiebreakers(candidates, marker)

    @staticmethod
    def _effective_aliases(marker: dict, subfamily: str | None) -> list[str]:
        names = [marker["name"]] + list(marker.get("aliases", []))
        if subfamily:
            extra = marker.get(f"aliases_{subfamily}")
            if extra:
                names.extend(extra)
        return names


def annotate_hmm_hits(
    records: list[SeqRecord], cfg: dict, *, force: bool = False
) -> None:
    """Batch-scan protein ``records`` against the configured HMM DB once and
    cache passing-annotated hits on each record's ``annotations["hmm_hits"]``.

    Scanning many proteins in one ``hmmscan`` call (then letting
    `HMMIdentifier` read the pre-computed hits per genome) is far cheaper than
    one scan per genome: hmmscan's start-up cost is paid once and identical AA
    sequences are de-duplicated inside `hmm.scan_proteins`.  Records already
    carrying the key are skipped unless ``force``.

    Mutates each record in place; relies on HMMER being on PATH and the
    database resolving (raises `HMMScanError` / `HMMDatabaseError` otherwise —
    HMM identification is authoritative when enabled, so failures surface).
    """
    to_scan = [
        r for r in records
        if force or _HMM_HITS_KEY not in (r.annotations or {})
    ]
    if not to_scan:
        return
    hcfg = cfg.get("hmm", {}) or {}
    ga_cutoffs = get_ga_cutoffs(cfg)
    default_evalue = float(hcfg.get("default_evalue", 1.0e-5))
    rel_len = float(hcfg.get("relative_length_cutoff", 0.5))
    use_ga = bool(hcfg.get("use_ga_when_available", True))

    # Index by a synthetic id so duplicate accessions can't collide as keys.
    idx_to_rec = {f"R{i:07d}": r for i, r in enumerate(to_scan)}
    queries = {k: str(r.seq) for k, r in idx_to_rec.items() if r.seq}
    raw = scan_proteins(queries, cfg) if queries else {}
    for k, r in idx_to_rec.items():
        hits = [dict(h) for h in raw.get(k, [])]
        for h in hits:
            h["passing"] = passes_cutoffs(
                h, ga_cutoffs, default_evalue, rel_len, use_ga,
            )
        r.annotations[_HMM_HITS_KEY] = hits


class HMMIdentifier:
    """HMMER-based marker identification behind the `MarkerIdentifier` Protocol
    (CONCAT_DESIGN §5.1/§8 — the "v2 swap").

    For each marker spec that declares ``hmms`` (a list of token strings), the
    HMM tier is **authoritative**: a protein satisfies the spec when its
    passing hits satisfy ANY token (tokens OR'd; multidomain ``A--B`` tokens
    require N-to-C order), and the **longest** satisfying protein is chosen —
    no name fallback for that spec.  Specs without ``hmms`` fall back to
    name+alias matching via `NameMatchIdentifier`, so mixed marker sets work.

    Hits are computed once per genome (or in bulk beforehand via
    `annotate_hmm_hits` / `prescan`) and cached on each record's annotations,
    so repeated `identify` calls over the same records do not rescan.
    """

    #: Lets the fetcher detect HMM mode (all-proteins fetch) without importing
    #: this class — ``getattr(identifier, "is_hmm", False)``.
    is_hmm = True

    def __init__(self, cfg: dict):
        self._cfg = cfg or {}
        hcfg = self._cfg.get("hmm", {}) or {}
        self._overlap_tol = int(hcfg.get("multidomain_overlap_tolerance", 30))
        self._name_fallback = NameMatchIdentifier()

    def prescan(self, proteins: list[SeqRecord]) -> None:
        """Batch-scan ``proteins`` once so later per-genome `identify` calls read
        cached hits instead of rescanning (one hmmscan per species, not per
        genome).  Safe to call repeatedly — already-annotated records are
        skipped."""
        annotate_hmm_hits(proteins, self._cfg)

    def identify(
        self,
        proteins: list[SeqRecord],
        marker_set: list[dict],
        species_lineage: list[dict] | None = None,
    ) -> dict[str, SeqRecord]:
        if not proteins:
            return {}
        annotate_hmm_hits(proteins, self._cfg)
        subfamily = _subfamily_from_lineage(species_lineage)
        result: dict[str, SeqRecord] = {}
        for marker in marker_set:
            tokens = list(marker.get("hmms") or [])
            if tokens:
                chosen = self._identify_one_hmm(proteins, tokens)
            else:
                # No HMM tokens for this marker → name+alias fallback.
                chosen = self._name_fallback._identify_one(
                    proteins, marker, subfamily,
                )
            if chosen is not None:
                result[marker["name"]] = chosen
        return result

    def _identify_one_hmm(
        self, proteins: list[SeqRecord], tokens: list[str]
    ) -> SeqRecord | None:
        satisfying = [
            p for p in proteins
            if _record_satisfies_any_token(
                p.annotations.get(_HMM_HITS_KEY) or [], tokens, self._overlap_tol,
            )
        ]
        if not satisfying:
            return None
        # Longest satisfying protein wins; lowest accession breaks ties.
        return sorted(satisfying, key=lambda p: (-len(p.seq), p.id or ""))[0]


def _record_satisfies_any_token(
    hits: list[dict], tokens: list[str], overlap_tol: int
) -> bool:
    for token in tokens:
        try:
            parsed = parse_hmm_token(token)
        except ValueError:
            continue
        if cds_satisfies_token(hits, parsed, overlap_tolerance=overlap_tol) is not None:
            return True
    return False


def make_identifier(cfg: dict) -> MarkerIdentifier:
    """Select the marker identifier from config.

    ``hmm.enabled: true`` → `HMMIdentifier` (HMM-authoritative for specs with
    ``hmms``, name fallback otherwise).  Otherwise the default
    `NameMatchIdentifier`.
    """
    hcfg = (cfg or {}).get("hmm", {}) or {}
    if hcfg.get("enabled"):
        return HMMIdentifier(cfg)
    return NameMatchIdentifier()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _subfamily_from_lineage(lineage: list[dict] | None) -> str | None:
    """Return the subfamily-rank entry name from an NCBI ranked lineage."""
    if not lineage:
        return None
    for entry in lineage:
        if entry.get("rank") == "subfamily":
            return entry.get("name")
    return None


def _record_matches_any(record: SeqRecord, names: list[str]) -> bool:
    haystack = _record_name_text(record).lower()
    return any(n.lower() in haystack for n in names if n)


def _record_name_text(record: SeqRecord) -> str:
    """Concatenate the protein-naming fields we match against.

    Uses the DEFINITION line plus any ``product``/``name`` qualifiers from
    Protein/CDS features.  Order is deliberate: description first so the
    most authoritative name dominates.
    """
    parts: list[str] = []
    desc = getattr(record, "description", None)
    if desc:
        parts.append(desc)
    for feat in getattr(record, "features", None) or []:
        if getattr(feat, "type", None) not in ("Protein", "CDS"):
            continue
        qualifiers = getattr(feat, "qualifiers", None) or {}
        for key in ("product", "name"):
            vals = qualifiers.get(key, [])
            parts.extend(v for v in vals if isinstance(v, str))
    return " ".join(parts)


def _apply_tiebreakers(candidates: list[SeqRecord], marker: dict) -> SeqRecord:
    if len(candidates) == 1:
        return candidates[0]

    # 1. locus_tag_hint regex
    hint = marker.get("locus_tag_hint")
    if hint:
        rx = re.compile(hint, re.IGNORECASE)
        hits = [c for c in candidates if rx.search(_record_locus_tag(c))]
        if hits:
            candidates = hits
            if len(candidates) == 1:
                return candidates[0]

    # 2. length_range midpoint proximity
    lr = marker.get("length_range")
    if lr:
        midpoint = (lr[0] + lr[1]) / 2.0
        candidates = sorted(
            candidates,
            key=lambda c: (abs(len(c.seq) - midpoint), c.id or ""),
        )
        return candidates[0]

    # 3. longest sequence; 4. lowest accession breaks ties
    return sorted(candidates, key=lambda c: (-len(c.seq), c.id or ""))[0]


def _record_locus_tag(record: SeqRecord) -> str:
    """Best-effort locus_tag extraction from annotations or feature qualifiers."""
    annotations = getattr(record, "annotations", None) or {}
    tag = annotations.get("locus_tag", "")
    if tag:
        return tag
    for feat in getattr(record, "features", None) or []:
        qualifiers = getattr(feat, "qualifiers", None) or {}
        tags = qualifiers.get("locus_tag", [])
        if tags:
            return tags[0]
    return ""
