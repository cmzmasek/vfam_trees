"""Marker identification for concatenated multi-protein mode.

A `MarkerIdentifier` maps a genome's protein records to its family's curated
marker set.  v1 (this module) does name + alias substring matching with
locus-tag and length tiebreakers.  v2 will swap in HMMER-based identification
behind the same Protocol — see CONCAT_DESIGN.md §5.1 and §8.

Marker specs are plain dicts (as delivered by the config layer); see
CONCAT_DESIGN.md §3.1 for the schema.
"""
from __future__ import annotations

import re
from typing import Protocol

from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)


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
