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

import hashlib
import re
from pathlib import Path
from typing import Protocol

from Bio.SeqRecord import SeqRecord

from .hmm import (
    cds_satisfies_token,
    get_ga_cutoffs,
    parse_hmm_token,
    passes_cutoffs,
    resolve_database_path,
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
    records: list[SeqRecord], cfg: dict, *, force: bool = False,
    ga_cutoffs: dict | None = None,
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
    *ga_cutoffs* (from `HMMIdentifier`, parsed once per run) saves re-reading
    the database; without it they are parsed here, and only when
    ``hmm.use_ga_when_available`` is on.
    """
    to_scan = [
        r for r in records
        if force or _HMM_HITS_KEY not in (r.annotations or {})
    ]
    if not to_scan:
        return
    hcfg = cfg.get("hmm", {}) or {}
    default_evalue = float(hcfg.get("default_evalue", 1.0e-5))
    rel_len = float(hcfg.get("relative_length_cutoff", 0.5))
    use_ga = bool(hcfg.get("use_ga_when_available", True))
    if ga_cutoffs is None:
        ga_cutoffs = get_ga_cutoffs(cfg) if use_ga else {}

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

    With ``hmm.best_hit_only`` a protein must also fit one of the spec's
    profiles better than any other profile in the database (`best_passing_target`),
    so a protein that only cross-hits a token — e.g. an orbivirus VP1 weakly
    hitting a palm-domain RdRp profile — is not taken as the marker.  Profiles
    in ``hmm.best_hit_also`` (other domains of the same protein) count as the
    spec's own for that test.
    """

    #: Lets the fetcher detect HMM mode (all-proteins fetch) without importing
    #: this class — ``getattr(identifier, "is_hmm", False)``.
    is_hmm = True

    def __init__(self, cfg: dict):
        self._cfg = cfg or {}
        hcfg = self._cfg.get("hmm", {}) or {}
        self._overlap_tol = int(hcfg.get("multidomain_overlap_tolerance", 30))
        self._best_hit_only = bool(hcfg.get("best_hit_only", False))
        self._best_hit_also = frozenset(hcfg.get("best_hit_also") or ())
        self._min_protein_length = hcfg.get("min_protein_length")
        self._use_ga = bool(hcfg.get("use_ga_when_available", True))
        self._ga_cutoffs: dict | None = None      # parsed once, when needed
        self._db_fingerprint: str | None = None
        self._name_fallback = NameMatchIdentifier()

    def _annotate(self, records: list[SeqRecord]) -> None:
        """`annotate_hmm_hits` with GA cutoffs parsed once per identifier,
        and only when there is something to scan."""
        if all(_HMM_HITS_KEY in (r.annotations or {}) for r in records):
            return
        if self._ga_cutoffs is None:
            self._ga_cutoffs = get_ga_cutoffs(self._cfg) if self._use_ga else {}
        annotate_hmm_hits(records, self._cfg, ga_cutoffs=self._ga_cutoffs)

    def database_fingerprint(self) -> str:
        """First 16 hex digits of the sha256 of the profile database's
        content (``hmm.database`` resolved; a directory set is its sorted
        concatenation).  Independent of install path and file times, so it
        changes exactly when a profile is added, removed or edited."""
        if self._db_fingerprint is None:
            hcfg = self._cfg.get("hmm", {}) or {}
            cache_dir = self._cfg.get("cache_dir")   # as hmm.runner resolves it
            db_path = resolve_database_path(
                hcfg.get("database"),
                cache_dir=Path(cache_dir).expanduser() if cache_dir else None,
            )
            self._db_fingerprint = hashlib.sha256(db_path.read_bytes()).hexdigest()[:16]
        return self._db_fingerprint

    @property
    def best_hit_only(self) -> bool:
        """Whether a marker protein must fit one of its spec's profiles best."""
        return self._best_hit_only

    @property
    def min_protein_length(self) -> int | None:
        """Shortest protein kept as a candidate (``hmm.min_protein_length``)."""
        return self._min_protein_length

    def best_hit_profiles(self, tokens: list[str]) -> set[str]:
        """Profiles a protein's best hit may be on for *tokens*' marker under
        ``best_hit_only``: the tokens' own profiles plus ``hmm.best_hit_also``."""
        return token_profile_names(tokens) | self._best_hit_also

    def crop(self, records: list[SeqRecord], tokens: list[str], padding: int) -> list[SeqRecord]:
        """Scan *records* (if not yet scanned) and cut each to its marker
        domain ± *padding* (`crop_to_marker_domain`)."""
        self._annotate(records)
        return [crop_to_marker_domain(r, tokens, padding) for r in records]

    def prescan(self, proteins: list[SeqRecord]) -> None:
        """Batch-scan ``proteins`` once so later per-genome `identify` calls read
        cached hits instead of rescanning (one hmmscan per species, not per
        genome).  Safe to call repeatedly — already-annotated records are
        skipped."""
        self._annotate(proteins)

    def identify(
        self,
        proteins: list[SeqRecord],
        marker_set: list[dict],
        species_lineage: list[dict] | None = None,
    ) -> dict[str, SeqRecord]:
        if not proteins:
            return {}
        self._annotate(proteins)
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
        if self._best_hit_only:
            names = self.best_hit_profiles(tokens)
            satisfying = [
                p for p in satisfying
                if best_passing_target(p.annotations.get(_HMM_HITS_KEY) or []) in names
            ]
        if not satisfying:
            return None
        # Longest satisfying protein wins; lowest accession breaks ties.
        return sorted(satisfying, key=lambda p: (-len(p.seq), p.id or ""))[0]


def token_profile_names(tokens: list[str]) -> set[str]:
    """Every profile name used in *tokens* (multidomain tokens split);
    malformed tokens are skipped (config validation reports them)."""
    names: set[str] = set()
    for token in tokens or []:
        try:
            names.update(parse_hmm_token(token))
        except ValueError:
            continue
    return names


def marker_hits(record: SeqRecord, tokens: list[str]) -> list[dict]:
    """*record*'s passing hits (cached by `annotate_hmm_hits`) on the profiles
    named in *tokens*: the hits `crop_to_marker_domain` cuts by.  Empty when
    the record was not scanned or has no such hit."""
    names = token_profile_names(tokens)
    return [
        h for h in (record.annotations.get(_HMM_HITS_KEY) or [])
        if h.get("passing") and h.get("target") in names
    ]


def best_passing_target(hits: list[dict]) -> str | None:
    """The profile that fits a protein best among its passing hits: highest
    full-sequence bit score, domain score breaking ties.  None without a
    passing hit."""
    passing = [h for h in hits or [] if h.get("passing")]
    if not passing:
        return None
    best = max(
        passing,
        key=lambda h: (h.get("full_score") or 0.0, h.get("dom_score") or 0.0),
    )
    return best.get("target")


def foreign_best_hit(record: SeqRecord, own_profiles: set[str]) -> str | None:
    """The best-fitting profile of *record* when it is not in *own_profiles*
    (a protein of another kind than the marker; see
    `HMMIdentifier.best_hit_profiles`), else None.  Reads hits cached by
    `annotate_hmm_hits`; None for a record never scanned or without a passing
    hit."""
    best = best_passing_target(record.annotations.get(_HMM_HITS_KEY) or [])
    if best is None or best in own_profiles:
        return None
    return best


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


# Passing marker hits further apart than this (aa) on one protein are not the
# same domain region.  On the RefSeq proteins of the 107 RdRp families the
# largest gap between neighbouring marker-domain hits was 383 aa (a
# phenuivirus L between its PA-C-like and RdRp domains; 99 % were <= 339).
MAX_DOMAIN_GAP = 500


def _marker_hit_chain(hits: list[dict]) -> list[dict]:
    """The run of *hits* around the best-scoring one in which neighbouring
    hits are at most `MAX_DOMAIN_GAP` residues apart."""
    best = max(hits, key=lambda h: (h.get("dom_score") or 0.0, h.get("full_score") or 0.0))
    ordered = sorted(hits, key=lambda h: (h["ali_from"], h["ali_to"]))
    chains: list[list[dict]] = []
    end = None
    for h in ordered:
        if end is None or h["ali_from"] - end > MAX_DOMAIN_GAP:
            chains.append([h])
            end = h["ali_to"]
        else:
            chains[-1].append(h)
            end = max(end, h["ali_to"])
    return next(c for c in chains if any(h is best for h in c))


def crop_to_marker_domain(
    record: SeqRecord, tokens: list[str], padding: int
) -> SeqRecord:
    """Return *record* cut to the span its marker hits cover, ± *padding* aa.

    The span runs from the first to the last residue aligned to the passing
    hits of profiles named in *tokens* (cached on ``annotations["hmm_hits"]``
    by `annotate_hmm_hits`) that form one region with the best-scoring hit:
    neighbouring hits at most `MAX_DOMAIN_GAP` apart.  A distant stray hit,
    or a second copy of the domain in a chimeric record, is left out (and
    noted in the COMMENT).  Used to keep only the conserved domain of a
    large multi-domain protein, e.g. the RdRp of a viral polyprotein: whole
    polyproteins would otherwise dominate the alignment with variable
    structural regions.

    Metadata survives: annotations are kept, the ``source`` and ``CDS``
    features keep their qualifiers (``/host``, ``/coded_by`` …) with their
    location reset to the new length, other features are dropped, and a
    COMMENT records which residues were kept.  Returns *record* itself when
    there is no matching hit or the span already covers the whole protein.
    """
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    hits = marker_hits(record, tokens)
    length = len(record.seq)
    if not hits or not length:
        return record
    chain = _marker_hit_chain(hits)
    n_distant = len(hits) - len(chain)
    if n_distant:
        log.info(
            "%s: %d marker hit(s) more than %d aa from the best one left out of "
            "the domain crop", record.id, n_distant, MAX_DOMAIN_GAP,
        )
    hits = chain
    start = max(0, min(h["ali_from"] for h in hits) - 1 - max(0, padding))
    end = min(length, max(h["ali_to"] for h in hits) + max(0, padding))
    if start == 0 and end == length:
        return record

    cropped = SeqRecord(
        record.seq[start:end],
        id=record.id,
        name=record.name,
        description=record.description,
    )
    cropped.annotations = {
        k: v for k, v in record.annotations.items() if k != _HMM_HITS_KEY
    }
    cropped.annotations.setdefault("molecule_type", "protein")
    note = (
        f"vfam_trees: residues {start + 1}-{end} of {length} "
        f"(marker domain {', '.join(sorted({h['target'] for h in hits}))} ± {padding} aa"
        + (f"; {n_distant} distant hit(s) left out" if n_distant else "")
        + ")"
    )
    previous = cropped.annotations.get("comment")
    cropped.annotations["comment"] = f"{previous}\n{note}" if previous else note
    cropped.features = [
        SeqFeature(FeatureLocation(0, end - start), type=f.type, qualifiers=dict(f.qualifiers))
        for f in record.features
        if f.type in ("source", "CDS")
    ]
    return cropped


def make_identifier(cfg: dict, threads: int | None = None) -> MarkerIdentifier:
    """Select the marker identifier from config.

    ``hmm.enabled: true`` → `HMMIdentifier` (HMM-authoritative for specs with
    ``hmms``, name fallback otherwise).  Otherwise the default
    `NameMatchIdentifier`.  *threads* (the run's thread count) is what
    hmmscan uses when ``hmm.threads`` is null.
    """
    hcfg = (cfg or {}).get("hmm", {}) or {}
    if hcfg.get("enabled"):
        if threads and not hcfg.get("threads"):
            cfg = {**cfg, "hmm": {**hcfg, "threads": int(threads)}}
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
