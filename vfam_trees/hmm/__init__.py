"""HMM-based marker-protein identification (HMMER hmmscan wrapper).

Vendored from the sister project ``repseq`` (``repseq/hmm/``) — pure stdlib +
HMMER subprocess, no new Python dependencies.  The token grammar and gating
logic (``parse_hmm_token``, ``cds_satisfies_token``, ``passes_cutoffs``) come
from repseq's ``clustering/marker.py`` selection layer.

Public API:
    - ``resolve_database_path``, ``ensure_pressed``, ``db_signature``,
      ``parse_ga_cutoffs``, ``profile_count`` — database management.
    - ``is_available``, ``scan`` — hmmscan invocation.
    - ``BUNDLED_HMMS_DIR`` — directory holding the bundled per-family sets.
    - ``passes_cutoffs``, ``coverage_of`` — per-hit GA/E-value + coverage gate.
    - ``parse_hmm_token``, ``cds_satisfies_token`` — multidomain token grammar.
    - ``HMMError``, ``HMMScanError``, ``HMMDatabaseError`` — exceptions.

See ``vfam_trees.markers.HMMIdentifier`` for how HMM hits feed marker
identification (the ``MarkerIdentifier`` Protocol swap point).
"""
from .database import (
    BUNDLED_DB_PATH,
    BUNDLED_HMMS_DIR,
    HMMPRESS_INDEX_SUFFIXES,
    db_signature,
    ensure_pressed,
    has_press_index,
    parse_ga_cutoffs,
    profile_count,
    resolve_database_path,
)
from .errors import HMMDatabaseError, HMMError, HMMScanError
from .hmmscan import is_available, scan
from .runner import (
    CACHE_SOURCE,
    cds_satisfies_token,
    coverage_of,
    get_ga_cutoffs,
    parse_hmm_token,
    passes_cutoffs,
    scan_proteins,
)

__all__ = [
    "BUNDLED_DB_PATH",
    "BUNDLED_HMMS_DIR",
    "CACHE_SOURCE",
    "HMMPRESS_INDEX_SUFFIXES",
    "HMMDatabaseError",
    "HMMError",
    "HMMScanError",
    "cds_satisfies_token",
    "coverage_of",
    "db_signature",
    "ensure_pressed",
    "get_ga_cutoffs",
    "has_press_index",
    "is_available",
    "parse_ga_cutoffs",
    "parse_hmm_token",
    "passes_cutoffs",
    "profile_count",
    "resolve_database_path",
    "scan",
    "scan_proteins",
]
