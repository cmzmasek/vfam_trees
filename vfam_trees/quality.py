"""Quality filtering of sequences."""
from __future__ import annotations

import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)

AMBIGUOUS_NUC = set("RYSWKMBDHVN")
AMBIGUOUS_AA = set("XBZJOU")

# Hard-coded absolute minimum lengths (applied even when user sets min_length explicitly)
MIN_LENGTH_NUC = 200
MIN_LENGTH_AA = 100

# Fractions of median tried in auto mode (most to least strict)
_AUTO_FRACTIONS = [0.5, 0.4, 0.3]


def filter_sequences(
    records: list[SeqRecord],
    seq_type: str,
    min_length: int | None,
    max_ambiguous: float,
    exclude_organisms: list[str] | None = None,
) -> tuple[list[SeqRecord], float | None, dict]:
    """Apply quality filters to a list of SeqRecords.

    Args:
        records: input sequences
        seq_type: 'nucleotide' or 'protein'
        min_length: minimum sequence length; None = auto (50% of median,
                    with fallback to 40%, 30% if no sequences pass)
        max_ambiguous: maximum fraction of ambiguous characters
        exclude_organisms: list of substrings matched case-insensitively
                           against the organism name in annotations

    Returns:
        Tuple of (filtered records, fraction_used, qc_stats).
        fraction_used is None when min_length was user-specified, or the
        median fraction (0.5, 0.4, or 0.3) that yielded >0 sequences.
        qc_stats is a dict with keys: n_excluded_organism, n_excluded_length,
        n_excluded_ambiguity, n_undefined.
    """
    empty_stats = {
        "n_excluded_organism": 0,
        "n_excluded_length": 0,
        "n_excluded_ambiguity": 0,
        "n_undefined": 0,
    }

    if not records:
        return [], None, empty_stats

    exclude_lower = [t.lower() for t in (exclude_organisms or [])]

    # Organism exclusion — applied before length filtering to keep median accurate
    passed_organism = []
    n_excluded_organism = 0
    for rec in records:
        organism = rec.annotations.get("organism", "").lower()
        if any(term in organism for term in exclude_lower):
            log.debug("Excluding organism: %s", rec.annotations.get("organism", ""))
            n_excluded_organism += 1
        else:
            passed_organism.append(rec)

    if n_excluded_organism:
        log.info("Excluded %d sequences matching organism exclusion list", n_excluded_organism)

    if not passed_organism:
        return [], None, {**empty_stats, "n_excluded_organism": n_excluded_organism}

    floor = MIN_LENGTH_NUC if seq_type == "nucleotide" else MIN_LENGTH_AA
    ambig_chars = AMBIGUOUS_NUC if seq_type == "nucleotide" else AMBIGUOUS_AA
    lengths = [len(r.seq) for r in passed_organism]
    median_len = statistics.median(lengths)

    if min_length is not None:
        # User-specified: apply once with floor enforcement
        effective = max(min_length, floor)
        if effective != min_length:
            log.info("min_length %d raised to hard floor %d", min_length, floor)
        else:
            log.info("Using min_length=%d", min_length)
        passed, n_short, n_ambig, n_undefined = _apply_length_filter(
            passed_organism, effective, ambig_chars, max_ambiguous
        )
        log.info(
            "Quality filter: %d passed, %d excluded organism, %d too short, "
            "%d too ambiguous, %d undefined sequence (from %d total)",
            len(passed), n_excluded_organism, n_short, n_ambig, n_undefined, len(records),
        )
        qc_stats = {
            "n_excluded_organism": n_excluded_organism,
            "n_excluded_length": n_short,
            "n_excluded_ambiguity": n_ambig,
            "n_undefined": n_undefined,
        }
        return passed, None, qc_stats

    # Auto mode: try fractions from most to least strict; stop when >0 sequences pass
    fraction_used = _AUTO_FRACTIONS[-1]
    passed: list[SeqRecord] = []
    n_short = n_ambig = n_undefined = 0
    for fraction in _AUTO_FRACTIONS:
        effective = max(int(median_len * fraction), floor)
        log.info(
            "Auto min_length set to %d (%.0f%% of median %d)",
            effective, fraction * 100, median_len,
        )
        passed, n_short, n_ambig, n_undefined = _apply_length_filter(
            passed_organism, effective, ambig_chars, max_ambiguous
        )
        log.info(
            "Quality filter: %d passed, %d excluded organism, %d too short, "
            "%d too ambiguous, %d undefined sequence (from %d total)",
            len(passed), n_excluded_organism, n_short, n_ambig, n_undefined, len(records),
        )
        fraction_used = fraction
        if passed:
            if fraction < 0.5:
                log.warning(
                    "Relaxed min_length threshold to %.0f%% of median (%d) to retain sequences",
                    fraction * 100, effective,
                )
            break

    qc_stats = {
        "n_excluded_organism": n_excluded_organism,
        "n_excluded_length": n_short,
        "n_excluded_ambiguity": n_ambig,
        "n_undefined": n_undefined,
    }
    return passed, fraction_used, qc_stats


def remove_length_outliers(
    records: list[SeqRecord],
    multiplier: float = 3.0,
) -> tuple[list[SeqRecord], int]:
    """Remove sequences longer than multiplier × median length.

    Args:
        records: input sequences
        multiplier: sequences longer than multiplier × median are removed

    Returns:
        (filtered records, number removed)
    """
    if len(records) < 2:
        return records, 0
    lengths = [len(r.seq) for r in records]
    median_len = statistics.median(lengths)
    cutoff = median_len * multiplier
    passed = [r for r in records if len(r.seq) <= cutoff]
    n_removed = len(records) - len(passed)
    if n_removed:
        log.info(
            "Removed %d length outlier(s) (>%.0fx median; cutoff %d bp/aa)",
            n_removed, multiplier, int(cutoff),
        )
    return passed, n_removed


def _apply_length_filter(
    records: list[SeqRecord],
    min_length: int,
    ambig_chars: set,
    max_ambiguous: float,
) -> tuple[list[SeqRecord], int, int, int]:
    """Filter records by length and ambiguity. Returns (passed, n_short, n_ambig, n_undefined)."""
    passed = []
    n_short = 0
    n_ambig = 0
    n_undefined = 0
    for rec in records:
        try:
            seq_str = str(rec.seq).upper()
        except Exception:
            n_undefined += 1
            continue
        if len(seq_str) < min_length:
            n_short += 1
            continue
        ambig_frac = sum(1 for c in seq_str if c in ambig_chars) / max(len(seq_str), 1)
        if ambig_frac > max_ambiguous:
            n_ambig += 1
            continue
        passed.append(rec)
    return passed, n_short, n_ambig, n_undefined


def deduplicate(records: list[SeqRecord]) -> list[SeqRecord]:
    """Remove exact duplicate accessions, keeping first occurrence."""
    seen = set()
    unique = []
    for rec in records:
        acc = rec.id
        if acc not in seen:
            seen.add(acc)
            unique.append(rec)
    removed = len(records) - len(unique)
    if removed:
        log.info("Removed %d duplicate accessions", removed)
    return unique


def write_fasta(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        SeqIO.write(records, f, "fasta")
    log.debug("Wrote %d sequences to %s", len(records), path)
