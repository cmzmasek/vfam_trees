"""Quality filtering of sequences."""
from __future__ import annotations

import math
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
                           against the ORGANISM, SOURCE, and DEFINITION fields

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

    # Organism exclusion — applied before length filtering to keep median accurate.
    # Terms are matched against ORGANISM, SOURCE, and DEFINITION (case-insensitive
    # substring), joined with a newline to prevent cross-field false matches.
    passed_organism = []
    n_excluded_organism = 0
    for rec in records:
        search_text = "\n".join([
            rec.annotations.get("organism", ""),
            rec.annotations.get("source", ""),
            rec.description,
        ]).lower()
        if any(term in search_text for term in exclude_lower):
            log.info("Excluding record: %s", rec.annotations.get("organism", ""))
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
        passed, n_short, n_ambig, n_undefined, pre_length_lengths = _apply_length_filter(
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
            "pre_length_lengths": pre_length_lengths,
        }
        return passed, None, qc_stats

    # Auto mode: try fractions from most to least strict; stop when >0 sequences pass
    fraction_used = _AUTO_FRACTIONS[-1]
    passed: list[SeqRecord] = []
    n_short = n_ambig = n_undefined = 0
    pre_length_lengths: list[int] = []
    for fraction in _AUTO_FRACTIONS:
        effective = max(int(median_len * fraction), floor)
        log.info(
            "Auto min_length set to %d (%.0f%% of median %d)",
            effective, fraction * 100, median_len,
        )
        passed, n_short, n_ambig, n_undefined, pre_length_lengths = _apply_length_filter(
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
        "pre_length_lengths": pre_length_lengths,
    }
    return passed, fraction_used, qc_stats


def log_mad_cutoffs(
    lengths: list[int],
    k: float,
) -> tuple[float | None, float | None, float, float]:
    """Compute length-outlier cutoffs from MAD on log-lengths.

    Why log-space: protein/CDS length is right-skewed and bounded at zero, so
    a symmetric ±k·σ window in log-space becomes a multiplicative window
    around the median — the natural shape for length data.

    Why MAD: the standard deviation of log-lengths would itself be inflated
    by the very outliers we want to detect.  MAD (median absolute deviation),
    scaled by 1.4826, is a robust estimator of σ for normal data.

    Returns ``(lo_cutoff, hi_cutoff, sigma_log, median_len)``.  Returns
    ``(None, None, 0.0, median_len)`` when ``k <= 0`` (filter disabled) or
    when MAD is exactly zero (data too uniform to characterize natural
    spread — filter degenerates to a no-op rather than fabricating one).
    """
    valid = [L for L in lengths if L > 0]
    if not valid:
        return None, None, 0.0, 0.0
    median_len = float(statistics.median(valid))
    if k <= 0:
        return None, None, 0.0, median_len
    log_lens = [math.log(L) for L in valid]
    med = statistics.median(log_lens)
    mad = statistics.median(abs(x - med) for x in log_lens)
    if mad <= 0:
        return None, None, 0.0, median_len
    sigma = 1.4826 * mad
    return (
        math.exp(med - k * sigma),
        math.exp(med + k * sigma),
        sigma,
        median_len,
    )


def remove_length_outliers(
    records: list[SeqRecord],
    k: float = 5.0,
    protected_ids: set[str] | None = None,
) -> tuple[list[SeqRecord], int, int]:
    """Remove sequences whose length is a MAD-on-log outlier.

    Args:
        records: input sequences
        k: number of MAD-units (in log-space, scaled by 1.4826) on either
           side of the log-median that defines the keep window.  A larger
           k is more lenient.  k <= 0 disables the filter.
        protected_ids: record IDs that must never be dropped even if flagged —
                       a warning is logged instead.

    Returns:
        (filtered records, n_long_removed, n_short_removed)
    """
    if len(records) < 2:
        return records, 0, 0
    protected = protected_ids or set()
    lengths = [len(r.seq) for r in records]
    lo_cutoff, hi_cutoff, sigma, median_len = log_mad_cutoffs(lengths, k)
    if lo_cutoff is None or hi_cutoff is None:
        return records, 0, 0

    passed: list[SeqRecord] = []
    n_long = 0
    n_short = 0
    for r in records:
        L = len(r.seq)
        too_long = L > hi_cutoff
        too_short = L < lo_cutoff
        if too_long or too_short:
            if r.id in protected:
                log.warning(
                    "RefSeq '%s' looks like a length outlier "
                    "(length=%d, median=%d, lo_cutoff=%.0f, hi_cutoff=%.0f) — "
                    "KEEPING (RefSeq protected)",
                    r.id, L, int(median_len), lo_cutoff, hi_cutoff,
                )
                passed.append(r)
            elif too_long:
                n_long += 1
            else:
                n_short += 1
        else:
            passed.append(r)

    if n_long or n_short:
        log.info(
            "Removed %d length outlier(s) for %d seqs "
            "(median=%d, sigma_log=%.3f, k=%.1f, "
            "keep window=[%.0f, %.0f]): %d long, %d short",
            n_long + n_short, len(records), int(median_len),
            sigma, k, lo_cutoff, hi_cutoff, n_long, n_short,
        )
    return passed, n_long, n_short


def _apply_length_filter(
    records: list[SeqRecord],
    min_length: int,
    ambig_chars: set,
    max_ambiguous: float,
) -> tuple[list[SeqRecord], int, int, int, list[int]]:
    """Filter records by ambiguity then length.

    Returns (passed, n_short, n_ambig, n_undefined, pre_length_lengths).
    pre_length_lengths: lengths of sequences that passed the ambiguity check,
    before the minimum-length check is applied.
    """
    passed = []
    n_short = 0
    n_ambig = 0
    n_undefined = 0
    pre_length_lengths: list[int] = []
    for rec in records:
        try:
            seq_str = str(rec.seq).upper()
        except Exception:
            n_undefined += 1
            continue
        ambig_frac = sum(1 for c in seq_str if c in ambig_chars) / max(len(seq_str), 1)
        if ambig_frac > max_ambiguous:
            n_ambig += 1
            continue
        pre_length_lengths.append(len(seq_str))
        if len(seq_str) < min_length:
            n_short += 1
            continue
        passed.append(rec)
    return passed, n_short, n_ambig, n_undefined, pre_length_lengths


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
