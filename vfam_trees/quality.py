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
    max_length: int | None = None,
    exclude_organisms: list[str] | None = None,
    manual_include_ids: set[str] | None = None,
    manual_exclude_ids: set[str] | None = None,
) -> tuple[list[SeqRecord], float | None, dict]:
    """Apply quality filters to a list of SeqRecords.

    Args:
        records: input sequences
        seq_type: 'nucleotide' or 'protein'
        min_length: minimum sequence length; None = auto (50% of median,
                    with fallback to 40%, 30% if no sequences pass)
        max_ambiguous: maximum fraction of ambiguous characters
        max_length: maximum sequence length; None = no upper limit
        exclude_organisms: list of substrings matched case-insensitively
                           against the ORGANISM, SOURCE, and DEFINITION fields
        manual_include_ids: accessions (exact match incl. version) that bypass
                            all QC and are appended to the passing set unchanged
        manual_exclude_ids: accessions dropped before any QC step

    Returns:
        Tuple of (filtered records, fraction_used, qc_stats).
        fraction_used is None when min_length was user-specified, or the
        median fraction (0.5, 0.4, or 0.3) that yielded >0 sequences.
        qc_stats is a dict with keys: n_excluded_organism, n_excluded_length,
        n_excluded_ambiguity, n_undefined, n_manual_include_bypassed,
        n_manual_exclude_dropped, manual_include_seen, manual_exclude_seen.
    """
    include_set = set(manual_include_ids or ())
    exclude_set = set(manual_exclude_ids or ())

    empty_stats = {
        "n_excluded_organism": 0,
        "n_excluded_length": 0,
        "n_excluded_ambiguity": 0,
        "n_undefined": 0,
        "n_manual_include_bypassed": 0,
        "n_manual_exclude_dropped": 0,
        "manual_include_seen": set(),
        "manual_exclude_seen": set(),
    }

    if not records:
        return [], None, empty_stats

    # Manual exclude/include passes — applied before any QC step.
    # Excludes win over includes (config validation prevents overlap, but be
    # defensive in case a caller passes overlapping sets directly).
    manual_include_seen: set[str] = set()
    manual_exclude_seen: set[str] = set()
    bypassed: list[SeqRecord] = []
    remaining: list[SeqRecord] = []
    n_manual_exclude_dropped = 0
    for rec in records:
        if rec.id in exclude_set:
            manual_exclude_seen.add(rec.id)
            n_manual_exclude_dropped += 1
            log.info("manual.exclude dropping %s", rec.id)
            continue
        if rec.id in include_set:
            manual_include_seen.add(rec.id)
            bypassed.append(rec)
            continue
        remaining.append(rec)

    n_manual_include_bypassed = len(bypassed)
    if n_manual_include_bypassed:
        log.info(
            "manual.include bypassed QC for %d record(s): %s",
            n_manual_include_bypassed,
            ", ".join(sorted(manual_include_seen)),
        )

    exclude_lower = [t.lower() for t in (exclude_organisms or [])]

    manual_stats = {
        "n_manual_include_bypassed": n_manual_include_bypassed,
        "n_manual_exclude_dropped": n_manual_exclude_dropped,
        "manual_include_seen": manual_include_seen,
        "manual_exclude_seen": manual_exclude_seen,
    }

    # Organism exclusion — applied before length filtering to keep median accurate.
    # Terms are matched against ORGANISM, SOURCE, and DEFINITION (case-insensitive
    # substring), joined with a newline to prevent cross-field false matches.
    passed_organism = []
    n_excluded_organism = 0
    for rec in remaining:
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
        # Bypassed records still go through even when nothing else survives QC.
        return list(bypassed), None, {
            **empty_stats,
            "n_excluded_organism": n_excluded_organism,
            **manual_stats,
        }

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
        passed, n_excl_len, n_ambig, n_undefined, pre_length_lengths = _apply_length_filter(
            passed_organism, effective, ambig_chars, max_ambiguous, max_length
        )
        log.info(
            "Quality filter: %d passed, %d excluded organism, %d excluded by length, "
            "%d too ambiguous, %d undefined sequence (from %d total)",
            len(passed), n_excluded_organism, n_excl_len, n_ambig, n_undefined, len(records),
        )
        qc_stats = {
            "n_excluded_organism": n_excluded_organism,
            "n_excluded_length": n_excl_len,
            "n_excluded_ambiguity": n_ambig,
            "n_undefined": n_undefined,
            "pre_length_lengths": pre_length_lengths,
            **manual_stats,
        }
        return list(bypassed) + passed, None, qc_stats

    # Auto mode: try fractions from most to least strict; stop when >0 sequences pass
    fraction_used = _AUTO_FRACTIONS[-1]
    passed: list[SeqRecord] = []
    n_excl_len = n_ambig = n_undefined = 0
    pre_length_lengths: list[int] = []
    for fraction in _AUTO_FRACTIONS:
        effective = max(int(median_len * fraction), floor)
        log.info(
            "Auto min_length set to %d (%.0f%% of median %d)",
            effective, fraction * 100, median_len,
        )
        passed, n_excl_len, n_ambig, n_undefined, pre_length_lengths = _apply_length_filter(
            passed_organism, effective, ambig_chars, max_ambiguous, max_length
        )
        log.info(
            "Quality filter: %d passed, %d excluded organism, %d excluded by length, "
            "%d too ambiguous, %d undefined sequence (from %d total)",
            len(passed), n_excluded_organism, n_excl_len, n_ambig, n_undefined, len(records),
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
        "n_excluded_length": n_excl_len,
        "n_excluded_ambiguity": n_ambig,
        "n_undefined": n_undefined,
        "pre_length_lengths": pre_length_lengths,
        **manual_stats,
    }
    return list(bypassed) + passed, fraction_used, qc_stats


def log_mad_cutoffs(
    lengths: list[int],
    k: float,
    min_lo_mult: float = 0.20,
    max_hi_mult: float = 5.0,
) -> tuple[float | None, float | None, float, float]:
    """Compute length-outlier cutoffs from MAD on log-lengths, with a
    hard-floor safety net.

    Why log-space: protein/CDS length is right-skewed and bounded at zero, so
    a symmetric ±k·σ window in log-space becomes a multiplicative window
    around the median — the natural shape for length data.

    Why MAD: the standard deviation of log-lengths would itself be inflated
    by the very outliers we want to detect.  MAD (median absolute deviation),
    scaled by 1.4826, is a robust estimator of σ for normal data.

    Why the floor: in tight families (sequences within a few percent of the
    median) MAD on log-lengths gives a very narrow keep window that would
    cut legitimate moderately-truncated variants.  The floor guarantees a
    minimum keep window of ``[min_lo_mult, max_hi_mult] × median``, taking
    the **union** of the MAD window and the floor.  Set either floor knob
    to 0 to disable that side.

    Returns ``(lo_cutoff, hi_cutoff, sigma_log, median_len)``.  Returns
    ``(None, None, 0.0, median_len)`` when both MAD and floor are disabled.
    """
    valid = [L for L in lengths if L > 0]
    if not valid:
        return None, None, 0.0, 0.0
    median_len = float(statistics.median(valid))

    floor_lo = median_len * min_lo_mult if min_lo_mult > 0 else None
    floor_hi = median_len * max_hi_mult if max_hi_mult > 0 else None

    sigma = 0.0
    mad_lo: float | None = None
    mad_hi: float | None = None
    if k > 0:
        log_lens = [math.log(L) for L in valid]
        med = statistics.median(log_lens)
        mad = statistics.median(abs(x - med) for x in log_lens)
        if mad > 0:
            sigma = 1.4826 * mad
            mad_lo = math.exp(med - k * sigma)
            mad_hi = math.exp(med + k * sigma)

    # Take the union: MAD widens the floor when the data warrants it, but
    # the floor never narrows a wider MAD window.
    if mad_lo is None:
        final_lo = floor_lo
    elif floor_lo is None:
        final_lo = mad_lo
    else:
        final_lo = min(mad_lo, floor_lo)

    if mad_hi is None:
        final_hi = floor_hi
    elif floor_hi is None:
        final_hi = mad_hi
    else:
        final_hi = max(mad_hi, floor_hi)

    return final_lo, final_hi, sigma, median_len


def remove_length_outliers(
    records: list[SeqRecord],
    k: float = 5.0,
    min_lo_mult: float = 0.20,
    max_hi_mult: float = 5.0,
    protected_ids: set[str] | None = None,
) -> tuple[list[SeqRecord], int, int]:
    """Remove sequences whose length is a MAD-on-log outlier.

    Args:
        records: input sequences
        k: number of MAD-units (in log-space, scaled by 1.4826) on either
           side of the log-median that defines the keep window.  A larger
           k is more lenient.  k <= 0 disables the MAD test.
        min_lo_mult: hard-floor lower bound — sequences ≥ ``min_lo_mult ×
                     median`` are never dropped, regardless of MAD.  0
                     disables the lower floor.
        max_hi_mult: hard-floor upper bound — sequences ≤ ``max_hi_mult ×
                     median`` are never dropped, regardless of MAD.  0
                     disables the upper floor.
        protected_ids: record IDs that must never be dropped even if flagged —
                       a warning is logged instead.

    Returns:
        (filtered records, n_long_removed, n_short_removed)
    """
    if len(records) < 2:
        return records, 0, 0
    protected = protected_ids or set()
    lengths = [len(r.seq) for r in records]
    lo_cutoff, hi_cutoff, sigma, median_len = log_mad_cutoffs(
        lengths, k, min_lo_mult=min_lo_mult, max_hi_mult=max_hi_mult,
    )
    if lo_cutoff is None and hi_cutoff is None:
        return records, 0, 0

    passed: list[SeqRecord] = []
    n_long = 0
    n_short = 0
    for r in records:
        L = len(r.seq)
        too_long = hi_cutoff is not None and L > hi_cutoff
        too_short = lo_cutoff is not None and L < lo_cutoff
        if too_long or too_short:
            if r.id in protected:
                log.warning(
                    "RefSeq '%s' looks like a length outlier "
                    "(length=%d, median=%d, lo_cutoff=%s, hi_cutoff=%s) — "
                    "KEEPING (RefSeq protected)",
                    r.id, L, int(median_len),
                    f"{lo_cutoff:.0f}" if lo_cutoff is not None else "disabled",
                    f"{hi_cutoff:.0f}" if hi_cutoff is not None else "disabled",
                )
                passed.append(r)
            elif too_long:
                n_long += 1
            else:
                n_short += 1
        else:
            passed.append(r)

    n_kept = len(passed)
    n_in = len(records)
    lo_str = f"{lo_cutoff:.0f}" if lo_cutoff is not None else "—"
    hi_str = f"{hi_cutoff:.0f}" if hi_cutoff is not None else "—"
    log.info(
        "Length-outlier filter: kept %d/%d sequences "
        "(dropped %d short + %d long). "
        "median length=%d aa/nt, keep window=[%s, %s] "
        "(MAD-on-log k=%.1f, floor=[%.2fx, %.2fx])",
        n_kept, n_in, n_short, n_long, int(median_len),
        lo_str, hi_str, k, min_lo_mult, max_hi_mult,
    )
    return passed, n_long, n_short


def _apply_length_filter(
    records: list[SeqRecord],
    min_length: int,
    ambig_chars: set,
    max_ambiguous: float,
    max_length: int | None = None,
) -> tuple[list[SeqRecord], int, int, int, list[int]]:
    """Filter records by ambiguity then length.

    Returns (passed, n_short, n_ambig, n_undefined, pre_length_lengths).
    pre_length_lengths: lengths of sequences that passed the ambiguity check,
    before the minimum-length check is applied.
    """
    passed = []
    n_excluded_length = 0
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
        L = len(seq_str)
        if L < min_length or (max_length is not None and L > max_length):
            n_excluded_length += 1
            continue
        passed.append(rec)
    return passed, n_excluded_length, n_ambig, n_undefined, pre_length_lengths


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
