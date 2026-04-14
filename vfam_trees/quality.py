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


def filter_sequences(
    records: list[SeqRecord],
    seq_type: str,
    min_length: int | None,
    max_ambiguous: float,
    exclude_organisms: list[str] | None = None,
) -> list[SeqRecord]:
    """Apply quality filters to a list of SeqRecords.

    Args:
        records: input sequences
        seq_type: 'nucleotide' or 'protein'
        min_length: minimum sequence length; None = auto (50% of median)
        max_ambiguous: maximum fraction of ambiguous characters
        exclude_organisms: list of substrings matched case-insensitively
                           against the organism name in annotations

    Returns:
        Filtered list of SeqRecords.
    """
    if not records:
        return []

    exclude_lower = [t.lower() for t in (exclude_organisms or [])]

    # Organism exclusion — applied before length filtering to keep median accurate
    passed_organism = []
    n_excluded = 0
    for rec in records:
        organism = rec.annotations.get("organism", "").lower()
        if any(term in organism for term in exclude_lower):
            log.debug("Excluding organism: %s", rec.annotations.get("organism", ""))
            n_excluded += 1
        else:
            passed_organism.append(rec)

    if n_excluded:
        log.info("Excluded %d sequences matching organism exclusion list", n_excluded)

    if not passed_organism:
        return []

    lengths = [len(r.seq) for r in passed_organism]
    median_len = statistics.median(lengths)

    if min_length is None:
        min_length = int(median_len * 0.5)
        log.info("Auto min_length set to %d (50%% of median %d)", min_length, median_len)
    else:
        log.info("Using min_length=%d", min_length)

    ambig_chars = AMBIGUOUS_NUC if seq_type == "nucleotide" else AMBIGUOUS_AA

    passed = []
    n_short = 0
    n_ambig = 0
    n_undefined = 0

    for rec in passed_organism:
        try:
            seq_str = str(rec.seq).upper()
        except Exception:
            # GenBank CONTIG records have no sequence content
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

    log.info(
        "Quality filter: %d passed, %d excluded organism, %d too short, "
        "%d too ambiguous, %d undefined sequence (from %d total)",
        len(passed), n_excluded, n_short, n_ambig, n_undefined, len(records),
    )
    return passed


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
