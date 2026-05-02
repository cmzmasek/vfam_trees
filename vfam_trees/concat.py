"""Concatenated multi-marker alignment helpers.

Pure helpers (per-marker length-outlier, concat builder, partition-file
writer) are testable without external tools.  The ``align_and_trim_markers``
orchestrator is a thin wrapper around the existing ``run_msa`` and
``run_trim`` and gets exercised in integration when the pipeline runs.

CONCAT_DESIGN.md §5.5–§5.6, §5.10.
"""
from __future__ import annotations

import re
import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Per-marker length-outlier (RefSeq exempt)
# CONCAT_DESIGN.md §5.5
# ---------------------------------------------------------------------------

def remove_per_marker_length_outliers(
    genomes: dict[str, dict[str, SeqRecord]],
    refseq_genome_ids: set[str],
    hi_mult: float = 3.0,
    lo_mult: float = 1.0 / 3.0,
) -> tuple[dict[str, dict[str, SeqRecord]], dict]:
    """Drop ``(genome, marker)`` cells whose protein length is a 2-sided
    outlier within that marker's length distribution across all genomes.

    RefSeq genomes (whose source-nuc accession appears in
    ``refseq_genome_ids``) are exempt: flagged outliers stay and a warning
    is logged.  The genome remains in the dataset; the dropped marker
    becomes a gap-padded block in the concat.

    Args:
        genomes:           {genome_id: {marker_name: SeqRecord}}
        refseq_genome_ids: source-nuc accessions immune to dropping
        hi_mult:           drop seqs longer than hi_mult × median (0 disables)
        lo_mult:           drop seqs shorter than lo_mult × median (0 disables)

    Returns:
        (updated_genomes, stats):
            updated_genomes — deep-copy with offending cells omitted
            stats — {n_long_dropped, n_short_dropped, per_marker_median, n_refseq_protected}
    """
    marker_lengths: dict[str, list[int]] = {}
    for mark_to_rec in genomes.values():
        for marker_name, rec in mark_to_rec.items():
            marker_lengths.setdefault(marker_name, []).append(len(rec.seq))

    medians: dict[str, float] = {}
    for marker_name, lengths in marker_lengths.items():
        if lengths:
            medians[marker_name] = float(statistics.median(lengths))

    n_long = 0
    n_short = 0
    n_refseq_protected = 0
    updated: dict[str, dict[str, SeqRecord]] = {}

    for genome_id, mark_to_rec in genomes.items():
        is_refseq = genome_id in refseq_genome_ids
        new_markers: dict[str, SeqRecord] = {}
        for marker_name, rec in mark_to_rec.items():
            median = medians.get(marker_name, 0.0)
            length = len(rec.seq)
            hi_cut = median * hi_mult if (hi_mult and hi_mult > 0 and median) else None
            lo_cut = median * lo_mult if (lo_mult and lo_mult > 0 and median) else None
            too_long = hi_cut is not None and length > hi_cut
            too_short = lo_cut is not None and length < lo_cut

            if too_long or too_short:
                if is_refseq:
                    log.warning(
                        "RefSeq genome %s: marker %r length=%d looks like a "
                        "length outlier (median=%.0f, hi_cutoff=%s, "
                        "lo_cutoff=%s) — KEEPING (RefSeq protected)",
                        genome_id, marker_name, length, median,
                        f"{hi_cut:.0f}" if hi_cut is not None else "disabled",
                        f"{lo_cut:.0f}" if lo_cut is not None else "disabled",
                    )
                    new_markers[marker_name] = rec
                    n_refseq_protected += 1
                elif too_long:
                    n_long += 1
                else:
                    n_short += 1
            else:
                new_markers[marker_name] = rec
        updated[genome_id] = new_markers

    if n_long or n_short:
        log.info(
            "Per-marker length-outlier removal: dropped %d long + %d short "
            "(genome × marker) cells (RefSeq-protected: %d)",
            n_long, n_short, n_refseq_protected,
        )

    stats = {
        "n_long_dropped":     n_long,
        "n_short_dropped":    n_short,
        "n_refseq_protected": n_refseq_protected,
        "per_marker_median":  medians,
    }
    return updated, stats


# ---------------------------------------------------------------------------
# Concatenation + partition map
# CONCAT_DESIGN.md §5.6
# ---------------------------------------------------------------------------

def concatenate_aligned_markers(
    aligned_per_marker: dict[str, dict[str, SeqRecord]],
    genome_ids: list[str],
    marker_order: list[str],
) -> tuple[dict[str, SeqRecord], dict[str, tuple[int, int]]]:
    """Stitch per-marker aligned MSAs into one concatenation per genome.

    Genomes lacking a given marker contribute an all-gap block of the
    marker's aligned column count, preserving alignment integrity.

    Args:
        aligned_per_marker: {marker_name: {genome_id: aligned SeqRecord}}.
                            All records within one marker must have equal
                            length (i.e. the per-marker MSA must be aligned).
        genome_ids:         universe of genome IDs to include in the concat.
                            Order is preserved in the output dict.
        marker_order:       order in which marker blocks are emitted.

    Returns:
        (concat, partition_map):
            concat:        {genome_id: concat SeqRecord}
            partition_map: {marker_name: (start_col, end_col)} 1-based,
                           inclusive — NEXUS-style coordinates suitable for
                           IQ-TREE ``-p partitions.nex``.
    """
    block_lengths: dict[str, int] = {}
    for marker in marker_order:
        per_genome = aligned_per_marker.get(marker, {})
        if not per_genome:
            block_lengths[marker] = 0
            continue
        lengths = {len(rec.seq) for rec in per_genome.values()}
        if len(lengths) > 1:
            raise ValueError(
                f"Aligned MSA for marker {marker!r} has inconsistent column "
                f"counts {lengths} — alignment may be corrupted."
            )
        block_lengths[marker] = lengths.pop()

    partition_map: dict[str, tuple[int, int]] = {}
    cursor = 1
    for marker in marker_order:
        block_len = block_lengths[marker]
        if block_len == 0:
            continue
        partition_map[marker] = (cursor, cursor + block_len - 1)
        cursor += block_len

    concat: dict[str, SeqRecord] = {}
    for genome_id in genome_ids:
        parts: list[str] = []
        for marker in marker_order:
            block_len = block_lengths[marker]
            if block_len == 0:
                continue
            rec = aligned_per_marker.get(marker, {}).get(genome_id)
            if rec is None:
                parts.append("-" * block_len)
            else:
                parts.append(str(rec.seq))
        concat[genome_id] = SeqRecord(Seq("".join(parts)), id=genome_id, description="")

    return concat, partition_map


# ---------------------------------------------------------------------------
# NEXUS partition file
# CONCAT_DESIGN.md §5.10
# ---------------------------------------------------------------------------

def write_partition_file_nexus(
    partition_map: dict[str, tuple[int, int]],
    output_path: Path,
) -> None:
    """Write a NEXUS partition file for IQ-TREE ``-p partitions.nex``.

    Marker names are sanitized to comply with NEXUS charset naming rules
    (alphanumerics + underscore, no whitespace).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#nexus", "begin sets;"]
    for marker, (start, end) in partition_map.items():
        safe = _safe_charset_name(marker)
        lines.append(f"    charset {safe} = {start}-{end};")
    lines.append("end;")
    output_path.write_text("\n".join(lines) + "\n")


def _safe_charset_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_") or "marker"


# ---------------------------------------------------------------------------
# Per-marker MSA + trim orchestrator
# CONCAT_DESIGN.md §5.6
# ---------------------------------------------------------------------------

def align_and_trim_markers(
    genomes: dict[str, dict[str, SeqRecord]],
    marker_order: list[str],
    work_dir: Path,
    msa_tool: str = "mafft",
    msa_options: str = "--auto",
    trim_enabled: bool = True,
    trim_tool: str = "trimal",
    trim_options: str = "-automated1",
    threads: int = 1,
) -> dict[str, dict[str, SeqRecord]]:
    """Per-marker MSA + trim, returning aligned records keyed by genome_id.

    For each marker, writes ``<work_dir>/<marker>/raw.fasta`` then runs
    MAFFT to ``aln.fasta``, then trimAl (if enabled) to ``aln.trim.fasta``.
    Markers with fewer than 2 sequences (across all genomes) are skipped
    and yield an empty per-marker dict — the concatenation step gap-pads
    accordingly.
    """
    from .msa import run_msa
    from .trim import run_trim

    aligned_per_marker: dict[str, dict[str, SeqRecord]] = {}

    for marker in marker_order:
        marker_dir = work_dir / _safe_charset_name(marker)
        marker_dir.mkdir(parents=True, exist_ok=True)

        records_for_marker: list[SeqRecord] = []
        for genome_id, mark_to_rec in genomes.items():
            rec = mark_to_rec.get(marker)
            if rec is None:
                continue
            records_for_marker.append(SeqRecord(rec.seq, id=genome_id, description=""))

        if len(records_for_marker) < 2:
            log.warning(
                "Marker %r has only %d genome(s) with this protein — "
                "skipping per-marker alignment (concat will gap-pad).",
                marker, len(records_for_marker),
            )
            aligned_per_marker[marker] = {}
            continue

        raw_fasta = marker_dir / "raw.fasta"
        with open(raw_fasta, "w") as f:
            SeqIO.write(records_for_marker, f, "fasta")

        aln_fasta = marker_dir / "aln.fasta"
        run_msa(raw_fasta, aln_fasta, tool=msa_tool, options=msa_options, threads=threads)

        if trim_enabled:
            trim_fasta = marker_dir / "aln.trim.fasta"
            run_trim(aln_fasta, trim_fasta, tool=trim_tool, options=trim_options)
            final = trim_fasta
        else:
            final = aln_fasta

        aligned_per_marker[marker] = {
            rec.id: rec for rec in SeqIO.parse(final, "fasta")
        }

    return aligned_per_marker
