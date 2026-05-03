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
from .subsample import (
    absorb_into_refseqs,
    adaptive_cluster_species,
    proportional_merge,
)

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Genome identity helpers — RefSeq detection + raw concat building
# CONCAT_DESIGN.md §5.7–§5.9
# ---------------------------------------------------------------------------

def is_refseq_genome(genome_id: str) -> bool:
    """Return True iff *genome_id* is a RefSeq nucleotide accession.

    Matches the NCBI RefSeq prefix convention: two uppercase letters
    followed by an underscore (NC_, NZ_, AC_, NW_, NT_, etc.).  Bare
    GenBank accessions (e.g. MN908947.3) return False.
    """
    acc = (genome_id or "").strip()
    return (
        len(acc) >= 3
        and acc[2] == "_"
        and acc[0:2].isalpha()
        and acc[0:2].isupper()
    )


def identify_refseq_genomes(
    species_genomes: dict[str, dict[str, dict[str, SeqRecord]]],
) -> set[str]:
    """Return the set of source-nuc accessions that are RefSeq genomes.

    Args:
        species_genomes: {species_name: {genome_id: {marker_name: SeqRecord}}}

    Returns:
        set of genome_ids (source-nuc accessions) matching the RefSeq prefix.
    """
    return {
        genome_id
        for genomes in species_genomes.values()
        for genome_id in genomes
        if is_refseq_genome(genome_id)
    }


def build_raw_concat(
    genome_proteins: dict[str, SeqRecord],
    marker_order: list[str],
    genome_id: str,
) -> SeqRecord:
    """Build an unaligned head-to-tail concatenation of one genome's markers.

    Used as the input sequence for genome-level clustering and RefSeq
    absorption.  Missing markers contribute nothing (no gap padding) — the
    concat is shorter for genomes with low coverage, which is fine for
    similarity-based clustering and biases absorption toward genomes that
    actually share most markers.

    Args:
        genome_proteins: {marker_name: SeqRecord} for one genome.
        marker_order:    canonical marker order from the family's spec.
        genome_id:       source-nuc accession; used as the SeqRecord id.

    Returns:
        SeqRecord with id=genome_id and seq = the head-to-tail concat.
    """
    parts: list[str] = []
    for marker in marker_order:
        rec = genome_proteins.get(marker)
        if rec is not None:
            parts.append(str(rec.seq))
    return SeqRecord(Seq("".join(parts)), id=genome_id, description="")


# ---------------------------------------------------------------------------
# Genome-level clustering + RefSeq absorption + cross-species merge
# CONCAT_DESIGN.md §5.7–§5.9.  Operates on the per-genome raw concat
# sequences built by ``build_raw_concat`` — MMseqs2 clustering only needs
# similarity, not alignment, so we work on unaligned head-to-tail concats
# to avoid running family-scale MAFFT before clustering reduces the set.
# ---------------------------------------------------------------------------

def cluster_and_merge_genomes(
    species_genomes: dict[str, dict[str, dict[str, SeqRecord]]],
    marker_order: list[str],
    target_n: int,
    max_reps_per_species: int,
    threshold_min: float,
    threshold_max: float,
    work_dir: Path,
    clustering_tool: str = "mmseqs2",
    refseq_absorption_enabled: bool = True,
    refseq_absorption_threshold: float = 0.99,
    seed: int = 42,
) -> tuple[set[str], dict]:
    """Reduce per-genome diversity to the target tree size.

    For each species:
      1. Build raw (unaligned, head-to-tail) concat per genome.
      2. RefSeq absorption on those concats — drops non-RefSeq genomes
         that are ≥ threshold identical to a RefSeq concat in the same
         species.
      3. Adaptive clustering on the absorbed concats; each species
         contributes at most ``max_reps_per_species`` representative
         genomes.

    Then a single cross-species proportional merge picks ``target_n``
    genomes overall, with RefSeq genomes preferred at the per-species
    quota stage.

    Returns:
        (selected_genome_ids, stats):
            selected_genome_ids — set of source-nuc accessions chosen.
            stats — diagnostic counters (n_refseq_absorbed, n_clusters_per_species,
                    cluster_thresholds_used).
    """
    refseq_ids = identify_refseq_genomes(species_genomes)

    species_reps: dict[str, list[SeqRecord]] = {}
    n_refseq_absorbed_total = 0
    thresholds_used: list[float] = []

    for sp_name, genomes in species_genomes.items():
        if not genomes:
            continue
        sp_safe = sp_name.replace(" ", "_").replace("/", "_")
        sp_work = work_dir / sp_safe

        # 1. Build raw concat per genome for this species
        concats = [
            build_raw_concat(genome_proteins, marker_order, genome_id)
            for genome_id, genome_proteins in genomes.items()
        ]

        # 2. RefSeq absorption (concat-level)
        if refseq_absorption_enabled and concats:
            concats, n_absorbed = absorb_into_refseqs(
                records=concats,
                refseq_ids=refseq_ids,
                threshold=refseq_absorption_threshold,
                seq_type="protein",
                work_dir=sp_work / "absorb",
                clustering_tool=clustering_tool,
            )
            if n_absorbed:
                log.info(
                    "RefSeq absorption (concat): %s — %d non-RefSeq genome(s) "
                    "absorbed at identity ≥ %.2f",
                    sp_name, n_absorbed, refseq_absorption_threshold,
                )
            n_refseq_absorbed_total += n_absorbed

        if not concats:
            continue

        # 3. Adaptive clustering — concats are protein-like sequences
        n_input = len(concats)
        reps, threshold_used = adaptive_cluster_species(
            records=concats,
            max_reps=max_reps_per_species,
            threshold_min=threshold_min,
            threshold_max=threshold_max,
            seq_type="protein",
            work_dir=sp_work / "cluster",
            clustering_tool=clustering_tool,
        )
        species_reps[sp_name] = reps
        thresholds_used.append(threshold_used)
        if n_input != len(reps):
            log.info(
                "Clustering %s: %d genomes → %d representatives at identity %.4f",
                sp_name, n_input, len(reps), threshold_used,
            )

    # Cross-species proportional merge — RefSeq priority preserved
    merged = proportional_merge(
        species_reps, target_n, seed=seed, priority_ids=refseq_ids,
    )
    selected = {rec.id for rec in merged}

    stats = {
        "n_refseq_absorbed":     n_refseq_absorbed_total,
        "n_species_with_reps":   len(species_reps),
        "n_total_reps":          sum(len(reps) for reps in species_reps.values()),
        "n_selected":            len(selected),
        "cluster_thresh_min":    round(min(thresholds_used), 4) if thresholds_used else "",
        "cluster_thresh_max":    round(max(thresholds_used), 4) if thresholds_used else "",
    }
    return selected, stats


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

    # Always log the per-marker median sample so users can spot annotation
    # issues even when no cells are dropped.
    if medians:
        median_str = ", ".join(
            f"{m}:{int(medians[m])}(n={len(marker_lengths[m])})"
            for m in medians
        )
        log.info(
            "Per-marker length-outlier scan: per-marker medians [%s] — "
            "dropped %d long + %d short (genome × marker) cells "
            "(RefSeq-protected: %d)",
            median_str, n_long, n_short, n_refseq_protected,
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
    n_markers = len(marker_order)

    for i, marker in enumerate(marker_order, start=1):
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
                "Marker [%d/%d] %r has only %d genome(s) with this protein — "
                "skipping per-marker alignment (concat will gap-pad).",
                i, n_markers, marker, len(records_for_marker),
            )
            aligned_per_marker[marker] = {}
            continue

        raw_fasta = marker_dir / "raw.fasta"
        with open(raw_fasta, "w") as f:
            SeqIO.write(records_for_marker, f, "fasta")

        aln_fasta = marker_dir / "aln.fasta"
        run_msa(raw_fasta, aln_fasta, tool=msa_tool, options=msa_options, threads=threads)

        # Count MAFFT output columns for the per-marker progress line.
        aln_cols_pre = next(
            (len(rec.seq) for rec in SeqIO.parse(aln_fasta, "fasta")), 0,
        )

        if trim_enabled:
            trim_fasta = marker_dir / "aln.trim.fasta"
            run_trim(aln_fasta, trim_fasta, tool=trim_tool, options=trim_options)
            final = trim_fasta
            aln_cols_post = next(
                (len(rec.seq) for rec in SeqIO.parse(final, "fasta")), 0,
            )
            log.info(
                "Marker [%d/%d] %r: %d genome(s) → MAFFT → %d cols → "
                "trimAl → %d cols",
                i, n_markers, marker, len(records_for_marker),
                aln_cols_pre, aln_cols_post,
            )
        else:
            final = aln_fasta
            log.info(
                "Marker [%d/%d] %r: %d genome(s) → MAFFT → %d cols (trim disabled)",
                i, n_markers, marker, len(records_for_marker), aln_cols_pre,
            )

        aligned_per_marker[marker] = {
            rec.id: rec for rec in SeqIO.parse(final, "fasta")
        }

    return aligned_per_marker
