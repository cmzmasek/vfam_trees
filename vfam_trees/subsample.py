"""Per-species adaptive clustering and proportional merging."""
from __future__ import annotations

import random
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)

BINARY_SEARCH_ITERATIONS = 7   # ~128 resolution steps over threshold range


def adaptive_cluster_species(
    records: list[SeqRecord],
    max_reps: int,
    threshold_min: float,
    threshold_max: float,
    seq_type: str,
    work_dir: Path,
    clustering_tool: str = "mmseqs2",
) -> tuple[list[SeqRecord], float]:
    """Find the tightest identity threshold that yields ≤ max_reps for one species.

    Uses binary search over [threshold_min, threshold_max].
    - If threshold_max already gives ≤ max_reps: use it (no clustering needed).
    - If threshold_min still gives > max_reps: use threshold_min anyway (best effort).
    - Otherwise: binary search to find highest threshold with ≤ max_reps.

    Args:
        records: quality-filtered SeqRecords for one species (short IDs)
        max_reps: maximum number of representatives desired
        threshold_min: lower bound for binary search (most aggressive clustering)
        threshold_max: upper bound for binary search (least aggressive clustering)
        seq_type: 'nucleotide' or 'protein'
        work_dir: temp directory for clustering runs
        clustering_tool: 'mmseqs2' (default) or 'cdhit'

    Returns:
        (representative_records, threshold_used)
    """
    if len(records) <= max_reps:
        log.debug("Species has %d sequences (≤ max_reps %d) — no clustering needed.",
                  len(records), max_reps)
        return records, threshold_max

    # Check if even the most aggressive threshold is sufficient
    reps_at_min = _cluster_at(records, threshold_min, seq_type,
                               work_dir / f"t{threshold_min:.2f}", clustering_tool)
    if len(reps_at_min) > max_reps:
        log.debug("threshold_min=%.2f still gives %d reps > max_reps %d — using threshold_min.",
                  threshold_min, len(reps_at_min), max_reps)
        return reps_at_min, threshold_min

    # Check if no clustering is needed at threshold_max
    reps_at_max = _cluster_at(records, threshold_max, seq_type,
                               work_dir / f"t{threshold_max:.2f}", clustering_tool)
    if len(reps_at_max) <= max_reps:
        log.debug("threshold_max=%.2f gives %d reps ≤ max_reps %d — using threshold_max.",
                  threshold_max, len(reps_at_max), max_reps)
        return reps_at_max, threshold_max

    # Binary search: find highest threshold where n_reps ≤ max_reps
    lo = threshold_min
    hi = threshold_max
    best_reps = reps_at_min
    best_threshold = threshold_min

    for _ in range(BINARY_SEARCH_ITERATIONS):
        mid = round((lo + hi) / 2, 4)
        reps = _cluster_at(records, mid, seq_type,
                            work_dir / f"t{mid:.4f}", clustering_tool)
        if len(reps) <= max_reps:
            best_reps = reps
            best_threshold = mid
            lo = mid   # can try higher (less aggressive)
        else:
            hi = mid   # need lower (more aggressive)

    log.debug(
        "Adaptive clustering: %d sequences → %d reps at threshold %.4f (max_reps=%d)",
        len(records), len(best_reps), best_threshold, max_reps,
    )
    return best_reps, best_threshold


def proportional_merge(
    species_reps: dict[str, list[SeqRecord]],
    target: int,
    seed: int = 42,
) -> list[SeqRecord]:
    """Select sequences proportionally across species to hit target count.

    Each species contributes slots proportional to its cluster count.
    Species with zero representatives are skipped.
    The total is trimmed or padded (using leftover slots) to hit exactly target.

    Args:
        species_reps: mapping of species_name → list of representative SeqRecords
        target: desired total number of sequences
        seed: random seed for reproducibility

    Returns:
        Flat list of selected SeqRecords (length ≤ target).
    """
    random.seed(seed)

    non_empty = {sp: reps for sp, reps in species_reps.items() if reps}
    if not non_empty:
        return []

    total_reps = sum(len(reps) for reps in non_empty.values())

    if total_reps <= target:
        log.info("Total representatives (%d) ≤ target (%d), using all.", total_reps, target)
        return [rec for reps in non_empty.values() for rec in reps]

    # Compute per-species quota proportional to cluster count
    quotas: dict[str, int] = {}
    for sp, reps in non_empty.items():
        quotas[sp] = max(1, round(len(reps) / total_reps * target))

    # Adjust total to exactly target (rounding may cause off-by-one)
    quota_total = sum(quotas.values())
    diff = target - quota_total
    if diff != 0:
        fractions = sorted(
            non_empty.keys(),
            key=lambda sp: (len(non_empty[sp]) / total_reps * target) % 1,
            reverse=True,
        )
        for sp in fractions:
            if diff == 0:
                break
            if diff > 0:
                quotas[sp] += 1
                diff -= 1
            elif quotas[sp] > 1:
                quotas[sp] -= 1
                diff += 1

    selected: list[SeqRecord] = []
    for sp, reps in non_empty.items():
        quota = quotas.get(sp, 0)
        if len(reps) <= quota:
            selected.extend(reps)
        else:
            selected.extend(random.sample(reps, quota))

    log.info(
        "Proportional merge: %d species, %d total reps → %d selected (target %d)",
        len(non_empty), total_reps, len(selected), target,
    )
    return selected


# ---------------------------------------------------------------------------
# Internal clustering helpers
# ---------------------------------------------------------------------------

def _cluster_at(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
    clustering_tool: str,
) -> list[SeqRecord]:
    """Cluster records at a given threshold and return representative SeqRecords."""
    tool = clustering_tool.lower().replace("-", "").replace("_", "")
    if tool == "mmseqs2":
        rep_ids = _mmseqs2_cluster(records, threshold, seq_type, work_dir)
    elif tool in ("cdhit", "cdhitest"):
        rep_ids = _cdhit_cluster(records, threshold, seq_type, work_dir)
    else:
        log.warning("Unknown clustering tool '%s', falling back to mmseqs2", clustering_tool)
        rep_ids = _mmseqs2_cluster(records, threshold, seq_type, work_dir)

    id_to_rec = {rec.id: rec for rec in records}
    return [id_to_rec[rid] for rid in rep_ids if rid in id_to_rec]


def _mmseqs2_cluster(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
) -> list[str]:
    work_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = work_dir / "input.fasta"
    output_prefix = work_dir / "output"
    tmp_dir = work_dir / "tmp"

    with open(input_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    dbtype = "2" if seq_type == "nucleotide" else "1"
    cmd = [
        "mmseqs", "easy-linclust",
        str(input_fasta), str(output_prefix), str(tmp_dir),
        "--min-seq-id", str(threshold),
        "--dbtype", dbtype,
        "--threads", "1",
        "-v", "0",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.warning("MMseqs2 failed at threshold %.4f, falling back to CD-HIT", threshold)
        return _cdhit_cluster(records, threshold, seq_type, work_dir)

    rep_fasta = Path(str(output_prefix) + "_rep_seq.fasta")
    if not rep_fasta.exists():
        log.warning("MMseqs2 produced no output at threshold %.4f, falling back to CD-HIT", threshold)
        return _cdhit_cluster(records, threshold, seq_type, work_dir)

    return [rec.id for rec in SeqIO.parse(rep_fasta, "fasta")]


def _cdhit_cluster(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
) -> list[str]:
    work_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = work_dir / "input.fasta"
    output_fasta = work_dir / "output.fasta"

    with open(input_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    if seq_type == "protein":
        cmd = ["cd-hit", "-i", str(input_fasta), "-o", str(output_fasta),
               "-c", str(threshold), "-n", _word_size_protein(threshold),
               "-T", "1", "-M", "0", "-d", "0"]
    else:
        cmd = ["cd-hit-est", "-i", str(input_fasta), "-o", str(output_fasta),
               "-c", str(threshold), "-n", _word_size_nuc(threshold),
               "-T", "1", "-M", "0", "-d", "0"]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not output_fasta.exists():
        log.error("CD-HIT failed at threshold %.4f — using all sequences", threshold)
        return [rec.id for rec in records]

    return [rec.id for rec in SeqIO.parse(output_fasta, "fasta")]


def _word_size_protein(threshold: float) -> str:
    if threshold >= 0.7:
        return "5"
    elif threshold >= 0.6:
        return "4"
    elif threshold >= 0.5:
        return "3"
    return "2"


def _word_size_nuc(threshold: float) -> str:
    if threshold >= 0.9:
        return "8"
    elif threshold >= 0.88:
        return "7"
    elif threshold >= 0.85:
        return "6"
    return "5"
