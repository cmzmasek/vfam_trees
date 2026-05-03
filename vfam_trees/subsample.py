"""Per-species adaptive clustering and proportional merging."""
from __future__ import annotations

import random
import re
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)

BINARY_SEARCH_ITERATIONS = 7   # ~128 resolution steps over threshold range


def absorb_into_refseqs(
    records: list[SeqRecord],
    refseq_ids: set[str],
    threshold: float,
    seq_type: str,
    work_dir: Path,
    clustering_tool: str = "mmseqs2",
) -> tuple[list[SeqRecord], int]:
    """Drop non-RefSeq sequences that are near-identical to a RefSeq in the same set.

    Clusters *records* at *threshold* identity; for each cluster that contains
    one or more RefSeq members, drops the non-RefSeq members (they are deemed
    near-identical to the RefSeq and therefore redundant).  All RefSeqs are
    always kept; clusters containing no RefSeq are passed through unchanged.

    A new MMseqs2/CD-HIT run is executed only when *records* contains at least
    one RefSeq; otherwise the function returns the input untouched.

    Args:
        records:    quality-filtered SeqRecords for one species (short IDs).
        refseq_ids: set of short IDs that are RefSeqs (immune from removal).
        threshold:  clustering identity (0–1); typical 0.99 = "near identical".
        seq_type:   'nucleotide' or 'protein'.
        work_dir:   temp directory for the clustering run.
        clustering_tool: 'mmseqs2' (default) or 'cdhit'.

    Returns:
        (kept_records, n_absorbed)
    """
    if len(records) < 2:
        return records, 0
    record_ids = {r.id for r in records}
    if not (record_ids & refseq_ids):
        return records, 0

    clusters = _cluster_membership(records, threshold, seq_type, work_dir, clustering_tool)
    keep_ids: set[str] = set()
    n_absorbed = 0
    for members in clusters:
        cluster_refseqs = members & refseq_ids
        if cluster_refseqs:
            keep_ids |= cluster_refseqs
            n_absorbed += len(members) - len(cluster_refseqs)
        else:
            keep_ids |= members

    if n_absorbed:
        log.info(
            "RefSeq absorption: %d non-RefSeq sequence(s) absorbed into "
            "RefSeq cluster reps at identity ≥ %.2f",
            n_absorbed, threshold,
        )
    return [r for r in records if r.id in keep_ids], n_absorbed


def _cluster_membership(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
    clustering_tool: str,
) -> list[set[str]]:
    """Cluster *records* and return the full membership of each cluster.

    Each item in the returned list is a set of record IDs belonging to one
    cluster (representative + members).
    """
    tool = clustering_tool.lower().replace("-", "").replace("_", "")
    if tool == "mmseqs2":
        return _mmseqs2_membership(records, threshold, seq_type, work_dir)
    if tool in ("cdhit", "cdhitest"):
        return _cdhit_membership(records, threshold, seq_type, work_dir)
    raise ValueError(
        f"Unsupported clustering tool: {clustering_tool!r}. "
        "Supported: 'mmseqs2', 'cdhit'."
    )


def _mmseqs2_membership(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
) -> list[set[str]]:
    """Run MMseqs2 easy-linclust and parse the cluster TSV for full membership."""
    _mmseqs2_cluster(records, threshold, seq_type, work_dir)
    cluster_tsv = work_dir / "output_cluster.tsv"
    if not cluster_tsv.exists():
        raise FileNotFoundError(
            f"MMseqs2 exited 0 but did not produce the expected cluster TSV "
            f"{cluster_tsv}. This likely indicates a change in MMseqs2's "
            "output-file naming conventions."
        )
    clusters: dict[str, set[str]] = {}
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters.setdefault(rep, set()).add(member)
    return list(clusters.values())


def _cdhit_membership(
    records: list[SeqRecord],
    threshold: float,
    seq_type: str,
    work_dir: Path,
) -> list[set[str]]:
    """Run CD-HIT and parse the .clstr file for full membership."""
    _cdhit_cluster(records, threshold, seq_type, work_dir)
    clstr_file = work_dir / "output.fasta.clstr"
    if not clstr_file.exists():
        raise FileNotFoundError(
            f"CD-HIT exited 0 but did not produce expected cluster file "
            f"{clstr_file}."
        )
    member_re = re.compile(r">([^\s]+?)\.\.\.")
    clusters: list[set[str]] = []
    current: set[str] = set()
    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = set()
            else:
                m = member_re.search(line)
                if m:
                    current.add(m.group(1))
    if current:
        clusters.append(current)
    return clusters


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
        log.info("Species has %d sequences (≤ max_reps %d) — no clustering needed.",
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
    priority_ids: set[str] | None = None,
) -> list[SeqRecord]:
    """Select sequences proportionally across species to hit target count.

    Each species contributes slots proportional to its cluster count.
    Species with zero representatives are skipped.
    The total is trimmed or padded (using leftover slots) to hit exactly target.

    Args:
        species_reps: mapping of species_name → list of representative SeqRecords
        target: desired total number of sequences
        seed: random seed for reproducibility
        priority_ids: optional set of record IDs to prefer when subsampling
                      within a species (e.g. RefSeq short IDs). When a
                      species' quota is smaller than its rep count, priority
                      records are taken first and the remainder is randomised.

    Returns:
        Flat list of selected SeqRecords (length ≤ target).
    """
    random.seed(seed)
    priority_ids = priority_ids or set()

    non_empty = {sp: reps for sp, reps in species_reps.items() if reps}
    if not non_empty:
        return []

    total_reps = sum(len(reps) for reps in non_empty.values())

    if total_reps <= target:
        log.info("Total representatives (%d) ≤ target (%d), using all.", total_reps, target)
        return [rec for reps in non_empty.values() for rec in reps]

    # When there are more species than slots, we cannot give every species a
    # slot.  Prefer species with more representatives (ties broken alphabetically
    # for reproducibility) and give each one exactly one slot.
    if len(non_empty) > target:
        sorted_species = sorted(
            non_empty.keys(),
            key=lambda s: (-len(non_empty[s]), s),
        )
        chosen_species = sorted_species[:target]
        selected: list[SeqRecord] = []
        n_priority_taken = 0
        n_priority_total = sum(
            1 for reps in non_empty.values() for r in reps if r.id in priority_ids
        )
        for sp in chosen_species:
            reps = non_empty[sp]
            prio = [r for r in reps if r.id in priority_ids]
            if prio:
                selected.append(random.choice(prio))
                n_priority_taken += 1
            else:
                selected.append(random.choice(reps))
        if priority_ids and n_priority_total:
            log.info(
                "Proportional merge: kept %d / %d priority record(s) (e.g. RefSeqs)",
                n_priority_taken, n_priority_total,
            )
        log.info(
            "Proportional merge: %d species (%d total reps) but target %d — "
            "selected 1 rep from %d largest species",
            len(non_empty), total_reps, target, len(chosen_species),
        )
        return selected

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
    n_priority_taken = 0
    n_priority_total = 0
    for sp, reps in non_empty.items():
        quota = quotas.get(sp, 0)
        if len(reps) <= quota:
            selected.extend(reps)
            n_priority_taken += sum(1 for r in reps if r.id in priority_ids)
            n_priority_total += sum(1 for r in reps if r.id in priority_ids)
        elif priority_ids:
            prio = [r for r in reps if r.id in priority_ids]
            rest = [r for r in reps if r.id not in priority_ids]
            n_priority_total += len(prio)
            if len(prio) >= quota:
                chosen = random.sample(prio, quota)
                selected.extend(chosen)
                n_priority_taken += quota
            else:
                selected.extend(prio)
                n_priority_taken += len(prio)
                remaining = quota - len(prio)
                selected.extend(random.sample(rest, remaining))
        else:
            selected.extend(random.sample(reps, quota))

    if priority_ids and n_priority_total:
        log.info(
            "Proportional merge: kept %d / %d priority record(s) (e.g. RefSeqs)",
            n_priority_taken, n_priority_total,
        )

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
        raise ValueError(
            f"Unsupported clustering tool: {clustering_tool!r}. "
            "Supported: 'mmseqs2', 'cdhit'."
        )

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
        log.error("MMseqs2 stderr:\n%s", (result.stderr or "")[-2000:])
        log.error("MMseqs2 stdout:\n%s", (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"MMseqs2 easy-linclust failed at threshold {threshold:.4f} "
            f"with exit code {result.returncode}. If MMseqs2 was upgraded "
            "recently, the CLI or default options may have changed — see "
            "stderr above. Consider setting clustering.tool: cdhit in the "
            "family config if this persists."
        )

    rep_fasta = Path(str(output_prefix) + "_rep_seq.fasta")
    if not rep_fasta.exists():
        raise FileNotFoundError(
            f"MMseqs2 exited 0 but did not produce the expected "
            f"representatives file {rep_fasta}. This likely indicates a "
            "change in MMseqs2's output-file naming conventions."
        )

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

    if result.returncode != 0:
        log.error("%s stderr:\n%s", cmd[0], (result.stderr or "")[-2000:])
        log.error("%s stdout:\n%s", cmd[0], (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"{cmd[0]} failed at threshold {threshold:.4f} with exit code "
            f"{result.returncode}. This may indicate a version "
            "incompatibility — see stderr above."
        )

    if not output_fasta.exists():
        raise FileNotFoundError(
            f"{cmd[0]} exited 0 but did not produce expected output "
            f"{output_fasta}. This may indicate a change in CD-HIT's output "
            "conventions."
        )

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
