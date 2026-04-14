"""Per-family summary TSV — appended after each family completes."""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree as BioTree

from .logger import get_logger

log = get_logger(__name__)

COLUMNS = [
    "family",
    "ncbi_taxid",
    "lineage",
    "molecule_region",
    "species_discovered",
    "species_with_seqs",
    # sequence length stats (post-QC)
    "seqlen_min",
    "seqlen_q1",
    "seqlen_median",
    "seqlen_mean",
    "seqlen_q3",
    "seqlen_max",
    "seqlen_iqr",
    # tree_500
    "tree500_leaves",
    "tree500_bs_min",
    "tree500_bs_q1",
    "tree500_bs_median",
    "tree500_bs_q3",
    "tree500_bs_max",
    "tree500_bs_iqr",
    "tree500_msa_length",
    "tree500_msa_gap_pct",
    # tree_100
    "tree100_leaves",
    "tree100_bs_min",
    "tree100_bs_q1",
    "tree100_bs_median",
    "tree100_bs_q3",
    "tree100_bs_max",
    "tree100_bs_iqr",
    "tree100_msa_length",
    "tree100_msa_gap_pct",
]


def compute_seqlen_stats(lengths: list[int]) -> dict:
    """Return summary statistics for a list of sequence lengths.

    Returns a dict with keys: min, q1, median, mean, q3, max, iqr.
    """
    if not lengths:
        return {k: "" for k in ("min", "q1", "median", "mean", "q3", "max", "iqr")}
    vals = sorted(lengths)
    mean = round(statistics.mean(vals), 1)
    if len(vals) >= 3:
        q1, median, q3 = statistics.quantiles(vals, n=4)
        iqr = q3 - q1
    else:
        q1 = median = q3 = statistics.median(vals)
        iqr = 0.0
    return {
        "min":    vals[0],
        "q1":     round(q1, 1),
        "median": round(median, 1),
        "mean":   mean,
        "q3":     round(q3, 1),
        "max":    vals[-1],
        "iqr":    round(iqr, 1),
    }


def compute_bootstrap_stats(tree: BioTree) -> dict:
    """Return bootstrap summary statistics for internal nodes of a tree.

    Only internal nodes with a non-None confidence value are included.
    Returns a dict with keys: min, q1, median, q3, max, iqr.
    All values are rounded to 1 decimal place.
    """
    vals = sorted(
        c.confidence
        for c in tree.find_clades()
        if not c.is_terminal() and c.confidence is not None
    )
    if not vals:
        return {k: "" for k in ("min", "q1", "median", "q3", "max", "iqr")}

    if len(vals) >= 3:
        q1, median, q3 = statistics.quantiles(vals, n=4)
        iqr = q3 - q1
    else:
        q1 = median = q3 = statistics.median(vals)
        iqr = 0.0
    return {
        "min":    round(vals[0], 1),
        "q1":     round(q1, 1),
        "median": round(median, 1),
        "q3":     round(q3, 1),
        "max":    round(vals[-1], 1),
        "iqr":    round(iqr, 1),
    }


def compute_msa_stats(msa_fasta: Path) -> dict:
    """Return alignment length and overall gap percentage for a FASTA MSA.

    Returns a dict with keys: length (int), gap_pct (float, 0–100, 1 d.p.).
    """
    records = list(SeqIO.parse(str(msa_fasta), "fasta"))
    if not records:
        return {"length": "", "gap_pct": ""}

    aln_len = len(records[0].seq)
    total_chars = aln_len * len(records)
    total_gaps = sum(str(r.seq).count("-") for r in records)
    gap_pct = round(100.0 * total_gaps / total_chars, 1) if total_chars else 0.0
    return {"length": aln_len, "gap_pct": gap_pct}


def format_molecule_region(seq_type: str, region: str, segment: str | None) -> str:
    """Build a human-readable molecule/region string for the summary."""
    mol = "protein" if seq_type == "protein" else "nucleotide"
    if segment:
        region_str = segment
    elif region == "whole_genome":
        region_str = "whole genome"
    else:
        region_str = f"gene: {region}"
    return f"{mol}, {region_str}"


def write_summary_row(summary_path: Path, row: dict) -> None:
    """Append one row to the summary TSV, writing the header if the file is new."""
    is_new = not summary_path.exists() or summary_path.stat().st_size == 0
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "a", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore"
        )
        if is_new:
            writer.writeheader()
        writer.writerow(row)
    log.info("Summary updated: %s", summary_path)


def build_summary_row(
    family: str,
    family_taxid: int | None,
    family_lineage: list[dict],
    seq_type: str,
    region: str,
    segment: str | None,
    n_species_discovered: int,
    n_species_with_seqs: int,
    seqlen_stats: dict,
    tree_stats: dict[str, dict],
) -> dict:
    """Assemble a summary row dict from collected pipeline stats.

    tree_stats should be keyed by label ("500", "100"), each value a dict with:
        leaves, bs (bootstrap stats dict), msa (msa stats dict).
    """
    lineage_str = "; ".join(e["name"] for e in family_lineage) if family_lineage else ""

    row: dict = {
        "family":             family,
        "ncbi_taxid":         family_taxid if family_taxid is not None else "",
        "lineage":            lineage_str,
        "molecule_region":    format_molecule_region(seq_type, region, segment),
        "species_discovered": n_species_discovered,
        "species_with_seqs":  n_species_with_seqs,
        "seqlen_min":         seqlen_stats.get("min", ""),
        "seqlen_q1":          seqlen_stats.get("q1", ""),
        "seqlen_median":      seqlen_stats.get("median", ""),
        "seqlen_mean":        seqlen_stats.get("mean", ""),
        "seqlen_q3":          seqlen_stats.get("q3", ""),
        "seqlen_max":         seqlen_stats.get("max", ""),
        "seqlen_iqr":         seqlen_stats.get("iqr", ""),
    }

    for label, prefix in (("500", "tree500"), ("100", "tree100")):
        stats = tree_stats.get(label, {})
        bs = stats.get("bs", {k: "" for k in ("min", "q1", "median", "q3", "max", "iqr")})
        msa = stats.get("msa", {"length": "", "gap_pct": ""})
        row[f"{prefix}_leaves"]     = stats.get("leaves", "")
        row[f"{prefix}_bs_min"]     = bs.get("min", "")
        row[f"{prefix}_bs_q1"]      = bs.get("q1", "")
        row[f"{prefix}_bs_median"]  = bs.get("median", "")
        row[f"{prefix}_bs_q3"]      = bs.get("q3", "")
        row[f"{prefix}_bs_max"]     = bs.get("max", "")
        row[f"{prefix}_bs_iqr"]     = bs.get("iqr", "")
        row[f"{prefix}_msa_length"] = msa.get("length", "")
        row[f"{prefix}_msa_gap_pct"] = msa.get("gap_pct", "")

    return row
