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
    "species_relaxed_threshold",
    "seqs_passing_qc",
    # QC exclusion breakdown
    "qc_excluded_organism",
    "qc_excluded_length",
    "qc_excluded_ambiguity",
    "qc_undefined",
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
    "tree500_support_type",
    "tree500_support_min",
    "tree500_support_q1",
    "tree500_support_median",
    "tree500_support_q3",
    "tree500_support_max",
    "tree500_support_iqr",
    "tree500_cluster_thresh_min",
    "tree500_cluster_thresh_max",
    "tree500_msa_length",
    "tree500_msa_length_pre_trim",
    "tree500_trim_tool",
    "tree500_trim_options",
    "tree500_msa_gap_pct",
    # tree_100
    "tree100_leaves",
    "tree100_support_type",
    "tree100_support_min",
    "tree100_support_q1",
    "tree100_support_median",
    "tree100_support_q3",
    "tree100_support_max",
    "tree100_support_iqr",
    "tree100_cluster_thresh_min",
    "tree100_cluster_thresh_max",
    "tree100_msa_length",
    "tree100_msa_length_pre_trim",
    "tree100_trim_tool",
    "tree100_trim_options",
    "tree100_msa_gap_pct",
    # diversity / outlier counts
    "tree500_n_outliers_removed",
    "tree100_n_outliers_removed",
    "tree500_n_length_outliers_long",
    "tree500_n_length_outliers_short",
    "tree100_n_length_outliers_long",
    "tree100_n_length_outliers_short",
    "tree500_n_genera",
    "tree100_n_genera",
    "tree500_n_subfamilies",
    "tree100_n_subfamilies",
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


def compute_support_stats(tree: BioTree) -> dict:
    """Return branch-support summary statistics for internal nodes of a tree.

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
    if aln_len == 0:
        return {"length": 0, "gap_pct": 0.0}
    total_chars = aln_len * len(records)
    total_gaps = sum(str(r.seq).count("-") for r in records)
    gap_pct = round(100.0 * total_gaps / total_chars, 1)
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
    n_species_relaxed: int = 0,
    total_seqs_qc: int = 0,
    qc_stats: dict | None = None,
) -> dict:
    """Assemble a summary row dict from collected pipeline stats.

    tree_stats should be keyed by label ("500", "100"), each value a dict with:
        leaves, support (branch-support stats dict), msa (msa stats dict).
    """
    lineage_str = "; ".join(e["name"] for e in family_lineage) if family_lineage else ""

    row: dict = {
        "family":             family,
        "ncbi_taxid":         family_taxid if family_taxid is not None else "",
        "lineage":            lineage_str,
        "molecule_region":    format_molecule_region(seq_type, region, segment),
        "species_discovered":        n_species_discovered,
        "species_with_seqs":         n_species_with_seqs,
        "species_relaxed_threshold": n_species_relaxed,
        "seqs_passing_qc":           total_seqs_qc,
        "qc_excluded_organism":      (qc_stats or {}).get("n_excluded_organism", ""),
        "qc_excluded_length":        (qc_stats or {}).get("n_excluded_length", ""),
        "qc_excluded_ambiguity":     (qc_stats or {}).get("n_excluded_ambiguity", ""),
        "qc_undefined":              (qc_stats or {}).get("n_undefined", ""),
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
        sh = stats.get("support", {k: "" for k in ("min", "q1", "median", "q3", "max", "iqr")})
        msa = stats.get("msa", {"length": "", "gap_pct": ""})
        row[f"{prefix}_leaves"]            = stats.get("leaves", "")
        row[f"{prefix}_support_type"]      = stats.get("support_type", "")
        row[f"{prefix}_support_min"]       = sh.get("min", "")
        row[f"{prefix}_support_q1"]        = sh.get("q1", "")
        row[f"{prefix}_support_median"]    = sh.get("median", "")
        row[f"{prefix}_support_q3"]        = sh.get("q3", "")
        row[f"{prefix}_support_max"]       = sh.get("max", "")
        row[f"{prefix}_support_iqr"]       = sh.get("iqr", "")
        row[f"{prefix}_cluster_thresh_min"]   = stats.get("cluster_thresh_min", "")
        row[f"{prefix}_cluster_thresh_max"]   = stats.get("cluster_thresh_max", "")
        row[f"{prefix}_msa_length"]           = msa.get("length", "")
        row[f"{prefix}_msa_length_pre_trim"]  = stats.get("msa_length_pre_trim", "")
        row[f"{prefix}_trim_tool"]            = stats.get("trim_tool", "")
        row[f"{prefix}_trim_options"]         = stats.get("trim_options", "")
        row[f"{prefix}_msa_gap_pct"]          = msa.get("gap_pct", "")
        row[f"{prefix}_seq_type"]             = stats.get("seq_type", "")
        row[f"{prefix}_msa_tool"]             = stats.get("msa_tool", "")
        row[f"{prefix}_msa_options"]          = stats.get("msa_options", "")
        row[f"{prefix}_tree_tool"]            = stats.get("tree_tool", "")
        row[f"{prefix}_tree_model"]           = stats.get("tree_model", "")
        row[f"{prefix}_tree_options"]         = stats.get("tree_options", "")
        row[f"{prefix}_n_outliers_removed"]      = stats.get("n_outliers_removed", "")
        row[f"{prefix}_n_length_outliers_long"]  = stats.get("n_length_outliers_long", "")
        row[f"{prefix}_n_length_outliers_short"] = stats.get("n_length_outliers_short", "")
        row[f"{prefix}_n_genera"]                = stats.get("n_genera", "")
        row[f"{prefix}_n_subfamilies"]           = stats.get("n_subfamilies", "")

    return row
