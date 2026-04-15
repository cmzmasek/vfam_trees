"""Per-family PDF report generation."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .logger import get_logger

log = get_logger(__name__)


def generate_family_report(
    family: str,
    output_pdf: Path,
    summary_row: dict,
    seq_lengths: list[int],
    tree_support: dict[str, list[float]],
    bio_trees: dict[str, Any],
) -> None:
    """Generate a per-family PDF report with stats and plots.

    Args:
        family: viral family name
        output_pdf: path for the output PDF
        summary_row: dict from build_summary_row (for the stats table)
        seq_lengths: list of sequence lengths passing QC
        tree_support: {"500": [sh_values], "100": [sh_values]}
        bio_trees: {"500": BioPython tree, "100": BioPython tree}
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from Bio import Phylo
    except ImportError as e:
        log.warning("PDF report skipped — matplotlib or biopython not available: %s", e)
        return

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    with PdfPages(str(output_pdf)) as pdf:

        # ------------------------------------------------------------------
        # Page 1: Summary statistics table
        # ------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")

        title = f"vfam_trees v{__version__}  —  {family}"
        fig.suptitle(title, fontsize=14, fontweight="bold", y=0.97)
        ax.text(0.5, 0.96, f"Generated {timestamp}", ha="center", va="top",
                transform=ax.transAxes, fontsize=9, color="gray")

        table_data = [
            ["NCBI taxid",          str(summary_row.get("ncbi_taxid", ""))],
            ["Lineage",             _wrap(str(summary_row.get("lineage", "")), 80)],
            ["Molecule / region",   str(summary_row.get("molecule_region", ""))],
            ["Species discovered",  str(summary_row.get("species_discovered", ""))],
            ["Species with seqs",   str(summary_row.get("species_with_seqs", ""))],
            ["Species relaxed QC",  str(summary_row.get("species_relaxed_threshold", ""))],
            ["Seqs passing QC",     str(summary_row.get("seqs_passing_qc", ""))],
            ["QC excl. organism",   str(summary_row.get("qc_excluded_organism", ""))],
            ["QC excl. length",     str(summary_row.get("qc_excluded_length", ""))],
            ["QC excl. ambiguity",  str(summary_row.get("qc_excluded_ambiguity", ""))],
            ["QC undefined seq",    str(summary_row.get("qc_undefined", ""))],
            ["Seq len min",         str(summary_row.get("seqlen_min", ""))],
            ["Seq len median",      str(summary_row.get("seqlen_median", ""))],
            ["Seq len max",         str(summary_row.get("seqlen_max", ""))],
            ["Seq len IQR",         str(summary_row.get("seqlen_iqr", ""))],
        ]
        for label in ("500", "100"):
            prefix = f"tree{label}"
            sup = "shlike" if label == "500" else "shalrt"
            table_data += [
                [f"Tree {label} — leaves",         str(summary_row.get(f"{prefix}_leaves", ""))],
                [f"Tree {label} — MSA length",     str(summary_row.get(f"{prefix}_msa_length", ""))],
                [f"Tree {label} — MSA gap %",      str(summary_row.get(f"{prefix}_msa_gap_pct", ""))],
                [f"Tree {label} — cluster thresh", f"{summary_row.get(f'{prefix}_cluster_thresh_min', '')}–{summary_row.get(f'{prefix}_cluster_thresh_max', '')}"],
                [f"Tree {label} — SH support median", str(summary_row.get(f"{prefix}_{sup}_median", ""))],
                [f"Tree {label} — SH support IQR",    str(summary_row.get(f"{prefix}_{sup}_iqr", ""))],
            ]

        tbl = ax.table(
            cellText=table_data,
            colLabels=["Parameter", "Value"],
            cellLoc="left",
            loc="upper center",
            bbox=[0.0, 0.0, 1.0, 0.92],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9)
        for (r, c), cell in tbl.get_celld().items():
            if r == 0:
                cell.set_facecolor("#3c6e9f")
                cell.set_text_props(color="white", fontweight="bold")
            elif r % 2 == 0:
                cell.set_facecolor("#f0f4f8")
            cell.set_edgecolor("#cccccc")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ------------------------------------------------------------------
        # Page 2: Sequence length histogram
        # ------------------------------------------------------------------
        if seq_lengths:
            fig, ax = plt.subplots(figsize=(9, 5))
            ax.hist(seq_lengths, bins=40, color="#3c6e9f", edgecolor="white", linewidth=0.5)
            ax.set_xlabel("Sequence length (bp/aa)", fontsize=11)
            ax.set_ylabel("Count", fontsize=11)
            ax.set_title(f"{family} — sequence length distribution (n={len(seq_lengths)})",
                         fontsize=12)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ------------------------------------------------------------------
        # Page 3: SH support histograms
        # ------------------------------------------------------------------
        has_500 = bool(tree_support.get("500"))
        has_100 = bool(tree_support.get("100"))
        if has_500 or has_100:
            n_panels = (1 if has_500 else 0) + (1 if has_100 else 0)
            fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels + 1, 5),
                                     squeeze=False)
            panel = 0
            colors = {"500": "#3c6e9f", "100": "#e07b39"}
            labels = {"500": "SH-like (FastTree)", "100": "SH-aLRT (IQ-TREE)"}
            for key in ("500", "100"):
                vals = tree_support.get(key, [])
                if not vals:
                    continue
                ax = axes[0][panel]
                ax.hist(vals, bins=20, range=(0, 100), color=colors[key],
                        edgecolor="white", linewidth=0.5)
                ax.set_xlabel("Support value", fontsize=11)
                ax.set_ylabel("Number of internal nodes", fontsize=11)
                ax.set_title(f"{family} tree_{key}\n{labels[key]} (n={len(vals)})",
                             fontsize=11)
                ax.set_xlim(0, 100)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                panel += 1
            fig.suptitle(f"{family} — branch support distributions", fontsize=12,
                         fontweight="bold")
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ------------------------------------------------------------------
        # Page 4: Tree_100 visualization
        # ------------------------------------------------------------------
        tree_100 = bio_trees.get("100")
        if tree_100 is not None:
            try:
                n_leaves = sum(1 for _ in tree_100.get_terminals())
                fig_h = max(6, n_leaves * 0.18)
                fig, ax = plt.subplots(figsize=(11, min(fig_h, 24)))
                Phylo.draw(tree_100, axes=ax, do_show=False)
                ax.set_title(f"{family} tree_100", fontsize=11, fontweight="bold")
                # Remove individual leaf labels if too many (>60) to keep readable
                if n_leaves > 60:
                    for label in ax.get_yticklabels():
                        label.set_fontsize(5)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
            except Exception as e:
                log.warning("Tree visualization skipped for %s: %s", family, e)

        # PDF metadata
        d = pdf.infodict()
        d["Title"] = f"vfam_trees report — {family}"
        d["Author"] = f"vfam_trees v{__version__}"
        d["CreationDate"] = datetime.now(timezone.utc)

    log.info("PDF report written to %s", output_pdf)


def _wrap(text: str, width: int) -> str:
    """Wrap a long string at word boundaries."""
    if len(text) <= width:
        return text
    words = text.split()
    lines = []
    current = ""
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(current)
    return "\n".join(lines)
