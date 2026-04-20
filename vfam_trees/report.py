"""Per-family PDF report generation."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .logger import get_logger
from .summary import compute_seqlen_stats

log = get_logger(__name__)


def generate_family_report(
    family: str,
    output_pdf: Path,
    summary_row: dict,
    seq_lengths: list[int],
    tree_seq_lengths: dict[str, list[int]] | None = None,
    tree_support: dict[str, list[float]] | None = None,
    bio_trees: dict[str, Any] | None = None,
    tree_leaf_colors: dict[str, dict] | None = None,
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
    if tree_seq_lengths is None:
        tree_seq_lengths = {}
    if tree_support is None:
        tree_support = {}
    if bio_trees is None:
        bio_trees = {}
    if tree_leaf_colors is None:
        tree_leaf_colors = {}

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
            tree_opts = str(summary_row.get(f"{prefix}_tree_options", "")).strip()
            msa_opts  = str(summary_row.get(f"{prefix}_msa_options", "")).strip()
            tsl = compute_seqlen_stats(tree_seq_lengths.get(label, []))
            table_data += [
                [f"Tree {label} — sequence type",    str(summary_row.get(f"{prefix}_seq_type", ""))],
                [f"Tree {label} — MSA tool",         str(summary_row.get(f"{prefix}_msa_tool", ""))],
                [f"Tree {label} — MSA options",      msa_opts or "—"],
                [f"Tree {label} — tree program",     str(summary_row.get(f"{prefix}_tree_tool", ""))],
                [f"Tree {label} — tree model",       str(summary_row.get(f"{prefix}_tree_model", ""))],
                [f"Tree {label} — tree options",     tree_opts or "—"],
                [f"Tree {label} — leaves",           str(summary_row.get(f"{prefix}_leaves", ""))],
                [f"Tree {label} — seq len min",      str(tsl.get("min", ""))],
                [f"Tree {label} — seq len median",   str(tsl.get("median", ""))],
                [f"Tree {label} — seq len max",      str(tsl.get("max", ""))],
                [f"Tree {label} — seq len IQR",      str(tsl.get("iqr", ""))],
                [f"Tree {label} — MSA length",       str(summary_row.get(f"{prefix}_msa_length", ""))],
                [f"Tree {label} — MSA gap %",        str(summary_row.get(f"{prefix}_msa_gap_pct", ""))],
                [f"Tree {label} — cluster thresh",   f"{summary_row.get(f'{prefix}_cluster_thresh_min', '')}–{summary_row.get(f'{prefix}_cluster_thresh_max', '')}"],
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
        # Pages 3a/3b: Sequence length histograms for tree_500 and tree_100
        # ------------------------------------------------------------------
        for label, color in (("500", "#3c6e9f"), ("100", "#e07b39")):
            lengths = tree_seq_lengths.get(label, [])
            if lengths:
                fig, ax = plt.subplots(figsize=(9, 5))
                ax.hist(lengths, bins=min(40, max(5, len(lengths) // 3)),
                        color=color, edgecolor="white", linewidth=0.5)
                ax.set_xlabel("Sequence length (bp/aa)", fontsize=11)
                ax.set_ylabel("Count", fontsize=11)
                ax.set_title(
                    f"{family} tree_{label} — sequence length distribution (n={len(lengths)})",
                    fontsize=12,
                )
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        # ------------------------------------------------------------------
        # Page 4: SH support histograms
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
            colors_100 = tree_leaf_colors.get("100", {})
            fig = _draw_tree_fig(
                tree_100, family, "100",
                label_colors=colors_100.get("display_to_color"),
                genus_to_color=colors_100.get("genus_to_color"),
                subfamily_to_genera=colors_100.get("subfamily_to_genera"),
            )
            if fig is not None:
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        # PDF metadata
        d = pdf.infodict()
        d["Title"] = f"vfam_trees report — {family}"
        d["Author"] = f"vfam_trees v{__version__}"
        d["CreationDate"] = datetime.now(timezone.utc)

    log.info("PDF report written to %s", output_pdf)


def save_tree_images(
    family: str,
    output_dir: Path,
    bio_trees: dict[str, Any],
    tree_leaf_colors: dict[str, dict] | None = None,
) -> None:
    """Save standalone PDF and PNG images of the tree_100 and tree_500 visualizations."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        log.warning("Tree image export skipped — matplotlib not available: %s", e)
        return

    if tree_leaf_colors is None:
        tree_leaf_colors = {}

    output_dir.mkdir(parents=True, exist_ok=True)

    for label in ("100", "500"):
        tree = bio_trees.get(label)
        if tree is None:
            continue
        colors = tree_leaf_colors.get(label, {})
        fig = _draw_tree_fig(
            tree, family, label,
            label_colors=colors.get("display_to_color"),
            genus_to_color=colors.get("genus_to_color"),
            subfamily_to_genera=colors.get("subfamily_to_genera"),
        )
        if fig is None:
            continue
        stem = f"{family}_tree_{label}"
        for path, kwargs in [
            (output_dir / f"{stem}.pdf", {}),
            (output_dir / f"{stem}.png", {"dpi": 150}),
        ]:
            try:
                fig.savefig(str(path), bbox_inches="tight", **kwargs)
                log.info("Tree image written to %s", path)
            except Exception as e:
                log.warning("Could not write %s: %s", path, e)
        plt.close(fig)


def _draw_tree_fig(
    tree,
    family: str,
    label: str = "100",
    label_colors: dict[str, str] | None = None,
    genus_to_color: dict[str, str] | None = None,
    subfamily_to_genera: dict[str, list[str]] | None = None,
):
    """Return a matplotlib Figure of the tree, or None on error."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.text as _mtext
        from Bio import Phylo
    except ImportError:
        return None

    try:
        n_leaves = sum(1 for _ in tree.get_terminals())
        fig_h = max(6, n_leaves * 0.18)

        has_legend = bool(genus_to_color)
        fig_w = 14 if has_legend else 11
        fig, ax = plt.subplots(figsize=(fig_w, min(fig_h, 24)))

        # BioPython's default label_func is str(clade), which truncates names
        # longer than 40 chars to "name[:37] + '...'".  Since our display names
        # include species|strain|accession|host and are almost always >40 chars,
        # we must override label_func to render the full name.
        Phylo.draw(
            tree, axes=ax, do_show=False,
            label_func=lambda c: c.name or "",
        )
        _thin_tree_lines(ax)
        ax.axis("off")
        ax.set_title(f"{family} tree_{label}", fontsize=11, fontweight="bold")
        font_size = max(4, min(8, int(200 / max(n_leaves, 1))))

        n_colored = 0
        for artist in ax.get_children():
            if not isinstance(artist, _mtext.Text):
                continue
            if artist is ax.title:
                continue
            artist.set_fontsize(font_size)
            if label_colors:
                name = artist.get_text().strip()
                color = label_colors.get(name)
                if color:
                    artist.set_color(color)
                    n_colored += 1

        log.debug(
            "tree_%s: colored %d / %d leaf labels in %s",
            label, n_colored, n_leaves, family,
        )

        if has_legend and genus_to_color:
            _draw_taxonomy_legend(ax, genus_to_color, subfamily_to_genera or {})

        return fig
    except Exception as e:
        log.warning("Tree visualization skipped for %s: %s", family, e)
        return None


def _thin_tree_lines(ax, scale: float = 0.5) -> None:
    """Halve (or scale) the line width of every Line2D artist on the axes.

    BioPython's Phylo.draw uses matplotlib's default linewidth (1.5 pt) for
    the tree branches; that reads heavy in dense trees.  Call immediately
    after Phylo.draw.
    """
    for line in ax.get_lines():
        line.set_linewidth(line.get_linewidth() * scale)


def _draw_taxonomy_legend(
    ax,
    genus_to_color: dict[str, str],
    subfamily_to_genera: dict[str, list[str]],
) -> None:
    """Add a genus/subfamily color legend to the axes."""
    import matplotlib.patches as mpatches

    handles = []
    # Build legend entries grouped by subfamily
    placed: set[str] = set()

    if subfamily_to_genera:
        for subfamily, genera in sorted(subfamily_to_genera.items()):
            # Subfamily header (bold, no patch)
            handles.append(mpatches.Patch(color="none", label=f"  {subfamily}"))
            for genus in sorted(genera):
                color = genus_to_color.get(genus, "#888888")
                handles.append(mpatches.Patch(color=color, label=f"    {genus}"))
                placed.add(genus)
    # Any remaining genera without a subfamily
    for genus, color in sorted(genus_to_color.items()):
        if genus not in placed:
            handles.append(mpatches.Patch(color=color, label=genus))

    if not handles:
        return

    n = len(handles)
    ncol = max(1, n // 30)
    legend = ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=6,
        frameon=True,
        framealpha=0.8,
        ncol=ncol,
        title="Genus (by subfamily)",
        title_fontsize=7,
        borderpad=0.5,
        handlelength=1.0,
        handletextpad=0.4,
    )


# Viral-group → (background color, legend label).
# Keyed by substring match against the lowercased NCBI lineage string.
# Iteration order is deliberate: more specific keys (kingdoms, phyla) are
# checked before broad realm keys so that, e.g., a Negarnaviricota virus is
# tagged "–ssRNA" instead of the catch-all "Riboviria" (+ssRNA/dsRNA).
# Multiple keys may map to the same (color, label) so new NCBI realm names
# (e.g. ICTV 2023's "Floreoviria" for ssDNA) can be added without losing
# coverage of older names.
_REALM_ENTRIES: list[tuple[str, str, str]] = [
    # (lineage substring, hex color, legend label)
    ("pararnavirae",    "#fff8e1", "RT viruses"),        # kingdom under Riboviria
    ("negarnaviricota", "#fce4ec", "–ssRNA"),            # phylum under Riboviria
    ("duplodnaviria",   "#e3f2fd", "dsDNA (tailed)"),    # realm
    ("varidnaviria",    "#dceefb", "dsDNA (non-tailed)"),# realm
    ("monodnaviria",    "#f3e5f5", "ssDNA"),             # realm (legacy)
    ("floreoviria",     "#f3e5f5", "ssDNA"),             # realm (ICTV 2023+)
    ("singelaviria",    "#f3e5f5", "ssDNA"),             # realm (ICTV 2023+)
    ("ribozyviria",     "#ede7f6", "deltavirus-like"),   # realm
    ("riboviria",       "#e8f5e9", "+ssRNA / dsRNA"),    # realm (catch-all last)
]
_DEFAULT_BG = "#f5f5f5"


def _match_realm(lineage_str: str) -> tuple[str, str | None]:
    """Return (bg_color, legend_label) for a lineage string.

    legend_label is None when no entry matched (caller should use defaults).
    """
    l = lineage_str.lower()
    for key, color, label in _REALM_ENTRIES:
        if key in l:
            return color, label
    return _DEFAULT_BG, None


def _read_summary_lineages(output_dir: Path) -> dict[str, str]:
    """Return {family: lineage_str} from summary.tsv if present."""
    import csv as _csv
    summary = output_dir / "summary.tsv"
    if not summary.exists():
        return {}
    result: dict[str, str] = {}
    try:
        with open(summary, newline="") as f:
            reader = _csv.DictReader(f, delimiter="\t")
            for row in reader:
                fam = row.get("family", "")
                lin = row.get("lineage", "")
                if fam:
                    result[fam] = lin
    except Exception:
        pass
    return result


def generate_overview_png(output_dir: Path, output_path: Path) -> None:
    """Thumbnail grid of every tree_100 found under output_dir.

    Scans output_dir for */*_tree_100.nwk files, draws each as a small
    topology-only thumbnail (no leaf labels), arranges them in a grid,
    and saves the result as a PNG.
    """
    try:
        import copy
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from Bio import Phylo
    except ImportError as e:
        log.warning("Overview PNG skipped — matplotlib/biopython not available: %s", e)
        return

    nwk_files = sorted(output_dir.glob("*/*_tree_100.nwk"))
    if not nwk_files:
        log.warning("No tree_100.nwk files found in %s — skipping overview PNG.", output_dir)
        return

    lineage_map = _read_summary_lineages(output_dir)

    entries: list[tuple[str, object]] = []
    for nwk_path in nwk_files:
        family = nwk_path.stem.replace("_tree_100", "")
        try:
            tree = next(iter(Phylo.parse(str(nwk_path), "newick")))
            entries.append((family, tree))
        except Exception as e:
            log.warning("Could not parse %s for overview: %s", nwk_path, e)

    if not entries:
        log.warning("No parseable tree_100 files found — skipping overview PNG.")
        return

    n = len(entries)
    n_cols = min(5, n)
    n_rows = (n + n_cols - 1) // n_cols

    thumb_w, thumb_h = 3.0, 3.5
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * thumb_w, n_rows * thumb_h))

    import numpy as np
    ax_flat = np.array(axes).flatten() if n > 1 else [axes]

    used_groups: dict[str, str] = {}  # label → color (for legend)

    def _strip_axis_decor(ax) -> None:
        # Hide spines/ticks/axis-labels without calling ax.axis("off"):
        # axison=False would skip the facecolor patch during rendering.
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax.set_xlabel("")
        ax.set_ylabel("")

    for idx, (family, tree) in enumerate(entries):
        ax = ax_flat[idx]
        lineage = lineage_map.get(family, "")
        bg_color, group_label = _match_realm(lineage)
        if group_label is not None:
            used_groups[group_label] = bg_color
        else:
            log.warning(
                "Overview PNG: no viral-group match for %s (lineage: %r) — "
                "using default background.",
                family, lineage,
            )
        ax.set_facecolor(bg_color)
        _strip_axis_decor(ax)
        try:
            tree_copy = copy.deepcopy(tree)
            for clade in tree_copy.get_terminals():
                clade.name = ""
            Phylo.draw(tree_copy, axes=ax, do_show=False)
            _thin_tree_lines(ax)
            for txt in ax.texts:
                txt.set_visible(False)
            # Phylo.draw re-enables axis decor and may reset the patch;
            # re-apply after it has drawn.
            ax.set_facecolor(bg_color)
            _strip_axis_decor(ax)
            ax.set_title(family, fontsize=6, pad=3, fontweight="bold")
        except Exception as e:
            ax.set_title(family, fontsize=6)
            log.warning("Could not draw overview thumbnail for %s: %s", family, e)

    # Hide unused axes completely (no background needed)
    for idx in range(len(entries), len(ax_flat)):
        ax_flat[idx].axis("off")

    fig.suptitle("tree_100 overview", fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()

    # Add viral-group color legend if lineage data was available
    if used_groups:
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(color=c, label=lbl) for lbl, c in sorted(used_groups.items())]
        fig.legend(handles=handles, loc="lower center", ncol=len(handles),
                   fontsize=7, frameon=True, title="Viral group",
                   title_fontsize=7, bbox_to_anchor=(0.5, -0.02))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        log.info("Overview PNG written to %s", output_path)
    except Exception as e:
        log.warning("Could not write overview PNG: %s", e)
    plt.close(fig)


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
