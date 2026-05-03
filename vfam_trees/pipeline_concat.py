"""Concatenated multi-marker pipeline runner.

Invoked by ``pipeline.run_family`` when ``sequence.region == "concatenated"``.
Per CONCAT_DESIGN.md:
  - §5.2: per-marker Entrez fetch with Policy A genome grouping
  - §5.5: per-marker length-outlier (RefSeq exempt)
  - §5.6: per-marker MAFFT + trimAl, concatenated alignment, partition map
  - §5.7–5.9: genome-level absorption, adaptive clustering, cross-species merge
  - §5.10: tree_500 single-model FastTree, tree_100 partitioned IQ-TREE

This first integration is intentionally minimal: it produces Newick and
PhyloXML trees plus summary/status rows, but defers tree-image rendering,
genus/subfamily coloring, and iterative branch-length outlier removal
to Phase 7.  The single-protein path in ``pipeline.run_family`` is
unchanged.
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
from typing import Any

from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord

from .branch_outliers import branch_length_stats, find_branch_length_outliers
from .colors import assign_leaf_colors
from .concat import (
    align_and_trim_markers,
    build_raw_concat,
    cluster_and_merge_genomes,
    concatenate_aligned_markers,
    identify_refseq_genomes,
    remove_per_marker_length_outliers,
    write_partition_file_nexus,
)
from .fetch import fetch_species_genomes, fetch_taxonomy_lineages
from .logger import get_logger
from .phyloxml_writer import write_phyloxml
from .report import generate_family_report, save_tree_icon, save_tree_images
from .summary import (
    build_status_row,
    build_summary_row,
    compute_msa_stats,
    compute_seqlen_stats,
    compute_support_stats,
    write_status_row,
    write_summary_row,
)
from .taxonomy import annotate_tree
from .tree import parse_iqtree_partition_models, run_tree, validate_newick

log = get_logger(__name__)


def run_family_concat(
    *,
    family: str,
    family_cfg: dict,
    family_taxid: int | None,
    family_lineage: list[dict],
    family_annotation: dict,
    family_dir: Path,
    work_dir: Path,
    output_dir: Path,
    species_list: list[dict],
    threads: int,
    summary_path: Path | None,
    status_path: Path | None,
    mark_done,
    mark_skipped,
) -> None:
    """Concat-mode entry point invoked from ``pipeline.run_family``.

    *mark_done* / *mark_skipped* are passed in by the caller so this module
    doesn't depend on private helpers from ``pipeline``.
    """
    seq_type = family_cfg["sequence"]["type"]
    region   = family_cfg["sequence"]["region"]
    segment  = family_cfg["sequence"].get("segment") or None  # always None in concat
    concat_cfg       = family_cfg["concatenation"]
    marker_set       = list(concat_cfg["proteins"])
    marker_order     = [m["name"] for m in marker_set]
    min_fraction     = float(concat_cfg.get("min_fraction", 0.7))
    quality_cfg      = family_cfg["quality"]
    clustering_cfg   = family_cfg["clustering"]
    clustering_tool  = clustering_cfg.get("tool", "mmseqs2")
    threshold_min    = float(clustering_cfg.get("threshold_min", 0.70))
    threshold_max    = float(clustering_cfg.get("threshold_max", 0.99))
    max_per_species  = int(family_cfg["download"]["max_per_species"])
    max_reps_500     = int(clustering_cfg.get("max_reps_500", 20))
    max_reps_100     = int(clustering_cfg.get("max_reps_100", 5))
    target_500       = int(family_cfg["targets"]["max_500"])
    target_100       = int(family_cfg["targets"]["max_100"])
    refseq_abs_cfg   = family_cfg.get("refseq_absorption", {}) or {}
    absorb_enabled   = bool(refseq_abs_cfg.get("enabled", True))
    absorb_threshold = float(refseq_abs_cfg.get("threshold", 0.99))
    length_outlier_cfg = family_cfg.get("length_outlier", {}) or {}
    hi_mult          = float(length_outlier_cfg.get("hi_mult", 3.0))
    lo_mult          = float(length_outlier_cfg.get("lo_mult", 1.0 / 3.0))

    log.info("=" * 60)
    log.info("Concatenated mode for %s: %d markers (%s)",
             family, len(marker_set), ", ".join(marker_order))
    log.info(
        "Downloading per-marker proteins for %d species "
        "(max %d non-RefSeq per species per marker, RefSeq uncapped) ...",
        len(species_list), max_per_species,
    )
    log.info("Sequence type: protein, concatenated (%d markers)", len(marker_set))
    log.info("=" * 60)

    # -------------------------------------------------------------------------
    # Per-species fetch
    # -------------------------------------------------------------------------
    species_genomes, species_lineages, fetch_stats = _fetch_all_species(
        family=family,
        species_list=species_list,
        marker_set=marker_set,
        work_dir=work_dir,
        max_per_species=max_per_species,
        min_fraction=min_fraction,
        exclude_organisms=quality_cfg.get("exclude_organisms"),
    )

    n_genomes_total = sum(len(g) for g in species_genomes.values())
    n_required_markers = max(1, int(round(min_fraction * len(marker_set))))
    log.info(
        "Concat fetch complete: %d proteins → %d source-nuc genome(s) found → "
        "%d kept after min_fraction (≥%d/%d markers); "
        "%d dropped (multi-accession / Policy A), "
        "%d dropped (low coverage), "
        "%d orphaned (no source-nuc accession) "
        "across %d species",
        fetch_stats["n_proteins_fetched"],
        fetch_stats["n_genomes_found"],
        n_genomes_total,
        n_required_markers, len(marker_set),
        fetch_stats["n_dropped_split_submission"],
        fetch_stats["n_dropped_min_fraction"],
        fetch_stats["n_orphaned_no_source"],
        len(species_genomes),
    )

    # Per-marker family-scale coverage summary — easy parity check vs the
    # marker preset's expected per-family yield.
    if species_genomes:
        per_marker_n = {m: 0 for m in marker_order}
        for genomes in species_genomes.values():
            for proteins in genomes.values():
                for m in proteins:
                    if m in per_marker_n:
                        per_marker_n[m] += 1
        coverage_str = ", ".join(f"{m}={per_marker_n[m]}" for m in marker_order)
        log.info("Per-marker coverage across %d genomes: %s",
                 n_genomes_total, coverage_str)

    # RefSeq detection — cross-family count for transparency.
    refseq_genome_ids = {
        gid for genomes in species_genomes.values() for gid in genomes
        if (len(gid) >= 3 and gid[2] == "_" and gid[:2].isalpha() and gid[:2].isupper())
    }
    if refseq_genome_ids:
        log.info(
            "Detected %d RefSeq genome(s) out of %d total — uncapped through "
            "clustering/merge and exempt from per-marker length-outlier and "
            "branch-length outlier removal",
            len(refseq_genome_ids), n_genomes_total,
        )

    if n_genomes_total < 4:
        reason = (f"too few genomes after concat fetch + min_fraction "
                  f"({n_genomes_total} < 4)")
        log.warning("%s — skipping.", reason)
        _emit_skip(family, family_taxid, family_lineage, family_annotation,
                   seq_type, region, segment, len(species_list),
                   len(species_genomes), reason,
                   summary_path, status_path,
                   family_dir, output_dir, mark_skipped)
        return

    # -------------------------------------------------------------------------
    # Per-target loop (tree_500 and tree_100)
    # -------------------------------------------------------------------------
    targets = [
        (target_500, max_reps_500, "500"),
        (target_100, max_reps_100, "100"),
    ]

    tree_stats: dict[str, dict] = {}
    for target_n, max_reps, label in targets:
        log.info("-" * 40)
        log.info("Building tree_%s (target=%d genomes, max_reps_per_species=%d)",
                 label, target_n, max_reps)
        target_stats = _run_target_concat(
            label=label,
            target_n=target_n,
            max_reps=max_reps,
            species_genomes=species_genomes,
            species_lineages=species_lineages,
            marker_set=marker_set,
            marker_order=marker_order,
            family=family,
            family_dir=family_dir,
            work_dir=work_dir,
            threads=threads,
            threshold_min=threshold_min,
            threshold_max=threshold_max,
            clustering_tool=clustering_tool,
            absorb_enabled=absorb_enabled,
            absorb_threshold=absorb_threshold,
            hi_mult=hi_mult,
            lo_mult=lo_mult,
            family_cfg=family_cfg,
        )
        tree_stats[label] = target_stats

    # -------------------------------------------------------------------------
    # Tree images + icon (PDF + PNG, rooted + unrooted) for concat trees.
    # Mirrors the single-protein output set so downstream consumers
    # (overview_tree_100.png, per-family report) work transparently.
    # -------------------------------------------------------------------------
    bio_trees = {
        label: stats["_display_tree"]
        for label, stats in tree_stats.items()
        if stats.get("_display_tree") is not None
    }
    tree_leaf_colors = {
        label: {
            "display_to_color":    stats.get("_display_to_color", {}),
            "genus_to_color":      stats.get("_genus_to_color", {}),
            "subfamily_to_genera": stats.get("_subfamily_to_genera", {}),
        }
        for label, stats in tree_stats.items()
    }
    visualization_cfg = (family_cfg.get("visualization") or {})
    branch_linewidth = float(visualization_cfg.get("branch_linewidth", 0.5))
    icon_cfg = (family_cfg.get("icon") or {})
    icon_size         = int(icon_cfg.get("size", 256))
    icon_bg_color     = icon_cfg.get("bg_color", "#EAF3F2")
    icon_branch_color = icon_cfg.get("branch_color", "#000000")

    if bio_trees:
        save_tree_images(
            family=family,
            output_dir=family_dir,
            bio_trees=bio_trees,
            tree_leaf_colors=tree_leaf_colors,
            branch_linewidth=branch_linewidth,
        )
        # Persist tree_100 leaf colors so generate_overview_png can read them back
        colors_100 = tree_leaf_colors.get("100", {}).get("display_to_color")
        if colors_100:
            with open(family_dir / f"{family}_colors_100.json", "w") as f:
                json.dump(colors_100, f)
        save_tree_icon(
            family=family,
            output_dir=family_dir,
            bio_trees=bio_trees,
            icon_size=icon_size,
            icon_bg_color=icon_bg_color,
            icon_branch_color=icon_branch_color,
            branch_linewidth=branch_linewidth,
            leaf_colors=tree_leaf_colors.get("100", {}).get("display_to_color"),
        )

    # -------------------------------------------------------------------------
    # Summary + status
    # -------------------------------------------------------------------------
    seq_lengths_all = [
        len(rec.seq)
        for genomes in species_genomes.values()
        for proteins in genomes.values()
        for rec in proteins.values()
    ]
    seqlen_stats = compute_seqlen_stats(seq_lengths_all)
    summary_row_for_report = build_summary_row(
        family=family,
        family_taxid=family_taxid,
        family_lineage=family_lineage,
        seq_type=seq_type,
        region=region,
        segment=segment,
        n_species_discovered=len(species_list),
        n_species_with_seqs=len(species_genomes),
        seqlen_stats=seqlen_stats,
        tree_stats=tree_stats,
        family_annotation=family_annotation,
    )
    if summary_path is not None:
        write_summary_row(summary_path, summary_row_for_report)
    # ---- Per-family report PDF (parity with single-protein, plus a
    #      per-marker coverage page — CONCAT_DESIGN.md §6.3) ----
    tree_seq_lengths = {
        label: stats.get("_seq_lengths", [])
        for label, stats in tree_stats.items()
        if stats.get("_seq_lengths")
    }
    tree_support = {
        label: stats.get("_support_vals", [])
        for label, stats in tree_stats.items()
        if stats.get("_support_vals")
    }
    marker_coverage_for_report = {
        label: stats.get("_marker_coverage_map", {})
        for label, stats in tree_stats.items()
        if stats.get("_marker_coverage_map")
    }
    try:
        generate_family_report(
            family=family,
            output_pdf=family_dir / f"{family}_report.pdf",
            summary_row=summary_row_for_report,
            seq_lengths=seq_lengths_all,
            tree_seq_lengths=tree_seq_lengths,
            tree_support=tree_support,
            bio_trees=bio_trees,
            tree_leaf_colors=tree_leaf_colors,
            branch_linewidth=branch_linewidth,
            marker_coverage=marker_coverage_for_report,
            concat_min_fraction=float(family_cfg["concatenation"].get("min_fraction", 0.7)),
        )
    except Exception as e:
        log.warning("Per-family report PDF skipped: %s", e)

    if status_path is not None:
        write_status_row(status_path, build_status_row(
            family=family,
            family_taxid=family_taxid,
            family_lineage=family_lineage,
            seq_type=seq_type,
            region=region,
            segment=segment,
            status="OK",
            family_annotation=family_annotation,
        ))
    mark_done(family_dir, family, output_dir)

    leaves_500 = tree_stats.get("500", {}).get("leaves", "?")
    leaves_100 = tree_stats.get("100", {}).get("leaves", "?")
    log.info("=" * 60)
    log.info(
        "Pipeline complete for %s  (concat: %d markers, "
        "tree_500: %s leaves, tree_100: %s leaves)",
        family, len(marker_set), leaves_500, leaves_100,
    )
    log.info("=" * 60)


# ---------------------------------------------------------------------------
# Per-target runner
# ---------------------------------------------------------------------------

def _run_target_concat(
    *,
    label: str,
    target_n: int,
    max_reps: int,
    species_genomes: dict[str, dict[str, dict[str, SeqRecord]]],
    species_lineages: dict[str, list[dict]],
    marker_set: list[dict],
    marker_order: list[str],
    family: str,
    family_dir: Path,
    work_dir: Path,
    threads: int,
    threshold_min: float,
    threshold_max: float,
    clustering_tool: str,
    absorb_enabled: bool,
    absorb_threshold: float,
    hi_mult: float,
    lo_mult: float,
    family_cfg: dict,
) -> dict:
    target_work = work_dir / f"target_{label}"
    target_work.mkdir(parents=True, exist_ok=True)

    # 1. Cluster + merge to pick selected genomes
    log.info(
        "Clustering genomes per species on concatenated sequences "
        "(%s, max_reps=%d, search range=%.2f–%.2f) ...",
        clustering_tool, max_reps, threshold_min, threshold_max,
    )
    selected_ids, sel_stats = cluster_and_merge_genomes(
        species_genomes=species_genomes,
        marker_order=marker_order,
        target_n=target_n,
        max_reps_per_species=max_reps,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        work_dir=target_work / "selection",
        clustering_tool=clustering_tool,
        refseq_absorption_enabled=absorb_enabled,
        refseq_absorption_threshold=absorb_threshold,
    )
    log.info(
        "Clustering complete (tree_%s): %d total representative genome(s) "
        "across %d species",
        label, sel_stats["n_total_reps"], sel_stats["n_species_with_reps"],
    )
    log.info(
        "Proportional merge: selected %d genomes for tree_%s (target=%d)",
        sel_stats["n_selected"], label, target_n,
    )

    if len(selected_ids) < 4:
        log.warning(
            "Only %d genome(s) selected for tree_%s — skipping (minimum 4).",
            len(selected_ids), label,
        )
        return {"leaves": 0, "skipped_reason": f"too few genomes ({len(selected_ids)})"}

    # 2. Build genome-level dicts limited to selected genomes
    #    (also collect per-genome species name + lineage for downstream tree annotation)
    selected_genomes: dict[str, dict[str, SeqRecord]] = {}
    genome_to_species: dict[str, str] = {}
    for sp_name, genomes in species_genomes.items():
        for genome_id, proteins in genomes.items():
            if genome_id in selected_ids:
                selected_genomes[genome_id] = proteins
                genome_to_species[genome_id] = sp_name

    # 3. Per-marker length-outlier (RefSeq exempt)
    refseq_ids = identify_refseq_genomes({"_": selected_genomes})  # bulk shape
    selected_genomes, lo_stats = remove_per_marker_length_outliers(
        selected_genomes, refseq_genome_ids=refseq_ids,
        hi_mult=hi_mult, lo_mult=lo_mult,
    )

    # 4-5. Iterative align + concat + tree, with branch-length outlier
    # removal between iterations (CONCAT_DESIGN.md §5.7 + READMEsec'8).
    # RefSeq genomes are exempt; the loop stops when no removable outliers
    # remain or max_iterations is reached.
    msa_cfg  = family_cfg[f"msa_{label}"]
    trim_cfg = family_cfg.get("msa_trim", {}) or {}
    msa_options  = msa_cfg.get("options_aa", "--auto")  # concat is always protein
    trim_enabled = bool(trim_cfg.get("enabled", True))
    trim_options = trim_cfg.get("options", "-automated1")

    outlier_cfg          = family_cfg.get("outlier_removal", {}) or {}
    do_outlier_removal   = bool(outlier_cfg.get("enabled", True))
    outlier_factor       = float(outlier_cfg.get("factor", 20.0))
    outlier_max_iter     = int(outlier_cfg.get("max_iterations", 3))
    outlier_min_genomes  = int(outlier_cfg.get("min_seqs", 40))

    tree_output_dir = target_work / "tree"
    tree_output_dir.mkdir(parents=True, exist_ok=True)
    tree_cfg  = family_cfg[f"tree_{label}"]
    tree_tool = tree_cfg.get("tool", "fasttree" if label == "500" else "iqtree")
    concat_fasta   = target_work / f"concat_{label}.fasta"
    partition_path = target_work / f"partitions_{label}.nex"

    n_outliers_removed = 0
    treefile: Path | None = None
    concat: dict[str, SeqRecord] = {}
    partition_map: dict[str, tuple[int, int]] = {}

    for iteration in range(outlier_max_iter + 1):
        log.info(
            "tree_%s iteration %d: aligning %d markers with MAFFT for "
            "%d genome(s) (options: %s)%s ...",
            label, iteration, len(marker_order), len(selected_genomes),
            msa_options or "(default)",
            f"; trimAl ({trim_options})" if trim_enabled else "; no trim",
        )
        aligned_per_marker = align_and_trim_markers(
            genomes=selected_genomes,
            marker_order=marker_order,
            work_dir=target_work / "msa",
            msa_tool=msa_cfg.get("tool", "mafft"),
            msa_options=msa_options,
            trim_enabled=trim_enabled,
            trim_options=trim_options,
            threads=threads,
        )

        concat, partition_map = concatenate_aligned_markers(
            aligned_per_marker=aligned_per_marker,
            genome_ids=list(selected_genomes.keys()),
            marker_order=marker_order,
        )
        if not concat:
            log.warning("tree_%s: concat alignment is empty — skipping.", label)
            return {"leaves": 0, "skipped_reason": "empty concat"}

        with open(concat_fasta, "w") as f:
            SeqIO.write(list(concat.values()), f, "fasta")
        write_partition_file_nexus(partition_map, partition_path)
        block_widths_str = ", ".join(
            f"{m}:{partition_map[m][1] - partition_map[m][0] + 1}"
            for m in marker_order if m in partition_map
        )
        total_cols = len(next(iter(concat.values())).seq)
        log.info(
            "Concatenated alignment for tree_%s: %d genomes × %d columns "
            "across %d markers; block widths [%s]",
            label, len(concat), total_cols, len(partition_map), block_widths_str,
        )
        log.info(
            "Wrote NEXUS partition file (%d charsets, %d columns) → %s",
            len(partition_map), total_cols, partition_path,
        )

        if label == "100" and tree_tool.lower().startswith("iqtree"):
            # Partitioned IQ-TREE: -p partitions.nex -m MFP + AA support flag
            tree_options = (
                f"-p {partition_path} "
                + (tree_cfg.get("options_aa", "-B 1000") or "")
            ).strip()
            model_aa = model_nuc = "MFP"
        else:
            tree_options = tree_cfg.get("options", "") or ""
            model_aa  = tree_cfg.get("model_aa", "LG+G")
            model_nuc = tree_cfg.get("model_nuc", "GTR+G")

        treefile = run_tree(
            alignment_fasta=concat_fasta,
            output_dir=tree_output_dir,
            prefix=f"tree_{label}",
            tool=tree_tool,
            seq_type="protein",
            model_nuc=model_nuc,
            model_aa=model_aa,
            options=tree_options,
            threads=threads,
        )
        validate_newick(treefile)

        # Branch-length outlier removal — RefSeq genomes are exempt.
        if not do_outlier_removal or iteration >= outlier_max_iter:
            break
        if len(selected_genomes) <= outlier_min_genomes:
            log.info(
                "tree_%s iteration %d: %d genomes ≤ min_seqs %d — stopping.",
                label, iteration, len(selected_genomes), outlier_min_genomes,
            )
            break

        all_outliers = find_branch_length_outliers(treefile, outlier_factor)
        if not all_outliers:
            log.info("tree_%s iteration %d: no branch-length outliers — stopping.",
                     label, iteration)
            break

        protected_outliers = all_outliers & refseq_ids
        removable_outliers = all_outliers - refseq_ids

        bl_stats = branch_length_stats(treefile)
        median_bl = bl_stats["median"]
        mad       = bl_stats["mad"]
        threshold = median_bl + outlier_factor * mad
        bl_map    = bl_stats["bl_map"]
        log.info(
            "tree_%s iteration %d: branch-length stats — median=%.5f, MAD=%.5f, "
            "threshold (median + %.1f×MAD)=%.5f",
            label, iteration, median_bl, mad, outlier_factor, threshold,
        )
        for oid in protected_outliers:
            bl = bl_map.get(oid, float("nan"))
            mads_above = (bl - median_bl) / mad if mad else float("nan")
            log.warning(
                "tree_%s iteration %d: RefSeq genome %r looks like a "
                "branch-length outlier — branch length %.5f (%.1f× MAD above "
                "median, threshold=%.5f) — KEEPING (RefSeq protected)",
                label, iteration, oid, bl, mads_above, threshold,
            )

        if not removable_outliers:
            log.info(
                "tree_%s iteration %d: only protected (RefSeq) outliers remain — stopping.",
                label, iteration,
            )
            break

        remaining = len(selected_genomes) - len(removable_outliers)
        if remaining < outlier_min_genomes:
            log.info(
                "tree_%s iteration %d: removing %d outlier(s) would leave %d genomes "
                "(min_seqs=%d) — stopping.",
                label, iteration, len(removable_outliers), remaining, outlier_min_genomes,
            )
            break

        for oid in removable_outliers:
            bl = bl_map.get(oid, float("nan"))
            mads_above = (bl - median_bl) / mad if mad else float("nan")
            log.info(
                "tree_%s iteration %d: removing outlier genome %r — branch length "
                "%.5f (%.1f× MAD above median, threshold=%.5f)",
                label, iteration, oid, bl, mads_above, threshold,
            )
        for oid in removable_outliers:
            selected_genomes.pop(oid, None)
        n_outliers_removed += len(removable_outliers)

        if len(selected_genomes) < 4:
            log.warning(
                "tree_%s: only %d genome(s) left after branch-length outlier removal — stopping.",
                label, len(selected_genomes),
            )
            break

    # 6. Build display names (one per genome) — used for visualization / PhyloXML.
    #    Format: <species_safe>|<accession>.  Spaces in species names are
    #    replaced with underscores so Newick display labels stay parseable.
    short_to_display: dict[str, str] = {}
    for gid in selected_genomes:
        sp_name = genome_to_species.get(gid, "")
        sp_safe = sp_name.replace(" ", "_") if sp_name else "unknown"
        short_to_display[gid] = f"{sp_safe}|{gid}"

    # 7. Annotate tree with LCA labels and root.
    metadata_for_annotation = [
        {
            "short_id":       gid,
            "lineage":        [e["name"] for e in species_lineages.get(genome_to_species.get(gid, ""), [])],
            "lineage_ranked": species_lineages.get(genome_to_species.get(gid, ""), []),
        }
        for gid in selected_genomes
    ]
    annotated_nwk = target_work / f"tree_{label}_annotated.nwk"
    bio_tree = annotate_tree(
        newick_path=treefile,
        id_map=short_to_display,
        metadata=metadata_for_annotation,
        output_nwk=annotated_nwk,
        lca_min_rank=family_cfg.get("taxonomy", {}).get("lca_min_rank", "none"),
    )

    # 8. Final Newick: deep copy bio_tree, rename terminals to display names,
    #    drop internal labels (matches single-protein convention).
    output_nwk = family_dir / f"{family}_tree_{label}.nwk"
    if bio_tree is not None:
        nwk_tree = copy.deepcopy(bio_tree)
        for clade in nwk_tree.find_clades():
            if clade.is_terminal():
                clade.name = short_to_display.get(clade.name, clade.name)
            else:
                clade.name = None
        Phylo.write(nwk_tree, str(output_nwk), "newick")
        log.info("Newick written to %s", output_nwk)

    # 9. Genus/subfamily leaf coloring — sel_records is a list of SeqRecord
    #    proxies (only .id is consulted), short_id_to_meta carries lineage.
    sel_records = list(concat.values())
    short_id_to_meta = {
        gid: {"lineage_ranked": species_lineages.get(genome_to_species.get(gid, ""), [])}
        for gid in selected_genomes
    }
    genus_inference = family_cfg.get("coloring", {}).get("genus_inference", "none")
    display_to_color, _short_to_color, genus_to_color, subfamily_to_genera = assign_leaf_colors(
        sel_records, short_id_to_meta, short_to_display, genus_inference=genus_inference,
    )
    n_subfamilies = sum(1 for k in subfamily_to_genera if k != "(unclassified)")

    # 10. Display-name copy of bio_tree for image rendering.  Stored on the
    #     target_stats so the post-loop run_family_concat can call
    #     save_tree_images / save_tree_icon.
    display_tree = None
    if bio_tree is not None:
        display_tree = copy.deepcopy(bio_tree)
        for clade in display_tree.find_clades():
            if clade.is_terminal():
                clade.name = short_to_display.get(clade.name, clade.name)

    # 11. PhyloXML — display names on leaves, source-nuc accession in <accession>.
    xml_path = family_dir / f"{family}_tree_{label}.xml"
    leaf_metadata = {
        gid: {"species":        genome_to_species.get(gid, ""),
              "accession":      gid,
              "seq_name":       short_to_display.get(gid, gid),
              "lineage_ranked": species_lineages.get(genome_to_species.get(gid, ""), [])}
        for gid in selected_genomes
    }
    write_phyloxml(
        newick_path=annotated_nwk,
        output_xml=xml_path,
        id_map=short_to_display,
        leaf_metadata=leaf_metadata,
        family=family,
        tree=bio_tree,
        phylogeny_name=f"{family} [concatenated|{len(marker_order)} markers]",
        phylogeny_detail=f"concatenated {len(marker_order)}-marker phylogeny "
                         f"(target_{label}, partitioned)" if label == "100"
                         else f"concatenated {len(marker_order)}-marker phylogeny (target_{label})",
        confidence_type="SH_aLRT" if label == "100" else "SH_like",
        leaf_colors=display_to_color,
    )

    # 8. Concat-specific stats: marker coverage, concat length, partition models
    n_markers_target = len(marker_set)
    n_markers_used   = sum(
        1 for m in marker_order
        if any(m in selected_genomes[gid] for gid in selected_genomes)
    )
    concat_length = len(next(iter(concat.values())).seq) if concat else 0
    marker_coverage = ",".join(
        f"{m}:{sum(1 for gid in selected_genomes if m in selected_genomes[gid])}"
        for m in marker_order
    )
    marker_models_str = ""
    if label == "100" and tree_tool.lower().startswith("iqtree"):
        # Read per-partition models from <prefix>.best_scheme.nex
        models_by_charset = parse_iqtree_partition_models(
            tree_output_dir / f"tree_{label}.iqtree"
        )
        if models_by_charset:
            # Map sanitized charset names back to canonical marker names
            from .concat import _safe_charset_name
            sanitized_to_marker = {_safe_charset_name(m): m for m in marker_order}
            marker_models_str = ", ".join(
                f"{sanitized_to_marker.get(cs, cs)}:{model}"
                for cs, model in models_by_charset.items()
            )
            log.info(
                "IQ-TREE ModelFinder per-partition (tree_%s): %s",
                label, marker_models_str,
            )
            # summary.tsv stores comma-only (no spaces) for grep-friendliness
            marker_models_str = marker_models_str.replace(", ", ",")
        else:
            log.warning(
                "tree_%s: could not parse per-partition models from "
                "<prefix>.best_scheme.nex — marker_models column will be empty",
                label,
            )

    target_stats: dict[str, Any] = {
        "leaves":               len(selected_genomes),
        "cluster_thresh_min":   sel_stats["cluster_thresh_min"],
        "cluster_thresh_max":   sel_stats["cluster_thresh_max"],
        "seq_type":             "protein",
        "msa_tool":             msa_cfg.get("tool", "mafft"),
        "msa_options":          msa_options,
        "trim_enabled":         trim_enabled,
        "trim_tool":            trim_cfg.get("tool", "trimal") if trim_enabled else "",
        "trim_options":         trim_options if trim_enabled else "",
        "tree_tool":            tree_tool,
        "tree_model":           model_aa if label == "100" else (tree_cfg.get("model_aa") or ""),
        "tree_options":         tree_options,
        "support_type":         "SH_aLRT" if label == "100" else "SH_like",
        "n_outliers_removed":   n_outliers_removed,
        "n_length_outliers_long":  lo_stats["n_long_dropped"],
        "n_length_outliers_short": lo_stats["n_short_dropped"],
        "n_refseq_absorbed":    sel_stats["n_refseq_absorbed"],
        "concat_n_markers_target": n_markers_target,
        "concat_n_markers_used":   n_markers_used,
        "concat_length":           concat_length,
        "marker_coverage":         marker_coverage,
        "marker_models":           marker_models_str,
        "n_genera":                len(genus_to_color),
        "n_subfamilies":           n_subfamilies,
        # Carried for post-loop image rendering (not written into summary.tsv)
        "_display_tree":           display_tree,
        "_display_to_color":       display_to_color,
        "_genus_to_color":         genus_to_color,
        "_subfamily_to_genera":    subfamily_to_genera,
    }
    if bio_tree is not None:
        target_stats["support"] = compute_support_stats(bio_tree)
        target_stats["_support_vals"] = [
            c.confidence for c in bio_tree.find_clades()
            if not c.is_terminal() and c.confidence is not None
        ]
    target_stats["msa"] = compute_msa_stats(concat_fasta)
    # Structured per-marker coverage map for the report PDF (underscored so
    # DictWriter(extrasaction='ignore') won't try to write it into summary.tsv).
    target_stats["_marker_coverage_map"] = {
        m: sum(1 for gid in selected_genomes if m in selected_genomes[gid])
        for m in marker_order
    }
    # Concat length per genome — used for the per-tree length histogram.
    target_stats["_seq_lengths"] = [len(rec.seq) for rec in concat.values()]
    return target_stats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fetch_all_species(
    *,
    family: str,
    species_list: list[dict],
    marker_set: list[dict],
    work_dir: Path,
    max_per_species: int,
    min_fraction: float,
    exclude_organisms: list[str] | None,
) -> tuple[dict[str, dict[str, dict[str, SeqRecord]]], dict[str, list[dict]], dict]:
    """Fetch genomes per species and return:
       - species_genomes: {species_name: {genome_id: {marker_name: SeqRecord}}}
       - species_lineages: {species_name: ranked NCBI lineage}
       - aggregated fetch stats
    """
    # Resolve species lineages once (one Entrez fetch instead of per-species)
    taxids = [str(sp["taxid"]) for sp in species_list]
    lineage_cache = work_dir / "concat_lineages.json"
    if lineage_cache.exists():
        with open(lineage_cache) as f:
            taxid_to_lineage = json.load(f)
        log.info("Reusing cached ranked taxonomy: %d taxa", len(taxid_to_lineage))
    else:
        log.info("Fetching ranked taxonomy from NCBI for %d species ...", len(taxids))
        taxid_to_lineage = fetch_taxonomy_lineages(taxids)
        with open(lineage_cache, "w") as f:
            json.dump(taxid_to_lineage, f)

    species_genomes: dict[str, dict[str, dict[str, SeqRecord]]] = {}
    species_lineages: dict[str, list[dict]] = {}

    agg = {
        "n_proteins_fetched":          0,
        "n_genomes_found":             0,
        "n_genomes_kept":              0,
        "n_dropped_min_fraction":      0,
        "n_dropped_split_submission":  0,
        "n_orphaned_no_source":        0,
    }

    # One-time warning that concat mode bypasses the global sequence cache.
    # Tracked in CONCAT_DESIGN.md as a deferred phase-10 task.
    log.warning(
        "Note: concat mode does not yet use the global sequence cache; "
        "every concat run hits NCBI for every species (CONCAT_DESIGN.md phase 10).",
    )

    for i, sp in enumerate(species_list, start=1):
        sp_name  = sp["name"]
        sp_taxid = sp["taxid"]
        sp_work  = work_dir / "species" / f"sp_{sp_taxid}_concat"
        sp_work.mkdir(parents=True, exist_ok=True)

        sp_lineage = taxid_to_lineage.get(str(sp_taxid), [])
        species_lineages[sp_name] = sp_lineage

        try:
            genomes, stats = fetch_species_genomes(
                taxid=sp_taxid,
                species_name=f"[{i}/{len(species_list)}] {sp_name}",
                species_lineage=sp_lineage,
                marker_set=marker_set,
                output_dir=sp_work,
                max_per_species=max_per_species,
                min_fraction=min_fraction,
                exclude_organisms=exclude_organisms,
            )
        except Exception as e:
            log.warning("[%d/%d] Concat fetch failed for %s: %s — skipping species.",
                        i, len(species_list), sp_name, e)
            continue

        if genomes:
            species_genomes[sp_name] = genomes
        for k in agg:
            agg[k] += stats.get(k, 0)

    return species_genomes, species_lineages, agg


def _emit_skip(
    family: str,
    family_taxid: int | None,
    family_lineage: list[dict],
    family_annotation: dict,
    seq_type: str,
    region: str,
    segment: str | None,
    n_species_discovered: int,
    n_species_with_seqs: int,
    reason: str,
    summary_path: Path | None,
    status_path: Path | None,
    family_dir: Path,
    output_dir: Path,
    mark_skipped,
) -> None:
    if summary_path is not None:
        row = build_summary_row(
            family=family,
            family_taxid=family_taxid,
            family_lineage=family_lineage,
            seq_type=seq_type,
            region=region,
            segment=segment,
            n_species_discovered=n_species_discovered,
            n_species_with_seqs=n_species_with_seqs,
            seqlen_stats={},
            tree_stats={},
            family_annotation=family_annotation,
        )
        write_summary_row(summary_path, row)
    if status_path is not None:
        write_status_row(status_path, build_status_row(
            family=family,
            family_taxid=family_taxid,
            family_lineage=family_lineage,
            seq_type=seq_type,
            region=region,
            segment=segment,
            status=reason,
            family_annotation=family_annotation,
        ))
    mark_skipped(family_dir, family, output_dir, reason)
