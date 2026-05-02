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

import csv
import json
from pathlib import Path
from typing import Any

from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord

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
from .tree import run_tree, validate_newick

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
    log.info("Concatenated mode: %s — %d markers", family, len(marker_set))
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
    log.info(
        "Concat fetch complete: %d genome(s) across %d species (proteins fetched: %d)",
        n_genomes_total, len(species_genomes), fetch_stats["n_proteins_fetched"],
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
    # Summary + status
    # -------------------------------------------------------------------------
    seqlen_stats = compute_seqlen_stats([
        len(rec.seq)
        for genomes in species_genomes.values()
        for proteins in genomes.values()
        for rec in proteins.values()
    ])
    if summary_path is not None:
        row = build_summary_row(
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
        write_summary_row(summary_path, row)
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
    log.info("tree_%s: selected %d genomes (from %d species)",
             label, sel_stats["n_selected"], sel_stats["n_species_with_reps"])

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

    # 4. Per-marker MSA + trim, then concatenate
    msa_cfg  = family_cfg[f"msa_{label}"]
    trim_cfg = family_cfg.get("msa_trim", {}) or {}
    msa_options = msa_cfg.get("options_aa", "--auto")  # concat is always protein
    trim_enabled = bool(trim_cfg.get("enabled", True))
    trim_options = trim_cfg.get("options", "-automated1")

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

    concat_fasta = target_work / f"concat_{label}.fasta"
    with open(concat_fasta, "w") as f:
        SeqIO.write(list(concat.values()), f, "fasta")

    partition_path = target_work / f"partitions_{label}.nex"
    write_partition_file_nexus(partition_map, partition_path)

    # 5. Tree inference
    tree_output_dir = target_work / "tree"
    tree_output_dir.mkdir(parents=True, exist_ok=True)
    tree_cfg = family_cfg[f"tree_{label}"]
    tree_tool = tree_cfg.get("tool", "fasttree" if label == "500" else "iqtree")

    if label == "100" and tree_tool.lower().startswith("iqtree"):
        # Partitioned IQ-TREE: -p partitions.nex -m MFP, plus existing AA support flag
        partition_options = (
            f"-p {partition_path} "
            + (tree_cfg.get("options_aa", "-B 1000") or "")
        ).strip()
        tree_options = partition_options
        model_aa = "MFP"
        model_nuc = "MFP"
    else:
        tree_options = tree_cfg.get("options", "") or ""
        model_aa = tree_cfg.get("model_aa", "LG+G")
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

    # 6. Annotate tree with LCA labels and root.
    #    id_map = {genome_id: genome_id} (concat leaves have human-readable IDs already).
    id_map = {gid: gid for gid in selected_genomes}
    metadata_for_annotation = [
        {
            "short_id":       gid,
            "lineage":        [e["name"] for e in species_lineages.get(genome_to_species.get(gid, ""), [])],
            "lineage_ranked": species_lineages.get(genome_to_species.get(gid, ""), []),
        }
        for gid in selected_genomes
    ]
    annotated_nwk = family_dir / f"{family}_tree_{label}.nwk"
    bio_tree = annotate_tree(
        newick_path=treefile,
        id_map=id_map,
        metadata=metadata_for_annotation,
        output_nwk=annotated_nwk,
        lca_min_rank=family_cfg.get("taxonomy", {}).get("lca_min_rank", "none"),
    )

    # 7. PhyloXML
    xml_path = family_dir / f"{family}_tree_{label}.xml"
    leaf_metadata = {
        gid: {"species": genome_to_species.get(gid, ""),
              "accession": gid,
              "seq_name": gid,
              "lineage_ranked": species_lineages.get(genome_to_species.get(gid, ""), [])}
        for gid in selected_genomes
    }
    write_phyloxml(
        newick_path=annotated_nwk,
        output_xml=xml_path,
        id_map=id_map,
        leaf_metadata=leaf_metadata,
        family=family,
        tree=bio_tree,
        phylogeny_name=f"{family} [concatenated|{len(marker_order)} markers]",
        phylogeny_detail=f"concatenated {len(marker_order)}-marker phylogeny "
                         f"(target_{label}, partitioned)" if label == "100"
                         else f"concatenated {len(marker_order)}-marker phylogeny (target_{label})",
        confidence_type="SH_aLRT" if label == "100" else "SH_like",
    )

    # 8. Stats
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
        "n_outliers_removed":   0,                                # Phase 7
        "n_length_outliers_long":  lo_stats["n_long_dropped"],
        "n_length_outliers_short": lo_stats["n_short_dropped"],
        "n_refseq_absorbed":    sel_stats["n_refseq_absorbed"],
        "n_genera":             "",
        "n_subfamilies":        "",
    }
    if bio_tree is not None:
        target_stats["support"] = compute_support_stats(bio_tree)
    target_stats["msa"] = compute_msa_stats(concat_fasta)
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
    else:
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

    for i, sp in enumerate(species_list, start=1):
        sp_name  = sp["name"]
        sp_taxid = sp["taxid"]
        sp_work  = work_dir / "species" / f"sp_{sp_taxid}_concat"
        sp_work.mkdir(parents=True, exist_ok=True)

        sp_lineage = taxid_to_lineage.get(str(sp_taxid), [])
        species_lineages[sp_name] = sp_lineage

        log.info("[%d/%d] Concat fetch: %s", i, len(species_list), sp_name)
        try:
            genomes, stats = fetch_species_genomes(
                taxid=sp_taxid,
                species_name=sp_name,
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
