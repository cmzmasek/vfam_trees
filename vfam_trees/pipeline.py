"""Per-family pipeline runner (called by Snakemake rules)."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import shutil
import statistics
from pathlib import Path

import yaml
from Bio import Phylo
from Bio.SeqRecord import SeqRecord

from .config import load_family_config, load_global_config
from .fetch import (
    configure_entrez, discover_species, get_family_taxid,
    fetch_species_sequences, parse_gb_records, extract_metadata,
    fetch_taxonomy_lineages,
)
from .summary import (
    compute_support_stats, compute_msa_stats, compute_seqlen_stats,
    build_summary_row, write_summary_row, load_family_annotations,
    build_status_row, write_status_row,
)
from .quality import filter_sequences, remove_length_outliers, deduplicate, write_fasta
from .rename import assign_short_ids, restore_fasta_names, restore_newick_names
from .subsample import absorb_into_refseqs, adaptive_cluster_species, proportional_merge
from .msa import run_msa, get_mafft_version, validate_msa
from .trim import run_trim, get_trimal_version
from .tree import (
    run_tree, get_tree_tool_version, validate_newick,
    parse_iqtree_best_model, is_model_finder_spec,
)
from .taxonomy import annotate_tree
from .phyloxml_writer import write_phyloxml
from .report import generate_family_report, save_tree_images, save_tree_icon, save_sequence_length_plot
from .colors import assign_leaf_colors
from .cache import SequenceCache
from .logger import setup_logger, get_logger

log = get_logger(__name__)


def _is_refseq_accession(accession: str) -> bool:
    """Return True if *accession* follows the NCBI RefSeq ``XX_`` prefix format.

    Examples: NC_045512.2, NZ_CP012345.1, YP_009724390.1 → True;
    MN908947.3, OP123456 → False.
    """
    acc = (accession or "").strip()
    return (
        len(acc) >= 3
        and acc[2] == "_"
        and acc[0:2].isalpha()
        and acc[0:2].isupper()
    )


def _compute_key(obj) -> str:
    """Return a short hex digest summarising *obj* (sort_keys for stability)."""
    payload = json.dumps(obj, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _resolve_tree_options(tree_cfg: dict, seq_type: str) -> str:
    """Return the tree-tool options string for *seq_type*.

    Prefers the per-sequence-type keys (``options_nuc`` / ``options_aa``)
    introduced in 1.0.15 and falls back to the legacy single ``options``
    key so pre-existing family YAMLs continue to work unchanged.
    """
    key = "options_aa" if seq_type == "protein" else "options_nuc"
    if key in tree_cfg:
        return tree_cfg.get(key, "") or ""
    return tree_cfg.get("options", "") or ""


def _support_type_for(tool: str, options: str) -> str:
    """Return a short label naming the branch-support measure produced by
    *tool* under *options* — one of 'SH_like', 'SH_aLRT', 'UFBoot',
    'SH_aLRT_UFBoot'.

    The IQ-TREE wrapper auto-injects ``-alrt 1000`` whenever no bootstrap
    flag is present, so plain or ``--fast`` runs yield SH-aLRT.
    """
    tool_norm = tool.lower().replace("-", "").replace("_", "")
    if tool_norm == "fasttree":
        return "SH_like"
    parts = options.split() if options else []
    has_ufboot = "-B" in parts or "--ufboot" in parts or "-bb" in parts
    has_alrt = "-alrt" in parts
    has_bootstrap = "-b" in parts
    if has_ufboot and has_alrt:
        return "SH_aLRT_UFBoot"
    if has_ufboot:
        return "UFBoot"
    if has_bootstrap:
        return "Bootstrap"
    return "SH_aLRT"


def _write_checkpoint(path: Path, key_obj: dict) -> None:
    """Write a checkpoint sidecar encoding the input key at completion time."""
    path.write_text(json.dumps({"key": _compute_key(key_obj)}))


def _check_checkpoint(path: Path, key_obj: dict) -> bool:
    """Return True when the checkpoint exists and its stored key matches."""
    if not path.exists():
        return False
    try:
        stored = json.loads(path.read_text())
    except Exception:
        return False
    return stored.get("key") == _compute_key(key_obj)


def run_family(
    family: str,
    global_config_path: Path,
    configs_dir: Path,
    output_dir: Path,
    threads: int = 1,
    log_level: str = "INFO",
    summary_path: Path | None = None,
    status_path: Path | None = None,
) -> None:
    """Execute the full pipeline for a single viral family."""
    global_cfg = load_global_config(global_config_path)
    configure_entrez(
        email=global_cfg["ncbi"]["email"],
        api_key=global_cfg["ncbi"].get("api_key") or None,
    )

    # External family-level annotations (Baltimore class, etc.) for the summary.
    # Path is taken from global config; default looks for the file next to it.
    ann_cfg = global_cfg.get("annotation_tsv")
    if ann_cfg:
        ann_path = Path(ann_cfg).expanduser()
        if not ann_path.is_absolute():
            ann_path = global_config_path.parent / ann_path
    else:
        ann_path = global_config_path.parent / "virus_families_annotation.tsv"
    family_annotation = load_family_annotations(ann_path).get(family.lower(), {})

    # Fetch family-level taxonomy (taxid + ranked lineage) for the summary.
    # Done first so the taxid can be embedded in the output directory name.
    family_taxid = get_family_taxid(family)
    if family_taxid is not None:
        family_lineage_map = fetch_taxonomy_lineages([family_taxid])
        family_lineage = family_lineage_map.get(str(family_taxid), [])
    else:
        family_lineage = []

    dir_name = f"{family}_{family_taxid}" if family_taxid is not None else family
    family_dir = output_dir / dir_name
    family_dir.mkdir(parents=True, exist_ok=True)

    log_file = family_dir / f"{family}.log"
    setup_logger(f"vfam_trees.{family}", log_file=log_file, level=log_level)
    log = get_logger(f"vfam_trees.{family}")

    log.info("=" * 60)
    log.info("Starting pipeline for %s", family)
    log.info("=" * 60)

    if family_taxid is None:
        log.warning("Could not resolve NCBI taxid for family %s", family)

    # Global sequence cache (optional — only active when cache.dir is set in global.yaml)
    cache_cfg  = global_cfg.get("cache") or {}
    cache_dir  = cache_cfg.get("dir") or None
    seq_cache: SequenceCache | None = None
    if cache_dir:
        ttl_days  = cache_cfg.get("ttl_days") or None
        seq_cache = SequenceCache(Path(cache_dir), ttl_days=ttl_days)
        cs = seq_cache.stats()
        log.info(
            "Global sequence cache: %s  (%d entries, %.1f MB)",
            seq_cache.cache_dir, cs["entries"], cs["size_mb"],
        )

    family_cfg, auto_generated = load_family_config(family, configs_dir, global_cfg)
    if auto_generated:
        log.warning("Using auto-generated config for %s — consider reviewing it.", family)

    _save_config_copy(family_cfg, family_dir / f"{family}.yaml")

    work_dir = family_dir / "_work"
    work_dir.mkdir(exist_ok=True)

    seq_type = family_cfg["sequence"]["type"]
    region = family_cfg["sequence"]["region"]
    segment = family_cfg["sequence"].get("segment") or None
    if segment:
        log.info("Segmented family detected — using segment: %s", segment)
    max_per_species = family_cfg["download"]["max_per_species"]
    clustering_tool = family_cfg["clustering"].get("tool", "mmseqs2")
    quality_cfg = family_cfg["quality"]

    # -------------------------------------------------------------------------
    # Step 1: Discover species
    # -------------------------------------------------------------------------
    species_cache = work_dir / "species_list.json"
    if species_cache.exists():
        with open(species_cache) as f:
            species_list = json.load(f)
        log.info("Reusing cached species list: %d species", len(species_list))
    else:
        species_list = discover_species(family)
        if not species_list:
            log.warning("No species found for %s — skipping.", family)
            _mark_skipped(family_dir, family, output_dir, "no species found in NCBI taxonomy")
            if summary_path is not None:
                row = build_summary_row(
                    family=family,
                    family_taxid=family_taxid,
                    family_lineage=family_lineage,
                    seq_type=seq_type,
                    region=region,
                    segment=segment,
                    n_species_discovered=0,
                    n_species_with_seqs=0,
                    seqlen_stats={},
                    tree_stats={},
                    n_species_relaxed=0,
                    total_seqs_qc=0,
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
                    status="no species found in NCBI taxonomy",
                    family_annotation=family_annotation,
                ))
            return
        with open(species_cache, "w") as f:
            json.dump(species_list, f, indent=2)
        log.info("Discovered %d species in %s", len(species_list), family)

    # -------------------------------------------------------------------------
    # Step 2: Per-species download + quality filter
    #         Results cached per species to avoid re-downloading
    # -------------------------------------------------------------------------
    species_data: dict[str, dict] = {}  # species_name → {records, metadata}
    species_pre_length_lengths: dict[str, list[int]] = {}  # species_name → lengths after ambiguity filter
    n_species_relaxed = 0  # species that needed a relaxed min_length threshold
    # Aggregate QC exclusion counts across all species
    family_qc_stats: dict[str, int] = {
        "n_excluded_organism": 0,
        "n_excluded_length": 0,
        "n_excluded_ambiguity": 0,
        "n_undefined": 0,
    }

    log.info("Downloading sequences for %d species (max %d per species) ...",
             len(species_list), max_per_species)
    from .summary import format_molecule_region
    log.info("Sequence type: %s", format_molecule_region(seq_type, region, segment))

    for i, sp in enumerate(species_list, start=1):
        sp_name = sp["name"]
        sp_taxid = sp["taxid"]
        sp_key = f"sp_{sp_taxid}"
        sp_work = work_dir / "species" / sp_key
        sp_work.mkdir(parents=True, exist_ok=True)

        gb_file = sp_work / "sequences.gb"
        db      = "nuccore" if seq_type == "nucleotide" else "protein"

        # --- Download with three-tier caching ---
        # Tier 1: local per-run cache (sp_work/sequences.gb)
        if gb_file.exists():
            log.info("[%d/%d] Local cache hit: %s", i, len(species_list), sp_name)

        # Tier 2: global shared cache
        elif seq_cache is not None:
            # Fast path: skip lock acquisition entirely for known-empty species.
            if seq_cache.get_empty(sp_taxid, db, region, segment, max_per_species):
                log.info("[%d/%d] Cached empty result: %s — skipping.", i, len(species_list), sp_name)
                continue

            with seq_cache.download_lock(sp_taxid, db, region, segment, max_per_species) as got_lock:
                # Re-check after acquiring the lock — a concurrent job may have
                # downloaded this species while we were waiting.
                cached = seq_cache.get(sp_taxid, db, region, segment, max_per_species)
                if cached is not None:
                    log.info("[%d/%d] Global cache hit: %s", i, len(species_list), sp_name)
                    shutil.copy2(cached, gb_file)
                elif seq_cache.get_empty(sp_taxid, db, region, segment, max_per_species):
                    log.info("[%d/%d] Cached empty result (concurrent): %s — skipping.",
                             i, len(species_list), sp_name)
                    continue
                else:
                    if not got_lock:
                        log.warning(
                            "[%d/%d] Proceeding without lock for %s", i, len(species_list), sp_name
                        )
                    log.info("[%d/%d] Downloading (will cache): %s", i, len(species_list), sp_name)
                    n = fetch_species_sequences(
                        taxid=sp_taxid,
                        species_name=sp_name,
                        seq_type=seq_type,
                        region=region,
                        output_gb=gb_file,
                        max_per_species=max_per_species,
                        exclude_organisms=quality_cfg.get("exclude_organisms"),
                        segment=segment,
                    )
                    if n == 0:
                        log.info("[%d/%d] No sequences found for %s — caching empty result.",
                                 i, len(species_list), sp_name)
                        seq_cache.store_empty(sp_taxid, db, region, segment, max_per_species, family=family)
                        continue
                    seq_cache.store(sp_taxid, db, region, segment, max_per_species, gb_file, n, family=family)

        # Tier 3: no global cache — plain download
        else:
            log.info("[%d/%d] Downloading: %s", i, len(species_list), sp_name)
            n = fetch_species_sequences(
                taxid=sp_taxid,
                species_name=sp_name,
                seq_type=seq_type,
                region=region,
                output_gb=gb_file,
                max_per_species=max_per_species,
                exclude_organisms=quality_cfg.get("exclude_organisms"),
                segment=segment,
            )
            if n == 0:
                log.info("[%d/%d] No sequences found for %s — skipping.",
                         i, len(species_list), sp_name)
                continue

        # Parse records
        raw_records: list[SeqRecord] = []
        raw_meta: list[dict] = []
        for rec in parse_gb_records(gb_file):
            raw_records.append(rec)
            raw_meta.append(extract_metadata(rec))

        if not raw_records:
            continue

        # Segment keyword validation for segmented families
        if segment:
            before = len(raw_records)
            raw_records = [r for r in raw_records if segment.lower() in r.description.lower()]
            n_missing_seg = before - len(raw_records)
            if n_missing_seg:
                log.warning(
                    "%d/%d record(s) for %s do not contain segment keyword '%s' "
                    "in title — excluding",
                    n_missing_seg, before, sp_name, segment,
                )
            if not raw_records:
                log.info("No records with segment '%s' for %s — skipping.", segment, sp_name)
                continue
            # Rebuild raw_meta to match filtered raw_records
            valid_ids = {r.id for r in raw_records}
            raw_meta = [m for m in raw_meta if m["accession"] in valid_ids]

        acc_to_raw_meta = {m["accession"]: m for m in raw_meta}
        raw_records = deduplicate(raw_records)
        raw_meta = [acc_to_raw_meta[rec.id] for rec in raw_records if rec.id in acc_to_raw_meta]
        raw_records = [rec for rec in raw_records if rec.id in acc_to_raw_meta]

        filtered, fraction_used, sp_qc_stats = filter_sequences(
            raw_records,
            seq_type=seq_type,
            min_length=quality_cfg["min_length"],
            max_ambiguous=quality_cfg["max_ambiguous"],
            exclude_organisms=quality_cfg.get("exclude_organisms"),
        )
        # Accumulate QC stats
        for k in family_qc_stats:
            family_qc_stats[k] += sp_qc_stats.get(k, 0)
        pre_lengths = sp_qc_stats.get("pre_length_lengths", [])
        if pre_lengths:
            species_pre_length_lengths[sp_name] = pre_lengths

        if not filtered:
            log.info("No sequences passed quality filter for %s", sp_name)
            continue
        if fraction_used is not None and fraction_used < 0.5:
            n_species_relaxed += 1

        acc_to_meta = {m["accession"]: m for m in raw_meta}
        filtered_meta = [acc_to_meta[r.id] for r in filtered if r.id in acc_to_meta]
        filtered = [r for r in filtered if r.id in acc_to_meta]

        species_data[sp_name] = {"records": filtered, "metadata": filtered_meta}

    # Sequence length violin plot (per species, after ambiguity filter, before length filter)
    if species_pre_length_lengths:
        save_sequence_length_plot(family, family_dir, species_pre_length_lengths)

    # Fetch ranked taxonomy lineages from NCBI (cached per family)
    taxid_cache = work_dir / "taxonomy_cache.json"
    if taxid_cache.exists():
        with open(taxid_cache) as f:
            taxid_to_lineage = json.load(f)
        log.info("Reusing cached ranked taxonomy: %d taxa", len(taxid_to_lineage))
    else:
        unique_taxids = {
            str(m["taxon_id"])
            for d in species_data.values()
            for m in d["metadata"]
            if m.get("taxon_id")
        }
        log.info("Fetching ranked taxonomy from NCBI for %d unique taxa ...", len(unique_taxids))
        taxid_to_lineage = fetch_taxonomy_lineages(unique_taxids)
        with open(taxid_cache, "w") as f:
            json.dump(taxid_to_lineage, f)

    # Attach ranked lineage to each metadata dict
    for data in species_data.values():
        for meta in data["metadata"]:
            ranked = taxid_to_lineage.get(str(meta.get("taxon_id", "")), [])
            meta["lineage_ranked"] = ranked
            if ranked:
                meta["lineage"] = [e["name"] for e in ranked]

    seq_lengths_all = [len(r.seq) for d in species_data.values() for r in d["records"]]
    seqlen_stats = compute_seqlen_stats(seq_lengths_all)
    total_seqs = sum(len(d["records"]) for d in species_data.values())

    if total_seqs < 4:
        log.warning(
            "Only %d sequence(s) passed quality filters for %s (minimum 4) — skipping.",
            total_seqs, family,
        )
        _mark_skipped(family_dir, family, output_dir, f"too few sequences after QC ({total_seqs} < 4)")
        if summary_path is not None:
            row = build_summary_row(
                family=family,
                family_taxid=family_taxid,
                family_lineage=family_lineage,
                seq_type=seq_type,
                region=region,
                segment=segment,
                n_species_discovered=len(species_list),
                n_species_with_seqs=len(species_data),
                seqlen_stats=seqlen_stats,
                tree_stats={},
                n_species_relaxed=n_species_relaxed,
                total_seqs_qc=total_seqs,
                qc_stats=family_qc_stats,
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
                status=f"too few sequences after QC ({total_seqs} < 4)",
                family_annotation=family_annotation,
            ))
        return
    log.info(
        "Download complete: %d sequences across %d / %d species passed quality filters",
        total_seqs, len(species_data), len(species_list),
    )

    # -------------------------------------------------------------------------
    # Step 3: Assign short IDs over the full set (once)
    # -------------------------------------------------------------------------
    all_records: list[SeqRecord] = []
    all_metadata: list[dict] = []
    for sp_name, data in species_data.items():
        all_records.extend(data["records"])
        all_metadata.extend(data["metadata"])

    id_map_path = family_dir / f"{family}_id_map.tsv"
    renamed_records, short_to_display = assign_short_ids(
        all_records, all_metadata, family, id_map_path
    )

    # Attach short_id back into metadata and rebuild species_data with short IDs
    for meta, renamed_rec in zip(all_metadata, renamed_records):
        meta["short_id"] = renamed_rec.id

    short_id_to_meta: dict[str, dict] = {
        rec.id: meta for rec, meta in zip(renamed_records, all_metadata)
    }

    # Identify RefSeq records so proportional_merge can prefer them when
    # subsampling clusters within each species.
    refseq_short_ids = {
        short_id for short_id, meta in short_id_to_meta.items()
        if _is_refseq_accession(meta.get("accession", ""))
    }
    if refseq_short_ids:
        log.info("Detected %d RefSeq record(s) — will prefer them during subsampling",
                 len(refseq_short_ids))

    # Rebuild species_data with renamed records
    renamed_iter = iter(renamed_records)
    for sp_name, data in species_data.items():
        n = len(data["records"])
        data["renamed"] = [next(renamed_iter) for _ in range(n)]

    # -------------------------------------------------------------------------
    # Steps 4–9: Run for each target size
    # -------------------------------------------------------------------------
    clustering_cfg = family_cfg["clustering"]
    targets = [
        (family_cfg["targets"]["max_500"], clustering_cfg["max_reps_500"], "500"),
        (family_cfg["targets"]["max_100"], clustering_cfg["max_reps_100"], "100"),
    ]

    tree_stats: dict[str, dict] = {}
    tree_support: dict[str, list[float]] = {}
    bio_trees: dict[str, object] = {}
    tree_seq_lengths: dict[str, list[int]] = {}
    tree_leaf_colors: dict[str, dict] = {}  # label → {display_to_color, genus_to_color, subfamily_to_genera}

    for target_n, max_reps, label in targets:
        log.info("-" * 40)
        log.info("Building tree_%s  (target=%d sequences, max_reps_per_species=%d)",
                 label, target_n, max_reps)
        stats, support_vals, bio_tree = _run_target(
            label=label,
            target_n=target_n,
            max_reps=max_reps,
            species_data=species_data,
            short_id_to_meta=short_id_to_meta,
            short_to_display=short_to_display,
            refseq_short_ids=refseq_short_ids,
            seq_type=seq_type,
            clustering_tool=clustering_tool,
            clustering_cfg=clustering_cfg,
            family_cfg=family_cfg,
            family=family,
            family_dir=family_dir,
            work_dir=work_dir,
            threads=threads,
            log=log,
        )
        tree_stats[label] = stats
        if stats.get("seq_lengths"):
            tree_seq_lengths[label] = stats["seq_lengths"]
        if stats.get("display_to_color"):
            tree_leaf_colors[label] = {
                "display_to_color":    stats["display_to_color"],
                "genus_to_color":      stats["genus_to_color"],
                "subfamily_to_genera": stats["subfamily_to_genera"],
            }
        if support_vals:
            tree_support[label] = support_vals
        if bio_tree is not None:
            # Store a display-name copy for the PDF report so leaf labels are
            # human-readable rather than short IDs.
            display_tree = copy.deepcopy(bio_tree)
            for clade in display_tree.find_clades():
                if clade.is_terminal():
                    clade.name = short_to_display.get(clade.name, clade.name)
            bio_trees[label] = display_tree

    log.info("=" * 60)
    log.info("Pipeline complete for %s", family)
    log.info("=" * 60)
    _mark_done(family_dir, family, output_dir)

    summary_row = build_summary_row(
        family=family,
        family_taxid=family_taxid,
        family_lineage=family_lineage,
        seq_type=seq_type,
        region=region,
        segment=segment,
        n_species_discovered=len(species_list),
        n_species_with_seqs=len(species_data),
        seqlen_stats=seqlen_stats,
        tree_stats=tree_stats,
        n_species_relaxed=n_species_relaxed,
        total_seqs_qc=total_seqs,
        qc_stats=family_qc_stats,
        family_annotation=family_annotation,
    )
    if summary_path is not None:
        write_summary_row(summary_path, summary_row)
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

    viz_cfg = global_cfg.get("visualization") or {}
    branch_linewidth = float(viz_cfg.get("branch_linewidth", 0.5))

    icon_cfg = global_cfg.get("icon") or {}
    icon_size = int(icon_cfg.get("size", 256))
    icon_bg_color = icon_cfg.get("bg_color", "#EAF3F2")
    icon_branch_color = icon_cfg.get("branch_color", "#000000")

    # PDF report
    generate_family_report(
        family=family,
        output_pdf=family_dir / f"{family}_report.pdf",
        summary_row=summary_row,
        seq_lengths=seq_lengths_all,
        tree_seq_lengths=tree_seq_lengths,
        tree_support=tree_support,
        bio_trees=bio_trees,
        tree_leaf_colors=tree_leaf_colors,
        branch_linewidth=branch_linewidth,
    )

    # Standalone tree images (PDF + PNG)
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
        colors_path = family_dir / f"{family}_colors_100.json"
        with open(colors_path, "w") as _f:
            json.dump(colors_100, _f)

    # Icon PNG (tree_100 topology only, square, no labels)
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


def _run_target(
    label: str,
    target_n: int,
    max_reps: int,
    species_data: dict,
    short_id_to_meta: dict[str, dict],
    short_to_display: dict[str, str],
    refseq_short_ids: set[str],
    seq_type: str,
    clustering_tool: str,
    clustering_cfg: dict,
    family_cfg: dict,
    family: str,
    family_dir: Path,
    work_dir: Path,
    threads: int,
    log,
) -> tuple[dict, list[float], object]:
    """Run clustering, MSA, tree inference, and annotation for one target size.

    Returns (target_stats, support_values, bio_tree).
    """
    target_work = work_dir / label
    target_work.mkdir(exist_ok=True)

    threshold_min = clustering_cfg.get("threshold_min", 0.70)
    threshold_max = clustering_cfg.get("threshold_max", 0.99)

    # RefSeq absorption (per-species, before adaptive clustering): drops
    # non-RefSeqs that are near-identical (≥ threshold) to a RefSeq.
    absorb_cfg = family_cfg.get("refseq_absorption", {})
    absorb_enabled = absorb_cfg.get("enabled", True)
    absorb_threshold = float(absorb_cfg.get("threshold", 0.99))
    n_refseq_absorbed = 0

    # Adaptive per-species clustering — track actual thresholds used
    log.info("Clustering sequences per species (%s, max_reps=%d, search range=%.2f–%.2f) ...",
             clustering_tool, max_reps, threshold_min, threshold_max)
    species_reps: dict[str, list[SeqRecord]] = {}
    thresholds_used: list[float] = []
    for sp_name, data in species_data.items():
        sp_safe = sp_name.replace(" ", "_").replace("/", "_")
        sp_work = target_work / "clustering" / sp_safe
        sp_records = data["renamed"]
        if absorb_enabled:
            sp_records, n_absorbed = absorb_into_refseqs(
                records=sp_records,
                refseq_ids=refseq_short_ids,
                threshold=absorb_threshold,
                seq_type=seq_type,
                work_dir=sp_work / "absorb",
                clustering_tool=clustering_tool,
            )
            if n_absorbed:
                log.info(
                    "RefSeq absorption: %s — %d non-RefSeq sequence(s) "
                    "absorbed at identity ≥ %.2f",
                    sp_name, n_absorbed, absorb_threshold,
                )
            n_refseq_absorbed += n_absorbed
        reps, threshold_used = adaptive_cluster_species(
            records=sp_records,
            max_reps=max_reps,
            threshold_min=threshold_min,
            threshold_max=threshold_max,
            seq_type=seq_type,
            work_dir=sp_work,
            clustering_tool=clustering_tool,
        )
        species_reps[sp_name] = reps
        thresholds_used.append(threshold_used)
    if absorb_enabled and n_refseq_absorbed:
        log.info(
            "RefSeq absorption total: %d non-RefSeq sequence(s) replaced "
            "by their near-identical RefSeq across %d species",
            n_refseq_absorbed, len(species_data),
        )

    cluster_thresh_min = round(min(thresholds_used), 4) if thresholds_used else ""
    cluster_thresh_max = round(max(thresholds_used), 4) if thresholds_used else ""

    total_reps = sum(len(r) for r in species_reps.values())
    log.info("Clustering complete: %d total representatives across %d species",
             total_reps, len(species_reps))

    # Proportional merge (prefers RefSeqs when subsampling within a species)
    sel_records = proportional_merge(
        species_reps, target_n, priority_ids=refseq_short_ids,
    )
    log.info("Proportional merge: selected %d sequences for tree_%s", len(sel_records), label)

    if len(sel_records) < 4:
        log.warning(
            "Only %d sequence(s) after proportional merge for tree_%s (minimum 4) — skipping.",
            len(sel_records), label,
        )
        return {}, [], None

    # Pre-MSA length-outlier removal (two-sided, config-driven)
    length_outlier_cfg = family_cfg.get("length_outlier", {})
    n_length_long = 0
    n_length_short = 0
    if length_outlier_cfg.get("enabled", True):
        hi_mult = float(length_outlier_cfg.get("hi_mult", 3.0))
        lo_mult = float(length_outlier_cfg.get("lo_mult", 0.333))
        sel_records, n_length_long, n_length_short = remove_length_outliers(
            sel_records, hi_mult=hi_mult, lo_mult=lo_mult,
            protected_ids=refseq_short_ids,
        )
        if n_length_long or n_length_short:
            log.info(
                "tree_%s: length-outlier removal dropped %d long + %d short "
                "(hi_mult=%.2f, lo_mult=%.2f)",
                label, n_length_long, n_length_short, hi_mult, lo_mult,
            )
    if len(sel_records) < 4:
        log.warning(
            "Only %d sequence(s) after length-outlier removal for tree_%s (minimum 4) — skipping.",
            len(sel_records), label,
        )
        return {}, [], None

    sel_metadata = [short_id_to_meta.get(r.id, {}) for r in sel_records]
    id_map = {r.id: short_to_display.get(r.id, r.id) for r in sel_records}

    # Write raw quality-controlled FASTA (short IDs → restored names)
    raw_short_fasta = target_work / f"sequences_raw_{label}_short.fasta"
    write_fasta(sel_records, raw_short_fasta)
    restore_fasta_names(raw_short_fasta, family_dir / f"{family}_sequences_raw_{label}.fasta", id_map)

    # Metadata TSV
    _write_metadata_tsv(sel_metadata, short_to_display, family_dir / f"{family}_metadata_{label}.tsv")

    msa_cfg = family_cfg[f"msa_{label}"]
    msa_options = (
        msa_cfg.get("options_aa", "--auto")
        if seq_type == "protein"
        else msa_cfg.get("options_nuc", "")
    )
    tree_cfg = family_cfg[f"tree_{label}"]
    tree_options = _resolve_tree_options(tree_cfg, seq_type)
    support_type = _support_type_for(tree_cfg["tool"], tree_options)
    trim_cfg = family_cfg.get("msa_trim", {})
    trim_enabled = bool(trim_cfg.get("enabled", True))
    trim_tool = trim_cfg.get("tool", "trimal")
    trim_options = trim_cfg.get("options", "-automated1")
    outlier_cfg = family_cfg.get("outlier_removal", {})
    do_outlier_removal = outlier_cfg.get("enabled", True)
    outlier_factor = float(outlier_cfg.get("factor", 20.0))
    max_iterations = int(outlier_cfg.get("max_iterations", 3))
    outlier_min_seqs = int(outlier_cfg.get("min_seqs", 40))

    msa_short_fasta = target_work / f"alignment_{label}_short.fasta"
    msa_trimmed_fasta = target_work / f"alignment_{label}_short.trim.fasta"
    # File fed to the tree step — trimmed when enabled, raw MAFFT output otherwise.
    tree_input_fasta = msa_trimmed_fasta if trim_enabled else msa_short_fasta
    msa_checkpoint = target_work / ".msa_done"
    tree_output_dir = target_work / "tree"
    expected_treefile = tree_output_dir / f"tree_{label}.nwk"
    tree_checkpoint = target_work / ".tree_done"
    treefile: Path | None = None
    n_outliers_removed = 0

    for iteration in range(max_iterations + 1):
        # ------------------------------------------------------------------
        # MSA (+ optional column trimming) — keyed on sequence set,
        # MAFFT options, and trim settings, so any change invalidates the cache.
        # ------------------------------------------------------------------
        msa_key = {
            "seq_ids":      sorted(r.id for r in sel_records),
            "msa_tool":     msa_cfg["tool"],
            "msa_options":  msa_options,
            "trim_enabled": trim_enabled,
            "trim_tool":    trim_tool if trim_enabled else None,
            "trim_options": trim_options if trim_enabled else None,
        }
        _needed_files = [msa_short_fasta, tree_input_fasta] if trim_enabled else [msa_short_fasta]
        if _check_checkpoint(msa_checkpoint, msa_key) and all(p.exists() for p in _needed_files):
            log.info("Reusing cached MSA%s for tree_%s (iteration %d)",
                     " + trim" if trim_enabled else "", label, iteration)
        else:
            log.info("Running MAFFT (%s) on %d sequences for tree_%s (iteration %d) ...",
                     msa_options, len(sel_records), label, iteration)
            write_fasta(sel_records, raw_short_fasta)
            run_msa(
                input_fasta=raw_short_fasta,
                output_fasta=msa_short_fasta,
                tool=msa_cfg["tool"],
                options=msa_options,
                threads=threads,
            )
            validate_msa(msa_short_fasta)
            log.info("MAFFT complete for tree_%s (iteration %d)", label, iteration)

            if trim_enabled:
                log.info("Running %s (%s) on tree_%s alignment (iteration %d) ...",
                         trim_tool, trim_options, label, iteration)
                run_trim(
                    input_fasta=msa_short_fasta,
                    output_fasta=msa_trimmed_fasta,
                    tool=trim_tool,
                    options=trim_options,
                )
                validate_msa(msa_trimmed_fasta)
                pre_len = compute_msa_stats(msa_short_fasta).get("length", 0) or 0
                post_len = compute_msa_stats(msa_trimmed_fasta).get("length", 0) or 0
                frac = (post_len / pre_len) if pre_len else 0.0
                log.info(
                    "Trim reduced tree_%s alignment from %d → %d columns (%.1f%% kept)",
                    label, pre_len, post_len, 100.0 * frac,
                )
            _write_checkpoint(msa_checkpoint, msa_key)

        # ------------------------------------------------------------------
        # Tree inference — keyed on tree-input MSA content + tool/model/options
        # (trimmed alignment when msa_trim is enabled)
        # ------------------------------------------------------------------
        tree_model_for_check = (
            tree_cfg.get("model_nuc") if seq_type == "nucleotide"
            else tree_cfg.get("model_aa")
        )
        tree_key = {
            "msa_hash": hashlib.sha256(tree_input_fasta.read_bytes()).hexdigest()[:16],
            "tool":     tree_cfg["tool"],
            "model":    tree_model_for_check,
            "options":  tree_options,
        }
        if _check_checkpoint(tree_checkpoint, tree_key) and expected_treefile.exists():
            log.info("Reusing cached tree for tree_%s (iteration %d)", label, iteration)
            treefile = expected_treefile
        else:
            log.info("Running %s (%s model) for tree_%s (iteration %d) ...",
                     tree_cfg["tool"], tree_model_for_check, label, iteration)
            treefile = run_tree(
                alignment_fasta=tree_input_fasta,
                output_dir=tree_output_dir,
                prefix=f"tree_{label}",
                tool=tree_cfg["tool"],
                seq_type=seq_type,
                model_nuc=tree_cfg.get("model_nuc", "GTR+G"),
                model_aa=tree_cfg.get("model_aa", "LG+G"),
                options=tree_options,
                threads=threads,
            )
            validate_newick(treefile)
            _write_checkpoint(tree_checkpoint, tree_key)
            log.info("Tree inference complete for tree_%s (iteration %d)", label, iteration)

        # ------------------------------------------------------------------
        # Post-tree branch-length outlier removal (iterative)
        # ------------------------------------------------------------------
        if not do_outlier_removal or iteration >= max_iterations:
            break

        if len(sel_records) <= outlier_min_seqs:
            log.info(
                "Skipping branch-length outlier removal for tree_%s (iteration %d): "
                "%d sequences ≤ min_seqs threshold %d.",
                label, iteration, len(sel_records), outlier_min_seqs,
            )
            break

        all_outlier_ids = _find_branch_length_outliers(treefile, outlier_factor)
        if not all_outlier_ids:
            log.info("No branch-length outliers found for tree_%s (iteration %d) — stopping.",
                     label, iteration)
            break

        # RefSeqs are never removed — flag them with a warning instead.
        protected_outliers = all_outlier_ids & refseq_short_ids
        outlier_ids = all_outlier_ids - refseq_short_ids

        # Stats header (shared by warning + removal logs)
        _bl_stats = _branch_length_stats(treefile)
        _median   = _bl_stats["median"]
        _mad      = _bl_stats["mad"]
        _threshold = _median + outlier_factor * _mad
        _bl_map   = _bl_stats["bl_map"]

        log.info(
            "tree_%s iteration %d: branch-length stats — median=%.5f, MAD=%.5f, "
            "threshold (median + %.1f×MAD)=%.5f",
            label, iteration, _median, _mad, outlier_factor, _threshold,
        )

        for oid in protected_outliers:
            display_name = short_to_display.get(oid, oid)
            bl = _bl_map.get(oid, float("nan"))
            mads_above = (bl - _median) / _mad if _mad else float("nan")
            log.warning(
                "tree_%s iteration %d: RefSeq '%s' looks like a branch-length "
                "outlier — branch length %.5f (%.1f× MAD above median, "
                "threshold=%.5f) — KEEPING (RefSeq protected)",
                label, iteration, display_name, bl, mads_above, _threshold,
            )

        if not outlier_ids:
            log.info(
                "Only protected (RefSeq) outliers remain for tree_%s "
                "(iteration %d) — stopping.",
                label, iteration,
            )
            break

        # Only remove outliers that still leave at least outlier_min_seqs sequences
        remaining_after = len(sel_records) - len(outlier_ids)
        if remaining_after < outlier_min_seqs:
            log.info(
                "Skipping branch-length outlier removal for tree_%s (iteration %d): "
                "removing %d outlier(s) would leave only %d sequences (min_seqs=%d).",
                label, iteration, len(outlier_ids), remaining_after, outlier_min_seqs,
            )
            break

        for oid in outlier_ids:
            display_name = short_to_display.get(oid, oid)
            bl = _bl_map.get(oid, float("nan"))
            mads_above = (bl - _median) / _mad if _mad else float("nan")
            log.info(
                "tree_%s iteration %d: removing outlier '%s' — branch length %.5f "
                "(%.1f× MAD above median, threshold=%.5f)",
                label, iteration, display_name, bl, mads_above, _threshold,
            )

        n_outliers_removed += len(outlier_ids)
        sel_records = [r for r in sel_records if r.id not in outlier_ids]
        if len(sel_records) < 4:
            log.warning(
                "Only %d sequence(s) after branch-length outlier removal for tree_%s — stopping.",
                len(sel_records), label,
            )
            break

        # Clear checkpoints so the next iteration re-runs MSA + tree
        msa_checkpoint.unlink(missing_ok=True)
        tree_checkpoint.unlink(missing_ok=True)
        # Update metadata/id_map for the reduced set
        sel_metadata = [short_id_to_meta.get(r.id, {}) for r in sel_records]
        id_map = {r.id: short_to_display.get(r.id, r.id) for r in sel_records}

    restore_fasta_names(tree_input_fasta, family_dir / f"{family}_alignment_{label}.fasta", id_map)

    # Taxonomy annotation + rooting
    log.info("Inferring internal node taxonomy and rooting tree_%s ...", label)
    annotated_nwk = target_work / f"tree_{label}_annotated.nwk"
    lca_min_rank = family_cfg.get("taxonomy", {}).get("lca_min_rank", "none")
    bio_tree = annotate_tree(
        newick_path=treefile,
        id_map=short_to_display,
        metadata=sel_metadata,
        output_nwk=annotated_nwk,
        lca_min_rank=lca_min_rank,
    )
    log.info("Taxonomy annotation complete for tree_%s", label)

    # Write final Newick from in-memory tree — display names on leaves, no internal labels
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
    else:
        source_nwk = annotated_nwk if annotated_nwk.exists() else treefile
        restore_newick_names(source_nwk, output_nwk, short_to_display)

    # Assign leaf colors by genus / subfamily
    genus_inference = family_cfg.get("coloring", {}).get("genus_inference", "none")
    display_to_color, short_to_color, genus_to_color, subfamily_to_genera = assign_leaf_colors(
        sel_records, short_id_to_meta, short_to_display, genus_inference=genus_inference
    )
    n_subfamilies = sum(1 for k in subfamily_to_genera if k != "(unclassified)")

    # Collect summary stats
    tree_model = tree_cfg.get("model_nuc") if seq_type == "nucleotide" else tree_cfg.get("model_aa")
    tool_norm = tree_cfg["tool"].lower().replace("-", "").replace("_", "")
    if tool_norm in ("iqtree", "iqtree2") and is_model_finder_spec(tree_model or ""):
        chosen = parse_iqtree_best_model(tree_output_dir / f"tree_{label}.iqtree")
        if chosen:
            log.info("IQ-TREE ModelFinder chose %s (requested: %s) for tree_%s",
                     chosen, tree_model, label)
            tree_model = chosen
        else:
            log.warning("Could not parse chosen model from IQ-TREE log for tree_%s", label)
    tree_seq_lengths = [len(r.seq) for r in sel_records]
    target_stats: dict = {
        "leaves":                len(sel_records),
        "cluster_thresh_min":    cluster_thresh_min,
        "cluster_thresh_max":    cluster_thresh_max,
        "seq_type":              seq_type,
        "msa_tool":              msa_cfg["tool"],
        "msa_options":           msa_options,
        "trim_enabled":          trim_enabled,
        "trim_tool":             trim_tool if trim_enabled else "",
        "trim_options":          trim_options if trim_enabled else "",
        "tree_tool":             tree_cfg["tool"],
        "tree_model":            tree_model or "",
        "tree_options":          tree_options,
        "support_type":          support_type,
        "seq_lengths":           tree_seq_lengths,
        "display_to_color":      display_to_color,
        "genus_to_color":        genus_to_color,
        "subfamily_to_genera":   subfamily_to_genera,
        "n_outliers_removed":    n_outliers_removed,
        "n_length_outliers_long":  n_length_long,
        "n_length_outliers_short": n_length_short,
        "n_refseq_absorbed":     n_refseq_absorbed,
        "n_genera":              len(genus_to_color),
        "n_subfamilies":         n_subfamilies,
    }
    if trim_enabled:
        pre_stats = compute_msa_stats(msa_short_fasta)
        target_stats["msa_length_pre_trim"] = pre_stats.get("length", "")
    support_vals: list[float] = []
    if bio_tree is not None:
        target_stats["support"] = compute_support_stats(bio_tree)
        support_vals = [
            c.confidence for c in bio_tree.find_clades()
            if not c.is_terminal() and c.confidence is not None
        ]
    target_stats["msa"] = compute_msa_stats(tree_input_fasta)

    # PhyloXML
    log.info("Writing PhyloXML for tree_%s ...", label)
    phylogeny_name, phylogeny_detail = _build_phylogeny_name(
        family=family,
        seq_type=seq_type,
        region=family_cfg["sequence"]["region"],
        segment=family_cfg["sequence"].get("segment") or None,
        target_n=target_n,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        msa_tool=msa_cfg["tool"],
        msa_version=get_mafft_version(),
        msa_options=msa_options,
        tree_tool=tree_cfg["tool"],
        tree_version=get_tree_tool_version(tree_cfg["tool"]),
        tree_model=tree_model,
        tree_options=tree_options,
    )
    confidence_type = support_type
    write_phyloxml(
        newick_path=annotated_nwk,
        output_xml=family_dir / f"{family}_tree_{label}.xml",
        id_map=short_to_display,
        leaf_metadata={r.id: short_id_to_meta.get(r.id, {}) for r in sel_records},
        family=family,
        tree=bio_tree,
        phylogeny_name=phylogeny_name,
        phylogeny_detail=phylogeny_detail,
        confidence_type=confidence_type,
        leaf_colors=short_to_color,
    )

    return target_stats, support_vals, bio_tree


def _branch_length_stats(treefile: Path) -> dict:
    """Return median, MAD, and per-leaf branch-length map for a Newick file."""
    from Bio import Phylo as _Phylo
    try:
        bio_tree = next(iter(_Phylo.parse(str(treefile), "newick")))
    except Exception:
        return {"median": 0.0, "mad": 0.0, "bl_map": {}}
    bl_map = {
        c.name: c.branch_length for c in bio_tree.get_terminals()
        if c.name and c.branch_length is not None
    }
    bls = [b for b in bl_map.values() if b > 0]
    if len(bls) < 3:
        return {"median": 0.0, "mad": 0.0, "bl_map": bl_map}
    median_bl = statistics.median(bls)
    mad = statistics.median([abs(x - median_bl) for x in bls])
    return {"median": median_bl, "mad": mad, "bl_map": bl_map}


def _find_branch_length_outliers(treefile: Path, factor: float) -> set[str]:
    """Return leaf IDs whose branch length exceeds median + factor × MAD."""
    stats = _branch_length_stats(treefile)
    median_bl, mad = stats["median"], stats["mad"]
    if median_bl == 0 or mad == 0:
        return set()
    threshold = median_bl + factor * mad
    return {
        name for name, bl in stats["bl_map"].items()
        if bl is not None and bl > threshold
    }



def _build_phylogeny_name(
    family: str,
    seq_type: str,
    region: str,
    segment: str | None,
    target_n: int,
    threshold_min: float,
    threshold_max: float,
    msa_tool: str,
    msa_version: str,
    msa_options: str,
    tree_tool: str,
    tree_version: str,
    tree_model: str,
    tree_options: str,
) -> tuple[str, str]:
    """Return (short_name, detail_str) for PhyloXML name and description."""
    mol = "protein" if seq_type == "protein" else "nucleotide"
    if segment:
        region_str = segment
    elif region == "whole_genome":
        region_str = "whole genome"
    else:
        region_str = f"gene: {region}"

    tree_name = tree_tool.upper().replace("FASTTREE", "FastTree").replace("IQTREE2", "IQ-TREE").replace("IQTREE", "IQ-TREE")

    short_name = f"{family} [{mol}|{region_str}|{msa_tool.upper()}|{tree_name} {tree_model}]"

    msa_detail = f"{msa_tool.upper()} {msa_version} {msa_options}".strip()
    tree_detail = f"{tree_name} {tree_version} {tree_model}"
    if tree_options:
        tree_detail += f" {tree_options}"
    detail_str = (
        f"{mol}, {region_str} "
        f"(target={target_n} seqs, cluster id {threshold_min:.2f}\u2013{threshold_max:.2f}) "
        f"| {msa_detail} | {tree_detail}"
    )

    return short_name, detail_str


def _write_metadata_tsv(
    meta_rows: list[dict],
    short_to_display: dict[str, str],
    output_path: Path,
) -> None:
    if not meta_rows:
        return
    cols = ["short_id", "display_name", "accession", "species", "strain", "host",
            "collection_date", "location", "taxon_id", "length"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for meta in meta_rows:
            row = dict(meta)
            row["display_name"] = short_to_display.get(row.get("short_id", ""), "")
            writer.writerow(row)


def _save_config_copy(cfg: dict, path: Path) -> None:
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(clean, f, default_flow_style=False, sort_keys=False)


def _mark_done(family_dir: Path, family: str, output_dir: Path) -> None:
    (family_dir / ".status.json").write_text(json.dumps({"family": family, "status": "done"}))
    (output_dir / f".done_{family}").write_text("")


def _mark_skipped(family_dir: Path, family: str, output_dir: Path, reason: str) -> None:
    (family_dir / ".status.json").write_text(
        json.dumps({"family": family, "status": "skipped", "reason": reason})
    )
    (output_dir / f".done_{family}").write_text("")
