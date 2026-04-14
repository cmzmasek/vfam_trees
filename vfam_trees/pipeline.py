"""Per-family pipeline runner (called by Snakemake rules)."""
from __future__ import annotations

import copy
import csv
import json
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
    compute_bootstrap_stats, compute_msa_stats, compute_seqlen_stats,
    build_summary_row, write_summary_row,
)
from .quality import filter_sequences, deduplicate, write_fasta
from .rename import assign_short_ids, restore_fasta_names, restore_newick_names
from .subsample import adaptive_cluster_species, proportional_merge
from .msa import run_msa, get_mafft_version
from .tree import run_tree, get_tree_tool_version
from .taxonomy import annotate_tree
from .phyloxml_writer import write_phyloxml
from .logger import setup_logger, get_logger

log = get_logger(__name__)


def run_family(
    family: str,
    global_config_path: Path,
    configs_dir: Path,
    output_dir: Path,
    threads: int = 1,
    log_level: str = "INFO",
    summary_path: Path | None = None,
) -> None:
    """Execute the full pipeline for a single viral family."""
    family_dir = output_dir / family
    family_dir.mkdir(parents=True, exist_ok=True)

    log_file = family_dir / f"{family}.log"
    setup_logger(f"vfam_trees.{family}", log_file=log_file, level=log_level)
    log = get_logger(f"vfam_trees.{family}")

    log.info("=" * 60)
    log.info("Starting pipeline for %s", family)
    log.info("=" * 60)

    global_cfg = load_global_config(global_config_path)
    configure_entrez(
        email=global_cfg["ncbi"]["email"],
        api_key=global_cfg["ncbi"].get("api_key") or None,
    )

    # Fetch family-level taxonomy (taxid + ranked lineage) for the summary
    family_taxid = get_family_taxid(family)
    if family_taxid is not None:
        family_lineage_map = fetch_taxonomy_lineages([family_taxid])
        family_lineage = family_lineage_map.get(str(family_taxid), [])
    else:
        log.warning("Could not resolve NCBI taxid for family %s", family)
        family_lineage = []

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
            _mark_skipped(family_dir, family, "no species found in NCBI taxonomy")
            return
        with open(species_cache, "w") as f:
            json.dump(species_list, f, indent=2)
        log.info("Discovered %d species in %s", len(species_list), family)

    # -------------------------------------------------------------------------
    # Step 2: Per-species download + quality filter
    #         Results cached per species to avoid re-downloading
    # -------------------------------------------------------------------------
    species_data: dict[str, dict] = {}  # species_name → {records, metadata}

    log.info("Downloading sequences for %d species (max %d per species) ...",
             len(species_list), max_per_species)

    for i, sp in enumerate(species_list, start=1):
        sp_name = sp["name"]
        sp_taxid = sp["taxid"]
        sp_key = f"sp_{sp_taxid}"
        sp_work = work_dir / "species" / sp_key
        sp_work.mkdir(parents=True, exist_ok=True)

        gb_file = sp_work / "sequences.gb"

        # Download (cached)
        if not gb_file.exists():
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
                log.info("[%d/%d] No sequences found for %s — skipping.", i, len(species_list), sp_name)
                continue
        else:
            log.info("[%d/%d] Using cached sequences: %s", i, len(species_list), sp_name)

        # Parse + quality filter
        raw_records: list[SeqRecord] = []
        raw_meta: list[dict] = []
        for rec in parse_gb_records(gb_file):
            raw_records.append(rec)
            raw_meta.append(extract_metadata(rec))

        if not raw_records:
            continue

        acc_to_raw_meta = {m["accession"]: m for m in raw_meta}
        raw_records = deduplicate(raw_records)
        raw_meta = [acc_to_raw_meta[rec.id] for rec in raw_records if rec.id in acc_to_raw_meta]
        raw_records = [rec for rec in raw_records if rec.id in acc_to_raw_meta]

        filtered = filter_sequences(
            raw_records,
            seq_type=seq_type,
            min_length=quality_cfg["min_length"],
            max_ambiguous=quality_cfg["max_ambiguous"],
            exclude_organisms=quality_cfg.get("exclude_organisms"),
        )
        if not filtered:
            log.debug("No sequences passed quality filter for %s", sp_name)
            continue

        acc_to_meta = {m["accession"]: m for m in raw_meta}
        filtered_meta = [acc_to_meta[r.id] for r in filtered if r.id in acc_to_meta]
        filtered = [r for r in filtered if r.id in acc_to_meta]

        species_data[sp_name] = {"records": filtered, "metadata": filtered_meta}

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

    # Attach ranked lineage to each metadata dict; prefer NCBI over GenBank's
    # unranked lineage because the names are used for LCA matching.
    for data in species_data.values():
        for meta in data["metadata"]:
            ranked = taxid_to_lineage.get(str(meta.get("taxon_id", "")), [])
            meta["lineage_ranked"] = ranked
            if ranked:
                meta["lineage"] = [e["name"] for e in ranked]

    seqlen_stats = compute_seqlen_stats(
        [len(r.seq) for d in species_data.values() for r in d["records"]]
    )

    total_seqs = sum(len(d["records"]) for d in species_data.values())

    if total_seqs < 5:
        log.warning(
            "Only %d sequence(s) passed quality filters for %s (minimum 5) — skipping.",
            total_seqs, family,
        )
        _mark_skipped(family_dir, family, f"too few sequences after QC ({total_seqs} < 5)")
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
    acc_to_short: dict[str, str] = {}
    for orig_rec, renamed_rec in zip(all_records, renamed_records):
        acc_to_short[orig_rec.id] = renamed_rec.id

    for meta, renamed_rec in zip(all_metadata, renamed_records):
        meta["short_id"] = renamed_rec.id

    short_id_to_meta: dict[str, dict] = {
        rec.id: meta for rec, meta in zip(renamed_records, all_metadata)
    }

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
    for target_n, max_reps, label in targets:
        log.info("-" * 40)
        log.info("Building tree_%s  (target=%d sequences, max_reps_per_species=%d)",
                 label, target_n, max_reps)
        tree_stats[label] = _run_target(
            label=label,
            target_n=target_n,
            max_reps=max_reps,
            species_data=species_data,
            short_id_to_meta=short_id_to_meta,
            short_to_display=short_to_display,
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

    log.info("=" * 60)
    log.info("Pipeline complete for %s", family)
    log.info("=" * 60)
    _mark_done(family_dir, family)

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
            tree_stats=tree_stats,
        )
        write_summary_row(summary_path, row)


def _run_target(
    label: str,
    target_n: int,
    max_reps: int,
    species_data: dict,
    short_id_to_meta: dict[str, dict],
    short_to_display: dict[str, str],
    seq_type: str,
    clustering_tool: str,
    clustering_cfg: dict,
    family_cfg: dict,
    family: str,
    family_dir: Path,
    work_dir: Path,
    threads: int,
    log,
) -> dict:
    target_work = work_dir / label
    target_work.mkdir(exist_ok=True)

    threshold_min = clustering_cfg.get("threshold_min", 0.70)
    threshold_max = clustering_cfg.get("threshold_max", 0.99)

    # Adaptive per-species clustering
    log.info("Clustering sequences per species (%s, max_reps=%d, search range=%.2f–%.2f) ...",
             clustering_tool, max_reps, threshold_min, threshold_max)
    species_reps: dict[str, list[SeqRecord]] = {}
    for sp_name, data in species_data.items():
        sp_safe = sp_name.replace(" ", "_").replace("/", "_")
        sp_work = target_work / "clustering" / sp_safe
        reps, threshold_used = adaptive_cluster_species(
            records=data["renamed"],
            max_reps=max_reps,
            threshold_min=threshold_min,
            threshold_max=threshold_max,
            seq_type=seq_type,
            work_dir=sp_work,
            clustering_tool=clustering_tool,
        )
        species_reps[sp_name] = reps

    total_reps = sum(len(r) for r in species_reps.values())
    log.info("Clustering complete: %d total representatives across %d species",
             total_reps, len(species_reps))

    # Proportional merge
    sel_records = proportional_merge(species_reps, target_n)
    log.info("Proportional merge: selected %d sequences for tree_%s", len(sel_records), label)

    if not sel_records:
        log.warning("No sequences selected for %s tree — skipping.", label)
        return {}

    sel_metadata = [short_id_to_meta.get(r.id, {}) for r in sel_records]
    id_map = {r.id: short_to_display.get(r.id, r.id) for r in sel_records}

    # Write raw quality-controlled FASTA (short IDs → restored names)
    raw_short_fasta = target_work / f"sequences_raw_{label}_short.fasta"
    write_fasta(sel_records, raw_short_fasta)
    restore_fasta_names(raw_short_fasta, family_dir / f"{family}_sequences_raw_{label}.fasta", id_map)

    # Metadata TSV
    _write_metadata_tsv(sel_metadata, short_to_display, family_dir / f"{family}_metadata_{label}.tsv")

    # MSA
    msa_short_fasta = target_work / f"alignment_{label}_short.fasta"
    msa_cfg = family_cfg["msa"]
    log.info("Running MAFFT (%s) on %d sequences for tree_%s ...",
             msa_cfg["options"], len(sel_records), label)
    run_msa(
        input_fasta=raw_short_fasta,
        output_fasta=msa_short_fasta,
        tool=msa_cfg["tool"],
        options=msa_cfg["options"],
        threads=threads,
    )
    log.info("MAFFT complete for tree_%s", label)
    restore_fasta_names(msa_short_fasta, family_dir / f"{family}_alignment_{label}.fasta", id_map)

    # Tree inference
    tree_cfg = family_cfg[f"tree_{label}"]
    log.info("Running %s (%s model) for tree_%s ...",
             tree_cfg["tool"],
             tree_cfg.get("model_nuc") if seq_type == "nucleotide" else tree_cfg.get("model_aa"),
             label)
    treefile = run_tree(
        alignment_fasta=msa_short_fasta,
        output_dir=target_work / "tree",
        prefix=f"tree_{label}",
        tool=tree_cfg["tool"],
        seq_type=seq_type,
        model_nuc=tree_cfg.get("model_nuc", "GTR+G"),
        model_aa=tree_cfg.get("model_aa", "WAG+G"),
        options=tree_cfg.get("options", ""),
        threads=threads,
    )

    log.info("Tree inference complete for tree_%s", label)

    # Taxonomy annotation + rooting
    log.info("Inferring internal node taxonomy and rooting tree_%s ...", label)
    annotated_nwk = target_work / f"tree_{label}_annotated.nwk"
    bio_tree = annotate_tree(
        newick_path=treefile,
        id_map=short_to_display,
        metadata=sel_metadata,
        output_nwk=annotated_nwk,
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
                clade.name = None  # taxonomy labels live in PhyloXML, not Newick
        Phylo.write(nwk_tree, str(output_nwk), "newick")
        log.info("Newick written to %s", output_nwk)
    else:
        source_nwk = annotated_nwk if annotated_nwk.exists() else treefile
        restore_newick_names(source_nwk, output_nwk, short_to_display)

    # Collect summary stats before writing outputs
    target_stats: dict = {"leaves": len(sel_records)}
    if bio_tree is not None:
        target_stats["bs"] = compute_bootstrap_stats(bio_tree)
    target_stats["msa"] = compute_msa_stats(msa_short_fasta)

    # PhyloXML — pass the in-memory tree to avoid lossy Newick roundtrip
    log.info("Writing PhyloXML for tree_%s ...", label)
    phylogeny_name = _build_phylogeny_name(
        family=family,
        seq_type=seq_type,
        region=family_cfg["sequence"]["region"],
        segment=family_cfg["sequence"].get("segment") or None,
        target_n=target_n,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        msa_tool=msa_cfg["tool"],
        msa_version=get_mafft_version(),
        msa_options=msa_cfg["options"],
        tree_tool=tree_cfg["tool"],
        tree_version=get_tree_tool_version(tree_cfg["tool"]),
        tree_model=tree_cfg.get("model_nuc") if seq_type == "nucleotide" else tree_cfg.get("model_aa"),
        tree_options=tree_cfg.get("options", ""),
    )
    write_phyloxml(
        newick_path=annotated_nwk,
        output_xml=family_dir / f"{family}_tree_{label}.xml",
        id_map=short_to_display,
        leaf_metadata={r.id: short_id_to_meta.get(r.id, {}) for r in sel_records},
        family=family,
        tree=bio_tree,
        phylogeny_name=phylogeny_name,
    )

    return target_stats


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
) -> str:
    mol = "protein" if seq_type == "protein" else "nucleotide"
    if segment:
        region_str = segment
    elif region == "whole_genome":
        region_str = "whole genome"
    else:
        region_str = f"gene: {region}"

    msa_str = f"{msa_tool.upper()} {msa_version} {msa_options}".strip()
    tree_name = tree_tool.upper().replace("FASTTREE", "FastTree").replace("IQTREE2", "IQ-TREE").replace("IQTREE", "IQ-TREE")
    tree_str = f"{tree_name} {tree_version} {tree_model}"
    if tree_options:
        tree_str += f" {tree_options}"

    return (
        f"{family} [{mol}, {region_str} "
        f"(target={target_n} seqs, cluster id {threshold_min:.2f}\u2013{threshold_max:.2f}) "
        f"| {msa_str} | {tree_str}]"
    )


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


def _mark_done(family_dir: Path, family: str) -> None:
    (family_dir / ".status.json").write_text(json.dumps({"family": family, "status": "done"}))


def _mark_skipped(family_dir: Path, family: str, reason: str) -> None:
    (family_dir / ".status.json").write_text(
        json.dumps({"family": family, "status": "skipped", "reason": reason})
    )
