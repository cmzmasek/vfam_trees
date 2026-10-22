"""vfam_trees command-line interface."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import click
import yaml

from . import __version__

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"
DEFAULT_GLOBAL_CFG = Path("config/global.yaml")
DEFAULT_CONFIGS_DIR = Path("configs")
DEFAULT_OUTPUT_DIR = Path("results")


@click.group()
@click.version_option(version=__version__, prog_name="vfam_trees")
def main():
    """vfam_trees — phylogenetic trees for viral families.

    Builds maximum-likelihood phylogenetic trees for one or more viral
    families using sequences downloaded from NCBI. For each family, two
    trees are produced: a broad 500-sequence tree (IQ-TREE, --fast for
    nucleotide, a 50-iteration search for protein) and a more collapsed
    100-sequence tree (IQ-TREE full search with UFBoot).

    \b
    Typical workflow:
      1. vfam_trees init [--full], then edit config/global.yaml
         (set NCBI email and API key)
      2. vfam_trees test          (verify all dependencies)
      3. vfam_trees init-configs -f families.txt   (review/edit per-family configs)
      4. vfam_trees run -f families.txt -j 4       (run the pipeline)
      5. vfam_trees status -f families.txt         (check progress)

    \b
    Examples:
      vfam_trees run -f families.txt
      vfam_trees run -f families.txt -j 4 -t 4 --verbose
      vfam_trees run -f families.txt --force
      vfam_trees status -f families.txt
      vfam_trees init-configs -f families.txt
      vfam_trees test
    """


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
@main.command()
@click.option(
    "--families", "-f",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Text file with one viral family name per line.",
)
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
    help="Path to global.yaml.",
)
@click.option(
    "--configs-dir", "-c",
    default=DEFAULT_CONFIGS_DIR,
    show_default=True,
    type=click.Path(path_type=Path),
    help="Directory containing per-family YAML configs.",
)
@click.option(
    "--output-dir", "-o",
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option(
    "--cores", "-j",
    default=1,
    show_default=True,
    type=int,
    help="Number of parallel family jobs (Snakemake --cores).",
)
@click.option(
    "--threads", "-t",
    default=1,
    show_default=True,
    type=int,
    help="Threads per family job (MAFFT, IQ-TREE).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Rerun families that were already processed.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Preview per-family config parameters without running the pipeline.",
)
@click.option(
    "--verbose", "log_level",
    flag_value="DEBUG",
    help="Verbose output (DEBUG level).",
)
@click.option(
    "--quiet", "log_level",
    flag_value="WARNING",
    help="Suppress informational output.",
)
@click.option(
    "--log-level",
    "log_level",
    default="INFO",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING"], case_sensitive=False),
    help="Logging verbosity.",
)
def run(
    families: Path,
    global_config: Path,
    configs_dir: Path,
    output_dir: Path,
    cores: int,
    threads: int,
    force: bool,
    dry_run: bool,
    log_level: str,
):
    """Run the full pipeline for one or more viral families.

    For each family the pipeline will:
    discover species (NCBI taxonomy) → download per species (RefSeq-priority) →
    quality filter → adaptive per-species clustering → proportional merge
    (RefSeqs preferred) → pre-MSA length-outlier removal →
    align (MAFFT) → trim alignment columns (trimAl) →
    infer trees (IQ-TREE: tree_500 --fast for nucleotide, -n 50 for protein;
    full search for tree_100) →
    iterative post-tree branch-length outlier removal (re-align + re-tree;
    rounds checked on quick --fast trees) →
    annotate internal nodes (LCA taxonomy) → write Newick, PhyloXML,
    Auspice v2 JSON (viewable at auspice.us / Nextclade / `auspice view`),
    PDF/PNG tree images, and per-family report.

    What is aligned depends on the family. Most RNA-virus families (107)
    are aligned on the RNA-dependent RNA polymerase (RdRp), the
    reverse-transcribing families (4) on the reverse transcriptase (RT) and
    the CRESS-DNA families (23) on the Rep protein: each genome's marker is
    found by HMMER profile homology against bundled Pfam profiles (the RdRp
    and RT cut to their domain), so these trees are protein trees. Large DNA
    viruses use curated marker proteins (single or concatenated), found by
    GenBank protein name or, where the family sets `hmm.enabled`, by
    profile homology. Other families use whole genomes (or one segment) as
    nucleotide. HMMER (hmmscan, hmmpress) must be on PATH for the default
    family list.

    MSA and tree-inference steps are content-hash-checkpointed, so reruns
    skip stages whose inputs and config are unchanged. A cross-family
    overview PNG (overview_tree_100.png) is generated automatically when
    the run completes.

    Already-processed families are skipped unless --force is used.
    Per-family configs are auto-generated in configs/ if missing.

    \b
    Examples:
      vfam_trees run -f families.txt
      vfam_trees run -f families.txt -j 4 -t 4
      vfam_trees run -f families.txt --dry-run
      vfam_trees run -f families.txt --force --verbose
      vfam_trees run -f families.txt -o /data/results -g /etc/vfam_global.yaml
    """
    family_list = _read_family_list(families)
    if not family_list:
        click.echo("No families found in input file.", err=True)
        sys.exit(1)

    if dry_run:
        _print_dry_run(family_list, global_config, configs_dir, output_dir, threads)
        sys.exit(0)

    # Checked here so a bad global.yaml stops the run before Snakemake starts
    # (each family job would fail on it); unknown keys are logged (WARNING).
    from .config import load_global_config
    try:
        global_cfg = load_global_config(global_config)
    except (FileNotFoundError, ValueError, yaml.YAMLError) as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    missing = _missing_tools(family_list, configs_dir, global_cfg)
    if missing:
        click.echo(
            "Error: this run needs tools that are not on PATH: "
            f"{', '.join(missing)}. Run `vfam_trees test` for details.",
            err=True,
        )
        sys.exit(1)

    click.echo(f"Families to process: {', '.join(family_list)}")

    if force:
        _clear_status(family_list, output_dir)

    snakemake_cfg = {
        "families": family_list,
        "global_config": str(global_config.resolve()),
        "configs_dir": str(configs_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
        "log_level": log_level.upper(),
        "threads": threads,
    }

    # Snakemake clamps a rule's per-job threads to --cores. To let
    # `-j N -t M` actually run N families × M threads each, --cores must be
    # at least N*M; otherwise rule run_family (threads: M) silently falls
    # back to --cores threads, and MAFFT/IQ-TREE only receive that many.
    total_cores = cores * threads
    cmd = [
        "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--cores", str(total_cores),
        "--config", *[f"{k}={_yaml_encode(v)}" for k, v in snakemake_cfg.items()],
        "--rerun-incomplete",
        "--printshellcmds",
    ]

    from .provenance import finish_run_record, start_run_record
    record = start_run_record(
        output_dir,
        argv=["vfam_trees", *sys.argv[1:]],
        families_file=families,
        families=family_list,
        options={"cores": cores, "threads": threads, "force": force,
                 "log_level": log_level.upper()},
        paths={"global_config": str(global_config.resolve()),
               "configs_dir": str(configs_dir.resolve()),
               "output_dir": str(output_dir.resolve())},
        global_cfg=global_cfg,
    )
    click.echo(f"Run record: {record}")

    click.echo(f"Running: {' '.join(cmd)}")
    exit_code = None
    try:
        result = subprocess.run(cmd)
        exit_code = result.returncode
    finally:
        finish_run_record(record, exit_code)

    if result.returncode == 0:
        from .report import generate_overview_png
        overview_path = output_dir / "overview_tree_100.png"
        click.echo(f"Generating overview PNG: {overview_path}")
        generate_overview_png(output_dir, overview_path)

    sys.exit(result.returncode)


def _missing_tools(family_list: list[str], configs_dir: Path, global_cfg: dict) -> list[str]:
    """External tools this run needs that are not on PATH: REQUIRED_TOOLS
    always, HMMER when a family identifies its marker by HMM (the RdRp and
    other HMM presets), FastTree when a tree block uses it.  Checked before
    Snakemake starts, so a missing tool stops the run at once instead of
    failing each family after its downloads."""
    import shutil
    from .config import preview_family_config

    uses_hmm = uses_fasttree = False
    for family in family_list:
        try:
            cfg = preview_family_config(family, configs_dir, global_cfg)
        except Exception:
            continue  # a broken family config is reported by that family's job
        uses_hmm |= bool((cfg.get("hmm") or {}).get("enabled"))
        uses_fasttree |= any(
            str((cfg.get(block) or {}).get("tool", "")).lower() == "fasttree"
            for block in ("tree_500", "tree_100")
        )
    needed = list(REQUIRED_TOOLS)
    if uses_hmm:
        needed += HMM_TOOLS
    if uses_fasttree:
        needed += FASTTREE_TOOLS
    return [tool for tool in needed if shutil.which(tool) is None]


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------
@main.command()
@click.option(
    "--families", "-f",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Text file with one viral family name per line.",
)
@click.option(
    "--output-dir", "-o",
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    type=click.Path(path_type=Path),
)
def status(families: Path, output_dir: Path):
    """Show processing status for each family.

    Reports one of: pending, in-progress (with current stage inferred
    from checkpoint files — downloading/QC, MSA, tree inference, or
    annotating), done, skipped (with reason), or error reading status.

    \b
    Examples:
      vfam_trees status -f families.txt
      vfam_trees status -f families.txt -o /data/results
    """
    family_list = _read_family_list(families)
    rows = []
    for family in family_list:
        status_file = _find_status_file(family, output_dir)
        if status_file is not None:
            try:
                state = json.loads(status_file.read_text())
                st = state.get("status", "unknown")
                reason = state.get("reason", "")
            except Exception:
                st, reason = "error reading status", ""
        else:
            family_dir = _find_family_dir(family, output_dir)
            if family_dir is not None:
                st = "in-progress"
                reason = _detect_stage(family_dir)
            else:
                st, reason = "pending", ""
        rows.append((family, st, reason))

    w = max(len(r[0]) for r in rows) + 2
    click.echo(f"{'Family':<{w}}  {'Status':<14}  Stage / Notes")
    click.echo("-" * 70)
    for family, st, reason in rows:
        click.echo(f"{family:<{w}}  {st:<14}  {reason}")


# ---------------------------------------------------------------------------
# init-configs
# ---------------------------------------------------------------------------
@main.command("init-configs")
@click.option(
    "--families", "-f",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Text file with one viral family name per line.",
)
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
)
@click.option(
    "--configs-dir", "-c",
    default=DEFAULT_CONFIGS_DIR,
    show_default=True,
    type=click.Path(path_type=Path),
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite existing configs with current defaults (user edits will be lost).",
)
def init_configs(families: Path, global_config: Path, configs_dir: Path, force: bool):
    """Generate default per-family configs without running the pipeline.

    Creates configs/{Family}.yaml for any family that does not already
    have one. Existing configs are left untouched unless --force is used.
    Run this before the pipeline to review and tune parameters (clustering
    thresholds, models, max sequences per species, etc.) before committing
    to a full run.

    \b
    Examples:
      vfam_trees init-configs -f families.txt
      vfam_trees init-configs -f families.txt --force
      vfam_trees init-configs -f families.txt -c /data/configs
    """
    from .config import load_global_config, load_family_config, make_minimal_global_cfg

    family_list = _read_family_list(families)
    if global_config.exists():
        global_cfg = load_global_config(global_config)
    else:
        click.echo(
            f"Notice: no global config found at {global_config} — "
            "using hardcoded defaults for all families."
        )
        global_cfg = make_minimal_global_cfg()

    for family in family_list:
        config_path = configs_dir / f"{family}.yaml"
        if force and config_path.exists():
            config_path.unlink()
        cfg, auto = load_family_config(family, configs_dir, global_cfg)
        if auto:
            click.echo(f"{'Regenerated' if force else 'Generated'}: {config_path}")
        else:
            click.echo(f"Exists:      {config_path}")


# ---------------------------------------------------------------------------
# init
# ---------------------------------------------------------------------------

GLOBAL_CONFIG_TEMPLATE = """\
# vfam_trees global configuration
# Edit this file before running the pipeline.

ncbi:
  email: "your.email@example.com"   # REQUIRED — your email address for NCBI Entrez API
  api_key: ""                        # Optional but recommended — get one at:
                                     # https://www.ncbi.nlm.nih.gov/account/

# The output folder and log level are command-line options of
# `vfam_trees run` (-o, --log-level), not settings in this file.

# Optional external family-annotation TSV.  Must have a 'family' column plus
# the extra columns to join into summary.tsv (currently: baltimore_class).
# null → virus_families_annotation.tsv next to this file; relative paths are
# resolved against this file's directory.  A missing file leaves the columns
# blank.
annotation_tsv: null

# Optional shared sequence download cache.
# Sequences are stored keyed by (taxid, db, region, segment, max_per_species)
# and reused across runs and families on the same machine — or on a shared
# filesystem across machines.  Remove or set dir to null to disable.
cache:
  dir: ~/.vfam_cache    # path to cache root; ~ is expanded automatically
  ttl_days: 90          # null = never expire; set to re-download after N days

# Visualization settings — applied to all tree PDF/PNG outputs.
visualization:
  branch_linewidth: 0.5        # branch line width in pt (BioPython default is ~1.5 pt)

# Per-family icon PNG (tree_100 topology only, no labels, square).
icon:
  size: 256                    # output size in pixels (square)
  bg_color: "#EAF3F2"          # background color (hex)
  branch_color: "#000000"      # branch color (hex; uniform, no genus coloring)

# Default values applied when per-family config is auto-generated.
# Override any of these here to change the default for all families.
defaults:
  download:
    max_per_species: 300

  # Family presets replace these for most families: RNA-virus families get
  # the HMM-identified RdRp protein, reverse-transcribing families the RT,
  # CRESS-DNA families the Rep, large DNA viruses marker proteins, and
  # segmented families their segment (README: "DNA virus families",
  # "RNA virus families (RdRp)", "Reverse-transcribing and CRESS-DNA
  # families").
  sequence:
    type: nucleotide    # nucleotide or protein
    region: whole_genome
    segment: null
    hmms: []            # HMM profile tokens of the marker (hmm.enabled); set by presets
    domain_padding: null  # HMM mode: keep the marker's profile hits ± this many
                          # residues (100 for RdRp and RT); null keeps whole proteins

  quality:
    min_length: null    # null = auto (50% of median sequence length)
    max_length: null    # null = no upper limit
    max_ambiguous: 0.01
    exclude_organisms:
      - synthetic construct
      - metagenome
      - "MAG:"
      - uncultured
      - unverified
      - vector
      - recombinant
      - patent

  clustering:
    tool: mmseqs2
    threshold_min: 0.70
    threshold_max: 0.99
    max_reps_500: 20
    max_reps_100: 5

  targets:
    max_500: 500
    max_100: 100

  # --adjustdirection (nucleotide only) reverse-complements GenBank records
  # deposited on the opposite strand; the pipeline strips MAFFT's "_R_"
  # name prefix and logs the flipped ids.
  msa_500:
    tool: mafft
    options_nuc: "--retree 2 --adjustdirection"
    options_aa: "--auto"

  msa_100:
    tool: mafft
    options_nuc: "--retree 3 --adjustdirection"
    options_aa: "--maxiterate 1000 --localpair"

  # Column trimming applied between MAFFT and tree inference (both tree_500
  # and tree_100).  Removes poorly-aligned positions that add noise to the
  # tree.  -gt 0.5 keeps columns with at most 50% gaps, for nucleotide and
  # protein alignments alike.  (-automated1, the default before 1.6.0, keeps
  # few columns on alignments with many short fragments.)
  msa_trim:
    enabled: true
    tool: trimal
    options: "-gt 0.5"

  # tree_500: IQ-TREE; the wrapper auto-adds -alrt 1000 for SH-aLRT support.
  # options_nuc applies to nucleotide runs, options_aa to protein runs
  # (marker proteins, RdRp / RT / Rep).  Nucleotide: --fast (same wall time
  # as FastTree, clearly better likelihood).  Protein: -n 50, the normal
  # search stopped after 50 iterations (--fast lands in a different local
  # optimum on each seed for weak-signal protein sets; ~12x --fast time).
  # FastTree remains available: set tool: fasttree AND options_nuc: "" /
  # options_aa: "" (IQ-TREE options such as --fast and -n are rejected at load).
  tree_500:
    tool: iqtree
    options_nuc: "--fast"
    options_aa: "-n 50"
    model_nuc: GTR+G
    model_aa: LG+G

  # tree_100: full IQ-TREE search with SH-aLRT and UFBoot (--bnni curbs
  # UFBoot overestimation under model violation).  Outputs carry the UFBoot
  # value; the SH-aLRT value stays in the IQ-TREE .treefile under _work.
  tree_100:
    tool: iqtree
    options_nuc: "-alrt 1000 -B 1000 --bnni"
    options_aa: "-alrt 1000 -B 1000 --bnni"
    model_nuc: GTR+G
    model_aa: TEST

  # Pre-MSA sequence-length outlier removal (applied after clustering +
  # proportional merge, before MAFFT).  Drops sequences whose length is a
  # MAD-on-log outlier — keep window = exp(median(log L) ± k·σ_log), with
  # σ_log = 1.4826 × MAD(log L).  Adapts to each family's natural length
  # variation.  A hard floor [min_lo_mult, max_hi_mult] × median is unioned
  # with the MAD window so tight families don't lose moderate truncations.
  length_outlier:
    enabled: true
    k: 5.0             # MAD-units on either side of log-median (0 disables MAD)
    min_lo_mult: 0.20  # never drop seqs ≥ this × median (0 disables lower floor)
    max_hi_mult: 5.0   # never drop seqs ≤ this × median (0 disables upper floor)

  # Post-tree branch-length outlier removal (iterative).
  # Drops leaves whose terminal branch length exceeds median + factor × MAD.
  # quick_check: each round is checked on a quick IQ-TREE --fast tree first;
  # the configured search runs once that finds nothing to remove (and its
  # tree is checked as before).  false: every round runs the full search.
  outlier_removal:
    enabled: true
    factor: 20.0
    max_iterations: 3
    min_seqs: 40
    quick_check: true

  # Leaf color inference strategy for sequences that have no formal genus rank
  # in their NCBI lineage.
  #   none    — unclassified sequences are drawn in grey (default, conservative).
  #   suffix  — treat any single-word name ending in "virus" as a genus proxy,
  #             regardless of NCBI's rank annotation (e.g. "Flavivirus" marked
  #             as "no rank" still gets its own color band).
  #   deepest — suffix match first; if that fails, use the deepest lineage entry
  #             above species level as a grouping key (maximally aggressive —
  #             useful for poorly annotated families).
  # Leaf label content — controls what appears in PhyloXML <name> elements,
  # PDF/PNG tree images, and FASTA display names.
  labeling:
    # Format string for leaf labels.  Use {placeholder} tokens from:
    #   {species}   NCBI species name
    #   {id}        accession number
    #   {host}      host organism
    #   {strain}    strain / isolate name
    #   {location}  collection country / region
    #   {year}      4-digit collection year (extracted from collection_date)
    #   {genus}     genus (when populated in metadata; empty otherwise)
    # Any literal text between tokens — including separators like "|" — is
    # reproduced verbatim.  Default: "{species}|{id}|{host}"
    format: "{species}|{id}|{host}"
    # Replace spaces in field values with underscores (keeps Newick parseable).
    replace_whitespace: true
    # When false (default): empty / "unknown" / "n/a" fields AND their
    # immediately preceding separator are dropped so no leading or consecutive
    # separators appear.  Set true to keep separators regardless.
    keep_separator_on_empty: false

  coloring:
    genus_inference: deepest

  # Minimum lineage depth a leaf must reach to participate in the LCA vote
  # used for internal-node taxonomy annotation.  Leaves whose lineage ends
  # above this rank are excluded from the vote — they remain in the tree but
  # do not blur internal-node annotations (e.g. a sequence filed only at
  # subfamily level will not drag a genus-level internal node up to subfamily).
  #   none      — all leaves vote (default, current behaviour).
  #   subfamily — exclude leaves whose lineage ends at subfamily or above.
  #   genus     — exclude leaves whose lineage does not reach genus level.
  #   species   — require species-level lineage to participate.
  taxonomy:
    lca_min_rank: species
"""


# Appended to GLOBAL_CONFIG_TEMPLATE by `vfam_trees init --full`, so the file
# lists every setting with its default.  Indented: these blocks continue the
# `defaults:` mapping that ends GLOBAL_CONFIG_TEMPLATE.
GLOBAL_CONFIG_TEMPLATE_ADVANCED = """\

  # ---------------------------------------------------------------------
  # Advanced settings.  Usually left alone or set per family; listed so
  # this file shows every default of this version.
  # ---------------------------------------------------------------------

  # Pre-clustering: drop non-RefSeq sequences ≥ threshold identical to a
  # RefSeq of the same species.  All RefSeqs are kept.
  refseq_absorption:
    enabled: true
    threshold: 0.99

  output:
    auspice_json: true           # Auspice v2 JSON next to each PhyloXML
    phyloxml_font_colors: false  # genus colour as style:font_color on PhyloXML tips
    archaeopteryx_html: true     # {prefix}_tree_{500,100}.html (Archaeopteryx.js page)

  # HMMER profile identification of marker proteins.  Off here; the RdRp,
  # RT, Rep, Poxviridae and Adenoviridae presets turn it on with their own
  # settings.
  hmm:
    enabled: false
    database: null               # bundled set name (e.g. RdRp), a directory of .hmm files, or one .hmm file
    default_evalue: 1.0e-05      # domain E-value gate for profiles without a curated GA
    use_ga_when_available: true  # prefer each profile's Pfam gathering threshold
    relative_length_cutoff: 0.5  # minimum fraction of the profile a hit must cover
    multidomain_overlap_tolerance: 30  # residues adjacent domains of an "A--B" token may overlap
    best_hit_only: false         # count a protein only if a marker profile is its best hit
    best_hit_also: []            # other profiles that count as the marker's own for best_hit_only
    min_protein_length: null     # proteins shorter than this (aa) are never candidates
    threads: null                # hmmscan --cpu; null → the run's thread count

  # Concatenated multi-marker mode (sequence.region: concatenated); the
  # Poxviridae, Orthoherpesviridae and other large-DNA-virus presets fill
  # in the markers.
  concatenation:
    proteins: []                 # marker specs: {name, aliases, hmms, length_range, ...}
    min_fraction: 0.7            # drop genomes with fewer than this fraction of markers
    source_nuc_min_length_frac: 0.3  # drop proteins whose source record is < this × the longest
    partition_tree_100: true     # one model per marker in tree_100
    partition_tree_500: false    # tree_500 is single-model

  # Curator overrides; normally set per family (see README).
  manual:
    include: []                  # accessions to force in (bypass QC)
    include_seq: []              # pasted sequences: {id, organism, sequence}
    include_fasta_files: []      # FASTA files injected like include_seq
    exclude: []                  # accessions to drop
    restrict_to_lineages: []     # only species under these taxa (names or taxids)
    name: ""                     # display name and output prefix; "" → family name

  # Forced-include sequences from an external database (Pathoplexus);
  # normally set per family (see README).
  additional_data_sources: []
"""


@main.command("init")
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
    help="Path to write the global config template.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite an existing config file.",
)
@click.option(
    "--full",
    is_flag=True,
    default=False,
    help="Also list the advanced settings (HMM, concatenation, output, "
         "manual), so the file shows every default of this version.",
)
def init(global_config: Path, force: bool, full: bool):
    """Generate a template global config file.

    Creates config/global.yaml (or the path given by --global-config)
    with the common default values and comments explaining each setting;
    --full lists every setting.  Edit the file — especially ncbi.email —
    before running the pipeline.

    Every run also records the settings it used: see run_records/ in the
    output folder and <Family>.yaml in each family folder.

    \b
    Examples:
      vfam_trees init
      vfam_trees init --full
      vfam_trees init --force
      vfam_trees init -g /path/to/my_global.yaml
    """
    if global_config.exists() and not force:
        click.echo(
            f"Config already exists: {global_config}\n"
            "Use --force to overwrite.",
            err=True,
        )
        sys.exit(1)

    global_config.parent.mkdir(parents=True, exist_ok=True)
    global_config.write_text(
        GLOBAL_CONFIG_TEMPLATE + (GLOBAL_CONFIG_TEMPLATE_ADVANCED if full else "")
    )
    click.echo(f"Created: {global_config}")
    click.echo("Edit the file and set ncbi.email before running the pipeline.")


# ---------------------------------------------------------------------------
# test
# ---------------------------------------------------------------------------

REQUIRED_TOOLS = ["mafft", "iqtree2", "mmseqs", "trimal"]
OPTIONAL_TOOLS = ["cd-hit", "cd-hit-est"]
FASTTREE_TOOLS = ["FastTree"]  # needed only when a tree block sets tool: fasttree
# HMMER: the RNA-virus RdRp presets and the Poxviridae / Adenoviridae presets
# identify their markers by HMM, so the default family list needs it.
HMM_TOOLS = ["hmmscan", "hmmpress"]
REQUIRED_PACKAGES = ["Bio", "click", "yaml", "snakemake", "requests"]


@main.command()
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
)
def test(global_config: Path):
    """Check that all dependencies and configuration are in place.

    Verifies:
      - Required Python packages (biopython, snakemake, taxopy, ...)
      - Required external tools (mafft, iqtree2, mmseqs, trimal)
      - Optional tools (cd-hit, cd-hit-est — only needed when clustering.tool: cdhit)
      - FastTree (only needed when tree_500.tool or tree_100.tool is fasttree)
      - HMMER tools (hmmscan, hmmpress — required by the HMM marker presets:
        the RNA-virus RdRp families, Poxviridae and Adenoviridae)
      - global.yaml exists and has ncbi.email set
      - Live connectivity to NCBI Entrez

    Exits with code 0 if all required checks pass, 1 otherwise.

    \b
    Examples:
      vfam_trees test
      vfam_trees test -g /path/to/global.yaml
    """
    import shutil
    ok = True

    click.echo("=== vfam_trees dependency check ===\n")

    # Python packages
    click.echo("Python packages:")
    for pkg in REQUIRED_PACKAGES:
        try:
            __import__(pkg)
            click.echo(f"  [OK]  {pkg}")
        except ImportError:
            click.echo(f"  [MISSING]  {pkg}")
            ok = False

    click.echo("")

    # Required external tools
    click.echo("Required external tools:")
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool)
        if path:
            click.echo(f"  [OK]  {tool}  ({path})")
        else:
            click.echo(f"  [MISSING]  {tool}")
            ok = False

    click.echo("")

    # Optional tools
    click.echo("Optional external tools (only used when clustering.tool: cdhit):")
    for tool in OPTIONAL_TOOLS:
        path = shutil.which(tool)
        if path:
            click.echo(f"  [OK]  {tool}  ({path})")
        else:
            click.echo(f"  [--]  {tool}  (not found — only needed if clustering.tool: cdhit)")

    click.echo("")

    # FastTree (optional since v1.5.0 — both trees default to IQ-TREE)
    click.echo("FastTree (only used when a tree block sets tool: fasttree):")
    for tool in FASTTREE_TOOLS:
        path = shutil.which(tool)
        if path:
            click.echo(f"  [OK]  {tool}  ({path})")
        else:
            click.echo(f"  [--]  {tool}  (not found — only needed if tool: fasttree)")

    click.echo("")

    # HMMER — required: the RdRp presets (most RNA-virus families) and the
    # Poxviridae / Adenoviridae presets identify markers by HMM.
    click.echo("HMMER tools (HMM marker presets: RNA-virus RdRp families, Poxviridae, Adenoviridae):")
    for tool in HMM_TOOLS:
        path = shutil.which(tool)
        if path:
            click.echo(f"  [OK]  {tool}  ({path})")
        else:
            click.echo(f"  [MISSING]  {tool}")
            ok = False

    click.echo("")

    # Global config
    click.echo("Configuration:")
    if not global_config.exists():
        click.echo(f"  [MISSING]  {global_config}")
        click.echo(f"             Run 'vfam_trees init' to generate a template.")
        ok = False
    else:
        click.echo(f"  [OK]  {global_config}")
        try:
            import yaml as _yaml
            with open(global_config) as f:
                cfg = _yaml.safe_load(f)
            email = cfg.get("ncbi", {}).get("email", "")
            api_key = cfg.get("ncbi", {}).get("api_key", "")
            if email:
                click.echo(f"  [OK]  ncbi.email = {email}")
            else:
                click.echo("  [MISSING]  ncbi.email is not set in global.yaml")
                ok = False
            if api_key:
                click.echo(f"  [OK]  ncbi.api_key is set")
            else:
                click.echo("  [--]  ncbi.api_key not set (rate-limited to 3 req/s)")
        except Exception as e:
            click.echo(f"  [ERROR]  Could not parse {global_config}: {e}")
            ok = False

    click.echo("")

    # NCBI connectivity
    click.echo("NCBI connectivity:")
    try:
        from Bio import Entrez
        Entrez.email = "test@test.com"
        handle = Entrez.esearch(db="taxonomy", term="Flaviviridae[Scientific Name]", retmax=1)
        result = Entrez.read(handle)
        handle.close()
        if result.get("IdList"):
            click.echo("  [OK]  NCBI Entrez reachable")
        else:
            click.echo("  [WARN]  NCBI reachable but returned no results (check connectivity)")
    except Exception as e:
        click.echo(f"  [FAIL]  Cannot reach NCBI Entrez: {e}")
        ok = False

    click.echo("")

    # Cache
    click.echo("Sequence cache:")
    try:
        import yaml as _yaml2
        with open(global_config) as f:
            _gcfg = _yaml2.safe_load(f)
        _cache_cfg = (_gcfg or {}).get("cache") or {}
        _cache_dir = _cache_cfg.get("dir") or None
        if _cache_dir:
            from pathlib import Path as _Path
            from .cache import SequenceCache as _SC
            _sc = _SC(_Path(_cache_dir))
            _st = _sc.stats()
            click.echo(f"  [OK]  {_sc.cache_dir}  ({_st['entries']} entries, {_st['size_mb']} MB)")
            _ttl = _cache_cfg.get("ttl_days")
            if _ttl:
                click.echo(f"  [OK]  ttl_days = {_ttl}")
            else:
                click.echo("  [--]  ttl_days not set (entries never expire)")
        else:
            click.echo("  [--]  cache.dir not set in global.yaml — caching disabled")
    except Exception:
        click.echo("  [--]  (could not read cache config)")

    click.echo("")
    if ok:
        click.echo("All checks passed.")
    else:
        click.echo("Some checks failed — see above.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# cache
# ---------------------------------------------------------------------------

@main.group()
def cache():
    """Manage the shared sequence download cache."""


def _load_cache(global_config: Path):
    """Return a SequenceCache instance from global config, or exit if unconfigured."""
    import yaml as _yaml
    from .cache import SequenceCache

    if not global_config.exists():
        click.echo(f"Global config not found: {global_config}", err=True)
        sys.exit(1)
    with open(global_config) as f:
        gcfg = _yaml.safe_load(f)
    cache_dir = ((gcfg or {}).get("cache") or {}).get("dir")
    if not cache_dir:
        click.echo("cache.dir is not set in global.yaml — caching is disabled.", err=True)
        sys.exit(1)
    return SequenceCache(Path(cache_dir).expanduser())


@cache.command("clear")
@click.argument("family", required=False)
@click.option(
    "--all", "clear_all",
    is_flag=True,
    default=False,
    help="Clear the entire cache (all families).",
)
@click.option(
    "--empty-only", "empty_only",
    is_flag=True,
    default=False,
    help=(
        "Remove only the negative-result (empty) sentinels, keeping downloaded "
        "sequence data.  Useful after a query fix when NCBI now returns results "
        "for species previously cached as empty."
    ),
)
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
)
@click.option(
    "--yes", "-y",
    is_flag=True,
    default=False,
    help="Skip confirmation prompt.",
)
def cache_clear(
    family: str | None,
    clear_all: bool,
    empty_only: bool,
    global_config: Path,
    yes: bool,
):
    """Clear cached sequences for a family or the entire cache.

    \b
    Examples:
      vfam_trees cache clear Hantaviridae
      vfam_trees cache clear Hantaviridae --empty-only
      vfam_trees cache clear --all
      vfam_trees cache clear --all --empty-only
      vfam_trees cache clear --all --yes
    """
    if not family and not clear_all:
        click.echo("Specify a FAMILY name or --all.", err=True)
        sys.exit(1)
    if family and clear_all:
        click.echo("Specify either a FAMILY name or --all, not both.", err=True)
        sys.exit(1)

    sc = _load_cache(global_config)

    if clear_all:
        if empty_only:
            if not yes:
                click.confirm(
                    "Delete all empty-result sentinels from the cache?", abort=True,
                )
            removed = sc.clear_empty_all()
            click.echo(f"Removed {removed} empty-result sentinel(s).")
        else:
            st = sc.stats()
            if not yes:
                click.confirm(
                    f"Delete all {st['entries']} cache entries ({st['size_mb']} MB) in {sc.cache_dir}?",
                    abort=True,
                )
            removed = sc.clear_all()
            click.echo(f"Removed {removed} cache entries.")
    else:
        if empty_only:
            if not yes:
                click.confirm(
                    f"Delete empty-result sentinels for {family}?", abort=True,
                )
            removed = sc.clear_empty_family(family)
            if removed:
                click.echo(f"Removed {removed} empty-result sentinel(s) for {family}.")
            else:
                click.echo(f"No empty-result sentinels found for {family}.")
        else:
            if not yes:
                click.confirm(f"Delete all cached sequences for {family}?", abort=True)
            removed = sc.clear_family(family)
            if removed:
                click.echo(f"Removed {removed} cache entries for {family}.")
            else:
                click.echo(f"No cache entries found for {family}.")


@cache.command("stats")
@click.option(
    "--global-config", "-g",
    default=DEFAULT_GLOBAL_CFG,
    show_default=True,
    type=click.Path(path_type=Path),
)
def cache_stats(global_config: Path):
    """Show cache size and entry count.

    \b
    Examples:
      vfam_trees cache stats
    """
    sc = _load_cache(global_config)
    st = sc.stats()
    click.echo(f"Cache directory : {st['cache_dir']}")
    click.echo(f"Entries         : {st['entries']}")
    click.echo(f"Total size      : {st['size_mb']} MB")


# ---------------------------------------------------------------------------
# overview
# ---------------------------------------------------------------------------

@main.command("overview")
@click.option(
    "--output-dir", "-o",
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    type=click.Path(path_type=Path),
    help="Results directory to scan for tree_100 files.",
)
@click.option(
    "--output", "-O",
    default=None,
    type=click.Path(path_type=Path),
    help="Output PNG path (default: <output-dir>/overview_tree_100.png).",
)
def overview(output_dir: Path, output: Path | None):
    """Generate a thumbnail grid PNG of all tree_100 trees.

    Scans the results directory for completed tree_100.nwk files,
    draws each as a small topology-only thumbnail (no leaf labels),
    and arranges them in a grid with the family name below each.

    This is run automatically at the end of 'vfam_trees run'. Use
    this command to regenerate the overview at any time.

    \b
    Examples:
      vfam_trees overview
      vfam_trees overview -o /data/results
      vfam_trees overview -O /data/results/my_overview.png
    """
    from .report import generate_overview_png
    out = output if output is not None else output_dir / "overview_tree_100.png"
    generate_overview_png(output_dir, Path(out))
    click.echo(f"Done: {out}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _read_family_list(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    return [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]


def _find_status_file(family: str, output_dir: Path) -> Path | None:
    """Return the .status.json path for a family, or None.

    Resolves the family's directory via _find_family_dir (which honours the
    .family marker, so manual.name-renamed directories are found too)."""
    family_dir = _find_family_dir(family, output_dir)
    if family_dir is not None:
        status = family_dir / ".status.json"
        if status.exists():
            return status
    return None


def _clear_status(families: list[str], output_dir: Path) -> None:
    for family in families:
        sentinel = output_dir / f".done_{family}"
        if sentinel.exists():
            sentinel.unlink()
            click.echo(f"Cleared status for {family}")

    # Drop any pre-existing rows for these families from summary.tsv and
    # status.tsv — both files are append-only at the per-family writer, so
    # without this --force re-runs would leave duplicate rows.
    families_set = set(families)
    for tsv_name in ("summary.tsv", "status.tsv"):
        path = output_dir / tsv_name
        n_dropped = _filter_tsv_rows(path, families_set)
        if n_dropped:
            click.echo(f"Dropped {n_dropped} stale row(s) from {tsv_name}")


def _filter_tsv_rows(path: Path, families_to_drop: set[str]) -> int:
    """Rewrite a TSV in place, removing rows whose 'family' is in the set.

    Returns the number of rows dropped.  No-op (and returns 0) if the file is
    missing, empty, or has no 'family' column.
    """
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "family" not in reader.fieldnames:
            return 0
        fieldnames = list(reader.fieldnames)
        kept_rows = []
        n_dropped = 0
        for row in reader:
            if row.get("family", "") in families_to_drop:
                n_dropped += 1
            else:
                kept_rows.append(row)
    if n_dropped == 0:
        return 0
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        for row in kept_rows:
            writer.writerow(row)
    return n_dropped


def _find_family_dir(family: str, output_dir: Path) -> Path | None:
    """Return the output directory for a family, or None.

    When ``manual.name`` is set the directory is named after the (sanitized)
    display name, so it cannot be located by the family name alone.  Each
    directory carries a ``.family`` marker recording the biological family that
    owns it; that is matched first.  If several match (e.g. a renamed re-run
    left a stale directory) the most recently modified one wins.  The name-glob
    fallbacks keep older result sets (written before the marker existed)
    discoverable.
    """
    if output_dir.is_dir():
        owned = [
            d for d in output_dir.iterdir()
            if d.is_dir() and _owner_family(d) == family
        ]
        if owned:
            return max(owned, key=lambda d: d.stat().st_mtime)
    plain = output_dir / family
    if plain.is_dir():
        return plain
    for d in sorted(output_dir.glob(f"{family}_*")):
        if d.is_dir():
            return d
    return None


def _owner_family(family_dir: Path) -> str | None:
    """Return the biological family recorded in the directory's .family marker."""
    marker = family_dir / ".family"
    if marker.is_file():
        try:
            return marker.read_text().strip() or None
        except OSError:
            return None
    return None


def _detect_stage(family_dir: Path) -> str:
    """Infer the current pipeline stage from work-dir checkpoint files."""
    work = family_dir / "_work"
    if not work.is_dir():
        return "initializing"
    checkpoints = [
        (work / "100" / ".tree_done", "annotating / writing outputs"),
        (work / "100" / ".msa_done",  "tree inference (tree_100)"),
        (work / "500" / ".tree_done", "MSA + clustering (tree_100)"),
        (work / "500" / ".msa_done",  "tree inference (tree_500)"),
    ]
    for path, stage in checkpoints:
        if path.exists():
            return stage
    if (work / "species_list.json").exists():
        return "downloading / QC"
    return "initializing"


def _print_dry_run(
    family_list: list[str],
    global_config: Path,
    configs_dir: Path,
    output_dir: Path,
    threads: int,
) -> None:
    """Print a per-family config preview table without running anything."""
    from .config import load_global_config, load_family_config

    click.echo(f"DRY RUN — vfam_trees pipeline preview")
    click.echo(f"Families file: {len(family_list)} families  |  Output: {output_dir}  |  Threads/job: {threads}")
    click.echo("")

    global_cfg = load_global_config(global_config)

    header = (
        f"{'Family':<22}  {'Config':<8}  {'Type':<12}  {'Region/segment':<18}  "
        f"{'max/sp':>6}  {'Tree-500':<18}  {'Tree-100':<18}"
    )
    click.echo(header)
    click.echo("─" * len(header))

    for family in family_list:
        cfg, auto = load_family_config(family, configs_dir, global_cfg)
        cfg_tag = "auto" if auto else "exists"
        seq_type = cfg["sequence"]["type"]
        region = cfg["sequence"].get("segment") or cfg["sequence"].get("region", "")
        max_sp = cfg["download"]["max_per_species"]

        t500 = cfg.get("tree_500", {})
        t100 = cfg.get("tree_100", {})
        model_key = "model_aa" if seq_type == "protein" else "model_nuc"
        tree500_str = f"{t500.get('tool','?')} {t500.get(model_key,'?')}"
        tree100_str = f"{t100.get('tool','?')} {t100.get(model_key,'?')}"

        # Mark already-completed families
        family_dir = _find_family_dir(family, output_dir)
        status_file = _find_status_file(family, output_dir)
        if status_file is not None:
            try:
                st = json.loads(status_file.read_text()).get("status", "")
                cfg_tag = f"{cfg_tag} [{st}]"
            except Exception:
                pass

        click.echo(
            f"{family:<22}  {cfg_tag:<8}  {seq_type:<12}  {region:<18}  "
            f"{max_sp:>6}  {tree500_str:<18}  {tree100_str:<18}"
        )


def _yaml_encode(value) -> str:
    """Encode a Python value for passing via snakemake --config key=value."""
    if isinstance(value, list):
        return json.dumps(value)
    return str(value)
