"""vfam_trees command-line interface."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import click
import yaml

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"
DEFAULT_GLOBAL_CFG = Path("config/global.yaml")
DEFAULT_CONFIGS_DIR = Path("configs")
DEFAULT_OUTPUT_DIR = Path("results")


@click.group()
@click.version_option()
def main():
    """vfam_trees — phylogenetic trees for viral families.

    Builds maximum-likelihood phylogenetic trees for one or more viral
    families using sequences downloaded from NCBI. For each family, two
    trees are produced: a broad 500-sequence tree (FastTree) and a more
    collapsed 100-sequence tree (IQ-TREE).

    \b
    Typical workflow:
      1. Edit config/global.yaml  (set NCBI email and API key)
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
    discover species (NCBI taxonomy) → download per species →
    quality filter → cluster per species → proportional merge →
    align (MAFFT) → infer trees (FastTree / IQ-TREE) →
    annotate internal nodes → write Newick + PhyloXML.

    Already-processed families are skipped unless --force is used.
    Per-family configs are auto-generated in configs/ if missing.

    \b
    Examples:
      vfam_trees run -f families.txt
      vfam_trees run -f families.txt -j 4 -t 4
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

    cmd = [
        "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--cores", str(cores),
        "--config", *[f"{k}={_yaml_encode(v)}" for k, v in snakemake_cfg.items()],
        "--rerun-incomplete",
        "--printshellcmds",
    ]

    click.echo(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)

    if result.returncode == 0:
        from .report import generate_overview_png
        overview_path = output_dir / "overview_tree_100.png"
        click.echo(f"Generating overview PNG: {overview_path}")
        generate_overview_png(output_dir, overview_path)

    sys.exit(result.returncode)


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

    Reports one of: pending, done, skipped (with reason), or
    error reading status.

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
    from .config import load_global_config, load_family_config

    family_list = _read_family_list(families)
    global_cfg = load_global_config(global_config)

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

output_dir: results
log_level: INFO         # DEBUG, INFO, WARNING

# Optional shared sequence download cache.
# Sequences are stored keyed by (taxid, db, region, segment, max_per_species)
# and reused across runs and families on the same machine — or on a shared
# filesystem across machines.  Remove or set dir to null to disable.
cache:
  dir: ~/.vfam_cache    # path to cache root; ~ is expanded automatically
  ttl_days: 90          # null = never expire; set to re-download after N days

# Default values applied when per-family config is auto-generated.
# Override any of these here to change the default for all families.
defaults:
  download:
    max_per_species: 200

  sequence:
    type: nucleotide    # nucleotide or protein
    region: whole_genome
    segment: null

  quality:
    min_length: null    # null = auto (50% of median sequence length)
    max_ambiguous: 0.01
    exclude_organisms:
      - synthetic construct
      - metagenome
      - MAG
      - uncultured
      - unverified
      - vector

  clustering:
    tool: mmseqs2
    threshold_min: 0.70
    threshold_max: 0.99
    max_reps_500: 20
    max_reps_100: 5

  targets:
    max_500: 500
    max_100: 100

  msa_500:
    tool: mafft
    options_nuc: "--6merpair --retree 1"
    options_aa: "--6merpair --retree 1"

  msa_100:
    tool: mafft
    options_nuc: "--retree 2"
    options_aa: "--auto"

  tree_500:
    tool: fasttree
    options: ""
    model_nuc: GTR+G
    model_aa: LG+G

  tree_100:
    tool: iqtree
    options: "--fast"
    model_nuc: GTR+G
    model_aa: TEST
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
def init(global_config: Path, force: bool):
    """Generate a template global config file.

    Creates config/global.yaml (or the path given by --global-config)
    with all default values and comments explaining each setting.
    Edit the file — especially ncbi.email — before running the pipeline.

    \b
    Examples:
      vfam_trees init
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
    global_config.write_text(GLOBAL_CONFIG_TEMPLATE)
    click.echo(f"Created: {global_config}")
    click.echo("Edit the file and set ncbi.email before running the pipeline.")


# ---------------------------------------------------------------------------
# test
# ---------------------------------------------------------------------------

REQUIRED_TOOLS = ["mafft", "iqtree2", "mmseqs", "FastTree"]
OPTIONAL_TOOLS = ["cd-hit", "cd-hit-est"]
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
      - Required external tools (mafft, iqtree2, mmseqs, FastTree)
      - Optional tools (cd-hit, cd-hit-est — only needed as fallback)
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
    click.echo("Optional external tools (fallback only):")
    for tool in OPTIONAL_TOOLS:
        path = shutil.which(tool)
        if path:
            click.echo(f"  [OK]  {tool}  ({path})")
        else:
            click.echo(f"  [--]  {tool}  (not found — only needed if clustering.tool: cdhit)")

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
def cache_clear(family: str | None, clear_all: bool, global_config: Path, yes: bool):
    """Clear cached sequences for a family or the entire cache.

    \b
    Examples:
      vfam_trees cache clear Asfarviridae
      vfam_trees cache clear --all
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
        st = sc.stats()
        if not yes:
            click.confirm(
                f"Delete all {st['entries']} cache entries ({st['size_mb']} MB) in {sc.cache_dir}?",
                abort=True,
            )
        removed = sc.clear_all()
        click.echo(f"Removed {removed} cache entries.")
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
    """Return the .status.json path for a family, searching both the plain
    name and the taxid-suffixed name (e.g. Asfarviridae_137992)."""
    plain = output_dir / family / ".status.json"
    if plain.exists():
        return plain
    for candidate in sorted(output_dir.glob(f"{family}_*/.status.json")):
        return candidate
    return None


def _clear_status(families: list[str], output_dir: Path) -> None:
    for family in families:
        sentinel = output_dir / f".done_{family}"
        if sentinel.exists():
            sentinel.unlink()
            click.echo(f"Cleared status for {family}")


def _find_family_dir(family: str, output_dir: Path) -> Path | None:
    """Return the output directory for a family (taxid-suffixed or plain), or None."""
    plain = output_dir / family
    if plain.is_dir():
        return plain
    for d in sorted(output_dir.glob(f"{family}_*")):
        if d.is_dir():
            return d
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
