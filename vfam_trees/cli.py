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
        status_file = output_dir / family / ".status.json"
        if status_file.exists():
            try:
                state = json.loads(status_file.read_text())
                st = state.get("status", "unknown")
                reason = state.get("reason", "")
            except Exception:
                st, reason = "error reading status", ""
        else:
            st, reason = "pending", ""
        rows.append((family, st, reason))

    w = max(len(r[0]) for r in rows) + 2
    click.echo(f"{'Family':<{w}}  {'Status':<12}  Notes")
    click.echo("-" * 60)
    for family, st, reason in rows:
        click.echo(f"{family:<{w}}  {st:<12}  {reason}")


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
def init_configs(families: Path, global_config: Path, configs_dir: Path):
    """Generate default per-family configs without running the pipeline.

    Creates configs/{Family}.yaml for any family that does not already
    have one. Existing configs are left untouched. Run this before the
    pipeline to review and tune parameters (clustering thresholds, models,
    max sequences per species, etc.) before committing to a full run.

    \b
    Examples:
      vfam_trees init-configs -f families.txt
      vfam_trees init-configs -f families.txt -c /data/configs
    """
    from .config import load_global_config, load_family_config

    family_list = _read_family_list(families)
    global_cfg = load_global_config(global_config)

    for family in family_list:
        cfg, auto = load_family_config(family, configs_dir, global_cfg)
        if auto:
            click.echo(f"Generated: {configs_dir / family}.yaml")
        else:
            click.echo(f"Exists:    {configs_dir / family}.yaml")


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
    if ok:
        click.echo("All checks passed.")
    else:
        click.echo("Some checks failed — see above.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _read_family_list(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    return [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]


def _clear_status(families: list[str], output_dir: Path) -> None:
    for family in families:
        status_file = output_dir / family / ".status.json"
        if status_file.exists():
            status_file.unlink()
            click.echo(f"Cleared status for {family}")


def _yaml_encode(value) -> str:
    """Encode a Python value for passing via snakemake --config key=value."""
    if isinstance(value, list):
        return json.dumps(value)
    return str(value)
