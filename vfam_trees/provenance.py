"""Run provenance: which version, tools and settings produced a results folder.

``vfam_trees run`` writes one run record per invocation into
``<output>/run_records/`` (``start_run_record`` / ``finish_run_record``), and
the pipeline writes each family's merged settings, with a header from
``family_config_header``, into the family folder.  Together they say what a
results folder was made with, without having to reconstruct a past version's
defaults.
"""
from __future__ import annotations

import platform
import re
import shlex
import shutil
import subprocess
from datetime import datetime
from importlib import metadata
from pathlib import Path

import yaml

from . import __version__
from .config import effective_global_config
from .logger import get_logger

log = get_logger(__name__)

RUN_RECORDS_DIR = "run_records"

#: Python distributions whose versions go into the run record.
_PACKAGES = ("biopython", "snakemake", "numpy", "matplotlib", "click", "PyYAML", "requests")

_RUN_RECORD_HEADER = """\
# vfam_trees run record.  Written by `vfam_trees run` when the run started;
# status, finished and exit_code are filled in when it ends (status stays
# "running" if the process was killed).
# global_config: every global setting the run used, i.e. the built-in
# defaults overlaid with global.yaml (ncbi.api_key is not recorded).  Family
# presets and family files are applied on top of `defaults` per family; the
# settings each family ran with are in <family folder>/<prefix>.yaml.

"""


def _probe(cmd: list[str], pattern: str) -> str:
    """First match group of *pattern* in a tool's output, or "unknown"."""
    if not shutil.which(cmd[0]):
        return "not found"
    try:
        res = subprocess.run(
            cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=30,
        )
        m = re.search(pattern, f"{res.stdout or ''}\n{res.stderr or ''}")
        return m.group(1) if m else "unknown"
    except Exception:
        return "unknown"


def tool_versions() -> dict[str, str]:
    """Versions of the external tools on PATH ("not found" when absent)."""
    from .msa import get_mafft_version
    from .tree import get_tree_tool_version
    from .trim import get_trimal_version

    def known(tool: str, probe) -> str:
        if not shutil.which(tool):
            return "not found"
        try:
            return probe() or "unknown"
        except Exception:
            return "unknown"

    return {
        "mafft": known("mafft", get_mafft_version),
        "trimal": known("trimal", get_trimal_version),
        "iqtree2": known("iqtree2", lambda: get_tree_tool_version("iqtree")),
        "FastTree": known("FastTree", lambda: get_tree_tool_version("fasttree")),
        "mmseqs": _probe(["mmseqs", "version"], r"^\s*(\S+)"),
        "cd-hit": _probe(["cd-hit", "-h"], r"CD-HIT version\s+(\S+)"),
        "hmmscan": _probe(["hmmscan", "-h"], r"HMMER\s+(\S+)"),
    }


def package_versions() -> dict[str, str]:
    """Python version and the versions of the main Python packages."""
    out = {"python": platform.python_version()}
    for dist in _PACKAGES:
        try:
            out[dist] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            out[dist] = "not installed"
    return out


def _platform() -> str:
    try:
        return platform.platform()
    except Exception:  # it runs `uname`; never let record-keeping stop a run
        return "unknown"


def _write_record(path: Path, record: dict) -> None:
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w") as f:
        f.write(_RUN_RECORD_HEADER)
        yaml.safe_dump(record, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    tmp.replace(path)


def start_run_record(
    output_dir: Path,
    *,
    argv: list[str],
    families_file: Path,
    families: list[str],
    options: dict,
    paths: dict,
    global_cfg: dict,
) -> Path:
    """Write ``<output_dir>/run_records/<start time>.yaml`` and return its path."""
    now = datetime.now()
    record = {
        "vfam_trees_version": __version__,
        "status": "running",
        "started": now.isoformat(timespec="seconds"),
        "finished": None,
        "exit_code": None,
        "command": shlex.join(argv),
        "working_directory": str(Path.cwd()),
        "families_file": str(Path(families_file).resolve()),
        "families": list(families),
        "options": dict(options),
        "paths": dict(paths),
        "platform": _platform(),
        "python_packages": package_versions(),
        "tools": tool_versions(),
        "global_config_unknown_keys": list(global_cfg.get("_unknown_keys") or []),
        "global_config": effective_global_config(global_cfg),
    }
    folder = Path(output_dir) / RUN_RECORDS_DIR
    folder.mkdir(parents=True, exist_ok=True)
    stem = now.strftime("%Y-%m-%dT%H-%M-%S")
    path = folder / f"{stem}.yaml"
    n = 2
    while path.exists():  # two runs started in the same second
        path = folder / f"{stem}_{n}.yaml"
        n += 1
    _write_record(path, record)
    return path


def finish_run_record(path: Path, exit_code: int | None) -> None:
    """Fill in status, finished and exit_code (None: interrupted)."""
    try:
        with open(path) as f:
            record = yaml.safe_load(f) or {}
        record["status"] = (
            "interrupted" if exit_code is None
            else "completed" if exit_code == 0 else "failed"
        )
        record["finished"] = datetime.now().isoformat(timespec="seconds")
        record["exit_code"] = exit_code
        _write_record(Path(path), record)
    except (OSError, yaml.YAMLError) as exc:
        log.warning("Could not complete the run record %s: %s", path, exc)


def family_config_header(
    *,
    family: str,
    started: datetime,
    threads: int,
    global_config_path: Path,
    family_config_path: Path,
    auto_generated: bool,
    unknown_keys: list[str],
) -> str:
    """Comment block at the top of ``<family folder>/<prefix>.yaml``."""
    source = (
        "auto-generated by this run" if auto_generated else "existing file"
    )
    lines = [
        f"# Settings vfam_trees {__version__} used for {family} "
        f"(started {started.isoformat(sep=' ', timespec='seconds')}).",
        "# Merged, lowest to highest priority: built-in defaults, the",
        f"# `defaults:` of {global_config_path},",
        "# the family preset (RdRp, DNA marker, segment), and",
        f"# {family_config_path} ({source}).",
        f"# Threads: {threads} (also hmmscan's, when hmm.threads is null).",
    ]
    if unknown_keys:
        lines.append(
            "# Not read by this version (no effect): " + ", ".join(unknown_keys) + "."
        )
    lines += [
        "# To repeat this analysis with the same version, use this file as",
        f"# the family's config ({family_config_path.name}).",
        "",
        "",
    ]
    return "\n".join(lines)
