"""MSA inference using MAFFT (or other configured tool)."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .logger import get_logger

log = get_logger(__name__)


def get_mafft_version() -> str:
    """Return the MAFFT version string, or 'unknown' if it cannot be determined."""
    try:
        result = subprocess.run(
            ["mafft", "--version"], capture_output=True, text=True
        )
        # MAFFT prints version to stderr: "v7.520 (2023/Apr/16)"
        output = (result.stderr or result.stdout).strip()
        m = re.search(r"v?(\d+\.\d+\S*)", output)
        return m.group(1) if m else output.split()[0] if output else "unknown"
    except Exception:
        return "unknown"


def run_msa(
    input_fasta: Path,
    output_fasta: Path,
    tool: str = "mafft",
    options: str = "--auto",
    threads: int = 1,
) -> None:
    """Run multiple sequence alignment.

    Args:
        input_fasta: unaligned sequences (short IDs)
        output_fasta: output aligned FASTA
        tool: alignment tool name ('mafft' supported)
        options: tool-specific options string
        threads: number of CPU threads
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if tool.lower() == "mafft":
        _run_mafft(input_fasta, output_fasta, options, threads)
    else:
        raise ValueError(f"Unsupported MSA tool: {tool}. Supported: mafft")

    log.info("MSA complete: %s", output_fasta)


def _run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    options: str,
    threads: int,
) -> None:
    cmd = ["mafft"]
    cmd += options.split()
    cmd += ["--thread", str(threads)]
    cmd += ["--out", str(output_fasta)]
    cmd += [str(input_fasta)]

    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("MAFFT stderr:\n%s", result.stderr)
        raise RuntimeError(f"MAFFT failed with exit code {result.returncode}")

    log.debug("MAFFT stdout: %s", result.stdout[:500] if result.stdout else "(none)")
