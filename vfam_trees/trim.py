"""Alignment column trimming — trimAl wrapper."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .logger import get_logger

log = get_logger(__name__)

# Keep columns with at most 50 % gaps.  Replaced -automated1 in v1.6.0:
# on alignments with many short PCR fragments -automated1 switches to its
# strict mode and keeps few columns (Retroviridae RT tree_500: 58 of 657),
# which cost genus resolution.  -gt 0.5 was never worse on 22 protein and
# 20 nucleotide alignments and better on 6 (CHANGELOG 1.6.0).
DEFAULT_TRIM_OPTIONS = "-gt 0.5"


def get_trimal_version() -> str:
    """Return the trimAl version string, or 'unknown' if it cannot be determined."""
    try:
        result = subprocess.run(
            ["trimal", "--version"], capture_output=True, text=True
        )
        output = (result.stdout or result.stderr).strip()
        m = re.search(r"(\d+\.\d+\S*)", output)
        if m:
            return m.group(1)
        return output.split()[0] if output else "unknown"
    except Exception:
        return "unknown"


def run_trim(
    input_fasta: Path,
    output_fasta: Path,
    tool: str = "trimal",
    options: str = DEFAULT_TRIM_OPTIONS,
) -> None:
    """Trim an MSA and write the column-filtered alignment to *output_fasta*.

    Args:
        input_fasta: MAFFT output alignment (short IDs)
        output_fasta: path to write the trimmed alignment
        tool: trimming tool ('trimal' supported)
        options: tool-specific options string (default DEFAULT_TRIM_OPTIONS)
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    tool_norm = tool.lower().strip()
    if tool_norm != "trimal":
        raise ValueError(
            f"Unsupported MSA trimming tool: {tool!r}. Supported: 'trimal'."
        )
    _run_trimal(input_fasta, output_fasta, options)
    log.debug("Trim complete: %s", output_fasta)


def _run_trimal(input_fasta: Path, output_fasta: Path, options: str) -> None:
    cmd = ["trimal", "-in", str(input_fasta), "-out", str(output_fasta)]
    cmd += options.split()

    log.info("Running trimAl (%s)", options)
    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("trimAl stderr:\n%s", (result.stderr or "")[-4000:])
        log.error("trimAl stdout:\n%s", (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"trimAl failed with exit code {result.returncode}. If trimAl "
            "was upgraded recently, its CLI flags may have changed — see "
            "stderr above."
        )

    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        log.error("trimAl stderr:\n%s", (result.stderr or "")[-4000:])
        if "empty" in (result.stderr or "").lower():
            raise RuntimeError(
                f"trimAl ({options}) kept no column of {input_fasta.name}: "
                "every column is too gappy for that threshold.  Loosen "
                "msa_trim.options (e.g. -gt 0.3) or set msa_trim.enabled: "
                "false in this family's config."
            )
        raise FileNotFoundError(
            f"trimAl exited 0 but produced no output at {output_fasta}. "
            "This may indicate a change in trimAl's output behavior — see "
            "stderr above."
        )
    log.info("trimAl complete")
