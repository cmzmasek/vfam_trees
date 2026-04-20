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


def validate_msa(fasta_path: Path) -> None:
    """Validate that a FASTA MSA file is non-empty and all sequences are the same length.

    Raises:
        FileNotFoundError: if the file does not exist or is empty
        ValueError: if fewer than 2 sequences, or sequences differ in length
    """
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        raise FileNotFoundError(f"MSA output is missing or empty: {fasta_path}")
    from Bio import SeqIO
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) < 2:
        raise ValueError(f"MSA contains fewer than 2 sequences: {fasta_path}")
    lengths = {len(r.seq) for r in records}
    if len(lengths) > 1:
        raise ValueError(
            f"MSA sequences have inconsistent lengths {lengths} — "
            f"alignment may be corrupted: {fasta_path}"
        )
    log.debug("MSA validated: %d sequences, alignment length %d", len(records), lengths.pop())


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
        log.error("MAFFT stderr:\n%s", (result.stderr or "")[-4000:])
        log.error("MAFFT stdout:\n%s", (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"MAFFT failed with exit code {result.returncode}. If MAFFT was "
            "upgraded recently, its CLI flags may have changed — see stderr "
            "above."
        )

    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        log.error("MAFFT stderr:\n%s", (result.stderr or "")[-4000:])
        raise FileNotFoundError(
            f"MAFFT exited 0 but produced no alignment at {output_fasta}. "
            "This may indicate a change in MAFFT's output behavior — see "
            "stderr above."
        )

    log.debug("MAFFT stdout: %s", result.stdout[:500] if result.stdout else "(none)")
