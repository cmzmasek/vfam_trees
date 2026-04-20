"""Phylogenetic tree inference — FastTree and IQ-TREE wrappers."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .logger import get_logger

log = get_logger(__name__)


def get_tree_tool_version(tool: str) -> str:
    """Return the version string for fasttree or iqtree, or 'unknown'."""
    tool_norm = tool.lower().replace("-", "").replace("_", "")
    try:
        if tool_norm == "fasttree":
            # FastTree prints version to stderr when run with no args
            result = subprocess.run(
                ["FastTree"], capture_output=True, text=True
            )
            output = result.stderr or ""
            m = re.search(r"FastTree\s+version\s+(\S+)", output, re.IGNORECASE)
            if not m:
                m = re.search(r"(\d+\.\d+\.\d+)", output)
            return m.group(1) if m else "unknown"
        elif tool_norm in ("iqtree", "iqtree2"):
            result = subprocess.run(
                ["iqtree2", "--version"], capture_output=True, text=True
            )
            output = result.stdout or result.stderr or ""
            m = re.search(r"version\s+(\d+\.\d+\S*)", output, re.IGNORECASE)
            return m.group(1) if m else "unknown"
    except Exception:
        pass
    return "unknown"


def parse_iqtree_best_model(log_path: Path) -> str | None:
    """Return the best-fit model string from an IQ-TREE `.iqtree` log, or None.

    Prefers the ModelFinder line ("Best-fit model according to BIC/AIC/AICc")
    and falls back to the "Model of substitution" line that IQ-TREE always
    writes.
    """
    if not log_path.exists():
        return None
    try:
        text = log_path.read_text()
    except Exception:
        return None
    m = re.search(r"^Best-fit model according to \S+:\s*(\S+)", text, re.MULTILINE)
    if m:
        return m.group(1).strip()
    m = re.search(r"^Model of substitution:\s*(\S+)", text, re.MULTILINE)
    if m:
        return m.group(1).strip()
    return None


def is_model_finder_spec(model: str) -> bool:
    """Return True if `model` asks IQ-TREE's ModelFinder to pick the model."""
    if not model:
        return True
    root = model.split("+", 1)[0].strip().upper()
    return root in {"TEST", "TESTONLY", "TESTNEW", "TESTNEWONLY", "MF", "MFP"}


def validate_newick(nwk_path: Path) -> None:
    """Validate that a Newick file is non-empty and parseable.

    Raises:
        FileNotFoundError: if the file does not exist or is empty
        ValueError: if the file cannot be parsed as a Newick tree
    """
    if not nwk_path.exists() or nwk_path.stat().st_size == 0:
        raise FileNotFoundError(f"Newick file is missing or empty: {nwk_path}")
    try:
        from Bio import Phylo
        trees = list(Phylo.parse(str(nwk_path), "newick"))
        if not trees:
            raise ValueError(f"No trees found in Newick file: {nwk_path}")
        log.debug("Newick validated: %d tree(s) parsed from %s", len(trees), nwk_path)
    except Exception as e:
        raise ValueError(f"Newick file failed to parse ({nwk_path}): {e}") from e


def run_tree(
    alignment_fasta: Path,
    output_dir: Path,
    prefix: str,
    tool: str,
    seq_type: str,
    model_nuc: str = "GTR+G",
    model_aa: str = "LG+G",
    options: str = "",
    threads: int = 1,
) -> Path:
    """Infer a phylogenetic tree.

    Args:
        alignment_fasta: aligned FASTA (short IDs)
        output_dir: directory for output files
        prefix: output file prefix (e.g. 'tree_500')
        tool: 'fasttree' or 'iqtree'
        seq_type: 'nucleotide' or 'protein'
        model_nuc: substitution model for nucleotide data
        model_aa: substitution model for amino acid data
        options: additional tool-specific options string
        threads: number of CPU threads (IQ-TREE only)

    Returns:
        Path to the output Newick (.nwk) file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_norm = tool.lower().replace("-", "").replace("_", "")

    if tool_norm == "fasttree":
        return _run_fasttree(alignment_fasta, output_dir, prefix, seq_type,
                             model_nuc, model_aa, options)
    elif tool_norm in ("iqtree", "iqtree2"):
        return _run_iqtree(alignment_fasta, output_dir, prefix, seq_type,
                           model_nuc, model_aa, options, threads)
    else:
        raise ValueError(f"Unsupported tree tool: {tool}. Supported: fasttree, iqtree")


# ---------------------------------------------------------------------------
# FastTree
# ---------------------------------------------------------------------------

def _fasttree_nuc_flags(model: str) -> list[str]:
    """Map a nucleotide model string to FastTree command-line flags.

    FastTree supports GTR (with ``-gtr``) or its default JC69. Any other
    substitution model falls back to GTR with a warning. ``+G`` / ``+GAMMA``
    suffixes enable the discrete-gamma rate model via ``-gamma``.
    """
    flags: list[str] = []
    root = model.split("+", 1)[0].strip().upper() if model else ""
    if root == "GTR" or not root:
        flags.append("-gtr")
    elif root in ("JC", "JC69"):
        pass  # FastTree default
    else:
        log.warning(
            "FastTree does not support nucleotide model %r — using GTR instead.",
            model,
        )
        flags.append("-gtr")
    upper = (model or "").upper()
    if "+G" in upper or "+GAMMA" in upper:
        flags.append("-gamma")
    return flags


def _fasttree_aa_flags(model: str) -> list[str]:
    """Map a protein model string to FastTree command-line flags.

    FastTree supports WAG (``-wag``), LG (``-lg``), or its default JTT.
    Unknown models fall back to WAG with a warning. ``+G`` / ``+GAMMA``
    suffixes enable discrete gamma via ``-gamma``.
    """
    flags: list[str] = []
    root = model.split("+", 1)[0].strip().upper() if model else ""
    if root == "LG":
        flags.append("-lg")
    elif root == "WAG":
        flags.append("-wag")
    elif root in ("", "JTT"):
        pass  # FastTree default
    else:
        log.warning(
            "FastTree does not support amino-acid model %r — using WAG instead.",
            model,
        )
        flags.append("-wag")
    upper = (model or "").upper()
    if "+G" in upper or "+GAMMA" in upper:
        flags.append("-gamma")
    return flags


def _run_fasttree(
    alignment_fasta: Path,
    output_dir: Path,
    prefix: str,
    seq_type: str,
    model_nuc: str,
    model_aa: str,
    options: str,
) -> Path:
    out_nwk = output_dir / f"{prefix}.nwk"

    cmd = ["FastTree"]
    if seq_type == "nucleotide":
        cmd += ["-nt"] + _fasttree_nuc_flags(model_nuc)
    else:
        cmd += _fasttree_aa_flags(model_aa)

    if options:
        cmd += options.split()

    cmd += [str(alignment_fasta)]

    log.debug("Running: %s > %s", " ".join(cmd), out_nwk)
    with open(out_nwk, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.error("FastTree stderr:\n%s", (result.stderr or "")[-4000:])
        raise RuntimeError(
            f"FastTree failed with exit code {result.returncode}. If "
            "FastTree was upgraded recently, its CLI flags may have changed "
            "— see stderr above."
        )

    if not out_nwk.exists() or out_nwk.stat().st_size == 0:
        log.error("FastTree stderr:\n%s", (result.stderr or "")[-4000:])
        raise FileNotFoundError(
            f"FastTree exited 0 but produced no Newick at {out_nwk}. "
            "This may indicate a change in FastTree's output behavior — see "
            "stderr above."
        )

    log.info("FastTree complete: %s", out_nwk)
    return out_nwk


# ---------------------------------------------------------------------------
# IQ-TREE
# ---------------------------------------------------------------------------

def _run_iqtree(
    alignment_fasta: Path,
    output_dir: Path,
    prefix: str,
    seq_type: str,
    model_nuc: str,
    model_aa: str,
    options: str,
    threads: int,
) -> Path:
    out_prefix = output_dir / prefix
    model = model_nuc if seq_type == "nucleotide" else model_aa

    # Ensure bootstrap support is always computed unless the user has explicitly
    # specified a bootstrap flag (-B ultrafast or -b standard).
    # Note: IQ-TREE's UFBoot (-B) is incompatible with --fast, so --fast is
    # removed automatically when bootstrap is injected.
    options_parts = options.split() if options else []
    if "-alrt" not in options_parts and "-B" not in options_parts and "-b" not in options_parts:
        options_parts = ["-alrt", "1000"] + options_parts
        log.debug("Auto-added -alrt 1000 for IQ-TREE SH-aLRT support")

    cmd = [
        "iqtree2",
        "-s", str(alignment_fasta),
        "--prefix", str(out_prefix),
        "-m", model,
        "-T", str(threads),
        "--redo",
    ] + options_parts

    log.debug("Running: %s", " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("IQ-TREE stderr:\n%s", (result.stderr or "")[-4000:])
        log.error("IQ-TREE stdout:\n%s", (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"IQ-TREE failed with exit code {result.returncode}. If "
            "IQ-TREE was upgraded recently, its CLI flags or model names "
            "may have changed — see stderr above."
        )

    treefile = Path(str(out_prefix) + ".treefile")
    if not treefile.exists():
        log.error("IQ-TREE stderr:\n%s", (result.stderr or "")[-4000:])
        raise FileNotFoundError(
            f"IQ-TREE exited 0 but did not produce the expected tree file "
            f"{treefile}. This may indicate a change in IQ-TREE's output "
            "conventions."
        )

    # Copy to a consistently named .nwk file
    out_nwk = output_dir / f"{prefix}.nwk"
    out_nwk.write_text(treefile.read_text())

    log.info("IQ-TREE complete: %s", out_nwk)
    return out_nwk
