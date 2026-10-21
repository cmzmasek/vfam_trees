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
            # stdin closed: with no arguments FastTree otherwise waits to read
            # an alignment from it, which hangs an interactive run.
            result = subprocess.run(
                ["FastTree"], capture_output=True, text=True,
                stdin=subprocess.DEVNULL,
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


# ---------------------------------------------------------------------------
# Branch-support flags — the single classifier every decision below uses
# ---------------------------------------------------------------------------

# IQ-TREE support measures, keyed by name: (option spellings, takes a value).
# Both the IQ-TREE 2 long forms and the legacy single-dash forms are listed.
_SUPPORT_FLAGS: dict[str, tuple[tuple[str, ...], bool]] = {
    "alrt":      (("-alrt", "--alrt"), True),
    "lbp":       (("-lbp", "--lbp"), True),
    "abayes":    (("-abayes", "--abayes"), False),
    "ufboot":    (("-B", "--ufboot", "-bb"), True),
    "bootstrap": (("-b", "--boot"), True),
}
# The order IQ-TREE joins several measures in one node label
# ("SH-aLRT / local bootstrap / aBayes / ultrafast bootstrap"); a standard
# bootstrap takes the last slot as UFBoot does (the two are exclusive).
_SUPPORT_ORDER = ("alrt", "lbp", "abayes", "ufboot", "bootstrap")
# --bnni only refines UFBoot, so it goes wherever UFBoot goes.
_UFBOOT_COMPANIONS = ("--bnni",)
# Measure kept when a label carries several, best first; and its name.
_KEPT_PREFERENCE = ("ufboot", "bootstrap", "alrt", "lbp", "abayes")
_SUPPORT_LABELS = {
    "ufboot": "UFBoot", "bootstrap": "Bootstrap", "alrt": "SH_aLRT",
    "lbp": "LBP", "abayes": "aBayes",
}


def support_measures(options_parts: list[str]) -> list[str]:
    """Return the support measures *options_parts* request, in label order."""
    present = {
        name for name, (flags, _) in _SUPPORT_FLAGS.items()
        if any(f in options_parts for f in flags)
    }
    return [m for m in _SUPPORT_ORDER if m in present]


def with_default_support(options_parts: list[str]) -> list[str]:
    """Return *options_parts* with ``-alrt 1000`` added when no SH-aLRT or
    bootstrap measure is requested, so every IQ-TREE tree carries support."""
    if {"alrt", "ufboot", "bootstrap"} & set(support_measures(options_parts)):
        return list(options_parts)
    return ["-alrt", "1000"] + list(options_parts)


def kept_support(options_parts: list[str]) -> tuple[str, int, int]:
    """Return ``(measure, index, count)`` for the value the outputs keep.

    IQ-TREE writes one value per requested measure, slash-joined; the
    wrapper keeps the bootstrap-type value when there is one, else SH-aLRT.
    *options_parts* must already include the default support
    (:func:`with_default_support`).
    """
    measures = support_measures(options_parts)
    if not measures:
        return "alrt", 0, 1
    kept = next(m for m in _KEPT_PREFERENCE if m in measures)
    return kept, measures.index(kept), len(measures)


def support_type_for(tool: str, options: str) -> str:
    """Return a short label naming the branch-support measure carried by the
    tree *tool* writes under *options* — 'SH_like' (FastTree), else one of
    'UFBoot', 'Bootstrap', 'SH_aLRT' (and 'LBP' / 'aBayes' only if a run
    somehow requests nothing else).

    The IQ-TREE wrapper adds ``-alrt 1000`` when no SH-aLRT or bootstrap is
    requested, and keeps one value per node when several are (see
    :func:`kept_support`); this returns the name of that value.
    """
    tool_norm = tool.lower().replace("-", "").replace("_", "")
    if tool_norm == "fasttree":
        return "SH_like"
    parts = with_default_support(options.split() if options else [])
    return _SUPPORT_LABELS[kept_support(parts)[0]]


def strip_ufboot_flags(options_parts: list[str]) -> list[str]:
    """Return *options_parts* without UFBoot (and ``--bnni``), keeping the
    other support measures; ``-alrt 1000`` is added if nothing is left.

    Used for the fallback when IQ-TREE aborts in its post-bootstrap
    consensus step.
    """
    ufboot_flags, _ = _SUPPORT_FLAGS["ufboot"]
    out: list[str] = []
    skip_next = False
    for part in options_parts:
        if skip_next:
            skip_next = False
            continue
        if part in ufboot_flags:
            skip_next = True
            continue
        if part in _UFBOOT_COMPANIONS:
            continue
        out.append(part)
    return with_default_support(out)


# ---------------------------------------------------------------------------
# Run record beside each tree
# ---------------------------------------------------------------------------

# Written by run_tree for every tool and removed before every run, so it
# always describes the tree next to it — cached trees included.
TREE_INFO_SUFFIX = ".tree_info.json"


def _tree_info_path(output_dir: Path, prefix: str) -> Path:
    return Path(output_dir) / f"{prefix}{TREE_INFO_SUFFIX}"


def read_tree_info(output_dir: Path, prefix: str) -> dict | None:
    """Return the run record written beside the tree, or None.

    Keys: ``tool``, ``options`` (the options the support values come from,
    in config terms), ``support_type`` (the measure the tree carries) and
    ``note`` (non-empty when IQ-TREE's UFBoot fallback was used).
    """
    import json
    try:
        info = json.loads(_tree_info_path(output_dir, prefix).read_text())
    except (OSError, ValueError):
        return None
    return info if isinstance(info, dict) and info.get("support_type") else None


def _write_tree_info(output_dir: Path, prefix: str, info: dict) -> None:
    import json
    _tree_info_path(output_dir, prefix).write_text(json.dumps(info, indent=2) + "\n")


# Signatures of an IQ-TREE abort (as opposed to a rejected command line).
_IQTREE_ABORT_MARKERS = ("CRASHES WITH SIGNAL", "Assertion")


def _looks_like_iqtree_abort(stderr: str) -> bool:
    return any(m in (stderr or "") for m in _IQTREE_ABORT_MARKERS)


# A node label of two or more slash-joined numbers, right after ")".  Leaf
# names never follow ")", so they are never touched.
_NUM = r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_MULTI_SUPPORT_RE = re.compile(rf"\)({_NUM}(?:/{_NUM})+)(?=[:,);\[])")


def keep_support_value(newick: str, index: int, count: int) -> str:
    """Reduce IQ-TREE's multi-value node labels (``97.3/0.98/100``) to one value.

    Downstream code (BioPython, the writers, summary statistics) carries one
    numeric confidence per node.  *index* is the position of the kept value
    among *count* (see :func:`kept_support`).  A label whose length differs
    from *count* keeps its last value when the last measure is the one kept,
    else its first.  Single-value labels are returned unchanged; the full
    labels stay in the IQ-TREE ``.treefile``.
    """
    def _pick(m: re.Match) -> str:
        values = m.group(1).split("/")
        if len(values) == count:
            return ")" + values[index]
        return ")" + (values[-1] if index == count - 1 else values[0])

    return _MULTI_SUPPORT_RE.sub(_pick, newick)


def parse_iqtree_partition_models(log_path: Path) -> dict[str, str]:
    """Return ``{charset_name: model}`` for a partitioned IQ-TREE run.

    Reads the sibling ``<prefix>.best_scheme.nex`` written by IQ-TREE 2
    when run with ``-p partitions.nex``.  Format:

    ::

        charpartition mymodels = LG+I+G4: DNApolymerase, WAG+G4: MCP;

    Returns an empty dict when the file is missing, unparseable, or the
    run wasn't partitioned (no charpartition line).
    """
    nex_path = log_path.with_suffix(".best_scheme.nex")
    if not nex_path.exists():
        return {}
    try:
        text = nex_path.read_text()
    except Exception:
        return {}
    m = re.search(r"charpartition\s+\w+\s*=\s*(.+?);", text, re.DOTALL | re.IGNORECASE)
    if not m:
        return {}
    body = m.group(1)
    result: dict[str, str] = {}
    for entry in body.split(","):
        if ":" not in entry:
            continue
        model_part, charset_part = entry.split(":", 1)
        model = model_part.strip()
        charset = charset_part.strip()
        if model and charset:
            result[charset] = model
    return result


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
    if tool_norm not in ("fasttree", "iqtree", "iqtree2"):
        raise ValueError(f"Unsupported tree tool: {tool}. Supported: fasttree, iqtree")

    # A record left by an earlier run (possibly another tool) must never be
    # read back for this tree, even if this run fails part-way.
    _tree_info_path(output_dir, prefix).unlink(missing_ok=True)

    if tool_norm == "fasttree":
        out_nwk = _run_fasttree(alignment_fasta, output_dir, prefix, seq_type,
                                model_nuc, model_aa, options)
        info = {"tool": "fasttree", "options": options or "",
                "support_type": "SH_like", "note": ""}
    else:
        out_nwk, info = _run_iqtree(alignment_fasta, output_dir, prefix, seq_type,
                                    model_nuc, model_aa, options, threads)
    _write_tree_info(output_dir, prefix, info)
    return out_nwk


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

    log.info("Running FastTree (%s)", " ".join(cmd[1:-1]))
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

    log.debug("FastTree complete: %s", out_nwk)
    log.info("FastTree complete")
    return out_nwk


# ---------------------------------------------------------------------------
# IQ-TREE
# ---------------------------------------------------------------------------

def _iqtree_cmd(
    alignment_fasta: Path,
    out_prefix: Path,
    model: str | None,
    seqtype_flag: str,
    threads: int,
    options_parts: list[str],
) -> list[str]:
    """Build an IQ-TREE command line; *model* None leaves ``-m`` out (the
    partition file then carries the models)."""
    cmd = ["iqtree2", "-s", str(alignment_fasta), "--prefix", str(out_prefix)]
    if model:
        cmd += ["-m", model]
    cmd += ["--seqtype", seqtype_flag, "-T", str(threads), "--redo"]
    return cmd + list(options_parts)


def _run_cmd(cmd: list[str]) -> subprocess.CompletedProcess:
    log.debug("Running: %s", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, capture_output=True, text=True)


def _option_value(parts: list[str], flag: str) -> str | None:
    """Return the value following *flag* in *parts*, or None."""
    if flag in parts:
        i = parts.index(flag)
        if i + 1 < len(parts):
            return parts[i + 1]
    return None


def _replace_option(parts: list[str], flag: str, value: str) -> list[str]:
    out = list(parts)
    out[out.index(flag) + 1] = value
    return out


def _logged_best_model(log_path: Path) -> str | None:
    """ModelFinder's pick from an IQ-TREE ``.log`` ("Best-fit model: X chosen
    according to BIC") — present even when the run aborted later."""
    try:
        text = log_path.read_text()
    except OSError:
        return None
    m = re.search(r"Best-fit model:\s*(\S+)\s+chosen according to", text)
    return m.group(1) if m else None


def _ufboot_rescore_parts(
    out_prefix: Path,
    model: str,
    fallback_parts: list[str],
) -> tuple[str | None, list[str]] | None:
    """Plan an SH-aLRT rescore of the ML tree an aborted UFBoot run left.

    IQ-TREE writes the final ML tree to ``.treefile`` (and ModelFinder's
    choices to ``.log`` / ``.best_scheme.nex``) before the consensus step
    that aborts, so scoring that fixed topology is all that is missing.
    Returns ``(model_or_None, options_parts)`` for the rescore, or None when
    the pieces are not there (the caller then re-runs the full search).
    """
    treefile = Path(str(out_prefix) + ".treefile")
    try:
        validate_newick(treefile)
    except (FileNotFoundError, ValueError):
        return None
    ml_tree = Path(str(out_prefix) + ".ml_tree.nwk")
    ml_tree.write_text(treefile.read_text())

    parts = [p for p in fallback_parts if p != "--fast"]
    rescore_model: str | None = model
    if _option_value(parts, "-p") is not None:
        scheme = Path(str(out_prefix) + ".best_scheme.nex")
        if scheme.exists():
            # Read from a copy: the rescore writes under the same prefix.
            scheme_copy = Path(str(out_prefix) + ".fallback_scheme.nex")
            scheme_copy.write_text(scheme.read_text())
            parts = _replace_option(parts, "-p", str(scheme_copy))
            rescore_model = None
        elif is_model_finder_spec(model):
            return None
    elif is_model_finder_spec(model):
        rescore_model = _logged_best_model(Path(str(out_prefix) + ".log"))
        if not rescore_model:
            return None
    return rescore_model, ["-t", str(ml_tree), "--tree-fix"] + parts


def _run_iqtree(
    alignment_fasta: Path,
    output_dir: Path,
    prefix: str,
    seq_type: str,
    model_nuc: str,
    model_aa: str,
    options: str,
    threads: int,
) -> tuple[Path, dict]:
    """Run IQ-TREE; return the Newick path and the run record for it."""
    out_prefix = output_dir / prefix
    model = model_nuc if seq_type == "nucleotide" else model_aa
    # Explicitly set the sequence type rather than relying on IQ-TREE's
    # auto-detection — protein alignments with high A/C/G/T content can be
    # mis-detected as nucleotide.
    seqtype_flag = "DNA" if seq_type == "nucleotide" else "AA"

    # Every tree carries support: -alrt 1000 is added when neither SH-aLRT
    # nor a bootstrap is requested.  --fast works with -alrt but not with
    # -B / UFBoot; that combination gets IQ-TREE's own error.
    configured_parts = options.split() if options else []
    run_parts = with_default_support(configured_parts)
    support_parts = run_parts      # the options the kept support comes from
    described = options or ""      # the same, in config terms, for provenance
    note = ""

    log.info("Running IQ-TREE (model=%s, %s)", model, " ".join(run_parts))
    result = _run_cmd(_iqtree_cmd(alignment_fasta, out_prefix, model, seqtype_flag,
                                  threads, run_parts))

    if (
        result.returncode != 0
        and "ufboot" in support_measures(run_parts)
        and _looks_like_iqtree_abort(result.stderr)
    ):
        # IQ-TREE 2.3.2 can abort in its post-bootstrap consensus step
        # ("Assertion `!findSplit(sp)' failed" in hashsplitset.cpp) on
        # particular taxon sets.  The tree file it leaves holds the final ML
        # tree but no support values, so rather than fail the family, score
        # that tree with the remaining measures (SH-aLRT); only if that is
        # impossible, repeat the search without UFBoot.
        log.debug("IQ-TREE stderr:\n%s", (result.stderr or "")[-4000:])
        support_parts = strip_ufboot_flags(run_parts)
        described = " ".join(strip_ufboot_flags(configured_parts))
        plan = _ufboot_rescore_parts(out_prefix, model, support_parts)
        if plan is not None:
            rescore_model, rescore_parts = plan
            log.warning(
                "IQ-TREE aborted after UFBoot (exit %d), most likely its "
                "bootstrap-consensus bug. Scoring the finished ML tree with %s "
                "instead (tree_info records the downgrade).",
                result.returncode, " ".join(support_parts),
            )
            note = ("UFBoot aborted in IQ-TREE's consensus step; ML tree from that "
                    "run scored with SH-aLRT on a fixed topology")
            result = _run_cmd(_iqtree_cmd(alignment_fasta, out_prefix, rescore_model,
                                          seqtype_flag, threads, rescore_parts))
            scheme = Path(str(out_prefix) + ".best_scheme.nex")
            scheme_copy = Path(str(out_prefix) + ".fallback_scheme.nex")
            if not scheme.exists() and scheme_copy.exists():
                scheme.write_text(scheme_copy.read_text())
        if plan is None or result.returncode != 0:
            log.warning(
                "IQ-TREE aborted after UFBoot (exit %d) and the ML tree could not be "
                "rescored; repeating the search with %s.",
                result.returncode, " ".join(support_parts),
            )
            note = ("UFBoot aborted in IQ-TREE's consensus step; search repeated "
                    "with SH-aLRT only")
            result = _run_cmd(_iqtree_cmd(alignment_fasta, out_prefix, model,
                                          seqtype_flag, threads, support_parts))

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

    # Copy to a consistently named .nwk file, keeping one support value per
    # node when IQ-TREE wrote several ("alrt/ufboot", "alrt/abayes/ufboot").
    out_nwk = output_dir / f"{prefix}.nwk"
    newick = treefile.read_text()
    kept, index, count = kept_support(support_parts)
    if count > 1:
        newick = keep_support_value(newick, index, count)
        log.debug("Kept the %s value of %d support values per node", kept, count)
    out_nwk.write_text(newick)

    log.debug("IQ-TREE complete: %s", out_nwk)
    log.info("IQ-TREE complete")
    info = {
        "tool": "iqtree",
        "options": described,
        "support_type": _SUPPORT_LABELS[kept],
        "note": note,
    }
    return out_nwk, info
