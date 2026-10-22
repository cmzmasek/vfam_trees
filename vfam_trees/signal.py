"""Phylogenetic-signal diagnostics for a finished alignment + tree.

Whole-genome nucleotide alignment of a divergent family can run out of signal:
trimAl keeps few columns, the columns it keeps are barely more similar than
unrelated sequence, and ML distances blow up to IQ-TREE's branch-length cap.
The tree then looks fine but its deep structure is noise.  These measures
make that visible per tree (in ``summary.tsv`` and the log) so such families
can be moved to a protein marker:

- ``core_frac`` — residues a sequence keeps in the trimmed alignment over its
  ungapped input length, median over sequences: how much of a typical
  sequence the tree is built from.
- ``identity_median`` / ``identity_p10`` — pairwise identity over the trimmed
  alignment (shared non-gap columns), on a deterministic sample of pairs.
- ``patristic_median`` / ``patristic_p90`` — pairwise tree distance in
  substitutions per site, on a deterministic sample of tip pairs; the
  classic saturation measure.
- ``n_branches_at_cap`` — branches at IQ-TREE's maximum length (10 by
  default), where the likelihood could not bound the length.

:func:`signal_verdict` turns them into ``ok`` / ``weak`` / ``saturated``.
Pure functions apart from reading the alignment and Newick files; stdlib +
NumPy + BioPython only.
"""
from __future__ import annotations

import random
import statistics
from pathlib import Path

import numpy as np

# IQ-TREE's default maximum branch length (-blmax 10); a branch at it means
# the optimiser hit the bound.  FastTree has no bound.
BRANCH_LENGTH_CAP = 10.0
_AT_CAP = 0.999 * BRANCH_LENGTH_CAP

# Pairs sampled for identity and patristic distances.  Fixed seed, so the
# numbers are reproducible run to run.
N_PAIRS = 2000
_SEED = 20261020

# Letters counted as residues for identity; anything else (gaps, N in a
# nucleotide alignment, IUPAC ambiguity codes, X, B/Z/J, "*", "?") is
# unknown.  In a protein alignment N is asparagine and counts.
_RESIDUES = {
    "nucleotide": np.frombuffer(b"ACGTU", dtype=np.uint8),
    "protein": np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", dtype=np.uint8),
}


def _read_alignment(path: Path) -> np.ndarray:
    """Return the alignment as an (n_seqs, n_cols) uint8 array, upper case."""
    from Bio import SeqIO
    rows = [str(r.seq).upper().encode("ascii", "replace") for r in SeqIO.parse(str(path), "fasta")]
    if not rows:
        return np.zeros((0, 0), dtype=np.uint8)
    width = max(len(r) for r in rows)
    arr = np.full((len(rows), width), ord("-"), dtype=np.uint8)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = np.frombuffer(r, dtype=np.uint8)
    return arr


def _sample_pairs(n: int, n_pairs: int = N_PAIRS) -> list[tuple[int, int]]:
    """All pairs when there are few, else a fixed-seed sample of distinct pairs."""
    total = n * (n - 1) // 2
    if total <= n_pairs:
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    rng = random.Random(_SEED)
    seen: set[tuple[int, int]] = set()
    while len(seen) < n_pairs:
        i, j = rng.randrange(n), rng.randrange(n)
        if i != j:
            seen.add((min(i, j), max(i, j)))
    return sorted(seen)


def alignment_identity(aln_path: Path, seq_type: str = "nucleotide") -> dict:
    """Median and 10th-percentile pairwise identity of an alignment.

    Identity is computed over columns where both sequences have a known
    residue — A/C/G/T/U for nucleotides, the 20 standard amino acids for
    proteins (`_RESIDUES`); pairs sharing fewer than 10 such columns are
    skipped.  Returns ``{"median": float|None, "p10": ...}``.
    """
    arr = _read_alignment(aln_path)
    if arr.shape[0] < 2 or arr.shape[1] == 0:
        return {"median": None, "p10": None}
    residues = _RESIDUES["protein" if seq_type == "protein" else "nucleotide"]
    is_res = np.isin(arr, residues)
    ids = []
    for i, j in _sample_pairs(arr.shape[0]):
        both = is_res[i] & is_res[j]
        n = int(both.sum())
        if n >= 10:
            ids.append(float((arr[i][both] == arr[j][both]).sum()) / n)
    if not ids:
        return {"median": None, "p10": None}
    ids.sort()
    return {"median": statistics.median(ids), "p10": ids[int(0.10 * (len(ids) - 1))]}


def _residue_counts(fasta_path: Path) -> dict[str, int]:
    """Residues (non-gap characters) per sequence id in *fasta_path*."""
    from Bio import SeqIO
    return {
        r.id: len(str(r.seq).replace("-", "").replace(".", ""))
        for r in SeqIO.parse(str(fasta_path), "fasta")
    }


def core_fraction(aln_trimmed: Path, input_fasta: Path) -> float | None:
    """Median over sequences of (residues kept in the trimmed alignment /
    ungapped input length).

    Sequences are matched by id; when the ids do not match (e.g. display
    names restored on one file only), the median kept residues are divided
    by the median input length instead.
    """
    kept = _residue_counts(aln_trimmed)
    full = _residue_counts(input_fasta)
    shared = [k for k in kept if full.get(k)]
    if shared:
        return statistics.median(kept[k] / full[k] for k in shared)
    lengths = [n for n in full.values() if n > 0]
    if not kept or not lengths:
        return None
    return statistics.median(kept.values()) / statistics.median(lengths)


def tree_distances(tree_path: Path) -> dict:
    """Patristic-distance summary and branch-length cap count for a Newick tree.

    Returns ``{"patristic_median", "patristic_p90", "n_branches_at_cap"}``;
    the distances are None when the tree has fewer than two tips.
    """
    from Bio import Phylo
    try:
        tree = Phylo.read(str(tree_path), "newick")
    except Exception:
        return {"patristic_median": None, "patristic_p90": None, "n_branches_at_cap": None}

    # Root-to-node depth and parent links, once; a pair's distance is
    # depth(a) + depth(b) - 2 * depth(lca).
    depth: dict[int, float] = {}
    parent: dict[int, object] = {}
    n_cap = 0
    stack = [(tree.root, 0.0, None)]
    while stack:
        clade, d, par = stack.pop()
        depth[id(clade)] = d
        parent[id(clade)] = par
        for child in clade.clades:
            bl = child.branch_length or 0.0
            if bl >= _AT_CAP:
                n_cap += 1
            stack.append((child, d + bl, clade))

    tips = tree.get_terminals()
    if len(tips) < 2:
        return {"patristic_median": None, "patristic_p90": None, "n_branches_at_cap": n_cap}

    def ancestors(clade) -> list:
        out = []
        while clade is not None:
            out.append(clade)
            clade = parent[id(clade)]
        return out

    lineage = {id(t): ancestors(t) for t in tips}
    dists = []
    for i, j in _sample_pairs(len(tips)):
        a, b = tips[i], tips[j]
        a_anc = {id(c) for c in lineage[id(a)]}
        lca = next(c for c in lineage[id(b)] if id(c) in a_anc)
        dists.append(depth[id(a)] + depth[id(b)] - 2.0 * depth[id(lca)])
    dists.sort()
    return {
        "patristic_median": statistics.median(dists),
        "patristic_p90": dists[int(0.90 * (len(dists) - 1))],
        "n_branches_at_cap": n_cap,
    }


# Verdict thresholds per sequence type.  The median tip-to-tip distance is
# the primary measure: past ~1.5 nucleotide substitutions per site the
# typical pair is saturated (under JC that is ~35% identity, close to what
# unrelated aligned sequence reaches), so deep relationships are not
# recoverable from the alignment.  Amino-acid data saturate far later.
# Calibrated on a 2026-10-20 survey of 127 families (CHANGELOG 1.6.0):
# Xinmoviridae, Hantaviridae (fish viruses at the cap) and most RNA-virus
# whole-genome alignments come out saturated; shallow families ok.
SIGNAL_THRESHOLDS = {
    "nucleotide": {"saturated": 1.5, "weak": 0.75, "identity_weak": 0.60, "core_weak": 0.25},
    "protein":    {"saturated": 4.0, "weak": 3.0,  "identity_weak": 0.30, "core_weak": 0.25},
}


def signal_verdict(seq_type: str, metrics: dict) -> tuple[str, list[str]]:
    """Return ``(verdict, reasons)`` for one tree's signal metrics.

    *metrics* holds ``patristic_median``, ``core_frac``, ``identity_median``
    and ``n_branches_at_cap`` (any may be None).

    - ``saturated``: median tip-to-tip distance above the saturated
      threshold, or a branch at IQ-TREE's cap in a tree whose median is
      past the weak threshold.
    - ``weak``: median past the weak threshold, median identity or core
      fraction below theirs, or any branch at the cap.
    - ``ok`` otherwise; ``""`` when nothing could be measured.
    """
    t = SIGNAL_THRESHOLDS["protein" if seq_type == "protein" else "nucleotide"]
    pat = metrics.get("patristic_median")
    core = metrics.get("core_frac")
    ident = metrics.get("identity_median")
    cap = metrics.get("n_branches_at_cap") or 0
    if pat is None and core is None and ident is None:
        return "", []

    cap_reason = f"{cap} branch(es) at the {BRANCH_LENGTH_CAP:g} subs/site cap"
    if pat is not None and pat > t["saturated"]:
        reasons = [f"median tip-to-tip distance {pat:.2f} subs/site > {t['saturated']}"]
        return "saturated", reasons + ([cap_reason] if cap else [])
    if cap and pat is not None and pat > t["weak"]:
        return "saturated", [cap_reason, f"median tip-to-tip distance {pat:.2f} subs/site"]

    reasons = []
    if pat is not None and pat > t["weak"]:
        reasons.append(f"median tip-to-tip distance {pat:.2f} subs/site > {t['weak']}")
    if ident is not None and ident < t["identity_weak"]:
        reasons.append(f"median identity {ident:.0%} < {t['identity_weak']:.0%}")
    if core is not None and core < t["core_weak"]:
        reasons.append(f"only {core:.0%} of a typical sequence kept after trimming")
    if cap:
        reasons.append(cap_reason)
    return ("weak" if reasons else "ok"), reasons


def measure(aln_trimmed: Path, input_fasta: Path | None, tree_path: Path, seq_type: str) -> dict:
    """All signal metrics for one tree, plus the verdict and its reasons.

    *aln_trimmed* is the alignment the tree was built from, *input_fasta* the
    unaligned sequences that went into it (None: no core fraction, e.g. a
    concatenation built from per-marker alignments), *tree_path* the tree.
    """
    ident = alignment_identity(aln_trimmed, seq_type)
    dist = tree_distances(tree_path)
    metrics = {
        "core_frac": core_fraction(aln_trimmed, input_fasta) if input_fasta else None,
        "identity_median": ident["median"],
        "identity_p10": ident["p10"],
        **dist,
    }
    verdict, reasons = signal_verdict(seq_type, metrics)
    metrics["verdict"] = verdict
    metrics["reasons"] = reasons
    return metrics


def summary_fields(metrics: dict) -> dict:
    """The summary.tsv values for *metrics* (rounded; '' when unmeasured)."""
    def r(key, nd):
        v = metrics.get(key)
        return round(v, nd) if isinstance(v, (int, float)) else ""
    cap = metrics.get("n_branches_at_cap")
    return {
        "signal":            metrics.get("verdict", ""),
        "patristic_median":  r("patristic_median", 3),
        "identity_median":   r("identity_median", 3),
        "core_frac":         r("core_frac", 3),
        "branches_at_cap":   cap if isinstance(cap, int) else "",
    }


def log_signal(log, label: str, seq_type: str, metrics: dict) -> None:
    """Log one tree's signal: INFO always, WARNING with advice when not ok."""
    f = summary_fields(metrics)
    log.info(
        "tree_%s signal: %s (median tip-to-tip %s subs/site, identity %s, "
        "core %s, branches at cap %s)",
        label, f["signal"] or "not measured", f["patristic_median"],
        f["identity_median"], f["core_frac"], f["branches_at_cap"],
    )
    verdict = metrics.get("verdict")
    if verdict not in ("weak", "saturated"):
        return
    advice = (
        "the deep structure of this tree is not reliable; a conserved protein "
        "marker (sequence.type: protein, e.g. the RdRp for RNA viruses) carries it far better"
        if seq_type != "protein" else
        "the deep structure of this tree is not reliable; check the marker choice"
    ) if verdict == "saturated" else "deep relationships may be poorly resolved"
    log.warning(
        "tree_%s phylogenetic signal is %s: %s — %s.",
        label, verdict.upper(), "; ".join(metrics.get("reasons") or []), advice,
    )
