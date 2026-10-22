"""Branch-length-outlier helpers shared by single-protein and concat pipelines.

Implements the median + factor × MAD criterion described in CONCAT_DESIGN.md
§5.5 / READMEsection on iterative outlier removal.  No I/O beyond reading the
Newick file; the tree builders are passed in.
"""
from __future__ import annotations

import statistics
from pathlib import Path
from typing import Callable

from Bio import Phylo

# Terminal branches at or below this length count as zero (identical or
# near-identical sequences) and are left out of the median / MAD.  FastTree
# writes such branches as 0.0, but IQ-TREE writes its minimum branch length,
# 1e-6, so a plain "> 0" test let them in.  1e-5 subs/site is under one
# substitution even for a 100 kb alignment, so no real branch is lost.
ZERO_BRANCH_LENGTH = 1e-5


def branch_length_stats(treefile: Path) -> dict:
    """Return median, MAD, and per-leaf branch-length map for a Newick file.

    ``mad`` is the median absolute deviation; together with the median it
    underpins the robust outlier criterion ``bl > median + factor × mad``.
    Returns zeros when the tree is unparseable or has fewer than 3 informative
    branches; ``bl_map`` is still populated with whatever leaves were read.
    """
    try:
        bio_tree = next(iter(Phylo.parse(str(treefile), "newick")))
    except Exception:
        return {"median": 0.0, "mad": 0.0, "bl_map": {}}
    bl_map = {
        c.name: c.branch_length for c in bio_tree.get_terminals()
        if c.name and c.branch_length is not None
    }
    bls = [b for b in bl_map.values() if b > ZERO_BRANCH_LENGTH]
    if len(bls) < 3:
        return {"median": 0.0, "mad": 0.0, "bl_map": bl_map}
    median_bl = statistics.median(bls)
    mad = statistics.median([abs(x - median_bl) for x in bls])
    return {"median": median_bl, "mad": mad, "bl_map": bl_map}


def find_branch_length_outliers(treefile: Path, factor: float) -> set[str]:
    """Return leaf IDs whose branch length exceeds ``median + factor × mad``.

    Returns an empty set when the tree has too few informative branches or
    when the median / MAD are zero (e.g., star tree, identical sequences).
    """
    stats = branch_length_stats(treefile)
    median_bl, mad = stats["median"], stats["mad"]
    if median_bl == 0 or mad == 0:
        return set()
    threshold = median_bl + factor * mad
    return {
        name for name, bl in stats["bl_map"].items()
        if bl is not None and bl > threshold
    }


# The pipelines stop removing outliers below this many leaves.
MIN_TREE_LEAVES = 4
# Log tag for decisions made on the quick tree.
QUICK_STAGE = " (quick tree)"


def select_outliers(
    treefile: Path,
    factor: float,
    protected: set[str],
    n_current: int,
    min_seqs: int,
    *,
    log,
    label: str,
    iteration: int,
    display: Callable[[str], str],
    noun: str = "sequence",
    stage: str = "",
) -> set[str]:
    """Return the leaves to remove from *treefile* as branch-length outliers.

    Shared by the single-marker and the concatenated pipeline, with the same
    log lines.  *protected* leaves (RefSeqs, forced-in records) are reported,
    never removed.  An empty set means stop: no outliers, only protected
    ones, or removing them would leave fewer than *min_seqs* of the
    *n_current* leaves.  *display* turns an id into a name for the log;
    *noun* names a leaf ("sequence" / "genome"); *stage* tags the tree the
    decision was made on (`QUICK_STAGE`).
    """
    tag = f"tree_{label} iteration {iteration}{stage}"
    found = find_branch_length_outliers(treefile, factor)
    if not found:
        log.info("%s: no branch-length outliers — stopping.", tag)
        return set()

    stats = branch_length_stats(treefile)
    median_bl, mad, bl_map = stats["median"], stats["mad"], stats["bl_map"]
    threshold = median_bl + factor * mad
    log.info(
        "%s: branch-length stats — median=%.5f, MAD=%.5f, "
        "threshold (median + %.1f×MAD)=%.5f",
        tag, median_bl, mad, factor, threshold,
    )

    def _describe(oid: str) -> tuple[float, float]:
        bl = bl_map.get(oid, float("nan"))
        return bl, ((bl - median_bl) / mad if mad else float("nan"))

    for oid in sorted(found & protected):
        bl, above = _describe(oid)
        log.warning(
            "%s: protected %s '%s' (RefSeq or forced in) looks like a "
            "branch-length outlier — branch length %.5f (%.1f× MAD above "
            "median, threshold=%.5f) — KEEPING",
            tag, noun, display(oid), bl, above, threshold,
        )
    removable = found - protected
    if not removable:
        log.info("%s: only protected outliers remain — stopping.", tag)
        return set()
    remaining = n_current - len(removable)
    if remaining < min_seqs:
        log.info(
            "%s: removing %d outlier(s) would leave only %d %ss (min_seqs=%d) — stopping.",
            tag, len(removable), remaining, noun, min_seqs,
        )
        return set()
    for oid in sorted(removable):
        bl, above = _describe(oid)
        log.info(
            "%s: removing outlier %s '%s' — branch length %.5f "
            "(%.1f× MAD above median, threshold=%.5f)",
            tag, noun, display(oid), bl, above, threshold,
        )
    return removable


def outlier_round(
    build_tree: Callable[[], Path],
    decide: Callable[[Path, str], set[str]],
    *,
    check: bool,
    n_current: int,
    build_quick_tree: Callable[[], Path] | None = None,
) -> tuple[Path | None, set[str]]:
    """One round of tree building and branch-length outlier removal.

    *build_tree* builds the configured tree; *decide(treefile, stage)* returns
    the ids to remove (`select_outliers`); *check* says whether this round
    looks for outliers at all.  With *build_quick_tree* (and *check*), the
    round first looks on that quick tree: when it finds outliers to remove,
    the configured search is skipped this round and ``(None, ids)`` is
    returned.  Otherwise the configured tree is built and checked as before,
    so the final tree is always the configured one, checked the same way.
    Returns ``(treefile or None, ids to remove)``.
    """
    if check and build_quick_tree is not None:
        remove = decide(build_quick_tree(), QUICK_STAGE)
        if remove and n_current - len(remove) >= MIN_TREE_LEAVES:
            return None, remove
    treefile = build_tree()
    if not check:
        return treefile, set()
    return treefile, decide(treefile, "")
