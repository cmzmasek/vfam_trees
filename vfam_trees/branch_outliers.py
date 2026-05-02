"""Branch-length-outlier helpers shared by single-protein and concat pipelines.

Implements the median + factor × MAD criterion described in CONCAT_DESIGN.md
§5.5 / READMEsection on iterative outlier removal.  Pure functions — no I/O
beyond reading the Newick file.
"""
from __future__ import annotations

import statistics
from pathlib import Path

from Bio import Phylo


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
    bls = [b for b in bl_map.values() if b > 0]
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
