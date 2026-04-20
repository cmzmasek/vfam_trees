"""Taxonomy annotation using lineages stored in GenBank records — no external DB needed."""
from __future__ import annotations

import copy
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioTree, Clade

from .logger import get_logger

log = get_logger(__name__)


def annotate_tree(
    newick_path: Path,
    id_map: dict[str, str],
    metadata: list[dict],
    output_nwk: Path,
) -> BioTree | None:
    """Load tree, root it optimally, annotate internal nodes with LCA taxonomy.

    Rooting strategy (in priority order):
      1. Taxonomy-guided rooting: find the root position that maximises
         mean LCA specificity across all internal nodes
      2. Midpoint rooting as final fallback

    LCA is derived from taxonomic lineages stored in GenBank record annotations.

    Args:
        newick_path: input Newick (short IDs)
        id_map: short_id → display_name
        metadata: list of metadata dicts (each must have 'short_id' and 'lineage')
        output_nwk: path to write annotated Newick

    Returns:
        Annotated BioPython Tree, or None on failure.
    """
    trees = list(Phylo.parse(str(newick_path), "newick"))
    if not trees:
        log.error("Could not parse tree from %s", newick_path)
        return None
    tree = trees[0]

    # Build short_id → lineage mapping
    short_id_to_lineage: dict[str, list[str]] = {
        meta["short_id"]: meta["lineage"]
        for meta in metadata
        if meta.get("short_id") and meta.get("lineage")
    }

    # Build authoritative name → rank map from NCBI ranked lineages
    name_to_rank: dict[str, str] = {}
    for meta in metadata:
        for entry in meta.get("lineage_ranked", []) or []:
            name = entry.get("name")
            rank = entry.get("rank", "")
            if name and rank and rank != "no rank" and name not in name_to_rank:
                name_to_rank[name] = rank

    # Root the tree
    _root_tree(tree, short_id_to_lineage)

    # Annotate internal nodes with LCA from lineages
    if short_id_to_lineage:
        _annotate_internal_nodes(tree.root, short_id_to_lineage, name_to_rank)
        _keep_deepest_labels(tree)
    else:
        log.warning("No lineage data available — skipping taxonomy annotation")

    _normalize_bootstrap(tree)

    output_nwk.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output_nwk), "newick")
    log.info("Annotated tree written to %s", output_nwk)
    return tree


# ---------------------------------------------------------------------------
# Rooting
# ---------------------------------------------------------------------------

def _root_tree(
    tree: BioTree,
    short_id_to_lineage: dict[str, list[str]],
) -> None:
    """Root the tree using the best available strategy."""
    # Priority 1: taxonomy-guided rooting
    if short_id_to_lineage:
        success = _taxonomy_guided_root(tree, short_id_to_lineage)
        if success:
            return

    # Priority 2: MAD rooting
    success = _mad_root(tree)
    if success:
        return

    # Priority 3: midpoint rooting (last resort)
    try:
        tree.root_at_midpoint()
        log.info("Rooted tree at midpoint (fallback)")
    except Exception as e:
        log.warning("Midpoint rooting failed (%s) — tree left unrooted", e)


def _taxonomy_guided_root(
    tree: BioTree,
    short_id_to_lineage: dict[str, list[str]],
) -> bool:
    """Try every branch as root, pick the one with highest mean LCA specificity.

    Specificity of an internal node = depth of its LCA taxon in the lineage
    (number of elements in the common lineage prefix). Higher = more specific.
    Mean specificity is averaged across all internal nodes with lineage data.

    Returns True if a root was applied, False if scoring failed.
    """
    log.info("Computing optimal root via taxonomy-guided search ...")

    # Collect all branches (non-root internal clades + all terminals)
    branches = list(tree.find_clades())
    if len(branches) < 3:
        return False

    best_score = -1.0
    best_clade = None
    best_branch_len = float("inf")

    for clade in branches:
        if clade is tree.root:
            continue
        score = _score_rooting(tree, clade, short_id_to_lineage)
        branch_len = clade.branch_length or 0.0
        if score > best_score or (score == best_score and branch_len < best_branch_len):
            best_score = score
            best_clade = clade
            best_branch_len = branch_len

    if best_clade is None:
        log.warning("Taxonomy-guided rooting: no valid root found")
        return False

    try:
        tree.root_with_outgroup(best_clade)
        log.info("Taxonomy-guided root applied (mean specificity score=%.3f)", best_score)
        return True
    except Exception as e:
        log.warning("Could not apply taxonomy-guided root (%s)", e)
        return False


def _score_rooting(
    tree: BioTree,
    clade: Clade,
    short_id_to_lineage: dict[str, list[str]],
) -> float:
    """Score a candidate root placement without permanently modifying the tree.

    Temporarily re-roots at `clade`, computes mean LCA specificity across all
    internal nodes, then restores the original root.

    Specificity of a node = length of the LCA lineage prefix of its leaves.
    """
    # Deep-copy the tree to avoid mutating the original
    try:
        target_leaves = frozenset(c.name for c in clade.get_terminals())
        tmp_tree = copy.deepcopy(tree)
        tmp_clade = _find_clade_by_leaves(tmp_tree, target_leaves)
        if tmp_clade is None:
            return -1.0

        tmp_tree.root_with_outgroup(tmp_clade)
        return _mean_lca_specificity(tmp_tree.root, short_id_to_lineage)
    except Exception:
        return -1.0


def _find_clade_by_leaves(tree: BioTree, target_leaves: frozenset) -> Clade | None:
    """Find a clade whose terminal set exactly matches target_leaves."""
    for clade in tree.find_clades():
        if frozenset(c.name for c in clade.get_terminals()) == target_leaves:
            return clade
    return None


def _mad_root(tree: BioTree) -> bool:
    """Root tree using Minimal Ancestor Deviation (Tria et al. 2017).

    For every branch, analytically finds the point along it that minimises
    the mean squared relative deviation ρ²  across all leaf pairs:

        ρ(i,j) = (d(i,r) - d(j,r)) / d(i,j)

    where d(i,r) is the path length from leaf i to candidate root r.
    The branch+position with the lowest mean ρ² is applied as the root.

    Returns True if a root was applied, False otherwise.
    """
    import collections

    terminals = tree.get_terminals()
    n = len(terminals)
    if n < 3:
        return False

    # ---- Build adjacency list (node id → [(neighbour, edge_len), ...]) ----
    adj: dict[int, list] = {}
    for clade in tree.find_clades():
        for child in clade.clades:
            edge = child.branch_length or 0.0
            adj.setdefault(id(clade), []).append((child, edge))
            adj.setdefault(id(child), []).append((clade, edge))

    # ---- Pairwise leaf distances via BFS from each terminal ----
    pdist: dict[tuple[str, str], float] = {}
    for start in terminals:
        queue: collections.deque = collections.deque([(start, 0.0)])
        visited: set[int] = {id(start)}
        while queue:
            node, d = queue.popleft()
            if node.is_terminal() and node is not start:
                pdist[(start.name, node.name)] = d
                pdist[(node.name, start.name)] = d
            for nbr, edge in adj.get(id(node), []):
                if id(nbr) not in visited:
                    visited.add(id(nbr))
                    queue.append((nbr, d + edge))

    # ---- Subtree distances: for each clade node, dist to each descendant leaf ----
    subdist: dict[int, dict[str, float]] = {}

    def _build_subdist(clade: Clade) -> dict[str, float]:
        if clade.is_terminal():
            subdist[id(clade)] = {clade.name: 0.0}
            return subdist[id(clade)]
        result: dict[str, float] = {}
        for child in clade.clades:
            edge = child.branch_length or 0.0
            for name, d in _build_subdist(child).items():
                result[name] = d + edge
        subdist[id(clade)] = result
        return result

    _build_subdist(tree.root)

    all_names = {t.name for t in terminals}
    n_pairs = n * (n - 1) / 2

    best_score = float("inf")
    best_clade: Clade | None = None
    best_x = 0.0

    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        L = clade.branch_length
        if not L or L <= 0:
            continue

        B = subdist[id(clade)]          # {name: dist_from_clade_root}
        A_names = all_names - B.keys()
        if not A_names:
            continue

        # a[i] = dist from leaf i to parent node u = d(i, j0) - L - b[j0]
        j0 = next(iter(B))
        b_j0 = B[j0]
        a = {ni: pdist[(ni, j0)] - L - b_j0 for ni in A_names}

        A_list = list(A_names)
        B_list = list(B)

        # Optimal x along branch: minimise Σ_cross ρ²  analytically
        sum_num = 0.0
        sum_den = 0.0
        for ni in A_list:
            for nj in B_list:
                d_ij = a[ni] + L + B[nj]
                if d_ij <= 0:
                    continue
                inv_d2 = 1.0 / (d_ij * d_ij)
                sum_num += (a[ni] - B[nj] - L) * inv_d2
                sum_den += inv_d2

        if sum_den == 0:
            continue

        x = max(0.0, min(L, -sum_num / (2.0 * sum_den)))

        # Compute MAD score at x
        mad = 0.0

        for ni in A_list:
            for nj in B_list:
                d_ij = a[ni] + L + B[nj]
                if d_ij <= 0:
                    continue
                rho = (a[ni] + x - B[nj] - (L - x)) / d_ij
                mad += rho * rho

        for p in range(len(A_list)):
            for q in range(p + 1, len(A_list)):
                d_ij = pdist.get((A_list[p], A_list[q]), 0.0)
                if d_ij <= 0:
                    continue
                mad += ((a[A_list[p]] - a[A_list[q]]) / d_ij) ** 2

        for p in range(len(B_list)):
            for q in range(p + 1, len(B_list)):
                d_ij = pdist.get((B_list[p], B_list[q]), 0.0)
                if d_ij <= 0:
                    continue
                mad += ((B[B_list[p]] - B[B_list[q]]) / d_ij) ** 2

        score = mad / n_pairs
        if score < best_score:
            best_score = score
            best_clade = clade
            best_x = x

    if best_clade is None:
        log.warning("MAD rooting: no valid branch found")
        return False

    B_leaves = frozenset(c.name for c in best_clade.get_terminals())
    orig_len = best_clade.branch_length
    try:
        tree.root_with_outgroup(best_clade)
        # Adjust branch lengths to reflect the optimal split point
        for child in tree.root.clades:
            child_leaves = frozenset(c.name for c in child.get_terminals())
            if child_leaves == B_leaves:
                child.branch_length = orig_len - best_x
            else:
                child.branch_length = best_x
        log.info("MAD rooting applied (score=%.6f)", best_score)
        return True
    except Exception as e:
        log.warning("Could not apply MAD root (%s)", e)
        return False


def _mean_lca_specificity(
    root: Clade,
    short_id_to_lineage: dict[str, list[str]],
) -> float:
    """Compute clade-size-weighted mean LCA lineage depth across internal nodes.

    Each internal node contributes `specificity × terminal_count`, so large
    root-ward clades dominate the score. This rewards rootings that get the
    deep taxonomic structure right rather than ones that carve out many tiny,
    highly-specific sub-clades.
    """
    weighted_sum = 0.0
    total_weight = 0
    for clade in root.find_clades():
        if clade.is_terminal():
            continue
        terminals = clade.get_terminals()
        leaf_lineages = [
            short_id_to_lineage[leaf.name]
            for leaf in terminals
            if leaf.name in short_id_to_lineage
        ]
        if not leaf_lineages or len(leaf_lineages) * 2 < len(terminals):
            continue
        lca_lineage = _find_lca_lineage(leaf_lineages)
        weight = len(terminals)
        weighted_sum += len(lca_lineage) * weight
        total_weight += weight

    return weighted_sum / total_weight if total_weight else 0.0


# ---------------------------------------------------------------------------
# LCA computation
# ---------------------------------------------------------------------------

# ICTV suffix → PhyloXML rank string.  Ordered longest-suffix first so that
# e.g. "viricetes" is matched before the shorter "virae" suffix.
_ICTV_RANK_SUFFIXES: list[tuple[str, str]] = [
    ("viricotina",  "subphylum"),
    ("viricetes",   "class"),
    ("viricota",    "phylum"),
    ("virinae",     "subfamily"),
    ("viridae",     "family"),
    ("virales",     "order"),
    ("virae",       "kingdom"),
    ("viria",       "realm"),
    ("vira",        "subrealm"),
]


def _infer_rank(name: str) -> str:
    """Infer PhyloXML rank string from an ICTV taxon name suffix.

    Returns an empty string when the rank cannot be determined. Single-word
    names ending in "virus" are treated as genus; multi-word names ending in
    "virus" are vernacular species names and left unranked.
    """
    lower = name.lower()
    for suffix, rank in _ICTV_RANK_SUFFIXES:
        if lower.endswith(suffix):
            return rank
    if lower.endswith("virus") and " " not in name.strip():
        return "genus"
    return ""


def _annotate_internal_nodes(
    clade: Clade,
    short_id_to_lineage: dict[str, list[str]],
    name_to_rank: dict[str, str] | None = None,
) -> None:
    """Recursively annotate internal nodes with LCA taxon from leaf lineages.

    Sets clade.name to the LCA taxon name and stores the rank from NCBI
    taxonomy (or the suffix-based fallback) in clade._taxonomy_rank.
    """
    if clade.is_terminal():
        return

    for child in clade.clades:
        _annotate_internal_nodes(child, short_id_to_lineage, name_to_rank)

    terminals = clade.get_terminals()
    leaf_lineages = [
        short_id_to_lineage[leaf.name]
        for leaf in terminals
        if leaf.name in short_id_to_lineage
    ]

    # Require ≥50% lineage coverage to avoid over-specific labels driven by
    # a small annotated minority.
    if not leaf_lineages or len(leaf_lineages) * 2 < len(terminals):
        return

    lca_lineage = _find_lca_lineage(leaf_lineages)
    if lca_lineage:
        taxon = lca_lineage[-1]
        clade.name = taxon
        rank = (name_to_rank or {}).get(taxon, "")
        clade._taxonomy_rank = rank or _infer_rank(taxon)


def _find_lca_lineage(lineages: list[list[str]]) -> list[str]:
    """Return the common prefix lineage across all input lineages."""
    if not lineages:
        return []
    common = lineages[0][:]
    for lin in lineages[1:]:
        trimmed = []
        for a, b in zip(common, lin):
            if a == b:
                trimmed.append(a)
            else:
                break
        common = trimmed
    return common


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def _keep_deepest_labels(tree: BioTree) -> None:
    """For each taxon label, keep it only at the most inclusive (crown) node.

    Sorting descending by terminal count ensures the largest clade (closest to
    root) is processed first for each label.  Any smaller sub-clade that shares
    the same LCA name is then cleared — it is redundant because the taxon is
    already annotated at its true crown node.

    This guarantees, for example, that the family name appears at the root of a
    single-family tree rather than being displaced by a smaller multi-genus clade
    that happens to share the same LCA.
    """
    seen_taxa: set[str] = set()
    # Primary key: descending terminal count (crown node first).
    # Secondary key: ascending preorder index for deterministic tie-breaking
    # when two clades share a label with the same terminal count.
    internal_nodes = []
    for idx, c in enumerate(tree.find_clades(order="preorder")):
        if not c.is_terminal() and c.name:
            internal_nodes.append((-len(c.get_terminals()), idx, c))
    internal_nodes.sort(key=lambda t: (t[0], t[1]))

    for _, _, clade in internal_nodes:
        taxon = clade.name or ""
        if not taxon:
            continue
        if taxon in seen_taxa:
            clade.name = ""
        else:
            seen_taxa.add(taxon)

    log.debug("Kept crown label for %d distinct taxa", len(seen_taxa))


def _normalize_bootstrap(tree: BioTree) -> None:
    """Scale confidence values to 0–100 if they appear to be in 0–1 range, then round to integer."""
    confidences = [
        c.confidence for c in tree.find_clades()
        if c.confidence is not None
    ]
    if not confidences:
        return
    if max(confidences) <= 1.0:
        for clade in tree.find_clades():
            if clade.confidence is not None:
                clade.confidence = round(clade.confidence * 100)
        log.debug("Bootstrap values scaled from [0,1] to [0,100]")
    else:
        for clade in tree.find_clades():
            if clade.confidence is not None:
                clade.confidence = round(clade.confidence)
