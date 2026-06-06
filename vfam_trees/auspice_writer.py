"""Generate Auspice v2 JSON output (Nextstrain interactive tree).

This is a *divergence* tree exporter: branch lengths are substitutions/site,
not time, so no temporal layout (``num_date``) is emitted.  Genus colors reuse
the same palette as the PDF/PNG/PhyloXML output by supplying an explicit
``scale`` on the ``genus`` coloring, so the "color by genus" view in Auspice
matches the other outputs exactly.

Single source of truth for Auspice export: ``write_auspice_json``.  Both the
single-protein pipeline (``pipeline.py``) and the concatenated-marker pipeline
(``pipeline_concat.py``) call it with the same data they already hand to
``write_phyloxml``.
"""
from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path

from Bio import Phylo

from . import __version__
from .logger import get_logger

log = get_logger(__name__)

_YEAR_RE = re.compile(r'\b(1[89]\d{2}|20\d{2})\b')

# Auspice v2 reserves ``node_attrs.accession`` as a plain string (not a
# ``{"value": ...}`` coloring object) and constrains it to this pattern.  A
# value that violates either rule fails schema validation for the *whole* tree,
# so we emit accession only when it matches.
_ACCESSION_RE = re.compile(r'^[0-9A-Za-z._-]+$')

# Placeholder species value emitted upstream by extract_metadata; treated as
# "no value" here, mirroring phyloxml_writer.
_PLACEHOLDERS = frozenset({"", "unknown", "n/a", "na", "none"})

# Categorical traits surfaced as colorings / filters when present on any leaf.
# (key, human title)
_CATEGORICAL_TRAITS = [
    ("genus",     "Genus"),
    ("subgenus",  "Subgenus"),
    ("subfamily", "Subfamily"),
    ("species",   "Species"),
    ("host",      "Host"),
    ("location",  "Location"),
    ("strain",    "Strain"),
    ("outbreak",  "Outbreak"),
]


def _extract_year(date_str: str) -> str:
    if not date_str:
        return ""
    m = _YEAR_RE.search(date_str)
    return m.group(1) if m else ""


def _clean(val) -> str:
    """Return a stripped string, or '' for placeholder / empty values."""
    s = str(val).strip() if val is not None else ""
    return "" if s.lower() in _PLACEHOLDERS else s


def _leaf_traits(meta: dict) -> dict[str, str]:
    """Extract categorical trait values for a leaf from its metadata dict."""
    taxa = {"genus": "", "subgenus": "", "subfamily": ""}
    for entry in meta.get("lineage_ranked", []) or []:
        if isinstance(entry, dict):
            rank = (entry.get("rank") or "").lower()
            if rank in taxa:
                taxa[rank] = entry.get("name", "")
    traits = {
        "genus":     _clean(taxa["genus"]),
        "subgenus":  _clean(taxa["subgenus"]),
        "subfamily": _clean(taxa["subfamily"]),
        "species":   _clean(meta.get("species", "")),
        "host":      _clean(meta.get("host", "")),
        "location":  _clean(meta.get("location", "")),
        "strain":    _clean(meta.get("strain", "")),
        "outbreak":  _clean(meta.get("outbreak", "")),
    }
    return {k: v for k, v in traits.items() if v}


def _format_support(conf) -> str:
    """Render a support value as a short string for a branch label."""
    try:
        f = float(conf)
    except (TypeError, ValueError):
        return str(conf)
    # FastTree emits 0–1, IQ-TREE 0–100; normalize sub-1 values to a percentage
    # so branch labels read consistently regardless of upstream tool.
    if f <= 1.0:
        f *= 100.0
    return f"{f:g}"


def write_auspice_json(
    output_json: Path,
    id_map: dict[str, str],
    leaf_metadata: dict[str, dict],
    family: str,
    tree=None,
    newick_path: Path | None = None,
    *,
    title: str | None = None,
    genus_to_color: dict[str, str] | None = None,
    confidence_type: str = "SH_aLRT",
    seq_type: str = "",
    description: str | None = None,
) -> None:
    """Write an Auspice v2 JSON tree.

    Args:
        output_json: destination ``*_auspice.json`` path
        id_map: short_id → display_name (used as the shown tip label)
        leaf_metadata: short_id → metadata dict (same dict given to phyloxml)
        family: biological family name (output / fallback title)
        tree: pre-parsed BioPython tree; parsed from ``newick_path`` if omitted
        newick_path: annotated Newick fallback when ``tree`` is None
        title: phylogeny title shown in Auspice (defaults to family)
        genus_to_color: genus → hex; emitted as the ``genus`` coloring scale so
            Auspice reproduces the PDF/PhyloXML palette
        confidence_type: branch-support label name (e.g. SH_aLRT, SH-like)
        seq_type: "protein" / "nucleotide" — recorded in the description only
        description: optional markdown shown in the Auspice sidebar
    """
    if tree is None:
        if newick_path is None:
            raise ValueError("write_auspice_json requires either tree or newick_path")
        tree = list(Phylo.parse(str(newick_path), "newick"))[0]
    tree.ladderize(reverse=True)

    # ------------------------------------------------------------------
    # Build the nested node structure, accumulating root→node divergence.
    # ------------------------------------------------------------------
    present_traits: set[str] = set()
    has_year = False
    has_display = False
    support_label = re.sub(r"[^0-9A-Za-z_]", "_", confidence_type) or "support"
    counter = {"n": 0}

    def build(clade, parent_div: float) -> dict:
        nonlocal has_year, has_display
        div = parent_div + (clade.branch_length or 0.0)
        node_attrs: dict = {"div": div}
        branch_attrs: dict = {}

        if clade.is_terminal():
            short_id = clade.name or ""
            name = short_id
            meta = leaf_metadata.get(short_id, {})
            for key, value in _leaf_traits(meta).items():
                node_attrs[key] = {"value": value}
                present_traits.add(key)
            year = _extract_year(meta.get("collection_date", ""))
            if year:
                node_attrs["year"] = {"value": int(year)}
                has_year = True
            accession = _clean(meta.get("accession", ""))
            if accession and _ACCESSION_RE.match(accession):
                # Reserved Auspice field: plain string, not a value-object.
                node_attrs["accession"] = accession
            display = id_map.get(short_id, short_id)
            if display:
                node_attrs["display_name"] = {"value": display}
                has_display = True
        else:
            counter["n"] += 1
            name = f"NODE_{counter['n']:07d}"
            labels: dict[str, str] = {}
            lca = (clade.name or "").strip()
            if lca:
                labels["clade"] = lca
            if clade.confidence is not None:
                labels[support_label] = _format_support(clade.confidence)
            if labels:
                branch_attrs["labels"] = labels

        node: dict = {"name": name, "node_attrs": node_attrs}
        if branch_attrs:
            node["branch_attrs"] = branch_attrs
        if not clade.is_terminal():
            node["children"] = [build(c, div) for c in clade.clades]
        return node

    tree_root = build(tree.root, 0.0)

    # ------------------------------------------------------------------
    # meta: colorings (only for traits actually present), filters, defaults.
    # ------------------------------------------------------------------
    colorings: list[dict] = []
    for key, title_ in _CATEGORICAL_TRAITS:
        if key not in present_traits:
            continue
        entry: dict = {"key": key, "title": title_, "type": "categorical"}
        if key == "genus" and genus_to_color:
            entry["scale"] = [[g, c] for g, c in sorted(genus_to_color.items())]
        colorings.append(entry)
    if has_year:
        colorings.append({"key": "year", "title": "Year", "type": "continuous"})

    color_by = "genus" if "genus" in present_traits else (
        colorings[0]["key"] if colorings else "div"
    )
    filter_keys = [k for k in ("genus", "subfamily", "host", "species", "outbreak")
                   if k in present_traits]
    if has_year:
        filter_keys.append("year")

    display_defaults: dict = {
        "color_by": color_by,
        "distance_measure": "div",
        "layout": "rect",
    }
    if has_display:
        display_defaults["tip_label"] = "display_name"

    seq_desc = {"protein": "amino-acid", "nucleotide": "nucleotide"}.get(seq_type, seq_type)
    meta: dict = {
        "title": title or family,
        "updated": date.today().isoformat(),
        "build_url": "",
        "description": description or (
            f"vfam_trees v{__version__} divergence tree for *{family}*"
            + (f" ({seq_desc})" if seq_desc else "")
            + ".  Branch lengths are substitutions per site; this tree is not "
              "time-resolved."
        ),
        "maintainers": [{"name": "vfam_trees"}],
        "panels": ["tree"],
        "colorings": colorings,
        "display_defaults": display_defaults,
    }
    if filter_keys:
        meta["filters"] = filter_keys

    doc = {"version": "v2", "meta": meta, "tree": tree_root}

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(doc, f, indent=2)
    log.debug("Auspice JSON written to %s", output_json)
