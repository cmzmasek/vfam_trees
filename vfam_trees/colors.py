"""Leaf color assignment by taxonomic rank (genus, grouped by subfamily)."""
from __future__ import annotations

import colorsys


def assign_leaf_colors(
    sel_records,
    short_id_to_meta: dict,
    short_to_display: dict,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, list[str]]]:
    """Assign colors to tree leaves by genus, grouped by subfamily where available.

    Uses Python's built-in colorsys (HLS) — no external dependencies.
    Subfamilies get evenly-spaced base hues; genera within a subfamily are
    distinguished by varying lightness.

    Returns:
        display_to_color   : display_name → hex color  (for PDF/PNG label_colors)
        short_to_color     : short_id → hex color      (for PhyloXML <color>)
        genus_to_color     : genus → hex color          (for legend)
        subfamily_to_genera: subfamily label → [genus, ...]  (for structured legend)
    """
    # ------------------------------------------------------------------
    # Step 1: extract genus / subfamily / subgenus per leaf from lineage
    # ------------------------------------------------------------------
    leaf_taxa: dict[str, dict[str, str]] = {}
    for r in sel_records:
        meta = short_id_to_meta.get(r.id, {})
        lineage = meta.get("lineage_ranked", [])
        taxa: dict[str, str] = {"subfamily": "", "genus": "", "subgenus": ""}
        for entry in lineage:
            if not isinstance(entry, dict):
                continue
            rank = (entry.get("rank") or "").lower()
            if rank in taxa:
                taxa[rank] = entry.get("name", "")
        leaf_taxa[r.id] = taxa

    # ------------------------------------------------------------------
    # Step 2: group genera by subfamily
    # ------------------------------------------------------------------
    subfamily_genera: dict[str, set[str]] = {}
    unclassified_genera: set[str] = set()

    for taxa in leaf_taxa.values():
        genus = taxa["genus"]
        if not genus:
            continue
        subfamily = taxa["subfamily"]
        if subfamily:
            subfamily_genera.setdefault(subfamily, set()).add(genus)
        else:
            unclassified_genera.add(genus)

    n_groups = len(subfamily_genera) + (1 if unclassified_genera else 0)
    if n_groups == 0:
        return {}, {}, {}, {}

    # ------------------------------------------------------------------
    # Step 3: assign HLS colors — one hue band per subfamily group
    # ------------------------------------------------------------------
    genus_to_color: dict[str, str] = {}
    subfamily_to_genera: dict[str, list[str]] = {}

    groups = sorted(subfamily_genera.keys())
    if unclassified_genera:
        groups.append("")  # sentinel for unclassified

    hue_step = 1.0 / n_groups

    for group_idx, subfamily in enumerate(groups):
        base_hue = group_idx * hue_step
        if subfamily:
            genera_in_group = sorted(subfamily_genera[subfamily])
            subfamily_to_genera[subfamily] = genera_in_group
        else:
            genera_in_group = sorted(unclassified_genera)
            subfamily_to_genera["(unclassified)"] = genera_in_group

        n = len(genera_in_group)
        no_subfamilies = len(subfamily_genera) == 0
        for g_idx, genus in enumerate(genera_in_group):
            if no_subfamilies:
                # No subfamily level: spread genera across the full hue wheel
                # so each genus gets a visually distinct color.
                hue = g_idx / max(n, 1)
                r_f, g_f, b_f = colorsys.hls_to_rgb(hue, 0.45, 0.80)
            else:
                # Within a subfamily hue band: vary lightness 0.35–0.60
                lightness = 0.35 + 0.25 * (g_idx / max(n - 1, 1)) if n > 1 else 0.45
                r_f, g_f, b_f = colorsys.hls_to_rgb(base_hue, lightness, 0.80)
            genus_to_color[genus] = "#{:02x}{:02x}{:02x}".format(
                int(r_f * 255), int(g_f * 255), int(b_f * 255)
            )

    # ------------------------------------------------------------------
    # Step 4: map short IDs and display names to colors
    # ------------------------------------------------------------------
    default_color = "#888888"
    display_to_color: dict[str, str] = {}
    short_to_color: dict[str, str] = {}

    for short_id, taxa in leaf_taxa.items():
        color = genus_to_color.get(taxa["genus"], default_color)
        display_to_color[short_to_display.get(short_id, short_id)] = color
        short_to_color[short_id] = color

    return display_to_color, short_to_color, genus_to_color, subfamily_to_genera
