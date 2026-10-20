"""Per-leaf fields derived from a metadata dict, shared by the tree writers.

Single source of truth for how the PhyloXML and Auspice writers turn a leaf's
metadata into: its rank-placed taxa (ICTV species, genus, ...), its country,
and its collection date as a decimal year.
"""
from __future__ import annotations

import calendar
import re
from datetime import date

# "unknown" is what extract_metadata emits for an absent qualifier; the rest
# are free-text placeholders and the INSDC missing-value vocabulary ("not
# collected", "missing: lab stock", ...), which GenBank now uses in
# /geo_loc_name, /host and /collection_date.
_PLACEHOLDERS = frozenset({
    "", "unknown", "n/a", "na", "none",
    "missing", "not applicable", "not collected", "not provided",
    "restricted access",
})


def clean_value(val) -> str:
    """A stripped string, or '' when the value is empty or a placeholder.

    Single source of truth for "no value" in the PhyloXML and Auspice writers.
    """
    s = str(val).strip() if val is not None else ""
    low = s.lower()
    if low in _PLACEHOLDERS or low.startswith("missing:"):
        return ""
    return s



# ---------------------------------------------------------------------------
# Taxa
# ---------------------------------------------------------------------------

def ranked_name(meta: dict, rank: str) -> str:
    """Name of the leaf's taxon at ``rank`` in its NCBI ranked lineage, or ''."""
    name = ""
    for entry in meta.get("lineage_ranked") or []:
        if isinstance(entry, dict) and (entry.get("rank") or "").lower() == rank:
            name = entry.get("name", "") or ""
    return clean_value(name)


def leaf_species(meta: dict) -> str:
    """The leaf's species, consistent with the internal LCA labels.

    When the leaf has a ranked lineage, this is its species-rank entry (the
    ICTV species, e.g. "Orthoebolavirus zairense"), and '' when the lineage
    has no species rank (an unclassified virus is left unplaced rather than
    given a pseudo-species such as "Sapovirus sp.").  Without a ranked lineage
    (taxonomy lookup failed, pasted sequence) it falls back to the organism.
    """
    if meta.get("lineage_ranked"):
        return ranked_name(meta, "species")
    return clean_value(meta.get("species", ""))


def own_taxon(meta: dict) -> dict:
    """NCBI's entry for the leaf's own taxid ({"name", "rank", "taxid"}), or {}.

    The last entry of a ranked lineage is the queried taxon itself, but only
    when its taxid matches the leaf's (in concat mode the lineage is the
    species', not the genome's own taxon).
    """
    taxon_id = str(meta.get("taxon_id") or "").strip()
    lineage = meta.get("lineage_ranked") or []
    if not taxon_id or not lineage:
        return {}
    last = lineage[-1]
    if isinstance(last, dict) and str(last.get("taxid") or "") == taxon_id:
        return last
    return {}


# ---------------------------------------------------------------------------
# Location
# ---------------------------------------------------------------------------

def location_fields(meta: dict) -> tuple[str, str]:
    """Return (country, detail) from an INSDC geo_loc_name ("Country: region").

    ``detail`` is the full string, and is '' when it says no more than the
    country, so a writer emits a separate Location only when it adds detail.
    """
    location = clean_value(meta.get("location", ""))
    if not location:
        return "", ""
    country = location.split(":", 1)[0].strip()
    detail = location if location != country else ""
    return country, detail


# ---------------------------------------------------------------------------
# Dates
# ---------------------------------------------------------------------------

_MONTHS = {m.lower(): i for i, m in enumerate(calendar.month_abbr) if m}

_ISO_RE = re.compile(r"^(\d{4})(?:-(\d{1,2})(?:-(\d{1,2}))?)?(?:T.*)?$")
_INSDC_RE = re.compile(r"^(?:(\d{1,2})-)?([A-Za-z]{3})-(\d{4})$")


def _interval(text: str) -> tuple[float, float] | None:
    """[start, end) in decimal years for one date, or None if unparsable."""
    text = text.strip()
    m = _ISO_RE.match(text)
    if m:
        year, month, day = int(m.group(1)), m.group(2), m.group(3)
        month = int(month) if month else 0
        day = int(day) if day else 0
    else:
        m = _INSDC_RE.match(text)
        if not m:
            return None
        month = _MONTHS.get(m.group(2).lower(), 0)
        if not month:
            return None
        year, day = int(m.group(3)), int(m.group(1) or 0)
    if not 1000 <= year <= 2999:
        return None
    n_days = 366 if calendar.isleap(year) else 365
    if not month:
        return float(year), float(year + 1)
    if not 1 <= month <= 12:
        return None
    first = date(year, month, 1).timetuple().tm_yday - 1
    if not day:
        length = calendar.monthrange(year, month)[1]
        return year + first / n_days, year + (first + length) / n_days
    try:
        doy = date(year, month, day).timetuple().tm_yday
    except ValueError:
        return None
    return year + (doy - 1) / n_days, year + doy / n_days


def decimal_year(collection_date: str) -> float | None:
    """Decimal year for a GenBank / ISO collection date, or None.

    Accepts "2014", "Mar-2014", "12-Mar-2014", "2014-03", "2014-03-12" (an
    ISO time suffix is ignored) and ranges "A/B".  The value is the midpoint
    of the stated interval, leap-year aware: a year is year + 0.5, a month
    its middle, a day its middle.  This is the convention both Archaeopteryx
    and Archaeopteryx.js use for dates read from tip names, so a date written
    here agrees with the one they would derive.
    """
    text = clean_value(collection_date)
    if not text:
        return None
    parts = text.split("/")
    if len(parts) > 2:
        return None
    bounds = [_interval(p) for p in parts]
    if any(b is None for b in bounds):
        return None
    start, end = bounds[0][0], bounds[-1][1]
    if end <= start:
        return None
    return (start + end) / 2
