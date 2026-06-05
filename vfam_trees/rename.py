"""Short ID assignment and restoration for sequences and trees."""
from __future__ import annotations

import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)

DEFAULT_LABEL_FORMAT = "{species}|{id}|{host}"

_FIELD_RE = re.compile(r'\{(\w+)\}')
_YEAR_RE  = re.compile(r'\b(1[89]\d{2}|20\d{2})\b')
_ABSENT   = frozenset({"", "unknown", "n/a"})

# Fields recognised in label format strings.
LABEL_FIELDS = frozenset({"species", "id", "host", "strain", "location", "year", "genus"})


def _normalize_field(val: str | None) -> str:
    """Return val, or empty string when it is None, blank, or a known placeholder."""
    if val is None:
        return ""
    s = str(val).strip()
    return "" if s.lower() in _ABSENT else s


def format_leaf_label(
    fmt: str,
    field_values: dict[str, str],
    replace_whitespace: bool = True,
    keep_separator_on_empty: bool = False,
) -> str:
    """Render a leaf-label format string substituting *field_values*.

    Supported placeholder names (``{name}``): ``species``, ``id``
    (accession), ``host``, ``strain``, ``location``, ``year``, ``genus``.
    Any literal text between placeholders — including separator characters
    such as ``|`` — is reproduced verbatim.

    *replace_whitespace* (default ``True``) — spaces in field values are
    replaced with underscores so Newick labels remain parseable.

    *keep_separator_on_empty* (default ``False``) — when a field value is
    empty (absent, ``"unknown"``, or ``"n/a"``, case-insensitive), its
    immediately preceding literal separator is also dropped so no
    consecutive or leading separators appear in the result.  Set ``True``
    to preserve separators regardless of field content.
    """
    parts = _FIELD_RE.split(fmt)
    n_fields = len(parts) // 2

    # Build (preceding_literal, resolved_value, original_index) triples.
    resolved: list[tuple[str, str, int]] = []
    for i in range(n_fields):
        preceding_lit = parts[2 * i]
        field_name    = parts[2 * i + 1]
        raw = _normalize_field(field_values.get(field_name, ""))
        val = raw.replace(" ", "_") if replace_whitespace else raw
        resolved.append((preceding_lit, val, i))
    trailing = parts[-1]

    if keep_separator_on_empty:
        return "".join(lit + val for lit, val, _ in resolved) + trailing

    # Drop empty fields and their preceding separators.
    non_empty = [(lit, val, orig_i) for lit, val, orig_i in resolved if val]
    if not non_empty:
        return ""
    result = ""
    for idx, (lit, val, orig_i) in enumerate(non_empty):
        if idx == 0 and orig_i > 0:
            # First non-empty field is not field-0: its preceding literal is
            # a separator from a skipped field — drop it.
            result += val
        else:
            result += lit + val
    return result + trailing


def _resolve_label_fields(meta: dict) -> dict[str, str]:
    """Extract all label-substitutable fields from a metadata dict.

    ``None``, blank, ``"unknown"``, and ``"n/a"`` (case-insensitive) all
    become the empty string so :func:`format_leaf_label` can cleanly omit
    absent fields when ``keep_separator_on_empty=False``.
    """
    raw_date = meta.get("collection_date", "") or ""
    year_m = _YEAR_RE.search(raw_date)
    return {
        "species":  _normalize_field(meta.get("species")),
        "id":       _normalize_field(meta.get("accession")),
        "host":     _normalize_field(meta.get("host")),
        "strain":   _normalize_field(meta.get("strain")),
        "location": _normalize_field(meta.get("location")),
        "year":     year_m.group(1) if year_m else "",
        "genus":    _normalize_field(meta.get("genus")),
    }


def canonical_leaf_label(
    species: str | None,
    accession: str | None,
    host: str | None,
) -> str:
    """Build a leaf label using the default ``{species}|{id}|{host}`` format.

    Delegates to :func:`format_leaf_label` with ``replace_whitespace=True``
    and ``keep_separator_on_empty=False``.  Absent / ``"unknown"`` / ``"n/a"``
    components are dropped rather than replaced with a placeholder.

    This is the single source of truth for the default label format used
    across single-protein and concat pipeline modes.
    """
    fields = {
        "species": _normalize_field(species),
        "id":      _normalize_field(accession),
        "host":    _normalize_field(host),
    }
    return format_leaf_label(DEFAULT_LABEL_FORMAT, fields)


def assign_short_ids(
    records: list[SeqRecord],
    metadata: list[dict],
    family: str,
    id_map_path: Path,
    label_format: str = DEFAULT_LABEL_FORMAT,
    replace_whitespace: bool = True,
    keep_separator_on_empty: bool = False,
) -> tuple[list[SeqRecord], dict[str, str]]:
    """Replace sequence IDs with short IDs and write the mapping table.

    Args:
        records: SeqRecords with original NCBI IDs
        metadata: list of metadata dicts (same order as records)
        family: family name used for prefix generation
        id_map_path: output path for id_map.tsv
        label_format: format string for leaf labels; default ``{species}|{id}|{host}``
        replace_whitespace: replace spaces with underscores in field values
        keep_separator_on_empty: keep separators adjacent to empty fields

    Returns:
        (renamed_records, short_to_display) where short_to_display maps
        short_id → full display name.
    """
    prefix = _family_prefix(family)
    short_to_display: dict[str, str] = {}
    renamed: list[SeqRecord] = []

    id_map_path.parent.mkdir(parents=True, exist_ok=True)

    with open(id_map_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["short_id", "accession", "display_name"])

        for i, (rec, meta) in enumerate(zip(records, metadata), start=1):
            short_id = f"{prefix}_{i:06d}"
            display_name = format_leaf_label(
                label_format,
                _resolve_label_fields(meta),
                replace_whitespace=replace_whitespace,
                keep_separator_on_empty=keep_separator_on_empty,
            )
            # Verbatim prefix for externally-sourced records (e.g. "PATHOPLEXUS_");
            # empty for everything else. Set by _inject_external_sequences.
            display_name = (meta.get("name_prefix") or "") + display_name
            short_to_display[short_id] = display_name

            writer.writerow([short_id, rec.id, display_name])

            new_rec = rec[:]
            new_rec.id = short_id
            new_rec.description = ""
            renamed.append(new_rec)

    log.info("Assigned %d short IDs", len(renamed))
    log.debug("id_map written to %s", id_map_path)
    return renamed, short_to_display


def load_id_map(id_map_path: Path) -> dict[str, str]:
    """Load short_id → display_name mapping from id_map.tsv."""
    mapping = {}
    with open(id_map_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["short_id"]] = row["display_name"]
    return mapping


def restore_fasta_names(input_fasta: Path, output_fasta: Path, id_map: dict[str, str]) -> None:
    """Replace short IDs with display names in a FASTA file."""
    records = list(SeqIO.parse(input_fasta, "fasta"))
    restored = []
    for rec in records:
        display = id_map.get(rec.id, rec.id)
        new_rec = rec[:]
        new_rec.id = display
        new_rec.description = ""
        restored.append(new_rec)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as f:
        SeqIO.write(restored, f, "fasta")
    log.debug("Restored names in %s → %s", input_fasta, output_fasta)


def _family_prefix(family: str) -> str:
    """Generate a short uppercase prefix from a family name."""
    consonants = [c for c in family.upper() if c.isalpha() and c not in "AEIOU"]
    return "".join(consonants[:4]) if len(consonants) >= 4 else family[:4].upper()


def _build_display_name(meta: dict) -> str:
    """Adaptor kept for test cross-checks — delegates to canonical_leaf_label."""
    return canonical_leaf_label(
        meta.get("species"),
        meta.get("accession"),
        meta.get("host"),
    )
