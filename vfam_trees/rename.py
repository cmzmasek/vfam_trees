"""Short ID assignment and restoration for sequences and trees."""
from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)


def assign_short_ids(
    records: list[SeqRecord],
    metadata: list[dict],
    family: str,
    id_map_path: Path,
) -> tuple[list[SeqRecord], dict[str, str]]:
    """Replace sequence IDs with short IDs and write the mapping table.

    Args:
        records: SeqRecords with original NCBI IDs
        metadata: list of metadata dicts (same order as records)
        family: family name used for prefix generation
        id_map_path: output path for id_map.tsv

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
            display_name = _build_display_name(meta)
            short_to_display[short_id] = display_name

            writer.writerow([short_id, rec.id, display_name])

            new_rec = rec[:]
            new_rec.id = short_id
            new_rec.description = ""
            renamed.append(new_rec)

    log.info("Assigned %d short IDs, map written to %s", len(renamed), id_map_path)
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


def restore_newick_names(input_nwk: Path, output_nwk: Path, id_map: dict[str, str]) -> None:
    """Replace short IDs with display names in a Newick file."""
    text = input_nwk.read_text()
    for short_id, display_name in id_map.items():
        text = text.replace(short_id, display_name)
    output_nwk.parent.mkdir(parents=True, exist_ok=True)
    output_nwk.write_text(text)
    log.debug("Restored names in %s → %s", input_nwk, output_nwk)


def _family_prefix(family: str) -> str:
    """Generate a short uppercase prefix from a family name."""
    consonants = [c for c in family.upper() if c.isalpha() and c not in "AEIOU"]
    return "".join(consonants[:4]) if len(consonants) >= 4 else family[:4].upper()


def _build_display_name(meta: dict) -> str:
    def _clean(s: str) -> str:
        return str(s).replace(" ", "_") if s else "unknown"

    parts = [
        _clean(meta.get("species")),
        _clean(meta.get("strain")),
        _clean(meta.get("accession")),
    ]
    host = _clean(meta.get("host"))
    if host != "unknown":
        parts.append(host)
    return "|".join(parts)
