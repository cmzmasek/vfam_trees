"""Fetch sequences and metadata from NCBI — per-species download."""
from __future__ import annotations

import time
from pathlib import Path
from typing import Generator

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

log = get_logger(__name__)

FETCH_BATCH_SIZE = 200
RETRY_DELAY = 5
MAX_RETRIES = 3


def configure_entrez(email: str, api_key: str | None = None) -> None:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


# ---------------------------------------------------------------------------
# Taxonomy discovery
# ---------------------------------------------------------------------------

def discover_species(family: str) -> list[dict]:
    """Return all species-rank taxa under a viral family from NCBI taxonomy.

    Args:
        family: ICTV family name (e.g. 'Flaviviridae')

    Returns:
        List of dicts with keys 'taxid' (int) and 'name' (str).
        Empty list if family not found or no species returned.
    """
    log.info("Discovering species in %s from NCBI taxonomy...", family)

    # Step 1: find the family's taxid
    family_taxid = _get_family_taxid(family)
    if family_taxid is None:
        log.warning("Could not find NCBI taxonomy entry for family: %s", family)
        return []

    log.debug("Family %s has taxid %d", family, family_taxid)

    # Step 2: get all species-rank descendants
    query = f"txid{family_taxid}[subtree] AND species[rank]"
    ids = _taxonomy_search(query, max_records=50000)
    if not ids:
        log.warning("No species found under %s (taxid %d)", family, family_taxid)
        return []

    # Step 3: fetch names for all taxids
    species = _fetch_taxon_names(ids)

    # Step 4: remove environmental samples and unclassified entries
    before = len(species)
    species = _filter_species(species, family)
    removed = before - len(species)
    if removed:
        log.info("Removed %d environmental/unclassified entries", removed)

    log.info("Found %d species in %s", len(species), family)
    return species


def get_family_taxid(family: str) -> int | None:
    """Return the NCBI taxonomy ID for a viral family name, or None if not found."""
    return _get_family_taxid(family)


def _get_family_taxid(family: str) -> int | None:
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.esearch(db="taxonomy", term=f'"{family}"[Scientific Name]')
            result = Entrez.read(handle)
            handle.close()
            ids = result.get("IdList", [])
            return int(ids[0]) if ids else None
        except Exception as e:
            log.warning("Taxonomy search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    return None


def _taxonomy_search(query: str, max_records: int) -> list[str]:
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.esearch(db="taxonomy", term=query, retmax=max_records)
            result = Entrez.read(handle)
            handle.close()
            return result.get("IdList", [])
        except Exception as e:
            log.warning("Taxonomy search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    return []


def _filter_species(species: list[dict], family: str) -> list[dict]:
    """Remove environmental samples and unclassified entries."""
    family_lower = family.lower()
    filtered = []
    for sp in species:
        name_lower = sp["name"].lower()
        if "environmental samples" in name_lower:
            log.debug("Skipping environmental sample: %s", sp["name"])
            continue
        if f"unclassified {family_lower}" in name_lower:
            log.debug("Skipping unclassified entry: %s", sp["name"])
            continue
        filtered.append(sp)
    return filtered


def fetch_taxonomy_lineages(taxids) -> dict[str, list[dict]]:
    """Fetch ranked lineages for a set of NCBI taxids.

    Returns {taxid: [{"name": str, "rank": str}, ...]} ordered root → tip,
    with the taxon itself appended as the final entry. Ranks are NCBI's
    authoritative values (including "no rank" / "clade").
    """
    ids = [str(t) for t in taxids if t]
    result: dict[str, list[dict]] = {}
    for batch in _batched(ids, 500):
        for attempt in range(MAX_RETRIES):
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                for rec in records:
                    taxid = str(rec["TaxId"])
                    lineage = [
                        {"name": str(e["ScientificName"]), "rank": str(e.get("Rank", "") or "")}
                        for e in rec.get("LineageEx", [])
                    ]
                    lineage.append({
                        "name": str(rec["ScientificName"]),
                        "rank": str(rec.get("Rank", "") or ""),
                    })
                    result[taxid] = lineage
                time.sleep(0.2)
                break
            except Exception as e:
                log.warning("Taxonomy lineage fetch attempt %d failed: %s", attempt + 1, e)
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY)
    return result


def _fetch_taxon_names(taxids: list[str]) -> list[dict]:
    """Fetch taxon names for a list of taxids, return list of {taxid, name}."""
    species = []
    for batch in _batched(taxids, 500):
        for attempt in range(MAX_RETRIES):
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                for rec in records:
                    species.append({
                        "taxid": int(rec["TaxId"]),
                        "name": rec["ScientificName"],
                    })
                time.sleep(0.2)
                break
            except Exception as e:
                log.warning("Taxon name fetch attempt %d failed: %s", attempt + 1, e)
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY)
    return species


# ---------------------------------------------------------------------------
# Per-species sequence download
# ---------------------------------------------------------------------------

def fetch_species_sequences(
    taxid: int,
    species_name: str,
    seq_type: str,
    region: str,
    output_gb: Path,
    max_per_species: int = 200,
    exclude_organisms: list[str] | None = None,
    segment: str | None = None,
) -> int:
    """Fetch GenBank records for a single species.

    RefSeqs are always included in full regardless of max_per_species.
    Remaining slots (up to max_per_species) are filled with non-RefSeq entries,
    with RefSeqs ordered first.

    Args:
        taxid: NCBI taxonomy ID for the species
        species_name: human-readable species name (for logging)
        seq_type: 'nucleotide' or 'protein'
        region: 'whole_genome' or a gene/segment name
        output_gb: path to write GenBank flat file
        max_per_species: cap on non-RefSeq sequences (RefSeqs are uncapped)

    Returns:
        Number of records fetched, or 0 if none found.
    """
    db = "nuccore" if seq_type == "nucleotide" else "protein"

    # Step 1: always fetch all RefSeq records (no cap)
    refseq_query = _build_species_query(taxid, seq_type, region, exclude_organisms, refseq_only=True, segment=segment)
    refseq_ids = _search_ids(db, refseq_query, max_records=10_000)
    refseq_set = set(refseq_ids)
    log.debug("%s: found %d RefSeq record(s)", species_name, len(refseq_ids))

    # Step 2: fill remaining slots with non-RefSeq entries
    n_remaining = max(0, max_per_species - len(refseq_ids))
    non_refseq_ids: list[str] = []
    if n_remaining > 0:
        all_query = _build_species_query(taxid, seq_type, region, exclude_organisms, segment=segment)
        all_ids = _search_ids(db, all_query, max_records=max_per_species)
        non_refseq_ids = [i for i in all_ids if i not in refseq_set][:n_remaining]

    final_ids = refseq_ids + non_refseq_ids
    if not final_ids:
        log.debug("No sequences found for %s", species_name)
        return 0

    log.debug("Fetching %s (taxid %d): %d RefSeq + %d other",
              species_name, taxid, len(refseq_ids), len(non_refseq_ids))

    output_gb.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with open(output_gb, "w") as out_f:
        for batch in _batched(final_ids, FETCH_BATCH_SIZE):
            data = _fetch_batch(db, batch)
            out_f.write(data)
            total += len(batch)

    log.debug("Fetched %d sequences for %s", total, species_name)
    return total


def _build_species_query(
    taxid: int,
    seq_type: str,
    region: str,
    exclude_organisms: list[str] | None = None,
    refseq_only: bool = False,
    segment: str | None = None,
) -> str:
    base = f"txid{taxid}[Organism:exp]"
    if segment:
        # Segmented virus: restrict to the specific segment by title
        base += f' AND "{segment}"[Title] AND "complete sequence"[Title]'
    elif region == "whole_genome":
        # Unsegmented whole-genome search: require complete genome/sequence in title
        if seq_type == "nucleotide":
            base += ' AND ("complete genome"[Title] OR "complete sequence"[Title])'
    else:
        # Marker gene (e.g. large DNA viruses): filter by gene name only —
        # individual gene records rarely carry "complete genome" in their title
        base += f' AND "{region}"[Gene]'
    if refseq_only:
        base += " AND refseq[filter]"
    base += " NOT patent[filter]"
    for term in (exclude_organisms or []):
        base += f' NOT "{term}"[Organism]'
    return base


# ---------------------------------------------------------------------------
# Parsing + metadata
# ---------------------------------------------------------------------------

def parse_gb_records(gb_file: Path) -> Generator[SeqRecord, None, None]:
    """Yield SeqRecord objects parsed from a GenBank flat file."""
    with open(gb_file) as f:
        yield from SeqIO.parse(f, "genbank")


def extract_metadata(record: SeqRecord) -> dict:
    """Extract metadata fields from a GenBank SeqRecord."""
    meta = {
        "accession": record.id,
        "seq_name": record.description or "",
        "species": "",
        "strain": "",
        "host": "unknown",
        "collection_date": "unknown",
        "location": "unknown",
        "taxon_id": "",
        "length": len(record.seq),
    }

    meta["species"] = record.annotations.get("organism", "").strip() or "unknown"
    meta["lineage"] = record.annotations.get("taxonomy", [])

    for feature in record.features:
        if feature.type == "source":
            q = feature.qualifiers
            meta["strain"] = _first(q.get("strain", [])) or _first(q.get("isolate", [])) or "unknown"
            meta["host"] = _first(q.get("host", [])) or _first(q.get("lab_host", [])) or "unknown"
            meta["collection_date"] = _first(q.get("collection_date", [])) or "unknown"
            meta["location"] = _first(q.get("country", [])) or "unknown"
            db_xref = _first(q.get("db_xref", []))
            if db_xref and db_xref.startswith("taxon:"):
                meta["taxon_id"] = db_xref.split("taxon:")[1]
            break

    return meta


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _search_ids(db: str, query: str, max_records: int) -> list[str]:
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.esearch(db=db, term=query, retmax=max_records)
            result = Entrez.read(handle)
            handle.close()
            return result["IdList"]
        except Exception as e:
            log.warning("Search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    return []


def _fetch_batch(db: str, ids: list[str]) -> str:
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db=db, id=",".join(ids), rettype="gb", retmode="text")
            data = handle.read()
            handle.close()
            time.sleep(0.34)
            return data
        except Exception as e:
            log.warning("Fetch attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    return ""


def _first(lst: list) -> str:
    return lst[0] if lst else ""


def _batched(items: list, size: int) -> Generator[list, None, None]:
    for i in range(0, len(items), size):
        yield items[i : i + size]
