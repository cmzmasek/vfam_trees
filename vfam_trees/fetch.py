"""Fetch sequences and metadata from NCBI — per-species download."""
from __future__ import annotations

import io
import re
import time
from pathlib import Path
from typing import Generator, TYPE_CHECKING

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from .logger import get_logger

if TYPE_CHECKING:
    from .markers import MarkerIdentifier

log = get_logger(__name__)

FETCH_BATCH_SIZE = 200
RETRY_DELAY = 10
MAX_RETRIES = 5


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
            with Entrez.esearch(db="taxonomy", term=f'"{family}"[Scientific Name]') as handle:
                result = Entrez.read(handle)
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
            with Entrez.esearch(db="taxonomy", term=query, retmax=max_records) as handle:
                result = Entrez.read(handle)
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
                with Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml") as handle:
                    records = Entrez.read(handle)
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
                with Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml") as handle:
                    records = Entrez.read(handle)
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
        # Segmented virus: restrict to the specific segment by title.
        # Accept any common "complete" title phrasing — per-segment records are
        # sometimes "complete cds" or "complete genome" rather than the strict
        # "complete sequence".
        base += (
            f' AND "{segment}"[Title]'
            ' AND ("complete sequence"[Title] OR "complete cds"[Title]'
            ' OR "complete genome"[Title])'
        )
    elif region == "whole_genome":
        # Unsegmented whole-genome search: require complete genome/sequence in title
        if seq_type == "nucleotide":
            base += ' AND ("complete genome"[Title] OR "complete sequence"[Title])'
    else:
        # Marker gene / marker protein for large DNA viruses.
        # nuccore: search [Gene] — gene-name annotation on nucleotide records.
        # protein: search [Protein Name] — descriptive names like "DNA polymerase"
        #          are stored there, not in [Gene] (which holds short symbols like UL30).
        #          Also include [Gene] as a fallback for gene-symbol-named markers
        #          like B646L that may appear in either field.
        if seq_type == "protein":
            base += (
                f' AND ("{region}"[Protein Name] OR "{region}"[Gene])'
            )
        else:
            base += (
                f' AND "{region}"[Gene]'
                ' NOT "complete genome"[Title]'
                ' NOT "complete sequence"[Title]'
            )
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
            with Entrez.esearch(db=db, term=query, retmax=max_records) as handle:
                result = Entrez.read(handle)
            return result["IdList"]
        except Exception as e:
            log.warning("Search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    return []


def _fetch_batch(db: str, ids: list[str]) -> str:
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.efetch(db=db, id=",".join(ids), rettype="gb", retmode="text") as handle:
                data = handle.read()
            time.sleep(0.34)
            # Warn if NCBI returned fewer records than requested
            n_returned = data.count("\nLOCUS ") + (1 if data.startswith("LOCUS ") else 0)
            if n_returned < len(ids):
                log.warning(
                    "NCBI returned %d record(s) for a batch of %d requested IDs "
                    "(partial response)",
                    n_returned, len(ids),
                )
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


# ---------------------------------------------------------------------------
# Concatenation mode — multi-marker per-genome fetching
# CONCAT_DESIGN.md §5.2 (fetcher contract).  Phase 3 of the rollout.
# ---------------------------------------------------------------------------

def _safe_marker_filename(marker_name: str) -> str:
    """Produce a filename-safe slug from a marker name."""
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", marker_name).strip("_")
    return slug or "marker"


_ACCESSION_RE = re.compile(r"\b([A-Z]{1,3}_?\d+(?:\.\d+)?)\b")


def _source_nuc_accession(record: SeqRecord) -> str:
    """Return the source nucleotide accession for a protein record, or ''.

    Tries the CDS/Protein feature ``coded_by`` qualifier first (the most
    reliable source under modern GenBank conventions), then falls back to
    parsing the ``db_source`` annotation.
    """
    for feat in getattr(record, "features", None) or []:
        qualifiers = getattr(feat, "qualifiers", None) or {}
        for cb in qualifiers.get("coded_by", []):
            m = _ACCESSION_RE.search(cb)
            if m:
                return m.group(1)
    annotations = getattr(record, "annotations", None) or {}
    db_source = annotations.get("db_source", "") or ""
    m = re.search(r"accession\s+([A-Z]{1,3}_?\d+(?:\.\d+)?)", db_source, re.IGNORECASE)
    if m:
        return m.group(1)
    return ""


def _extract_isolate(record: SeqRecord) -> str:
    """Best-effort isolate / strain qualifier from the protein's source feature."""
    for feat in getattr(record, "features", None) or []:
        if getattr(feat, "type", None) != "source":
            continue
        qualifiers = getattr(feat, "qualifiers", None) or {}
        for key in ("isolate", "strain"):
            vals = qualifiers.get(key, [])
            if vals and isinstance(vals[0], str) and vals[0].strip():
                return vals[0].strip()
    return ""


def _count_split_submissions(bucket: dict[str, dict[str, list[SeqRecord]]]) -> int:
    """Count distinct isolate names that appear in 2+ source-nuc accessions.

    Diagnostic only: under Policy A each source-nuc accession is its own
    genome, so an isolate that's been split across multiple GenBank
    submissions surfaces as multiple low-coverage genomes.  This counter
    is a proxy for "how many specimens did Policy A's strict accession
    grouping cost us"; the value motivates whether to revisit Policy B
    in a future release.
    """
    isolate_to_accs: dict[str, set[str]] = {}
    for src_acc, marker_to_records in bucket.items():
        # One isolate per source_acc is enough — pick the first non-empty one.
        for records in marker_to_records.values():
            picked = False
            for rec in records:
                isolate = _extract_isolate(rec)
                if isolate:
                    isolate_to_accs.setdefault(isolate, set()).add(src_acc)
                    picked = True
                    break
            if picked:
                break
    return sum(1 for accs in isolate_to_accs.values() if len(accs) > 1)


def group_proteins_by_genome(
    proteins_by_marker: dict[str, list[SeqRecord]],
    marker_set: list[dict],
    species_lineage: list[dict] | None,
    min_fraction: float,
    identifier: "MarkerIdentifier | None" = None,
) -> tuple[dict[str, dict[str, SeqRecord]], dict]:
    """Group fetched proteins into genomes (Policy A) and apply min-fraction drop.

    Pure function — no network calls.  Tests exercise this directly.

    Args:
        proteins_by_marker: {marker_name: [SeqRecord]} from per-marker queries.
        marker_set:         the family's curated marker spec list.
        species_lineage:    NCBI ranked lineage; passed through to identifier
                            for subfamily-aware alias resolution.
        min_fraction:       genomes covering fewer than ceil(min_fraction × N)
                            markers are dropped.
        identifier:         MarkerIdentifier; defaults to NameMatchIdentifier().

    Returns:
        (genomes, stats):
            genomes: {source_nuc_accession: {marker_name: SeqRecord}}
            stats:   diagnostic counters (n_genomes_found, n_genomes_kept,
                     n_dropped_min_fraction, n_dropped_split_submission).
    """
    if identifier is None:
        from .markers import NameMatchIdentifier
        identifier = NameMatchIdentifier()

    # Group raw protein records by source nucleotide accession + marker name.
    bucket: dict[str, dict[str, list[SeqRecord]]] = {}
    n_orphaned_no_source = 0
    for marker_name, records in proteins_by_marker.items():
        for rec in records:
            src = _source_nuc_accession(rec)
            if not src:
                n_orphaned_no_source += 1
                continue
            bucket.setdefault(src, {}).setdefault(marker_name, []).append(rec)

    n_genomes_found = len(bucket)
    n_required_markers = max(1, int(round(min_fraction * len(marker_set))))

    genomes: dict[str, dict[str, SeqRecord]] = {}
    n_dropped_min_fraction = 0

    for src_acc, marker_to_records in bucket.items():
        # Flatten this genome's candidate proteins, then run the identifier
        # (which applies tiebreakers per marker against the marker_set).
        all_candidates: list[SeqRecord] = []
        for records in marker_to_records.values():
            all_candidates.extend(records)
        chosen = identifier.identify(all_candidates, marker_set, species_lineage)
        if len(chosen) >= n_required_markers:
            genomes[src_acc] = chosen
        else:
            n_dropped_min_fraction += 1

    stats = {
        "n_genomes_found":            n_genomes_found,
        "n_genomes_kept":             len(genomes),
        "n_dropped_min_fraction":     n_dropped_min_fraction,
        "n_dropped_split_submission": _count_split_submissions(bucket),
        "n_orphaned_no_source":       n_orphaned_no_source,
        "n_required_markers":         n_required_markers,
    }
    return genomes, stats


def fetch_species_genomes(
    taxid: int,
    species_name: str,
    species_lineage: list[dict],
    marker_set: list[dict],
    output_dir: Path,
    max_per_species: int = 200,
    min_fraction: float = 0.7,
    exclude_organisms: list[str] | None = None,
    identifier: "MarkerIdentifier | None" = None,
) -> tuple[dict[str, dict[str, SeqRecord]], dict]:
    """Concatenation-mode counterpart to ``fetch_species_sequences``.

    Issues one Entrez query per marker (RefSeq matches uncapped, non-RefSeq
    capped at ``max_per_species``), groups returned proteins by source
    nucleotide accession (Policy A), applies the per-genome paralog
    tiebreaker via ``identifier``, and drops genomes covering fewer than
    ``min_fraction × N`` markers.

    Returns:
        (genomes, stats):
            genomes: {source_nuc_accession: {marker_name: SeqRecord}}
            stats:   diagnostic counters (see group_proteins_by_genome).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    proteins_by_marker: dict[str, list[SeqRecord]] = {}
    n_proteins_fetched = 0

    for marker in marker_set:
        marker_name = marker["name"]
        marker_proteins = _fetch_marker_proteins(
            taxid=taxid,
            marker=marker,
            species_lineage=species_lineage,
            output_gb=output_dir / f"{_safe_marker_filename(marker_name)}.gb",
            max_per_species=max_per_species,
            exclude_organisms=exclude_organisms,
            species_name=species_name,
        )
        proteins_by_marker[marker_name] = marker_proteins
        n_proteins_fetched += len(marker_proteins)

    # Flag markers that returned 0 hits while sibling markers found something —
    # often a sign of annotation drift (alias list miss).  Quiet when the
    # species had no hits at all, since that's just an absent virus.
    n_markers_found = sum(1 for v in proteins_by_marker.values() if v)
    if n_markers_found > 0:
        for marker_name, recs in proteins_by_marker.items():
            if not recs:
                log.info(
                    "%s [marker=%s]: 0 protein records found — marker missing "
                    "for this species (sibling markers had hits)",
                    species_name, marker_name,
                )

    genomes, stats = group_proteins_by_genome(
        proteins_by_marker=proteins_by_marker,
        marker_set=marker_set,
        species_lineage=species_lineage,
        min_fraction=min_fraction,
        identifier=identifier,
    )
    stats["n_proteins_fetched"] = n_proteins_fetched
    # Count RefSeq genomes within this species' kept set (uses inline
    # accession-prefix detection — concat.is_refseq_genome would create an
    # import cycle; this stays in fetch.py).
    n_refseq_kept = sum(
        1 for gid in genomes
        if len(gid) >= 3 and gid[2] == "_" and gid[:2].isalpha() and gid[:2].isupper()
    )
    stats["n_refseq_kept"] = n_refseq_kept
    log.info(
        "%s: fetched %d proteins across %d markers; "
        "%d genome(s) kept (%d RefSeq), %d dropped (min_fraction)",
        species_name, n_proteins_fetched, len(marker_set),
        stats["n_genomes_kept"], n_refseq_kept,
        stats["n_dropped_min_fraction"],
    )
    return genomes, stats


def _fetch_marker_proteins(
    taxid: int,
    marker: dict,
    species_lineage: list[dict],
    output_gb: Path,
    max_per_species: int,
    exclude_organisms: list[str] | None,
    species_name: str = "",
) -> list[SeqRecord]:
    """Fetch all protein records for one marker within one species.

    Always retrieves RefSeq matches in full; tops up with non-RefSeq matches
    up to ``max_per_species``.  Writes a per-marker GenBank flat file and
    returns the parsed records.

    ``species_name`` is used only for log messages — the actual fetch is keyed
    by ``taxid``.
    """
    db = "protein"
    refseq_query = _build_marker_query(taxid, marker, species_lineage,
                                       exclude_organisms, refseq_only=True)
    refseq_ids = _search_ids(db, refseq_query, max_records=10_000)
    refseq_set = set(refseq_ids)

    n_remaining = max(0, max_per_species - len(refseq_ids))
    non_refseq_ids: list[str] = []
    if n_remaining > 0:
        all_query = _build_marker_query(taxid, marker, species_lineage,
                                        exclude_organisms)
        all_ids = _search_ids(db, all_query, max_records=max_per_species)
        non_refseq_ids = [i for i in all_ids if i not in refseq_set][:n_remaining]

    final_ids = refseq_ids + non_refseq_ids
    log.debug(
        "%s [marker=%s]: %d RefSeq + %d non-RefSeq protein records",
        species_name or f"taxid {taxid}", marker.get("name", "?"),
        len(refseq_ids), len(non_refseq_ids),
    )
    if not final_ids:
        return []

    output_gb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_gb, "w") as out_f:
        for batch in _batched(final_ids, FETCH_BATCH_SIZE):
            data = _fetch_batch(db, batch)
            out_f.write(data)

    return list(SeqIO.parse(output_gb, "genbank"))


def _build_marker_query(
    taxid: int,
    marker: dict,
    species_lineage: list[dict] | None,
    exclude_organisms: list[str] | None,
    refseq_only: bool = False,
) -> str:
    """Build an Entrez protein query for one marker within one species.

    Aliases (including any subfamily-specific aliases applicable to this
    species) are OR-combined as `[Protein Name]` clauses.
    """
    from .markers import _subfamily_from_lineage  # local import to avoid cycle
    subfamily = _subfamily_from_lineage(species_lineage)
    names = [marker["name"]] + list(marker.get("aliases", []))
    if subfamily:
        names.extend(marker.get(f"aliases_{subfamily}", []))

    or_clause = " OR ".join(f'"{n}"[Protein Name]' for n in names if n)
    base = f"txid{taxid}[Organism:exp] AND ({or_clause})"
    if refseq_only:
        base += " AND refseq[filter]"
    base += " NOT patent[filter]"
    for term in (exclude_organisms or []):
        base += f' NOT "{term}"[Organism]'
    return base
