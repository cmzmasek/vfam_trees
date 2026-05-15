"""Fetch sequences and metadata from NCBI — per-species download."""
from __future__ import annotations

import io
import re
import time
from pathlib import Path
from typing import Callable, Generator, Iterable, TYPE_CHECKING

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

    log.info("Family %s has taxid %d", family, family_taxid)

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
    last_exc: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.esearch(db="taxonomy", term=f'"{family}"[Scientific Name]') as handle:
                result = Entrez.read(handle)
            ids = result.get("IdList", [])
            return int(ids[0]) if ids else None
        except Exception as e:
            last_exc = e
            log.warning("Taxonomy search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    log.error(
        "Taxonomy search for family %r failed after %d attempts: %s",
        family, MAX_RETRIES, last_exc,
    )
    return None


def _taxonomy_search(query: str, max_records: int) -> list[str]:
    last_exc: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.esearch(db="taxonomy", term=query, retmax=max_records) as handle:
                result = Entrez.read(handle)
            return result.get("IdList", [])
        except Exception as e:
            last_exc = e
            log.warning("Taxonomy search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    log.error(
        "Taxonomy search %r failed after %d attempts — returning empty result: %s",
        query, MAX_RETRIES, last_exc,
    )
    return []


def fetch_species_taxids_under_lineages(entries: list[str]) -> dict[str, set[int]]:
    """For each entry (taxid digit-string or scientific name), return the
    set of species-rank descendant taxids under that taxon in NCBI taxonomy.

    Entries can be at any rank — species, genus, subfamily, family, etc.
    The returned mapping ``{entry: {species_taxids}}`` is empty for an entry
    whose name does not resolve, or whose taxon has no species-rank
    descendants.  Callers can use the empty sets to warn per-entry.
    """
    result: dict[str, set[int]] = {}
    for entry in entries:
        if entry.isdigit():
            taxid: int | None = int(entry)
        else:
            taxid = _get_family_taxid(entry)
            if taxid is None:
                result[entry] = set()
                continue
        ids = _taxonomy_search(
            f"txid{taxid}[subtree] AND species[rank]", max_records=50000,
        )
        result[entry] = {int(i) for i in ids}
    return result


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
    n_failed_batches = 0
    for batch in _batched(ids, 500):
        last_exc: Exception | None = None
        succeeded = False
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
                succeeded = True
                break
            except Exception as e:
                last_exc = e
                log.warning("Taxonomy lineage fetch attempt %d failed: %s", attempt + 1, e)
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY)
        if not succeeded:
            n_failed_batches += 1
            log.error(
                "Taxonomy lineage fetch for %d taxid(s) failed after %d attempts: %s — "
                "those lineages will be missing from downstream LCA voting",
                len(batch), MAX_RETRIES, last_exc,
            )
    if n_failed_batches:
        log.error(
            "fetch_taxonomy_lineages: %d batch(es) lost — %d/%d taxids resolved",
            n_failed_batches, len(result), len(ids),
        )
    return result


def _fetch_taxon_names(taxids: list[str]) -> list[dict]:
    """Fetch taxon names for a list of taxids, return list of {taxid, name}."""
    species = []
    n_failed_batches = 0
    for batch in _batched(taxids, 500):
        last_exc: Exception | None = None
        succeeded = False
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
                succeeded = True
                break
            except Exception as e:
                last_exc = e
                log.warning("Taxon name fetch attempt %d failed: %s", attempt + 1, e)
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY)
        if not succeeded:
            n_failed_batches += 1
            log.error(
                "Taxon name fetch for %d taxid(s) failed after %d attempts: %s — "
                "those species will be missing from the discovered list",
                len(batch), MAX_RETRIES, last_exc,
            )
    if n_failed_batches:
        log.error(
            "_fetch_taxon_names: %d batch(es) lost — %d/%d taxids resolved",
            n_failed_batches, len(species), len(taxids),
        )
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


def fetch_accessions_directly(accessions: Iterable[str]) -> list[SeqRecord]:
    """Fetch a set of accessions from NCBI, auto-detecting nuccore vs protein db.

    Accessions that match the nuccore shape are fetched from nuccore; all
    others are assumed to be protein accessions and fetched from the protein
    db.  Used to actively retrieve manual.include IDs that were not returned
    by the normal per-species download — for example a nuccore accession
    specified on a protein-mode run, or an accession from a species outside
    the family's NCBI taxonomy subtree.

    Returns parsed SeqRecords for every accession NCBI successfully returns.
    Unfetchable accessions are silently omitted; callers should compare the
    input and output sets to detect gaps.
    """
    accs = list(accessions)
    nuccore_accs = [a for a in accs if _is_nuccore_accession(a)]
    protein_accs = [a for a in accs if not _is_nuccore_accession(a)]

    records: list[SeqRecord] = []
    for db, db_accs in [("nuccore", nuccore_accs), ("protein", protein_accs)]:
        if not db_accs:
            continue
        data = _fetch_batch(db, db_accs)
        if data.strip():
            records.extend(SeqIO.parse(io.StringIO(data), "genbank"))
    return records


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
        # nuccore: search [Title] OR [Gene].  [Gene] is sparsely annotated for
        #   many viral families and uses symbol-level names (e.g. "GPC", "G1");
        #   [Title] reliably carries the descriptive name from the FASTA defline
        #   (e.g. "glycoprotein g1") so it must be the primary field.
        # protein: search [Protein Name] OR [Title] OR [Gene].
        #   [Protein Name] should hold descriptive names like "DNA polymerase",
        #   but a 2026-05 audit found this field is sparsely populated for
        #   many viral families (Marseilleviridae, Pithoviridae, Mimiviridae,
        #   Pandoraviridae, Iridoviridae, Phycodnaviridae) — the actual
        #   protein name is reliably in [Title] (the FASTA defline) instead.
        #   [Gene] catches gene-symbol-named markers like UL30 / B646L that
        #   may appear in either field.
        if seq_type == "protein":
            base += (
                f' AND ("{region}"[Protein Name]'
                f' OR "{region}"[Title]'
                f' OR "{region}"[Gene])'
            )
        else:
            base += (
                f' AND ("{region}"[Title]'
                f' OR "{region}"[Gene])'
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
    last_exc: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.esearch(db=db, term=query, retmax=max_records) as handle:
                result = Entrez.read(handle)
            return result["IdList"]
        except Exception as e:
            last_exc = e
            log.warning("Search attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    log.error(
        "Entrez search on %s failed after %d attempts — returning empty result. "
        "Query: %r. Last error: %s",
        db, MAX_RETRIES, query, last_exc,
    )
    return []


# Counts a record by matching the LOCUS keyword anchored at the start of a
# line.  Substring counting against "\nLOCUS " misses the first record and
# can also be inflated by "LOCUS " appearing inside REFERENCE / COMMENT blocks
# of long records.
_LOCUS_RE = re.compile(r"^LOCUS\s", re.MULTILINE)

# Strict nuccore-accession shape used to validate UIDs before sending them to
# Entrez esummary.  See _ACCESSION_RE for the rationale on the digit minimum.
# The three alternatives reflect real NCBI conventions:
#   - [A-Z]{1,4}_\d+       any underscore-separated (RefSeq-style; protein
#                          prefixes are rejected separately below)
#   - [A-Z]{1,2}\d+        bare GenBank nucleotide (1+5 or 2+6 forms)
#   - [A-Z]{4}\d+          WGS nucleotide (4+8/9/10 forms)
# Note that 3-letter, no-underscore accessions (e.g. ``ABO61246``) are
# *protein* GenBank accessions and are deliberately excluded — they have
# the same shape as a nuccore accession but resolve to the protein db,
# which makes esummary return an "Otherdb db=protein" batch-aborting error.
_NUCCORE_ACCESSION_RE = re.compile(
    r"^(?:[A-Z]{1,4}_\d{3,}|[A-Z]{1,2}\d{3,}|[A-Z]{4}\d{3,})(?:\.\d+)?$"
)

# Known RefSeq *protein* prefixes — any accession with one of these prefixes
# is the protein record's own accession, not a source nucleotide.  Sending
# one to nuccore esummary makes NCBI return an "Otherdb db=protein" error
# that aborts the whole batch.  These leak into _source_nuc_accession via
# the db_source fallback (typical RefSeq protein db_source is
# "REFSEQ: accession YP_009047263.1").
_PROTEIN_REFSEQ_PREFIXES = ("NP_", "XP_", "YP_", "AP_", "WP_", "ZP_", "ELP_")


def _is_nuccore_accession(acc: str) -> bool:
    """Return True if ``acc`` has a valid nuccore-accession shape and is not
    a known RefSeq protein prefix."""
    if not _NUCCORE_ACCESSION_RE.match(acc):
        return False
    if acc.startswith(_PROTEIN_REFSEQ_PREFIXES):
        return False
    return True


def fetch_nuc_lengths(accessions: Iterable[str]) -> dict[str, int]:
    """Return ``{accession: length_bp}`` for nucleotide accessions via esummary.

    Used in concat mode to identify partial single-gene submissions: each
    such submission has its own short source-nuc accession, and length-based
    filtering against the longest source-nuc in the family lets the pipeline
    keep only proteins drawn from genome-scale submissions.

    Failures or missing entries are silently omitted from the result so
    callers can decide whether to skip filtering when the lookup is
    incomplete.  Versioned accessions (``NC_002617.1``) and unversioned
    accessions (``NC_002617``) both work — esummary returns
    ``AccessionVersion`` and we key on the input string when present.
    """
    # Dedupe and drop anything that doesn't look like a nuccore accession —
    # NCBI rejects the *entire* batch if any single UID is malformed
    # (e.g. a "Q8" fragment that leaked in from a free-text qualifier),
    # which would otherwise erase the length map for every well-formed
    # accession in the same batch.
    accs: list[str] = []
    seen: set[str] = set()
    skipped: list[str] = []
    for a in accessions:
        if not a or a in seen:
            continue
        seen.add(a)
        if _is_nuccore_accession(a):
            accs.append(a)
        else:
            skipped.append(a)
    if skipped:
        log.debug(
            "fetch_nuc_lengths: skipping %d malformed accession(s): %s",
            len(skipped), skipped[:5],
        )
    result: dict[str, int] = {}
    for batch in _batched(accs, 200):
        result.update(_esummary_lengths_with_recovery(batch))
        time.sleep(0.34)  # courtesy throttle, matching efetch path
    return result


def _parse_summaries(summaries, batch_set: set[str]) -> dict[str, int]:
    """Pull (accession → length) entries out of an esummary response."""
    out: dict[str, int] = {}
    for s in summaries:
        try:
            length = int(s.get("Length", 0))
        except (TypeError, ValueError):
            continue
        if length <= 0:
            continue
        for key in (s.get("AccessionVersion"), s.get("Caption")):
            if key and key in batch_set:
                out[key] = length
                break
    return out


# Substrings that mark *deterministic* NCBI esummary errors — retrying won't
# help, the only recovery is to remove the offending UID from the batch.
_DETERMINISTIC_ESUMMARY_ERRORS = ("Otherdb", "Invalid uid")


def _esummary_lengths_with_recovery(batch: list[str]) -> dict[str, int]:
    """Run esummary on one batch with retry+split error recovery.

    A single bad UID (e.g. one that resolves into the protein db, or a
    malformed fragment that slipped past the shape filter) makes NCBI
    reject the entire batch.  Retrying changes nothing for these
    deterministic errors, so on detection we fall through to a binary
    split: each half retries independently, recursively isolating the
    bad UID.  Singletons that still fail are logged at debug and
    omitted from the length map (the caller treats missing entries
    as "unknown length, do not filter").
    """
    if not batch:
        return {}
    last_exc: Exception | None = None
    summaries = None
    deterministic = False
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.esummary(db="nuccore", id=",".join(batch)) as handle:
                summaries = Entrez.read(handle)
            break
        except Exception as e:
            last_exc = e
            err_str = str(e)
            if any(tok in err_str for tok in _DETERMINISTIC_ESUMMARY_ERRORS):
                deterministic = True
                break
            log.warning("nuc-length esummary attempt %d failed: %s",
                        attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    if summaries is not None:
        return _parse_summaries(summaries, set(batch))
    if len(batch) == 1:
        log.debug(
            "nuc-length: dropping unresolvable accession %r (%s)",
            batch[0], last_exc,
        )
        return {}
    if not deterministic:
        log.error(
            "nuc-length esummary for %d accession(s) failed after %d "
            "attempts — splitting batch and retrying. Last error: %s",
            len(batch), MAX_RETRIES, last_exc,
        )
    mid = len(batch) // 2
    out = _esummary_lengths_with_recovery(batch[:mid])
    out.update(_esummary_lengths_with_recovery(batch[mid:]))
    return out


def _fetch_batch(db: str, ids: list[str]) -> str:
    last_exc: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.efetch(db=db, id=",".join(ids), rettype="gb", retmode="text") as handle:
                data = handle.read()
            time.sleep(0.34)
            # Warn if NCBI returned fewer records than requested.
            n_returned = len(_LOCUS_RE.findall(data))
            if n_returned < len(ids):
                log.warning(
                    "NCBI returned %d record(s) for a batch of %d requested IDs "
                    "(partial response)",
                    n_returned, len(ids),
                )
            return data
        except Exception as e:
            last_exc = e
            log.warning("Fetch attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    log.error(
        "Entrez efetch on %s failed after %d attempts for %d ID(s) — "
        "returning empty payload. Last error: %s",
        db, MAX_RETRIES, len(ids), last_exc,
    )
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


# Nuccore-shape accession matcher used inside _source_nuc_accession's text
# scans (coded_by qualifiers, db_source).  Mirrors _NUCCORE_ACCESSION_RE but
# anchored on word boundaries instead of full-string anchors.  Excludes the
# 3-letter no-underscore shape (GenBank protein) — see _NUCCORE_ACCESSION_RE
# for rationale.
_ACCESSION_RE = re.compile(
    r"\b((?:[A-Z]{1,4}_\d{3,}|[A-Z]{1,2}\d{3,}|[A-Z]{4}\d{3,})(?:\.\d+)?)\b"
)

# Stricter matcher for the db_source fallback path: only known RefSeq
# *nucleotide* prefixes are accepted.  An unrestricted match here would also
# pull in the protein record's own home (RefSeq protein YP_/NP_, UniProt
# P12345, GenBank protein ABO61246).
_DB_SOURCE_NUC_RE = re.compile(
    r"accession\s+((?:NC|NG|NM|NR|NS|NT|NW|NZ|AC|XM|XR)_\d{3,}(?:\.\d+)?)",
    re.IGNORECASE,
)


def _source_nuc_accession(record: SeqRecord) -> str:
    """Return the source nucleotide accession for a protein record, or ''.

    Tries the CDS/Protein feature ``coded_by`` qualifier first (the most
    reliable source under modern GenBank conventions), then falls back to
    parsing the ``db_source`` annotation — but the fallback is restricted
    to known RefSeq *nucleotide* prefixes only.  ``db_source`` typically
    describes the *protein record's own* home (``YP_…``, ``NP_…``,
    UniProt ``P12345``, GenBank ``ABO61246``), and an unrestricted shape
    match there silently hands those protein accessions to nuccore
    esummary downstream — making NCBI return "Otherdb db=protein" errors
    that abort the entire batch.
    """
    for feat in getattr(record, "features", None) or []:
        qualifiers = getattr(feat, "qualifiers", None) or {}
        for cb in qualifiers.get("coded_by", []):
            for m in _ACCESSION_RE.finditer(cb):
                cand = m.group(1)
                if not cand.startswith(_PROTEIN_REFSEQ_PREFIXES):
                    return cand
    annotations = getattr(record, "annotations", None) or {}
    db_source = annotations.get("db_source", "") or ""
    m = _DB_SOURCE_NUC_RE.search(db_source)
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
    source_nuc_min_length_frac: float = 0.0,
    nuc_length_lookup: Callable[[Iterable[str]], dict[str, int]] | None = None,
) -> tuple[dict[str, dict[str, SeqRecord]], dict]:
    """Group fetched proteins into genomes (Policy A) and apply min-fraction drop.

    Args:
        proteins_by_marker: {marker_name: [SeqRecord]} from per-marker queries.
        marker_set:         the family's curated marker spec list.
        species_lineage:    NCBI ranked lineage; passed through to identifier
                            for subfamily-aware alias resolution.
        min_fraction:       genomes covering fewer than ceil(min_fraction × N)
                            markers are dropped.
        identifier:         MarkerIdentifier; defaults to NameMatchIdentifier().
        source_nuc_min_length_frac:  if > 0, drop source-nuc accessions whose
                            parent nucleotide record is shorter than this
                            fraction of the longest source-nuc parent in the
                            current bucket set.  Filters single-gene partial
                            submissions whose entire "genome" is one CDS.
                            0 disables (default — pure / no network).
        nuc_length_lookup:  callable that maps an iterable of accessions to
                            ``{accession: length_bp}``.  Required when
                            ``source_nuc_min_length_frac > 0``.  Defaults to
                            :func:`fetch_nuc_lengths` (an Entrez esummary call).
                            Tests inject a stub.

    Returns:
        (genomes, stats):
            genomes: {source_nuc_accession: {marker_name: SeqRecord}}
            stats:   diagnostic counters (n_genomes_found, n_genomes_kept,
                     n_dropped_min_fraction, n_dropped_split_submission,
                     n_dropped_short_source_nuc).

    The function is pure when ``source_nuc_min_length_frac <= 0`` or a
    non-network ``nuc_length_lookup`` is injected; the default network path
    only fires when the caller opts in.
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

    # Source-nuc length filter — drops single-gene partial submissions whose
    # parent nucleotide record is much shorter than the longest in this set.
    n_dropped_short_source_nuc = 0
    if source_nuc_min_length_frac > 0 and bucket:
        if nuc_length_lookup is None:
            nuc_length_lookup = fetch_nuc_lengths
        try:
            lengths = nuc_length_lookup(set(bucket.keys()))
        except Exception as exc:
            log.warning(
                "Source-nuc length lookup failed (%s) — skipping length filter.", exc,
            )
            lengths = {}
        if lengths:
            longest = max(lengths.values())
            min_len = source_nuc_min_length_frac * longest
            short_accs = [acc for acc in bucket if lengths.get(acc, 0) < min_len]
            for acc in short_accs:
                bucket.pop(acc, None)
            n_dropped_short_source_nuc = len(short_accs)
            if short_accs:
                log.info(
                    "Source-nuc length filter dropped %d / %d accession(s) shorter "
                    "than %.0f bp (< %.0f%% of longest %d bp); kept %d genome(s) "
                    "for marker assembly",
                    len(short_accs), len(short_accs) + len(bucket),
                    min_len, source_nuc_min_length_frac * 100, longest, len(bucket),
                )

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
        "n_genomes_found":              n_genomes_found,
        "n_genomes_kept":               len(genomes),
        "n_dropped_min_fraction":       n_dropped_min_fraction,
        "n_dropped_split_submission":   _count_split_submissions(bucket),
        "n_dropped_short_source_nuc":   n_dropped_short_source_nuc,
        "n_orphaned_no_source":         n_orphaned_no_source,
        "n_required_markers":           n_required_markers,
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
    source_nuc_min_length_frac: float = 0.0,
    nuc_length_lookup: Callable[[Iterable[str]], dict[str, int]] | None = None,
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
        source_nuc_min_length_frac=source_nuc_min_length_frac,
        nuc_length_lookup=nuc_length_lookup,
    )
    from .concat import is_refseq_genome  # local to keep fetch standalone-importable
    stats["n_proteins_fetched"] = n_proteins_fetched
    n_refseq_kept = sum(1 for gid in genomes if is_refseq_genome(gid))
    stats["n_refseq_kept"] = n_refseq_kept
    log.info(
        "%s: fetched %d proteins across %d markers; "
        "%d genome(s) kept (%d RefSeq), %d dropped (min_fraction), "
        "%d dropped (short source-nuc)",
        species_name, n_proteins_fetched, len(marker_set),
        stats["n_genomes_kept"], n_refseq_kept,
        stats["n_dropped_min_fraction"],
        stats.get("n_dropped_short_source_nuc", 0),
    )
    return genomes, stats


def load_proteins_from_marker_dir(
    marker_dir: Path,
    marker_set: list[dict],
) -> tuple[dict[str, list[SeqRecord]], int]:
    """Parse cached per-marker GenBank files into ``proteins_by_marker``.

    Counterpart to the per-marker fetch loop in ``fetch_species_genomes`` —
    used when the global sequence cache hits and we want to skip the network
    call.  Looks for ``<safe_marker>.gb`` files keyed by the same naming
    scheme the fetcher uses.

    Returns:
        (proteins_by_marker, n_total_proteins)
    """
    proteins_by_marker: dict[str, list[SeqRecord]] = {}
    n_total = 0
    for marker in marker_set:
        marker_name = marker["name"]
        gb_path = marker_dir / f"{_safe_marker_filename(marker_name)}.gb"
        if gb_path.exists() and gb_path.stat().st_size > 0:
            recs = list(SeqIO.parse(gb_path, "genbank"))
        else:
            recs = []
        proteins_by_marker[marker_name] = recs
        n_total += len(recs)
    return proteins_by_marker, n_total


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
    log.info(
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

    # Search [Protein Name] OR [Title] OR [Gene]:
    #   [Protein Name] should hold descriptive names like "DNA polymerase",
    #     but a 2026-05 audit found this field is sparsely populated for
    #     several viral families (Marseilleviridae, Pithoviridae, Mimiviridae,
    #     Pandoraviridae, Iridoviridae, Phycodnaviridae) — the actual name
    #     is reliably in [Title] (the FASTA defline) instead.
    #   [Gene] catches gene-symbol-named markers (B646L, UL30, ...) that may
    #     live in either field.
    or_terms = []
    for n in names:
        if not n:
            continue
        or_terms.append(f'"{n}"[Protein Name]')
        or_terms.append(f'"{n}"[Title]')
        or_terms.append(f'"{n}"[Gene]')
    or_clause = " OR ".join(or_terms)
    base = f"txid{taxid}[Organism:exp] AND ({or_clause})"
    if refseq_only:
        base += " AND refseq[filter]"
    base += " NOT patent[filter]"
    for term in (exclude_organisms or []):
        base += f' NOT "{term}"[Organism]'
    return base
