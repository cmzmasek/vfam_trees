"""Fetch forced-include sequences from external pathogen databases.

Currently only **Pathoplexus** (https://pathoplexus.org) is supported, queried
through its GenSpectrum LAPIS API.  Sequences returned here are injected into
the pipeline like ``manual.include`` records: they bypass QC and clustering and
are protected from outlier removal (see ``pipeline._inject_external_sequences``).

Pathoplexus organises data into a fixed set of curated *organisms* (each its own
LAPIS engine at ``https://lapis.pathoplexus.org/{organism}``), addressed by slug
— not by taxid — so a genus-level tree fans out to several organism entries
(e.g. Orthoebolavirus -> ``ebola-zaire`` + ``ebola-sudan`` + ``ebola-bdbv``).

Single source of truth for external fetching.  Uses only the standard library
(``urllib`` + ``json``) so it adds no compiled dependencies.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request

from .logger import get_logger

log = get_logger(__name__)

PATHOPLEXUS_BASE = "https://lapis.pathoplexus.org"

# The curated organism slugs Pathoplexus exposes (verified against the live
# instance).  Used for config validation and the per-organism LAPIS base URL.
KNOWN_PATHOPLEXUS_ORGANISMS = frozenset({
    "andv", "cchf", "dengue", "ebola-bdbv", "ebola-sudan", "ebola-zaire",
    "hmpv", "marburg", "measles", "mpox", "rsv-a", "rsv-b", "west-nile",
    "yellow-fever",
})

# Safety cap: these sequences bypass subsampling, so an unbounded outbreak
# query could dominate a tree.  Per-entry ``max_seqs`` overrides this.
DEFAULT_MAX_SEQS = 200

# Map config limiter keys -> LAPIS query parameters.  Note that the raw
# ``sampleCollectionDate`` field is a string (often partial) and is NOT
# range-queryable; the date-typed bounds ``sampleCollectionDateRange{Lower,Upper}``
# are, via the ``From``/``To`` suffix.
_FILTER_PARAMS = {
    "country":   "geoLocCountry",                      # exact match
    "host":      "hostNameScientific",                 # exact match
    "date_from": "sampleCollectionDateRangeLowerFrom",  # ISO YYYY-MM-DD
    "date_to":   "sampleCollectionDateRangeUpperTo",    # ISO YYYY-MM-DD
}

_USER_AGENT = "vfam_trees (https://github.com/; viral family tree pipeline)"
_HTTP_TIMEOUT = 60
_MAX_RETRIES = 5
_RETRY_DELAY = 10


def _http_get(url: str) -> bytes:
    """GET *url* with retries.  4xx fails fast (client error); 5xx/network retry."""
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    last_exc: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=_HTTP_TIMEOUT) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", "replace")[:500]
            if 400 <= e.code < 500:
                raise RuntimeError(
                    f"Pathoplexus/LAPIS request failed (HTTP {e.code}) for {url}: {body}"
                ) from e
            last_exc = e
            log.warning("LAPIS request attempt %d failed (HTTP %d): %s",
                        attempt + 1, e.code, body)
        except urllib.error.URLError as e:
            last_exc = e
            log.warning("LAPIS request attempt %d failed: %s", attempt + 1, e)
        if attempt < _MAX_RETRIES - 1:
            time.sleep(_RETRY_DELAY)
    raise RuntimeError(
        f"Pathoplexus/LAPIS request to {url} failed after {_MAX_RETRIES} attempts: {last_exc}"
    )


def _build_filters(entry: dict) -> dict[str, str]:
    """Translate a config entry's limiters into LAPIS query parameters."""
    params: dict[str, str] = {}
    for cfg_key, lapis_key in _FILTER_PARAMS.items():
        val = entry.get(cfg_key)
        if val:
            params[lapis_key] = str(val)
    return params


def _parse_fasta(text: str) -> dict[str, str]:
    """Parse FASTA text into {header_id: sequence}.  Header id is the first
    whitespace-delimited token after '>'."""
    out: dict[str, str] = {}
    cur: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur is not None:
                out[cur] = "".join(chunks)
            header = line[1:].strip()
            cur = header.split()[0] if header else ""
            chunks = []
        elif cur is not None:
            chunks.append(line.strip())
    if cur is not None:
        out[cur] = "".join(chunks)
    return out


def fetch_pathoplexus(entry: dict, *, ncbi_accessions: set[str]) -> list[dict]:
    """Fetch nucleotide sequences for one ``additional_data_sources`` entry.

    Returns a list of injection dicts with keys ``id, sequence, species, host,
    location, collection_date, strain, taxon_id, source`` — carrying *real*
    metadata (unlike pasted sequences which have none).  Records whose
    underlying INSDC accession is already in *ncbi_accessions* are dropped when
    ``dedup_vs_ncbi`` is set (default).  Raises ``RuntimeError`` on API failure.
    """
    organism = entry["organism"]
    base = f"{PATHOPLEXUS_BASE}/{organism}/sample"
    filters = _build_filters(entry)
    # Pathoplexus keeps every revision of an accession; restrict to the current
    # record so we don't inject several near-identical versions of one genome
    # (LATEST_VERSION already excludes revoked records).
    filters["versionStatus"] = "LATEST_VERSION"
    max_seqs = int(entry.get("max_seqs", DEFAULT_MAX_SEQS))
    dedup = entry.get("dedup_vs_ncbi", True)

    # Preflight count so the run log shows how many records match the filters.
    agg_url = f"{base}/aggregated?{urllib.parse.urlencode(filters)}"
    count = json.loads(_http_get(agg_url))["data"][0]["count"]
    log.info("Pathoplexus[%s]: %d sequence(s) match filters %s",
             organism, count, filters or "{}")
    if count == 0:
        return []
    if count > max_seqs:
        log.info("Pathoplexus[%s]: capping at max_seqs=%d (of %d matching)",
                 organism, max_seqs, count)

    # Identical filters + limit + orderBy on both endpoints so the metadata and
    # sequence pulls return the same records; join on accessionVersion.
    common = dict(filters)
    common["limit"] = str(max_seqs)
    common["orderBy"] = "accessionVersion"  # deterministic / reproducible

    details_url = f"{base}/details?{urllib.parse.urlencode({**common, 'dataFormat': 'json'})}"
    records = json.loads(_http_get(details_url))["data"]
    meta_by_id = {r["accessionVersion"]: r for r in records}

    seq_url = (f"{base}/unalignedNucleotideSequences?"
               f"{urllib.parse.urlencode({**common, 'dataFormat': 'fasta'})}")
    seqs = _parse_fasta(_http_get(seq_url).decode("utf-8"))

    norm_ncbi = {a.split(".")[0] for a in ncbi_accessions}
    out: list[dict] = []
    n_dedup = 0
    for acc, seq in seqs.items():
        rec = meta_by_id.get(acc)
        if rec is None or not seq:
            continue
        insdc = (rec.get("insdcAccessionBase") or "").strip()
        if dedup and insdc and insdc in norm_ncbi:
            n_dedup += 1
            continue
        date = (rec.get("sampleCollectionDate")
                or rec.get("sampleCollectionDateRangeUpper") or "")
        out.append({
            "id": acc,
            "sequence": seq.upper(),
            "species": (rec.get("ncbiVirusName") or "").strip(),
            "host": (rec.get("hostNameScientific") or "").strip(),
            "location": (rec.get("geoLocCountry") or "").strip(),
            "collection_date": str(date).strip(),
            "strain": "",
            "taxon_id": str(rec.get("ncbiVirusTaxId") or "").strip(),
            "source": "pathoplexus",
        })
    if n_dedup:
        log.info("Pathoplexus[%s]: skipped %d record(s) already present via NCBI "
                 "(insdcAccessionBase match)", organism, n_dedup)
    log.info("Pathoplexus[%s]: returning %d sequence(s)", organism, len(out))
    return out


# Registry of external sources.  Adding a source = one entry here + its fetcher.
SOURCE_FETCHERS = {
    "pathoplexus": fetch_pathoplexus,
}
