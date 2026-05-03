"""Global file-based cache for per-species GenBank sequence downloads.

Cache layout
------------
{cache_dir}/
  {db}/                          nuccore or protein
    txid{taxid}_{region}_{segment}_{max}/
      sequences.gb               GenBank flat file (single-protein / whole-genome)
      manifest.json              download metadata
      _no_results                sentinel: NCBI returned 0 sequences (JSON timestamp)
      .lock                      present only during active download

    txid{taxid}_concat_marker{hash8}_{max}/
      <safe_marker>.gb           per-marker GenBank flat file (concat mode)
      manifest.json              includes marker_set_hash + marker filenames
      _no_results, .lock         same semantics as single-protein

The cache key encodes every parameter that affects what NCBI returns, so
changing region, segment, max_per_species, or — in concat mode — any element
of the marker_set spec automatically triggers a fresh download.  Negative
results (0 sequences found) are also cached so that repeated runs do not
re-query NCBI for species with no data.  Sentinels use the same TTL as
positive entries.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import time
from contextlib import contextmanager
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Generator

from .logger import get_logger

log = get_logger(__name__)

# Stale-lock threshold: a lock file older than this is assumed to belong to a
# crashed process and is silently overridden.
LOCK_TIMEOUT_SECONDS = 1800   # 30 minutes

# How long to sleep between lock-poll attempts, and the maximum number of polls
# before giving up waiting on a concurrent download.
LOCK_POLL_INTERVAL  = 5       # seconds
LOCK_POLL_MAX       = 120     # 10 minutes total wait


def marker_set_hash(marker_set: list[dict]) -> str:
    """Stable 8-character hex hash over a concat-mode marker_set spec.

    Hash inputs include every field that affects what NCBI returns (name,
    aliases, subfamily-aware aliases, length_range, locus_tag_hint), so any
    edit to the curated preset invalidates the cache for that family.
    """
    canonical: list[dict] = []
    for marker in marker_set:
        canonical.append({
            "name":     marker.get("name", ""),
            "aliases":  list(marker.get("aliases") or []),
            "subfam":   {k: list(v) for k, v in marker.items()
                         if k.startswith("aliases_") and isinstance(v, list)},
            "length":   list(marker.get("length_range") or []) or None,
            "locus":    marker.get("locus_tag_hint") or "",
        })
    payload = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:8]


class SequenceCache:
    """File-based cache for per-species GenBank downloads.

    Parameters
    ----------
    cache_dir:
        Root directory for the cache (e.g. ``~/.vfam_cache``).
        Created automatically if it does not exist.
    ttl_days:
        Maximum age of a valid cache entry in days.
        ``None`` means entries never expire.
    """

    def __init__(self, cache_dir: Path, ttl_days: int | None = None) -> None:
        self.cache_dir = Path(cache_dir).expanduser().resolve()
        self.ttl_days = ttl_days
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------

    def get(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> Path | None:
        """Return the path to a valid cached .gb file, or ``None`` on miss.

        A cache entry is invalid if:
        - the ``sequences.gb`` or ``manifest.json`` file is missing
        - the ``sequences.gb`` file is empty
        - the entry is older than ``ttl_days`` (when set)
        """
        entry_dir = self._entry_dir(taxid, db, region, segment, max_per_species)
        gb_file      = entry_dir / "sequences.gb"
        manifest_file = entry_dir / "manifest.json"

        if not gb_file.exists() or not manifest_file.exists():
            return None

        if gb_file.stat().st_size == 0:
            log.debug("Cache entry is empty — treating as miss: %s", gb_file)
            return None

        if self.ttl_days is not None:
            try:
                manifest = json.loads(manifest_file.read_text())
                downloaded = datetime.fromisoformat(manifest["downloaded"])
                age = datetime.now(timezone.utc) - downloaded
                if age > timedelta(days=self.ttl_days):
                    log.debug(
                        "Cache entry expired (age %.1f d > ttl %d d): %s",
                        age.total_seconds() / 86400,
                        self.ttl_days,
                        entry_dir,
                    )
                    return None
            except Exception as exc:
                log.debug("Cannot read cache manifest %s: %s — treating as miss", manifest_file, exc)
                return None

        return gb_file

    def get_empty(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> bool:
        """Return True if a valid negative-result sentinel exists for this key.

        A sentinel is invalid if it is older than ``ttl_days`` (when set).
        """
        sentinel = self._entry_dir(taxid, db, region, segment, max_per_species) / "_no_results"
        if not sentinel.exists():
            return False
        if self.ttl_days is not None:
            try:
                data = json.loads(sentinel.read_text())
                downloaded = datetime.fromisoformat(data["downloaded"])
                age = datetime.now(timezone.utc) - downloaded
                if age > timedelta(days=self.ttl_days):
                    log.debug(
                        "Negative cache entry expired (age %.1f d > ttl %d d): %s",
                        age.total_seconds() / 86400,
                        self.ttl_days,
                        sentinel,
                    )
                    return False
            except Exception as exc:
                log.debug("Cannot read negative sentinel %s: %s — treating as miss", sentinel, exc)
                return False
        return True

    def store_empty(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
        family: str = "",
    ) -> None:
        """Write a negative-result sentinel for this key."""
        entry_dir = self._entry_dir(taxid, db, region, segment, max_per_species)
        entry_dir.mkdir(parents=True, exist_ok=True)
        sentinel = entry_dir / "_no_results"
        payload = {
            "downloaded": datetime.now(timezone.utc).isoformat(),
            "taxid":      taxid,
            "family":     family,
        }
        try:
            sentinel.write_text(json.dumps(payload))
            log.debug("Stored negative result in cache: %s", sentinel)
        except Exception as exc:
            log.warning("Could not write negative cache sentinel (%s)", exc)

    def store(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
        source_gb: Path,
        n_records: int,
        query: str = "",
        family: str = "",
    ) -> Path:
        """Copy *source_gb* into the cache and write a manifest.

        Returns the path to the cached file.  Any error during the copy
        is caught and logged as a warning so the pipeline can continue.
        """
        entry_dir = self._entry_dir(taxid, db, region, segment, max_per_species)
        entry_dir.mkdir(parents=True, exist_ok=True)
        cached_gb = entry_dir / "sequences.gb"

        try:
            shutil.copy2(source_gb, cached_gb)
        except Exception as exc:
            log.warning("Could not write to sequence cache (%s) — skipping cache store", exc)
            return source_gb

        manifest = {
            "taxid":           taxid,
            "db":              db,
            "region":          region,
            "segment":         segment,
            "max_per_species": max_per_species,
            "downloaded":      datetime.now(timezone.utc).isoformat(),
            "n_records":       n_records,
            "query":           query,
            "family":          family,
        }
        try:
            (entry_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
        except Exception as exc:
            log.warning("Could not write cache manifest (%s)", exc)

        log.debug("Stored %d record(s) in cache: %s", n_records, cached_gb)
        return cached_gb

    @contextmanager
    def download_lock(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> Generator[bool, None, None]:
        """Context manager that acquires a per-entry advisory lock file.

        Yields ``True`` if the lock was acquired (this process should
        download), or ``False`` if the wait timed out (caller should
        re-check the cache and fall back to downloading without a lock).

        Stale locks older than ``LOCK_TIMEOUT_SECONDS`` are automatically
        removed so a crashed process never blocks the pipeline indefinitely.
        """
        entry_dir = self._entry_dir(taxid, db, region, segment, max_per_species)
        with self._acquire_lock(entry_dir, lock_label=f"taxid {taxid}") as acquired:
            yield acquired

    @contextmanager
    def _acquire_lock(
        self,
        entry_dir: Path,
        lock_label: str,
    ) -> Generator[bool, None, None]:
        """Shared advisory-lock implementation used by both single-protein and
        concat-mode lock contexts.  ``lock_label`` appears in info/warning logs
        so both flavours surface in `<family>.log` with the right context.
        """
        entry_dir.mkdir(parents=True, exist_ok=True)
        lock_path = entry_dir / ".lock"
        acquired  = False

        for poll in range(LOCK_POLL_MAX):
            try:
                # open("x") is atomic on POSIX: raises FileExistsError if present
                with lock_path.open("x") as fh:
                    fh.write(datetime.now(timezone.utc).isoformat())
                acquired = True
                break
            except FileExistsError:
                try:
                    age = time.time() - lock_path.stat().st_mtime
                    if age > LOCK_TIMEOUT_SECONDS:
                        log.warning(
                            "Removing stale download lock (%.0f s old): %s",
                            age, lock_path,
                        )
                        lock_path.unlink(missing_ok=True)
                        continue  # retry immediately
                except OSError:
                    pass

                if poll == 0:
                    log.info(
                        "Another job is downloading %s — waiting (max %d s) ...",
                        lock_label, LOCK_POLL_MAX * LOCK_POLL_INTERVAL,
                    )
                time.sleep(LOCK_POLL_INTERVAL)

        if not acquired:
            log.warning(
                "Timed out waiting for download lock on %s — proceeding anyway",
                lock_label,
            )

        try:
            yield acquired
        finally:
            if acquired:
                lock_path.unlink(missing_ok=True)

    def invalidate(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> bool:
        """Delete a single cache entry.  Returns ``True`` if it existed."""
        entry_dir = self._entry_dir(taxid, db, region, segment, max_per_species)
        if entry_dir.exists():
            shutil.rmtree(entry_dir)
            log.info("Invalidated cache entry: %s", entry_dir)
            return True
        return False

    def clear_family(self, family: str) -> int:
        """Delete all cache entries whose manifest or sentinel lists *family*.

        Returns the number of entries removed.
        """
        to_remove = []
        for manifest_path in self.cache_dir.rglob("manifest.json"):
            try:
                manifest = json.loads(manifest_path.read_text())
            except Exception:
                continue
            if manifest.get("family") == family:
                to_remove.append(manifest_path.parent)
        for sentinel_path in self.cache_dir.rglob("_no_results"):
            try:
                data = json.loads(sentinel_path.read_text())
            except Exception:
                continue
            if data.get("family") == family:
                entry_dir = sentinel_path.parent
                if entry_dir not in to_remove:
                    to_remove.append(entry_dir)
        for entry_dir in to_remove:
            shutil.rmtree(entry_dir)
            log.info("Cleared cache entry for %s: %s", family, entry_dir)
        return len(to_remove)

    def clear_all(self) -> int:
        """Delete every entry in the cache.  Returns the number removed."""
        removed = 0
        seen: set[Path] = set()
        for marker in list(self.cache_dir.rglob("manifest.json")) + \
                      list(self.cache_dir.rglob("_no_results")):
            entry_dir = marker.parent
            if entry_dir not in seen and entry_dir.exists():
                seen.add(entry_dir)
                shutil.rmtree(entry_dir)
                removed += 1
        return removed

    def stats(self) -> dict:
        """Return basic cache statistics (entry count, negative count, total size).

        Counts both single-protein ``sequences.gb`` and concat-mode per-marker
        ``<name>.gb`` files via one tree walk.
        """
        total         = 0
        n_empty       = 0
        size_bytes    = 0
        for gb in self.cache_dir.rglob("*.gb"):
            total      += 1
            size_bytes += gb.stat().st_size
        for _ in self.cache_dir.rglob("_no_results"):
            n_empty += 1
        return {
            "entries":         total,
            "empty_entries":   n_empty,
            "size_mb":         round(size_bytes / 1_048_576, 1),
            "cache_dir":       str(self.cache_dir),
        }

    # ------------------------------------------------------------------
    # Concat mode — multi-marker per-species cache.  Entry key is
    # (taxid, marker_set_hash, max_per_species); same TTL and lock
    # semantics as the single-protein methods above.
    # ------------------------------------------------------------------

    def get_concat(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
    ) -> Path | None:
        """Return the entry directory containing per-marker .gb files, or None.

        Same TTL semantics as ``get``.  Caller iterates the directory's
        ``*.gb`` files (each named after the safe-marker filename used by
        the fetcher) to load the cached records.
        """
        entry_dir = self._concat_entry_dir(taxid, marker_set_hash_str, max_per_species)
        manifest_file = entry_dir / "manifest.json"
        if not manifest_file.exists():
            return None
        try:
            manifest = json.loads(manifest_file.read_text())
        except Exception as exc:
            log.debug("Cannot read concat cache manifest %s: %s — treating as miss",
                      manifest_file, exc)
            return None
        if not manifest.get("marker_files"):
            return None
        if self.ttl_days is not None:
            try:
                downloaded = datetime.fromisoformat(manifest["downloaded"])
            except Exception:
                return None
            age = datetime.now(timezone.utc) - downloaded
            if age > timedelta(days=self.ttl_days):
                log.debug(
                    "Concat cache entry expired (age %.1f d > ttl %d d): %s",
                    age.total_seconds() / 86400, self.ttl_days, entry_dir,
                )
                return None
        return entry_dir

    def get_empty_concat(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
    ) -> bool:
        """Return True if a valid negative-result sentinel exists for this concat key."""
        sentinel = self._concat_entry_dir(
            taxid, marker_set_hash_str, max_per_species,
        ) / "_no_results"
        if not sentinel.exists():
            return False
        if self.ttl_days is not None:
            try:
                data = json.loads(sentinel.read_text())
                downloaded = datetime.fromisoformat(data["downloaded"])
                age = datetime.now(timezone.utc) - downloaded
                if age > timedelta(days=self.ttl_days):
                    return False
            except Exception:
                return False
        return True

    def store_concat(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
        marker_dir: Path,
        family: str = "",
    ) -> Path:
        """Copy every ``*.gb`` from *marker_dir* into the cache entry.

        Writes a manifest.json carrying the marker-set hash and the list of
        per-marker filenames.  Errors during copy are logged but do not
        propagate — the pipeline always has its working copy in *marker_dir*.
        """
        entry_dir = self._concat_entry_dir(taxid, marker_set_hash_str, max_per_species)
        entry_dir.mkdir(parents=True, exist_ok=True)
        marker_files: list[str] = []
        for gb in sorted(marker_dir.glob("*.gb")):
            try:
                shutil.copy2(gb, entry_dir / gb.name)
                marker_files.append(gb.name)
            except Exception as exc:
                log.warning("Could not copy %s into concat cache (%s)", gb, exc)
        manifest = {
            "taxid":              taxid,
            "marker_set_hash":    marker_set_hash_str,
            "max_per_species":    max_per_species,
            "marker_files":       marker_files,
            "downloaded":         datetime.now(timezone.utc).isoformat(),
            "family":             family,
        }
        try:
            (entry_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
        except Exception as exc:
            log.warning("Could not write concat cache manifest (%s)", exc)
        log.debug("Stored %d marker file(s) in concat cache: %s",
                  len(marker_files), entry_dir)
        return entry_dir

    def store_empty_concat(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
        family: str = "",
    ) -> None:
        """Write a negative-result sentinel for this concat key."""
        entry_dir = self._concat_entry_dir(taxid, marker_set_hash_str, max_per_species)
        entry_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "downloaded":      datetime.now(timezone.utc).isoformat(),
            "taxid":           taxid,
            "marker_set_hash": marker_set_hash_str,
            "family":          family,
        }
        try:
            (entry_dir / "_no_results").write_text(json.dumps(payload))
        except Exception as exc:
            log.warning("Could not write concat negative sentinel (%s)", exc)

    @contextmanager
    def download_lock_concat(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
    ) -> Generator[bool, None, None]:
        """Concat-mode counterpart to ``download_lock`` — same advisory-lock
        semantics, scoped to the concat cache entry."""
        entry_dir = self._concat_entry_dir(taxid, marker_set_hash_str, max_per_species)
        label = f"concat-mode taxid {taxid} (marker hash {marker_set_hash_str})"
        with self._acquire_lock(entry_dir, lock_label=label) as acquired:
            yield acquired

    def _concat_entry_dir(
        self,
        taxid: int,
        marker_set_hash_str: str,
        max_per_species: int,
    ) -> Path:
        # Layout: protein/txid<taxid>_concat_marker<hash8>_<max>/
        key = f"txid{taxid}_concat_marker{marker_set_hash_str}_{max_per_species}"
        return self.cache_dir / "protein" / key

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _cache_key(
        self,
        taxid: int,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> str:
        seg = segment.replace(" ", "_") if segment else "none"
        reg = region.replace(" ", "_").replace("/", "_")
        return f"txid{taxid}_{reg}_{seg}_{max_per_species}"

    def _entry_dir(
        self,
        taxid: int,
        db: str,
        region: str,
        segment: str | None,
        max_per_species: int,
    ) -> Path:
        return self.cache_dir / db / self._cache_key(taxid, region, segment, max_per_species)
