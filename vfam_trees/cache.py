"""Global file-based cache for per-species GenBank sequence downloads.

Cache layout
------------
{cache_dir}/
  {db}/                          nuccore or protein
    txid{taxid}_{region}_{segment}_{max}/
      sequences.gb               GenBank flat file
      manifest.json              download metadata
      .lock                      present only during active download

The cache key encodes every parameter that affects what NCBI returns, so
changing region, segment, or max_per_species automatically triggers a fresh
download.
"""
from __future__ import annotations

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
                        "Another job is downloading taxid %d — waiting (max %d s) ...",
                        taxid, LOCK_POLL_MAX * LOCK_POLL_INTERVAL,
                    )
                time.sleep(LOCK_POLL_INTERVAL)

        if not acquired:
            log.warning(
                "Timed out waiting for download lock on taxid %d — proceeding anyway",
                taxid,
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
        """Delete all cache entries whose manifest lists *family*.

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
        for entry_dir in to_remove:
            shutil.rmtree(entry_dir)
            log.info("Cleared cache entry for %s: %s", family, entry_dir)
        return len(to_remove)

    def clear_all(self) -> int:
        """Delete every entry in the cache.  Returns the number removed."""
        removed = 0
        for manifest_path in list(self.cache_dir.rglob("manifest.json")):
            entry_dir = manifest_path.parent
            if entry_dir.exists():
                shutil.rmtree(entry_dir)
                removed += 1
        return removed

    def stats(self) -> dict:
        """Return basic cache statistics (entry count and total size)."""
        total      = 0
        size_bytes = 0
        for gb in self.cache_dir.rglob("sequences.gb"):
            total      += 1
            size_bytes += gb.stat().st_size
        return {
            "entries":  total,
            "size_mb":  round(size_bytes / 1_048_576, 1),
            "cache_dir": str(self.cache_dir),
        }

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
