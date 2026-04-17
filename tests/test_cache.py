"""Tests for vfam_trees.cache."""
import json
import tempfile
from pathlib import Path

import pytest

from vfam_trees.cache import SequenceCache


@pytest.fixture
def cache_dir(tmp_path):
    return tmp_path / "cache"


@pytest.fixture
def sc(cache_dir):
    return SequenceCache(cache_dir)


def _make_gb(tmp_path, content="LOCUS test\n"):
    p = tmp_path / "seq.gb"
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# store / get round-trip
# ---------------------------------------------------------------------------

def test_store_and_get(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1234, "nuccore", "whole_genome", None, 200, gb, 5,
             query="test query", family="Flaviviridae")
    result = sc.get(1234, "nuccore", "whole_genome", None, 200)
    assert result is not None
    assert result.exists()


def test_manifest_contains_family(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1234, "nuccore", "B646L", None, 200, gb, 3, family="Asfarviridae")
    entry_dir = sc._entry_dir(1234, "nuccore", "B646L", None, 200)
    manifest = json.loads((entry_dir / "manifest.json").read_text())
    assert manifest["family"] == "Asfarviridae"


def test_get_returns_none_for_missing(sc):
    assert sc.get(9999, "nuccore", "whole_genome", None, 200) is None


def test_get_returns_none_for_empty_file(sc, tmp_path):
    gb = _make_gb(tmp_path, content="")
    sc.store(1234, "nuccore", "whole_genome", None, 200, gb, 0)
    assert sc.get(1234, "nuccore", "whole_genome", None, 200) is None


# ---------------------------------------------------------------------------
# clear_family
# ---------------------------------------------------------------------------

def test_clear_family_removes_matching_entries(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "B646L", None, 200, gb, 5, family="Asfarviridae")
    sc.store(1002, "nuccore", "B646L", None, 200, gb, 5, family="Asfarviridae")
    sc.store(2001, "nuccore", "whole_genome", None, 200, gb, 5, family="Flaviviridae")

    removed = sc.clear_family("Asfarviridae")
    assert removed == 2
    assert sc.get(1001, "nuccore", "B646L", None, 200) is None
    assert sc.get(1002, "nuccore", "B646L", None, 200) is None
    # Flaviviridae entry untouched
    assert sc.get(2001, "nuccore", "whole_genome", None, 200) is not None


def test_clear_family_returns_zero_when_none_found(sc):
    assert sc.clear_family("Nonexistentviridae") == 0


# ---------------------------------------------------------------------------
# clear_all
# ---------------------------------------------------------------------------

def test_clear_all(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5, family="FamilyA")
    sc.store(2001, "nuccore", "whole_genome", None, 200, gb, 5, family="FamilyB")
    removed = sc.clear_all()
    assert removed == 2
    assert sc.stats()["entries"] == 0


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

def test_stats_empty_cache(sc):
    st = sc.stats()
    assert st["entries"] == 0
    assert st["size_mb"] == 0.0


def test_stats_counts_entries(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5)
    sc.store(1002, "nuccore", "whole_genome", None, 200, gb, 5)
    assert sc.stats()["entries"] == 2


# ---------------------------------------------------------------------------
# TTL
# ---------------------------------------------------------------------------

def test_ttl_expiry(tmp_path):
    from datetime import datetime, timezone, timedelta
    sc = SequenceCache(tmp_path / "cache", ttl_days=1)
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5)

    # Manually backdate the manifest
    entry_dir = sc._entry_dir(1001, "nuccore", "whole_genome", None, 200)
    manifest_path = entry_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    old_time = (datetime.now(timezone.utc) - timedelta(days=2)).isoformat()
    manifest["downloaded"] = old_time
    manifest_path.write_text(json.dumps(manifest))

    assert sc.get(1001, "nuccore", "whole_genome", None, 200) is None


def test_ttl_not_expired(tmp_path):
    sc = SequenceCache(tmp_path / "cache", ttl_days=30)
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5)
    assert sc.get(1001, "nuccore", "whole_genome", None, 200) is not None
