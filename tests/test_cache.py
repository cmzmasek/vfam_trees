"""Tests for vfam_trees.cache."""
import copy
import json
import tempfile
from pathlib import Path

import pytest

from vfam_trees.cache import SequenceCache, marker_set_hash


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


# ---------------------------------------------------------------------------
# store_empty / get_empty
# ---------------------------------------------------------------------------

def test_store_empty_and_get_empty(sc):
    sc.store_empty(5001, "nuccore", "whole_genome", None, 200, family="Testviridae")
    assert sc.get_empty(5001, "nuccore", "whole_genome", None, 200) is True


def test_get_empty_returns_false_when_not_stored(sc):
    assert sc.get_empty(9999, "nuccore", "whole_genome", None, 200) is False


def test_store_empty_sentinel_contains_family(sc):
    sc.store_empty(5001, "nuccore", "whole_genome", None, 200, family="Flaviviridae")
    sentinel = sc._entry_dir(5001, "nuccore", "whole_genome", None, 200) / "_no_results"
    data = json.loads(sentinel.read_text())
    assert data["family"] == "Flaviviridae"
    assert "downloaded" in data


def test_empty_ttl_expiry(tmp_path):
    from datetime import datetime, timezone, timedelta
    sc = SequenceCache(tmp_path / "cache", ttl_days=1)
    sc.store_empty(5001, "nuccore", "whole_genome", None, 200)
    sentinel = sc._entry_dir(5001, "nuccore", "whole_genome", None, 200) / "_no_results"
    data = json.loads(sentinel.read_text())
    old_time = (datetime.now(timezone.utc) - timedelta(days=2)).isoformat()
    data["downloaded"] = old_time
    sentinel.write_text(json.dumps(data))
    assert sc.get_empty(5001, "nuccore", "whole_genome", None, 200) is False


def test_empty_ttl_not_expired(tmp_path):
    sc = SequenceCache(tmp_path / "cache", ttl_days=30)
    sc.store_empty(5001, "nuccore", "whole_genome", None, 200)
    assert sc.get_empty(5001, "nuccore", "whole_genome", None, 200) is True


def test_clear_family_removes_empty_entries(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5, family="Flaviviridae")
    sc.store_empty(1002, "nuccore", "whole_genome", None, 200, family="Flaviviridae")
    sc.store_empty(2001, "nuccore", "whole_genome", None, 200, family="Herpesviridae")
    removed = sc.clear_family("Flaviviridae")
    assert removed == 2
    assert sc.get_empty(1002, "nuccore", "whole_genome", None, 200) is False
    assert sc.get_empty(2001, "nuccore", "whole_genome", None, 200) is True


def test_clear_all_removes_empty_entries(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5)
    sc.store_empty(2001, "nuccore", "whole_genome", None, 200)
    removed = sc.clear_all()
    assert removed == 2
    assert sc.stats()["entries"] == 0
    assert sc.stats()["empty_entries"] == 0


def test_stats_counts_empty_entries(sc, tmp_path):
    gb = _make_gb(tmp_path)
    sc.store(1001, "nuccore", "whole_genome", None, 200, gb, 5)
    sc.store_empty(2001, "nuccore", "whole_genome", None, 200)
    sc.store_empty(3001, "nuccore", "whole_genome", None, 200)
    st = sc.stats()
    assert st["entries"] == 1
    assert st["empty_entries"] == 2


# ---------------------------------------------------------------------------
# marker_set_hash — stable concat-mode cache key extension
# ---------------------------------------------------------------------------

class TestMarkerSetHash:
    def _spec(self):
        return [
            {"name": "DNA polymerase",
             "aliases": ["DNA-directed DNA polymerase"],
             "locus_tag_hint": r"polB"},
            {"name": "MCP",
             "aliases": [],
             "length_range": [400, 800]},
        ]

    def test_deterministic(self):
        h1 = marker_set_hash(self._spec())
        h2 = marker_set_hash(self._spec())
        assert h1 == h2
        assert len(h1) == 8 and all(c in "0123456789abcdef" for c in h1)

    def test_alias_change_invalidates(self):
        s1 = self._spec()
        s2 = copy.deepcopy(s1)
        s2[0]["aliases"].append("DNA pol")
        assert marker_set_hash(s1) != marker_set_hash(s2)

    def test_locus_tag_hint_change_invalidates(self):
        s1 = self._spec()
        s2 = copy.deepcopy(s1)
        s2[0]["locus_tag_hint"] = r"polB|G1211R"
        assert marker_set_hash(s1) != marker_set_hash(s2)

    def test_length_range_change_invalidates(self):
        s1 = self._spec()
        s2 = copy.deepcopy(s1)
        s2[1]["length_range"] = [400, 900]
        assert marker_set_hash(s1) != marker_set_hash(s2)

    def test_subfamily_alias_change_invalidates(self):
        s1 = self._spec()
        s2 = copy.deepcopy(s1)
        s2[0]["aliases_Entomopoxvirinae"] = ["DNA polymerase B"]
        assert marker_set_hash(s1) != marker_set_hash(s2)

    def test_marker_order_matters(self):
        # Stability for canonical input order — reordering marker entries
        # corresponds to a different concat layout, so the hash should change.
        s1 = self._spec()
        s2 = list(reversed(self._spec()))
        assert marker_set_hash(s1) != marker_set_hash(s2)

    def test_unrelated_field_does_not_affect_hash(self):
        # Hash inputs are limited to fields that affect what NCBI returns.
        # Adding an irrelevant key (e.g. a comment) should not invalidate.
        s1 = self._spec()
        s2 = copy.deepcopy(s1)
        s2[0]["__comment"] = "added by curator"
        assert marker_set_hash(s1) == marker_set_hash(s2)


# ---------------------------------------------------------------------------
# Concat cache round-trip
# ---------------------------------------------------------------------------

def _make_marker_dir(tmp_path, marker_files: dict[str, str]) -> Path:
    """Build a sp_work-style directory with one .gb file per marker."""
    d = tmp_path / "sp_work_concat"
    d.mkdir(parents=True, exist_ok=True)
    for fname, content in marker_files.items():
        (d / fname).write_text(content)
    return d


class TestConcatCache:
    def test_store_then_get_round_trip(self, sc, tmp_path):
        marker_dir = _make_marker_dir(tmp_path, {
            "DNA_polymerase.gb": "LOCUS YP_1\nproduct: DNA polymerase\n//\n",
            "MCP.gb":            "LOCUS YP_2\nproduct: MCP\n//\n",
        })
        sc.store_concat(taxid=12345, marker_set_hash_str="abc12345",
                        max_per_species=300, marker_dir=marker_dir,
                        family="Poxviridae")
        cached = sc.get_concat(12345, "abc12345", 300)
        assert cached is not None
        assert (cached / "DNA_polymerase.gb").exists()
        assert (cached / "MCP.gb").exists()
        assert (cached / "manifest.json").exists()
        manifest = json.loads((cached / "manifest.json").read_text())
        assert manifest["marker_set_hash"] == "abc12345"
        assert manifest["family"] == "Poxviridae"
        assert set(manifest["marker_files"]) == {"DNA_polymerase.gb", "MCP.gb"}

    def test_get_returns_none_for_missing(self, sc):
        assert sc.get_concat(99999, "deadbeef", 300) is None

    def test_get_returns_none_when_no_gb_files(self, sc, tmp_path):
        # Manifest exists but no .gb files (corrupted state)
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        sc.store_concat(taxid=1, marker_set_hash_str="abc12345",
                        max_per_species=300, marker_dir=empty_dir,
                        family="X")
        # Now no .gb files were copied; get_concat must return None
        assert sc.get_concat(1, "abc12345", 300) is None

    def test_store_empty_then_get_empty(self, sc):
        sc.store_empty_concat(7777, "abc12345", 300, family="Asfarviridae")
        assert sc.get_empty_concat(7777, "abc12345", 300) is True
        assert sc.get_empty_concat(7777, "abc12345", 200) is False  # different max
        assert sc.get_empty_concat(7777, "different", 300) is False  # different hash

    def test_marker_hash_isolates_caches(self, sc, tmp_path):
        # Two runs of the same taxid with different marker sets must NOT collide.
        d1 = _make_marker_dir(tmp_path / "v1", {"polB.gb": "LOCUS A\n//\n"})
        sc.store_concat(taxid=42, marker_set_hash_str="hash1111",
                        max_per_species=300, marker_dir=d1, family="X")
        # Different marker hash → cache miss
        assert sc.get_concat(42, "hash2222", 300) is None
        # Original hash still hits
        assert sc.get_concat(42, "hash1111", 300) is not None

    def test_concat_cache_does_not_collide_with_single_protein(self, sc, tmp_path):
        # Single-protein entry for taxid 42, then concat entry for the same
        # taxid — both must coexist independently.
        gb = _make_gb(tmp_path)
        sc.store(42, "protein", "DNA polymerase", None, 300, gb, 5,
                 family="X")
        d = _make_marker_dir(tmp_path, {"polB.gb": "LOCUS A\n//\n"})
        sc.store_concat(taxid=42, marker_set_hash_str="abc12345",
                        max_per_species=300, marker_dir=d, family="X")
        # Both retrievable
        assert sc.get(42, "protein", "DNA polymerase", None, 300) is not None
        assert sc.get_concat(42, "abc12345", 300) is not None

    def test_concat_cache_in_protein_subdirectory(self, sc, tmp_path):
        # Layout invariant: concat entries live under cache_dir/protein/
        d = _make_marker_dir(tmp_path, {"polB.gb": "LOCUS A\n//\n"})
        sc.store_concat(taxid=42, marker_set_hash_str="abc12345",
                        max_per_species=300, marker_dir=d, family="X")
        cached = sc.get_concat(42, "abc12345", 300)
        assert cached is not None
        # Path includes /protein/ segment
        assert "protein" in cached.parts

    def test_stats_counts_concat_marker_files(self, sc, tmp_path):
        d = _make_marker_dir(tmp_path, {
            "polB.gb": "LOCUS A\n//\n",
            "MCP.gb":  "LOCUS B\n//\n",
            "hel.gb":  "LOCUS C\n//\n",
        })
        sc.store_concat(taxid=42, marker_set_hash_str="abc12345",
                        max_per_species=300, marker_dir=d, family="X")
        st = sc.stats()
        # 3 per-marker .gb files counted as 3 entries
        assert st["entries"] >= 3
