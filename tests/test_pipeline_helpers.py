"""Tests for pipeline helper functions."""
import json
from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees.pipeline import _mark_done, _mark_skipped, _inject_pasted_sequences, _filter_species_by_include_list


@pytest.fixture
def dirs(tmp_path):
    family_dir = tmp_path / "Flaviviridae_11050"
    family_dir.mkdir()
    output_dir = tmp_path
    return family_dir, output_dir


def test_mark_done_writes_status(dirs):
    family_dir, output_dir = dirs
    _mark_done(family_dir, "Flaviviridae", output_dir)
    status = json.loads((family_dir / ".status.json").read_text())
    assert status["status"] == "done"
    assert status["family"] == "Flaviviridae"


def test_mark_done_writes_sentinel(dirs):
    family_dir, output_dir = dirs
    _mark_done(family_dir, "Flaviviridae", output_dir)
    assert (output_dir / ".done_Flaviviridae").exists()


def test_mark_skipped_writes_status(dirs):
    family_dir, output_dir = dirs
    _mark_skipped(family_dir, "Flaviviridae", output_dir, "no species found")
    status = json.loads((family_dir / ".status.json").read_text())
    assert status["status"] == "skipped"
    assert "no species" in status["reason"]


def test_mark_skipped_writes_sentinel(dirs):
    family_dir, output_dir = dirs
    _mark_skipped(family_dir, "Flaviviridae", output_dir, "too few sequences")
    assert (output_dir / ".done_Flaviviridae").exists()


# ---------------------------------------------------------------------------
# _inject_pasted_sequences — manual.include_fasta injection
# ---------------------------------------------------------------------------

def _make_fetched_species(name: str, accession: str, seq: str = "ACGT") -> dict:
    """Build a species_data bucket as the fetch loop would."""
    rec = SeqRecord(Seq(seq), id=accession, description=f"{name} reference")
    meta = {
        "accession": accession,
        "seq_name": f"{name} reference",
        "species": name,
        "strain": "unknown",
        "host": "unknown",
        "collection_date": "unknown",
        "location": "unknown",
        "taxon_id": "12345",
        "length": len(seq),
        "lineage": ["Viruses"],
    }
    return {"records": [rec], "metadata": [meta]}


class TestInjectPastedSequences:
    def test_empty_entries_is_noop(self):
        species_data = {"Foo virus": _make_fetched_species("Foo virus", "NC_001.1")}
        manual_include = set()
        n = _inject_pasted_sequences(species_data, manual_include, [], "Test")
        assert n == 0
        assert list(species_data.keys()) == ["Foo virus"]
        assert manual_include == set()

    def test_single_entry_creates_new_species_bucket(self):
        species_data: dict = {}
        manual_include: set = set()
        n = _inject_pasted_sequences(
            species_data, manual_include,
            [{"id": "MySeq1", "organism": "Novel virus", "sequence": "ACGTACGT"}],
            "Test",
        )
        assert n == 1
        assert "Novel virus" in species_data
        bucket = species_data["Novel virus"]
        assert len(bucket["records"]) == 1
        rec = bucket["records"][0]
        assert rec.id == "MySeq1"
        assert rec.description == "Novel virus"
        assert str(rec.seq) == "ACGTACGT"
        assert "MySeq1" in manual_include

    def test_metadata_fields_populated(self):
        species_data: dict = {}
        _inject_pasted_sequences(
            species_data, set(),
            [{"id": "MySeq1", "organism": "Novel virus", "sequence": "ACGTACGT"}],
            "Test",
        )
        meta = species_data["Novel virus"]["metadata"][0]
        assert meta["accession"] == "MySeq1"
        assert meta["species"] == "Novel virus"
        assert meta["seq_name"] == "Novel virus"
        assert meta["taxon_id"] == ""
        assert meta["lineage"] == []
        assert meta["length"] == 8
        assert meta["strain"] == "unknown"
        assert meta["host"] == "unknown"

    def test_entry_joins_existing_species_bucket(self):
        species_data = {"Foo virus": _make_fetched_species("Foo virus", "NC_001.1")}
        _inject_pasted_sequences(
            species_data, set(),
            [{"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGT"}],
            "Test",
        )
        bucket = species_data["Foo virus"]
        assert len(bucket["records"]) == 2
        assert {r.id for r in bucket["records"]} == {"NC_001.1", "MySeq1"}
        assert list(species_data.keys()) == ["Foo virus"]

    def test_multiple_entries_same_organism_share_bucket(self):
        species_data: dict = {}
        _inject_pasted_sequences(
            species_data, set(),
            [
                {"id": "MySeq1", "organism": "Novel virus", "sequence": "ACGT"},
                {"id": "MySeq2", "organism": "Novel virus", "sequence": "TGCA"},
            ],
            "Test",
        )
        assert list(species_data.keys()) == ["Novel virus"]
        assert len(species_data["Novel virus"]["records"]) == 2

    def test_multiple_entries_different_organisms_separate_buckets(self):
        species_data: dict = {}
        manual_include: set = set()
        n = _inject_pasted_sequences(
            species_data, manual_include,
            [
                {"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGT"},
                {"id": "MySeq2", "organism": "Bar virus", "sequence": "TGCA"},
            ],
            "Test",
        )
        assert n == 2
        assert set(species_data.keys()) == {"Foo virus", "Bar virus"}
        assert manual_include == {"MySeq1", "MySeq2"}

    def test_id_collision_with_fetched_accession_raises(self):
        species_data = {"Foo virus": _make_fetched_species("Foo virus", "NC_001.1")}
        with pytest.raises(ValueError, match=r"collides with an accession returned by NCBI"):
            _inject_pasted_sequences(
                species_data, set(),
                [{"id": "NC_001.1", "organism": "Foo virus", "sequence": "ACGT"}],
                "Flaviviridae",
            )

    def test_id_collision_error_names_family_and_id(self):
        species_data = {"Foo virus": _make_fetched_species("Foo virus", "NC_001.1")}
        with pytest.raises(ValueError) as exc:
            _inject_pasted_sequences(
                species_data, set(),
                [{"id": "NC_001.1", "organism": "Foo virus", "sequence": "ACGT"}],
                "Flaviviridae",
            )
        msg = str(exc.value)
        assert "Flaviviridae" in msg
        assert "NC_001.1" in msg

    def test_collision_check_runs_before_injection(self):
        """If any entry collides, no entries should be injected (fail-fast)."""
        species_data = {"Foo virus": _make_fetched_species("Foo virus", "NC_001.1")}
        manual_include: set = set()
        with pytest.raises(ValueError):
            _inject_pasted_sequences(
                species_data, manual_include,
                [
                    {"id": "NC_001.1", "organism": "Foo virus", "sequence": "ACGT"},
                    {"id": "MySeq2",   "organism": "Bar virus", "sequence": "TGCA"},
                ],
                "Test",
            )
        # The first entry hit the collision before MySeq2 was processed.
        assert "MySeq2" not in manual_include
        assert "Bar virus" not in species_data

    def test_manual_include_ids_accumulates_across_calls(self):
        species_data: dict = {}
        manual_include = {"pre-existing-id"}
        _inject_pasted_sequences(
            species_data, manual_include,
            [{"id": "MySeq1", "organism": "Foo virus", "sequence": "ACGT"}],
            "Test",
        )
        assert manual_include == {"pre-existing-id", "MySeq1"}


# ---------------------------------------------------------------------------
# _filter_species_by_include_list
# ---------------------------------------------------------------------------

def _make_species(name: str, taxid: int) -> dict:
    return {"name": name, "taxid": taxid}


class _NullLog:
    """Minimal logger stub that records warning calls."""
    def __init__(self):
        self.warnings: list[str] = []

    def info(self, *a, **kw):
        pass

    def warning(self, msg, *args, **kw):
        self.warnings.append(msg % args if args else msg)


class TestFilterSpeciesByIncludeList:
    def _log(self):
        return _NullLog()

    def test_name_match_case_insensitive(self):
        species = [_make_species("Dengue virus", 11103)]
        log = self._log()
        result = _filter_species_by_include_list(species, ["dengue virus"], "Flaviviridae", log)
        assert len(result) == 1
        assert result[0]["taxid"] == 11103

    def test_name_match_exact_case(self):
        species = [_make_species("Dengue virus", 11103)]
        result = _filter_species_by_include_list(species, ["Dengue virus"], "Flaviviridae", self._log())
        assert len(result) == 1

    def test_taxid_match_by_digit_string(self):
        species = [_make_species("Dengue virus", 11103)]
        result = _filter_species_by_include_list(species, ["11103"], "Flaviviridae", self._log())
        assert len(result) == 1
        assert result[0]["name"] == "Dengue virus"

    def test_unmatched_entry_produces_warning(self):
        species = [_make_species("Dengue virus", 11103)]
        log = self._log()
        result = _filter_species_by_include_list(species, ["Zika virus"], "Flaviviridae", log)
        assert result == []
        assert any("Zika virus" in w for w in log.warnings)

    def test_matched_and_unmatched_mixed(self):
        species = [
            _make_species("Dengue virus", 11103),
            _make_species("Zika virus", 64320),
        ]
        log = self._log()
        result = _filter_species_by_include_list(species, ["Dengue virus", "West Nile virus"], "F", log)
        assert len(result) == 1
        assert result[0]["taxid"] == 11103
        assert any("West Nile virus" in w for w in log.warnings)

    def test_name_and_taxid_both_accepted(self):
        species = [
            _make_species("Dengue virus", 11103),
            _make_species("Zika virus", 64320),
        ]
        result = _filter_species_by_include_list(
            species, ["Dengue virus", "64320"], "Flaviviridae", self._log()
        )
        assert len(result) == 2
        assert {sp["taxid"] for sp in result} == {11103, 64320}

    def test_duplicate_taxid_via_name_and_taxid_deduped(self):
        species = [_make_species("Dengue virus", 11103)]
        result = _filter_species_by_include_list(
            species, ["Dengue virus", "11103"], "Flaviviridae", self._log()
        )
        assert len(result) == 1

    def test_order_follows_include_list_not_species_list(self):
        species = [
            _make_species("Alpha virus", 1),
            _make_species("Beta virus", 2),
            _make_species("Gamma virus", 3),
        ]
        result = _filter_species_by_include_list(
            species, ["Gamma virus", "Alpha virus"], "F", self._log()
        )
        assert [sp["name"] for sp in result] == ["Gamma virus", "Alpha virus"]

    def test_empty_include_list_returns_empty(self):
        species = [_make_species("Dengue virus", 11103)]
        result = _filter_species_by_include_list(species, [], "Flaviviridae", self._log())
        assert result == []

    def test_empty_species_list_all_unmatched_warns(self):
        log = self._log()
        result = _filter_species_by_include_list([], ["Dengue virus"], "Flaviviridae", log)
        assert result == []
        assert log.warnings
