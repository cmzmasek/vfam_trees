"""Tests for pipeline helper functions."""
import json
import logging
from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees.pipeline import (
    _mark_done, _mark_skipped,
    _inject_pasted_sequences, _load_fasta_file_entries,
    _filter_species_by_lineages,
)


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
# _inject_pasted_sequences — manual.include_seq / include_fasta_files injection
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

    def test_empty_organism_buckets_by_id_and_keeps_species_empty(self):
        # Single-token FASTA headers yield an empty organism.  meta["species"]
        # must stay empty (so the leaf-label formatter omits it instead of
        # duplicating the id), and the bucket key falls back to the id so
        # distinct entries are not merged into one "" bucket.
        species_data: dict = {}
        _inject_pasted_sequences(
            species_data, set(),
            [
                {"id": "SeqA", "organism": "", "sequence": "ACGT"},
                {"id": "SeqB", "organism": "", "sequence": "TGCA"},
            ],
            "Test",
        )
        assert set(species_data.keys()) == {"SeqA", "SeqB"}
        meta_a = species_data["SeqA"]["metadata"][0]
        assert meta_a["species"] == ""
        assert meta_a["seq_name"] == "SeqA"

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


# ---------------------------------------------------------------------------
# _load_fasta_file_entries — FASTA file parsing for include_fasta_files
# ---------------------------------------------------------------------------

class TestLoadFastaFileEntries:
    def _write_fasta(self, tmp_path, content: str, name: str = "seqs.fasta"):
        p = tmp_path / name
        p.write_text(content)
        return p

    def test_id_and_description_parsed(self, tmp_path):
        p = self._write_fasta(tmp_path, ">NC_001.1 Foo virus complete genome\nACGTACGT\n")
        entries = _load_fasta_file_entries(p, "Test")
        assert len(entries) == 1
        assert entries[0]["id"] == "NC_001.1"
        assert entries[0]["organism"] == "Foo virus complete genome"
        assert entries[0]["sequence"] == "ACGTACGT"

    def test_id_only_header_leaves_organism_empty(self, tmp_path):
        # A single-token header has no organism remainder.  The id must not
        # be reused as the organism — doing so duplicates it in the leaf
        # label (id and species would both render the same string).
        p = self._write_fasta(tmp_path, ">MySeq1\nACGT\n")
        entries = _load_fasta_file_entries(p, "Test")
        assert entries[0]["id"] == "MySeq1"
        assert entries[0]["organism"] == ""

    def test_sequence_uppercased(self, tmp_path):
        p = self._write_fasta(tmp_path, ">seq1 Foo\nacgtacgt\n")
        entries = _load_fasta_file_entries(p, "Test")
        assert entries[0]["sequence"] == "ACGTACGT"

    def test_multiple_records(self, tmp_path):
        p = self._write_fasta(tmp_path,
            ">seq1 Foo virus\nACGT\n>seq2 Bar virus\nTGCA\n")
        entries = _load_fasta_file_entries(p, "Test")
        assert len(entries) == 2
        assert entries[0]["id"] == "seq1"
        assert entries[1]["id"] == "seq2"

    def test_empty_file_returns_empty_list(self, tmp_path):
        p = self._write_fasta(tmp_path, "")
        entries = _load_fasta_file_entries(p, "Test")
        assert entries == []

    def test_missing_file_raises_value_error(self, tmp_path):
        p = tmp_path / "nonexistent.fasta"
        with pytest.raises(ValueError, match="could not read"):
            _load_fasta_file_entries(p, "Test")

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
# _filter_species_by_lineages
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


def _resolver_from(mapping: dict[str, set[int]]):
    """Build a fake taxid_resolver callback from a fixed mapping.

    Entries absent from *mapping* resolve to an empty set (simulating
    a name that didn't resolve in NCBI or a taxon with no descendants).
    """
    def _resolve(entries):
        return {e: mapping.get(e, set()) for e in entries}
    return _resolve


class TestFilterSpeciesByLineages:
    def _log(self):
        return _NullLog()

    def test_species_entry_resolves_to_single_taxid(self):
        species = [_make_species("Dengue virus", 11103)]
        resolver = _resolver_from({"Dengue virus": {11103}})
        result = _filter_species_by_lineages(
            species, ["Dengue virus"], "Flaviviridae", self._log(),
            taxid_resolver=resolver,
        )
        assert [sp["taxid"] for sp in result] == [11103]

    def test_genus_entry_pulls_in_all_species_under_genus(self):
        species = [
            _make_species("Dengue virus", 11103),
            _make_species("Zika virus", 64320),
            _make_species("Hantaan orthohantavirus", 11595),
        ]
        # "Flavivirus" genus subtree contains Dengue and Zika; not Hantaan.
        resolver = _resolver_from({"Flavivirus": {11103, 64320, 99999}})
        result = _filter_species_by_lineages(
            species, ["Flavivirus"], "Flaviviridae", self._log(),
            taxid_resolver=resolver,
        )
        assert {sp["taxid"] for sp in result} == {11103, 64320}

    def test_taxid_digit_string_entry_accepted(self):
        species = [_make_species("Dengue virus", 11103)]
        resolver = _resolver_from({"11103": {11103}})
        result = _filter_species_by_lineages(
            species, ["11103"], "Flaviviridae", self._log(),
            taxid_resolver=resolver,
        )
        assert len(result) == 1

    def test_unresolved_entry_produces_warning(self):
        species = [_make_species("Dengue virus", 11103)]
        log = self._log()
        resolver = _resolver_from({})  # "Typo virus" resolves to empty
        result = _filter_species_by_lineages(
            species, ["Typo virus"], "Flaviviridae", log,
            taxid_resolver=resolver,
        )
        assert result == []
        assert any("Typo virus" in w for w in log.warnings)

    def test_entry_outside_family_produces_warning(self):
        species = [_make_species("Dengue virus", 11103)]
        log = self._log()
        # "Mammalia" resolves to many species, none in Flaviviridae.
        resolver = _resolver_from({"Mammalia": {9606, 10090, 9913}})
        result = _filter_species_by_lineages(
            species, ["Mammalia"], "Flaviviridae", log,
            taxid_resolver=resolver,
        )
        assert result == []
        assert any("Mammalia" in w for w in log.warnings)

    def test_matched_and_unmatched_mixed(self):
        species = [_make_species("Dengue virus", 11103)]
        log = self._log()
        resolver = _resolver_from({
            "Dengue virus": {11103},
            "Mammalia": {9606, 10090},
        })
        result = _filter_species_by_lineages(
            species, ["Dengue virus", "Mammalia"], "Flaviviridae", log,
            taxid_resolver=resolver,
        )
        assert [sp["taxid"] for sp in result] == [11103]
        assert any("Mammalia" in w for w in log.warnings)

    def test_overlapping_lineages_deduplicate(self):
        species = [
            _make_species("Dengue virus", 11103),
            _make_species("Zika virus", 64320),
        ]
        resolver = _resolver_from({
            "Flavivirus": {11103, 64320},
            "Dengue virus": {11103},  # subset of Flavivirus
        })
        result = _filter_species_by_lineages(
            species, ["Flavivirus", "Dengue virus"], "Flaviviridae", self._log(),
            taxid_resolver=resolver,
        )
        assert {sp["taxid"] for sp in result} == {11103, 64320}

    def test_result_preserves_species_list_order(self):
        species = [
            _make_species("Alpha virus", 1),
            _make_species("Beta virus", 2),
            _make_species("Gamma virus", 3),
        ]
        resolver = _resolver_from({"Gamma virus": {3}, "Alpha virus": {1}})
        result = _filter_species_by_lineages(
            species, ["Gamma virus", "Alpha virus"], "F", self._log(),
            taxid_resolver=resolver,
        )
        # Discovered-list order, not user-list order.
        assert [sp["name"] for sp in result] == ["Alpha virus", "Gamma virus"]

    def test_empty_lineages_list_returns_empty(self):
        species = [_make_species("Dengue virus", 11103)]
        resolver = _resolver_from({})
        result = _filter_species_by_lineages(
            species, [], "Flaviviridae", self._log(),
            taxid_resolver=resolver,
        )
        assert result == []

    def test_empty_species_list_all_unmatched_warns(self):
        log = self._log()
        resolver = _resolver_from({"Dengue virus": {11103}})
        result = _filter_species_by_lineages(
            [], ["Dengue virus"], "Flaviviridae", log,
            taxid_resolver=resolver,
        )
        assert result == []
        assert log.warnings


# ---------------------------------------------------------------------------
# Config warning: nucleotide + named region + no segment
# ---------------------------------------------------------------------------

_MINIMAL_FAMILY_CFG = {
    "sequence": {"type": "nucleotide", "region": "glycoprotein g1", "segment": None},
    "download": {"max_per_species": 10},
    "quality": {"min_length": None, "max_ambiguous": 0.01, "exclude_organisms": []},
    "clustering": {
        "tool": "mmseqs2", "threshold_min": 0.7, "threshold_max": 0.99,
        "max_reps_500": 5, "max_reps_100": 2,
    },
    "targets": {"max_500": 50, "max_100": 10},
    "manual": {},
}


def test_nuc_named_region_no_segment_config_warning(tmp_path, caplog, monkeypatch):
    """Regression: run_family must emit a WARNING when seq_type=nucleotide,
    region is not whole_genome, and segment is null.  This combination silently
    returned zero sequences before v1.2.35 because the nuccore query searched
    only [Gene], which is sparsely populated."""
    import vfam_trees.pipeline as pipeline_mod

    monkeypatch.setattr(pipeline_mod, "load_global_config",
                        lambda p: {"ncbi": {"email": "t@t.com", "api_key": None}})
    monkeypatch.setattr(pipeline_mod, "configure_entrez", lambda email, api_key: None)
    monkeypatch.setattr(pipeline_mod, "load_family_annotations", lambda p: {})
    monkeypatch.setattr(pipeline_mod, "get_family_taxid", lambda f: None)
    monkeypatch.setattr(pipeline_mod, "load_family_config",
                        lambda f, d, g: (_MINIMAL_FAMILY_CFG, False))
    monkeypatch.setattr(pipeline_mod, "_save_config_copy", lambda cfg, p: None)
    monkeypatch.setattr(pipeline_mod, "discover_species", lambda f: [])

    with caplog.at_level(logging.WARNING):
        pipeline_mod.run_family(
            family="TestVirus",
            global_config_path=tmp_path / "global.yaml",
            configs_dir=tmp_path / "configs",
            output_dir=tmp_path / "results",
        )

    warning_texts = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("glycoprotein g1" in m for m in warning_texts), (
        f"Expected config warning about 'glycoprotein g1', got: {warning_texts}"
    )


def test_whole_genome_nucleotide_no_warning(tmp_path, caplog, monkeypatch):
    """whole_genome nucleotide runs must not trigger the niche-config warning."""
    import vfam_trees.pipeline as pipeline_mod

    cfg = {**_MINIMAL_FAMILY_CFG, "sequence": {"type": "nucleotide", "region": "whole_genome", "segment": None}}
    monkeypatch.setattr(pipeline_mod, "load_global_config",
                        lambda p: {"ncbi": {"email": "t@t.com", "api_key": None}})
    monkeypatch.setattr(pipeline_mod, "configure_entrez", lambda email, api_key: None)
    monkeypatch.setattr(pipeline_mod, "load_family_annotations", lambda p: {})
    monkeypatch.setattr(pipeline_mod, "get_family_taxid", lambda f: None)
    monkeypatch.setattr(pipeline_mod, "load_family_config", lambda f, d, g: (cfg, False))
    monkeypatch.setattr(pipeline_mod, "_save_config_copy", lambda cfg, p: None)
    monkeypatch.setattr(pipeline_mod, "discover_species", lambda f: [])

    with caplog.at_level(logging.WARNING):
        pipeline_mod.run_family(
            family="TestVirus",
            global_config_path=tmp_path / "global.yaml",
            configs_dir=tmp_path / "configs",
            output_dir=tmp_path / "results",
        )

    warning_texts = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
    assert not any("niche configuration" in m for m in warning_texts)


def test_protein_named_region_no_warning(tmp_path, caplog, monkeypatch):
    """protein + named region is the normal path and must not trigger the warning."""
    import vfam_trees.pipeline as pipeline_mod

    cfg = {**_MINIMAL_FAMILY_CFG, "sequence": {"type": "protein", "region": "glycoprotein g1", "segment": None}}
    monkeypatch.setattr(pipeline_mod, "load_global_config",
                        lambda p: {"ncbi": {"email": "t@t.com", "api_key": None}})
    monkeypatch.setattr(pipeline_mod, "configure_entrez", lambda email, api_key: None)
    monkeypatch.setattr(pipeline_mod, "load_family_annotations", lambda p: {})
    monkeypatch.setattr(pipeline_mod, "get_family_taxid", lambda f: None)
    monkeypatch.setattr(pipeline_mod, "load_family_config", lambda f, d, g: (cfg, False))
    monkeypatch.setattr(pipeline_mod, "_save_config_copy", lambda cfg, p: None)
    monkeypatch.setattr(pipeline_mod, "discover_species", lambda f: [])

    with caplog.at_level(logging.WARNING):
        pipeline_mod.run_family(
            family="TestVirus",
            global_config_path=tmp_path / "global.yaml",
            configs_dir=tmp_path / "configs",
            output_dir=tmp_path / "results",
        )

    warning_texts = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
    assert not any("niche configuration" in m for m in warning_texts)
