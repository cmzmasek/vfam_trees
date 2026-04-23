"""Tests for vfam_trees.rename — short ID assignment and name restoration."""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees.rename import (
    _build_display_name,
    _family_prefix,
    assign_short_ids,
    load_id_map,
    restore_fasta_names,
    restore_newick_names,
)


def _rec(acc, seq="ATGC"):
    return SeqRecord(Seq(seq), id=acc, description="raw desc")


def _meta(species="Dengue virus", strain="DEN1", accession="NC_001477", host="Homo sapiens"):
    return {"species": species, "strain": strain, "accession": accession, "host": host}


# ---------------------------------------------------------------------------
# _family_prefix
# ---------------------------------------------------------------------------

class TestFamilyPrefix:
    def test_four_consonants(self):
        # "Flaviviridae" consonants: F, L, V, V, R, D → "FLVV"
        assert _family_prefix("Flaviviridae") == "FLVV"

    def test_first_four_consonants_only(self):
        # "Coronaviridae" consonants: C, R, N, V, R, D → "CRNV"
        assert _family_prefix("Coronaviridae") == "CRNV"

    def test_too_few_consonants_uses_first_four_chars(self):
        # Hypothetical family with <4 consonants: falls back to first 4 uppercase
        assert _family_prefix("Abc") == "ABC"

    def test_vowel_heavy_name_fallback(self):
        # "Aeiou" has 0 consonants → fallback to first 4 chars uppercased
        assert _family_prefix("Aeiou") == "AEIO"


# ---------------------------------------------------------------------------
# _build_display_name
# ---------------------------------------------------------------------------

class TestBuildDisplayName:
    def test_full_metadata(self):
        m = _meta()
        assert _build_display_name(m) == "Dengue_virus|DEN1|NC_001477|Homo_sapiens"

    def test_host_unknown_omitted(self):
        m = _meta(host="unknown")
        assert _build_display_name(m) == "Dengue_virus|DEN1|NC_001477"

    def test_missing_fields_replaced_by_unknown(self):
        m = {"accession": "ACC1"}
        assert _build_display_name(m) == "unknown|unknown|ACC1"

    def test_whitespace_replaced_with_underscore(self):
        m = _meta(species="Zika virus", host="Aedes aegypti")
        name = _build_display_name(m)
        assert "Zika_virus" in name
        assert "Aedes_aegypti" in name
        assert " " not in name


# ---------------------------------------------------------------------------
# assign_short_ids
# ---------------------------------------------------------------------------

class TestAssignShortIds:
    def test_short_ids_sequential_with_family_prefix(self, tmp_path):
        records = [_rec("ACC1"), _rec("ACC2"), _rec("ACC3")]
        metas = [_meta(accession="ACC1"), _meta(accession="ACC2"), _meta(accession="ACC3")]
        id_map_path = tmp_path / "id_map.tsv"
        renamed, short_to_display = assign_short_ids(records, metas, "Flaviviridae", id_map_path)
        ids = [r.id for r in renamed]
        assert ids == ["FLVV_000001", "FLVV_000002", "FLVV_000003"]
        assert len(short_to_display) == 3

    def test_short_ids_zero_padded_to_6_digits(self, tmp_path):
        records = [_rec(f"ACC{i}") for i in range(5)]
        metas = [_meta(accession=f"ACC{i}") for i in range(5)]
        id_map_path = tmp_path / "id_map.tsv"
        renamed, _ = assign_short_ids(records, metas, "Flaviviridae", id_map_path)
        for i, r in enumerate(renamed, start=1):
            assert r.id == f"FLVV_{i:06d}"

    def test_id_map_tsv_written_with_header(self, tmp_path):
        records = [_rec("ACC1")]
        metas = [_meta(accession="ACC1")]
        id_map_path = tmp_path / "id_map.tsv"
        assign_short_ids(records, metas, "Flaviviridae", id_map_path)
        rows = id_map_path.read_text().splitlines()
        assert rows[0] == "short_id\taccession\tdisplay_name"
        assert rows[1].startswith("FLVV_000001\tACC1\t")

    def test_original_sequence_preserved_on_renamed_record(self, tmp_path):
        rec = _rec("ACC1", seq="ATGCATGC")
        records = [rec]
        metas = [_meta(accession="ACC1")]
        id_map_path = tmp_path / "id_map.tsv"
        renamed, _ = assign_short_ids(records, metas, "Flaviviridae", id_map_path)
        assert str(renamed[0].seq) == "ATGCATGC"

    def test_description_cleared_on_renamed_record(self, tmp_path):
        rec = _rec("ACC1")
        rec.description = "some long description"
        records = [rec]
        metas = [_meta(accession="ACC1")]
        id_map_path = tmp_path / "id_map.tsv"
        renamed, _ = assign_short_ids(records, metas, "Flaviviridae", id_map_path)
        assert renamed[0].description == ""


# ---------------------------------------------------------------------------
# load_id_map (round-trip)
# ---------------------------------------------------------------------------

def test_load_id_map_roundtrip(tmp_path):
    records = [_rec("ACC1"), _rec("ACC2")]
    metas = [_meta(accession="ACC1"), _meta(accession="ACC2")]
    id_map_path = tmp_path / "id_map.tsv"
    _, short_to_display = assign_short_ids(records, metas, "Flaviviridae", id_map_path)
    loaded = load_id_map(id_map_path)
    assert loaded == short_to_display


# ---------------------------------------------------------------------------
# restore_newick_names
# ---------------------------------------------------------------------------

def test_restore_newick_names_replaces_short_ids(tmp_path):
    nwk_in = tmp_path / "in.nwk"
    nwk_out = tmp_path / "out.nwk"
    nwk_in.write_text("((FLVV_000001:0.1,FLVV_000002:0.2):0.05,FLVV_000003:0.3);")
    id_map = {
        "FLVV_000001": "Dengue|s|A",
        "FLVV_000002": "Zika|s|B",
        "FLVV_000003": "YFV|s|C",
    }
    restore_newick_names(nwk_in, nwk_out, id_map)
    out = nwk_out.read_text()
    assert "Dengue|s|A" in out
    assert "Zika|s|B" in out
    assert "YFV|s|C" in out
    assert "FLVV_" not in out


# ---------------------------------------------------------------------------
# restore_fasta_names
# ---------------------------------------------------------------------------

def test_restore_fasta_names_replaces_ids(tmp_path):
    in_fasta = tmp_path / "in.fasta"
    out_fasta = tmp_path / "out.fasta"
    in_fasta.write_text(">FLVV_000001\nATGC\n>FLVV_000002\nGGCC\n")
    id_map = {"FLVV_000001": "Dengue|s|A", "FLVV_000002": "Zika|s|B"}
    restore_fasta_names(in_fasta, out_fasta, id_map)
    out = out_fasta.read_text()
    assert ">Dengue|s|A" in out
    assert ">Zika|s|B" in out
    assert "FLVV_" not in out
