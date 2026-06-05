"""Tests for vfam_trees.rename — short ID assignment and name restoration.

Canonical regression test for the leaf-label naming convention used across
the entire pipeline (sequence FASTAs, PhyloXML ``<name>`` elements, Newick
leaves, and rendered tree-image labels).

Label format (default):

    ``{species}|{id}|{host}``

    — absent / "unknown" / "n/a" components (case-insensitive) are dropped
      and their preceding separator is suppressed so no leading or consecutive
      separators appear (``keep_separator_on_empty=False``).

Custom formats use the same ``{placeholder}`` syntax with any literal
separators; the same absent-field rule applies unless
``keep_separator_on_empty=True``.

If you change the leaf naming convention, update the docstring of
``format_leaf_label`` in ``vfam_trees/rename.py`` AND update / add the
relevant tests here.
"""
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees.rename import (
    DEFAULT_LABEL_FORMAT,
    LABEL_FIELDS,
    _build_display_name,
    _family_prefix,
    _normalize_field,
    _resolve_label_fields,
    assign_short_ids,
    canonical_leaf_label,
    format_leaf_label,
    load_id_map,
    restore_fasta_names,
)


def _rec(acc, seq="ATGC"):
    return SeqRecord(Seq(seq), id=acc, description="raw desc")


def _meta(species="Dengue virus", accession="NC_001477", host="Homo sapiens"):
    return {"species": species, "accession": accession, "host": host}


# ---------------------------------------------------------------------------
# _family_prefix
# ---------------------------------------------------------------------------

class TestFamilyPrefix:
    def test_four_consonants(self):
        assert _family_prefix("Flaviviridae") == "FLVV"

    def test_first_four_consonants_only(self):
        assert _family_prefix("Coronaviridae") == "CRNV"

    def test_too_few_consonants_uses_first_four_chars(self):
        assert _family_prefix("Abc") == "ABC"

    def test_vowel_heavy_name_fallback(self):
        assert _family_prefix("Aeiou") == "AEIO"


# ---------------------------------------------------------------------------
# _normalize_field
# ---------------------------------------------------------------------------

class TestNormalizeField:
    @pytest.mark.parametrize("absent", [None, "", " ", "  ", "unknown", "Unknown",
                                         "UNKNOWN", "n/a", "N/A"])
    def test_absent_values_become_empty_string(self, absent):
        assert _normalize_field(absent) == ""

    def test_real_value_preserved(self):
        assert _normalize_field("Dengue virus") == "Dengue virus"

    def test_strips_surrounding_whitespace(self):
        assert _normalize_field("  Dengue  ") == "Dengue"


# ---------------------------------------------------------------------------
# _resolve_label_fields
# ---------------------------------------------------------------------------

class TestResolveLabelFields:
    def test_all_standard_fields_present(self):
        meta = {
            "species": "Dengue virus", "accession": "NC_001477",
            "host": "Homo sapiens", "strain": "DEN-1",
            "location": "Thailand", "collection_date": "2005-07-14",
        }
        f = _resolve_label_fields(meta)
        assert f["species"] == "Dengue virus"
        assert f["id"]      == "NC_001477"
        assert f["host"]    == "Homo sapiens"
        assert f["strain"]  == "DEN-1"
        assert f["location"] == "Thailand"
        assert f["year"]    == "2005"
        assert f["genus"]   == ""

    def test_year_extracted_from_full_date(self):
        assert _resolve_label_fields({"collection_date": "2020-03-15"})["year"] == "2020"

    def test_year_extracted_from_year_only(self):
        assert _resolve_label_fields({"collection_date": "1999"})["year"] == "1999"

    def test_year_empty_when_no_date(self):
        assert _resolve_label_fields({})["year"] == ""

    def test_year_empty_when_date_unknown(self):
        assert _resolve_label_fields({"collection_date": "unknown"})["year"] == ""

    def test_absent_fields_normalized_to_empty(self):
        f = _resolve_label_fields({"accession": "ACC1", "host": "unknown"})
        assert f["species"] == ""
        assert f["host"]    == ""

    def test_all_label_fields_are_present_in_output(self):
        f = _resolve_label_fields({})
        assert set(f.keys()) >= LABEL_FIELDS


# ---------------------------------------------------------------------------
# format_leaf_label — the core engine
# ---------------------------------------------------------------------------

class TestFormatLeafLabel:
    """Lock down format_leaf_label rendering for all combinations."""

    def test_default_format_all_fields_present(self):
        f = {"species": "Dengue virus", "id": "NC_001477", "host": "Homo sapiens"}
        assert format_leaf_label(DEFAULT_LABEL_FORMAT, f) == "Dengue_virus|NC_001477|Homo_sapiens"

    def test_default_format_host_absent_two_components(self):
        f = {"species": "Dengue virus", "id": "NC_001477", "host": ""}
        assert format_leaf_label(DEFAULT_LABEL_FORMAT, f) == "Dengue_virus|NC_001477"

    def test_custom_format_species_strain_id(self):
        f = {"species": "Dengue virus", "strain": "DEN-1", "id": "NC_001477"}
        assert format_leaf_label("{species}|{strain}|{id}", f) == "Dengue_virus|DEN-1|NC_001477"

    def test_custom_separator_not_pipe(self):
        f = {"species": "Dengue virus", "id": "NC_001477"}
        assert format_leaf_label("{species}_{id}", f) == "Dengue_virus_NC_001477"

    def test_empty_middle_field_drops_its_separator(self):
        f = {"species": "Dengue virus", "strain": "", "id": "NC_001477"}
        out = format_leaf_label("{species}|{strain}|{id}", f)
        assert out == "Dengue_virus|NC_001477"
        assert "||" not in out

    def test_empty_first_field_drops_leading_separator_of_second(self):
        f = {"strain": "", "species": "Dengue virus", "id": "NC_001477"}
        out = format_leaf_label("{strain}|{species}|{id}", f)
        assert out == "Dengue_virus|NC_001477"
        assert out.startswith("|") is False

    def test_all_fields_empty_returns_empty_string(self):
        assert format_leaf_label("{species}|{id}", {"species": "", "id": ""}) == ""

    def test_keep_separator_on_empty_true_preserves_all_seps(self):
        f = {"species": "Dengue virus", "strain": "", "id": "NC_001477"}
        out = format_leaf_label("{species}|{strain}|{id}", f, keep_separator_on_empty=True)
        assert out == "Dengue_virus||NC_001477"

    def test_replace_whitespace_true_spaces_become_underscores(self):
        f = {"species": "Dengue virus", "id": "NC 001"}
        assert " " not in format_leaf_label("{species}|{id}", f, replace_whitespace=True)

    def test_replace_whitespace_false_spaces_kept(self):
        f = {"species": "Dengue virus", "id": "NC_001"}
        out = format_leaf_label("{species}|{id}", f, replace_whitespace=False)
        assert "Dengue virus" in out

    def test_year_placeholder_uses_year_value(self):
        f = {"species": "Dengue virus", "id": "NC_001", "year": "2020"}
        assert format_leaf_label("{species}|{id}|{year}", f) == "Dengue_virus|NC_001|2020"

    def test_location_placeholder(self):
        f = {"species": "Sp", "id": "ACC", "location": "Thailand"}
        assert format_leaf_label("{species}|{id}|{location}", f) == "Sp|ACC|Thailand"

    def test_trailing_literal_always_appended(self):
        f = {"species": "Sp", "id": "ACC"}
        out = format_leaf_label("{species}|{id}_v1", f)
        assert out.endswith("_v1")

    def test_unknown_field_name_renders_empty(self):
        f = {"species": "Sp", "id": "ACC"}
        out = format_leaf_label("{species}|{bogus}|{id}", f)
        # {bogus} is absent → empty → its separator dropped
        assert out == "Sp|ACC"

    @pytest.mark.parametrize("absent", ["", "unknown", "Unknown", "n/a", "N/A"])
    def test_absent_values_drop_field_and_separator(self, absent):
        f = {"species": "Sp", "strain": absent, "id": "ACC"}
        out = format_leaf_label("{species}|{strain}|{id}", f)
        assert out == "Sp|ACC"


# ---------------------------------------------------------------------------
# canonical_leaf_label — default format with normalized absent values
# ---------------------------------------------------------------------------

class TestCanonicalLeafLabel:
    """Lock down the canonical default-format label behaviour."""

    def test_with_host_three_components_in_order(self):
        assert canonical_leaf_label("Dengue virus", "NC_001477", "Homo sapiens") \
            == "Dengue_virus|NC_001477|Homo_sapiens"

    def test_without_host_two_components_in_order(self):
        assert canonical_leaf_label("Dengue virus", "NC_001477", None) \
            == "Dengue_virus|NC_001477"

    @pytest.mark.parametrize("absent_host", [
        None, "", " ", "  ", "\t", "unknown", "Unknown", "UNKNOWN", " unknown ",
        "n/a", "N/A",
    ])
    def test_host_absent_variants_all_omit_host_component(self, absent_host):
        out = canonical_leaf_label("X virus", "ACC1", absent_host)
        assert out == "X_virus|ACC1"
        assert out.count("|") == 1

    def test_host_present_always_yields_three_components(self):
        out = canonical_leaf_label("X virus", "ACC1", "Bat")
        assert out.count("|") == 2
        assert out.split("|") == ["X_virus", "ACC1", "Bat"]

    def test_spaces_replaced_with_underscores_everywhere(self):
        out = canonical_leaf_label("Foot and mouth virus", "AB 12", "Bos taurus")
        assert " " not in out
        assert out == "Foot_and_mouth_virus|AB_12|Bos_taurus"

    def test_missing_species_omitted(self):
        # Absent species → just accession (no "unknown|" prefix)
        assert canonical_leaf_label(None, "ACC1", None) == "ACC1"
        assert canonical_leaf_label("", "ACC1", None) == "ACC1"

    def test_missing_accession_omitted(self):
        assert canonical_leaf_label("Sp", None, None) == "Sp"
        assert canonical_leaf_label("Sp", "", None) == "Sp"

    def test_all_missing_returns_empty(self):
        assert canonical_leaf_label(None, None, None) == ""

    def test_host_with_only_punctuation_and_spaces_kept(self):
        assert canonical_leaf_label("Sp", "ACC1", "Aedes spp.") \
            == "Sp|ACC1|Aedes_spp."


# ---------------------------------------------------------------------------
# _build_display_name (single-protein adaptor)
# ---------------------------------------------------------------------------

class TestBuildDisplayName:
    def test_full_metadata(self):
        m = _meta()
        assert _build_display_name(m) == "Dengue_virus|NC_001477|Homo_sapiens"

    def test_host_unknown_omitted(self):
        assert _build_display_name(_meta(host="unknown")) == "Dengue_virus|NC_001477"

    def test_host_empty_omitted(self):
        assert _build_display_name(_meta(host="")) == "Dengue_virus|NC_001477"

    def test_host_none_omitted(self):
        assert _build_display_name(_meta(host=None)) == "Dengue_virus|NC_001477"

    def test_missing_fields_omitted(self):
        m = {"accession": "ACC1"}
        assert _build_display_name(m) == "ACC1"

    def test_whitespace_replaced_with_underscore(self):
        name = _build_display_name(_meta(species="Zika virus", host="Aedes aegypti"))
        assert "Zika_virus" in name
        assert "Aedes_aegypti" in name
        assert " " not in name

    def test_adaptor_matches_canonical_for_same_inputs(self):
        for species in ("Dengue virus", "X", "", None):
            for accession in ("NC_001", "", None):
                for host in ("Homo sapiens", "", None, "unknown"):
                    meta = {"species": species, "accession": accession, "host": host}
                    assert _build_display_name(meta) \
                        == canonical_leaf_label(species, accession, host)


# ---------------------------------------------------------------------------
# Concat display-name builder — must match single-protein for same inputs
# ---------------------------------------------------------------------------

class TestConcatDisplayNamesMatchSingleProtein:
    """Both modes must produce the exact same canonical label for identical inputs."""

    @staticmethod
    def _genome(host=None, strain=None, location=None, date=None):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("M"), id="prot1", description="")
        quals = {}
        if host     is not None: quals["host"]            = [host]
        if strain   is not None: quals["strain"]          = [strain]
        if location is not None: quals["country"]         = [location]
        if date     is not None: quals["collection_date"] = [date]
        rec.features = [SeqFeature(FeatureLocation(0, 1), type="source", qualifiers=quals)]
        rec.annotations = {}
        return {"polB": rec}

    def test_default_format_with_host(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome(host="Homo sapiens")}
        out = _build_concat_display_names(selected, {"NC_001": "Dengue virus"})
        assert out["NC_001"] == "Dengue_virus|NC_001|Homo_sapiens"

    def test_default_format_without_host(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome()}
        out = _build_concat_display_names(selected, {"NC_001": "Dengue virus"})
        assert out["NC_001"] == "Dengue_virus|NC_001"

    @pytest.mark.parametrize("absent", ["", "unknown", "Unknown", " "])
    def test_absent_host_variants_drop_host_component(self, absent):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome(host=absent)}
        out = _build_concat_display_names(selected, {"NC_001": "Sp"})
        assert out["NC_001"] == "Sp|NC_001"
        assert out["NC_001"].count("|") == 1

    def test_empty_marker_dict_yields_label_without_host(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        out = _build_concat_display_names({"NC_x": {}}, {"NC_x": "X virus"})
        assert out["NC_x"] == "X_virus|NC_x"

    def test_no_species_yields_just_accession(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        out = _build_concat_display_names({"NC_x": self._genome()}, {})
        assert out["NC_x"] == "NC_x"

    def test_custom_format_strain_included(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome(host="Homo sapiens", strain="DEN-1")}
        out = _build_concat_display_names(
            selected, {"NC_001": "Dengue virus"},
            label_format="{species}|{strain}|{id}|{host}",
        )
        assert out["NC_001"] == "Dengue_virus|DEN-1|NC_001|Homo_sapiens"

    def test_custom_format_year_and_location(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome(location="Thailand", date="2005-07-14")}
        out = _build_concat_display_names(
            selected, {"NC_001": "Dengue virus"},
            label_format="{species}|{id}|{location}|{year}",
        )
        assert out["NC_001"] == "Dengue_virus|NC_001|Thailand|2005"

    def test_custom_format_empty_field_dropped(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome()}  # no strain
        out = _build_concat_display_names(
            selected, {"NC_001": "Dengue virus"},
            label_format="{species}|{strain}|{id}",
        )
        assert out["NC_001"] == "Dengue_virus|NC_001"

    def test_keep_separator_on_empty_true(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome()}  # no strain
        out = _build_concat_display_names(
            selected, {"NC_001": "Dengue virus"},
            label_format="{species}|{strain}|{id}",
            keep_separator_on_empty=True,
        )
        assert out["NC_001"] == "Dengue_virus||NC_001"


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
        id_map_path = tmp_path / "id_map.tsv"
        renamed, _ = assign_short_ids([rec], [_meta(accession="ACC1")], "Flaviviridae", id_map_path)
        assert str(renamed[0].seq) == "ATGCATGC"

    def test_description_cleared_on_renamed_record(self, tmp_path):
        rec = _rec("ACC1")
        rec.description = "some long description"
        id_map_path = tmp_path / "id_map.tsv"
        renamed, _ = assign_short_ids([rec], [_meta(accession="ACC1")], "Flaviviridae", id_map_path)
        assert renamed[0].description == ""

    def test_custom_format_strain_in_display_name(self, tmp_path):
        meta = {**_meta(accession="NC_001477", host=None), "strain": "DEN-1"}
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("NC_001477")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{strain}|{id}",
        )
        assert list(s2d.values()) == ["Dengue_virus|DEN-1|NC_001477"]

    def test_custom_format_year_extracted(self, tmp_path):
        meta = {**_meta(accession="NC_001477", host=None), "collection_date": "2005-07-14"}
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("NC_001477")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{id}|{year}",
        )
        assert list(s2d.values()) == ["Dengue_virus|NC_001477|2005"]

    def test_default_format_matches_canonical_leaf_label(self, tmp_path):
        meta = _meta(accession="NC_001477")
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids([_rec("NC_001477")], [meta], "Flaviviridae", id_map_path)
        expected = canonical_leaf_label(meta["species"], meta["accession"], meta["host"])
        assert list(s2d.values()) == [expected]

    def test_name_prefix_prepended_verbatim(self, tmp_path):
        meta = {**_meta(accession="PP_1", host=None), "name_prefix": "PATHOPLEXUS_"}
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("PP_1")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{id}",
        )
        assert list(s2d.values()) == ["PATHOPLEXUS_Dengue_virus|PP_1"]

    def test_absent_name_prefix_leaves_label_unchanged(self, tmp_path):
        meta = _meta(accession="NC_001477", host=None)  # no name_prefix key
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("NC_001477")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{id}",
        )
        assert list(s2d.values()) == ["Dengue_virus|NC_001477"]

    def test_empty_field_not_in_label_by_default(self, tmp_path):
        meta = {**_meta(accession="NC_001477", host=None), "strain": "unknown"}
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("NC_001477")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{strain}|{id}",
        )
        label = list(s2d.values())[0]
        assert label == "Dengue_virus|NC_001477"
        assert label.count("|") == 1

    def test_keep_separator_on_empty_produces_double_pipe(self, tmp_path):
        meta = {**_meta(accession="NC_001477", host=None), "strain": ""}
        id_map_path = tmp_path / "id_map.tsv"
        _, s2d = assign_short_ids(
            [_rec("NC_001477")], [meta], "Flaviviridae", id_map_path,
            label_format="{species}|{strain}|{id}",
            keep_separator_on_empty=True,
        )
        assert list(s2d.values()) == ["Dengue_virus||NC_001477"]


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


# ---------------------------------------------------------------------------
# End-to-end naming-consistency tests
# ---------------------------------------------------------------------------

class TestCanonicalNameFlowsThroughFastaOutput:
    def test_with_host(self, tmp_path):
        rec = _rec("NC_001")
        meta = {"species": "Dengue virus", "accession": "NC_001", "host": "Homo sapiens"}
        id_map_path = tmp_path / "id_map.tsv"
        renamed, mapping = assign_short_ids([rec], [meta], "Flaviviridae", id_map_path)
        assert list(mapping.values()) == ["Dengue_virus|NC_001|Homo_sapiens"]
        in_fa = tmp_path / "short.fasta"
        in_fa.write_text(f">{renamed[0].id}\nATGC\n")
        out_fa = tmp_path / "display.fasta"
        restore_fasta_names(in_fa, out_fa, mapping)
        assert ">Dengue_virus|NC_001|Homo_sapiens" in out_fa.read_text()

    def test_without_host(self, tmp_path):
        rec = _rec("NC_002")
        meta = {"species": "Zika virus", "accession": "NC_002", "host": None}
        id_map_path = tmp_path / "id_map.tsv"
        renamed, mapping = assign_short_ids([rec], [meta], "Flaviviridae", id_map_path)
        assert list(mapping.values()) == ["Zika_virus|NC_002"]
        in_fa = tmp_path / "short.fasta"
        in_fa.write_text(f">{renamed[0].id}\nATGC\n")
        out_fa = tmp_path / "display.fasta"
        restore_fasta_names(in_fa, out_fa, mapping)
        assert ">Zika_virus|NC_002\n" in out_fa.read_text()


class TestCanonicalNameFlowsThroughPhyloXML:
    def _xml_leaf_names(self, xml_path):
        import xml.etree.ElementTree as ET
        ns = {"px": "http://www.phyloxml.org"}
        tree = ET.parse(str(xml_path))
        names = []
        for clade in tree.iter("{http://www.phyloxml.org}clade"):
            children = clade.findall("px:clade", ns)
            if not children:
                name_el = clade.find("px:name", ns)
                if name_el is not None and name_el.text:
                    names.append(name_el.text)
        return names

    def test_single_protein_phyloxml_name_with_host(self, tmp_path):
        from vfam_trees.phyloxml_writer import write_phyloxml
        canonical = canonical_leaf_label("Dengue virus", "NC_001", "Homo sapiens")
        nwk = tmp_path / "tree.nwk"
        nwk.write_text("(SHORT_001:0.1,SHORT_002:0.1);")
        xml_path = tmp_path / "tree.xml"
        write_phyloxml(
            newick_path=nwk, output_xml=xml_path,
            id_map={"SHORT_001": canonical, "SHORT_002": "Other_sp|ACC2"},
            leaf_metadata={"SHORT_001": {"accession": "NC_001"}, "SHORT_002": {"accession": "ACC2"}},
            family="Flaviviridae", tree=None, phylogeny_name="Flaviviridae",
            phylogeny_detail="test", confidence_type="SH_like",
            leaf_colors={}, aligned_seqs={}, seq_type="protein",
        )
        names = self._xml_leaf_names(xml_path)
        assert canonical in names
        assert "Dengue_virus|NC_001|Homo_sapiens" in names

    def test_single_protein_phyloxml_name_without_host(self, tmp_path):
        from vfam_trees.phyloxml_writer import write_phyloxml
        canonical = canonical_leaf_label("Zika virus", "NC_002", None)
        nwk = tmp_path / "tree.nwk"
        nwk.write_text("(SHORT_001:0.1,SHORT_002:0.1);")
        xml_path = tmp_path / "tree.xml"
        write_phyloxml(
            newick_path=nwk, output_xml=xml_path,
            id_map={"SHORT_001": canonical, "SHORT_002": "Other|ACC2"},
            leaf_metadata={"SHORT_001": {"accession": "NC_002"}, "SHORT_002": {}},
            family="Flaviviridae", tree=None, phylogeny_name="Flaviviridae",
            phylogeny_detail="test", confidence_type="SH_like",
            leaf_colors={}, aligned_seqs={}, seq_type="protein",
        )
        assert canonical in self._xml_leaf_names(xml_path)
        assert canonical == "Zika_virus|NC_002"
