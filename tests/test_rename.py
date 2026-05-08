"""Tests for vfam_trees.rename — short ID assignment and name restoration.

This file is the canonical regression test for the leaf-label naming
convention used across the entire pipeline (sequence FASTAs, PhyloXML
``<name>`` elements, Newick leaves, and rendered tree-image labels).

The naming format is locked down here:

    no host known  → ``<species>|<accession>``
    host known     → ``<species>|<accession>|<host>``

Single-protein and concat modes both delegate to :func:`canonical_leaf_label`
so identical ``(species, accession, host)`` inputs yield byte-identical
labels — the cross-mode consistency check below enforces that explicitly.

If you change the leaf naming convention, update the docstring of
``canonical_leaf_label`` in ``vfam_trees/rename.py`` AND update / add the
relevant tests here.  The format has regressed before; tests in this file
exist specifically to catch that.
"""
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees.rename import (
    _build_display_name,
    _family_prefix,
    assign_short_ids,
    canonical_leaf_label,
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
# canonical_leaf_label — single source of truth for leaf naming
# ---------------------------------------------------------------------------

class TestCanonicalLeafLabel:
    """Lock down the canonical leaf-label format.

    The naming convention is intentionally narrow:
      - exactly two pipe-separated components when host is absent
      - exactly three pipe-separated components when host is present
      - components are: species, accession, host (in that order)
      - spaces in any component become underscores
      - host treated as absent when None / "" / whitespace / "unknown"
    """

    def test_with_host_three_components_in_order(self):
        assert canonical_leaf_label("Dengue virus", "NC_001477", "Homo sapiens") \
            == "Dengue_virus|NC_001477|Homo_sapiens"

    def test_without_host_two_components_in_order(self):
        assert canonical_leaf_label("Dengue virus", "NC_001477", None) \
            == "Dengue_virus|NC_001477"

    @pytest.mark.parametrize("absent_host", [
        None, "", " ", "  ", "\t", "unknown", "Unknown", "UNKNOWN", " unknown ",
    ])
    def test_host_absent_variants_all_omit_host_component(self, absent_host):
        out = canonical_leaf_label("X virus", "ACC1", absent_host)
        assert out == "X_virus|ACC1"
        assert out.count("|") == 1   # exactly two components

    def test_host_present_always_yields_three_components(self):
        out = canonical_leaf_label("X virus", "ACC1", "Bat")
        assert out.count("|") == 2   # exactly three components
        assert out.split("|") == ["X_virus", "ACC1", "Bat"]

    def test_spaces_replaced_with_underscores_everywhere(self):
        out = canonical_leaf_label("Foot and mouth virus", "AB 12", "Bos taurus")
        assert " " not in out
        assert out == "Foot_and_mouth_virus|AB_12|Bos_taurus"

    def test_missing_species_becomes_unknown(self):
        assert canonical_leaf_label(None, "ACC1", None) == "unknown|ACC1"
        assert canonical_leaf_label("", "ACC1", None) == "unknown|ACC1"

    def test_missing_accession_becomes_unknown(self):
        assert canonical_leaf_label("Sp", None, None) == "Sp|unknown"
        assert canonical_leaf_label("Sp", "", None) == "Sp|unknown"

    def test_all_missing(self):
        assert canonical_leaf_label(None, None, None) == "unknown|unknown"

    def test_host_with_only_punctuation_and_spaces_kept(self):
        # Real edge case: hosts like "Aedes spp." should be kept verbatim
        # (only None/empty/whitespace/"unknown" are dropped).
        assert canonical_leaf_label("Sp", "ACC1", "Aedes spp.") \
            == "Sp|ACC1|Aedes_spp."

    def test_label_components_never_contain_pipes_in_inputs(self):
        # If species or host already contains a pipe (very rare), the label
        # would silently fragment.  Document this limitation: the canonical
        # builder does not sanitize internal pipes.  If this ever matters,
        # update both the function and this test deliberately.
        out = canonical_leaf_label("Weird|name", "ACC", None)
        assert out == "Weird|name|ACC"   # caller responsibility


# ---------------------------------------------------------------------------
# _build_display_name (single-protein adaptor) — must match canonical_leaf_label
# ---------------------------------------------------------------------------

class TestBuildDisplayName:
    def test_full_metadata(self):
        m = _meta()
        assert _build_display_name(m) == "Dengue_virus|NC_001477|Homo_sapiens"

    def test_host_unknown_omitted(self):
        m = _meta(host="unknown")
        assert _build_display_name(m) == "Dengue_virus|NC_001477"

    def test_host_empty_omitted(self):
        m = _meta(host="")
        assert _build_display_name(m) == "Dengue_virus|NC_001477"

    def test_host_none_omitted(self):
        m = _meta(host=None)
        assert _build_display_name(m) == "Dengue_virus|NC_001477"

    def test_missing_fields_replaced_by_unknown(self):
        m = {"accession": "ACC1"}
        assert _build_display_name(m) == "unknown|ACC1"

    def test_strain_field_ignored(self):
        # Strain is no longer part of the display name even if present in meta.
        m = _meta()
        m["strain"] = "DEN1"
        assert "DEN1" not in _build_display_name(m)

    def test_whitespace_replaced_with_underscore(self):
        m = _meta(species="Zika virus", host="Aedes aegypti")
        name = _build_display_name(m)
        assert "Zika_virus" in name
        assert "Aedes_aegypti" in name
        assert " " not in name

    def test_only_three_meta_keys_consulted(self):
        # Adaptor must consult only species / accession / host.  Any other
        # field (strain, country, year, ...) must NOT influence the label.
        baseline = _meta()
        contaminated = {
            **baseline,
            "strain":          "DEN1",
            "country":         "USA",
            "collection_year": "2024",
            "extra":           "should_not_appear",
        }
        assert _build_display_name(baseline) == _build_display_name(contaminated)

    def test_adaptor_matches_canonical_for_same_inputs(self):
        # Cross-check: feeding the same (species, accession, host) into
        # _build_display_name and canonical_leaf_label must yield the same
        # string.  This catches anyone re-introducing inline construction
        # in the single-protein adaptor.
        for species in ("Dengue virus", "X", "", None):
            for accession in ("NC_001", "", None):
                for host in ("Homo sapiens", "", None, "unknown"):
                    meta = {"species": species, "accession": accession, "host": host}
                    assert _build_display_name(meta) \
                        == canonical_leaf_label(species, accession, host)


# ---------------------------------------------------------------------------
# Concat display-name builder — must produce identical labels to single-protein
# ---------------------------------------------------------------------------

class TestConcatDisplayNamesMatchSingleProtein:
    """Both modes must produce the exact same canonical label for identical
    (species, accession, host) inputs.

    This test imports the concat helper directly and feeds it minimal inputs.
    The concat helper extracts host from a SeqRecord's source feature; we
    construct minimal SeqRecords here to exercise that path.
    """

    @staticmethod
    def _genome_with_host(host_value: str | None):
        """Build a one-marker genome dict whose source feature carries
        ``/host=<host_value>`` (or no host qualifier when None)."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("M"), id="prot1", description="")
        quals = {}
        if host_value is not None:
            quals["host"] = [host_value]
        rec.features = [SeqFeature(FeatureLocation(0, 1), type="source", qualifiers=quals)]
        rec.annotations = {}
        return {"polB": rec}

    def test_with_host_matches_single_protein_label(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome_with_host("Homo sapiens")}
        species_map = {"NC_001": "Dengue virus"}
        out = _build_concat_display_names(selected, species_map)
        meta = {"species": "Dengue virus", "accession": "NC_001", "host": "Homo sapiens"}
        assert out["NC_001"] == _build_display_name(meta) \
            == "Dengue_virus|NC_001|Homo_sapiens"

    def test_without_host_matches_single_protein_label(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome_with_host(None)}
        species_map = {"NC_001": "Dengue virus"}
        out = _build_concat_display_names(selected, species_map)
        meta = {"species": "Dengue virus", "accession": "NC_001", "host": None}
        assert out["NC_001"] == _build_display_name(meta) == "Dengue_virus|NC_001"

    @pytest.mark.parametrize("absent", ["", "unknown", "Unknown", " "])
    def test_absent_host_variants_drop_host_component(self, absent):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        selected = {"NC_001": self._genome_with_host(absent)}
        out = _build_concat_display_names(selected, {"NC_001": "Sp"})
        assert out["NC_001"] == "Sp|NC_001"
        assert out["NC_001"].count("|") == 1

    def test_empty_marker_dict_still_yields_label_without_host(self):
        # Defensive: a genome with zero markers should still get a label
        # (no host available → two-component form).
        from vfam_trees.pipeline_concat import _build_concat_display_names
        out = _build_concat_display_names({"NC_x": {}}, {"NC_x": "X virus"})
        assert out["NC_x"] == "X_virus|NC_x"

    def test_no_species_yields_unknown_species(self):
        from vfam_trees.pipeline_concat import _build_concat_display_names
        out = _build_concat_display_names({"NC_x": self._genome_with_host(None)}, {})
        assert out["NC_x"] == "unknown|NC_x"


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
#
# These exercise the full path from canonical_leaf_label → through whatever
# output writer the user sees.  Adding a new output channel?  Add a test
# here verifying the canonical name flows through unchanged.
# ---------------------------------------------------------------------------

class TestCanonicalNameFlowsThroughFastaOutput:
    """Sequence FASTAs must contain the canonical label as the record ID."""

    def test_with_host(self, tmp_path):
        rec = _rec("NC_001")
        meta = {"species": "Dengue virus", "accession": "NC_001", "host": "Homo sapiens"}
        id_map_path = tmp_path / "id_map.tsv"
        renamed, mapping = assign_short_ids([rec], [meta], "Flaviviridae", id_map_path)

        # 1. id_map carries the canonical name
        assert list(mapping.values()) == ["Dengue_virus|NC_001|Homo_sapiens"]

        # 2. Round-trip through restore_fasta_names produces it as the FASTA ID
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
    """PhyloXML ``<name>`` for external nodes must be the canonical label."""

    def _xml_leaf_names(self, xml_path):
        import xml.etree.ElementTree as ET
        ns = {"px": "http://www.phyloxml.org"}
        tree = ET.parse(str(xml_path))
        names = []
        for clade in tree.iter("{http://www.phyloxml.org}clade"):
            # External node = no nested <clade> children
            children = clade.findall("px:clade", ns)
            if not children:
                name_el = clade.find("px:name", ns)
                if name_el is not None and name_el.text:
                    names.append(name_el.text)
        return names

    def test_single_protein_phyloxml_name_with_host(self, tmp_path):
        from vfam_trees.phyloxml_writer import write_phyloxml
        # canonical label that the rest of the pipeline would produce
        canonical = canonical_leaf_label("Dengue virus", "NC_001", "Homo sapiens")
        nwk = tmp_path / "tree.nwk"
        nwk.write_text("(SHORT_001:0.1,SHORT_002:0.1);")
        xml_path = tmp_path / "tree.xml"
        write_phyloxml(
            newick_path=nwk,
            output_xml=xml_path,
            id_map={"SHORT_001": canonical, "SHORT_002": "Other_sp|ACC2"},
            leaf_metadata={
                "SHORT_001": {"accession": "NC_001"},
                "SHORT_002": {"accession": "ACC2"},
            },
            family="Flaviviridae",
            tree=None,
            phylogeny_name="Flaviviridae",
            phylogeny_detail="test",
            confidence_type="SH_like",
            leaf_colors={},
            aligned_seqs={},
            seq_type="protein",
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
            newick_path=nwk,
            output_xml=xml_path,
            id_map={"SHORT_001": canonical, "SHORT_002": "Other|ACC2"},
            leaf_metadata={"SHORT_001": {"accession": "NC_002"}, "SHORT_002": {}},
            family="Flaviviridae",
            tree=None,
            phylogeny_name="Flaviviridae",
            phylogeny_detail="test",
            confidence_type="SH_like",
            leaf_colors={},
            aligned_seqs={},
            seq_type="protein",
        )
        assert canonical in self._xml_leaf_names(xml_path)
        # exactly two pipe-separated components, no host
        assert canonical == "Zika_virus|NC_002"
