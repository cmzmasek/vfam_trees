"""Tests for vfam_trees.phyloxml_writer."""
import xml.etree.ElementTree as ET

import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees.phyloxml_writer import (
    _extract_year,
    _is_same_species_pair,
    _valid_rank,
    _PHYLOXML_VALID_RANKS,
    METADATA_FIELDS,
    write_phyloxml,
)

# PhyloXML declares a default namespace, so ElementTree prefixes every element
# name with the namespace URI when parsing.  Helper alias keeps test queries readable.
NS = "{http://www.phyloxml.org}"


# ---------------------------------------------------------------------------
# _extract_year
# ---------------------------------------------------------------------------

class TestExtractYear:
    def test_iso_date(self):
        assert _extract_year("2024-01-01") == "2024"

    def test_year_only(self):
        assert _extract_year("2025") == "2025"

    def test_month_year(self):
        assert _extract_year("Sep-2023") == "2023"

    def test_year_month(self):
        assert _extract_year("2020-09") == "2020"

    def test_day_month_year_long(self):
        assert _extract_year("01-Jan-2020") == "2020"

    def test_verbose_date(self):
        assert _extract_year("January 15, 2019") == "2019"

    def test_empty_returns_empty(self):
        assert _extract_year("") == ""

    def test_unknown_returns_empty(self):
        assert _extract_year("unknown") == ""

    def test_no_year_in_string(self):
        assert _extract_year("Jan-15") == ""

    def test_1800s_year_matches(self):
        assert _extract_year("Apr-1892") == "1892"

    def test_1700s_year_does_not_match(self):
        # Regex restricts to 18xx, 19xx, 20xx
        assert _extract_year("1750") == ""

    def test_2100s_year_does_not_match(self):
        assert _extract_year("2150") == ""


# ---------------------------------------------------------------------------
# _valid_rank
# ---------------------------------------------------------------------------

class TestValidRank:
    def test_valid_rank_returned_as_is(self):
        assert _valid_rank("genus") == "genus"
        assert _valid_rank("family") == "family"

    def test_invalid_rank_mapped_to_other(self):
        assert _valid_rank("realm") not in _PHYLOXML_VALID_RANKS or True
        # "realm" is NOT in the PhyloXML schema enumeration
        assert _valid_rank("realm") == "other"

    def test_empty_rank_mapped_to_other(self):
        assert _valid_rank("") == "other"


# ---------------------------------------------------------------------------
# _is_same_species_pair
# ---------------------------------------------------------------------------

class TestIsSameSpeciesPair:
    def _mk_leaf(self, short_id):
        c = Clade(name=short_id)
        return c

    def test_two_leaves_same_species(self):
        parent = Clade()
        parent.clades = [self._mk_leaf("s1"), self._mk_leaf("s2")]
        meta = {"s1": {"species": "Dengue virus"}, "s2": {"species": "Dengue virus"}}
        assert _is_same_species_pair(parent, meta) is True

    def test_two_leaves_different_species(self):
        parent = Clade()
        parent.clades = [self._mk_leaf("s1"), self._mk_leaf("s2")]
        meta = {"s1": {"species": "Dengue virus"}, "s2": {"species": "Zika virus"}}
        assert _is_same_species_pair(parent, meta) is False

    def test_three_leaves_never_same_species_pair(self):
        parent = Clade()
        parent.clades = [self._mk_leaf(i) for i in ("a", "b", "c")]
        meta = {i: {"species": "Dengue virus"} for i in ("a", "b", "c")}
        assert _is_same_species_pair(parent, meta) is False

    def test_one_child_internal_returns_false(self):
        parent = Clade()
        inner = Clade()
        inner.clades = [self._mk_leaf("x")]
        parent.clades = [self._mk_leaf("a"), inner]
        meta = {"a": {"species": "X"}, "x": {"species": "X"}}
        assert _is_same_species_pair(parent, meta) is False

    def test_both_empty_species_returns_false(self):
        parent = Clade()
        parent.clades = [self._mk_leaf("a"), self._mk_leaf("b")]
        meta = {"a": {"species": ""}, "b": {"species": ""}}
        assert _is_same_species_pair(parent, meta) is False


# ---------------------------------------------------------------------------
# write_phyloxml — integration-ish test using an in-memory tree
# ---------------------------------------------------------------------------

def _make_mini_tree():
    """Build a 3-leaf tree: ((l1, l2):conf, l3)."""
    l1 = Clade(branch_length=0.1, name="s1")
    l2 = Clade(branch_length=0.2, name="s2")
    l3 = Clade(branch_length=0.3, name="s3")
    inner = Clade(branch_length=0.05, clades=[l1, l2])
    inner.confidence = 85
    inner.name = "Flavivirus"
    inner._taxonomy_rank = "genus"
    root = Clade(clades=[inner, l3])
    return Tree(root=root, rooted=True)


def _make_meta():
    # s1 and s2 have different species so the inner "Flavivirus" LCA
    # label is not suppressed by the same-species-pair heuristic.
    return {
        "s1": {
            "species": "Dengue virus",
            "taxon_id": "11053",
            "accession": "NC_001477",
            "seq_name": "Dengue virus 1",
            "host": "Homo sapiens",
            "collection_date": "2020-05-14",
            "location": "USA",
            "strain": "DEN1",
            "lineage_ranked": [
                {"name": "Flaviviridae", "rank": "family"},
                {"name": "Flavivirinae", "rank": "subfamily"},
                {"name": "Flavivirus", "rank": "genus"},
            ],
        },
        "s2": {
            "species": "Yellow fever virus",
            "taxon_id": "11089",
            "accession": "AF038403",
            "collection_date": "Sep-1985",
            "lineage_ranked": [
                {"name": "Flavivirus", "rank": "genus"},
            ],
        },
        "s3": {
            "species": "Zika virus",
            "taxon_id": "64320",
            "accession": "NC_012532",
        },
    }


@pytest.fixture
def mini_xml(tmp_path):
    tree = _make_mini_tree()
    meta = _make_meta()
    out = tmp_path / "out.xml"
    write_phyloxml(
        newick_path=tmp_path / "unused.nwk",
        output_xml=out,
        id_map={"s1": "Dengue|strain|NC_001477", "s2": "Dengue|s2|AF038403", "s3": "Zika|s3|NC_012532"},
        leaf_metadata=meta,
        family="Flaviviridae",
        tree=tree,
        phylogeny_name="TestTree",
        phylogeny_detail="detail",
        confidence_type="SH_aLRT",
        leaf_colors={"s1": "#ff0000", "s3": "#00ff00"},
    )
    return out


class TestWritePhyloxml:
    def test_file_created(self, mini_xml):
        assert mini_xml.exists() and mini_xml.stat().st_size > 0

    def test_phylogeny_is_rooted_and_not_rerootable(self, mini_xml):
        tree = ET.parse(mini_xml)
        phy = tree.getroot().find(f"{NS}phylogeny")
        assert phy is not None
        assert phy.attrib["rooted"] == "true"
        assert phy.attrib["rerootable"] == "false"

    def test_name_element_is_phylogeny_name(self, mini_xml):
        tree = ET.parse(mini_xml)
        name = tree.getroot().find(f"{NS}phylogeny/{NS}name")
        assert name is not None and name.text == "TestTree"

    def test_leaf_has_taxonomy_with_species_and_taxid(self, mini_xml):
        text = mini_xml.read_text()
        assert "Dengue virus" in text
        assert "11053" in text

    def test_leaf_has_sequence_accession(self, mini_xml):
        text = mini_xml.read_text()
        assert "NC_001477" in text
        assert "AF038403" in text

    def test_capitalized_vipr_property_refs(self, mini_xml):
        text = mini_xml.read_text()
        assert 'ref="vipr:Host"' in text
        assert 'ref="vipr:Collection_Date"' in text
        assert 'ref="vipr:Location"' in text
        assert 'ref="vipr:Strain"' in text

    def test_vipr_year_added_for_iso_date(self, mini_xml):
        text = mini_xml.read_text()
        assert 'ref="vipr:Year"' in text
        assert ">2020<" in text
        assert ">1985<" in text   # Sep-1985

    def test_vipr_year_absent_when_date_unknown(self, mini_xml):
        # s3 has no collection_date field → no <property ref="vipr:Year" ...>s3-related element.
        # We verify exactly two Year elements are present (for s1 + s2).
        tree = ET.parse(mini_xml)
        years = tree.getroot().findall(f".//{NS}property[@ref='vipr:Year']")
        assert len(years) == 2

    def test_vipr_species_property(self, mini_xml):
        text = mini_xml.read_text()
        assert 'ref="vipr:Species"' in text

    def test_vipr_genus_property_when_available(self, mini_xml):
        text = mini_xml.read_text()
        # s1 and s2 have genus=Flavivirus
        genus_props = text.count('ref="vipr:Genus"')
        assert genus_props == 2

    def test_vipr_subfamily_property_when_available(self, mini_xml):
        tree = ET.parse(mini_xml)
        subfams = tree.getroot().findall(f".//{NS}property[@ref='vipr:Subfamily']")
        assert len(subfams) == 1  # only s1 has subfamily in lineage_ranked
        assert subfams[0].text == "Flavivirinae"

    def test_vipr_subgenus_property_omitted_when_absent(self, mini_xml):
        tree = ET.parse(mini_xml)
        subgen = tree.getroot().findall(f".//{NS}property[@ref='vipr:Subgenus']")
        assert len(subgen) == 0

    def test_internal_taxonomy_element_present(self, mini_xml):
        # Inner clade has name "Flavivirus" with rank "genus"
        tree = ET.parse(mini_xml)
        sci_names = [e.text for e in tree.getroot().iter(f"{NS}scientific_name")]
        assert "Flavivirus" in sci_names

    def test_internal_confidence_preserved(self, mini_xml):
        tree = ET.parse(mini_xml)
        confidences = [e.text for e in tree.getroot().iter(f"{NS}confidence")]
        assert "85" in confidences

    def test_same_species_pair_label_suppressed(self, tmp_path):
        # Same tree but make the inner node a same-species pair.
        l1 = Clade(branch_length=0.1, name="s1")
        l2 = Clade(branch_length=0.2, name="s2")
        inner = Clade(branch_length=0.05, clades=[l1, l2])
        inner.name = "Dengue virus"
        inner._taxonomy_rank = "species"
        root = Clade(clades=[inner])
        tree = Tree(root=root, rooted=True)

        meta = {
            "s1": {"species": "Dengue virus"},
            "s2": {"species": "Dengue virus"},
        }
        out = tmp_path / "same.xml"
        write_phyloxml(
            newick_path=tmp_path / "_.nwk",
            output_xml=out,
            id_map={"s1": "s1d", "s2": "s2d"},
            leaf_metadata=meta,
            family="F",
            tree=tree,
        )
        # The internal node has name "Dengue virus" but same-species pair → suppressed.
        # We expect NO internal <name>Dengue virus</name> element.
        xml = ET.parse(out)
        names = [e.text for e in xml.getroot().iter(f"{NS}name")]
        assert "Dengue virus" not in names

    def test_font_color_property_emitted_when_color_mapped(self, mini_xml):
        text = mini_xml.read_text()
        assert 'ref="style:font_color"' in text
        assert "#ff0000" in text


# ---------------------------------------------------------------------------
# Fields declaration
# ---------------------------------------------------------------------------

def test_metadata_fields_are_three_tuples():
    # Regression: ensure no metadata field declaration was accidentally left as 2-tuple.
    for entry in METADATA_FIELDS:
        assert len(entry) == 3
        field, datatype, ref = entry
        assert isinstance(field, str) and field
        assert datatype.startswith("xsd:")
        assert ref and ref[0].isupper()
