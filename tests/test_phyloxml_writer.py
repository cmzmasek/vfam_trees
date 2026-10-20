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

    def test_same_ictv_species_different_organism(self):
        # Leaves show vipr:Species = ICTV species, so a pair sharing it is
        # redundant even when the GenBank organism strings differ.
        lineage = [{"name": "Norovirus", "rank": "genus"},
                   {"name": "Norovirus norwalkense", "rank": "species"}]
        parent = Clade()
        parent.clades = [self._mk_leaf("a"), self._mk_leaf("b")]
        meta = {"a": {"species": "Norovirus GII.4", "lineage_ranked": lineage},
                "b": {"species": "Norovirus GII.17", "lineage_ranked": lineage}}
        assert _is_same_species_pair(parent, meta) is True

    def test_placeholder_species_is_not_a_match(self):
        parent = Clade()
        parent.clades = [self._mk_leaf("a"), self._mk_leaf("b")]
        meta = {"a": {"species": "unknown"}, "b": {"species": "unknown"}}
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
        assert 'ref="vipr:Country"' in text
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

    def test_vipr_outbreak_property_emitted_only_for_tagged_leaf(self, tmp_path):
        # Only s1 carries an outbreak tag; s2/s3 must not get the property.
        meta = _make_meta()
        meta["s1"]["outbreak"] = "Bdbv-2026"
        out = tmp_path / "outbreak.xml"
        write_phyloxml(
            newick_path=tmp_path / "unused.nwk",
            output_xml=out,
            id_map={"s1": "x", "s2": "y", "s3": "z"},
            leaf_metadata=meta,
            family="Flaviviridae",
            tree=_make_mini_tree(),
        )
        props = ET.parse(out).getroot().findall(
            f".//{NS}property[@ref='vipr:outbreak']"
        )
        assert len(props) == 1
        assert props[0].text == "Bdbv-2026"
        assert props[0].attrib["datatype"] == "xsd:string"

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

    def test_no_mol_seq_by_default(self, mini_xml):
        text = mini_xml.read_text()
        assert "mol_seq" not in text

    def test_no_type_attr_on_sequence_by_default(self, mini_xml):
        xml = ET.parse(mini_xml)
        for seq_el in xml.getroot().iter(f"{NS}sequence"):
            assert "type" not in seq_el.attrib


class TestMolSeq:
    """Aligned-sequence embedding via aligned_seqs + seq_type."""

    def _write(self, tmp_path, aligned_seqs=None, seq_type=""):
        tree = _make_mini_tree()
        meta = _make_meta()
        out = tmp_path / "out.xml"
        write_phyloxml(
            newick_path=tmp_path / "unused.nwk",
            output_xml=out,
            id_map={"s1": "Leaf1", "s2": "Leaf2", "s3": "Leaf3"},
            leaf_metadata=meta,
            family="Flaviviridae",
            tree=tree,
            aligned_seqs=aligned_seqs,
            seq_type=seq_type,
        )
        return out

    def test_mol_seq_present_for_each_leaf(self, tmp_path):
        seqs = {"s1": "ACDEF-GH", "s2": "ACDE--GH", "s3": "ACDEK-GH"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="protein")
        xml = ET.parse(out)
        mol_seqs = list(xml.getroot().iter(f"{NS}mol_seq"))
        assert len(mol_seqs) == 3
        texts = {e.text for e in mol_seqs}
        assert "ACDEF-GH" in texts
        assert "ACDEK-GH" in texts

    def test_mol_seq_is_aligned_attribute_true(self, tmp_path):
        seqs = {"s1": "ACDEF", "s2": "ACDE-", "s3": "ACDE-"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="protein")
        xml = ET.parse(out)
        for mol_el in xml.getroot().iter(f"{NS}mol_seq"):
            assert mol_el.attrib.get("is_aligned") == "true"

    def test_sequence_type_protein_for_protein(self, tmp_path):
        seqs = {"s1": "ACDEF", "s2": "ACDE-", "s3": "ACDE-"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="protein")
        xml = ET.parse(out)
        for seq_el in xml.getroot().iter(f"{NS}sequence"):
            assert seq_el.attrib.get("type") == "protein"

    def test_sequence_type_dna_for_nucleotide(self, tmp_path):
        seqs = {"s1": "ATCG--", "s2": "ATCGAT", "s3": "ATCG--"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="nucleotide")
        xml = ET.parse(out)
        for seq_el in xml.getroot().iter(f"{NS}sequence"):
            assert seq_el.attrib.get("type") == "dna"

    def test_mol_seq_absent_when_aligned_seqs_none(self, tmp_path):
        out = self._write(tmp_path, aligned_seqs=None, seq_type="protein")
        text = out.read_text()
        assert "mol_seq" not in text

    def test_mol_seq_absent_for_leaf_missing_from_dict(self, tmp_path):
        # Only s1 has an aligned sequence; s2 and s3 should not get mol_seq.
        seqs = {"s1": "ACDEF"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="protein")
        xml = ET.parse(out)
        mol_seqs = list(xml.getroot().iter(f"{NS}mol_seq"))
        assert len(mol_seqs) == 1
        assert mol_seqs[0].text == "ACDEF"

    def test_accession_and_mol_seq_coexist_in_same_sequence_el(self, tmp_path):
        seqs = {"s1": "ACDEF", "s2": "ACDE-", "s3": "ACDE-"}
        out = self._write(tmp_path, aligned_seqs=seqs, seq_type="protein")
        xml = ET.parse(out)
        # Each <sequence> that has an accession should also have a mol_seq child
        for seq_el in xml.getroot().iter(f"{NS}sequence"):
            acc = seq_el.find(f"{NS}accession")
            mol = seq_el.find(f"{NS}mol_seq")
            if acc is not None:
                assert mol is not None, "accession without mol_seq"


# ---------------------------------------------------------------------------
# Fields declaration
# ---------------------------------------------------------------------------

def test_metadata_fields_are_three_tuples():
    # Regression: ensure no metadata field declaration was accidentally left as 2-tuple.
    # Ref names are capitalised (vipr:Host, vipr:Location, ...) except "outbreak",
    # which is intentionally lowercase to mirror the additional_data_sources field.
    _LOWERCASE_REFS = {"outbreak"}
    for entry in METADATA_FIELDS:
        assert len(entry) == 3
        field, datatype, ref = entry
        assert isinstance(field, str) and field
        assert datatype.startswith("xsd:")
        assert ref and (ref[0].isupper() or ref in _LOWERCASE_REFS)


# ---------------------------------------------------------------------------
# Archaeopteryx / Archaeopteryx.js integration
# ---------------------------------------------------------------------------

class TestViewerIntegration:
    """What the Archaeopteryx viewers read: unit, tip dates, taxa, sources."""

    def _meta(self):
        meta = _make_meta()
        # s1: own taxon = species node (taxid matches the last lineage entry)
        meta["s1"]["lineage_ranked"] = [
            {"name": "Flaviviridae", "rank": "family", "taxid": "11050"},
            {"name": "Orthoflavivirus", "rank": "genus", "taxid": "3044782"},
            {"name": "Orthoflavivirus denguei", "rank": "species", "taxid": "3052464"},
            {"name": "Dengue virus", "rank": "no rank", "taxid": "11053"},
        ]
        meta["s1"]["location"] = "USA: Texas, Austin"
        # s2: Pathoplexus record
        meta["s2"]["source"] = "pathoplexus"
        meta["s2"]["accession"] = "PP_00711S5.2"
        # s3: pasted sequence, no lineage
        meta["s3"]["source"] = "user"
        return meta

    def _write(self, tmp_path, tree=None):
        out = tmp_path / "out.xml"
        write_phyloxml(
            newick_path=tmp_path / "unused.nwk", output_xml=out,
            id_map={"s1": "a", "s2": "b", "s3": "c"}, leaf_metadata=self._meta(),
            family="Flaviviridae", tree=tree or _make_mini_tree(),
        )
        return ET.parse(out).getroot()

    def _leaf(self, root, name):
        for clade in root.iter(f"{NS}clade"):
            n = clade.find(f"{NS}name")
            if n is not None and n.text == name and clade.find(f"{NS}clade") is None:
                return clade
        raise AssertionError(name)

    def _props(self, clade):
        return {p.attrib["ref"]: p.text for p in clade.findall(f"{NS}property")}

    def test_branch_length_unit_and_attribute_order(self, tmp_path):
        phy = self._write(tmp_path).find(f"{NS}phylogeny")
        assert phy.attrib["branch_length_unit"] == "subs/site"
        assert list(phy.attrib) == ["rooted", "branch_length_unit", "rerootable"]

    def test_tip_date_decimal_year_with_raw_desc(self, tmp_path):
        root = self._write(tmp_path)
        d1 = self._leaf(root, "a").find(f"{NS}date")
        assert d1.attrib["unit"] == "year"
        assert d1.find(f"{NS}desc").text == "2020-05-14"
        # 14 May 2020 = day 135 of a leap year, mid-day
        assert d1.find(f"{NS}value").text == f"{2020 + 134.5 / 366:.5f}"
        assert d1.find(f"{NS}minimum") is None and d1.find(f"{NS}maximum") is None
        assert self._leaf(root, "b").find(f"{NS}date") is not None   # Sep-1985
        assert self._leaf(root, "c").find(f"{NS}date") is None       # no date

    def test_date_sits_between_sequence_and_properties(self, tmp_path):
        tags = [c.tag.replace(NS, "") for c in self._leaf(self._write(tmp_path), "a")]
        assert tags.index("sequence") < tags.index("date") < tags.index("property")

    def test_leaf_taxonomy_is_own_ncbi_taxon(self, tmp_path):
        tax = self._leaf(self._write(tmp_path), "a").find(f"{NS}taxonomy")
        assert tax.find(f"{NS}id").text == "11053"
        assert tax.find(f"{NS}scientific_name").text == "Dengue virus"
        assert tax.find(f"{NS}rank") is None          # NCBI "no rank"

    def test_leaf_taxonomy_rank_written_when_ncbi_rank_is_valid(self, tmp_path):
        meta = self._meta()
        meta["s1"]["taxon_id"] = "3052464"
        meta["s1"]["lineage_ranked"] = meta["s1"]["lineage_ranked"][:3]
        out = tmp_path / "out.xml"
        write_phyloxml(newick_path=tmp_path / "x.nwk", output_xml=out,
                       id_map={"s1": "a", "s2": "b", "s3": "c"}, leaf_metadata=meta,
                       family="Flaviviridae", tree=_make_mini_tree())
        tax = self._leaf(ET.parse(out).getroot(), "a").find(f"{NS}taxonomy")
        assert tax.find(f"{NS}scientific_name").text == "Orthoflavivirus denguei"
        assert tax.find(f"{NS}rank").text == "species"

    def test_vipr_species_is_ictv_species(self, tmp_path):
        root = self._write(tmp_path)
        assert self._props(self._leaf(root, "a"))["vipr:Species"] == "Orthoflavivirus denguei"
        # s2: ranked lineage without a species rank -> no pseudo-species
        assert "vipr:Species" not in self._props(self._leaf(root, "b"))
        # s3: no ranked lineage -> organism fallback
        assert self._props(self._leaf(root, "c"))["vipr:Species"] == "Zika virus"

    def test_country_and_location_detail(self, tmp_path):
        props = self._props(self._leaf(self._write(tmp_path), "a"))
        assert props["vipr:Country"] == "USA"
        assert props["vipr:Location"] == "USA: Texas, Austin"

    def test_accession_sources(self, tmp_path):
        root = self._write(tmp_path)
        acc = {n: self._leaf(root, n).find(f"{NS}sequence/{NS}accession")
               for n in ("a", "b", "c")}
        assert acc["a"].attrib["source"] == "ncbi"
        assert acc["b"].attrib["source"] == "pathoplexus"
        assert acc["c"].attrib["source"] == "user"
        uri = self._leaf(root, "b").find(f"{NS}sequence/{NS}uri")
        assert uri.text == "https://pathoplexus.org/seq/PP_00711S5.2"
        assert self._leaf(root, "a").find(f"{NS}sequence/{NS}uri") is None
        assert self._leaf(root, "c").find(f"{NS}sequence/{NS}uri") is None

    def test_no_font_color_without_leaf_colors(self, tmp_path):
        root = self._write(tmp_path)
        assert all("style:font_color" not in self._props(c) for c in root.iter(f"{NS}clade"))

    def test_internal_ncbi_no_rank_is_not_inferred_from_suffix(self, tmp_path):
        tree = _make_mini_tree()
        inner = tree.root.clades[0]
        inner.name = "unclassified Flaviviridae"
        inner._taxonomy_rank = ""          # NCBI: no rank
        inner._taxonomy_id = "38144"
        tax = self._write(tmp_path, tree=tree).find(f"{NS}phylogeny/{NS}clade/{NS}clade/{NS}taxonomy")
        assert tax.find(f"{NS}scientific_name").text == "unclassified Flaviviridae"
        assert tax.find(f"{NS}rank") is None
        assert tax.find(f"{NS}id").text == "38144"
        assert tax.find(f"{NS}id").attrib["provider"] == "ncbi"

    def test_placeholder_values_not_written(self, tmp_path):
        meta = self._meta()
        meta["s1"]["host"] = "n/a"
        meta["s1"]["strain"] = "not provided"
        meta["s1"]["location"] = "missing: control sample"
        meta["s1"]["collection_date"] = "not collected"
        out = tmp_path / "out.xml"
        write_phyloxml(
            newick_path=tmp_path / "unused.nwk", output_xml=out,
            id_map={"s1": "a", "s2": "b", "s3": "c"}, leaf_metadata=meta,
            family="Flaviviridae", tree=_make_mini_tree(),
        )
        leaf = self._leaf(ET.parse(out).getroot(), "a")
        props = self._props(leaf)
        for ref in ("vipr:Host", "vipr:Strain", "vipr:Country", "vipr:Location",
                    "vipr:Collection_Date", "vipr:Year"):
            assert ref not in props
        assert leaf.find(f"{NS}date") is None
