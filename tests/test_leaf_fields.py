"""Tests for vfam_trees.leaf_fields."""
import pytest

from vfam_trees.leaf_fields import (
    clean_value, decimal_year, leaf_species, location_fields, own_taxon, ranked_name,
)


class TestCleanValue:
    @pytest.mark.parametrize("val", [
        None, "", "  ", "unknown", "N/A", "none", "missing", "Not Collected",
        "not provided", "not applicable", "restricted access",
        "missing: control sample", "Missing: lab stock",
    ])
    def test_placeholders(self, val):
        assert clean_value(val) == ""

    def test_real_values_kept(self):
        assert clean_value("  Homo sapiens ") == "Homo sapiens"
        assert clean_value("Missingham") == "Missingham"


class TestDecimalYear:
    """Midpoint of the stated interval, as Archaeopteryx's TipDateExtractor."""

    def test_year_only(self):
        assert decimal_year("2014") == pytest.approx(2014.5)

    def test_month_year_is_mid_month(self):
        # March 2014 starts after 59 days and lasts 31.
        assert decimal_year("Mar-2014") == pytest.approx(2014 + (59 + 15.5) / 365)
        assert decimal_year("2014-03") == decimal_year("Mar-2014")

    def test_day_is_mid_day(self):
        # 12 March 2014 is day 71.
        assert decimal_year("12-Mar-2014") == pytest.approx(2014 + 70.5 / 365)
        assert decimal_year("2014-03-12") == decimal_year("12-Mar-2014")
        assert decimal_year("2014-03-12T10:30:00Z") == decimal_year("2014-03-12")

    def test_leap_year(self):
        assert decimal_year("2016-02-29") == pytest.approx(2016 + 59.5 / 366)

    def test_case_insensitive_month(self):
        assert decimal_year("mar-2014") == decimal_year("Mar-2014")

    def test_range_is_midpoint(self):
        assert decimal_year("2014/2015") == pytest.approx(2015.0)
        assert decimal_year("2014-03-01/2014-03-31") == pytest.approx(decimal_year("Mar-2014"))

    @pytest.mark.parametrize("bad", [
        "", "unknown", "missing", "2014-13", "31-Feb-2014", "Foo-2014",
        "2014/2013", "1/2/3", "circa 2014", "14", None,
    ])
    def test_unparsable(self, bad):
        assert decimal_year(bad) is None


class TestLocationFields:
    def test_country_with_detail(self):
        assert location_fields({"location": "USA: Texas, Austin"}) == ("USA", "USA: Texas, Austin")

    def test_country_only(self):
        assert location_fields({"location": "Guinea"}) == ("Guinea", "")

    def test_placeholder(self):
        assert location_fields({"location": "unknown"}) == ("", "")
        assert location_fields({"location": "missing: control sample"}) == ("", "")
        assert location_fields({"location": "not collected"}) == ("", "")
        assert location_fields({}) == ("", "")


_LINEAGE = [
    {"name": "Filoviridae", "rank": "family", "taxid": "11266"},
    {"name": "Orthoebolavirus", "rank": "genus", "taxid": "3044781"},
    {"name": "Orthoebolavirus zairense", "rank": "species", "taxid": "3052462"},
    {"name": "Zaire ebolavirus", "rank": "no rank", "taxid": "186538"},
]


class TestTaxa:
    def test_ranked_name(self):
        meta = {"lineage_ranked": _LINEAGE}
        assert ranked_name(meta, "genus") == "Orthoebolavirus"
        assert ranked_name(meta, "subfamily") == ""

    def test_leaf_species_from_lineage(self):
        meta = {"species": "Zaire ebolavirus", "lineage_ranked": _LINEAGE}
        assert leaf_species(meta) == "Orthoebolavirus zairense"

    def test_leaf_species_absent_without_species_rank(self):
        meta = {"species": "Sapovirus sp.", "lineage_ranked": _LINEAGE[:2]}
        assert leaf_species(meta) == ""

    def test_leaf_species_falls_back_to_organism(self):
        assert leaf_species({"species": "my virus", "lineage_ranked": []}) == "my virus"
        assert leaf_species({"species": "unknown"}) == ""

    def test_own_taxon_requires_matching_taxid(self):
        assert own_taxon({"taxon_id": "186538", "lineage_ranked": _LINEAGE})["name"] == "Zaire ebolavirus"
        assert own_taxon({"taxon_id": "999", "lineage_ranked": _LINEAGE}) == {}
        # entries without taxids (lineages cached before taxids were stored)
        old = [{"name": e["name"], "rank": e["rank"]} for e in _LINEAGE]
        assert own_taxon({"taxon_id": "186538", "lineage_ranked": old}) == {}
