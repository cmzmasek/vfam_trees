"""Tests for the Auspice v2 JSON exporter (vfam_trees.auspice_writer)."""
import json
from io import StringIO

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from vfam_trees.auspice_writer import (
    write_auspice_json,
    _extract_year,
    _clean,
    _leaf_traits,
    _format_support,
)
from vfam_trees.config import DEFAULT_FAMILY_CONFIG


# --------------------------------------------------------------------------
# Fixtures / helpers
# --------------------------------------------------------------------------
def _build_tree():
    """((A:0.1,B:0.2)Ebolavirus:0.05,C:0.3) with support on the internal node."""
    a = Clade(branch_length=0.1, name="A")
    b = Clade(branch_length=0.2, name="B")
    internal = Clade(branch_length=0.05, name="Ebolavirus", clades=[a, b])
    internal.confidence = 95.0
    c = Clade(branch_length=0.3, name="C")
    root = Clade(clades=[internal, c])
    return Tree(root=root)


_META = {
    "A": {
        "species": "Zaire ebolavirus",
        "accession": "NC_002549.1",
        "host": "Homo sapiens",
        "location": "Guinea",
        "collection_date": "2014-05-01",
        "lineage_ranked": [
            {"rank": "subfamily", "name": "Orthovirinae"},
            {"rank": "genus", "name": "Ebolavirus"},
        ],
    },
    "B": {
        "species": "Sudan ebolavirus",
        "accession": "NC_006432.1",
        "host": "Homo sapiens",
        "collection_date": "2000",
        "lineage_ranked": [
            {"rank": "subfamily", "name": "Orthovirinae"},
            {"rank": "genus", "name": "Ebolavirus"},
        ],
    },
    "C": {
        # placeholder species + no lineage → minimal traits
        "species": "unknown",
        "accession": "NC_004161.1",
        "host": "",
        "lineage_ranked": [],
    },
}

_ID_MAP = {
    "A": "Zaire_ebolavirus|NC_002549.1|human",
    "B": "Sudan_ebolavirus|NC_006432.1|human",
    "C": "NC_004161.1",
}

_GENUS_COLORS = {"Ebolavirus": "#ff0000"}


def _load(tmp_path):
    out = tmp_path / "Filoviridae_tree_100_auspice.json"
    write_auspice_json(
        output_json=out,
        id_map=_ID_MAP,
        leaf_metadata=_META,
        family="Filoviridae",
        tree=_build_tree(),
        title="Filoviridae (test)",
        genus_to_color=_GENUS_COLORS,
        confidence_type="SH-aLRT",
        seq_type="nucleotide",
    )
    with open(out) as f:
        return json.load(f)


def _all_nodes(node):
    yield node
    for c in node.get("children", []):
        yield from _all_nodes(c)


def _leaves(node):
    return [n for n in _all_nodes(node) if "children" not in n]


# --------------------------------------------------------------------------
# Top-level document
# --------------------------------------------------------------------------
def test_version_is_v2(tmp_path):
    assert _load(tmp_path)["version"] == "v2"


def test_meta_has_required_keys(tmp_path):
    meta = _load(tmp_path)["meta"]
    for key in ("title", "updated", "panels", "colorings", "display_defaults"):
        assert key in meta
    assert meta["panels"] == ["tree"]
    assert meta["title"] == "Filoviridae (test)"


def test_valid_json_round_trips(tmp_path):
    # _load already json.load()s it; a second dump/parse confirms it's clean.
    doc = _load(tmp_path)
    json.loads(json.dumps(doc))


# --------------------------------------------------------------------------
# Tree structure & divergence
# --------------------------------------------------------------------------
def test_leaf_count_matches_input(tmp_path):
    assert len(_leaves(_load(tmp_path)["tree"])) == 3


def test_root_div_is_zero(tmp_path):
    assert _load(tmp_path)["tree"]["node_attrs"]["div"] == 0.0


def test_div_is_cumulative(tmp_path):
    root = _load(tmp_path)["tree"]
    by_name = {n["name"]: n for n in _all_nodes(root)}
    # A sits below internal(0.05) then A(0.1) → 0.15
    assert by_name["A"]["node_attrs"]["div"] == pytest.approx(0.15)
    assert by_name["B"]["node_attrs"]["div"] == pytest.approx(0.25)
    assert by_name["C"]["node_attrs"]["div"] == pytest.approx(0.30)


def test_internal_nodes_get_synthetic_names(tmp_path):
    root = _load(tmp_path)["tree"]
    internals = [n for n in _all_nodes(root) if "children" in n]
    assert all(n["name"].startswith("NODE_") for n in internals)


def test_leaf_names_are_short_ids(tmp_path):
    names = {n["name"] for n in _leaves(_load(tmp_path)["tree"])}
    assert names == {"A", "B", "C"}


# --------------------------------------------------------------------------
# Branch attributes: support + clade labels
# --------------------------------------------------------------------------
def test_support_emitted_as_branch_label(tmp_path):
    root = _load(tmp_path)["tree"]
    labelled = [
        n for n in _all_nodes(root)
        if n.get("branch_attrs", {}).get("labels", {}).get("SH_aLRT")
    ]
    assert labelled
    assert labelled[0]["branch_attrs"]["labels"]["SH_aLRT"] == "95"


def test_clade_label_from_internal_name(tmp_path):
    root = _load(tmp_path)["tree"]
    clades = [
        n["branch_attrs"]["labels"]["clade"]
        for n in _all_nodes(root)
        if n.get("branch_attrs", {}).get("labels", {}).get("clade")
    ]
    assert "Ebolavirus" in clades


# --------------------------------------------------------------------------
# Node attributes / traits
# --------------------------------------------------------------------------
def test_leaf_traits_emitted_as_value_objects(tmp_path):
    a = {n["name"]: n for n in _leaves(_load(tmp_path)["tree"])}["A"]
    assert a["node_attrs"]["genus"] == {"value": "Ebolavirus"}
    assert a["node_attrs"]["species"] == {"value": "Zaire ebolavirus"}
    assert a["node_attrs"]["host"] == {"value": "Homo sapiens"}


def test_placeholder_species_skipped(tmp_path):
    c = {n["name"]: n for n in _leaves(_load(tmp_path)["tree"])}["C"]
    assert "species" not in c["node_attrs"]
    assert "host" not in c["node_attrs"]


def test_year_extracted_as_int(tmp_path):
    a = {n["name"]: n for n in _leaves(_load(tmp_path)["tree"])}["A"]
    assert a["node_attrs"]["year"] == {"value": 2014}


def test_display_name_attr_present(tmp_path):
    a = {n["name"]: n for n in _leaves(_load(tmp_path)["tree"])}["A"]
    assert a["node_attrs"]["display_name"]["value"] == "Zaire_ebolavirus|NC_002549.1|human"


# --------------------------------------------------------------------------
# Colorings / filters / defaults
# --------------------------------------------------------------------------
def test_genus_coloring_carries_scale(tmp_path):
    colorings = _load(tmp_path)["meta"]["colorings"]
    genus = next(c for c in colorings if c["key"] == "genus")
    assert genus["type"] == "categorical"
    assert genus["scale"] == [["Ebolavirus", "#ff0000"]]


def test_year_coloring_is_continuous(tmp_path):
    colorings = _load(tmp_path)["meta"]["colorings"]
    year = next(c for c in colorings if c["key"] == "year")
    assert year["type"] == "continuous"


def test_absent_trait_has_no_coloring(tmp_path):
    # No leaf carries a 'strain' value → no strain coloring.
    keys = {c["key"] for c in _load(tmp_path)["meta"]["colorings"]}
    assert "strain" not in keys


def test_display_defaults(tmp_path):
    dd = _load(tmp_path)["meta"]["display_defaults"]
    assert dd["color_by"] == "genus"
    assert dd["distance_measure"] == "div"
    assert dd["tip_label"] == "display_name"


def test_filters_restricted_to_present_traits(tmp_path):
    filters = _load(tmp_path)["meta"]["filters"]
    assert "genus" in filters
    assert "year" in filters
    assert "strain" not in filters


# --------------------------------------------------------------------------
# Parsing fallback + edge cases
# --------------------------------------------------------------------------
def test_parses_from_newick_path(tmp_path):
    nwk = tmp_path / "t.nwk"
    Phylo.write(_build_tree(), str(nwk), "newick")
    out = tmp_path / "out_auspice.json"
    write_auspice_json(
        output_json=out,
        id_map=_ID_MAP,
        leaf_metadata=_META,
        family="Filoviridae",
        newick_path=nwk,
    )
    doc = json.loads(out.read_text())
    assert len(_leaves(doc["tree"])) == 3


def test_requires_tree_or_newick(tmp_path):
    with pytest.raises(ValueError):
        write_auspice_json(
            output_json=tmp_path / "x.json",
            id_map={}, leaf_metadata={}, family="X",
        )


def test_no_genus_colors_omits_scale(tmp_path):
    out = tmp_path / "ns_auspice.json"
    write_auspice_json(
        output_json=out, id_map=_ID_MAP, leaf_metadata=_META,
        family="Filoviridae", tree=_build_tree(),
    )
    colorings = json.loads(out.read_text())["meta"]["colorings"]
    genus = next(c for c in colorings if c["key"] == "genus")
    assert "scale" not in genus


# --------------------------------------------------------------------------
# Helper units
# --------------------------------------------------------------------------
def test_extract_year():
    assert _extract_year("2014-05-01") == "2014"
    assert _extract_year("collected 1998 in DRC") == "1998"
    assert _extract_year("") == ""
    assert _extract_year("no date here") == ""


def test_clean_drops_placeholders():
    assert _clean("Homo sapiens") == "Homo sapiens"
    assert _clean("unknown") == ""
    assert _clean("  N/A ") == ""
    assert _clean(None) == ""


def test_format_support_normalizes_fraction():
    assert _format_support(0.95) == "95"
    assert _format_support(95) == "95"
    assert _format_support(98.5) == "98.5"


def test_leaf_traits_extracts_lineage_ranks():
    traits = _leaf_traits(_META["A"])
    assert traits["genus"] == "Ebolavirus"
    assert traits["subfamily"] == "Orthovirinae"
    assert "host" in traits


# --------------------------------------------------------------------------
# Config wiring
# --------------------------------------------------------------------------
def test_default_config_enables_auspice():
    assert DEFAULT_FAMILY_CONFIG["output"]["auspice_json"] is True
