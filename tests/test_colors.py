"""Tests for vfam_trees.colors — genus inference and leaf color assignment."""
import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees.colors import _infer_genus_from_lineage, _SPECIES_OR_BELOW, assign_leaf_colors


# ---------------------------------------------------------------------------
# _infer_genus_from_lineage
# ---------------------------------------------------------------------------

def _entry(name, rank):
    return {"name": name, "rank": rank}


def test_none_mode_returns_empty():
    lin = [_entry("Flavivirus", "no rank")]
    assert _infer_genus_from_lineage(lin, "none") == ""


def test_suffix_single_word_virus_matches():
    lin = [
        _entry("Flaviviridae", "family"),
        _entry("Flavivirus", "no rank"),
        _entry("Dengue virus", "species"),
    ]
    assert _infer_genus_from_lineage(lin, "suffix") == "Flavivirus"


def test_suffix_multiword_virus_not_matched():
    lin = [
        _entry("Flaviviridae", "family"),
        _entry("Dengue virus", "no rank"),
    ]
    assert _infer_genus_from_lineage(lin, "suffix") == ""


def test_suffix_already_genus_rank_still_matched():
    # If NCBI correctly marks it genus, suffix still finds it (safety net).
    lin = [_entry("Alphacoronavirus", "genus")]
    assert _infer_genus_from_lineage(lin, "suffix") == "Alphacoronavirus"


def test_suffix_no_match_returns_empty():
    lin = [
        _entry("Flaviviridae", "family"),
        _entry("Orthoflavi", "no rank"),
    ]
    assert _infer_genus_from_lineage(lin, "suffix") == ""


def test_deepest_falls_back_to_subfamily_when_no_suffix():
    lin = [
        _entry("Coronaviridae", "family"),
        _entry("Orthocoronavirinae", "subfamily"),
        _entry("Betacoronavirus SARS", "species"),
    ]
    result = _infer_genus_from_lineage(lin, "deepest")
    # Should pick subfamily (deepest entry above species)
    assert result == "Orthocoronavirinae"


def test_deepest_prefers_suffix_over_deeper_non_suffix():
    lin = [
        _entry("Coronaviridae", "family"),
        _entry("Orthocoronavirinae", "subfamily"),
        _entry("Betacoronavirus", "no rank"),   # suffix match, even if not deepest rank
        _entry("SARS-CoV-2", "species"),
    ]
    result = _infer_genus_from_lineage(lin, "deepest")
    # suffix wins because it's tried first
    assert result == "Betacoronavirus"


def test_deepest_skips_species_or_below():
    lin = [
        _entry("Flaviviridae", "family"),
        _entry("Dengue virus", "species"),
        _entry("Dengue virus type 1", "subspecies"),
    ]
    result = _infer_genus_from_lineage(lin, "deepest")
    assert result == "Flaviviridae"


def test_empty_lineage_returns_empty_for_all_modes():
    for mode in ("none", "suffix", "deepest"):
        assert _infer_genus_from_lineage([], mode) == ""


# ---------------------------------------------------------------------------
# assign_leaf_colors — genus_inference integration
# ---------------------------------------------------------------------------

def _make_record(short_id):
    r = SeqRecord(Seq("ATGC"), id=short_id)
    return r


def _make_meta(genus="", subfamily="", lineage_ranked=None):
    return {
        "lineage_ranked": lineage_ranked or [],
        "genus": genus,
        "subfamily": subfamily,
    }


def test_genus_inference_none_leaves_unclassified_grey():
    rec = _make_record("seq1")
    meta = _make_meta(lineage_ranked=[_entry("Flavivirus", "no rank")])
    d2c, _, g2c, _ = assign_leaf_colors(
        [rec],
        {"seq1": meta},
        {"seq1": "Dengue_virus|strain|ACC1"},
        genus_inference="none",
    )
    # No genus → all return empty dicts (no groups → early return)
    assert d2c == {}


def test_genus_inference_suffix_colors_inferred_genus():
    rec = _make_record("seq1")
    meta = _make_meta(lineage_ranked=[
        _entry("Flaviviridae", "family"),
        _entry("Flavivirus", "no rank"),
    ])
    d2c, _, g2c, _ = assign_leaf_colors(
        [rec],
        {"seq1": meta},
        {"seq1": "Dengue_virus|strain|ACC1"},
        genus_inference="suffix",
    )
    assert "Flavivirus" in g2c
    assert "Dengue_virus|strain|ACC1" in d2c
    assert d2c["Dengue_virus|strain|ACC1"] != "#888888"


def test_genus_inference_deepest_uses_subfamily_fallback():
    rec = _make_record("seq1")
    meta = _make_meta(lineage_ranked=[
        _entry("Coronaviridae", "family"),
        _entry("Orthocoronavirinae", "subfamily"),
        _entry("SARS-CoV-2", "species"),
    ])
    d2c, _, g2c, _ = assign_leaf_colors(
        [rec],
        {"seq1": meta},
        {"seq1": "SARSCoV2|wuhan|MN908947"},
        genus_inference="deepest",
    )
    assert "Orthocoronavirinae" in g2c


def test_single_subfamily_genera_span_full_hue_wheel():
    import colorsys
    # Four genera, all within one subfamily (Coronaviridae-like scenario).
    genera = ["Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus", "Deltacoronavirus"]
    records = [_make_record(f"seq{i}") for i in range(4)]
    metas = {
        f"seq{i}": _make_meta(lineage_ranked=[
            _entry("Coronaviridae", "family"),
            _entry("Orthocoronavirinae", "subfamily"),
            _entry(g, "genus"),
        ])
        for i, g in enumerate(genera)
    }
    display = {f"seq{i}": f"{g}|strain|ACC{i}" for i, g in enumerate(genera)}
    _, _, g2c, _ = assign_leaf_colors(records, metas, display)
    assert len(g2c) == 4
    hues = []
    for hex_c in g2c.values():
        r, g, b = (int(hex_c.lstrip("#")[i:i+2], 16) / 255 for i in (0, 2, 4))
        h, _, _ = colorsys.rgb_to_hls(r, g, b)
        hues.append(h)
    # Hues should span at least 0.5 of the wheel (full spread, not a narrow band)
    hue_range = max(hues) - min(hues)
    assert hue_range >= 0.5


def test_two_subfamilies_still_use_banding():
    import colorsys
    # Two subfamilies — hue banding should apply (each gets its own hue zone).
    records = [_make_record(f"seq{i}") for i in range(2)]
    metas = {
        "seq0": _make_meta(lineage_ranked=[
            _entry("Flaviviridae", "family"),
            _entry("Flavivirinae", "subfamily"),
            _entry("Flavivirus", "genus"),
        ]),
        "seq1": _make_meta(lineage_ranked=[
            _entry("Flaviviridae", "family"),
            _entry("Pegivirinae", "subfamily"),
            _entry("Pegivirus", "genus"),
        ]),
    }
    display = {"seq0": "Dengue|s|A", "seq1": "Pegivirus|s|B"}
    _, _, g2c, s2g = assign_leaf_colors(records, metas, display)
    assert len(g2c) == 2
    assert len(s2g) == 2   # two distinct subfamily groups


def test_formal_genus_rank_takes_precedence_over_inference():
    rec = _make_record("seq1")
    meta = _make_meta(lineage_ranked=[
        _entry("Flaviviridae", "family"),
        _entry("Orthoflavivirus", "genus"),
        _entry("Flavivirus_like", "no rank"),
    ])
    d2c, _, g2c, _ = assign_leaf_colors(
        [rec],
        {"seq1": meta},
        {"seq1": "Dengue|strain|ACC"},
        genus_inference="suffix",
    )
    # Formal genus should be picked up in the main loop, not suffix inference
    assert "Orthoflavivirus" in g2c
    assert "Flavivirus_like" not in g2c
