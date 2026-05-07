"""Tests for vfam_trees.pipeline_concat — concat dispatch + skip path.

The full concat run path requires MMseqs2 / MAFFT / IQ-TREE / Entrez and is
exercised in integration when the pipeline runs.  These tests verify the
public interface and the skip-path with stubbed fetch — no external tools
required.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees import pipeline_concat
from vfam_trees.pipeline_concat import _build_concat_leaf_metadata


# ---------------------------------------------------------------------------
# Module surface
# ---------------------------------------------------------------------------

def test_run_family_concat_is_exported():
    assert hasattr(pipeline_concat, "run_family_concat")
    assert callable(pipeline_concat.run_family_concat)


# ---------------------------------------------------------------------------
# Skip path: no genomes pass min_fraction → cleanly emits skip row.
# ---------------------------------------------------------------------------

def _minimal_family_cfg() -> dict:
    return {
        "download": {"max_per_species": 100},
        "sequence": {"type": "protein", "region": "concatenated", "segment": None},
        "quality": {
            "min_length": None,
            "max_ambiguous": 0.01,
            "exclude_organisms": [],
        },
        "clustering": {
            "tool": "mmseqs2",
            "threshold_min": 0.7,
            "threshold_max": 0.99,
            "max_reps_500": 20,
            "max_reps_100": 5,
        },
        "targets":   {"max_500": 500, "max_100": 100},
        "msa_500":   {"tool": "mafft", "options_aa": "--6merpair --retree 1"},
        "msa_100":   {"tool": "mafft", "options_aa": "--auto"},
        "msa_trim":  {"enabled": True, "tool": "trimal", "options": "-automated1"},
        "tree_500":  {"tool": "fasttree", "options": "", "model_aa": "LG+G"},
        "tree_100":  {"tool": "iqtree", "options_aa": "-B 1000", "model_aa": "TEST"},
        "refseq_absorption": {"enabled": True, "threshold": 0.99},
        "length_outlier":    {"enabled": True, "hi_mult": 3.0, "lo_mult": 1.0 / 3.0},
        "outlier_removal":   {"enabled": True, "factor": 20.0,
                              "max_iterations": 3, "min_seqs": 40},
        "taxonomy":          {"lca_min_rank": "none"},
        "concatenation": {
            "proteins": [
                {"name": "DNA polymerase", "aliases": []},
                {"name": "MCP",            "aliases": []},
            ],
            "min_fraction": 0.7,
            "partition_tree_100": True,
            "partition_tree_500": False,
        },
    }


def test_skip_path_when_fetch_returns_no_genomes(tmp_path, monkeypatch):
    """When fetch_species_genomes returns no genomes for any species, the
    runner emits the 'too few genomes' skip and writes summary/status rows.
    """
    # Stub the heavy I/O — Entrez + per-marker fetch
    def fake_fetch_species_genomes(*args, **kwargs):
        return {}, {"n_proteins_fetched": 0, "n_genomes_found": 0,
                    "n_genomes_kept": 0, "n_dropped_min_fraction": 0,
                    "n_dropped_split_submission": 0, "n_orphaned_no_source": 0}
    def fake_fetch_taxonomy_lineages(taxids):
        return {}
    monkeypatch.setattr(pipeline_concat, "fetch_species_genomes", fake_fetch_species_genomes)
    monkeypatch.setattr(pipeline_concat, "fetch_taxonomy_lineages", fake_fetch_taxonomy_lineages)

    # Capture mark_skipped / mark_done invocations
    mark_skipped_calls: list[tuple] = []
    mark_done_calls: list[tuple] = []
    def fake_mark_skipped(family_dir, family, output_dir, reason):
        mark_skipped_calls.append((family, reason))
    def fake_mark_done(family_dir, family, output_dir):
        mark_done_calls.append((family,))

    family_dir = tmp_path / "Poxviridae_10240"
    family_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    summary_path = output_dir / "summary.tsv"
    status_path  = output_dir / "status.tsv"

    species_list = [
        {"name": "Test virus 1", "taxid": 100001},
        {"name": "Test virus 2", "taxid": 100002},
    ]

    pipeline_concat.run_family_concat(
        family="Poxviridae",
        family_cfg=_minimal_family_cfg(),
        family_taxid=10240,
        family_lineage=[{"rank": "family", "name": "Poxviridae"}],
        family_annotation={"baltimore_class": "I"},
        family_dir=family_dir,
        work_dir=work_dir,
        output_dir=output_dir,
        species_list=species_list,
        threads=1,
        summary_path=summary_path,
        status_path=status_path,
        mark_done=fake_mark_done,
        mark_skipped=fake_mark_skipped,
    )

    # Skip path must have run
    assert len(mark_skipped_calls) == 1
    fam, reason = mark_skipped_calls[0]
    assert fam == "Poxviridae"
    assert "too few genomes" in reason
    assert mark_done_calls == []

    # Summary + status rows written
    assert summary_path.exists()
    assert status_path.exists()
    summary_text = summary_path.read_text()
    status_text = status_path.read_text()
    assert "Poxviridae" in summary_text
    assert "Poxviridae" in status_text
    assert "too few genomes" in status_text


# ---------------------------------------------------------------------------
# _build_concat_leaf_metadata — regression tests for PhyloXML metadata
# enrichment (host / collection_date / location / strain / taxon_id pulled
# from the protein record source feature so concat-mode XML carries the same
# vipr: properties as single-protein mode).
# ---------------------------------------------------------------------------

def _make_protein_record(
    accession: str,
    *,
    host: str | None = None,
    country: str | None = None,
    collection_date: str | None = None,
    strain: str | None = None,
    taxon_id: str | None = None,
    organism: str = "Test virus",
) -> SeqRecord:
    """Build a minimal protein SeqRecord with a populated source feature."""
    qualifiers: dict[str, list[str]] = {}
    if host:
        qualifiers["host"] = [host]
    if country:
        qualifiers["country"] = [country]
    if collection_date:
        qualifiers["collection_date"] = [collection_date]
    if strain:
        qualifiers["strain"] = [strain]
    if taxon_id:
        qualifiers["db_xref"] = [f"taxon:{taxon_id}"]

    rec = SeqRecord(Seq("MKQ"), id=accession, description=f"protein from {accession}")
    rec.annotations = {"organism": organism, "taxonomy": []}
    rec.features = [
        SeqFeature(FeatureLocation(0, 3), type="source", qualifiers=qualifiers),
    ]
    return rec


class TestBuildConcatLeafMetadata:
    def test_basic_fields_always_populated(self):
        rec = _make_protein_record("YP_001")
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_001": {"DNA polymerase": rec}},
            genome_to_species={"NC_001": "Test virus"},
            short_to_display={"NC_001": "Test_virus|NC_001"},
            species_lineages={"Test virus": [{"name": "Genusname", "rank": "genus"}]},
        )
        assert "NC_001" in meta
        assert meta["NC_001"]["species"] == "Test virus"
        assert meta["NC_001"]["accession"] == "NC_001"
        assert meta["NC_001"]["seq_name"] == "Test_virus|NC_001"
        assert meta["NC_001"]["lineage_ranked"] == [{"name": "Genusname", "rank": "genus"}]

    def test_source_feature_fields_extracted(self):
        rec = _make_protein_record(
            "YP_002",
            host="Homo sapiens",
            country="USA",
            collection_date="2020-05-14",
            strain="DEN1",
            taxon_id="11053",
        )
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_002": {"DNA polymerase": rec}},
            genome_to_species={"NC_002": "Dengue virus"},
            short_to_display={"NC_002": "Dengue|NC_002"},
            species_lineages={},
        )
        m = meta["NC_002"]
        assert m["host"] == "Homo sapiens"
        assert m["location"] == "USA"
        assert m["collection_date"] == "2020-05-14"
        assert m["strain"] == "DEN1"
        assert m["taxon_id"] == "11053"

    def test_unknown_placeholder_skipped(self):
        # No source-feature qualifiers → extract_metadata returns "unknown"
        # placeholders; the helper must drop them so absent fields stay absent.
        rec = _make_protein_record("YP_003")
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_003": {"DNA polymerase": rec}},
            genome_to_species={"NC_003": "Test virus"},
            short_to_display={"NC_003": "Test|NC_003"},
            species_lineages={},
        )
        m = meta["NC_003"]
        for absent in ("host", "collection_date", "location", "strain"):
            assert absent not in m

    def test_partial_source_feature_only_keeps_present_fields(self):
        rec = _make_protein_record("YP_004", host="Homo sapiens", country="Japan")
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_004": {"DNA polymerase": rec}},
            genome_to_species={"NC_004": "Test virus"},
            short_to_display={"NC_004": "Test|NC_004"},
            species_lineages={},
        )
        m = meta["NC_004"]
        assert m["host"] == "Homo sapiens"
        assert m["location"] == "Japan"
        assert "collection_date" not in m
        assert "strain" not in m

    def test_uses_first_marker_record_when_multiple(self):
        # Two markers per genome — the helper picks one (any), source feature is
        # the same on both so the result is identical regardless of order.
        rec_a = _make_protein_record("YP_a", host="Homo sapiens")
        rec_b = _make_protein_record("YP_b", host="Homo sapiens")
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_005": {"DNA polymerase": rec_a, "MCP": rec_b}},
            genome_to_species={"NC_005": "Test virus"},
            short_to_display={"NC_005": "Test|NC_005"},
            species_lineages={},
        )
        assert meta["NC_005"]["host"] == "Homo sapiens"

    def test_empty_marker_records_skips_extraction(self):
        # Defensive: a genome with no protein records (shouldn't happen in
        # practice but the helper must not crash).
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_006": {}},
            genome_to_species={"NC_006": "Test virus"},
            short_to_display={"NC_006": "Test|NC_006"},
            species_lineages={},
        )
        m = meta["NC_006"]
        assert m["accession"] == "NC_006"
        assert "host" not in m
        assert "location" not in m

    def test_keyed_by_short_id_not_display_name(self):
        # Critical wiring: keys must be short_id (gid) since the PhyloXML
        # writer looks up leaf_metadata via clade.name (= short_id).
        rec = _make_protein_record("YP_007", host="Homo sapiens")
        meta = _build_concat_leaf_metadata(
            selected_genomes={"NC_007": {"DNA polymerase": rec}},
            genome_to_species={"NC_007": "Test virus"},
            short_to_display={"NC_007": "Test_virus|NC_007"},
            species_lineages={},
        )
        assert "NC_007" in meta
        assert "Test_virus|NC_007" not in meta

    def test_multiple_genomes_each_get_own_metadata(self):
        rec1 = _make_protein_record("YP_x", host="Homo sapiens", strain="A")
        rec2 = _make_protein_record("YP_y", host="Mus musculus", strain="B")
        meta = _build_concat_leaf_metadata(
            selected_genomes={
                "NC_x": {"DNA polymerase": rec1},
                "NC_y": {"DNA polymerase": rec2},
            },
            genome_to_species={"NC_x": "Virus X", "NC_y": "Virus Y"},
            short_to_display={"NC_x": "X|NC_x", "NC_y": "Y|NC_y"},
            species_lineages={},
        )
        assert meta["NC_x"]["host"] == "Homo sapiens"
        assert meta["NC_x"]["strain"] == "A"
        assert meta["NC_y"]["host"] == "Mus musculus"
        assert meta["NC_y"]["strain"] == "B"


# ---------------------------------------------------------------------------
# PhyloXML output — end-to-end regression: exercise _build_concat_leaf_metadata
# through write_phyloxml so we verify both the metadata enrichment AND that
# leaf_colors is keyed by short_id (not display name) for the font_color
# property to fire.
# ---------------------------------------------------------------------------

class TestConcatPhyloxmlIntegration:
    def _write_xml(self, tmp_path):
        from Bio.Phylo.BaseTree import Clade, Tree
        from vfam_trees.phyloxml_writer import write_phyloxml

        rec = _make_protein_record(
            "YP_001",
            host="Homo sapiens",
            country="USA",
            collection_date="2020-05-14",
            strain="DEN1",
            taxon_id="11053",
        )
        selected_genomes = {"NC_001": {"DNA polymerase": rec}}
        short_to_display = {"NC_001": "Dengue_virus|NC_001"}
        leaf_metadata = _build_concat_leaf_metadata(
            selected_genomes=selected_genomes,
            genome_to_species={"NC_001": "Dengue virus"},
            short_to_display=short_to_display,
            species_lineages={"Dengue virus": [{"name": "Flavivirus", "rank": "genus"}]},
        )
        # Match the wiring used in pipeline_concat._run_target_concat: the
        # leaf_colors dict is keyed by short_id (gid), NOT display name.
        short_to_color = {"NC_001": "#ff0000"}

        leaf = Clade(branch_length=0.1, name="NC_001")
        root = Clade(clades=[leaf])
        tree = Tree(root=root, rooted=True)

        out = tmp_path / "concat.xml"
        write_phyloxml(
            newick_path=tmp_path / "_.nwk",
            output_xml=out,
            id_map=short_to_display,
            leaf_metadata=leaf_metadata,
            family="Flaviviridae",
            tree=tree,
            phylogeny_name="ConcatTest",
            phylogeny_detail="concatenated 2-marker phylogeny [DNA polymerase, MCP] (target_500)",
            leaf_colors=short_to_color,
        )
        return out

    def test_vipr_host_property_present(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="vipr:Host"' in xml
        assert "Homo sapiens" in xml

    def test_vipr_collection_date_property_present(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="vipr:Collection_Date"' in xml
        assert "2020-05-14" in xml

    def test_vipr_location_property_present(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="vipr:Location"' in xml
        assert "USA" in xml

    def test_vipr_strain_property_present(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="vipr:Strain"' in xml
        assert "DEN1" in xml

    def test_vipr_year_property_present(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="vipr:Year"' in xml
        assert ">2020<" in xml

    def test_style_font_color_property_present(self, tmp_path):
        # Critical regression: leaf_colors must be keyed by short_id (gid),
        # otherwise the writer's lookup (by clade.name = short_id) misses
        # and the font_color property is silently omitted.
        xml = self._write_xml(tmp_path).read_text()
        assert 'ref="style:font_color"' in xml
        assert "#ff0000" in xml

    def test_taxon_id_emitted_as_taxonomy_id(self, tmp_path):
        xml = self._write_xml(tmp_path).read_text()
        assert "11053" in xml
