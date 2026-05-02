"""Tests for vfam_trees.fetch — query building and concat-mode helpers
(no network calls)."""
import pytest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from vfam_trees.fetch import (
    _build_marker_query,
    _build_species_query,
    _count_split_submissions,
    _extract_isolate,
    _safe_marker_filename,
    _source_nuc_accession,
    group_proteins_by_genome,
)


def test_whole_genome_nucleotide():
    q = _build_species_query(12345, "nucleotide", "whole_genome")
    assert "complete genome" in q or "complete sequence" in q
    assert "[Gene]" not in q


def test_whole_genome_protein():
    q = _build_species_query(12345, "protein", "whole_genome")
    assert "complete genome" not in q
    assert "[Gene]" not in q


def test_segment_query():
    q = _build_species_query(12345, "nucleotide", "whole_genome", segment="RNA1")
    assert '"RNA1"[Title]' in q
    assert "complete sequence" in q


def test_marker_gene_query():
    q = _build_species_query(12345, "nucleotide", "B646L")
    assert '"B646L"[Gene]' in q


def test_marker_gene_excludes_complete_genome():
    q = _build_species_query(12345, "nucleotide", "B646L")
    assert 'NOT "complete genome"[Title]' in q
    assert 'NOT "complete sequence"[Title]' in q


def test_marker_gene_hexon_excludes_complete():
    q = _build_species_query(12345, "nucleotide", "hexon")
    assert '"hexon"[Gene]' in q
    assert 'NOT "complete genome"[Title]' in q


def test_refseq_only_adds_filter():
    q = _build_species_query(12345, "nucleotide", "whole_genome", refseq_only=True)
    assert "refseq[filter]" in q


def test_patent_always_excluded():
    for region in ("whole_genome", "B646L"):
        q = _build_species_query(12345, "nucleotide", region)
        assert "patent[filter]" in q


def test_exclude_organisms():
    q = _build_species_query(12345, "nucleotide", "whole_genome",
                             exclude_organisms=["metagenome", "synthetic construct"])
    assert '"metagenome"[Organism]' in q
    assert '"synthetic construct"[Organism]' in q


def test_taxid_in_query():
    q = _build_species_query(99999, "nucleotide", "whole_genome")
    assert "txid99999" in q


def test_protein_marker_uses_protein_name_field():
    q = _build_species_query(12345, "protein", "DNA polymerase")
    assert '"DNA polymerase"[Protein Name]' in q


def test_protein_marker_also_includes_gene_fallback():
    q = _build_species_query(12345, "protein", "B646L")
    assert '"B646L"[Protein Name]' in q
    assert '"B646L"[Gene]' in q


def test_protein_marker_does_not_use_nucleotide_gene_only():
    q = _build_species_query(12345, "protein", "DNA polymerase")
    # [Gene] alone (without [Protein Name]) would miss most protein records
    assert '"DNA polymerase"[Protein Name]' in q


# ---------------------------------------------------------------------------
# Concat mode — helpers and group_proteins_by_genome
# ---------------------------------------------------------------------------


def _protein(acc: str, description: str = "DNA polymerase",
             length: int = 500, source_acc: str | None = None,
             isolate: str | None = None,
             db_source: str | None = None,
             locus_tag: str | None = None) -> SeqRecord:
    """Build a minimal protein-like SeqRecord for concat-mode helper tests."""
    rec = SeqRecord(Seq("M" * length), id=acc, description=description)
    rec.annotations = {"organism": "Test virus"}
    if db_source:
        rec.annotations["db_source"] = db_source
    if locus_tag:
        rec.annotations["locus_tag"] = locus_tag
    rec.features = []
    if source_acc:
        cds = SeqFeature(FeatureLocation(0, length), type="CDS")
        cds.qualifiers = {"coded_by": [f"{source_acc}:1..{length * 3}"]}
        rec.features.append(cds)
    if isolate:
        src = SeqFeature(FeatureLocation(0, length), type="source")
        src.qualifiers = {"isolate": [isolate]}
        rec.features.append(src)
    return rec


# ---- _safe_marker_filename ----

class TestSafeMarkerFilename:
    def test_basic(self):
        assert _safe_marker_filename("DNA polymerase") == "DNA_polymerase"

    def test_special_characters_replaced(self):
        assert _safe_marker_filename("DNA-directed RNA polymerase 147 kDa") == \
            "DNA-directed_RNA_polymerase_147_kDa"

    def test_empty_falls_back(self):
        assert _safe_marker_filename("") == "marker"


# ---- _source_nuc_accession ----

class TestSourceNucAccession:
    def test_coded_by_simple(self):
        rec = _protein("YP_1", source_acc="NC_001234.1")
        assert _source_nuc_accession(rec) == "NC_001234.1"

    def test_coded_by_complement(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        cds = SeqFeature(FeatureLocation(0, 100), type="CDS")
        cds.qualifiers = {"coded_by": ["complement(NC_001234.1:1..300)"]}
        rec.features = [cds]
        rec.annotations = {}
        assert _source_nuc_accession(rec) == "NC_001234.1"

    def test_db_source_fallback(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        rec.annotations = {"db_source": "REFSEQ: accession NC_999999.1"}
        rec.features = []
        assert _source_nuc_accession(rec) == "NC_999999.1"

    def test_unknown_returns_empty(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        rec.annotations = {}
        rec.features = []
        assert _source_nuc_accession(rec) == ""

    def test_coded_by_preferred_over_db_source(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        rec.annotations = {"db_source": "REFSEQ: accession NC_999999.1"}
        cds = SeqFeature(FeatureLocation(0, 100), type="CDS")
        cds.qualifiers = {"coded_by": ["NC_001234.1:1..300"]}
        rec.features = [cds]
        assert _source_nuc_accession(rec) == "NC_001234.1"


# ---- _extract_isolate ----

class TestExtractIsolate:
    def test_extracts_from_source_feature(self):
        rec = _protein("YP_1", isolate="ISO-42")
        assert _extract_isolate(rec) == "ISO-42"

    def test_falls_back_to_strain(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        src = SeqFeature(FeatureLocation(0, 100), type="source")
        src.qualifiers = {"strain": ["X-strain"]}
        rec.features = [src]
        rec.annotations = {}
        assert _extract_isolate(rec) == "X-strain"

    def test_no_source_feature_returns_empty(self):
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        rec.features = []
        rec.annotations = {}
        assert _extract_isolate(rec) == ""


# ---- _count_split_submissions ----

class TestCountSplitSubmissions:
    def test_no_shared_isolates_returns_zero(self):
        bucket = {
            "NC_001": {"DNA polymerase": [_protein("YP_1", source_acc="NC_001", isolate="A")]},
            "NC_002": {"MCP":            [_protein("YP_2", source_acc="NC_002", isolate="B")]},
        }
        assert _count_split_submissions(bucket) == 0

    def test_isolate_in_two_accessions_counts(self):
        # Same isolate "X" submitted as two separate GenBank records.
        bucket = {
            "NC_001": {"DNA polymerase": [_protein("YP_1", source_acc="NC_001", isolate="X")]},
            "NC_002": {"MCP":            [_protein("YP_2", source_acc="NC_002", isolate="X")]},
        }
        assert _count_split_submissions(bucket) == 1

    def test_isolate_only_in_one_accession_doesnt_count(self):
        # Same isolate appears multiple times within one accession — that's fine.
        bucket = {
            "NC_001": {
                "DNA polymerase": [_protein("YP_1", source_acc="NC_001", isolate="X")],
                "MCP":            [_protein("YP_2", source_acc="NC_001", isolate="X")],
            },
        }
        assert _count_split_submissions(bucket) == 0

    def test_records_without_isolate_are_skipped(self):
        bucket = {
            "NC_001": {"DNA polymerase": [_protein("YP_1", source_acc="NC_001")]},  # no isolate
            "NC_002": {"MCP":            [_protein("YP_2", source_acc="NC_002")]},
        }
        assert _count_split_submissions(bucket) == 0


# ---- group_proteins_by_genome ----

_POLB = {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase"]}
_MCP  = {"name": "major capsid protein", "aliases": ["MCP"]}
_HEL  = {"name": "DNA helicase", "aliases": []}


class TestGroupProteinsByGenome:
    def test_basic_grouping_keeps_complete_genome(self):
        proteins_by_marker = {
            "DNA polymerase":         [_protein("YP_1", "DNA polymerase", source_acc="NC_001")],
            "major capsid protein":   [_protein("YP_2", "major capsid protein", source_acc="NC_001")],
            "DNA helicase":           [_protein("YP_3", "DNA helicase", source_acc="NC_001")],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
        )
        assert "NC_001" in genomes
        assert set(genomes["NC_001"].keys()) == {
            "DNA polymerase", "major capsid protein", "DNA helicase"
        }
        assert stats["n_genomes_kept"] == 1
        assert stats["n_dropped_min_fraction"] == 0

    def test_two_genomes_separated_by_source_acc(self):
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "DNA polymerase", source_acc="NC_001"),
                _protein("YP_2", "DNA polymerase", source_acc="NC_002"),
            ],
            "major capsid protein": [
                _protein("YP_3", "major capsid protein", source_acc="NC_001"),
                _protein("YP_4", "major capsid protein", source_acc="NC_002"),
            ],
            "DNA helicase": [
                _protein("YP_5", "DNA helicase", source_acc="NC_001"),
                _protein("YP_6", "DNA helicase", source_acc="NC_002"),
            ],
        }
        genomes, _ = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
        )
        assert set(genomes.keys()) == {"NC_001", "NC_002"}
        for g in genomes.values():
            assert len(g) == 3

    def test_drops_low_coverage_genomes(self):
        # 3-marker set with min_fraction=0.7 → require ceil(2.1) = 2 markers.
        # NC_002 has only 1 marker → dropped.  NC_001 has 3 → kept.
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "DNA polymerase", source_acc="NC_001"),
                _protein("YP_4", "DNA polymerase", source_acc="NC_002"),
            ],
            "major capsid protein": [
                _protein("YP_2", "major capsid protein", source_acc="NC_001"),
            ],
            "DNA helicase": [
                _protein("YP_3", "DNA helicase", source_acc="NC_001"),
            ],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
        )
        assert set(genomes.keys()) == {"NC_001"}
        assert stats["n_genomes_kept"] == 1
        assert stats["n_dropped_min_fraction"] == 1
        assert stats["n_required_markers"] == 2

    def test_orphaned_proteins_counted(self):
        # Protein with no source_acc is excluded from grouping.
        rec = SeqRecord(Seq("M" * 100), id="YP_orphan", description="DNA polymerase")
        rec.annotations = {}
        rec.features = []
        proteins_by_marker = {
            "DNA polymerase": [rec],
            "major capsid protein": [],
            "DNA helicase": [],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
        )
        assert genomes == {}
        assert stats["n_orphaned_no_source"] == 1

    def test_paralog_tiebreaker_applied(self):
        # Two candidates for DNA polymerase in the same genome — tiebreaker picks
        # the one matching locus_tag_hint.
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_a", "DNA polymerase", source_acc="NC_001", locus_tag="A24R"),
                _protein("YP_b", "DNA polymerase", source_acc="NC_001", locus_tag="E9L"),
            ],
        }
        genomes, _ = group_proteins_by_genome(
            proteins_by_marker, [marker], None, min_fraction=1.0,
        )
        assert genomes["NC_001"]["DNA polymerase"].id == "YP_b"

    def test_subfamily_aware_aliases_passed_through(self):
        marker = {
            "name": "DNA polymerase",
            "aliases": [],
            "aliases_Entomopoxvirinae": ["polB-spec"],
        }
        # Only matches via the Entomopoxvirinae alias
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "polB-spec [some virus]", source_acc="NC_001"),
            ],
        }
        lineage = [{"rank": "subfamily", "name": "Entomopoxvirinae"}]
        genomes, _ = group_proteins_by_genome(
            proteins_by_marker, [marker], lineage, min_fraction=1.0,
        )
        assert "NC_001" in genomes

    def test_split_submission_detected(self):
        # NC_001 and NC_002 share isolate "ISO-42" but each carries only 1 marker.
        # Both fail min_fraction=0.7 (need 2 of 3 markers) and get dropped, but
        # the diagnostic counter detects the shared isolate.
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "DNA polymerase", source_acc="NC_001", isolate="ISO-42"),
            ],
            "major capsid protein": [
                _protein("YP_2", "major capsid protein", source_acc="NC_002", isolate="ISO-42"),
            ],
            "DNA helicase": [],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
        )
        assert stats["n_genomes_kept"] == 0
        assert stats["n_dropped_min_fraction"] == 2
        assert stats["n_dropped_split_submission"] == 1


# ---- _build_marker_query ----

class TestBuildMarkerQuery:
    def test_basic_query_uses_protein_name(self):
        q = _build_marker_query(
            taxid=12345,
            marker={"name": "DNA polymerase", "aliases": []},
            species_lineage=None,
            exclude_organisms=None,
        )
        assert "txid12345[Organism:exp]" in q
        assert '"DNA polymerase"[Protein Name]' in q

    def test_aliases_or_combined(self):
        q = _build_marker_query(
            taxid=12345,
            marker={"name": "DNA polymerase", "aliases": ["DNA pol", "pol B"]},
            species_lineage=None,
            exclude_organisms=None,
        )
        assert '"DNA polymerase"[Protein Name]' in q
        assert '"DNA pol"[Protein Name]' in q
        assert '"pol B"[Protein Name]' in q

    def test_subfamily_aliases_unioned_when_lineage_matches(self):
        marker = {
            "name": "DNA polymerase",
            "aliases": [],
            "aliases_Entomopoxvirinae": ["DNA polymerase B"],
        }
        lineage = [{"rank": "subfamily", "name": "Entomopoxvirinae"}]
        q = _build_marker_query(12345, marker, lineage, None)
        assert '"DNA polymerase B"[Protein Name]' in q

    def test_subfamily_aliases_skipped_when_lineage_other(self):
        marker = {
            "name": "DNA polymerase",
            "aliases": [],
            "aliases_Entomopoxvirinae": ["DNA polymerase B"],
        }
        lineage = [{"rank": "subfamily", "name": "Chordopoxvirinae"}]
        q = _build_marker_query(12345, marker, lineage, None)
        assert '"DNA polymerase B"[Protein Name]' not in q

    def test_refseq_filter(self):
        q = _build_marker_query(
            12345, {"name": "X", "aliases": []}, None, None, refseq_only=True,
        )
        assert "refseq[filter]" in q

    def test_excludes_terms(self):
        q = _build_marker_query(
            12345, {"name": "X", "aliases": []}, None, ["synthetic construct"],
        )
        assert 'NOT "synthetic construct"[Organism]' in q

    def test_patent_always_excluded(self):
        q = _build_marker_query(12345, {"name": "X", "aliases": []}, None, None)
        assert "NOT patent[filter]" in q
