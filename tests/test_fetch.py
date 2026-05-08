"""Tests for vfam_trees.fetch — query building and concat-mode helpers
(no network calls)."""
import pytest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from vfam_trees.fetch import (
    _LOCUS_RE,
    _build_marker_query,
    _build_species_query,
    _count_split_submissions,
    _extract_isolate,
    _safe_marker_filename,
    _source_nuc_accession,
    fetch_nuc_lengths,
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


# ---- _build_marker_query (concat-mode) ----

def test_concat_marker_query_includes_gene_fallback():
    """Regression for v1.2.12: concat-mode protein queries must search
    [Gene] as well as [Protein Name] so gene-symbol-named markers
    (B646L, UL30, ...) match records that only carry the symbol in [Gene]."""
    marker = {"name": "B646L", "aliases": ["p72", "major capsid protein p72"]}
    q = _build_marker_query(taxid=137992, marker=marker, species_lineage=None,
                            exclude_organisms=None)
    for n in ("B646L", "p72", "major capsid protein p72"):
        assert f'"{n}"[Protein Name]' in q
        assert f'"{n}"[Gene]' in q


def test_concat_marker_query_subfamily_aliases_get_gene_fallback():
    """Subfamily-specific aliases (e.g. aliases_Entomopoxvirinae) must also
    be expanded with the [Gene] fallback."""
    marker = {
        "name": "DNA polymerase",
        "aliases": ["DNA-directed DNA polymerase"],
        "aliases_Entomopoxvirinae": ["EsV-1-110"],
    }
    species_lineage = [
        {"rank": "subfamily", "name": "Entomopoxvirinae"},
        {"rank": "species",   "name": "Entomopoxvirus X"},
    ]
    q = _build_marker_query(taxid=12345, marker=marker,
                            species_lineage=species_lineage,
                            exclude_organisms=None)
    assert '"EsV-1-110"[Protein Name]' in q
    assert '"EsV-1-110"[Gene]' in q


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

    def test_short_digit_fragment_not_extracted(self):
        # Real-world: a free-text coded_by fragment like "Q8" must not be
        # picked up as a source-nuc accession — NCBI's esummary then rejects
        # the whole batch with "Invalid uid Q8".  Real nuccore accessions have
        # at least 4 digits, so the regex requires that.
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        cds = SeqFeature(FeatureLocation(0, 100), type="CDS")
        cds.qualifiers = {"coded_by": ["join(Q8, NC_001234.1:1..300)"]}
        rec.features = [cds]
        rec.annotations = {}
        assert _source_nuc_accession(rec) == "NC_001234.1"


# ---- fetch_nuc_lengths input validation ----

class TestFetchNucLengthsInputValidation:
    def test_malformed_uids_filtered_before_esummary(self, monkeypatch):
        from vfam_trees import fetch as fetch_mod

        captured_batches: list[list[str]] = []

        class _FakeHandle:
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False

        def fake_esummary(db, id):
            captured_batches.append(id.split(","))
            return _FakeHandle()

        def fake_read(handle):
            return []

        monkeypatch.setattr(fetch_mod.Entrez, "esummary", fake_esummary)
        monkeypatch.setattr(fetch_mod.Entrez, "read", fake_read)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        # Mix of well-formed and malformed accessions
        fetch_mod.fetch_nuc_lengths(
            ["NC_001234.1", "Q8", "AB123456", "", "X1", "NC_001234.1"]
        )
        assert captured_batches == [["NC_001234.1", "AB123456"]]


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

    def test_source_nuc_length_filter_drops_short_partial_submission(self):
        # Mimic the Asfarviridae pattern: one complete genome (NC_001, big
        # parent nuc) shares the family with a partial single-gene
        # submission (PG_001, tiny parent nuc) that happens to also hit the
        # marker query.  Without the filter, the partial submission appears
        # as a 1-marker "genome" and dilutes the dataset; with the filter
        # it gets dropped before min_fraction.
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "DNA polymerase", source_acc="NC_001"),
                _protein("YP_4", "DNA polymerase", source_acc="PG_001"),
            ],
            "major capsid protein": [
                _protein("YP_2", "major capsid protein", source_acc="NC_001"),
            ],
            "DNA helicase": [
                _protein("YP_3", "DNA helicase", source_acc="NC_001"),
            ],
        }
        # Complete genome 200 kb, partial submission 2 kb — partial < 30% of 200kb.
        lengths = {"NC_001": 200_000, "PG_001": 2_000}
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
            source_nuc_min_length_frac=0.3,
            nuc_length_lookup=lambda accs: {a: lengths[a] for a in accs if a in lengths},
        )
        assert set(genomes.keys()) == {"NC_001"}
        assert stats["n_dropped_short_source_nuc"] == 1
        assert stats["n_genomes_found"] == 1   # post-filter count
        assert stats["n_dropped_min_fraction"] == 0

    def test_source_nuc_length_filter_disabled_by_default(self):
        # frac=0 means no length lookup is invoked at all; behaviour matches
        # the legacy pure path.
        sentinel = {"called": False}
        def boom(_):
            sentinel["called"] = True
            raise AssertionError("should not be called when frac is 0")
        proteins_by_marker = {
            "DNA polymerase":         [_protein("YP_1", "DNA polymerase", source_acc="NC_001")],
            "major capsid protein":   [_protein("YP_2", "major capsid protein", source_acc="NC_001")],
            "DNA helicase":           [_protein("YP_3", "DNA helicase", source_acc="NC_001")],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB, _MCP, _HEL], None, min_fraction=0.7,
            source_nuc_min_length_frac=0.0,
            nuc_length_lookup=boom,
        )
        assert sentinel["called"] is False
        assert "NC_001" in genomes
        assert stats["n_dropped_short_source_nuc"] == 0

    def test_source_nuc_length_filter_keeps_all_when_lookup_returns_empty(self):
        # If the network call returns nothing (e.g. offline / Entrez hiccup),
        # the filter logs a warning and lets every bucket through unchanged.
        proteins_by_marker = {
            "DNA polymerase": [
                _protein("YP_1", "DNA polymerase", source_acc="NC_001"),
                _protein("YP_2", "DNA polymerase", source_acc="PG_001"),
            ],
        }
        genomes, stats = group_proteins_by_genome(
            proteins_by_marker, [_POLB], None, min_fraction=1.0,
            source_nuc_min_length_frac=0.3,
            nuc_length_lookup=lambda accs: {},  # empty = lookup failed
        )
        assert set(genomes.keys()) == {"NC_001", "PG_001"}
        assert stats["n_dropped_short_source_nuc"] == 0


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


# ---------------------------------------------------------------------------
# _LOCUS_RE — regex used by _fetch_batch to count returned GenBank records.
# Regression for v1.2.12: previous implementation used substring counting
# against a leading-newline "LOCUS " marker which missed the first record
# and could be inflated by "LOCUS " appearing inside REFERENCE / COMMENT blocks.
# ---------------------------------------------------------------------------

class TestLocusRegex:
    def test_counts_first_record(self):
        # First record has no leading newline — a "\\nLOCUS " substring
        # search would miss it; the anchored regex must catch it.
        data = "LOCUS       NC_001234           500 bp\n//\n"
        assert len(_LOCUS_RE.findall(data)) == 1

    def test_counts_multiple_records(self):
        data = (
            "LOCUS       NC_001234           500 bp\n//\n"
            "LOCUS       MN567890           600 bp\n//\n"
            "LOCUS       KX111111           700 bp\n//\n"
        )
        assert len(_LOCUS_RE.findall(data)) == 3

    def test_does_not_count_locus_in_comment_block(self):
        # "LOCUS " mid-record (e.g. inside a COMMENT or REFERENCE) must not
        # be counted as a record.  The regex anchors at line start.
        data = (
            "LOCUS       NC_001234           500 bp\n"
            "COMMENT     This entry references LOCUS NC_999999 in another DB.\n"
            "//\n"
        )
        assert len(_LOCUS_RE.findall(data)) == 1

    def test_empty_string_is_zero(self):
        assert len(_LOCUS_RE.findall("")) == 0

