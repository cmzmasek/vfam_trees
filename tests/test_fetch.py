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
    fetch_accessions_directly,
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


def test_nuc_marker_also_searches_title():
    """Regression: nucleotide region searches must include [Title] so that
    hantavirus / other records annotated by name in the title (not the [Gene]
    field) are found when seq_type is 'nucleotide' with a non-whole_genome
    region (e.g. 'glycoprotein g1')."""
    q = _build_species_query(12345, "nucleotide", "glycoprotein g1")
    assert '"glycoprotein g1"[Title]' in q
    assert '"glycoprotein g1"[Gene]' in q


def test_nuc_marker_title_and_gene_or_combined():
    """Title and Gene must be OR-ed inside a parenthesised clause, not added
    as independent AND clauses."""
    q = _build_species_query(12345, "nucleotide", "RNA polymerase")
    # Both fields must appear and must be within the same OR clause
    title_term = '"RNA polymerase"[Title]'
    gene_term = '"RNA polymerase"[Gene]'
    assert title_term in q
    assert gene_term in q
    # The OR must sit between them — find both sides of OR
    ti = q.index(title_term)
    gi = q.index(gene_term)
    between = q[min(ti, gi):max(ti, gi) + len(gene_term if gi > ti else title_term)]
    assert " OR " in between


def test_marker_gene_excludes_complete_genome():
    q = _build_species_query(12345, "nucleotide", "B646L")
    assert 'NOT "complete genome"[Title]' in q
    assert 'NOT "complete sequence"[Title]' in q


def test_marker_gene_hexon_excludes_complete():
    q = _build_species_query(12345, "nucleotide", "hexon")
    assert '"hexon"[Gene]' in q
    assert '"hexon"[Title]' in q
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


def test_protein_marker_includes_title_field():
    """Regression for v1.2.30: NCBI's [Protein Name] index is sparsely
    populated for several viral families (Marseilleviridae, Pithoviridae,
    Mimiviridae, ...).  The query must also search [Title] so the FASTA
    defline catches records that aren't indexed under [Protein Name]."""
    q = _build_species_query(12345, "protein", "DNA polymerase")
    assert '"DNA polymerase"[Title]' in q


def test_protein_marker_also_includes_gene_fallback():
    q = _build_species_query(12345, "protein", "B646L")
    assert '"B646L"[Protein Name]' in q
    assert '"B646L"[Title]' in q
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


def test_concat_marker_query_includes_title_field():
    """Regression for v1.2.30: every alias must also be searched in [Title]
    to catch records where NCBI didn't index [Protein Name]."""
    marker = {"name": "DNA polymerase", "aliases": ["polymerase B"]}
    q = _build_marker_query(taxid=944644, marker=marker, species_lineage=None,
                            exclude_organisms=None)
    for n in ("DNA polymerase", "polymerase B"):
        assert f'"{n}"[Title]' in q


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

    def test_protein_refseq_prefix_in_db_source_rejected(self):
        # RefSeq protein records have db_source describing the protein's
        # *own* accession (YP_, NP_, XP_, etc.) — not the source nucleotide.
        # Picking it up and feeding it to nuccore esummary makes NCBI return
        # an "Otherdb db=protein" error that aborts the whole batch.
        rec = SeqRecord(Seq("M" * 100), id="YP_009047263.1")
        rec.annotations = {"db_source": "REFSEQ: accession YP_009047263.1"}
        rec.features = []
        assert _source_nuc_accession(rec) == ""

    def test_protein_refseq_prefix_in_coded_by_rejected(self):
        # Defensive: if a protein accession somehow appears inside a coded_by
        # qualifier, fall through to the next candidate rather than returning it.
        rec = SeqRecord(Seq("M" * 100), id="YP_1")
        cds = SeqFeature(FeatureLocation(0, 100), type="CDS")
        cds.qualifiers = {"coded_by": ["YP_009047263.1, NC_001234.1:1..300"]}
        rec.features = [cds]
        rec.annotations = {}
        assert _source_nuc_accession(rec) == "NC_001234.1"

    def test_three_letter_genbank_protein_in_db_source_rejected(self):
        # 3-letter, no-underscore accessions like "ABO61246.1" are GenBank
        # *protein* accessions, not nucleotide.  They surface in the
        # db_source of GenBank protein records as the protein's own home
        # accession, and would otherwise abort a nuccore esummary batch
        # with an "Otherdb db=protein" error.
        rec = SeqRecord(Seq("M" * 100), id="ABO61246.1")
        rec.annotations = {"db_source": "accession ABO61246.1"}
        rec.features = []
        assert _source_nuc_accession(rec) == ""

    def test_three_letter_genbank_protein_in_coded_by_rejected(self):
        rec = SeqRecord(Seq("M" * 100), id="ABO61246.1")
        cds = SeqFeature(FeatureLocation(0, 100), type="CDS")
        # Defensive — fall through to the real source nuc later in the qualifier.
        cds.qualifiers = {"coded_by": ["ABO61246.1, NC_001234.1:1..300"]}
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

    def test_protein_refseq_prefixes_filtered(self, monkeypatch):
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

        monkeypatch.setattr(fetch_mod.Entrez, "esummary", fake_esummary)
        monkeypatch.setattr(fetch_mod.Entrez, "read", lambda h: [])
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        fetch_mod.fetch_nuc_lengths(
            ["YP_009047263.1", "NP_009123456", "XP_555.1",
             "AP_001.1", "WP_002.1", "ZP_003.1", "ELP_004.1",
             "NC_001234.1", "AB123456"]
        )
        assert captured_batches == [["NC_001234.1", "AB123456"]]

    def test_three_letter_genbank_protein_filtered(self, monkeypatch):
        # 3-letter, no-underscore GenBank protein accessions (e.g. ABO61246.1)
        # match the nuccore-shape regex by length but resolve in the protein
        # database — must be filtered before esummary.
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

        monkeypatch.setattr(fetch_mod.Entrez, "esummary", fake_esummary)
        monkeypatch.setattr(fetch_mod.Entrez, "read", lambda h: [])
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        fetch_mod.fetch_nuc_lengths(
            ["ABO61246.1", "AAC54321", "BAB12345.1",  # 3-letter GenBank protein
             "NC_001234.1", "AB123456"]               # real nuc accessions
        )
        assert captured_batches == [["NC_001234.1", "AB123456"]]

    def test_otherdb_error_triggers_binary_split_recovery(self, monkeypatch):
        # When NCBI returns a deterministic "Otherdb" error for a batch with
        # one bad UID, retrying changes nothing; the recovery is to split
        # the batch and isolate the bad UID, preserving lengths for the rest.
        from vfam_trees import fetch as fetch_mod

        # Shape: 1 letter + 5 digits.  Passes the nuccore-shape filter (it's
        # the same shape as a 1+5 GenBank nuc) but in real data this form is
        # often a UniProt accession that resolves to the protein db.
        bad_uid = "P12345"

        # Track every esummary call's input batch.
        calls: list[list[str]] = []

        class _FakeHandle:
            def __init__(self, summaries):
                self._summaries = summaries
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False

        def fake_esummary(db, id):
            ids = id.split(",")
            calls.append(ids)
            if bad_uid in ids:
                raise RuntimeError(
                    f'Otherdb uid="42" db="protein" term="{bad_uid}"'
                )
            # Build summaries for each good ID (length proportional to index).
            summaries = [
                {"AccessionVersion": acc, "Caption": acc.split(".")[0],
                 "Length": str(1000 * (i + 1))}
                for i, acc in enumerate(ids)
            ]
            return _FakeHandle(summaries)

        def fake_read(handle):
            return handle._summaries

        monkeypatch.setattr(fetch_mod.Entrez, "esummary", fake_esummary)
        monkeypatch.setattr(fetch_mod.Entrez, "read", fake_read)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        good = ["NC_000001.1", "NC_000002.1", "NC_000003.1"]
        result = fetch_mod.fetch_nuc_lengths(good + [bad_uid])

        # All three good accessions should have lengths recovered after splitting.
        assert set(result.keys()) == set(good)
        assert all(result[a] > 0 for a in good)
        # The first call hits Otherdb (bad in the batch); subsequent splits
        # isolate the bad UID without retrying the original 5 times.
        assert calls[0] == good + [bad_uid]
        # Bad UID ends up alone in some later batch and fails — that's expected;
        # the important thing is we don't retry the full batch 5x first.
        assert sum(1 for c in calls if c == (good + [bad_uid])) == 1

    def test_db_source_fallback_only_accepts_refseq_nuc_prefixes(self):
        from vfam_trees.fetch import _source_nuc_accession

        # GenBank protein accession in db_source — REJECTED.
        rec = SeqRecord(Seq("M" * 100), id="ABO61246.1")
        rec.annotations = {"db_source": "accession ABO61246.1"}
        rec.features = []
        assert _source_nuc_accession(rec) == ""

        # UniProt-style accession (1 letter + 5 digits) — REJECTED, even
        # though it matches the nuccore-shape regex.
        rec.annotations = {"db_source": "UniProtKB: accession P12345"}
        assert _source_nuc_accession(rec) == ""

        # RefSeq nuc accession in db_source — accepted.
        rec.annotations = {"db_source": "REFSEQ: accession NC_999999.1"}
        assert _source_nuc_accession(rec) == "NC_999999.1"


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


# ---------------------------------------------------------------------------
# fetch_accessions_directly — auto-routing nuccore vs protein
# ---------------------------------------------------------------------------

# Minimal GenBank flat-file stub for two records: one nuccore, one protein.
_NUC_GB = """\
LOCUS       KC567260                 100 bp    RNA     linear   VRL 01-JAN-2020
DEFINITION  Orthohantavirus andesense glycoprotein G1 gene.
ACCESSION   KC567260
VERSION     KC567260.1
ORGANISM    Orthohantavirus andesense
            Viruses; Riboviria; Orthornavirae; Negarnaviricota;
            Ellioviricetes; Bunyavirales; Hantaviridae.
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="Orthohantavirus andesense"
                     /mol_type="genomic RNA"
                     /db_xref="taxon:1980456"
     gene            1..100
                     /gene="G1"
ORIGIN
        1 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg
       61 atgatgatga tgatgatgat gatgatgatg
//
"""

_PROT_GB = """\
LOCUS       YP_009999999             50 aa     linear   VRL 01-JAN-2020
DEFINITION  glycoprotein G1 [Test virus].
ACCESSION   YP_009999999
VERSION     YP_009999999.1
ORGANISM    Test virus
            Viruses; Riboviria.
FEATURES             Location/Qualifiers
     source          1..50
                     /organism="Test virus"
                     /db_xref="taxon:99999"
ORIGIN
        1 mmmmmmmmmm mmmmmmmmmm mmmmmmmmmm mmmmmmmmmm mmmmmmmmmm
//
"""


class TestFetchAccessionsdirectly:
    """Regression: fetch_accessions_directly auto-routes to the correct NCBI db
    and returns parsed SeqRecords without hitting the network."""

    def _make_efetch(self, monkeypatch, nuc_data: str, prot_data: str):
        """Patch Entrez.efetch to return canned data per db, and disable sleep."""
        import vfam_trees.fetch as fetch_mod

        class _FakeHandle:
            def __init__(self, text):
                self._text = text
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
            def read(self):
                return self._text

        def fake_efetch(db, id, rettype, retmode):
            return _FakeHandle(nuc_data if db == "nuccore" else prot_data)

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", fake_efetch)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

    def test_nuccore_accession_routed_to_nuccore(self, monkeypatch):
        """A 2+6 nuccore accession (e.g. KC567260.1) must be fetched from nuccore."""
        import vfam_trees.fetch as fetch_mod
        calls: list[str] = []

        class _FakeHandle:
            def __init__(self, text):
                self._text = text
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
            def read(self):
                return self._text

        def fake_efetch(db, id, rettype, retmode):
            calls.append(db)
            return _FakeHandle(_NUC_GB if db == "nuccore" else "")

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", fake_efetch)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        recs = fetch_accessions_directly(["KC567260.1"])
        assert "nuccore" in calls
        assert "protein" not in calls
        assert len(recs) == 1
        assert recs[0].id == "KC567260.1"

    def test_protein_accession_routed_to_protein(self, monkeypatch):
        """A RefSeq protein accession (YP_…) must be fetched from the protein db."""
        import vfam_trees.fetch as fetch_mod
        calls: list[str] = []

        class _FakeHandle:
            def __init__(self, text):
                self._text = text
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
            def read(self):
                return self._text

        def fake_efetch(db, id, rettype, retmode):
            calls.append(db)
            return _FakeHandle(_PROT_GB if db == "protein" else "")

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", fake_efetch)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        recs = fetch_accessions_directly(["YP_009999999.1"])
        assert "protein" in calls
        assert "nuccore" not in calls
        assert len(recs) == 1

    def test_mixed_accessions_split_across_databases(self, monkeypatch):
        """Nuccore and protein accessions in one call must be routed to separate dbs."""
        import vfam_trees.fetch as fetch_mod
        calls: list[str] = []

        class _FakeHandle:
            def __init__(self, text):
                self._text = text
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
            def read(self):
                return self._text

        def fake_efetch(db, id, rettype, retmode):
            calls.append(db)
            return _FakeHandle(_NUC_GB if db == "nuccore" else _PROT_GB)

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", fake_efetch)
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        recs = fetch_accessions_directly(["KC567260.1", "YP_009999999.1"])
        assert "nuccore" in calls
        assert "protein" in calls
        assert len(recs) == 2

    def test_empty_ncbi_response_returns_empty_list(self, monkeypatch):
        """When NCBI returns nothing, the function must return [] without error."""
        import vfam_trees.fetch as fetch_mod

        class _FakeHandle:
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
            def read(self):
                return ""

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", lambda **_: _FakeHandle())
        monkeypatch.setattr(fetch_mod.time, "sleep", lambda *_: None)

        recs = fetch_accessions_directly(["KC567260.1"])
        assert recs == []

    def test_empty_input_returns_empty_list(self, monkeypatch):
        """Empty accession set must short-circuit without any network call."""
        import vfam_trees.fetch as fetch_mod
        called = {"n": 0}

        def boom(*a, **kw):
            called["n"] += 1

        monkeypatch.setattr(fetch_mod.Entrez, "efetch", boom)
        recs = fetch_accessions_directly([])
        assert recs == []
        assert called["n"] == 0


# ---------------------------------------------------------------------------
# HMM fetch mode (all-proteins-per-species → HMM marker assignment)
# ---------------------------------------------------------------------------

def _hmm_hit(target, *, ali_from=1, ali_to=100):
    return {"target": target, "passing": True, "dom_evalue": 1e-30,
            "ali_from": ali_from, "ali_to": ali_to}


def _hprot(acc, source_acc, length, hits, description="protein"):
    """Protein record with a coded_by source genome AND pre-computed HMM hits
    (so the real HMMIdentifier runs without hmmscan)."""
    rec = _protein(acc, description=description, length=length, source_acc=source_acc)
    rec.annotations["hmm_hits"] = hits
    return rec


class TestNoNameQuery:
    def test_all_proteins_query_has_no_name_clause(self):
        from vfam_trees.fetch import _build_all_proteins_query
        q = _build_all_proteins_query(10244, None)
        assert q.startswith("txid10244[Organism:exp]")
        assert "[Protein Name]" not in q and "[Title]" not in q
        assert "NOT patent[filter]" in q

    def test_refseq_and_exclude(self):
        from vfam_trees.fetch import _build_all_proteins_query
        q = _build_all_proteins_query(1, ["Homo sapiens"], refseq_only=True)
        assert "refseq[filter]" in q
        assert 'NOT "Homo sapiens"[Organism]' in q


class TestCapGenomes:
    def test_caps_refseq_first(self):
        from vfam_trees.fetch import _cap_genomes_refseq_first
        g = {"NC_1.1": {}, "AY_2.1": {}, "AY_3.1": {}, "AY_1.1": {}}
        kept = _cap_genomes_refseq_first(g, 2)
        assert "NC_1.1" in kept          # RefSeq always kept
        assert len(kept) == 2

    def test_no_cap_when_zero_or_under(self):
        from vfam_trees.fetch import _cap_genomes_refseq_first
        g = {"NC_1.1": {}, "AY_2.1": {}}
        assert len(_cap_genomes_refseq_first(g, 0)) == 2
        assert len(_cap_genomes_refseq_first(g, 5)) == 2


class TestHmmFetch:
    def test_fetch_species_genomes_hmm(self, monkeypatch, tmp_path):
        from vfam_trees import fetch as fetch_mod
        from vfam_trees.markers import HMMIdentifier
        recs = [
            _hprot("P1", "NC_001.1", 900, [_hmm_hit("Adeno_hexon")], "hexon"),
            _hprot("P2", "NC_001.1", 400, [_hmm_hit("Adeno_penton")], "penton"),
            _hprot("P3", "AB_002.1", 880, [_hmm_hit("Adeno_hexon")], "hexon"),
        ]
        monkeypatch.setattr(fetch_mod, "_fetch_all_species_proteins",
                            lambda *a, **k: recs)
        ident = HMMIdentifier({"hmm": {"enabled": True}})
        marker_set = [{"name": "hexon", "hmms": ["Adeno_hexon"]},
                      {"name": "penton", "hmms": ["Adeno_penton"]}]
        genomes, stats = fetch_mod._fetch_species_genomes_hmm(
            taxid=1, species_name="sp", species_lineage=None,
            marker_set=marker_set, output_dir=tmp_path, max_per_species=200,
            min_fraction=0.5, exclude_organisms=None, identifier=ident,
            source_nuc_min_length_frac=0.0, nuc_length_lookup=None,
        )
        assert set(genomes) == {"NC_001.1", "AB_002.1"}
        assert genomes["NC_001.1"]["hexon"].id == "P1"
        assert genomes["NC_001.1"]["penton"].id == "P2"
        assert "penton" not in genomes["AB_002.1"]
        assert stats["n_proteins_fetched"] == 3
        assert stats["n_genomes_kept"] == 2

    def test_genome_cap_applied(self, monkeypatch, tmp_path):
        from vfam_trees import fetch as fetch_mod
        from vfam_trees.markers import HMMIdentifier
        recs = [
            _hprot("P1", "NC_001.1", 900, [_hmm_hit("Adeno_hexon")], "hexon"),
            _hprot("P3", "AY_002.1", 880, [_hmm_hit("Adeno_hexon")], "hexon"),
            _hprot("P4", "AY_003.1", 870, [_hmm_hit("Adeno_hexon")], "hexon"),
        ]
        monkeypatch.setattr(fetch_mod, "_fetch_all_species_proteins",
                            lambda *a, **k: recs)
        ident = HMMIdentifier({"hmm": {"enabled": True}})
        genomes, stats = fetch_mod._fetch_species_genomes_hmm(
            taxid=1, species_name="sp", species_lineage=None,
            marker_set=[{"name": "hexon", "hmms": ["Adeno_hexon"]}],
            output_dir=tmp_path, max_per_species=2, min_fraction=1.0,
            exclude_organisms=None, identifier=ident,
            source_nuc_min_length_frac=0.0, nuc_length_lookup=None,
        )
        assert len(genomes) == 2               # capped from 3
        assert "NC_001.1" in genomes           # RefSeq survives the cap
        assert stats["n_dropped_genome_cap"] == 1

    def test_fetch_species_marker_hmm_writes_selected(self, monkeypatch, tmp_path):
        from Bio import SeqIO
        from vfam_trees import fetch as fetch_mod
        from vfam_trees.markers import HMMIdentifier
        recs = [
            _hprot("P1", "NC_001.1", 900, [_hmm_hit("Adeno_hexon")], "hexon"),
            _hprot("P2", "NC_001.1", 400, [_hmm_hit("Adeno_penton")], "penton"),
            _hprot("P3", "AB_002.1", 880, [_hmm_hit("Adeno_hexon")], "hexon"),
        ]
        monkeypatch.setattr(fetch_mod, "_fetch_all_species_proteins",
                            lambda *a, **k: recs)
        ident = HMMIdentifier({"hmm": {"enabled": True}})
        out = tmp_path / "out.gb"
        n = fetch_mod._fetch_species_marker_hmm(
            taxid=1, species_name="sp",
            marker={"name": "hexon", "hmms": ["Adeno_hexon"]},
            output_gb=out, max_per_species=200, exclude_organisms=None,
            identifier=ident, species_lineage=None,
        )
        assert n == 2  # hexon selected in both genomes
        parsed = list(SeqIO.parse(out, "genbank"))
        assert {r.id for r in parsed} == {"P1", "P3"}

    def test_fetch_species_sequences_routes_to_hmm(self, monkeypatch, tmp_path):
        from vfam_trees import fetch as fetch_mod

        captured = {}

        def fake_marker_hmm(**kwargs):
            captured.update(kwargs)
            return 7

        monkeypatch.setattr(fetch_mod, "_fetch_species_marker_hmm", fake_marker_hmm)

        class _StubHmm:
            is_hmm = True

        n = fetch_mod.fetch_species_sequences(
            taxid=1, species_name="s", seq_type="protein", region="hexon",
            output_gb=tmp_path / "x.gb", hmm_identifier=_StubHmm(),
            hmm_marker={"name": "hexon", "hmms": ["Adeno_hexon--Adeno_hexon_C"]},
        )
        assert n == 7
        assert captured["marker"]["name"] == "hexon"

    def test_name_mode_unaffected_when_no_hmm_identifier(self, monkeypatch, tmp_path):
        # Without an HMM identifier the classic name-query path is used.
        from vfam_trees import fetch as fetch_mod
        monkeypatch.setattr(fetch_mod, "_search_ids", lambda *a, **k: [])
        n = fetch_mod.fetch_species_sequences(
            taxid=1, species_name="s", seq_type="protein", region="hexon",
            output_gb=tmp_path / "x.gb",
        )
        assert n == 0

