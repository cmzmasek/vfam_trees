"""Tests for vfam_trees.quality."""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees.quality import filter_sequences, remove_length_outliers


def _rec(seq, taxid="12345", acc="ACC1", organism="Test virus",
         source="", description=""):
    r = SeqRecord(Seq(seq), id=acc, description=description)
    r.annotations["organism"] = organism
    r.annotations["source"] = source
    r.features = []
    r.dbxrefs = [f"taxon:{taxid}"]
    return r


# ---------------------------------------------------------------------------
# filter_sequences
# ---------------------------------------------------------------------------

GOOD_SEQ = "ATCGATCG" * 40   # 320 bp — above MIN_LENGTH_NUC (200)


def test_passes_good_sequence():
    records = [_rec(GOOD_SEQ)]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=["synthetic construct"],
    )
    assert len(passed) == 1


def test_excludes_by_organism():
    records = [
        _rec(GOOD_SEQ, organism="synthetic construct"),
        _rec(GOOD_SEQ, organism="Real virus"),
    ]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=["synthetic construct"],
    )
    assert len(passed) == 1
    assert stats["n_excluded_organism"] == 1


def test_excludes_by_source():
    records = [
        _rec(GOOD_SEQ, organism="Real virus", source="metagenome"),
        _rec(GOOD_SEQ, organism="Real virus", source="Homo sapiens"),
    ]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=["metagenome"],
    )
    assert len(passed) == 1
    assert stats["n_excluded_organism"] == 1


def test_excludes_by_definition():
    records = [
        _rec(GOOD_SEQ, organism="Real virus",
             description="MAG: virus isolate xyz, complete genome"),
        _rec(GOOD_SEQ, organism="Real virus",
             description="Virus strain ABC, complete genome"),
    ]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=["MAG:"],
    )
    assert len(passed) == 1
    assert stats["n_excluded_organism"] == 1


def test_mag_colon_not_matched_by_unrelated_name():
    # "Magnolia" does not contain "MAG:" so should pass
    records = [_rec(GOOD_SEQ, organism="Magnolia virus 1")]
    passed, _, _ = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=["MAG:"],
    )
    assert len(passed) == 1


def test_excludes_by_ambiguity():
    records = [_rec("N" * 320)]  # above floor, but 100% ambiguous
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    assert len(passed) == 0
    assert stats["n_excluded_ambiguity"] == 1


def test_excludes_by_min_length():
    records = [_rec(GOOD_SEQ)]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=10000,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    assert len(passed) == 0
    assert stats["n_excluded_length"] == 1


def test_accepts_when_no_min_length():
    records = [_rec(GOOD_SEQ)]
    passed, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    assert len(passed) == 1


def test_pre_length_lengths_present_in_stats():
    records = [_rec(GOOD_SEQ)]
    _, _, stats = filter_sequences(
        records, seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    assert "pre_length_lengths" in stats
    assert stats["pre_length_lengths"] == [len(GOOD_SEQ)]


def test_pre_length_lengths_excludes_ambiguous():
    ambig = _rec("N" * 320, acc="AMB")
    clean = _rec(GOOD_SEQ, acc="OK")
    _, _, stats = filter_sequences(
        [ambig, clean], seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    # Ambiguous sequence must not appear in pre_length_lengths
    assert len(stats["pre_length_lengths"]) == 1
    assert stats["pre_length_lengths"][0] == len(GOOD_SEQ)


def test_pre_length_lengths_includes_too_short():
    short = _rec("A" * 210, acc="SHORT")   # above floor (200) but below min_length
    good  = _rec(GOOD_SEQ, acc="GOOD")
    _, _, stats = filter_sequences(
        [short, good], seq_type="nucleotide", min_length=1000,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    # Both passed ambiguity; both should appear in pre_length_lengths
    assert sorted(stats["pre_length_lengths"]) == sorted([210, len(GOOD_SEQ)])
    assert stats["n_excluded_length"] == 2


def test_ambiguity_filter_before_length_filter():
    # A sequence that is both too ambiguous AND too short must be counted
    # as ambiguity-excluded, not length-excluded (ambiguity check runs first).
    both_bad = _rec("N" * 50, acc="BAD")   # < 200 bp floor AND 100% ambiguous
    _, _, stats = filter_sequences(
        [both_bad], seq_type="nucleotide", min_length=None,
        max_ambiguous=0.01, exclude_organisms=[],
    )
    assert stats["n_excluded_ambiguity"] == 1
    assert stats["n_excluded_length"] == 0


# ---------------------------------------------------------------------------
# remove_length_outliers
# ---------------------------------------------------------------------------

def test_remove_length_outliers_keeps_normal():
    records = [_rec("A" * 300), _rec("A" * 310), _rec("A" * 290)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == 3


def test_remove_length_outliers_drops_extreme_long():
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 30000)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 1
    assert n_short == 0
    assert len(kept) == 10


def test_remove_length_outliers_drops_short():
    # median = 300; lo_mult default = 1/3 → cutoff ~100, so a 50 bp seq is dropped
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 50)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 0
    assert n_short == 1
    assert len(kept) == 10


def test_remove_length_outliers_two_sided():
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 30000), _rec("A" * 50)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 1
    assert n_short == 1
    assert len(kept) == 10


def test_remove_length_outliers_lo_disabled():
    # lo_mult=0 restores the old upper-bound-only behaviour
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 50)]
    kept, n_long, n_short = remove_length_outliers(records, lo_mult=0)
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == 11


def test_remove_length_outliers_protects_refseq_long(caplog):
    # Extreme-long sequence flagged but protected → kept, warning logged
    normal = [_rec("A" * 300, acc=f"N{i}") for i in range(10)]
    refseq = _rec("A" * 30000, acc="NC_000001")
    records = normal + [refseq]
    with caplog.at_level("WARNING"):
        kept, n_long, n_short = remove_length_outliers(
            records, protected_ids={"NC_000001"},
        )
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == 11
    assert "NC_000001" in {r.id for r in kept}
    assert any("NC_000001" in m and "protected" in m for m in caplog.messages)


def test_remove_length_outliers_protects_refseq_short(caplog):
    normal = [_rec("A" * 300, acc=f"N{i}") for i in range(10)]
    refseq = _rec("A" * 50, acc="NC_000002")
    records = normal + [refseq]
    with caplog.at_level("WARNING"):
        kept, n_long, n_short = remove_length_outliers(
            records, protected_ids={"NC_000002"},
        )
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == 11
    assert "NC_000002" in {r.id for r in kept}
    assert any("NC_000002" in m and "protected" in m for m in caplog.messages)


def test_remove_length_outliers_non_protected_still_dropped():
    # Mix a protected RefSeq and an unprotected outlier — only the latter drops
    normal = [_rec("A" * 300, acc=f"N{i}") for i in range(10)]
    refseq = _rec("A" * 30000, acc="NC_000001")
    unprotected = _rec("A" * 30000, acc="UNPROT1")
    records = normal + [refseq, unprotected]
    kept, n_long, n_short = remove_length_outliers(
        records, protected_ids={"NC_000001"},
    )
    assert n_long == 1
    assert n_short == 0
    kept_ids = {r.id for r in kept}
    assert "NC_000001" in kept_ids
    assert "UNPROT1" not in kept_ids
