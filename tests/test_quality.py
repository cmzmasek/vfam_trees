"""Tests for vfam_trees.quality."""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from vfam_trees.quality import filter_sequences, remove_length_outliers


def _rec(seq, taxid="12345", acc="ACC1", organism="Test virus"):
    r = SeqRecord(Seq(seq), id=acc, description="")
    r.annotations["organism"] = organism
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


# ---------------------------------------------------------------------------
# remove_length_outliers
# ---------------------------------------------------------------------------

def test_remove_length_outliers_keeps_normal():
    records = [_rec("A" * 300), _rec("A" * 310), _rec("A" * 290)]
    kept, n_removed = remove_length_outliers(records)
    assert n_removed == 0
    assert len(kept) == 3


def test_remove_length_outliers_drops_extreme():
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 30000)]
    kept, n_removed = remove_length_outliers(records)
    assert n_removed == 1
    assert len(kept) == 10
