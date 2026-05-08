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

# Realistic small variation around 300 → MAD on log-lengths is non-zero,
# so the filter has signal to work with.  Pure [300]*N synthetic data
# yields MAD = 0 and the filter correctly degenerates to a no-op.
_NORMAL_LENGTHS = (290, 295, 298, 300, 300, 302, 305, 308, 310, 312)


def test_remove_length_outliers_keeps_normal():
    records = [_rec("A" * L) for L in _NORMAL_LENGTHS]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == len(_NORMAL_LENGTHS)


def test_remove_length_outliers_drops_extreme_long():
    records = [_rec("A" * L) for L in _NORMAL_LENGTHS] + [_rec("A" * 30000)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 1
    assert n_short == 0
    assert len(kept) == len(_NORMAL_LENGTHS)


def test_remove_length_outliers_drops_short():
    records = [_rec("A" * L) for L in _NORMAL_LENGTHS] + [_rec("A" * 30)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 0
    assert n_short == 1
    assert len(kept) == len(_NORMAL_LENGTHS)


def test_remove_length_outliers_two_sided():
    records = (
        [_rec("A" * L) for L in _NORMAL_LENGTHS]
        + [_rec("A" * 30000), _rec("A" * 30)]
    )
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 1
    assert n_short == 1
    assert len(kept) == len(_NORMAL_LENGTHS)


def test_remove_length_outliers_disabled():
    # k=0 disables the filter entirely
    records = [_rec("A" * L) for L in _NORMAL_LENGTHS] + [_rec("A" * 30)]
    kept, n_long, n_short = remove_length_outliers(records, k=0)
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == len(_NORMAL_LENGTHS) + 1


def test_remove_length_outliers_zero_mad_is_noop():
    # All lengths identical → MAD = 0 → no characterization of natural
    # spread possible → filter is a no-op (does not fabricate a fallback).
    records = [_rec("A" * 300)] * 10 + [_rec("A" * 30000)]
    kept, n_long, n_short = remove_length_outliers(records)
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == 11


def test_remove_length_outliers_larger_k_more_lenient():
    # Borderline-short sequence is dropped at k=5 but kept at k=20.
    borderline = (240,)
    records = [_rec("A" * L) for L in _NORMAL_LENGTHS + borderline]
    _, _, n_short_k5 = remove_length_outliers(records, k=5.0)
    _, _, n_short_k20 = remove_length_outliers(records, k=20.0)
    assert n_short_k5 >= n_short_k20


def test_remove_length_outliers_protects_refseq_long(caplog):
    normal = [_rec("A" * L, acc=f"N{i}") for i, L in enumerate(_NORMAL_LENGTHS)]
    refseq = _rec("A" * 30000, acc="NC_000001")
    records = normal + [refseq]
    with caplog.at_level("WARNING"):
        kept, n_long, n_short = remove_length_outliers(
            records, protected_ids={"NC_000001"},
        )
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == len(normal) + 1
    assert "NC_000001" in {r.id for r in kept}
    assert any("NC_000001" in m and "protected" in m for m in caplog.messages)


def test_remove_length_outliers_protects_refseq_short(caplog):
    normal = [_rec("A" * L, acc=f"N{i}") for i, L in enumerate(_NORMAL_LENGTHS)]
    refseq = _rec("A" * 30, acc="NC_000002")
    records = normal + [refseq]
    with caplog.at_level("WARNING"):
        kept, n_long, n_short = remove_length_outliers(
            records, protected_ids={"NC_000002"},
        )
    assert n_long == 0
    assert n_short == 0
    assert len(kept) == len(normal) + 1
    assert "NC_000002" in {r.id for r in kept}
    assert any("NC_000002" in m and "protected" in m for m in caplog.messages)


def test_remove_length_outliers_non_protected_still_dropped():
    normal = [_rec("A" * L, acc=f"N{i}") for i, L in enumerate(_NORMAL_LENGTHS)]
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
