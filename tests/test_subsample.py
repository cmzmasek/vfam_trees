"""Tests for vfam_trees.subsample.proportional_merge."""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees.subsample import proportional_merge


def _rec(rec_id):
    return SeqRecord(Seq("ATGC"), id=rec_id)


def _reps(sp_name, n):
    return [_rec(f"{sp_name}_{i}") for i in range(n)]


# ---------------------------------------------------------------------------
# Basic behaviour
# ---------------------------------------------------------------------------

def test_empty_input_returns_empty():
    assert proportional_merge({}, 10) == []


def test_all_empty_species_returns_empty():
    assert proportional_merge({"A": [], "B": []}, 10) == []


def test_total_reps_leq_target_uses_all():
    species = {"A": _reps("A", 3), "B": _reps("B", 2)}
    out = proportional_merge(species, 10)
    assert len(out) == 5


# ---------------------------------------------------------------------------
# Proportional allocation
# ---------------------------------------------------------------------------

def test_proportional_allocation_hits_target():
    species = {"A": _reps("A", 50), "B": _reps("B", 30), "C": _reps("C", 20)}
    out = proportional_merge(species, 10)
    assert len(out) == 10


def test_proportional_allocation_preserves_ratios():
    species = {"A": _reps("A", 80), "B": _reps("B", 20)}
    out = proportional_merge(species, 10)
    assert len(out) == 10
    a_count = sum(1 for r in out if r.id.startswith("A_"))
    b_count = sum(1 for r in out if r.id.startswith("B_"))
    # ~8:2 ratio (allow ±1 for rounding)
    assert 7 <= a_count <= 9
    assert 1 <= b_count <= 3


# ---------------------------------------------------------------------------
# n_species > target: pathological case that the old code couldn't fix
# ---------------------------------------------------------------------------

def test_hits_target_when_more_species_than_target():
    # 20 species × 1 rep each, target 5. Old code returned 20; fix returns 5.
    species = {f"SP{i:02d}": _reps(f"SP{i:02d}", 1) for i in range(20)}
    out = proportional_merge(species, 5)
    assert len(out) == 5


def test_hits_target_when_every_species_at_quota_floor():
    # 10 species × 1 rep, target 4
    species = {f"SP{i}": _reps(f"SP{i}", 1) for i in range(10)}
    out = proportional_merge(species, 4)
    assert len(out) == 4


def test_prefers_species_with_more_reps_when_oversubscribed():
    # 5 species, rep counts 10/5/3/2/1, target 3 → top 3 by rep count win
    species = {
        "big":   _reps("big",   10),
        "mid":   _reps("mid",   5),
        "small": _reps("small", 3),
        "tiny":  _reps("tiny",  2),
        "one":   _reps("one",   1),
    }
    out = proportional_merge(species, 3)
    assert len(out) == 3
    chosen_prefixes = {r.id.split("_")[0] for r in out}
    assert chosen_prefixes == {"big", "mid", "small"}


# ---------------------------------------------------------------------------
# Priority IDs (RefSeqs)
# ---------------------------------------------------------------------------

def test_priority_ids_preferred_when_subsampling_species():
    # Two species, each oversubscribed; RefSeq rep should survive.
    species = {
        "A": [_rec("A_refseq"), _rec("A_2"), _rec("A_3"), _rec("A_4")],
        "B": [_rec("B_1"), _rec("B_2"), _rec("B_3"), _rec("B_4")],
    }
    out = proportional_merge(species, 2, priority_ids={"A_refseq"})
    ids = {r.id for r in out}
    assert "A_refseq" in ids


def test_priority_ids_survive_when_oversubscribed_at_species_level():
    # 10 species × 2 reps each, target 3 — only 3 species pass through.
    # RefSeq rep in chosen species should be retained.
    species = {}
    for i in range(10):
        species[f"SP{i}"] = [_rec(f"SP{i}_refseq"), _rec(f"SP{i}_2")]
    priority = {f"SP{i}_refseq" for i in range(10)}
    out = proportional_merge(species, 3, priority_ids=priority)
    assert len(out) == 3
    # Each chosen species contributes its RefSeq (ties broken by name)
    assert all(r.id.endswith("_refseq") for r in out)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def test_same_seed_gives_same_result():
    species = {f"SP{i}": _reps(f"SP{i}", 5) for i in range(6)}
    a = proportional_merge(species, 10, seed=123)
    b = proportional_merge(species, 10, seed=123)
    assert [r.id for r in a] == [r.id for r in b]
