"""Tests for vfam_trees.subsample.proportional_merge."""
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vfam_trees import subsample
from vfam_trees.subsample import absorb_into_refseqs, proportional_merge


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


# ---------------------------------------------------------------------------
# absorb_into_refseqs
# ---------------------------------------------------------------------------

def _mock_membership(monkeypatch, clusters: list[set[str]]):
    """Patch _cluster_membership to return *clusters* verbatim."""
    def _fake(records, threshold, seq_type, work_dir, clustering_tool):
        return clusters
    monkeypatch.setattr(subsample, "_cluster_membership", _fake)


def test_absorb_no_records_returns_empty(tmp_path):
    kept, n = absorb_into_refseqs([], {"NC_1"}, 0.99, "nucleotide", tmp_path)
    assert kept == []
    assert n == 0


def test_absorb_single_record_unchanged(tmp_path):
    rec = _rec("NC_1")
    kept, n = absorb_into_refseqs([rec], {"NC_1"}, 0.99, "nucleotide", tmp_path)
    assert kept == [rec]
    assert n == 0


def test_absorb_no_refseqs_in_records_skips_clustering(tmp_path, monkeypatch):
    # Should never call _cluster_membership when no RefSeqs are present
    def _boom(*a, **k):
        raise AssertionError("_cluster_membership should not be called")
    monkeypatch.setattr(subsample, "_cluster_membership", _boom)
    records = [_rec("X1"), _rec("X2"), _rec("X3")]
    kept, n = absorb_into_refseqs(records, {"NC_1"}, 0.99, "nucleotide", tmp_path)
    assert [r.id for r in kept] == ["X1", "X2", "X3"]
    assert n == 0


def test_absorb_drops_non_refseq_in_refseq_cluster(tmp_path, monkeypatch):
    # NC_1 (RefSeq) clusters with X1, X2 (non-RefSeqs) → drop X1, X2
    records = [_rec("NC_1"), _rec("X1"), _rec("X2"), _rec("Y1")]
    _mock_membership(monkeypatch, [{"NC_1", "X1", "X2"}, {"Y1"}])
    kept, n = absorb_into_refseqs(records, {"NC_1"}, 0.99, "nucleotide", tmp_path)
    kept_ids = {r.id for r in kept}
    assert kept_ids == {"NC_1", "Y1"}
    assert n == 2


def test_absorb_keeps_all_refseqs_in_one_cluster(tmp_path, monkeypatch):
    # Two RefSeqs near-identical to each other → both kept
    records = [_rec("NC_1"), _rec("NC_2"), _rec("X1")]
    _mock_membership(monkeypatch, [{"NC_1", "NC_2", "X1"}])
    kept, n = absorb_into_refseqs(
        records, {"NC_1", "NC_2"}, 0.99, "nucleotide", tmp_path,
    )
    kept_ids = {r.id for r in kept}
    assert kept_ids == {"NC_1", "NC_2"}
    assert n == 1


def test_absorb_passes_through_non_refseq_clusters(tmp_path, monkeypatch):
    # Cluster {X1, X2} has no RefSeq → both kept; cluster {NC_1, Y1} → Y1 absorbed
    records = [_rec("NC_1"), _rec("X1"), _rec("X2"), _rec("Y1")]
    _mock_membership(monkeypatch, [{"NC_1", "Y1"}, {"X1", "X2"}])
    kept, n = absorb_into_refseqs(records, {"NC_1"}, 0.99, "nucleotide", tmp_path)
    kept_ids = {r.id for r in kept}
    assert kept_ids == {"NC_1", "X1", "X2"}
    assert n == 1


def test_absorb_refseq_alone_in_cluster_unchanged(tmp_path, monkeypatch):
    records = [_rec("NC_1"), _rec("X1"), _rec("X2")]
    _mock_membership(monkeypatch, [{"NC_1"}, {"X1", "X2"}])
    kept, n = absorb_into_refseqs(records, {"NC_1"}, 0.99, "nucleotide", tmp_path)
    kept_ids = {r.id for r in kept}
    assert kept_ids == {"NC_1", "X1", "X2"}
    assert n == 0
