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
    out, stats = proportional_merge({}, 10)
    assert out == []
    assert stats["n_species_dropped_at_cap"] == 0


def test_all_empty_species_returns_empty():
    out, stats = proportional_merge({"A": [], "B": []}, 10)
    assert out == []
    assert stats["n_species_dropped_at_cap"] == 0


def test_total_reps_leq_target_uses_all():
    species = {"A": _reps("A", 3), "B": _reps("B", 2)}
    out, stats = proportional_merge(species, 10)
    assert len(out) == 5
    assert stats["n_species_dropped_at_cap"] == 0


# ---------------------------------------------------------------------------
# Proportional allocation
# ---------------------------------------------------------------------------

def test_proportional_allocation_hits_target():
    species = {"A": _reps("A", 50), "B": _reps("B", 30), "C": _reps("C", 20)}
    out, _ = proportional_merge(species, 10)
    assert len(out) == 10


def test_proportional_allocation_preserves_ratios():
    species = {"A": _reps("A", 80), "B": _reps("B", 20)}
    out, _ = proportional_merge(species, 10)
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
    out, stats = proportional_merge(species, 5)
    assert len(out) == 5
    assert stats["n_species_dropped_at_cap"] == 15
    assert stats["n_refseq_species_dropped_at_cap"] == 0


def test_hits_target_when_every_species_at_quota_floor():
    # 10 species × 1 rep, target 4
    species = {f"SP{i}": _reps(f"SP{i}", 1) for i in range(10)}
    out, stats = proportional_merge(species, 4)
    assert len(out) == 4
    assert stats["n_species_dropped_at_cap"] == 6


def test_prefers_species_with_more_reps_when_oversubscribed():
    # 5 species, rep counts 10/5/3/2/1, target 3 → top 3 by rep count win
    species = {
        "big":   _reps("big",   10),
        "mid":   _reps("mid",   5),
        "small": _reps("small", 3),
        "tiny":  _reps("tiny",  2),
        "one":   _reps("one",   1),
    }
    out, _ = proportional_merge(species, 3)
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
    out, _ = proportional_merge(species, 2, priority_ids={"A_refseq"})
    ids = {r.id for r in out}
    assert "A_refseq" in ids


def test_priority_ids_survive_when_oversubscribed_at_species_level():
    # 10 species × 2 reps each, target 3 — only 3 species pass through.
    # RefSeq rep in chosen species should be retained.
    species = {}
    for i in range(10):
        species[f"SP{i}"] = [_rec(f"SP{i}_refseq"), _rec(f"SP{i}_2")]
    priority = {f"SP{i}_refseq" for i in range(10)}
    out, stats = proportional_merge(species, 3, priority_ids=priority)
    assert len(out) == 3
    # Each chosen species contributes its RefSeq (ties broken by name)
    assert all(r.id.endswith("_refseq") for r in out)
    # All 10 species had a RefSeq → 7 RefSeq-bearing species got dropped at cap.
    assert stats["n_species_dropped_at_cap"] == 7
    assert stats["n_refseq_species_dropped_at_cap"] == 7


# ---------------------------------------------------------------------------
# v1.2.14: RefSeq-bearing species kept first when species count > target
# ---------------------------------------------------------------------------

def test_refseq_bearing_species_kept_first_when_capped():
    # 10 species; only 3 have RefSeqs.  Target 5.  All 3 RefSeq-bearing
    # species must make the cut, regardless of rep counts.  The 5th slot
    # is the largest non-RefSeq species.
    species = {}
    refseq_ids = set()
    # 3 RefSeq-bearing species, each with a single sequence.
    for i in range(3):
        sp = f"REF{i}"
        species[sp] = _reps(sp, 1)
        refseq_ids.add(species[sp][0].id)
    # 7 non-RefSeq species with progressively more reps.
    for i in range(7):
        sp = f"NRF{i}"
        species[sp] = _reps(sp, i + 2)  # 2..8 reps
    out, stats = proportional_merge(species, 5, priority_ids=refseq_ids)
    assert len(out) == 5
    # All 3 RefSeq-bearing species selected.
    chosen_prefixes = {r.id.rsplit("_", 1)[0] for r in out}
    assert {"REF0", "REF1", "REF2"} <= chosen_prefixes
    # 2 of the 7 NRF species filled the remaining slots — the rep-count
    # winners (NRF6 with 8 reps, NRF5 with 7 reps).
    assert "NRF6" in chosen_prefixes
    assert "NRF5" in chosen_prefixes
    # 5 species dropped, 0 of which were RefSeq-bearing.
    assert stats["n_species_dropped_at_cap"] == 5
    assert stats["n_refseq_species_dropped_at_cap"] == 0


def test_dropped_refseq_species_counted_when_more_refseqs_than_target():
    # 110 species, 105 with RefSeqs, target 100 → 10 dropped, 5 of which
    # had a RefSeq.  Surfaces in stats["n_refseq_species_dropped_at_cap"].
    species = {}
    refseq_ids = set()
    for i in range(105):
        sp = f"REF{i:03d}"
        species[sp] = _reps(sp, 1)
        refseq_ids.add(species[sp][0].id)
    for i in range(5):
        sp = f"NRF{i}"
        species[sp] = _reps(sp, 1)
    out, stats = proportional_merge(species, 100, priority_ids=refseq_ids)
    assert len(out) == 100
    assert stats["n_species_dropped_at_cap"] == 10
    # 5 RefSeq species dropped (alphabetical tie-break: REF100..REF104)
    assert stats["n_refseq_species_dropped_at_cap"] == 5


def test_warning_emitted_when_species_dropped_at_cap(caplog):
    import logging
    species = {f"SP{i:02d}": _reps(f"SP{i:02d}", 1) for i in range(10)}
    with caplog.at_level(logging.WARNING, logger="vfam_trees.subsample"):
        proportional_merge(species, 5)
    assert any("species exceeds target" in r.getMessage() for r in caplog.records)


def test_no_warning_when_species_fit_within_target(caplog):
    import logging
    species = {f"SP{i:02d}": _reps(f"SP{i:02d}", 1) for i in range(3)}
    with caplog.at_level(logging.WARNING, logger="vfam_trees.subsample"):
        out, stats = proportional_merge(species, 100)
    assert stats["n_species_dropped_at_cap"] == 0
    assert not any("species exceeds target" in r.getMessage() for r in caplog.records)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def test_same_seed_gives_same_result():
    species = {f"SP{i}": _reps(f"SP{i}", 5) for i in range(6)}
    a, _ = proportional_merge(species, 10, seed=123)
    b, _ = proportional_merge(species, 10, seed=123)
    assert [r.id for r in a] == [r.id for r in b]


def test_returns_tuple_of_list_and_dict():
    """API shape: ensure callers always unpack (records, stats)."""
    out = proportional_merge({"A": _reps("A", 2)}, 10)
    assert isinstance(out, tuple)
    assert len(out) == 2
    records, stats = out
    assert isinstance(records, list)
    assert isinstance(stats, dict)
    for key in ("n_species_dropped_at_cap", "n_refseq_species_dropped_at_cap",
                "n_priority_kept", "n_priority_total"):
        assert key in stats


# ---------------------------------------------------------------------------
# absorb_into_refseqs
# ---------------------------------------------------------------------------

def _mock_membership(monkeypatch, clusters: list[set[str]]):
    """Patch _cluster_membership to return *clusters* verbatim."""
    def _fake(records, threshold, seq_type, work_dir, clustering_tool, *args, **kwargs):
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


# ---------------------------------------------------------------------------
# Threads plumbing — regression for v1.2.11.  Earlier code hard-coded
# `--threads 1` in the MMseqs2 invocation and `-T 1` in the CD-HIT
# invocation, so MAFFT/IQ-TREE got the requested threads but clustering
# always ran single-threaded.
# ---------------------------------------------------------------------------

class _CapturedRun:
    """subprocess.run stand-in: captures argv, simulates a clean exit, and
    writes the minimum output files the wrappers expect."""

    def __init__(self):
        self.cmds: list[list[str]] = []

    def __call__(self, cmd, *args, **kwargs):
        self.cmds.append(list(cmd))
        # Synthesize the output the parser expects so the wrapper does not
        # raise a FileNotFoundError after subprocess.run "completes".
        if cmd[0] == "mmseqs":
            out_prefix = Path(cmd[3])  # easy-linclust positional 2 → output prefix
            (out_prefix.parent).mkdir(parents=True, exist_ok=True)
            (Path(str(out_prefix) + "_rep_seq.fasta")).write_text(">A\nATGC\n")
            (out_prefix.parent / "output_cluster.tsv").write_text("A\tA\n")
        elif cmd[0] in ("cd-hit", "cd-hit-est"):
            i = cmd.index("-o")
            out = Path(cmd[i + 1])
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(">A\nATGC\n")
            (out.parent / "output.fasta.clstr").write_text(">Cluster 0\n0\t4nt, >A...\n")

        class _R:
            returncode = 0
            stdout = ""
            stderr = ""
        return _R()


def _arg_after(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


class TestClusteringThreadsPlumbing:
    def test_mmseqs2_receives_threads_flag(self, tmp_path, monkeypatch):
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        subsample._mmseqs2_cluster(
            [_rec("A")], 0.99, "nucleotide", tmp_path, threads=8,
        )
        assert run.cmds, "mmseqs invocation never happened"
        assert _arg_after(run.cmds[0], "--threads") == "8"

    def test_mmseqs2_default_threads_is_one(self, tmp_path, monkeypatch):
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        subsample._mmseqs2_cluster([_rec("A")], 0.99, "nucleotide", tmp_path)
        assert _arg_after(run.cmds[0], "--threads") == "1"

    def test_cdhit_receives_threads_flag_protein(self, tmp_path, monkeypatch):
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        subsample._cdhit_cluster(
            [SeqRecord(Seq("MMMM"), id="A")], 0.9, "protein", tmp_path, threads=4,
        )
        assert run.cmds[0][0] == "cd-hit"
        assert _arg_after(run.cmds[0], "-T") == "4"

    def test_cdhit_receives_threads_flag_nucleotide(self, tmp_path, monkeypatch):
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        subsample._cdhit_cluster(
            [_rec("A")], 0.95, "nucleotide", tmp_path, threads=6,
        )
        assert run.cmds[0][0] == "cd-hit-est"
        assert _arg_after(run.cmds[0], "-T") == "6"

    def test_threads_flow_through_adaptive_cluster_species(
        self, tmp_path, monkeypatch
    ):
        """End-to-end: adaptive_cluster_species(threads=N) must reach the
        actual mmseqs argv."""
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        records = [SeqRecord(Seq("ATGC" * 20), id=f"r{i}") for i in range(10)]
        subsample.adaptive_cluster_species(
            records=records,
            max_reps=3,
            threshold_min=0.7,
            threshold_max=0.99,
            seq_type="nucleotide",
            work_dir=tmp_path,
            clustering_tool="mmseqs2",
            threads=8,
        )
        assert run.cmds, "mmseqs was never invoked"
        # Every mmseqs invocation in the binary search must use --threads 8.
        for cmd in run.cmds:
            assert _arg_after(cmd, "--threads") == "8", cmd

    def test_threads_flow_through_absorb_into_refseqs(
        self, tmp_path, monkeypatch
    ):
        run = _CapturedRun()
        monkeypatch.setattr(subsample.subprocess, "run", run)
        records = [_rec("NC_1"), _rec("X1"), _rec("X2")]
        absorb_into_refseqs(
            records=records,
            refseq_ids={"NC_1"},
            threshold=0.99,
            seq_type="nucleotide",
            work_dir=tmp_path,
            clustering_tool="mmseqs2",
            threads=4,
        )
        assert run.cmds, "mmseqs was never invoked"
        assert _arg_after(run.cmds[0], "--threads") == "4"


# ---------------------------------------------------------------------------
# _cluster_at: protected_ids force-inclusion (regression for v1.2.31)
#
# Bug: adaptive_cluster_species delegates representative selection to the
# external clustering tool (mmseqs2/cd-hit).  That tool has no concept of
# "protected" sequences, so a manual.include record could be assigned to a
# cluster whose representative was a different sequence — silently dropping
# it before proportional_merge ever sees it.  This happened only for the
# 100-tree (smaller max_reps → more aggressive threshold → more collisions)
# while the 500-tree was unaffected.
# ---------------------------------------------------------------------------

class TestClusterAtProtectedIds:
    def test_protected_id_force_added_when_tool_omits_it(self, tmp_path, monkeypatch):
        """Core regression: protected record absent from tool output is added back."""
        records = [_rec("A"), _rec("MANUAL_INCLUDE"), _rec("C")]
        monkeypatch.setattr(subsample, "_mmseqs2_cluster",
                            lambda *a, **k: ["A", "C"])
        result = subsample._cluster_at(
            records, 0.8, "nucleotide", tmp_path, "mmseqs2",
            protected_ids={"MANUAL_INCLUDE"},
        )
        assert {r.id for r in result} == {"A", "C", "MANUAL_INCLUDE"}

    def test_protected_id_not_duplicated_when_already_chosen(self, tmp_path, monkeypatch):
        records = [_rec("A"), _rec("MANUAL_INCLUDE")]
        monkeypatch.setattr(subsample, "_mmseqs2_cluster",
                            lambda *a, **k: ["A", "MANUAL_INCLUDE"])
        result = subsample._cluster_at(
            records, 0.8, "nucleotide", tmp_path, "mmseqs2",
            protected_ids={"MANUAL_INCLUDE"},
        )
        assert [r.id for r in result].count("MANUAL_INCLUDE") == 1

    def test_no_protected_ids_result_unchanged(self, tmp_path, monkeypatch):
        """Backwards compat: omitting protected_ids behaves exactly as before."""
        records = [_rec("A"), _rec("B"), _rec("C")]
        monkeypatch.setattr(subsample, "_mmseqs2_cluster",
                            lambda *a, **k: ["A"])
        result = subsample._cluster_at(
            records, 0.8, "nucleotide", tmp_path, "mmseqs2",
        )
        assert [r.id for r in result] == ["A"]

    def test_multiple_protected_ids_all_forced_in(self, tmp_path, monkeypatch):
        records = [_rec("A"), _rec("P1"), _rec("P2"), _rec("B")]
        monkeypatch.setattr(subsample, "_mmseqs2_cluster",
                            lambda *a, **k: ["A"])
        result = subsample._cluster_at(
            records, 0.8, "nucleotide", tmp_path, "mmseqs2",
            protected_ids={"P1", "P2"},
        )
        assert {r.id for r in result} >= {"A", "P1", "P2"}

    def test_protected_id_not_in_input_is_silently_ignored(self, tmp_path, monkeypatch):
        """A protected ID that isn't even in records cannot be injected."""
        records = [_rec("A"), _rec("B")]
        monkeypatch.setattr(subsample, "_mmseqs2_cluster",
                            lambda *a, **k: ["A"])
        result = subsample._cluster_at(
            records, 0.8, "nucleotide", tmp_path, "mmseqs2",
            protected_ids={"GHOST"},
        )
        assert {r.id for r in result} == {"A"}


# ---------------------------------------------------------------------------
# adaptive_cluster_species: protected_ids survive aggressive clustering
# ---------------------------------------------------------------------------

class TestAdaptiveClusterSpeciesProtectedIds:
    def _records(self, n: int, protected_id: str = "MANUAL_INCLUDE") -> list[SeqRecord]:
        recs = [SeqRecord(Seq("ATGCATGC"), id=f"r{i}") for i in range(n)]
        recs.append(SeqRecord(Seq("ATGCATGC"), id=protected_id))
        return recs

    def test_protected_id_survives_when_tool_omits_it(self, tmp_path, monkeypatch):
        """Regression: manual.include seq present in 500-tree but missing from
        100-tree because smaller max_reps triggers tighter clustering where the
        tool picks a different cluster representative."""
        records = self._records(8)

        # Simulate clustering tool that never picks MANUAL_INCLUDE as a rep.
        monkeypatch.setattr(
            subsample, "_mmseqs2_cluster",
            lambda recs, *a, **k: [r.id for r in recs if r.id != "MANUAL_INCLUDE"][:3],
        )

        result, _ = subsample.adaptive_cluster_species(
            records=records,
            max_reps=3,
            threshold_min=0.70,
            threshold_max=0.99,
            seq_type="nucleotide",
            work_dir=tmp_path,
            protected_ids={"MANUAL_INCLUDE"},
        )
        assert "MANUAL_INCLUDE" in {r.id for r in result}, (
            "manual.include sequence was dropped by adaptive_cluster_species — "
            "regression of v1.2.31 bug"
        )

    def test_without_protected_ids_manual_include_can_be_dropped(
        self, tmp_path, monkeypatch
    ):
        """Control: without protected_ids, the same clustering tool omitting
        MANUAL_INCLUDE causes it to be absent — confirming the pre-fix behaviour
        that motivated the fix."""
        records = self._records(8)

        monkeypatch.setattr(
            subsample, "_mmseqs2_cluster",
            lambda recs, *a, **k: [r.id for r in recs if r.id != "MANUAL_INCLUDE"][:3],
        )

        result, _ = subsample.adaptive_cluster_species(
            records=records,
            max_reps=3,
            threshold_min=0.70,
            threshold_max=0.99,
            seq_type="nucleotide",
            work_dir=tmp_path,
            # no protected_ids → pre-fix behaviour
        )
        assert "MANUAL_INCLUDE" not in {r.id for r in result}

    def test_protected_ids_forwarded_to_every_cluster_at_call(
        self, tmp_path, monkeypatch
    ):
        """Plumbing: every _cluster_at call inside binary search receives protected_ids."""
        received: list[set | None] = []

        real_cluster_at = subsample._cluster_at

        def _spy(recs, threshold, seq_type, work_dir, tool, threads=1, protected_ids=None):
            received.append(protected_ids)
            return real_cluster_at(
                recs, threshold, seq_type, work_dir, tool, threads,
                protected_ids=protected_ids,
            )

        monkeypatch.setattr(subsample, "_cluster_at", _spy)
        monkeypatch.setattr(
            subsample, "_mmseqs2_cluster",
            lambda recs, *a, **k: [r.id for r in recs][:3],
        )

        records = [SeqRecord(Seq("ATGCATGC"), id=f"r{i}") for i in range(8)]
        subsample.adaptive_cluster_species(
            records=records,
            max_reps=3,
            threshold_min=0.70,
            threshold_max=0.99,
            seq_type="nucleotide",
            work_dir=tmp_path,
            protected_ids={"MANUAL_INCLUDE"},
        )
        assert received, "_cluster_at was never called"
        assert all(p == {"MANUAL_INCLUDE"} for p in received), (
            f"Some _cluster_at calls did not receive protected_ids: {received}"
        )
