"""Tests for vfam_trees.signal (phylogenetic-signal diagnostics)."""
import logging

import pytest

from vfam_trees import signal
from vfam_trees.signal import (
    alignment_identity,
    core_fraction,
    log_signal,
    measure,
    signal_verdict,
    summary_fields,
    tree_distances,
)


def _fasta(path, recs):
    path.write_text("".join(f">{i}\n{s}\n" for i, s in recs))
    return path


class TestAlignmentIdentity:
    def test_identical_rows(self, tmp_path):
        aln = _fasta(tmp_path / "a.fa", [("a", "ACGTACGTACGT"), ("b", "ACGTACGTACGT")])
        assert alignment_identity(aln) == {"median": 1.0, "p10": 1.0}

    def test_gaps_and_n_are_ignored(self, tmp_path):
        # 12 shared columns after dropping gaps / N; 3 mismatches → 0.75
        aln = _fasta(tmp_path / "a.fa", [("a", "ACGTACGTACGT--NN"),
                                         ("b", "ACGTACGTATTA" + "AAAA")])
        assert alignment_identity(aln)["median"] == pytest.approx(9 / 12)

    def test_protein_asparagine_counts(self, tmp_path):
        # N is asparagine in a protein: 13 of 23 shared residues identical.
        aln = _fasta(tmp_path / "a.fa", [("a", "NNNNNNNNNNKLMAAAAAAAAAA"),
                                         ("b", "NNNNNNNNNNKLMCCCCCCCCCC")])
        assert alignment_identity(aln, "protein")["median"] == pytest.approx(13 / 23)
        # Read as nucleotides, only the 10 A/C columns are residues.
        assert alignment_identity(aln, "nucleotide")["median"] == pytest.approx(0.0)

    def test_protein_x_and_ambiguity_codes_ignored(self, tmp_path):
        aln = _fasta(tmp_path / "a.fa", [("a", "MKLVXBZAAAAAAAAAA"),
                                         ("b", "MKLVAAAAAAAAAAAAA")])
        assert alignment_identity(aln, "protein")["median"] == pytest.approx(1.0)

    def test_nucleotide_iupac_codes_ignored(self, tmp_path):
        # R / Y are ambiguity codes, not mismatching residues.
        aln = _fasta(tmp_path / "a.fa", [("a", "ACGTACGTACGTRY"), ("b", "ACGTACGTACGTAC")])
        assert alignment_identity(aln, "nucleotide")["median"] == pytest.approx(1.0)

    def test_measure_passes_the_sequence_type(self, tmp_path, monkeypatch):
        seen = []
        monkeypatch.setattr(signal, "alignment_identity",
                            lambda p, t="nucleotide": seen.append(t) or {"median": None, "p10": None})
        monkeypatch.setattr(signal, "tree_distances", lambda p: {
            "patristic_median": None, "patristic_p90": None, "n_branches_at_cap": None})
        signal.measure(tmp_path / "a.fa", None, tmp_path / "t.nwk", "protein")
        assert seen == ["protein"]

    def test_pairs_with_too_little_overlap_are_skipped(self, tmp_path):
        aln = _fasta(tmp_path / "a.fa", [("a", "ACGTA-----"), ("b", "-----ACGTA")])
        assert alignment_identity(aln) == {"median": None, "p10": None}

    def test_sampling_is_deterministic(self, tmp_path):
        import random
        rng = random.Random(1)
        recs = [(f"s{i}", "".join(rng.choice("ACGT") for _ in range(60))) for i in range(80)]
        aln = _fasta(tmp_path / "a.fa", recs)
        assert alignment_identity(aln) == alignment_identity(aln)
        assert len(signal._sample_pairs(80)) == signal.N_PAIRS   # 3160 pairs > N_PAIRS


class TestCoreFraction:
    def test_per_sequence_by_id(self, tmp_path):
        raw = _fasta(tmp_path / "raw.fa", [("a", "A" * 100), ("b", "C" * 200), ("c", "G" * 100)])
        aln = _fasta(tmp_path / "aln.fa", [("a", "A" * 50 + "--"), ("b", "C" * 50 + "--"),
                                           ("c", "G" * 25 + "-" * 27)])
        # kept/full: 0.5, 0.25, 0.25 → median 0.25
        assert core_fraction(aln, raw) == pytest.approx(0.25)

    def test_falls_back_to_medians_when_ids_differ(self, tmp_path):
        raw = _fasta(tmp_path / "raw.fa", [("x", "A" * 100), ("y", "A" * 100)])
        aln = _fasta(tmp_path / "aln.fa", [("a", "A" * 40), ("b", "A" * 40)])
        assert core_fraction(aln, raw) == pytest.approx(0.4)

    def test_gappy_long_alignment_is_not_above_one(self, tmp_path):
        # A trimmed alignment far longer than the sequences (mostly gaps)
        # used to give a "core" above 1; per-sequence residues cannot.
        raw = _fasta(tmp_path / "raw.fa", [("a", "A" * 100), ("b", "C" * 100)])
        aln = _fasta(tmp_path / "aln.fa", [("a", "A" * 90 + "-" * 500), ("b", "-" * 500 + "C" * 90)])
        assert core_fraction(aln, raw) == pytest.approx(0.9)


class TestTreeDistances:
    def test_patristic_and_cap(self, tmp_path):
        nwk = tmp_path / "t.nwk"
        # A–B = 0.3, C–D = 0.7, cross pairs = 0.1+0.2+0.05+0.3 / … ; cap branch on E
        nwk.write_text("(((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.05):0.0,E:10.0);")
        d = tree_distances(nwk)
        assert d["n_branches_at_cap"] == 1
        # 10 pairs; the 4 pairs with E are > 10, the median of all 10 is the
        # mean of the 5th and 6th smallest distances.
        assert d["patristic_median"] > 0.3
        assert d["patristic_p90"] > 10.0

    def test_unparseable_tree(self, tmp_path):
        nwk = tmp_path / "t.nwk"
        nwk.write_text("not a tree((")
        d = tree_distances(nwk)
        assert d["patristic_median"] is None

    def test_matches_biopython_distance(self, tmp_path):
        from Bio import Phylo
        nwk = tmp_path / "t.nwk"
        nwk.write_text("((A:0.1,B:0.2):0.3,((C:0.05,D:0.15):0.2,E:0.6):0.1);")
        tree = Phylo.read(str(nwk), "newick")
        tips = tree.get_terminals()
        import statistics
        expected = statistics.median(
            tree.distance(a, b) for i, a in enumerate(tips) for b in tips[i + 1:])
        assert tree_distances(nwk)["patristic_median"] == pytest.approx(expected)


class TestVerdict:
    def test_nucleotide_levels(self):
        assert signal_verdict("nucleotide", {"patristic_median": 2.0})[0] == "saturated"
        assert signal_verdict("nucleotide", {"patristic_median": 1.0})[0] == "weak"
        assert signal_verdict("nucleotide", {"patristic_median": 0.4, "identity_median": 0.8,
                                             "core_frac": 0.9})[0] == "ok"

    def test_cap_with_weak_median_is_saturated(self):
        v, reasons = signal_verdict("nucleotide", {"patristic_median": 1.0, "n_branches_at_cap": 2})
        assert v == "saturated" and any("cap" in r for r in reasons)

    def test_cap_alone_is_weak(self):
        assert signal_verdict("nucleotide", {"patristic_median": 0.3, "n_branches_at_cap": 1})[0] == "weak"

    def test_low_identity_or_core_is_weak(self):
        assert signal_verdict("nucleotide", {"patristic_median": 0.3, "identity_median": 0.5})[0] == "weak"
        assert signal_verdict("nucleotide", {"patristic_median": 0.3, "core_frac": 0.1})[0] == "weak"

    def test_protein_thresholds_are_higher(self):
        assert signal_verdict("protein", {"patristic_median": 2.0, "identity_median": 0.5})[0] == "ok"
        assert signal_verdict("protein", {"patristic_median": 3.5})[0] == "weak"
        assert signal_verdict("protein", {"patristic_median": 4.5})[0] == "saturated"

    def test_nothing_measured(self):
        assert signal_verdict("nucleotide", {}) == ("", [])


class TestReporting:
    def test_summary_fields_round_and_blank(self):
        f = summary_fields({"verdict": "weak", "patristic_median": 1.234567,
                            "identity_median": None, "core_frac": 0.5, "n_branches_at_cap": 0})
        assert f == {"signal": "weak", "patristic_median": 1.235, "identity_median": "",
                     "core_frac": 0.5, "branches_at_cap": 0}
        assert summary_fields({})["signal"] == ""

    def test_log_warns_only_when_not_ok(self, caplog):
        log = logging.getLogger("test_signal")
        with caplog.at_level(logging.INFO, logger="test_signal"):
            log_signal(log, "100", "nucleotide", {"verdict": "ok", "patristic_median": 0.2})
            assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
            log_signal(log, "100", "nucleotide",
                       {"verdict": "saturated", "reasons": ["median tip-to-tip distance 2.00"]})
        warn = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warn) == 1 and "SATURATED" in warn[0].getMessage()
        assert "protein marker" in warn[0].getMessage()


def test_measure_end_to_end(tmp_path):
    raw = _fasta(tmp_path / "raw.fa", [(c, "ACGTACGTACGTACGTACGT") for c in "ABCD"])
    aln = _fasta(tmp_path / "aln.fa", [(c, "ACGTACGTACGTACGTACGT") for c in "ABCD"])
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((A:0.01,B:0.01):0.01,(C:0.01,D:0.01):0.01);")
    m = measure(aln, raw, nwk, "nucleotide")
    assert m["verdict"] == "ok"
    assert m["core_frac"] == pytest.approx(1.0) and m["identity_median"] == 1.0
    assert measure(aln, None, nwk, "protein")["core_frac"] is None
