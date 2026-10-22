"""Tests for vfam_trees.msa helpers (no MAFFT calls; MAFFT is scripted)."""
from pathlib import Path

from Bio import SeqIO

import vfam_trees.msa as msa_mod
from vfam_trees.msa import (
    MAFFT_REVERSED_PREFIX,
    read_flipped_ids,
    reverse_complement,
    run_msa,
    strip_mafft_direction_prefix,
)


def _ids(path: Path) -> list[str]:
    return [r.id for r in SeqIO.parse(str(path), "fasta")]


def _seqs(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq) for r in SeqIO.parse(str(path), "fasta")}


class TestReverseComplement:
    def test_iupac_case_and_gaps(self):
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("aacg-t") == "a-cgtt"
        assert reverse_complement("RYKMBVDHNSW") == "WSNDHBVKMRY"
        assert reverse_complement("--AC.") == ".GT--"

    def test_involution(self):
        s = "acgtRYKMbvdhNSW--n"
        assert reverse_complement(reverse_complement(s)) == s


class TestStripMafftDirectionPrefix:
    def test_strips_prefix_and_reports_ids(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(
            ">S000001\nACGT--\n"
            ">_R_S000002 some description\nAC--GT\n"
            ">S000003\n--ACGT\n"
            ">_R_S000004\nACGTAC\n"
        )
        flipped = strip_mafft_direction_prefix(aln)
        assert flipped == ["S000002", "S000004"]
        text = aln.read_text()
        assert MAFFT_REVERSED_PREFIX not in text
        assert ">S000002 some description\n" in text
        assert ">S000004\n" in text
        # Sequence rows are left exactly as MAFFT wrote them.
        assert _seqs(aln)["S000002"] == "AC--GT"
        assert _seqs(aln)["S000004"] == "ACGTAC"

    def test_no_prefix_leaves_file_untouched(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        original = ">S000001\nACGT\n>S000002\nACGA\n"
        aln.write_text(original)
        mtime = aln.stat().st_mtime_ns
        assert strip_mafft_direction_prefix(aln) == []
        assert aln.read_text() == original
        assert aln.stat().st_mtime_ns == mtime

    def test_prefix_inside_name_is_not_a_marker(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">S_R_001\nACGT\n")
        assert strip_mafft_direction_prefix(aln) == []
        assert aln.read_text() == ">S_R_001\nACGT\n"


# A strand-specific test sequence and its reverse complement.
FWD = "AACCGGTTTAGC"
REV = reverse_complement(FWD)


class _FakeMafft:
    """Plays MAFFT --adjustdirection: rows on the other strand from the first
    record are reverse-complemented and marked with _R_, as MAFFT does."""

    def __init__(self):
        self.calls = 0

    def __call__(self, input_fasta, output_fasta, options, threads):
        self.calls += 1
        recs = list(SeqIO.parse(str(input_fasta), "fasta"))
        lead = str(recs[0].seq)
        with open(output_fasta, "w") as fh:
            for r in recs:
                s = str(r.seq)
                if s != lead:
                    fh.write(f">{MAFFT_REVERSED_PREFIX}{r.id}\n{reverse_complement(s)}\n")
                else:
                    fh.write(f">{r.id}\n{s}\n")
        if "--adjustdirection" not in options:
            return []
        return strip_mafft_direction_prefix(output_fasta)


def _write_input(path: Path, strands: dict[str, str]) -> None:
    path.write_text("".join(
        f">{i}\n{FWD if s == 'F' else REV}\n" for i, s in strands.items()
    ))


def _run(tmp_path, monkeypatch, strands):
    fake = _FakeMafft()
    monkeypatch.setattr(msa_mod, "_run_mafft", fake)
    inp = tmp_path / "raw.fasta"
    _write_input(inp, strands)
    out = tmp_path / "aln.fasta"
    flipped = run_msa(inp, out, tool="mafft", options="--retree 2 --adjustdirection")
    return fake, out, flipped


class TestMajorityStrand:
    def test_minority_leader_is_fixed_without_second_mafft_run(self, tmp_path, monkeypatch):
        fake, out, flipped = _run(tmp_path, monkeypatch,
                                  {"S1": "R", "S2": "F", "S3": "F", "S4": "R", "S5": "F"})
        assert fake.calls == 1
        assert flipped == ["S1", "S4"]                # the true minority
        assert _ids(out) == ["S1", "S2", "S3", "S4", "S5"]
        assert set(_seqs(out).values()) == {FWD}      # all on the majority strand
        assert not list(tmp_path.glob("*reoriented*"))

    def test_majority_leader_kept_as_is(self, tmp_path, monkeypatch):
        _, out, flipped = _run(tmp_path, monkeypatch, {"S1": "F", "S2": "F", "S3": "R"})
        assert flipped == ["S3"]
        assert set(_seqs(out).values()) == {FWD}

    def test_tie_goes_to_strand_of_smallest_id_whatever_the_order(self, tmp_path, monkeypatch):
        results = []
        for order in (["S3", "S1", "S4", "S2"], ["S1", "S2", "S3", "S4"], ["S4", "S2", "S1", "S3"]):
            strands = {i: ("F" if i in ("S1", "S2") else "R") for i in order}
            sub = tmp_path / "_".join(order)
            sub.mkdir()
            _, out, flipped = _run(sub, monkeypatch, strands)
            results.append((sorted(flipped), _seqs(out)["S1"]))
        assert all(r == (["S3", "S4"], FWD) for r in results)

    def test_no_adjustdirection_no_flips(self, tmp_path, monkeypatch):
        fake = _FakeMafft()
        monkeypatch.setattr(msa_mod, "_run_mafft", fake)
        inp = tmp_path / "raw.fasta"
        _write_input(inp, {"S1": "F", "S2": "F"})
        assert run_msa(inp, tmp_path / "aln.fasta", tool="mafft", options="--retree 2") == []
        assert fake.calls == 1


class TestFlippedRecord:
    def test_recorded_for_a_cached_alignment_and_cleared_on_rerun(self, tmp_path, monkeypatch):
        _, out, flipped = _run(tmp_path, monkeypatch, {"S1": "F", "S2": "R", "S3": "F"})
        assert read_flipped_ids(out) == flipped == ["S2"]
        # A later run of the same alignment path with nothing flipped clears it.
        _write_input(tmp_path / "raw.fasta", {"S1": "F", "S2": "F", "S3": "F"})
        run_msa(tmp_path / "raw.fasta", out, tool="mafft", options="--adjustdirection")
        assert read_flipped_ids(out) == []

    def test_missing_record_is_empty(self, tmp_path):
        assert read_flipped_ids(tmp_path / "nothing.fasta") == []


class TestProteinInput:
    """run_msa(seq_type="protein") aligns a copy with residues MAFFT rejects
    written as X; the caller's FASTA is left as deposited."""

    def _capture(self, monkeypatch):
        seen = {}

        def fake(input_fasta, output_fasta, options, threads):
            seen["input"] = input_fasta
            seen["seqs"] = _seqs(input_fasta)
            output_fasta.write_text(input_fasta.read_text())
            return []
        monkeypatch.setattr(msa_mod, "_run_mafft", fake)
        return seen

    def test_rejected_residues_written_as_x_in_a_copy(self, tmp_path, monkeypatch, caplog):
        import logging
        seen = self._capture(monkeypatch)
        inp = tmp_path / "raw.fasta"
        inp.write_text(">p1\nMKOLAUV\n>p2\nMKLAV*\n")
        with caplog.at_level(logging.INFO, logger="vfam_trees.msa"):
            run_msa(inp, tmp_path / "aln.fasta", options="--auto", seq_type="protein")
        assert seen["input"] != inp
        assert seen["seqs"] == {"p1": "MKXLAXV", "p2": "MKLAV"}
        assert inp.read_text() == ">p1\nMKOLAUV\n>p2\nMKLAV*\n"     # untouched
        assert "1 protein sequence(s)" in caplog.text and "p1" in caplog.text

    def test_clean_input_aligned_as_is(self, tmp_path, monkeypatch, caplog):
        import logging
        seen = self._capture(monkeypatch)
        inp = tmp_path / "raw.fasta"
        inp.write_text(">p1\nMKLAV\n>p2\nMKLAV\n")
        with caplog.at_level(logging.INFO, logger="vfam_trees.msa"):
            run_msa(inp, tmp_path / "aln.fasta", options="--auto", seq_type="protein")
        assert seen["input"] == inp
        assert "cannot align" not in caplog.text

    def test_nucleotide_input_never_touched(self, tmp_path, monkeypatch):
        seen = self._capture(monkeypatch)
        inp = tmp_path / "raw.fasta"
        inp.write_text(">n1\nACGTU\n")
        run_msa(inp, tmp_path / "aln.fasta", options="--auto")
        assert seen["input"] == inp


class TestMafftThreads:
    """Threaded iterative refinement is nondeterministic: --threadit 0 by default."""

    def _cmd(self, tmp_path, monkeypatch, options, threads=8):
        import subprocess
        seen = {}

        def fake_run(cmd, **kw):
            seen["cmd"] = cmd
            out = Path(cmd[cmd.index("--out") + 1])
            out.write_text(">a\nMKV\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")
        monkeypatch.setattr(msa_mod.subprocess, "run", fake_run)
        inp = tmp_path / "in.fasta"
        inp.write_text(">a\nMKV\n")
        msa_mod._run_mafft(inp, tmp_path / "out.fasta", options, threads)
        return seen["cmd"]

    def test_refinement_on_one_thread_by_default(self, tmp_path, monkeypatch):
        cmd = self._cmd(tmp_path, monkeypatch, "--maxiterate 1000 --localpair")
        i = cmd.index("--thread")
        assert cmd[i:i + 4] == ["--thread", "8", "--threadit", "0"]

    def test_configured_threadit_wins(self, tmp_path, monkeypatch):
        cmd = self._cmd(tmp_path, monkeypatch, "--auto --threadit 4")
        assert cmd.count("--threadit") == 1
        assert cmd[cmd.index("--threadit") + 1] == "4"

    def test_thread_options_helper(self):
        assert msa_mod.mafft_thread_options("--retree 2", 1) == ["--thread", "1", "--threadit", "0"]
        assert msa_mod.mafft_thread_options("--threadit=2 --auto", 3) == ["--thread", "3"]
