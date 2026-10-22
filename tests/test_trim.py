"""Tests for vfam_trees.trim (trimAl wrapper)."""
import shutil
import subprocess
from pathlib import Path

import pytest

from vfam_trees import trim
from vfam_trees.config import DEFAULT_FAMILY_CONFIG
from vfam_trees.trim import DEFAULT_TRIM_OPTIONS, run_trim


def test_family_default_is_the_trim_default():
    assert DEFAULT_FAMILY_CONFIG["msa_trim"]["options"] == DEFAULT_TRIM_OPTIONS
    assert DEFAULT_TRIM_OPTIONS == "-gt 0.5"


def test_options_reach_trimal_as_separate_arguments(tmp_path, monkeypatch):
    """'-gt 0.5' as one argv entry makes trimAl fail; it must be two."""
    seen = {}

    def fake_run(cmd, capture_output, text):
        seen["cmd"] = cmd
        Path(cmd[cmd.index("-out") + 1]).write_text(">a\nAC\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(trim.subprocess, "run", fake_run)
    run_trim(tmp_path / "in.fasta", tmp_path / "out.fasta")
    assert seen["cmd"][-2:] == ["-gt", "0.5"]


def test_empty_output_names_the_setting(tmp_path, monkeypatch):
    def fake_run(cmd, capture_output, text):
        return subprocess.CompletedProcess(
            cmd, 0, "",
            "WARNING: Output alignment has not been generated. It is empty.\n")

    monkeypatch.setattr(trim.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="kept no column.*msa_trim"):
        run_trim(tmp_path / "in.fasta", tmp_path / "out.fasta")


def test_unexplained_missing_output_still_reported(tmp_path, monkeypatch):
    def fake_run(cmd, capture_output, text):
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(trim.subprocess, "run", fake_run)
    with pytest.raises(FileNotFoundError):
        run_trim(tmp_path / "in.fasta", tmp_path / "out.fasta")


@pytest.mark.skipif(shutil.which("trimal") is None, reason="trimAl not installed")
class TestRealTrimal:
    def test_default_drops_columns_with_more_than_half_gaps(self, tmp_path):
        src = tmp_path / "in.fasta"
        # column 5 is a gap in all 3 rows, column 4 in 2 of 3 → both go;
        # column 3 is a gap in 1 of 3 → kept.
        src.write_text(">a\nACGT-A\n>b\nACG--A\n>c\nAC---A\n")
        out = tmp_path / "out.fasta"
        run_trim(src, out)
        seqs = [l for l in out.read_text().splitlines() if not l.startswith(">")]
        assert seqs == ["ACGA", "ACGA", "AC-A"]

    def test_nothing_left_raises_clear_error(self, tmp_path):
        src = tmp_path / "in.fasta"
        src.write_text(">a\nAC----\n>b\n--GT--\n>c\n----AC\n")
        with pytest.raises(RuntimeError, match="kept no column"):
            run_trim(src, tmp_path / "out.fasta")
