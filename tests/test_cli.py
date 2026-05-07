"""Tests for vfam_trees.cli — argument plumbing into Snakemake.

Covers the v1.2.10 fix where `--cores` is now passed as `cores * threads`
so Snakemake does not silently clamp the per-rule `threads:` declaration
down to `--cores` (which previously made `-j 1 -t 8` run single-threaded).
"""
from pathlib import Path

import pytest
from click.testing import CliRunner

import vfam_trees.cli as cli_mod
from vfam_trees.cli import main


def _families_file(tmp_path: Path) -> Path:
    f = tmp_path / "families.txt"
    f.write_text("Poxviridae\nFiloviridae\n")
    return f


class _FakeCompleted:
    """Minimal stand-in for subprocess.CompletedProcess."""

    def __init__(self, returncode: int = 0):
        self.returncode = returncode


@pytest.fixture
def captured_cmd(monkeypatch):
    """Patch subprocess.run + sys.exit so the CLI's snakemake invocation
    is intercepted and the constructed argv can be inspected."""
    captured: dict[str, list] = {}

    def _fake_run(cmd, *args, **kwargs):
        captured["cmd"] = list(cmd)
        return _FakeCompleted(returncode=0)

    monkeypatch.setattr(cli_mod.subprocess, "run", _fake_run)
    # Suppress the post-run overview PNG step (returncode==0 path triggers it).
    import vfam_trees.report as report_mod
    monkeypatch.setattr(report_mod, "generate_overview_png",
                        lambda *a, **k: None)
    return captured


def _cores_arg(cmd: list[str]) -> int:
    """Extract the integer that follows --cores in the constructed argv."""
    i = cmd.index("--cores")
    return int(cmd[i + 1])


class TestSnakemakeCoresPlumbing:
    def test_j1_t8_passes_cores_8(self, tmp_path, captured_cmd):
        """`-j 1 -t 8` must yield --cores 8 so the run_family rule can
        actually use 8 threads (Snakemake clamps to --cores)."""
        runner = CliRunner()
        result = runner.invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"),
            "-j", "1", "-t", "8",
        ])
        assert result.exit_code == 0, result.output
        assert _cores_arg(captured_cmd["cmd"]) == 8

    def test_j4_t4_passes_cores_16(self, tmp_path, captured_cmd):
        runner = CliRunner()
        result = runner.invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"),
            "-j", "4", "-t", "4",
        ])
        assert result.exit_code == 0, result.output
        assert _cores_arg(captured_cmd["cmd"]) == 16

    def test_j2_t1_passes_cores_2(self, tmp_path, captured_cmd):
        # Single-threaded jobs but parallel families — total cores = 2.
        runner = CliRunner()
        result = runner.invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"),
            "-j", "2", "-t", "1",
        ])
        assert result.exit_code == 0, result.output
        assert _cores_arg(captured_cmd["cmd"]) == 2

    def test_threads_value_propagated_to_snakemake_config(
        self, tmp_path, captured_cmd
    ):
        """The per-job thread count must also reach the Snakefile via --config
        threads=N (the Snakefile uses this for `rule run_family: threads:`)."""
        runner = CliRunner()
        result = runner.invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"),
            "-j", "1", "-t", "8",
        ])
        assert result.exit_code == 0, result.output
        cmd = captured_cmd["cmd"]
        # Locate `threads=8` among the --config tokens.
        config_tokens = [
            tok for tok in cmd if isinstance(tok, str) and tok.startswith("threads=")
        ]
        assert config_tokens == ["threads=8"]

    def test_default_cores_and_threads_are_one(self, tmp_path, captured_cmd):
        # Defaults are -j 1 -t 1 → --cores 1 (regression: no accidental
        # multiplier on the default path).
        runner = CliRunner()
        result = runner.invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"),
        ])
        assert result.exit_code == 0, result.output
        assert _cores_arg(captured_cmd["cmd"]) == 1
