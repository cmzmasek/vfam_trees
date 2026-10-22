"""Tests for vfam_trees.cli — argument plumbing into Snakemake.

Covers the v1.2.10 fix where `--cores` is now passed as `cores * threads`
so Snakemake does not silently clamp the per-rule `threads:` declaration
down to `--cores` (which previously made `-j 1 -t 8` run single-threaded).
"""
import csv
import json
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
def captured_cmd(monkeypatch, tmp_path):
    """Patch subprocess.run + sys.exit so the CLI's snakemake invocation
    is intercepted and the constructed argv can be inspected.  Runs in
    tmp_path with a valid config/global.yaml (run checks it before
    Snakemake starts); tool probing for the run record is stubbed."""
    captured: dict[str, list] = {}
    monkeypatch.chdir(tmp_path)
    (tmp_path / "config").mkdir(exist_ok=True)
    (tmp_path / "config" / "global.yaml").write_text(
        'ncbi:\n  email: "test@example.org"\n  api_key: "SECRET123"\n'
    )
    import vfam_trees.provenance as prov_mod
    monkeypatch.setattr(prov_mod, "tool_versions", lambda: {"mafft": "7.0"})
    monkeypatch.setattr(cli_mod, "_missing_tools", lambda *a: [])

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


# ---------------------------------------------------------------------------
# _filter_tsv_rows + _clear_status — --force should drop stale rows from
# summary.tsv and status.tsv to prevent duplicate per-family rows on re-run.
# Regression for v1.2.12.
# ---------------------------------------------------------------------------

class TestFilterTsvRows:
    def _write_summary(self, path: Path, families: list[str]) -> None:
        with open(path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["family", "ncbi_taxid", "lineage"])
            for fam in families:
                writer.writerow([fam, "0", "Viruses"])

    def test_drops_matching_rows(self, tmp_path):
        from vfam_trees.cli import _filter_tsv_rows
        path = tmp_path / "summary.tsv"
        self._write_summary(path, ["Poxviridae", "Filoviridae", "Reoviridae"])
        n = _filter_tsv_rows(path, {"Poxviridae", "Reoviridae"})
        assert n == 2
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            kept = [row["family"] for row in reader]
        assert kept == ["Filoviridae"]

    def test_no_op_when_no_match(self, tmp_path):
        from vfam_trees.cli import _filter_tsv_rows
        path = tmp_path / "summary.tsv"
        self._write_summary(path, ["Filoviridae"])
        original = path.read_bytes()
        n = _filter_tsv_rows(path, {"Poxviridae"})
        assert n == 0
        # File unchanged when nothing matches (avoid spurious mtime update)
        assert path.read_bytes() == original

    def test_missing_file_no_op(self, tmp_path):
        from vfam_trees.cli import _filter_tsv_rows
        n = _filter_tsv_rows(tmp_path / "nope.tsv", {"Poxviridae"})
        assert n == 0

    def test_empty_file_no_op(self, tmp_path):
        from vfam_trees.cli import _filter_tsv_rows
        path = tmp_path / "summary.tsv"
        path.write_text("")
        assert _filter_tsv_rows(path, {"Poxviridae"}) == 0

    def test_clear_status_drops_summary_rows(self, tmp_path):
        # End-to-end: running --force on a family must remove its stale
        # summary.tsv row, otherwise the re-run produces duplicates.
        from vfam_trees.cli import _clear_status
        out = tmp_path / "out"
        out.mkdir()
        (out / ".done_Poxviridae").write_text("")
        self._write_summary(out / "summary.tsv", ["Poxviridae", "Filoviridae"])
        self._write_summary(out / "status.tsv", ["Poxviridae", "Filoviridae"])
        _clear_status(["Poxviridae"], out)
        # Sentinel removed
        assert not (out / ".done_Poxviridae").exists()
        # Summary + status: only Filoviridae left
        for name in ("summary.tsv", "status.tsv"):
            with open(out / name, newline="") as f:
                kept = [row["family"] for row in csv.DictReader(f, delimiter="\t")]
            assert kept == ["Filoviridae"], name


class TestFindFamilyDir:
    """Directory discovery via the .family ownership marker."""

    def _make_dir(self, output_dir, name, owner=None, status=None):
        d = output_dir / name
        d.mkdir(parents=True)
        if owner is not None:
            (d / ".family").write_text(owner)
        if status is not None:
            (d / ".status.json").write_text(json.dumps({"family": owner, "status": status}))
        return d

    def test_finds_taxid_suffixed_dir_via_marker(self, tmp_path):
        d = self._make_dir(tmp_path, "Filoviridae_11266", owner="Filoviridae")
        assert cli_mod._find_family_dir("Filoviridae", tmp_path) == d

    def test_finds_renamed_dir_via_marker(self, tmp_path):
        # Directory named after manual.name, not the biological family.
        d = self._make_dir(tmp_path, "Ebola_3052317", owner="Orthoebolavirus")
        assert cli_mod._find_family_dir("Orthoebolavirus", tmp_path) == d

    def test_glob_fallback_for_markerless_legacy_dir(self, tmp_path):
        d = self._make_dir(tmp_path, "Filoviridae_11266", owner=None)
        assert cli_mod._find_family_dir("Filoviridae", tmp_path) == d

    def test_plain_name_fallback(self, tmp_path):
        d = self._make_dir(tmp_path, "Filoviridae", owner=None)
        assert cli_mod._find_family_dir("Filoviridae", tmp_path) == d

    def test_returns_none_when_absent(self, tmp_path):
        assert cli_mod._find_family_dir("Nope", tmp_path) is None

    def test_marker_takes_priority_and_prefers_newest(self, tmp_path):
        import os, time
        old = self._make_dir(tmp_path, "OldName_99", owner="Hantaviridae")
        new = self._make_dir(tmp_path, "NewName_99", owner="Hantaviridae")
        # Make `new` clearly more recent.
        t = time.time()
        os.utime(old, (t - 100, t - 100))
        os.utime(new, (t, t))
        assert cli_mod._find_family_dir("Hantaviridae", tmp_path) == new

    def test_owner_family_reads_marker(self, tmp_path):
        d = self._make_dir(tmp_path, "X_1", owner="RealFamily")
        assert cli_mod._owner_family(d) == "RealFamily"
        assert cli_mod._owner_family(tmp_path / "missing") is None

    def test_find_status_file_via_marker(self, tmp_path):
        d = self._make_dir(tmp_path, "Ebola_3052317", owner="Orthoebolavirus", status="done")
        assert cli_mod._find_status_file("Orthoebolavirus", tmp_path) == d / ".status.json"

    def test_find_status_file_none_when_no_status(self, tmp_path):
        self._make_dir(tmp_path, "Ebola_3052317", owner="Orthoebolavirus")
        assert cli_mod._find_status_file("Orthoebolavirus", tmp_path) is None


# ---------------------------------------------------------------------------
# `vfam_trees test`: FastTree is optional since both trees default to IQ-TREE.
# ---------------------------------------------------------------------------

class TestToolCheckFastTreeOptional:
    def test_fasttree_not_required(self):
        assert "FastTree" not in cli_mod.REQUIRED_TOOLS
        assert "FastTree" in cli_mod.FASTTREE_TOOLS

    def test_missing_fasttree_does_not_fail_the_tool_section(self, tmp_path, monkeypatch):
        import shutil
        real_which = shutil.which
        monkeypatch.setattr(shutil, "which",
                            lambda name: None if name == "FastTree" else (real_which(name) or f"/fake/{name}"))
        # Python packages and NCBI connectivity are outside this test; only the
        # tool sections are checked for FastTree's status line.
        result = CliRunner().invoke(main, ["test", "-g", str(tmp_path / "missing.yaml")])
        out = result.output
        tools = out.split("Configuration:")[0]
        assert "[MISSING]  FastTree" not in tools
        assert "[--]  FastTree" in tools


# ---------------------------------------------------------------------------
# init --full, run records, global config checks
# ---------------------------------------------------------------------------

def _flat_keys(d: dict, prefix: str = "") -> set[str]:
    out = set()
    for k, v in d.items():
        if isinstance(v, dict) and v:
            out |= _flat_keys(v, f"{prefix}{k}.")
        else:
            out.add(prefix + k)
    return out


class TestInitTemplates:
    def test_full_template_lists_every_default(self, tmp_path):
        """`init --full` must show every setting with this version's default:
        the defaults block equals DEFAULT_FAMILY_CONFIG exactly, and the
        global-only keys are exactly GLOBAL_SETTINGS_DEFAULTS'."""
        import yaml
        from vfam_trees.config import DEFAULT_FAMILY_CONFIG, GLOBAL_SETTINGS_DEFAULTS
        out = tmp_path / "global.yaml"
        result = CliRunner().invoke(main, ["init", "--full", "-g", str(out)])
        assert result.exit_code == 0, result.output
        cfg = yaml.safe_load(out.read_text())
        assert cfg["defaults"] == DEFAULT_FAMILY_CONFIG
        top = {k: v for k, v in cfg.items() if k != "defaults"}
        assert _flat_keys(top) == _flat_keys(GLOBAL_SETTINGS_DEFAULTS)

    def test_short_template_values_match_the_defaults(self, tmp_path):
        import yaml
        from vfam_trees.config import DEFAULT_FAMILY_CONFIG
        out = tmp_path / "global.yaml"
        result = CliRunner().invoke(main, ["init", "-g", str(out)])
        assert result.exit_code == 0, result.output
        defaults = yaml.safe_load(out.read_text())["defaults"]
        for block, values in defaults.items():
            for key, value in values.items():
                assert DEFAULT_FAMILY_CONFIG[block][key] == value, f"{block}.{key}"
        assert "hmm" not in defaults  # advanced blocks only with --full

    def test_full_template_changes_no_family_config(self, tmp_path, caplog):
        """Pinning every default in global.yaml is a no-op: each preset
        family's config equals the one built from the built-in defaults, and
        the file raises no unknown-key warning."""
        import logging
        from vfam_trees import config as c
        out = tmp_path / "global.yaml"
        CliRunner().invoke(main, ["init", "--full", "-g", str(out)])
        out.write_text(out.read_text().replace("your.email@example.com", "a@b.org"))
        with caplog.at_level(logging.WARNING):
            full = c.load_global_config(out)
        assert full["_unknown_keys"] == []
        assert "Unknown config key" not in caplog.text
        families = (set(c.DNA_FAMILIES) | set(c.CONCATENATION_FAMILIES)
                    | set(c.SEGMENTED_FAMILIES) | set(c.HMM_MARKER_FAMILIES) | {"Geminiviridae"})
        for fam in sorted(families):
            a = c._generate_default_family_config(fam, full)
            b = c._generate_default_family_config(fam, c.make_minimal_global_cfg())
            assert a == b, fam


class TestRunRecord:
    def _run(self, tmp_path, extra=()):
        return CliRunner().invoke(main, [
            "run", "-f", str(_families_file(tmp_path)),
            "-o", str(tmp_path / "out"), "-j", "1", "-t", "4", *extra,
        ])

    def test_record_written_and_completed(self, tmp_path, captured_cmd):
        import yaml
        from vfam_trees import __version__
        from vfam_trees.config import DEFAULT_FAMILY_CONFIG
        result = self._run(tmp_path)
        assert result.exit_code == 0, result.output
        records = list((tmp_path / "out" / "run_records").glob("*.yaml"))
        assert len(records) == 1
        text = records[0].read_text()
        assert text.startswith("# vfam_trees run record")
        rec = yaml.safe_load(text)
        assert rec["vfam_trees_version"] == __version__
        assert rec["status"] == "completed" and rec["exit_code"] == 0
        assert rec["finished"] is not None
        assert rec["families"] == ["Poxviridae", "Filoviridae"]
        assert rec["options"] == {"cores": 1, "threads": 4, "force": False,
                                  "log_level": "INFO"}
        assert rec["tools"] == {"mafft": "7.0"}
        assert rec["global_config"]["defaults"] == DEFAULT_FAMILY_CONFIG
        assert rec["global_config"]["icon"]["size"] == 256
        # The API key never lands in the results folder.
        assert "SECRET123" not in text
        assert rec["global_config"]["ncbi"]["email"] == "test@example.org"

    def test_failed_run_is_recorded(self, tmp_path, captured_cmd, monkeypatch):
        import yaml
        monkeypatch.setattr(
            cli_mod.subprocess, "run",
            lambda cmd, *a, **k: _FakeCompleted(returncode=3 if cmd[0] == "snakemake" else 0))
        result = self._run(tmp_path)
        assert result.exit_code == 3
        rec = yaml.safe_load(
            next((tmp_path / "out" / "run_records").glob("*.yaml")).read_text())
        assert rec["status"] == "failed" and rec["exit_code"] == 3

    def test_interrupted_run_is_recorded(self, tmp_path, captured_cmd, monkeypatch):
        import yaml

        def _interrupt(cmd, *a, **k):
            if cmd[0] == "snakemake":  # Ctrl-C while Snakemake runs
                raise KeyboardInterrupt
            return _FakeCompleted()
        monkeypatch.setattr(cli_mod.subprocess, "run", _interrupt)
        self._run(tmp_path)
        rec = yaml.safe_load(
            next((tmp_path / "out" / "run_records").glob("*.yaml")).read_text())
        assert rec["status"] == "interrupted" and rec["exit_code"] is None

    def test_two_runs_in_the_same_second_keep_both_records(self, tmp_path, captured_cmd):
        assert self._run(tmp_path).exit_code == 0
        assert self._run(tmp_path).exit_code == 0
        assert len(list((tmp_path / "out" / "run_records").glob("*.yaml"))) == 2

    def test_unknown_global_keys_recorded(self, tmp_path, captured_cmd):
        import yaml
        (tmp_path / "config" / "global.yaml").write_text(
            'ncbi:\n  email: "test@example.org"\noutput_dir: results\n'
            'defaults:\n  length_outlier:\n    hi_mult: 3.0\n'
        )
        result = self._run(tmp_path)
        assert result.exit_code == 0, result.output
        rec = yaml.safe_load(
            next((tmp_path / "out" / "run_records").glob("*.yaml")).read_text())
        assert rec["global_config_unknown_keys"] == [
            "output_dir", "defaults.length_outlier.hi_mult"]

    def test_bad_global_config_stops_before_snakemake(self, tmp_path, captured_cmd):
        (tmp_path / "config" / "global.yaml").write_text("ncbi:\n  email: \"\"\n")
        result = self._run(tmp_path)
        assert result.exit_code == 1
        assert "ncbi.email must be set" in result.output
        assert "cmd" not in captured_cmd
        assert not (tmp_path / "out" / "run_records").exists()

    def test_dry_run_writes_no_record(self, tmp_path, captured_cmd):
        self._run(tmp_path, extra=("--dry-run",))
        assert not (tmp_path / "out" / "run_records").exists()


class TestToolPreflight:
    def _which(self, absent):
        return lambda tool: None if tool in absent else f"/usr/bin/{tool}"

    def test_hmmer_needed_only_by_hmm_families(self, tmp_path, monkeypatch):
        import shutil
        from vfam_trees.config import make_minimal_global_cfg
        monkeypatch.setattr(shutil, "which", self._which({"hmmscan", "hmmpress"}))
        g = make_minimal_global_cfg()
        # Togaviridae: RdRp preset (HMM); Kolmioviridae: whole-genome nucleotide.
        assert cli_mod._missing_tools(["Togaviridae"], tmp_path, g) == ["hmmscan", "hmmpress"]
        assert cli_mod._missing_tools(["Kolmioviridae"], tmp_path, g) == []
        assert not list(tmp_path.iterdir())  # nothing written

    def test_family_file_opting_out_needs_no_hmmer(self, tmp_path, monkeypatch):
        import shutil
        from vfam_trees.config import make_minimal_global_cfg
        monkeypatch.setattr(shutil, "which", self._which({"hmmscan", "hmmpress"}))
        (tmp_path / "Togaviridae.yaml").write_text(
            "sequence:\n  type: nucleotide\n  region: whole_genome\nhmm:\n  enabled: false\n")
        assert cli_mod._missing_tools(["Togaviridae"], tmp_path, make_minimal_global_cfg()) == []

    def test_required_and_fasttree(self, tmp_path, monkeypatch):
        import shutil
        from vfam_trees.config import make_minimal_global_cfg
        monkeypatch.setattr(shutil, "which", self._which({"mafft", "FastTree"}))
        (tmp_path / "Geminiviridae.yaml").write_text(
            'tree_500:\n  tool: fasttree\n  options_nuc: ""\n  options_aa: ""\n')
        assert cli_mod._missing_tools(["Geminiviridae"], tmp_path, make_minimal_global_cfg()) == [
            "mafft", "FastTree"]

    def test_run_stops_before_snakemake(self, tmp_path, captured_cmd, monkeypatch):
        monkeypatch.setattr(cli_mod, "_missing_tools", lambda *a: ["hmmscan", "hmmpress"])
        result = CliRunner().invoke(main, [
            "run", "-f", str(_families_file(tmp_path)), "-o", str(tmp_path / "out")])
        assert result.exit_code == 1
        assert "needs tools that are not on PATH: hmmscan, hmmpress" in result.output
        assert "cmd" not in captured_cmd
