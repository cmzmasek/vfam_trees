"""Tests for vfam_trees.provenance — run records and family settings headers."""
import logging
from datetime import datetime
from pathlib import Path

import yaml

from vfam_trees import __version__
from vfam_trees import provenance as prov


def _start(tmp_path, monkeypatch, global_cfg=None):
    monkeypatch.setattr(prov, "tool_versions", lambda: {"mafft": "7.520"})
    return prov.start_run_record(
        tmp_path / "out",
        argv=["vfam_trees", "run", "-f", "fams.txt"],
        families_file=tmp_path / "fams.txt",
        families=["Togaviridae"],
        options={"cores": 1, "threads": 8, "force": False, "log_level": "INFO"},
        paths={"output_dir": str(tmp_path / "out")},
        global_cfg=global_cfg or {"ncbi": {"email": "a@b.org", "api_key": "KEY"}},
    )


class TestRunRecord:
    def test_contents(self, tmp_path, monkeypatch):
        path = _start(tmp_path, monkeypatch)
        assert path.parent == tmp_path / "out" / "run_records"
        rec = yaml.safe_load(path.read_text())
        assert rec["vfam_trees_version"] == __version__
        assert rec["status"] == "running" and rec["finished"] is None
        assert rec["command"] == "vfam_trees run -f fams.txt"
        assert rec["tools"] == {"mafft": "7.520"}
        assert rec["python_packages"]["python"]
        assert rec["python_packages"]["biopython"] != "not installed"
        assert rec["global_config"]["ncbi"]["api_key"] == "<set, not recorded>"
        assert "KEY" not in path.read_text()
        assert not list(path.parent.glob("*.part"))

    def test_finish(self, tmp_path, monkeypatch):
        path = _start(tmp_path, monkeypatch)
        prov.finish_run_record(path, 0)
        text = path.read_text()
        assert text.startswith("# vfam_trees run record")  # header kept
        rec = yaml.safe_load(text)
        assert (rec["status"], rec["exit_code"]) == ("completed", 0)
        assert datetime.fromisoformat(rec["finished"]) >= datetime.fromisoformat(rec["started"])

    def test_finish_missing_file_only_warns(self, tmp_path, caplog):
        with caplog.at_level(logging.WARNING):
            prov.finish_run_record(tmp_path / "gone.yaml", 0)
        assert "Could not complete the run record" in caplog.text

    def test_missing_tool(self, monkeypatch):
        monkeypatch.setattr(prov.shutil, "which", lambda name: None)
        assert set(prov.tool_versions().values()) == {"not found"}


class TestFamilyHeader:
    def test_header(self):
        h = prov.family_config_header(
            family="Togaviridae", started=datetime(2026, 10, 21, 10, 9, 19), threads=8,
            global_config_path=Path("/run/config/global.yaml"),
            family_config_path=Path("/run/configs/Togaviridae.yaml"),
            auto_generated=False, unknown_keys=[],
        )
        assert h.startswith(f"# Settings vfam_trees {__version__} used for Togaviridae "
                            "(started 2026-10-21 10:09:19).")
        assert "/run/configs/Togaviridae.yaml (existing file)" in h
        assert "Not read by this version" not in h
        assert h.endswith("\n\n")
        assert all(line.startswith("#") for line in h.strip().splitlines())
