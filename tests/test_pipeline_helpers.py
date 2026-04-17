"""Tests for pipeline helper functions."""
import json
from pathlib import Path

import pytest

from vfam_trees.pipeline import _mark_done, _mark_skipped


@pytest.fixture
def dirs(tmp_path):
    family_dir = tmp_path / "Flaviviridae_11050"
    family_dir.mkdir()
    output_dir = tmp_path
    return family_dir, output_dir


def test_mark_done_writes_status(dirs):
    family_dir, output_dir = dirs
    _mark_done(family_dir, "Flaviviridae", output_dir)
    status = json.loads((family_dir / ".status.json").read_text())
    assert status["status"] == "done"
    assert status["family"] == "Flaviviridae"


def test_mark_done_writes_sentinel(dirs):
    family_dir, output_dir = dirs
    _mark_done(family_dir, "Flaviviridae", output_dir)
    assert (output_dir / ".done_Flaviviridae").exists()


def test_mark_skipped_writes_status(dirs):
    family_dir, output_dir = dirs
    _mark_skipped(family_dir, "Flaviviridae", output_dir, "no species found")
    status = json.loads((family_dir / ".status.json").read_text())
    assert status["status"] == "skipped"
    assert "no species" in status["reason"]


def test_mark_skipped_writes_sentinel(dirs):
    family_dir, output_dir = dirs
    _mark_skipped(family_dir, "Flaviviridae", output_dir, "too few sequences")
    assert (output_dir / ".done_Flaviviridae").exists()
