"""Tests for vfam_trees.pipeline_concat — concat dispatch + skip path.

The full concat run path requires MMseqs2 / MAFFT / IQ-TREE / Entrez and is
exercised in integration when the pipeline runs.  These tests verify the
public interface and the skip-path with stubbed fetch — no external tools
required.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from vfam_trees import pipeline_concat


# ---------------------------------------------------------------------------
# Module surface
# ---------------------------------------------------------------------------

def test_run_family_concat_is_exported():
    assert hasattr(pipeline_concat, "run_family_concat")
    assert callable(pipeline_concat.run_family_concat)


# ---------------------------------------------------------------------------
# Skip path: no genomes pass min_fraction → cleanly emits skip row.
# ---------------------------------------------------------------------------

def _minimal_family_cfg() -> dict:
    return {
        "download": {"max_per_species": 100},
        "sequence": {"type": "protein", "region": "concatenated", "segment": None},
        "quality": {
            "min_length": None,
            "max_ambiguous": 0.01,
            "exclude_organisms": [],
        },
        "clustering": {
            "tool": "mmseqs2",
            "threshold_min": 0.7,
            "threshold_max": 0.99,
            "max_reps_500": 20,
            "max_reps_100": 5,
        },
        "targets":   {"max_500": 500, "max_100": 100},
        "msa_500":   {"tool": "mafft", "options_aa": "--6merpair --retree 1"},
        "msa_100":   {"tool": "mafft", "options_aa": "--auto"},
        "msa_trim":  {"enabled": True, "tool": "trimal", "options": "-automated1"},
        "tree_500":  {"tool": "fasttree", "options": "", "model_aa": "LG+G"},
        "tree_100":  {"tool": "iqtree", "options_aa": "-B 1000", "model_aa": "TEST"},
        "refseq_absorption": {"enabled": True, "threshold": 0.99},
        "length_outlier":    {"enabled": True, "hi_mult": 3.0, "lo_mult": 1.0 / 3.0},
        "outlier_removal":   {"enabled": True, "factor": 20.0,
                              "max_iterations": 3, "min_seqs": 40},
        "taxonomy":          {"lca_min_rank": "none"},
        "concatenation": {
            "proteins": [
                {"name": "DNA polymerase", "aliases": []},
                {"name": "MCP",            "aliases": []},
            ],
            "min_fraction": 0.7,
            "partition_tree_100": True,
            "partition_tree_500": False,
        },
    }


def test_skip_path_when_fetch_returns_no_genomes(tmp_path, monkeypatch):
    """When fetch_species_genomes returns no genomes for any species, the
    runner emits the 'too few genomes' skip and writes summary/status rows.
    """
    # Stub the heavy I/O — Entrez + per-marker fetch
    def fake_fetch_species_genomes(*args, **kwargs):
        return {}, {"n_proteins_fetched": 0, "n_genomes_found": 0,
                    "n_genomes_kept": 0, "n_dropped_min_fraction": 0,
                    "n_dropped_split_submission": 0, "n_orphaned_no_source": 0}
    def fake_fetch_taxonomy_lineages(taxids):
        return {}
    monkeypatch.setattr(pipeline_concat, "fetch_species_genomes", fake_fetch_species_genomes)
    monkeypatch.setattr(pipeline_concat, "fetch_taxonomy_lineages", fake_fetch_taxonomy_lineages)

    # Capture mark_skipped / mark_done invocations
    mark_skipped_calls: list[tuple] = []
    mark_done_calls: list[tuple] = []
    def fake_mark_skipped(family_dir, family, output_dir, reason):
        mark_skipped_calls.append((family, reason))
    def fake_mark_done(family_dir, family, output_dir):
        mark_done_calls.append((family,))

    family_dir = tmp_path / "Poxviridae_10240"
    family_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    summary_path = output_dir / "summary.tsv"
    status_path  = output_dir / "status.tsv"

    species_list = [
        {"name": "Test virus 1", "taxid": 100001},
        {"name": "Test virus 2", "taxid": 100002},
    ]

    pipeline_concat.run_family_concat(
        family="Poxviridae",
        family_cfg=_minimal_family_cfg(),
        family_taxid=10240,
        family_lineage=[{"rank": "family", "name": "Poxviridae"}],
        family_annotation={"baltimore_class": "I"},
        family_dir=family_dir,
        work_dir=work_dir,
        output_dir=output_dir,
        species_list=species_list,
        threads=1,
        summary_path=summary_path,
        status_path=status_path,
        mark_done=fake_mark_done,
        mark_skipped=fake_mark_skipped,
    )

    # Skip path must have run
    assert len(mark_skipped_calls) == 1
    fam, reason = mark_skipped_calls[0]
    assert fam == "Poxviridae"
    assert "too few genomes" in reason
    assert mark_done_calls == []

    # Summary + status rows written
    assert summary_path.exists()
    assert status_path.exists()
    summary_text = summary_path.read_text()
    status_text = status_path.read_text()
    assert "Poxviridae" in summary_text
    assert "Poxviridae" in status_text
    assert "too few genomes" in status_text
