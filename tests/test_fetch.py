"""Tests for vfam_trees.fetch — query building (no network calls)."""
import pytest

from vfam_trees.fetch import _build_species_query


def test_whole_genome_nucleotide():
    q = _build_species_query(12345, "nucleotide", "whole_genome")
    assert "complete genome" in q or "complete sequence" in q
    assert "[Gene]" not in q


def test_whole_genome_protein():
    q = _build_species_query(12345, "protein", "whole_genome")
    assert "complete genome" not in q
    assert "[Gene]" not in q


def test_segment_query():
    q = _build_species_query(12345, "nucleotide", "whole_genome", segment="RNA1")
    assert '"RNA1"[Title]' in q
    assert "complete sequence" in q


def test_marker_gene_query():
    q = _build_species_query(12345, "nucleotide", "B646L")
    assert '"B646L"[Gene]' in q


def test_marker_gene_excludes_complete_genome():
    q = _build_species_query(12345, "nucleotide", "B646L")
    assert 'NOT "complete genome"[Title]' in q
    assert 'NOT "complete sequence"[Title]' in q


def test_marker_gene_hexon_excludes_complete():
    q = _build_species_query(12345, "nucleotide", "hexon")
    assert '"hexon"[Gene]' in q
    assert 'NOT "complete genome"[Title]' in q


def test_refseq_only_adds_filter():
    q = _build_species_query(12345, "nucleotide", "whole_genome", refseq_only=True)
    assert "refseq[filter]" in q


def test_patent_always_excluded():
    for region in ("whole_genome", "B646L"):
        q = _build_species_query(12345, "nucleotide", region)
        assert "patent[filter]" in q


def test_exclude_organisms():
    q = _build_species_query(12345, "nucleotide", "whole_genome",
                             exclude_organisms=["metagenome", "synthetic construct"])
    assert '"metagenome"[Organism]' in q
    assert '"synthetic construct"[Organism]' in q


def test_taxid_in_query():
    q = _build_species_query(99999, "nucleotide", "whole_genome")
    assert "txid99999" in q


def test_protein_marker_uses_protein_name_field():
    q = _build_species_query(12345, "protein", "DNA polymerase")
    assert '"DNA polymerase"[Protein Name]' in q


def test_protein_marker_also_includes_gene_fallback():
    q = _build_species_query(12345, "protein", "B646L")
    assert '"B646L"[Protein Name]' in q
    assert '"B646L"[Gene]' in q


def test_protein_marker_does_not_use_nucleotide_gene_only():
    q = _build_species_query(12345, "protein", "DNA polymerase")
    # [Gene] alone (without [Protein Name]) would miss most protein records
    assert '"DNA polymerase"[Protein Name]' in q
