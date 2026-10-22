"""Shared test setup.

The translation fallback and the species record count ask NCBI when a
species has no marker, which many pipeline and fetch tests produce on
purpose.  They are stubbed here (no records), so no test reaches the
network; a test of the helpers themselves opts out with
``@pytest.mark.real_ncbi_helpers`` and stubs the Entrez layer below them.
"""
import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "real_ncbi_helpers: keep fetch._fetch_species_nucleotide_records and "
        "fetch.count_species_records unstubbed",
    )


@pytest.fixture(autouse=True)
def _no_ncbi_fallback(request, monkeypatch):
    if request.node.get_closest_marker("real_ncbi_helpers"):
        return
    import vfam_trees.fetch as fetch_mod
    monkeypatch.setattr(fetch_mod, "_fetch_species_nucleotide_records", lambda *a, **k: [])
    monkeypatch.setattr(fetch_mod, "count_species_records", lambda *a, **k: 0)
