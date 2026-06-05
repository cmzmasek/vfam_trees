"""Tests for the external data-source fetchers (Pathoplexus / LAPIS).

All HTTP is mocked by patching ``datasources._http_get`` — no network access.
"""
import json

import pytest

from vfam_trees import datasources as ds


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

class TestBuildFilters:
    def test_maps_all_limiters_to_lapis_params(self):
        params = ds._build_filters({
            "country": "Uganda",
            "host": "Homo sapiens",
            "date_from": "2024-01-01",
            "date_to": "2024-12-31",
        })
        assert params == {
            "geoLocCountry": "Uganda",
            "hostNameScientific": "Homo sapiens",
            "sampleCollectionDateRangeLowerFrom": "2024-01-01",
            "sampleCollectionDateRangeUpperTo": "2024-12-31",
        }

    def test_omits_absent_and_falsy_limiters(self):
        assert ds._build_filters({"country": None, "host": ""}) == {}
        assert ds._build_filters({}) == {}

    def test_date_uses_range_bound_fields_not_raw_string(self):
        # the raw sampleCollectionDate field is not range-queryable
        params = ds._build_filters({"date_from": "2024-01-01"})
        assert "sampleCollectionDate" not in params
        assert params == {"sampleCollectionDateRangeLowerFrom": "2024-01-01"}


class TestParseFasta:
    def test_parses_multi_record_and_joins_wrapped_lines(self):
        out = ds._parse_fasta(">PP_1 some desc\nACGT\nTTGG\n>PP_2\nAAAA\n")
        assert out == {"PP_1": "ACGTTTGG", "PP_2": "AAAA"}

    def test_header_id_is_first_token(self):
        out = ds._parse_fasta(">PP_000LAAX.1 Zaire ebolavirus|x\nACGT")
        assert list(out) == ["PP_000LAAX.1"]


# ---------------------------------------------------------------------------
# fetch_pathoplexus (mocked transport)
# ---------------------------------------------------------------------------

def _wire_responses(monkeypatch, *, count, details, fasta):
    """Patch _http_get to answer aggregated / details / sequences by URL."""
    calls = []

    def fake_get(url):
        calls.append(url)
        if "/aggregated" in url:
            return json.dumps({"data": [{"count": count}], "info": {}}).encode()
        if "/details" in url:
            return json.dumps({"data": details, "info": {}}).encode()
        if "/unalignedNucleotideSequences" in url:
            return fasta.encode()
        raise AssertionError(f"unexpected URL: {url}")

    monkeypatch.setattr(ds, "_http_get", fake_get)
    return calls


def _rec(acc, *, insdc="", taxid="186538", name="Orthoebolavirus zairense",
         host="Homo sapiens", country="Uganda", date="2024-03-15"):
    return {
        "accessionVersion": acc,
        "insdcAccessionBase": insdc,
        "ncbiVirusTaxId": taxid,
        "ncbiVirusName": name,
        "hostNameScientific": host,
        "geoLocCountry": country,
        "sampleCollectionDate": date,
    }


class TestFetchPathoplexus:
    def test_returns_joined_entries_with_real_metadata(self, monkeypatch):
        _wire_responses(
            monkeypatch,
            count=2,
            details=[_rec("PP_1"), _rec("PP_2", country="Kenya")],
            fasta=">PP_1\nACGT\n>PP_2\nTTGG\n",
        )
        out = ds.fetch_pathoplexus(
            {"organism": "ebola-zaire"}, ncbi_accessions=set()
        )
        assert {e["id"] for e in out} == {"PP_1", "PP_2"}
        e1 = next(e for e in out if e["id"] == "PP_1")
        assert e1["sequence"] == "ACGT"
        assert e1["species"] == "Orthoebolavirus zairense"
        assert e1["host"] == "Homo sapiens"
        assert e1["location"] == "Uganda"
        assert e1["collection_date"] == "2024-03-15"
        assert e1["taxon_id"] == "186538"
        assert e1["source"] == "pathoplexus"

    def test_zero_count_short_circuits(self, monkeypatch):
        calls = _wire_responses(monkeypatch, count=0, details=[], fasta="")
        out = ds.fetch_pathoplexus({"organism": "cchf"}, ncbi_accessions=set())
        assert out == []
        # only the aggregated preflight was issued
        assert len(calls) == 1 and "/aggregated" in calls[0]

    def test_dedup_vs_ncbi_drops_insdc_mirrored_records(self, monkeypatch):
        _wire_responses(
            monkeypatch,
            count=2,
            details=[_rec("PP_1", insdc="MN908947"), _rec("PP_2", insdc="OQ123456")],
            fasta=">PP_1\nACGT\n>PP_2\nTTGG\n",
        )
        # NCBI set carries a versioned accession; comparison strips the version
        out = ds.fetch_pathoplexus(
            {"organism": "ebola-zaire", "dedup_vs_ncbi": True},
            ncbi_accessions={"MN908947.1"},
        )
        assert {e["id"] for e in out} == {"PP_2"}

    def test_dedup_disabled_keeps_mirrored_records(self, monkeypatch):
        _wire_responses(
            monkeypatch,
            count=1,
            details=[_rec("PP_1", insdc="MN908947")],
            fasta=">PP_1\nACGT\n",
        )
        out = ds.fetch_pathoplexus(
            {"organism": "ebola-zaire", "dedup_vs_ncbi": False},
            ncbi_accessions={"MN908947.1"},
        )
        assert {e["id"] for e in out} == {"PP_1"}

    def test_max_seqs_passed_as_limit(self, monkeypatch):
        calls = _wire_responses(
            monkeypatch, count=999, details=[_rec("PP_1")], fasta=">PP_1\nACGT\n",
        )
        ds.fetch_pathoplexus(
            {"organism": "ebola-zaire", "max_seqs": 25}, ncbi_accessions=set()
        )
        detail_url = next(u for u in calls if "/details" in u)
        assert "limit=25" in detail_url
        assert "orderBy=accessionVersion" in detail_url

    def test_sequence_without_metadata_is_skipped(self, monkeypatch):
        _wire_responses(
            monkeypatch,
            count=2,
            details=[_rec("PP_1")],  # PP_2 absent from metadata
            fasta=">PP_1\nACGT\n>PP_2\nTTGG\n",
        )
        out = ds.fetch_pathoplexus({"organism": "ebola-zaire"}, ncbi_accessions=set())
        assert {e["id"] for e in out} == {"PP_1"}

    def test_collection_date_falls_back_to_range_upper(self, monkeypatch):
        rec = _rec("PP_1", date="")
        rec["sampleCollectionDateRangeUpper"] = "2006-02-18"
        _wire_responses(monkeypatch, count=1, details=[rec], fasta=">PP_1\nACGT\n")
        out = ds.fetch_pathoplexus({"organism": "ebola-zaire"}, ncbi_accessions=set())
        assert out[0]["collection_date"] == "2006-02-18"

    def test_filters_are_sent_to_all_three_endpoints(self, monkeypatch):
        calls = _wire_responses(
            monkeypatch, count=1, details=[_rec("PP_1")], fasta=">PP_1\nACGT\n",
        )
        ds.fetch_pathoplexus(
            {"organism": "ebola-zaire", "country": "Uganda"}, ncbi_accessions=set()
        )
        assert all("geoLocCountry=Uganda" in u for u in calls)

    def test_restricts_to_latest_version(self, monkeypatch):
        # Pathoplexus keeps every revision; we must request only the current one.
        calls = _wire_responses(
            monkeypatch, count=1, details=[_rec("PP_1")], fasta=">PP_1\nACGT\n",
        )
        ds.fetch_pathoplexus({"organism": "ebola-zaire"}, ncbi_accessions=set())
        assert all("versionStatus=LATEST_VERSION" in u for u in calls)

    def test_missing_taxid_uses_organism_fallback(self, monkeypatch):
        # Fresh outbreak records lack ncbiVirusName/ncbiVirusTaxId.
        _wire_responses(
            monkeypatch, count=1,
            details=[_rec("PP_1", taxid="", name="")], fasta=">PP_1\nACGT\n",
        )
        out = ds.fetch_pathoplexus({"organism": "ebola-bdbv"}, ncbi_accessions=set())
        assert out[0]["taxon_id"] == "565995"              # Bundibugyo fallback
        assert out[0]["species"] == "Bundibugyo ebolavirus"

    def test_present_taxid_not_overridden_by_fallback(self, monkeypatch):
        _wire_responses(
            monkeypatch, count=1,
            details=[_rec("PP_1", taxid="186538", name="Zaire ebolavirus")],
            fasta=">PP_1\nACGT\n",
        )
        out = ds.fetch_pathoplexus({"organism": "ebola-bdbv"}, ncbi_accessions=set())
        # the record's own NCBI taxid/name win over the organism fallback
        assert out[0]["taxon_id"] == "186538"
        assert out[0]["species"] == "Zaire ebolavirus"


class TestRegistry:
    def test_pathoplexus_registered(self):
        assert ds.SOURCE_FETCHERS["pathoplexus"] is ds.fetch_pathoplexus

    def test_known_organisms_count(self):
        assert len(ds.KNOWN_PATHOPLEXUS_ORGANISMS) == 14
        assert "ebola-zaire" in ds.KNOWN_PATHOPLEXUS_ORGANISMS

    def test_taxon_map_covers_all_organisms(self):
        # KNOWN_PATHOPLEXUS_ORGANISMS is derived from the taxon map's keys
        assert set(ds.PATHOPLEXUS_ORGANISM_TAXON) == set(ds.KNOWN_PATHOPLEXUS_ORGANISMS)
        for name, taxid in ds.PATHOPLEXUS_ORGANISM_TAXON.values():
            assert name and taxid.isdigit()


class TestHttpGet:
    def test_4xx_fails_fast_without_retry(self, monkeypatch):
        import io
        import urllib.error

        attempts = {"n": 0}

        def boom(req, timeout):
            attempts["n"] += 1
            raise urllib.error.HTTPError(
                req.full_url, 400, "Bad Request", {}, io.BytesIO(b"bad param")
            )

        monkeypatch.setattr(ds.urllib.request, "urlopen", boom)
        with pytest.raises(RuntimeError, match="HTTP 400"):
            ds._http_get("https://lapis.pathoplexus.org/ebola-zaire/sample/aggregated")
        assert attempts["n"] == 1  # no retry on client error
