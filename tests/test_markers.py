"""Tests for vfam_trees.markers — NameMatchIdentifier (v1)."""
from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from vfam_trees.markers import NameMatchIdentifier


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rec(acc: str, description: str, length: int = 500,
         locus_tag: str | None = None, product: str | None = None) -> SeqRecord:
    """Build a minimal protein-like SeqRecord for matcher input."""
    rec = SeqRecord(Seq("M" * length), id=acc, description=description)
    rec.annotations = {"organism": "Test virus"}
    if locus_tag:
        rec.annotations["locus_tag"] = locus_tag
    if product:
        feat = SeqFeature(FeatureLocation(0, length), type="Protein")
        feat.qualifiers = {"product": [product]}
        rec.features = [feat]
    else:
        rec.features = []
    return rec


# Common fixture: a tiny marker set
_POLB = {"name": "DNA polymerase", "aliases": ["DNA-directed DNA polymerase", "DNA pol"]}
_MCP  = {"name": "major capsid protein", "aliases": ["MCP"]}


# ---------------------------------------------------------------------------
# Basic name + alias matching
# ---------------------------------------------------------------------------

class TestBasicMatching:
    def setup_method(self):
        self.id = NameMatchIdentifier()

    def test_name_match(self):
        proteins = [_rec("YP_1", "DNA polymerase [Vaccinia virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert "DNA polymerase" in out
        assert out["DNA polymerase"].id == "YP_1"

    def test_alias_match(self):
        proteins = [_rec("YP_1", "DNA pol [Vaccinia virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_case_insensitive(self):
        proteins = [_rec("YP_1", "DNA POLYMERASE [some virus]")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_missing_marker_omitted_from_result(self):
        proteins = [_rec("YP_1", "DNA polymerase")]
        out = self.id.identify(proteins, [_POLB, _MCP])
        assert "DNA polymerase" in out
        assert "major capsid protein" not in out

    def test_empty_proteins_returns_empty_dict(self):
        out = self.id.identify([], [_POLB, _MCP])
        assert out == {}

    def test_empty_marker_set_returns_empty_dict(self):
        out = self.id.identify([_rec("YP_1", "DNA polymerase")], [])
        assert out == {}

    def test_multiple_markers_all_found(self):
        proteins = [
            _rec("YP_1", "DNA polymerase"),
            _rec("YP_2", "major capsid protein"),
        ]
        out = self.id.identify(proteins, [_POLB, _MCP])
        assert out["DNA polymerase"].id == "YP_1"
        assert out["major capsid protein"].id == "YP_2"

    def test_match_against_protein_feature_product(self):
        # Description is generic; product qualifier carries the marker name.
        proteins = [_rec("YP_1", "hypothetical protein", product="DNA polymerase")]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"


# ---------------------------------------------------------------------------
# Subfamily-aware aliases
# ---------------------------------------------------------------------------

class TestSubfamilyAwareAliases:
    def setup_method(self):
        self.id = NameMatchIdentifier()
        self.marker = {
            "name": "DNA polymerase",
            "aliases": ["DNA pol"],
            "aliases_Entomopoxvirinae": ["DNA polymerase B"],
        }

    def test_subfamily_alias_used_when_lineage_matches(self):
        proteins = [_rec("YP_1", "DNA polymerase B [Amsacta entomopoxvirus]")]
        lineage = [
            {"rank": "family", "name": "Poxviridae"},
            {"rank": "subfamily", "name": "Entomopoxvirinae"},
            {"rank": "genus", "name": "Alphaentomopoxvirus"},
        ]
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert out["DNA polymerase"].id == "YP_1"

    def test_subfamily_alias_ignored_when_lineage_is_other_subfamily(self):
        # Chordopoxvirinae lineage → entomopox-specific alias not added.
        # The protein only matches "DNA polymerase B", which is entomopox-specific.
        # Note: "DNA polymerase B" still contains "DNA pol" substring — so it
        # would still match via the base alias.  Use a description that does
        # NOT contain any base alias to test true exclusion.
        proteins = [_rec("YP_1", "polymerase variant Z")]  # no base alias hit
        lineage = [
            {"rank": "subfamily", "name": "Chordopoxvirinae"},
        ]
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert "DNA polymerase" not in out

    def test_subfamily_alias_ignored_when_no_lineage_provided(self):
        # No lineage → only base aliases used.
        proteins = [_rec("YP_1", "polymerase variant Z")]
        out = self.id.identify(proteins, [self.marker], species_lineage=None)
        assert "DNA polymerase" not in out

    def test_subfamily_alias_ignored_when_lineage_has_no_subfamily_rank(self):
        proteins = [_rec("YP_1", "polymerase variant Z")]
        lineage = [{"rank": "family", "name": "Poxviridae"}]  # no subfamily entry
        out = self.id.identify(proteins, [self.marker], species_lineage=lineage)
        assert "DNA polymerase" not in out


# ---------------------------------------------------------------------------
# Tiebreakers
# ---------------------------------------------------------------------------

class TestTiebreakers:
    def setup_method(self):
        self.id = NameMatchIdentifier()

    def test_locus_tag_hint_wins(self):
        # Two candidates both match the name; only one has the matching locus_tag.
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        proteins = [
            _rec("YP_1", "DNA polymerase", locus_tag="A24R"),
            _rec("YP_2", "DNA polymerase", locus_tag="E9L"),
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_locus_tag_hint_regex_alternation(self):
        marker = {**_POLB, "locus_tag_hint": r"E9L|polB"}
        proteins = [
            _rec("YP_1", "DNA polymerase", locus_tag="random"),
            _rec("YP_2", "DNA polymerase", locus_tag="some-polB-thing"),
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_locus_tag_hint_no_match_falls_through(self):
        # No candidate has the hint locus tag; fall through to next tiebreaker.
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        proteins = [
            _rec("YP_1", "DNA polymerase", length=900),
            _rec("YP_2", "DNA polymerase", length=1100),  # longest wins
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_2"

    def test_length_range_midpoint(self):
        marker = {**_POLB, "length_range": [800, 1200]}  # midpoint = 1000
        proteins = [
            _rec("YP_1", "DNA polymerase", length=950),   # |950-1000| = 50
            _rec("YP_2", "DNA polymerase", length=1100),  # |1100-1000| = 100
            _rec("YP_3", "DNA polymerase", length=995),   # |995-1000| = 5  ← winner
        ]
        out = self.id.identify(proteins, [marker])
        assert out["DNA polymerase"].id == "YP_3"

    def test_longest_wins_when_no_hints(self):
        proteins = [
            _rec("YP_1", "DNA polymerase", length=500),
            _rec("YP_2", "DNA polymerase", length=1000),
            _rec("YP_3", "DNA polymerase", length=800),
        ]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_2"

    def test_lowest_accession_breaks_length_ties(self):
        proteins = [
            _rec("YP_3", "DNA polymerase", length=900),
            _rec("YP_1", "DNA polymerase", length=900),
            _rec("YP_2", "DNA polymerase", length=900),
        ]
        out = self.id.identify(proteins, [_POLB])
        assert out["DNA polymerase"].id == "YP_1"

    def test_tiebreaker_priority_locus_tag_then_length(self):
        # Locus-tag hint defined AND length range defined.  Locus-tag must win.
        marker = {
            **_POLB,
            "locus_tag_hint": r"E9L",
            "length_range": [900, 1100],   # midpoint = 1000
        }
        proteins = [
            # Hint match but length far from midpoint (700)
            _rec("YP_1", "DNA polymerase", length=700, locus_tag="E9L"),
            # No hint match but length right at midpoint (1000)
            _rec("YP_2", "DNA polymerase", length=1000, locus_tag="random"),
        ]
        out = self.id.identify(proteins, [marker])
        # Locus-tag hint is priority 1; YP_1 wins despite being further from midpoint.
        assert out["DNA polymerase"].id == "YP_1"

    def test_locus_tag_via_feature_qualifier(self):
        # When annotations don't carry locus_tag, fall back to feature qualifier.
        rec = _rec("YP_1", "DNA polymerase")
        feat = SeqFeature(FeatureLocation(0, 500), type="CDS")
        feat.qualifiers = {"locus_tag": ["E9L"]}
        rec.features = [feat]
        marker = {**_POLB, "locus_tag_hint": r"E9L"}
        out = self.id.identify([rec], [marker])
        assert out["DNA polymerase"].id == "YP_1"


# ---------------------------------------------------------------------------
# HMMIdentifier — driven by PRE-ANNOTATED hits (no live hmmscan needed)
# ---------------------------------------------------------------------------

from vfam_trees.markers import HMMIdentifier, make_identifier


def _hit(target: str, *, passing: bool = True, dom_evalue: float = 1e-30,
         ali_from: int = 1, ali_to: int = 100) -> dict:
    return {
        "target": target, "passing": passing, "dom_evalue": dom_evalue,
        "ali_from": ali_from, "ali_to": ali_to,
    }


def _hrec(acc: str, length: int, hits: list[dict],
          description: str = "") -> SeqRecord:
    """Protein-like record carrying pre-computed HMM hits (skips hmmscan)."""
    rec = SeqRecord(Seq("M" * length), id=acc, description=description or acc)
    rec.annotations = {"hmm_hits": hits}
    rec.features = []
    return rec


class TestMakeIdentifier:
    def test_disabled_returns_namematch(self):
        assert type(make_identifier({"hmm": {"enabled": False}})).__name__ == \
            "NameMatchIdentifier"

    def test_missing_block_returns_namematch(self):
        assert type(make_identifier({})).__name__ == "NameMatchIdentifier"

    def test_enabled_returns_hmm(self):
        assert type(make_identifier({"hmm": {"enabled": True}})).__name__ == \
            "HMMIdentifier"


class TestHMMIdentifier:
    ident = HMMIdentifier({"hmm": {"enabled": True}})

    def test_single_token_picks_satisfying(self):
        hexon = _hrec("A1", 900, [_hit("Adeno_hexon")])
        penton = _hrec("A2", 500, [_hit("Adeno_penton")])
        out = self.ident.identify(
            [hexon, penton],
            [{"name": "hexon", "hmms": ["Adeno_hexon"]},
             {"name": "penton", "hmms": ["Adeno_penton"]}],
        )
        assert out["hexon"].id == "A1"
        assert out["penton"].id == "A2"

    def test_longest_satisfying_wins(self):
        short = _hrec("S", 300, [_hit("Adeno_hexon")])
        long_ = _hrec("L", 950, [_hit("Adeno_hexon")])
        out = self.ident.identify(
            [short, long_], [{"name": "hexon", "hmms": ["Adeno_hexon"]}],
        )
        assert out["hexon"].id == "L"

    def test_non_passing_hit_ignored(self):
        rec = _hrec("A1", 900, [_hit("Adeno_hexon", passing=False)])
        out = self.ident.identify(
            [rec], [{"name": "hexon", "hmms": ["Adeno_hexon"]}],
        )
        assert "hexon" not in out  # HMM-authoritative: no passing hit → dropped

    def test_hmm_authoritative_no_name_fallback(self):
        # Record name says "hexon" but carries no passing HMM hit → still dropped
        # (a spec with hmms is HMM-authoritative; aliases are NOT consulted).
        rec = _hrec("A1", 900, [], description="adenovirus hexon protein")
        out = self.ident.identify(
            [rec],
            [{"name": "hexon", "aliases": ["hexon"], "hmms": ["Adeno_hexon"]}],
        )
        assert "hexon" not in out

    def test_multidomain_token_requires_order(self):
        # A--B satisfied only when A is N-terminal to B on the protein.
        good = _hrec("G", 800, [
            _hit("DomA", ali_from=10, ali_to=200),
            _hit("DomB", ali_from=300, ali_to=600),
        ])
        bad = _hrec("Bd", 800, [
            _hit("DomB", ali_from=10, ali_to=200),
            _hit("DomA", ali_from=300, ali_to=600),
        ])
        spec = [{"name": "fusion", "hmms": ["DomA--DomB"]}]
        assert self.ident.identify([good], spec)["fusion"].id == "G"
        assert "fusion" not in self.ident.identify([bad], spec)

    def test_or_tokens_alternative_architectures(self):
        # Two tokens in one spec are OR'd: either architecture satisfies.
        alpha = _hrec("alpha", 600, [_hit("CoV_S1"), _hit("CoV_S2",
                                                          ali_from=300, ali_to=590)])
        spec = [{"name": "spike", "hmms": ["CoV_S1--CoV_S2", "bCoV_RBD--CoV_S2"]}]
        out = self.ident.identify([alpha], spec)
        assert out["spike"].id == "alpha"

    def test_spec_without_hmms_uses_name_fallback(self):
        rec = _hrec("B1", 300, [], description="DNA polymerase catalytic subunit")
        out = self.ident.identify(
            [rec], [{"name": "pol", "aliases": ["DNA polymerase"]}],
        )
        assert out["pol"].id == "B1"

    def test_empty_proteins_returns_empty(self):
        assert self.ident.identify([], [{"name": "x", "hmms": ["X"]}]) == {}


# ---------------------------------------------------------------------------
# crop_to_marker_domain — keep only the marker domain of a long protein
# (e.g. the RdRp of a viral polyprotein), with metadata intact.
# ---------------------------------------------------------------------------

from vfam_trees.markers import crop_to_marker_domain


def _polyprotein(n=1000, hits=None):
    rec = SeqRecord(Seq("M" + "ACDEFGHIKL" * ((n - 1) // 10) + "A" * ((n - 1) % 10)),
                    id="YP_000001.1", name="YP_000001",
                    description="polyprotein [Testvirus alpha]")
    rec.annotations = {"molecule_type": "protein", "organism": "Testvirus alpha",
                       "taxonomy": ["Viruses", "Riboviria"], "db_source": "REFSEQ: accession NC_000001.1",
                       "hmm_hits": hits or []}
    rec.features = [
        SeqFeature(FeatureLocation(0, n), type="source",
                   qualifiers={"organism": ["Testvirus alpha"], "host": ["Homo sapiens"],
                               "collection_date": ["2019"], "geo_loc_name": ["Peru: Lima"],
                               "db_xref": ["taxon:12345"]}),
        SeqFeature(FeatureLocation(0, n), type="Protein", qualifiers={"product": ["polyprotein"]}),
        SeqFeature(FeatureLocation(600, 900), type="Region", qualifiers={"region_name": ["RdRP_1"]}),
        SeqFeature(FeatureLocation(0, n), type="CDS",
                   qualifiers={"coded_by": ["NC_000001.1:100..3102"]}),
    ]
    return rec


def _crop_hit(target, a, b, passing=True):
    return {"target": target, "ali_from": a, "ali_to": b, "passing": passing}


class TestCropToMarkerDomain:
    def test_crops_to_hit_span_plus_padding(self):
        rec = _polyprotein(1000, [_crop_hit("RdRP_1", 601, 900)])
        out = crop_to_marker_domain(rec, ["RdRP_1"], 50)
        assert len(out.seq) == 400                      # residues 551–950
        assert str(out.seq) == str(rec.seq[550:950])
        assert "residues 551-950 of 1000" in out.annotations["comment"]

    def test_span_covers_all_matching_hits_and_clamps(self):
        rec = _polyprotein(1000, [_crop_hit("Bunya_RdRp", 20, 400), _crop_hit("L_thumb_ring_vir", 700, 990)])
        out = crop_to_marker_domain(rec, ["Bunya_RdRp", "L_thumb_ring_vir"], 50)
        assert len(out.seq) == 1000 and out is rec      # 0..1000 after clamping: unchanged

    def test_non_passing_and_other_profiles_are_ignored(self):
        rec = _polyprotein(1000, [_crop_hit("RdRP_1", 601, 900, passing=False), _crop_hit("Other", 1, 100)])
        assert crop_to_marker_domain(rec, ["RdRP_1"], 50) is rec

    def test_distant_stray_hit_left_out(self):
        # The real palm hit near residue 1800 plus a weak stray one near 300:
        # crop around the best hit only, not most of the polyprotein.
        rec = _polyprotein(2200, [
            {**_crop_hit("RdRP_1", 1701, 2000), "dom_score": 300.0},
            {**_crop_hit("RdRP_1", 301, 380), "dom_score": 25.0},
        ])
        out = crop_to_marker_domain(rec, ["RdRP_1"], 100)
        assert str(out.seq) == str(rec.seq[1600:2100])
        assert "1 distant hit(s) left out" in out.annotations["comment"]

    def test_second_domain_copy_left_out(self):
        rec = _polyprotein(3000, [
            {**_crop_hit("RdRP_1", 401, 800), "dom_score": 200.0},
            {**_crop_hit("RdRP_1", 2001, 2400), "dom_score": 350.0},
        ])
        out = crop_to_marker_domain(rec, ["RdRP_1"], 0)
        assert str(out.seq) == str(rec.seq[2000:2400])

    def test_neighbouring_domains_stay_together(self):
        # Bunyavirus L: PA-C-like and RdRp domains ~380 aa apart are one region.
        rec = _polyprotein(2000, [
            {**_crop_hit("L_PA-C-like", 101, 400), "dom_score": 90.0},
            {**_crop_hit("Bunya_RdRp", 781, 1500), "dom_score": 800.0},
        ])
        out = crop_to_marker_domain(rec, ["Bunya_RdRp", "L_PA-C-like"], 0)
        assert str(out.seq) == str(rec.seq[100:1500])

    def test_multidomain_tokens_contribute_their_names(self):
        rec = _polyprotein(1000, [_crop_hit("A", 301, 400), _crop_hit("B", 501, 600)])
        out = crop_to_marker_domain(rec, ["A--B"], 0)
        assert len(out.seq) == 300

    def test_metadata_survives_a_genbank_round_trip(self, tmp_path):
        from Bio import SeqIO
        from vfam_trees.fetch import _source_nuc_accession, extract_metadata
        rec = _polyprotein(1000, [_crop_hit("RdRP_1", 601, 900)])
        out = crop_to_marker_domain(rec, ["RdRP_1"], 50)
        assert "hmm_hits" not in out.annotations
        assert {f.type for f in out.features} == {"source", "CDS"}
        path = tmp_path / "x.gb"
        with open(path, "w") as fh:
            SeqIO.write([out], fh, "genbank")
        back = next(SeqIO.parse(str(path), "genbank"))
        assert len(back.seq) == 400
        meta = extract_metadata(back)
        assert meta["species"] == "Testvirus alpha"
        assert meta["host"] == "Homo sapiens"
        assert meta["location"] == "Peru: Lima"
        assert meta["collection_date"] == "2019"
        assert meta["taxon_id"] == "12345"
        assert meta["length"] == 400
        assert _source_nuc_accession(back) == "NC_000001.1"


# ---------------------------------------------------------------------------
# hmm.best_hit_only — a marker protein must fit one of its spec's profiles
# better than any other profile in the database (foreign-RdRp filter).
# ---------------------------------------------------------------------------

from vfam_trees.markers import best_passing_target, foreign_best_hit, token_profile_names


def _scored(target, score, *, passing=True, dom_score=None, ali_from=1, ali_to=300):
    return {"target": target, "passing": passing, "full_score": score,
            "dom_score": score if dom_score is None else dom_score,
            "dom_evalue": 1e-30, "ali_from": ali_from, "ali_to": ali_to}


class TestBestPassingTarget:
    def test_highest_full_score_wins(self):
        assert best_passing_target([_scored("RdRP_1", 60), _scored("Orbi_VP1", 2200)]) == "Orbi_VP1"

    def test_non_passing_hits_are_ignored(self):
        hits = [_scored("Orbi_VP1", 2200, passing=False), _scored("RdRP_1", 60)]
        assert best_passing_target(hits) == "RdRP_1"

    def test_domain_score_breaks_ties(self):
        hits = [_scored("A", 100, dom_score=40), _scored("B", 100, dom_score=90)]
        assert best_passing_target(hits) == "B"

    def test_none_without_passing_hit(self):
        assert best_passing_target([]) is None
        assert best_passing_target([_scored("A", 10, passing=False)]) is None


class TestTokenProfileNames:
    def test_splits_multidomain_and_skips_malformed(self):
        assert token_profile_names(["A--B", "C", "D--"]) == {"A", "B", "C"}
        assert token_profile_names([]) == set()


class TestForeignBestHit:
    def test_foreign_profile_reported(self):
        rec = _hrec("L", 2257, [_scored("Bunya_RdRp", 790), _scored("L_thumb_ring_vir", 120)])
        assert foreign_best_hit(rec, {"RdRP_1", "RdRP_2"}) == "Bunya_RdRp"

    def test_own_profile_is_not_foreign(self):
        rec = _hrec("nsP4", 2500, [_scored("RdRP_2", 500), _scored("RdRP_1", 90)])
        assert foreign_best_hit(rec, {"RdRP_1", "RdRP_2"}) is None

    def test_unscanned_or_hitless_record_is_not_foreign(self):
        rec = SeqRecord(Seq("M" * 50), id="u")
        assert foreign_best_hit(rec, {"RdRP_1"}) is None
        assert foreign_best_hit(_hrec("n", 50, []), {"RdRP_1"}) is None


class TestBestHitProfiles:
    def test_tokens_split_and_best_hit_also_added(self):
        ident = HMMIdentifier({"hmm": {"enabled": True, "best_hit_also": ["L_thumb_ring_vir"]}})
        assert ident.best_hit_profiles(["A--B", "Bunya_RdRp"]) == {
            "A", "B", "Bunya_RdRp", "L_thumb_ring_vir"}

    def test_sibling_domain_best_hit_is_kept(self):
        # A divergent hantavirus L whose best full-sequence score is the L
        # thumb-ring domain: the family's own kind, not a foreign virus.
        spec = [{"name": "RdRp", "hmms": ["Bunya_RdRp", "Flu_PB1"]}]
        rec = _hrec("L", 2150, [_scored("L_thumb_ring_vir", 420), _scored("Bunya_RdRp", 380)])
        with_also = HMMIdentifier({"hmm": {"enabled": True, "best_hit_only": True,
                                           "best_hit_also": ["L_PA-C-like", "L_thumb_ring_vir"]}})
        without = HMMIdentifier({"hmm": {"enabled": True, "best_hit_only": True}})
        assert with_also.identify([rec], spec)["RdRp"].id == "L"
        assert "RdRp" not in without.identify([rec], spec)
        assert foreign_best_hit(rec, with_also.best_hit_profiles(spec[0]["hmms"])) is None


class TestIdentifierSettings:
    def test_min_protein_length_exposed(self):
        assert HMMIdentifier({"hmm": {"enabled": True, "min_protein_length": 100}}).min_protein_length == 100
        assert HMMIdentifier({"hmm": {"enabled": True}}).min_protein_length is None

    def test_run_threads_used_when_hmm_threads_null(self):
        ident = make_identifier({"hmm": {"enabled": True, "threads": None}}, threads=8)
        assert ident._cfg["hmm"]["threads"] == 8

    def test_explicit_hmm_threads_win(self):
        cfg = {"hmm": {"enabled": True, "threads": 2}}
        assert make_identifier(cfg, threads=8)._cfg["hmm"]["threads"] == 2
        assert cfg["hmm"]["threads"] == 2   # caller's dict untouched


class TestHMMIdentifierBestHitOnly:
    spec = [{"name": "RdRp", "hmms": ["RdRP_1", "RdRP_2", "RdRP_3", "RdRP_4"]}]
    strict = HMMIdentifier({"hmm": {"enabled": True, "best_hit_only": True}})
    loose = HMMIdentifier({"hmm": {"enabled": True}})

    def test_flag_exposed(self):
        assert self.strict.best_hit_only is True
        assert self.loose.best_hit_only is False

    def test_cross_hitting_foreign_protein_rejected_only_when_on(self):
        # An orbivirus VP1 filed under an alphavirus species: weak RdRP_1
        # cross-hit, but Orbi_VP1 fits far better.
        vp1 = _hrec("VP1", 1295, [_scored("Orbi_VP1", 2290), _scored("RdRP_1", 55)])
        assert "RdRp" not in self.strict.identify([vp1], self.spec)
        assert self.loose.identify([vp1], self.spec)["RdRp"].id == "VP1"

    def test_own_rdrp_kept_when_shorter_than_a_foreign_one(self):
        # Longest-satisfying applies among the proteins that pass the
        # best-hit check, not among all token-satisfying proteins.
        foreign = _hrec("F", 2200, [_scored("Bunya_RdRp", 700), _scored("RdRP_4", 40)])
        own = _hrec("O", 1500, [_scored("RdRP_2", 480), _scored("RdRP_1", 70)])
        assert self.strict.identify([foreign, own], self.spec)["RdRp"].id == "O"
        assert self.loose.identify([foreign, own], self.spec)["RdRp"].id == "F"


class TestReviewRound3:
    """Code review 2026-10-21 (second v1.6.0 round): shared hit test, GA cutoffs
    parsed once and only when used, profile-set fingerprint."""

    def test_marker_hits(self):
        from vfam_trees.markers import marker_hits
        rec = _hrec("P", 900, [_hit("Bunya_RdRp"), _hit("L_thumb_ring_vir"),
                               _hit("Bunya_RdRp", passing=False)])
        assert [h["target"] for h in marker_hits(rec, ["Bunya_RdRp"])] == ["Bunya_RdRp"]
        assert marker_hits(_hrec("Q", 900, []), ["Bunya_RdRp"]) == []
        assert marker_hits(SeqRecord(Seq("M" * 10), id="R"), ["Bunya_RdRp"]) == []

    def _unscanned(self, n=2):
        return [SeqRecord(Seq("M" * 50), id=f"U{i}") for i in range(n)]

    def test_ga_cutoffs_skipped_when_unused(self, monkeypatch):
        import vfam_trees.markers as m
        calls = []
        monkeypatch.setattr(m, "get_ga_cutoffs", lambda cfg: calls.append(1) or {})
        monkeypatch.setattr(m, "scan_proteins", lambda q, cfg: {})
        ident = HMMIdentifier({"hmm": {"enabled": True, "use_ga_when_available": False}})
        ident.prescan(self._unscanned())
        ident.prescan(self._unscanned())
        assert calls == []

    def test_ga_cutoffs_parsed_once_per_identifier(self, monkeypatch):
        import vfam_trees.markers as m
        calls = []
        monkeypatch.setattr(m, "get_ga_cutoffs", lambda cfg: calls.append(1) or {"X": 25.0})
        monkeypatch.setattr(m, "scan_proteins", lambda q, cfg: {})
        ident = HMMIdentifier({"hmm": {"enabled": True, "use_ga_when_available": True}})
        ident.prescan(self._unscanned())
        ident.prescan(self._unscanned())
        ident.prescan([_hrec("done", 10, [])])   # nothing to scan: no parse needed
        assert calls == [1]

    def test_database_fingerprint_follows_content(self, tmp_path):
        a = tmp_path / "a"; a.mkdir()
        b = tmp_path / "b"; b.mkdir()
        (a / "P1.hmm").write_text("HMMER3/f\nNAME  P1\n//\n")
        (b / "P1.hmm").write_text("HMMER3/f\nNAME  P1\n//\n")
        cfg = lambda d: {"hmm": {"enabled": True, "database": str(d)},
                         "cache_dir": str(tmp_path / "cache")}
        fa = HMMIdentifier(cfg(a)).database_fingerprint()
        assert len(fa) == 16
        assert HMMIdentifier(cfg(b)).database_fingerprint() == fa   # same content, other place
        (b / "P2.hmm").write_text("HMMER3/f\nNAME  P2\n//\n")          # a profile added
        assert HMMIdentifier(cfg(b)).database_fingerprint() != fa

    def test_bundled_rdrp_fingerprint_is_stable(self):
        from vfam_trees.config import rdrp_preset
        cfg = {"hmm": rdrp_preset("Togaviridae")["hmm"]}
        one = HMMIdentifier(cfg).database_fingerprint()
        assert one == HMMIdentifier(cfg).database_fingerprint()
