"""Tests for vfam_trees.tree helpers (no external tool calls)."""
from pathlib import Path

import pytest

from vfam_trees.tree import (
    _fasttree_aa_flags,
    _fasttree_nuc_flags,
    is_model_finder_spec,
    parse_iqtree_best_model,
    parse_iqtree_partition_models,
    support_type_for,
    validate_newick,
)


# ---------------------------------------------------------------------------
# is_model_finder_spec
# ---------------------------------------------------------------------------

class TestIsModelFinderSpec:
    def test_empty_is_model_finder(self):
        assert is_model_finder_spec("") is True

    def test_mf_is_model_finder(self):
        assert is_model_finder_spec("MF") is True

    def test_mfp_is_model_finder(self):
        assert is_model_finder_spec("MFP") is True

    def test_test_is_model_finder(self):
        assert is_model_finder_spec("TEST") is True
        assert is_model_finder_spec("TESTONLY") is True

    def test_concrete_model_not_model_finder(self):
        assert is_model_finder_spec("GTR+G") is False
        assert is_model_finder_spec("LG+G") is False
        assert is_model_finder_spec("JC69") is False

    def test_case_insensitive(self):
        assert is_model_finder_spec("mf") is True
        assert is_model_finder_spec("mfp") is True
        assert is_model_finder_spec("test") is True

    def test_model_with_rate_suffix(self):
        assert is_model_finder_spec("MFP+G") is True
        assert is_model_finder_spec("GTR+I+G") is False


# ---------------------------------------------------------------------------
# parse_iqtree_best_model
# ---------------------------------------------------------------------------

class TestParseIqtreeBestModel:
    def test_missing_file_returns_none(self, tmp_path):
        assert parse_iqtree_best_model(tmp_path / "nope.iqtree") is None

    def test_best_fit_line_preferred(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text(
            "Some header\n"
            "Best-fit model according to BIC: GTR+I+G4\n"
            "Other stuff\n"
            "Model of substitution: IGNORED\n"
        )
        assert parse_iqtree_best_model(log) == "GTR+I+G4"

    def test_fallback_to_model_of_substitution(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text("Model of substitution: LG+F+G\n")
        assert parse_iqtree_best_model(log) == "LG+F+G"

    def test_no_match_returns_none(self, tmp_path):
        log = tmp_path / "log.iqtree"
        log.write_text("some unrelated log content\n")
        assert parse_iqtree_best_model(log) is None


# ---------------------------------------------------------------------------
# _fasttree_nuc_flags / _fasttree_aa_flags
# ---------------------------------------------------------------------------

class TestFasttreeNucFlags:
    def test_gtr_adds_gtr_flag(self):
        assert _fasttree_nuc_flags("GTR") == ["-gtr"]

    def test_gtr_plus_gamma_adds_both(self):
        flags = _fasttree_nuc_flags("GTR+G")
        assert "-gtr" in flags and "-gamma" in flags

    def test_jc_gives_default_no_flags(self):
        assert _fasttree_nuc_flags("JC") == []
        assert _fasttree_nuc_flags("JC69") == []

    def test_jc_plus_gamma(self):
        flags = _fasttree_nuc_flags("JC+G")
        assert "-gamma" in flags
        assert "-gtr" not in flags

    def test_unsupported_falls_back_to_gtr(self):
        # HKY is not supported → warn + use GTR
        flags = _fasttree_nuc_flags("HKY")
        assert "-gtr" in flags

    def test_empty_model_defaults_to_gtr(self):
        assert "-gtr" in _fasttree_nuc_flags("")


class TestFasttreeAaFlags:
    def test_lg_adds_lg(self):
        assert _fasttree_aa_flags("LG") == ["-lg"]

    def test_lg_plus_gamma(self):
        flags = _fasttree_aa_flags("LG+G")
        assert "-lg" in flags and "-gamma" in flags

    def test_wag_adds_wag(self):
        assert _fasttree_aa_flags("WAG") == ["-wag"]

    def test_jtt_is_default_no_flag(self):
        assert _fasttree_aa_flags("JTT") == []

    def test_empty_model_defaults_to_jtt(self):
        assert _fasttree_aa_flags("") == []

    def test_unsupported_falls_back_to_wag(self):
        flags = _fasttree_aa_flags("DAYHOFF")
        assert "-wag" in flags


# ---------------------------------------------------------------------------
# validate_newick
# ---------------------------------------------------------------------------

class TestValidateNewick:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            validate_newick(tmp_path / "missing.nwk")

    def test_empty_file_raises(self, tmp_path):
        empty = tmp_path / "e.nwk"
        empty.write_text("")
        with pytest.raises(FileNotFoundError):
            validate_newick(empty)

    def test_invalid_newick_raises(self, tmp_path):
        # Unbalanced parentheses make the Newick parser raise.
        bad = tmp_path / "bad.nwk"
        bad.write_text("(A:0.1,B:0.2,(C:0.3")
        with pytest.raises(ValueError):
            validate_newick(bad)

    def test_valid_newick_passes(self, tmp_path):
        good = tmp_path / "good.nwk"
        good.write_text("(A:0.1,B:0.2,(C:0.3,D:0.4):0.05);")
        validate_newick(good)  # should not raise


# ---------------------------------------------------------------------------
# parse_iqtree_partition_models — reads sibling .best_scheme.nex
# ---------------------------------------------------------------------------

class TestParseIqtreePartitionModels:
    """Verifies the NEXUS charpartition parser used to extract per-partition
    models from IQ-TREE's -p output."""

    def _write_nex(self, tmp_path, name, content):
        path = tmp_path / f"{name}.iqtree"
        # Sibling that the parser actually reads
        nex = tmp_path / f"{name}.best_scheme.nex"
        nex.write_text(content)
        return path

    def test_basic_two_partitions(self, tmp_path):
        log = self._write_nex(tmp_path, "tree_100",
            """#nexus
begin sets;
charset DNA_polymerase = 1-1234;
charset MCP = 1235-2456;
charpartition mymodels = LG+I+G4: DNA_polymerase, WAG+G4: MCP;
end;
""")
        out = parse_iqtree_partition_models(log)
        assert out == {"DNA_polymerase": "LG+I+G4", "MCP": "WAG+G4"}

    def test_extra_whitespace_handled(self, tmp_path):
        log = self._write_nex(tmp_path, "tree_100",
            """#nexus
begin sets;
charset polB = 1-100;
charpartition  scheme   =    LG+I+G4 :  polB ,  WAG+G4 : MCP ;
end;
""")
        out = parse_iqtree_partition_models(log)
        assert out == {"polB": "LG+I+G4", "MCP": "WAG+G4"}

    def test_missing_nex_returns_empty(self, tmp_path):
        log = tmp_path / "tree_100.iqtree"
        log.write_text("just an iqtree report")
        # No best_scheme.nex sibling
        assert parse_iqtree_partition_models(log) == {}

    def test_no_charpartition_returns_empty(self, tmp_path):
        log = self._write_nex(tmp_path, "tree_100",
            """#nexus
begin sets;
charset polB = 1-100;
end;
""")
        assert parse_iqtree_partition_models(log) == {}

    def test_multi_marker_partitioned_run(self, tmp_path):
        log = self._write_nex(tmp_path, "tree_100",
            """#nexus
begin sets;
charset polB = 1-1000;
charset MCP = 1001-1500;
charset hel = 1501-2000;
charset ATPase = 2001-2300;
charpartition mymodels = LG+I+G4: polB, WAG+G4: MCP, JTT+G: hel, LG+F+I+G: ATPase;
end;
""")
        out = parse_iqtree_partition_models(log)
        assert out == {
            "polB": "LG+I+G4",
            "MCP": "WAG+G4",
            "hel": "JTT+G",
            "ATPase": "LG+F+I+G",
        }


# ---------------------------------------------------------------------------
# support_type_for — branch-support label inferred from tool + options.
# Regression for v1.2.12: concat-mode tree_100 was hardcoding "SH_aLRT" even
# when the actual IQ-TREE invocation used "-B 1000" (UFBoot).
# ---------------------------------------------------------------------------

class TestSupportTypeFor:
    def test_fasttree_is_sh_like(self):
        assert support_type_for("fasttree", "") == "SH_like"
        assert support_type_for("fasttree", "-gtr -gamma") == "SH_like"

    def test_iqtree_default_is_sh_alrt(self):
        # Wrapper auto-injects -alrt 1000 when no support flag is set.
        assert support_type_for("iqtree", "") == "SH_aLRT"
        assert support_type_for("iqtree", "--fast") == "SH_aLRT"

    def test_iqtree_with_ufboot_returns_ufboot(self):
        # `-B 1000` is the concat-mode protein default — must not be labelled
        # SH_aLRT in summary.tsv / PhyloXML.
        assert support_type_for("iqtree", "-B 1000") == "UFBoot"

    def test_iqtree_with_ufboot_and_alrt(self):
        assert support_type_for("iqtree", "-B 1000 -alrt 1000") == "SH_aLRT_UFBoot"

    def test_iqtree_with_standard_bootstrap(self):
        assert support_type_for("iqtree", "-b 100") == "Bootstrap"

    def test_iqtree2_alias(self):
        # tool name canonicalization: "iqtree", "iqtree2", "iq-tree" all OK.
        assert support_type_for("iqtree2", "") == "SH_aLRT"
        assert support_type_for("iq-tree", "") == "SH_aLRT"


# ---------------------------------------------------------------------------
# IQ-TREE invocation — must pass --seqtype DNA|AA explicitly so protein
# alignments with high A/C/G/T content are not misdetected as nucleotide.
# Regression for v1.2.12.
# ---------------------------------------------------------------------------

class TestIqtreeSeqtype:
    @staticmethod
    def _run_iqtree_and_capture(monkeypatch, tmp_path, seq_type):
        from vfam_trees import tree as tree_mod

        captured = {}

        class _Result:
            returncode = 0
            stdout = ""
            stderr = ""

        def fake_run(cmd, *a, **k):
            captured["cmd"] = list(cmd)
            # Synthesize the output the wrapper expects so it doesn't raise.
            i = cmd.index("--prefix")
            prefix = cmd[i + 1]
            Path(prefix).parent.mkdir(parents=True, exist_ok=True)
            Path(prefix + ".treefile").write_text("(A:1,B:1);")
            return _Result()

        monkeypatch.setattr(tree_mod.subprocess, "run", fake_run)
        aln = tmp_path / "aln.fasta"
        aln.write_text(">A\nATGC\n>B\nATGC\n")
        tree_mod._run_iqtree(
            alignment_fasta=aln,
            output_dir=tmp_path,
            prefix="tree_test",
            seq_type=seq_type,
            model_nuc="GTR+G",
            model_aa="LG+G",
            options="",
            threads=1,
        )
        return captured["cmd"]

    def test_nucleotide_passes_dna_seqtype(self, tmp_path, monkeypatch):
        cmd = self._run_iqtree_and_capture(monkeypatch, tmp_path, "nucleotide")
        i = cmd.index("--seqtype")
        assert cmd[i + 1] == "DNA"

    def test_protein_passes_aa_seqtype(self, tmp_path, monkeypatch):
        cmd = self._run_iqtree_and_capture(monkeypatch, tmp_path, "protein")
        i = cmd.index("--seqtype")
        assert cmd[i + 1] == "AA"

