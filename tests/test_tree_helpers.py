"""Tests for vfam_trees.tree helpers (no external tool calls)."""
from pathlib import Path

import pytest

from vfam_trees.tree import (
    TREE_INFO_SUFFIX,
    _fasttree_aa_flags,
    _fasttree_nuc_flags,
    is_model_finder_spec,
    keep_support_value,
    kept_support,
    parse_iqtree_best_model,
    parse_iqtree_partition_models,
    read_tree_info,
    strip_ufboot_flags,
    support_measures,
    support_type_for,
    validate_newick,
    with_default_support,
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

    def test_iqtree_with_ufboot_and_alrt_carries_ufboot(self):
        # v1.5.0 tree_100 default: IQ-TREE writes "alrt/ufboot" labels and the
        # wrapper keeps the UFBoot value, so that is what the outputs carry.
        assert support_type_for("iqtree", "-B 1000 -alrt 1000") == "UFBoot"
        assert support_type_for("iqtree", "-alrt 1000 -B 1000 --bnni") == "UFBoot"

    def test_iqtree_with_alrt_and_standard_bootstrap(self):
        assert support_type_for("iqtree", "-alrt 1000 -b 100") == "Bootstrap"

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



# ---------------------------------------------------------------------------
# Support-flag classifier: every IQ-TREE spelling, IQ-TREE's label order.
# ---------------------------------------------------------------------------

class TestSupportClassifier:
    def test_long_and_legacy_spellings(self):
        assert support_measures("--alrt 1000 -B 1000".split()) == ["alrt", "ufboot"]
        assert support_measures("-alrt 1000 --ufboot 1000".split()) == ["alrt", "ufboot"]
        assert support_measures("-bb 1000".split()) == ["ufboot"]
        assert support_measures("--boot 100".split()) == ["bootstrap"]
        assert support_measures("--abayes --lbp 1000".split()) == ["lbp", "abayes"]

    def test_label_order_matches_iqtree(self):
        # IQ-TREE 2.3.2: "SH-aLRT / local bootstrap / aBayes / ultrafast bootstrap"
        parts = "-B 1000 --abayes --lbp 1000 --alrt 1000".split()
        assert support_measures(parts) == ["alrt", "lbp", "abayes", "ufboot"]

    def test_default_alrt_added_only_without_alrt_or_bootstrap(self):
        assert with_default_support(["--fast"]) == ["-alrt", "1000", "--fast"]
        assert with_default_support(["--abayes"]) == ["-alrt", "1000", "--abayes"]
        assert with_default_support(["--alrt", "1000"]) == ["--alrt", "1000"]
        assert with_default_support(["-B", "1000"]) == ["-B", "1000"]
        assert with_default_support(["-b", "100"]) == ["-b", "100"]

    def test_kept_support_prefers_bootstrap_then_alrt(self):
        assert kept_support("-alrt 1000 -B 1000 --bnni".split()) == ("ufboot", 1, 2)
        assert kept_support("-alrt 1000 --abayes -B 1000".split()) == ("ufboot", 2, 3)
        assert kept_support("-alrt 1000 --abayes".split()) == ("alrt", 0, 2)
        assert kept_support(["-alrt", "1000"]) == ("alrt", 0, 1)

    def test_support_type_for_aliases(self):
        assert support_type_for("iqtree", "--alrt 1000 -B 1000") == "UFBoot"
        assert support_type_for("iqtree", "-alrt 1000 --abayes -B 1000") == "UFBoot"
        assert support_type_for("iqtree", "--abayes") == "SH_aLRT"
        assert support_type_for("iqtree", "--boot 100") == "Bootstrap"


# ---------------------------------------------------------------------------
# keep_support_value — IQ-TREE's slash-joined node labels reduced to one value
# ---------------------------------------------------------------------------

class TestKeepSupportValue:
    def test_two_values_keep_bootstrap(self):
        nwk = "((A:0.1,B:0.2)96.2/94:0.05,(C:0.3,D:0.4)100/100:0.01)88.5/72;"
        assert keep_support_value(nwk, 1, 2) == (
            "((A:0.1,B:0.2)94:0.05,(C:0.3,D:0.4)100:0.01)72;"
        )

    def test_three_values_from_alrt_abayes_ufboot(self):
        # IQ-TREE 2.3.2 with -alrt 1000 -abayes -B 1000 writes ")64.4/0.595/61:".
        nwk = "((A:0.1,B:0.2)64.4/0.595/61:0.05,C:0.3);"
        assert keep_support_value(nwk, 2, 3) == "((A:0.1,B:0.2)61:0.05,C:0.3);"

    def test_four_values_and_first_kept(self):
        nwk = "((A:0.1,B:0.2)100/100/1/100:0.05,C:0.3);"
        assert keep_support_value(nwk, 3, 4) == "((A:0.1,B:0.2)100:0.05,C:0.3);"
        assert keep_support_value("((A:1,B:1)64.4/0.595:1,C:1);", 0, 2) == "((A:1,B:1)64.4:1,C:1);"

    def test_single_labels_unchanged(self):
        nwk = "((A:0.1,B:0.2)94:0.05,(C:0.3,D:0.4)100:0.01);"
        assert keep_support_value(nwk, 1, 2) == nwk

    def test_leaf_names_with_slash_untouched(self):
        nwk = "((A/1:0.1,B:0.2)96.2/94:0.05,C:0.3);"
        assert keep_support_value(nwk, 1, 2) == "((A/1:0.1,B:0.2)94:0.05,C:0.3);"

    def test_unexpected_value_count_falls_back_by_position(self):
        nwk = "((A:1,B:1)10/20/30:1,C:1);"
        assert keep_support_value(nwk, 1, 2) == "((A:1,B:1)30:1,C:1);"   # last kept
        assert keep_support_value(nwk, 0, 2) == "((A:1,B:1)10:1,C:1);"   # first kept

    def test_result_parses_with_confidences(self, tmp_path):
        from Bio import Phylo
        nwk = "((A:0.1,B:0.2)64.4/0.595/61:0.05,(C:0.3,D:0.4)100/1/100:0.01);"
        path = tmp_path / "t.nwk"
        path.write_text(keep_support_value(nwk, 2, 3))
        tree = Phylo.read(str(path), "newick")
        confs = sorted(c.confidence for c in tree.find_clades() if c.confidence is not None)
        assert confs == [61.0, 100.0]


# ---------------------------------------------------------------------------
# strip_ufboot_flags — options for the fallback when UFBoot aborts
# ---------------------------------------------------------------------------

class TestStripUfbootFlags:
    def test_drops_ufboot_and_bnni_keeps_alrt(self):
        assert strip_ufboot_flags("-alrt 1000 -B 1000 --bnni".split()) == ["-alrt", "1000"]

    def test_adds_alrt_when_only_ufboot_was_requested(self):
        for flag in ("-B", "-bb", "--ufboot"):
            assert strip_ufboot_flags([flag, "1000"]) == ["-alrt", "1000"]

    def test_keeps_unrelated_options_and_other_measures(self):
        assert strip_ufboot_flags("-B 1000 -ntmax 4".split()) == ["-alrt", "1000", "-ntmax", "4"]
        assert strip_ufboot_flags("--alrt 1000 --abayes -B 1000".split()) == ["--alrt", "1000", "--abayes"]

    def test_standard_bootstrap_is_not_touched(self):
        assert strip_ufboot_flags(["-b", "100"]) == ["-b", "100"]

    def test_fallback_options_yield_sh_alrt_label(self):
        fallback = strip_ufboot_flags("-alrt 1000 -B 1000 --bnni".split())
        assert support_type_for("iqtree", " ".join(fallback)) == "SH_aLRT"


# ---------------------------------------------------------------------------
# Run record beside each tree, and the UFBoot-abort fallback (no real tools:
# subprocess.run is replaced by a scripted IQ-TREE).
# ---------------------------------------------------------------------------

class _Proc:
    def __init__(self, returncode=0, stderr=""):
        self.returncode, self.stdout, self.stderr = returncode, "", stderr


_ABORT = ("ERROR: hashsplitset.cpp:60: void SplitIntMap::insertSplit(Split *, int): "
          "Assertion `!findSplit(sp)' failed.\n*** IQ-TREE CRASHES WITH SIGNAL ABORTED\n")
_ML_TREE = "((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.01);"


def _arg(cmd, flag):
    return cmd[cmd.index(flag) + 1] if flag in cmd else None


class _ScriptedIqtree:
    """Plays IQ-TREE: a UFBoot run aborts after writing the ML tree."""

    def __init__(self, abort=True, log_text="", scheme=False, leave_treefile=True):
        self.calls, self.abort, self.log_text = [], abort, log_text
        self.scheme, self.leave_treefile = scheme, leave_treefile

    def __call__(self, cmd, *a, **k):
        cmd = [str(c) for c in cmd]
        self.calls.append(cmd)
        prefix = _arg(cmd, "--prefix")
        ufboot = "-B" in cmd
        if ufboot and self.abort:
            if self.leave_treefile:
                Path(prefix + ".treefile").write_text(_ML_TREE)
            Path(prefix + ".log").write_text(self.log_text)
            if self.scheme:
                Path(prefix + ".best_scheme.nex").write_text(
                    "#nexus\nbegin sets;\n  charpartition mymodels = LG+G4: m1, WAG+G4: m2;\nend;\n")
            return _Proc(-6, _ABORT)
        labels = "100/98" if ufboot else "97"
        Path(prefix + ".treefile").write_text(
            f"((A:0.1,B:0.2){labels}:0.05,(C:0.3,D:0.4){labels}:0.01);")
        return _Proc(0)


def _run(monkeypatch, tmp_path, fake, options="-alrt 1000 -B 1000 --bnni", model="GTR+G"):
    from vfam_trees import tree as tree_mod
    monkeypatch.setattr(tree_mod.subprocess, "run", fake)
    aln = tmp_path / "aln.fasta"
    aln.write_text(">A\nATGC\n>B\nATGC\n>C\nATGA\n>D\nATGG\n")
    nwk = tree_mod.run_tree(aln, tmp_path, "tree_100", tool="iqtree", seq_type="nucleotide",
                            model_nuc=model, options=options)
    return nwk, read_tree_info(tmp_path, "tree_100")


class TestTreeInfo:
    def test_missing_or_bad_record_is_none(self, tmp_path):
        assert read_tree_info(tmp_path, "tree_100") is None
        (tmp_path / f"tree_100{TREE_INFO_SUFFIX}").write_text("not json")
        assert read_tree_info(tmp_path, "tree_100") is None

    def test_normal_ufboot_run(self, tmp_path, monkeypatch):
        nwk, info = _run(monkeypatch, tmp_path, _ScriptedIqtree(abort=False))
        assert info == {"tool": "iqtree", "options": "-alrt 1000 -B 1000 --bnni",
                        "support_type": "UFBoot", "note": ""}
        assert "100/98" not in nwk.read_text() and ")98:" in nwk.read_text()

    def test_fasttree_run_replaces_a_stale_iqtree_record(self, tmp_path, monkeypatch):
        from vfam_trees import tree as tree_mod
        (tmp_path / f"tree_500{TREE_INFO_SUFFIX}").write_text(
            '{"tool": "iqtree", "options": "--fast", "support_type": "SH_aLRT", "note": ""}')

        def fake_fasttree(alignment_fasta, output_dir, prefix, *a):
            out = output_dir / f"{prefix}.nwk"
            out.write_text("((A:0.1,B:0.2)0.97:0.05,C:0.3);")
            return out

        monkeypatch.setattr(tree_mod, "_run_fasttree", fake_fasttree)
        tree_mod.run_tree(tmp_path / "aln.fasta", tmp_path, "tree_500", tool="fasttree",
                          seq_type="nucleotide")
        assert read_tree_info(tmp_path, "tree_500")["support_type"] == "SH_like"

    def test_failed_run_leaves_no_record(self, tmp_path, monkeypatch):
        from vfam_trees import tree as tree_mod
        (tmp_path / f"tree_100{TREE_INFO_SUFFIX}").write_text(
            '{"tool": "iqtree", "options": "", "support_type": "UFBoot", "note": ""}')
        monkeypatch.setattr(tree_mod.subprocess, "run", lambda *a, **k: _Proc(2, "bad option"))
        aln = tmp_path / "aln.fasta"
        aln.write_text(">A\nAC\n>B\nAC\n")
        with pytest.raises(RuntimeError):
            tree_mod.run_tree(aln, tmp_path, "tree_100", tool="iqtree", seq_type="nucleotide")
        assert read_tree_info(tmp_path, "tree_100") is None


class TestUfbootAbortFallback:
    def test_rescores_fixed_ml_tree_with_alrt(self, tmp_path, monkeypatch):
        fake = _ScriptedIqtree()
        nwk, info = _run(monkeypatch, tmp_path, fake)
        assert len(fake.calls) == 2
        rescore = fake.calls[1]
        assert "--tree-fix" in rescore and _arg(rescore, "-t").endswith(".ml_tree.nwk")
        assert "-B" not in rescore and "--bnni" not in rescore and "-alrt" in rescore
        assert _arg(rescore, "-m") == "GTR+G"
        assert info["support_type"] == "SH_aLRT"
        assert info["options"] == "-alrt 1000"
        assert "fixed topology" in info["note"]
        assert ")97:" in nwk.read_text()

    def test_modelfinder_pick_read_from_the_aborted_log(self, tmp_path, monkeypatch):
        fake = _ScriptedIqtree(log_text="Best-fit model: GTR+F+I+G4 chosen according to BIC\n")
        _, info = _run(monkeypatch, tmp_path, fake, model="TEST")
        assert _arg(fake.calls[1], "-m") == "GTR+F+I+G4"
        assert info["support_type"] == "SH_aLRT"

    def test_partitioned_run_rescored_from_best_scheme_without_m(self, tmp_path, monkeypatch):
        parts = tmp_path / "partitions_100.nex"
        parts.write_text("#nexus\n")
        fake = _ScriptedIqtree(scheme=True)
        _, info = _run(monkeypatch, tmp_path, fake, model="MFP",
                       options=f"-p {parts} -alrt 1000 -B 1000 --bnni")
        rescore = fake.calls[1]
        assert "-m" not in rescore
        assert _arg(rescore, "-p").endswith(".fallback_scheme.nex")
        assert (tmp_path / "tree_100.best_scheme.nex").exists()
        assert info["options"] == f"-p {parts} -alrt 1000"

    def test_full_search_without_ufboot_when_no_ml_tree(self, tmp_path, monkeypatch):
        fake = _ScriptedIqtree(leave_treefile=False)
        _, info = _run(monkeypatch, tmp_path, fake)
        assert len(fake.calls) == 2
        rerun = fake.calls[1]
        assert "--tree-fix" not in rerun and "-B" not in rerun and "-alrt" in rerun
        assert info["support_type"] == "SH_aLRT"
        assert "search repeated" in info["note"]

    def test_modelfinder_without_logged_pick_repeats_search(self, tmp_path, monkeypatch):
        fake = _ScriptedIqtree(log_text="no model line\n")
        _run(monkeypatch, tmp_path, fake, model="MFP")
        assert "--tree-fix" not in fake.calls[1] and _arg(fake.calls[1], "-m") == "MFP"

    def test_ordinary_failure_is_not_retried(self, tmp_path, monkeypatch):
        calls = []

        def fake(cmd, *a, **k):
            calls.append(cmd)
            return _Proc(2, "ERROR: Unknown option -X")

        with pytest.raises(RuntimeError):
            _run(monkeypatch, tmp_path, fake)
        assert len(calls) == 1
