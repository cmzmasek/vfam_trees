# Changelog

All notable user-facing changes are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The full per-commit history is in `git log`; this file summarizes only
changes that affect outputs, configuration, CLI behaviour, or other
user-visible surface area.

## [1.6.3] — 2026-10-28

Mimiviridae and Marseilleviridae change from name-based to HMM-identified
markers. Delete their family configs written by an earlier run so the new
presets apply.

### Changed
- **Mimiviridae and Marseilleviridae: markers identified by HMM.** Their
  concatenation presets found the marker proteins by GenBank name, and NCBI
  names the RNA polymerase subunits and the packaging ATPase in many ways, so
  many genomes fell short of `min_fraction`. Both now use six markers of the
  bundled `NCLDV` set: DNA polymerase B, RNA polymerase largest and second
  subunits, D5-like primase-helicase, mRNA capping enzyme and topoisomerase II
  (`relative_length_cutoff` 0.30, as for Asfarviridae; at 0.50 HMMER's split
  of the polymerase domain loses it in more than half of the Mimiviridae
  genomes). The major capsid protein and VLTF-3 are left out: Mimiviridae
  genomes carry 2–5 capsid paralogs and 18 of 68 two VLTF-3-like proteins,
  and the pipeline takes the longest hit per marker, which would mix
  paralogs. The packaging ATPase (Pfam `Pox_A32`) is not used either; it has
  no hit in the klosneuviruses and most marseilleviruses. Run end to end
  (2026-10-27; Mimiviridae checked again on 2026-10-28 with the two changes
  below), compared with the 2026-10-22 run:
  - Marseilleviridae: 37 genomes of 20 species, trees of 30 / 26 tips
    (before: 23 / 23 tips of 18 species); Noumeavirus and Golden
    Marseillevirus are new. The records of *Marseillevirus marseillevirus*
    and Golden Marseillevirus have no protein with a DNA polymerase B domain,
    so they go in with 5 of the 6 markers.
  - Mimiviridae: 69 genomes of 45 species, trees of 59 / 59 tips of 45
    species (before: 42 / 41 tips of 29 species). New: the four
    chlorarachniophyte viruses, Cotonvirus, Indivirus, *Catovirus
    klosterense*, six moumouviruses and four megaviruses. Lost: Megavirus
    Me1-2 (a draft genome with its markers on different contigs). All ICTV
    genera with two or more tips are monophyletic, *Tupanvirus* among them
    (it was not before), except that two records NCBI files as "Mimivirus
    sp." (MT663335, PQ010613) sit inside *Megavirus*. Megamimivirinae and
    Klosneuvirinae are monophyletic; Aliimimivirinae is too, apart from the
    unclassified Chlorella virus XW01, which sits beside the
    chlorarachniophyte viruses.
  - Signal `ok` on all four trees; PhyloXML valid.
  - Mimiviridae takes longer: all proteins of 151 species (~92,000) are
    downloaded instead of six name searches (~2 h), and the partitioned
    tree_100 search on ~60 genomes takes ~1.5 h; the family runs ~3.5–4 h
    (1 h before). Marseilleviridae: 27 min.
- **Mimiviridae: no branch-length outlier removal** (`outlier_removal.enabled:
  false` in its preset). Most Mimiviridae genomes are near-identical
  Megamimivirinae (median tip branch 0.005 subs/site), so the median + 20 × MAD
  threshold fell to ~0.10 and real members of the divergent Klosneuvirinae
  (0.2–0.5) were removed: Klosneuvirus KNV1, Catovirus CTV1, Indivirus ILV1,
  *Catovirus klosterense* and Chlorella virus XW01. Five RefSeq genomes of
  small groups (0.18–0.87) were over the threshold too and stayed only
  because RefSeq genomes are protected.
- **Mimiviridae: Megavirus haemonchus182742 left out** (`manual.exclude:
  [BK066332.1]` in its preset). The record is a TPA assembled from a
  *Haemonchus contortus* (nematode) genome project. With outlier removal
  off it entered the trees, did not group with *Megavirus*, and sat alone
  on a long branch (0.8 subs/site) at the base of Megamimivirinae.

## [1.6.2] — 2026-10-27

Two families of the 2026-10-22 run that had no tree get one, and a family
that finishes without a tree says so. Delete the Asfarviridae and
Pithoviridae configs written by an earlier run so the new presets apply.

### Changed
- **A family that finishes without a tree says so.** When each tree target
  ends with fewer than the 4 sequences a tree needs (after per-species
  clustering or length-outlier removal; in concatenated mode, genomes), the
  `status.tsv` value is `no tree: <reason>` instead of `OK`, e.g. `no tree:
  too few sequences after per-species clustering (2 < 4)`, the family log
  ends with a WARNING, and `vfam_trees status` shows `no tree`. In the
  2026-10-22 run 12 families were `OK` without a tree. A family with one of
  its two trees stays `OK`.

### Fixed
- **Asfarviridae and Pithoviridae had no tree.** Their concatenation presets
  found the marker proteins by GenBank name, and NCBI names most of them
  differently: African swine fever virus matched only its DNA polymerase and
  p72, so all 46 full-length genomes failed `min_fraction`; pithoviruses
  kept 2 genomes. Both now identify the markers by HMM, from a new bundled
  set `NCLDV` of nine Pfam profiles (accessions checked against InterPro):
  DNA polymerase B (`DNA_pol_B`), major capsid protein (`Capsid_NCLDV`,
  Asfarviridae only; pithoviruses have none Pfam detects), RNA polymerase
  largest and second subunits (`RNA_pol_Rpb1_2`, `RNA_pol_Rpb2_6`), D5-like
  primase-helicase (`D5_N`), VLTF-3 (`Pox_VLTF3`), mRNA capping enzyme
  (`mRNA_G-N7_MeTrfase` or `mRNA_cap_enzyme`) and topoisomerase II
  (`DNA_topoisoIV`). The profiles were chosen from the Pfam domains UniProt
  gives these viruses' proteins and tested on the pipeline's own download.
  `relative_length_cutoff` is 0.30 for both: HMMER splits their polymerase
  domain at an insertion into hits of 0.23–0.47 of the profile. Run end to
  end (2026-10-27): Asfarviridae keeps 35 genomes (33 ASFV, 2 Abalone
  asfa-like virus) with all 8 markers each, trees of 12 / 7 tips; Pithoviridae
  keeps the 18 complete genomes of 7 species, 13 / 13 tips, cedratviruses
  and pithoviruses in separate clades; signal `ok` on all four trees. The
  metagenome-assembled pithoviruses (LCPAC…) are split over contigs with at
  most 4 markers each and stay out.

## [1.6.1] — 2026-10-26

Two fixes found by the full re-run with 1.6.0. Mamonoviridae changes from
whole-genome nucleotide to a DNA polymerase tree; delete its family config
from an earlier run so the new preset applies.

### Fixed
- **A GenBank record without a sequence stopped the family.** NCBI files
  some WGS project records ("… whole genome shotgun sequencing project",
  a list of contig accessions) and CON-division records under virus taxa;
  Biopython parses them with an undefined sequence, and the translation
  fallback crashed on one (`UndefinedSequenceError`; Nodaviridae, from
  `CABPSU000000000.1` under Sander lucioperca betanodavirus, a pikeperch
  genome project). Such records are now left out with a log line in the
  fallback's download, in the translation, and wherever downloaded or
  `manual.include` records are read.
- **Mamonoviridae ran as whole-genome nucleotide.** Its presets were filed
  under *Medusaviridae*, a proposed name NCBI never used, so neither the
  DNA polymerase preset nor the NCLDV concatenation fallback applied, and
  MAFFT spent hours on four 360–560 kb genomes. Mamonoviridae now has its
  own preset: the family-B DNA polymerase identified by HMM (bundled
  `DNA_pol_B`, PF00136) and cut to the domain ± 100 aa. Two of the four
  genomes split the polymerase gene into two proteins (Medusavirus
  'euryale' between the exonuclease and polymerase parts, Furtivovirus B10
  inside the polymerase domain), and the spliced polymerases of
  *M. medusae* and *M. sthenus* carry ~680 N-terminal residues the others
  lack, so a name-based pick would align different pieces.
  Run end to end (2026-10-26): one polymerase per genome, cut to 404–663
  aa, a 4 × 663 alignment, signal `ok` on both trees, 2 minutes. A family
  config written by an earlier run pins whole-genome nucleotide; delete it
  so the preset applies.

## [1.6.0] — 2026-10-22

Phylogenetic-signal diagnostics for every tree, and 134 families moved from
whole-genome nucleotide alignment to an HMM-identified marker protein: the
RdRp for 107 RNA-virus families, the RT for 4 reverse-transcribing families
and Rep for 23 CRESS-DNA and satellite families. Also: a translation
fallback and a report for species without the marker, run records and
`init --full` for provenance, a 50-iteration search for protein tree_500,
outlier rounds checked on quick trees, trimAl `-gt 0.5`, and repeatable
MAFFT and IQ-TREE runs. Only templates and presets changed; per-family YAML
files written by earlier versions keep nucleotide analysis (a warning says
so), so use a fresh configs directory — and a fresh output directory, or
each species is downloaded again (see Fixed). A `global.yaml` that sets the
`msa_trim` or `tree_500` blocks keeps its old values; regenerate it with
`vfam_trees init --full` or delete those blocks.

### Added
- **Phylogenetic-signal check** (`signal.py`) on every final tree: median
  tip-to-tip distance (subs/site), median pairwise identity of the trimmed
  alignment, core fraction (residues a typical sequence keeps after trimming)
  and branches at IQ-TREE's 10 subs/site cap, with a verdict `ok` / `weak` /
  `saturated`. New `summary.tsv` columns `tree{500,100}_signal`,
  `_patristic_median`, `_identity_median`, `_core_frac`, `_branches_at_cap`;
  rows in the report PDF; an INFO line per tree and a WARNING with reasons
  when a tree is not `ok`. Nucleotide: saturated above a median of 1.5
  subs/site (or a branch at the cap with a median above 0.75), weak above
  0.75 or below 60 % identity or 25 % core. Protein: 4.0 / 3.0 / 30 %.
- **RdRp protein marker for 107 RNA-virus families** (`RDRP_FAMILIES`,
  kingdom Orthornavirae). The RdRp is identified by HMM against a bundled set
  of 23 Pfam RdRp profiles (`vfam_trees/data/hmms/RdRp/`, accessions verified
  against InterPro, recipe in `scripts/build_bundled_hmms.sh`), whatever
  GenBank calls it, and cut to the RdRp domain ± 100 residues. Hit rule for
  the set: domain E-value ≤ 1e-3, ≥ 20 % profile coverage, no Pfam GA.
- **Only the family's own kind of RdRp.** Each family's marker tokens are
  the profiles its own RefSeq RdRps hit (`RDRP_FAMILY_PROFILES`; the four
  broad palm-domain profiles RdRP_1–RdRP_4 always together), and the new
  `hmm.best_hit_only` (default `false`, `true` in the RdRp preset) counts a
  protein only when one of them fits it better than any other profile in the
  database. This keeps out other viruses filed under the family's taxa: NCBI
  had about 20 orbivirus and bunyavirus species (taxids 345xxxx, 2025
  deposits) under *Alphavirus*; their VP1 and L proteins fit Orbi_VP1 /
  Bunya_RdRp best and had joined the Togaviridae trees (25 of 92 tree_100
  tips and 15 of 162 tree_500 tips in the first validation run). Whole-genome
  nucleotide mode took them too. Each species affected gets a WARNING in the
  family log naming the profile and the accessions left out. Profiles for
  other domains of the same polymerase count as the family's own
  (`hmm.best_hit_also`, set by the preset from `RDRP_SIBLING_DOMAINS`): 6 of
  149 Peribunyaviridae RefSeq L proteins fit the L thumb-ring profile best,
  and Hantaviridae or Nairoviridae, whose RefSeqs never hit it, would
  otherwise lose a divergent L as "foreign".
- **RT protein marker for the reverse-transcribing families and Rep for the
  CRESS-DNA families** (`HMM_MARKER_FAMILIES`, `RT_FAMILY_PROFILES`,
  `REP_FAMILIES`; bundled Pfam sets `RT` and `Rep`, accessions verified
  against InterPro 2026-10-21). Same machinery, hit rule, length floor,
  opt-out and load-time checks as the RdRp. RefSeq coverage (same test):
  - RT (`RVT_1`): Caulimoviridae 108/108, Hepadnaviridae 27/27,
    Metaviridae 2/2, Retroviridae 64/89 — the 25 others have no RT gene
    (transforming viruses with an onc gene in place of *pol*, endogenous-
    retrovirus fragments) and are not in the RT tree. Cut to the RT domain
    ± 100 (~310–370 aa); `RVT_thumb` / `RVT_connect` count only for the
    best-hit rule, `RVT_2` (Ty1/copia) marks such RTs as foreign.
    Hepadnaviridae was whole-genome nucleotide.
  - Rep (`Viral_Rep`, `RNA_helicase`, `Gemini_AL1`, `Gemini_AL1_M`):
    ≥ 94 % of RefSeq viruses in 22 *Cressdnaviricota* families and
    Alphasatellitidae; whole proteins (~270–470 aa), every profile selects,
    no best-hit rule (the profiles separate geminivirus-like from
    circovirus-like Reps, not families). Circoviridae used a name-based Rep.
    Smacoviridae keeps it (Pfam: 38/88); Metaxyviridae (no RefSeq proteins),
    Tolecusatellitidae (no Rep) and Pseudoviridae (no RefSeq proteins) stay
    as they were.
  In the 2026-10-20 signal survey 22 of the 24 of these families with a
  result were saturated (17) or weak (5) as whole-genome nucleotide (tree_100
  median tip-to-tip distance 0.9–2.7 subs/site); the other two are small
  (Oomyviridae: 5 near-identical genomes; Redondoviridae: 23).
- **Translation fallback for HMM markers.** When no annotated protein of a
  species carries the marker, its nucleotide records (≤ 100, complete
  genomes first) are translated in six frames and run through the same
  profiles, best-hit rule and domain cut (`fetch.translated_orfs`,
  `_translated_marker_genomes`). Recovers genomes whose marker NCBI never
  annotated: on Retroviridae the galago prosimian foamy virus (12 kb genome,
  no pol CDS), bat spumavirus 1 (no proteins at all), the *Xenopus* ERV Xen1,
  and SHIV (70,422 proteins, the first 6000 all env/tat/rev/vpu clones).
  The record takes its genome's accession and metadata and carries a COMMENT.
- **Species without the marker are reported** in every pipeline mode: a
  WARNING in the family log, `<Family>_species_without_marker.tsv`, and
  `summary.tsv` columns `species_without_marker` / `species_marker_translated`
  (also in the report PDF). Before, a species with NCBI records but none of
  the marker or region was logged like one with no data at all ("No
  sequences found"). In concatenated mode a species whose genomes carry some
  of the markers, but none as many as `min_fraction` asks for, is listed as
  `too few markers` (Marseilleviridae: Noumeavirus, 916 NCBI proteins, has 3
  of the 6 name-matched markers). Costs one record-count query per such species outside
  HMM mode. In the whole-genome Retroviridae tree_100, 17 of 99 species had no
  RT protein; 5 of them come back by translation, 11 are acutely transforming
  or defective viruses with no pol gene.
- `sequence.domain_padding` (single-marker HMM mode; in the template, default
  `null`: whole protein): keep only the marker's hits ± that many residues.
  The span is the run of hits around the best-scoring one; a hit more than
  500 residues from its neighbours (`markers.MAX_DOMAIN_GAP`; the largest gap
  between domains of one RefSeq RdRp was 383) is left out, so a stray hit or
  a second copy in a chimeric record cannot bring back the polyprotein.
  Cropped records keep their metadata; a GenBank COMMENT records the
  residues used.
- `hmm.min_protein_length` (default `null`; 100 in the RdRp preset): proteins
  shorter than this are left out of the candidates in HMM mode. At 100 aa (the
  protein QC minimum) no tree changes. The floor is applied after download: an
  NCBI length clause (`[SLEN]`) would have saved 2–14 % of the records but
  made every protein search 4–6 s slower.
- **Run records.** Each `vfam_trees run` writes
  `<output>/run_records/<start time>.yaml`: version, command line, options,
  paths, Python / package / external-tool versions (MAFFT, trimAl, IQ-TREE,
  FastTree, mmseqs, CD-HIT, HMMER), status and exit code (filled in when the
  run ends), and every global setting the run used (built-in defaults
  overlaid with `global.yaml`; the NCBI API key is not recorded). The
  `<Family>.yaml` each family folder already received, its merged settings,
  now starts with a header naming the version, start time, thread count, the
  files merged and any keys the version did not read.
- **`vfam_trees init --full`** writes a global config listing every setting
  with its default and a comment (the plain template leaves out 29 of the 75
  per-family settings: the `hmm`, `output`, `manual`, `concatenation` and
  `refseq_absorption` blocks and a few keys). A test keeps it equal to the
  built-in defaults. `annotation_tsv`, `sequence.hmms`,
  `sequence.domain_padding` and `quality.max_length` are now in the plain
  template too; `sequence.hmms: []` is now in the built-in defaults (the
  presets set it; no behaviour change).

### Changed
- **tree_500 on protein alignments is a 50-iteration search (`-n 50`)
  instead of `--fast`.** On weak-signal protein sets `--fast` lands in a
  different local optimum on every seed: on 498 Tombusviridae RdRps its five
  seeds spanned −104,656 to −104,863 log-likelihood and each tree had 6–13
  splits with SH-aLRT ≥ 95 that a full search contradicts. `-n 50` (IQ-TREE's
  normal search, stopped after 50 iterations) reaches −104,483 with 4 such
  splits, as close to the best full search as a second full search is, at
  ~12× the `--fast` time (a full search is 35–75×). `--fast --runs 5` gains
  nothing. Nucleotide tree_500 keeps `--fast`. Benchmark below ("tree_500
  search effort"). `-n` is rejected under `tool: fasttree`, where it means
  something else (the number of alignments in the input).
- **Branch-length outlier rounds are checked on a quick tree first**
  (`outlier_removal.quick_check`, default true). Each round builds an
  IQ-TREE `--fast` tree (GTR+G / LG+G, no ModelFinder) and removes the
  outliers found on it; the configured search runs once a quick tree shows
  nothing more to remove, and its tree is checked exactly as before, so the
  final tree is always the configured one. Before, tree_100's UFBoot search
  (10–25 min for a protein family) repeated in every round — four times for
  Togaviridae — and tree_500's would have with `-n 50`. Measured on
  Togaviridae with its downloads reused: 92 → 52 min for the family; tree_100
  one UFBoot search after three quick trees instead of four searches (54 →
  12 min); tree_500 8 min (two `-n 50` searches, since the first configured
  tree still had two outliers; 1 min with `--fast` before). The final trees
  share 145 of 147 and 84 of 88 tips with the old ones, with the same signal
  verdicts and median support. Not used when the
  configured tree is FastTree or IQ-TREE `--fast` already. Both pipelines
  share the round logic (`branch_outliers.outlier_round` /
  `select_outliers`); the concatenated pipeline's protected-outlier and
  removal log lines now read like the single-marker ones.
- **Column trimming keeps columns with at most 50 % gaps (trimAl `-gt
  0.5`) instead of `-automated1`**, for nucleotide and protein alike. On
  alignments with many short PCR fragments `-automated1` switches to its
  strict mode and keeps few columns (Retroviridae RT tree_500: 58 of 657;
  Circoviridae 104 of 1660), and genera fall apart. On 22 protein and 20
  nucleotide alignments `-gt 0.5` grouped as many ICTV genera as
  `-automated1` in every one and more in six (Retroviridae tree_500 9 of 9
  instead of 4 of 9); tree inference takes ~1.2× as long for protein and
  ~1.4× for nucleotide. Benchmark below ("Trimming"). The signal check
  reads the trimmed alignment and its tree, so its verdicts move towards
  `saturated` (more variable columns are kept); the thresholds are
  unchanged. A trim that keeps no column now stops with an error naming
  `msa_trim.options`, instead of one blaming trimAl's output behaviour.
- HMMER (`hmmscan`, `hmmpress`) is required by the default family list and
  `vfam_trees test` reports it as such; `vfam_trees run --help` describes
  the RdRp default.
- Segmented RNA families on the RdRp preset no longer use a segment keyword:
  the HMM finds the RdRp on whichever segment carries it (and the log no
  longer announces one).
- **Opting out of the RdRp preset is all or nothing.** A family file that
  sets `sequence.type` other than protein, `sequence.region` other than
  `RdRp`, or `hmm.enabled: false` gets no part of the preset: it runs on the
  template under the file's own values (whole-genome nucleotide, with the
  usual segment for segmented families). Before, a file setting only
  `sequence.type: nucleotide` (e.g. one written before v1.3) kept the
  preset's HMM protein fetch and lost its segment: protein records then
  failed nucleotide QC, and for Hantaviridae every species was skipped. The
  load warning now names the settings that apply.
- **Half-edited HMM configs stop at load** with the lines to fix:
  `hmm.enabled: true` with `sequence.type: nucleotide` (e.g. a v1.6.0 RdRp
  file with only its type changed, or a pre-v1.3 Adenoviridae file under
  the hexon HMM preset) would fetch marker proteins and reject them all in
  nucleotide quality control, and a `sequence.segment` with an HMM protein
  marker would drop every record (protein titles name no segment). The
  opt-out warning now reports the segment that will really be used: a
  v1.6.0-generated file's `segment: null` wins over the usual segment, and
  the warning says that all segments then go into one alignment and which
  line to add.
- **The non-RdRp L domains no longer select or cut a protein.**
  `L_PA-C-like` and `L_thumb_ring_vir` were marker tokens for
  Peribunyaviridae, Phenuiviridae and Cruliviridae (their RefSeq L proteins
  hit them), so a partial L covering only such a domain would have been
  chosen and cut to it, with no RdRp columns. They now only count for the
  best-hit rule (`RDRP_NON_CORE_DOMAINS`). Every RefSeq L protein of the
  three families also hits `Bunya_RdRp`; their cut shrinks from a median of
  ~1300 to ~900 aa, the RdRp domain ± 100, as for Hantaviridae.
- `vfam_trees run` checks before Snakemake starts that the external tools
  the listed families need are on PATH (MAFFT, IQ-TREE, mmseqs, trimAl;
  HMMER for HMM families; FastTree when a tree block uses it). Before, a
  run without HMMER downloaded up to 6000 proteins for each RdRp family
  and then failed it, family after family.
- **Config keys this version does not read are reported** (WARNING), at any
  depth and for `global.yaml` as well as family files; before, only unknown
  top-level keys of family files were. Such keys had silently no effect: the
  `length_outlier.hi_mult` / `lo_mult` of v1.2.19 and earlier, for instance,
  survive in older `global.yaml` files and from there in every generated
  family file, while the filter uses `k`, `min_lo_mult` and `max_hi_mult`.
  The pipeline repeats the warning in each family log.
- `output_dir` and `log_level` are no longer in the `global.yaml` template:
  nothing ever read them (the output folder and log level are `vfam_trees
  run` options, `-o` and `--log-level`); an existing file gets a warning
  saying so.
- `vfam_trees run` reads and checks `global.yaml` before starting Snakemake,
  so a missing file or missing NCBI email stops the run at once instead of
  failing every family job.

### Fixed
- **A lost NCBI connection stalled a run forever.** Biopython opens its
  connections without a timeout, so a connection that died while the
  computer slept or moved to another network left the read waiting with no
  end: a validation run sat for two hours on one, at 0 % CPU. A connection
  that sends nothing for 180 s now fails (`fetch.NCBI_SOCKET_TIMEOUT`, set
  by `configure_entrez` unless the process already has a socket timeout),
  and the usual 5 retries apply. Sleep still pauses a run; keep the computer
  awake for long runs (macOS: `caffeinate -is vfam_trees run …`).
- **Repeated runs gave different alignments and trees.** MAFFT's threaded
  iterative refinement (FFT-NS-i, which `--auto` picks for a few hundred
  protein sequences, and L-INS-i, the tree_100 protein mode) returns a
  different alignment on every run: 100 Tombusviridae RdRps, `--thread 8`,
  three runs, three alignments of 1366–1400 columns. trimAl `-automated1`
  can magnify that: the same input kept 939 columns in one validation run and
  441 in the next. IQ-TREE drew a new random seed on every run, too: two
  `--fast` runs on the same Tombusviridae tree_500 alignment differed in 366
  of 495 splits. MAFFT now runs with `--threadit 0` (only the refinement on
  one thread; the alignment is then the same on every run and for any
  `--threads`, checked with 1, 2, 4 and 8), and IQ-TREE with `-seed 42` (the
  same tree and support values on every run at the same `--threads`,
  checked for `--fast` and for the tree_100 search with ModelFinder and
  UFBoot). Options that set `--threadit` or `-seed` themselves win. Cost:
  the refinement stage is no longer parallel, about 2× MAFFT time for
  L-INS-i (100 sequences: 12 → 22 s) and 4× for FFT-NS-i (500 sequences:
  28 → 110 s), a few minutes per protein family; the nucleotide modes have
  no refinement stage.
- **A tree could have more sequences than its target.** When many species
  have only one or two representatives, each is still given one slot, and the
  rounded shares then add up to more than the target. The correction took
  back at most one slot per species in a single pass, which was not enough:
  Genomoviridae (419 species, 686 representatives) got 550 sequences in
  tree_500. It now repeats until the total is exact, and never gives a
  species more slots than it has representatives (which left a tree short
  when rounding fell below the target). Selections that already hit the
  target are unchanged (checked on 1,200 random cases). Present since the
  proportional merge was introduced; tree_100 could overshoot the same way.
- **Nucleotide sequences in protein trees.** `additional_data_sources`
  (Pathoplexus genomes) are refused at config load for a family analysed as
  protein, with a message saying how to switch it to nucleotide; pasted
  sequences (`include_seq`, `include_fasta_files`) of the other molecule type
  are a hard error (`quality.sequence_type_mismatch`: ≥ 90 % A/C/G/T/U/N in a
  protein family, > 10 % non-nucleotide letters in a nucleotide one); and a
  `manual.include` accession of the other NCBI database is left out with a
  warning instead of being fetched. A/C/G/T/N are amino-acid letters too, so
  before this such a genome went into the RdRp alignment as a protected tip.
  Protein records given by hand are now also cut to the marker domain like
  fetched ones.
- **Downloads reused under other settings.** The per-species download kept in
  an output folder (`_work/species/*/sequences.gb`) was reused whatever the
  sequence settings, so `--force` on a v1.5.0 output folder aligned its
  nucleotide genomes as protein without a word. It is now reused only when
  type, region, segment, HMM marker settings, `max_per_species` and
  excluded organisms match (recorded beside it in `sequences.key`), and
  downloaded again otherwise. The alignment checkpoint now also covers the
  sequences themselves, not only their (per-run) short ids.
- **Concatenated-marker families ignored `visualization` and `icon`** in
  `global.yaml` (Poxviridae and the other large-DNA-virus presets): their
  pipeline read the settings from the family config, which never has them,
  so tree images always used the built-in line width and icon colours. The
  fallback values now live in one place (`config.GLOBAL_SETTINGS_DEFAULTS`)
  for both pipelines.
- **Re-runs of an older output folder** (`--force`): the per-species
  download check now also covers the content of the HMM profile set (an
  edited or added profile changes which protein is chosen and how it is
  cut), and the folder's taxonomy cache is topped up with the lineages of
  taxa it lacks. Before, new taxa after a re-download got no lineage (grey,
  no genus, no LCA vote).
- **`manual.include` accessions of no NCBI shape** (Pathoplexus ids, which
  the code expects there) were reported as "protein accessions … left out"
  in nucleotide families, as were 6-letter WGS accessions and lower-case
  entries. Only accessions whose shape clearly belongs to the other
  database are left out now, case-insensitively; 6-letter WGS accessions
  are also fetched from nuccore.
- **hmmscan ran on one thread**: `hmm.threads: null` now means the run's
  `--threads`, as documented.
- **Parallel families raced on the shared profile database.** All RdRp
  families use one combined `RdRp` database; with `-j` > 1, two families
  could build or `hmmpress` it at once and scan half-written indexes
  (reproduced: errors in 5 of 5 trials with 12 parallel workers, none with
  the fix). Building and pressing now happen under a file lock.
- **Protein identity in the signal check** dropped `N` (asparagine) as if it
  were an unknown base; identity now counts A/C/G/T/U for nucleotides and the
  20 standard amino acids for proteins, ignoring ambiguity codes and X.
- **Protein residues MAFFT cannot align.** GenBank translations can carry O
  (pyrrolysine, and how some records write an amber read-through codon), U,
  J or an internal stop; MAFFT aborts on them ("Alphabet 'O' is unknown") and
  IQ-TREE would too. Found in the RdRp validation: one tombusvirus RdRp of
  500 wrote its read-through position as O and failed the family. Such
  residues are now written as X (terminal stops and gap symbols dropped) in
  the copy `msa.run_msa` aligns (`seq_type="protein"`), in the single-marker
  and concatenation paths alike; the raw FASTA keeps the deposited residues,
  and the log names the records changed. Quality control already counted
  them as ambiguous, so which records pass is unchanged.

### Survey behind the change (2026-10-20)
Earlier results of 127 families (April and October runs) were re-measured
with the new signal check. Of the 117 families still aligned as whole-genome
nucleotide in v1.5.0, tree_100 was `saturated` for 84, `weak` for 24 and `ok`
for 9 (all with 4–23 tips). Median pairwise identity of their trimmed
alignments was 0.35–0.5 for most, close to what unrelated nucleotide sequence
reaches.

RdRp coverage was then measured on each family's RefSeq proteins (hmmscan
with the pipeline's cutoff code). With Pfam GA and 50 % coverage (the
pipeline default), 78 of 138 Orthornavirae families had ≥ 70 % of their
RefSeq viruses covered; with the rule above, 107. An audit of the hits gained
found them RdRps under other names (potyvirus NIb, alphavirus nsP4, bromovirus
2a, reovirus lambda-3, tombusvirus p86/p88, …) except about 4 in several
thousand. The 31 families below 70 % stay nucleotide; they are listed with
their coverage in `config.py`.

### Validation (2026-10-20/21)
Five families run end to end with fresh configs: tree_100 and tree_500
verdicts of the v1.5.0 whole-genome nucleotide trees (from the survey) against
the v1.6.0 RdRp trees. Medians are subs/site; nucleotide and protein
distances are not on the same scale, which is why the thresholds differ.

| Family | v1.5.0 nucleotide tree_100 / tree_500 | v1.6.0 RdRp tree_100 (median, identity, core) | v1.6.0 RdRp tree_500 |
|---|---|---|---|
| Xinmoviridae | saturated (2.01, core 0.29) / saturated | ok (2.28, 0.40, 0.64) | ok |
| Hantaviridae | saturated (2.07, 3 at cap) / saturated | weak (1.13, 0.73, 1.00; 1 at cap) | weak (1 at cap) |
| Tombusviridae | saturated (2.49, core 0.11) / saturated | ok (2.04, 0.42, 0.72) | ok (2.87) |
| Partitiviridae | saturated (3.05) / saturated | weak (3.66, 0.26, 0.54) | weak (3.37) |
| Togaviridae | saturated (2.21, 3 at cap) / saturated | weak (0.72, 0.78, 0.84; 3 at cap) | weak (3 at cap) |

The v1.6.0 columns are from the rerun after the first round of review
fixes (2026-10-21), before the fixed MAFFT and IQ-TREE seeds below. The
first run's verdicts were the same except Hantaviridae tree_500 (ok, no
branch at the cap) and Tombusviridae tree_500 (weak, 3.04). Both differences
are run-to-run variation: Tombusviridae tree_100 had the same 100 input
sequences in both runs, yet trimAl kept 939 columns of the first alignment
and 441 of the second (median 2.97 → 2.04, core 0.99 → 0.72); see "Repeated
runs gave different alignments and trees" under Fixed. Hantaviridae also
had one L protein cut differently (Buffalo Bayou virus QFQ60707.1, 1297 →
909 aa, now its RdRp domain ± 100 like the rest).

- Hantaviridae's capped branch is Wenling hagfish virus, a RefSeq and the
  most divergent hantavirus. Togaviridae's three lead to unclassified fungal
  and termite "alpha-like" viruses NCBI files under the family; the misfiled
  orbi- and bunyavirus species are all gone (every selected protein fits
  RdRP_2, the alphavirus nsP4 profile, best). Partitiviridae is deep at the
  RdRp level too.
- All 10 PhyloXML files validate against phyloXML 1.20; no IQ-TREE UFBoot
  abort. Wall time per family 11–93 min with 8 threads, of which tree_100
  (ModelFinder `TEST`, UFBoot, up to four rounds of outlier removal) took
  7–52 min.
- Trimming was benchmarked on the eight alignments of Hantaviridae,
  Partitiviridae, Tombusviridae and Xinmoviridae (trimAl `-automated1`,
  `-gappyout`, `-gt 0.5`, none; score: ICTV genera monophyletic among
  classified tips). The modes agreed in 7 of 8, `-gt 0.5` better in the
  eighth; with the RT/Rep and nucleotide sets added, `-gt 0.5` became the
  default (see "Trimming" below).

### RT and Rep validation (2026-10-21/22)
Seven families run end to end with fresh configs. "Without marker": species
with NCBI records but no marker protein, even by translation; "translated":
species found only by translating their genomes.

| Family | Marker | tree_500 / tree_100 tips | Signal 500 / 100 | Without marker | Translated |
|---|---|---|---|---|---|
| Caulimoviridae | RT | 489 / 100 | ok / ok | 9 | 3 |
| Hepadnaviridae | RT | 187 / 98 | ok / ok | 4 | 0 |
| Retroviridae | RT | 491 / 100 | weak / weak | 87 | 54 |
| Circoviridae | Rep | 500 / 100 | ok / ok | 55 | 48 |
| Genomoviridae | Rep | 549 / 100 | ok / ok | 9 | 25 |
| Geminiviridae | Rep | 500 / 100 | ok / ok | 123 | 28 |
| Alphasatellitidae | Rep | 494 / 100 | ok / ok | 4 | 6 |

- Genomoviridae's 549 is what led to the `proportional_merge` fix (Fixed,
  above); it ran before the fix.
- Species in tree_500, old whole-genome tree → new: Retroviridae 144 → 195
  (41 by translation), Circoviridae 393 → 500, Caulimoviridae 155 → 197,
  Hepadnaviridae 52 → 57, Alphasatellitidae 182 → 192, Genomoviridae 471 →
  421. Species of the old trees lost for want of the marker: Retroviridae 5
  (Abelson MLV, avian myelocytomatosis virus, murine osteosarcoma virus,
  LNras*SN, a murine AIDS provirus: no pol), Geminiviridae 10, Circoviridae
  6, Caulimoviridae 1, the others none.
- Species without the marker are records that cannot carry it: coat-protein
  PCR fragments and DNA-B components (most of Geminiviridae's 123), Rep PCR
  fragments shorter than 100 aa, capsid-only records, PCV2-derived agents,
  retroviruses without pol, and one caulimovirus genome (PX904244.1) with no
  RT even among ORFs of ≥ 30 aa.
- All 14 PhyloXML files validate; no ERROR in any log. The quick outlier
  check (Retroviridae, the only family run after it was added): tree_500
  removed 9 sequences over three quick rounds of ~30 s, then one `-n 50`
  search (8 min); tree_100's quick tree was clean, one UFBoot search
  (15 min).
- Marseilleviridae (concatenated, name-based markers) ran end to end too:
  23 genomes, so no outlier rounds; tree_500 `-n 50` in 45 s. 14 species are
  missing for lack of markers, among them Noumeavirus, Golden Marseillevirus
  and three LCMAC genomes (400–900 NCBI proteins each), whose packaging
  ATPase and RNA polymerase subunits NCBI names differently.
- Retroviridae is weak because of trimming: trimAl `-automated1` keeps 58 of
  657 tree_500 columns (Circoviridae 104 of 1660, Genomoviridae 189 of 1217);
  the larger sets have many short PCR fragments. This led to the trimming
  change (Changed, above): `-gt 0.5` keeps 328 Retroviridae columns.

### Trimming (2026-10-22)
The untrimmed alignments of the validation runs (RdRp, RT, Rep) and of the
v1.4.0 run (whole-genome nucleotide), trimmed with trimAl `-automated1` and
with `-gt 0.5`; IQ-TREE `--fast -alrt 1000` (LG+G4 / GTR+G4), 2 threads.
"Genera": ICTV genera with ≥ 2 classified tips that form one clade, of those
testable. Columns: untrimmed → `-automated1` / `-gt 0.5`.

| Protein | Tree | Columns | Genera | Median SH-aLRT |
|---|---|---|---|---|
| Hantaviridae (RdRp) | 500 | 984 → 922 / 923 | 0 / 0 of 5 | 86.6 / 87.1 |
| | 100 | 963 → 922 / 924 | 0 / 0 of 5 | 96.0 / 95.4 |
| Partitiviridae (RdRp) | 500 | 1518 → 255 / 416 | 3 / 3 of 5 | 88.8 / 93.8 |
| | 100 | 876 → 267 / 424 | 4 / 4 of 4 | 91.3 / 92.8 |
| Tombusviridae (RdRp) | 500 | 2361 → 192 / 468 | 9 / **10** of 15 | 87.3 / 93.3 |
| | 100 | 1357 → 939 / 596 | 13 / 13 of 14 | 97.3 / 94.9 |
| Xinmoviridae (RdRp) | 500 | 1336 → 678 / 1000 | 5 / 5 of 8 | 92.0 / 93.4 |
| | 100 | 1329 → 665 / 1012 | 5 / 5 of 8 | 96.7 / 97.9 |
| Alphasatellitidae (Rep) | 500 | 603 → 221 / 315 | 5 / **7** of 13 | 78.4 / 85.2 |
| | 100 | 383 → 225 / 311 | 5 / **6** of 10 | 86.1 / 91.0 |
| Caulimoviridae (RT) | 500 | 442 → 247 / 339 | 3 / 3 of 8 | 81.3 / 85.0 |
| | 100 | 427 → 271 / 360 | 5 / 5 of 6 | 87.2 / 93.2 |
| Circoviridae (Rep) | 500 | 1660 → 104 / 269 | 0 / 0 of 2 | 80.2 / 86.1 |
| | 100 | 804 → 206 / 281 | 2 / 2 of 2 | 82.7 / 85.2 |
| Geminiviridae (Rep) | 500 | 829 → 347 / 361 | 6 / 6 of 8 | 87.4 / 88.2 |
| | 100 | 468 → 347 / 361 | 1 / 1 of 2 | 90.1 / 85.1 |
| Genomoviridae (Rep) | 500 | 1217 → 189 / 323 | 1 / 1 of 7 | 85.8 / 92.0 |
| | 100 | 581 → 235 / 324 | 2 / 2 of 5 | 89.7 / 94.3 |
| Hepadnaviridae (RT) | 500 | 443 → 169 / 325 | 3 / 3 of 4 | 79.8 / 87.2 |
| | 100 | 453 → 170 / 327 | 2 / 2 of 3 | 84.9 / 91.8 |
| Retroviridae (RT) | 500 | 657 → 58 / 328 | 4 / **9** of 9 | 62.3 / 81.0 |
| | 100 | 474 → 304 / 353 | 5 / **7** of 8 | 92.0 / 92.1 |

| Nucleotide | Tree | Columns | Genera | Median SH-aLRT |
|---|---|---|---|---|
| Alphasatellitidae | 500 | 5359 → 983 / 1282 | 7 / 7 of 13 | 93.4 / 96.4 |
| | 100 | 2897 → 492 / 1231 | 11 / 11 of 11 | 89.9 / 98.3 |
| Deltaflexiviridae | 500 | 12134 → 3362 / 7786 | none testable | 99.8 / 100 |
| | 100 | 11817 → 3466 / 7724 | none testable | 98.1 / 100 |
| Hantaviridae | 500 | 9411 → 6506 / 6534 | 1 / 1 of 4 | 100 / 100 |
| | 100 | 8604 → 6370 / 6535 | 1 / 1 of 4 | 100 / 100 |
| Hepaciviridae | 500 | 28181 → 6128 / 8945 | 0 / 0 of 2 | 99.9 / 100 |
| | 100 | 24512 → 2949 / 9027 | 2 / 2 of 2 | 99.7 / 100 |
| Hepadnaviridae | 500 | 8306 → 2157 / 3097 | 3 / 3 of 4 | 97.1 / 99.2 |
| | 100 | 7360 → 2174 / 3113 | 3 / 3 of 4 | 98.9 / 100 |
| Hepeviridae | 500 | 16211 → 5550 / 6961 | 1 / 1 of 6 | 100 / 100 |
| | 100 | 14797 → 3714 / 7092 | 1 / 1 of 6 | 100 / 100 |
| Pneumoviridae | 500 | 19834 → 11264 / 14382 | 2 / 2 of 2 | 100 / 100 |
| | 100 | 19703 → 11634 / 14293 | 2 / 2 of 2 | 100 / 100 |
| Solemoviridae | 500 | 17821 → 4055 / 5298 | 1 / 1 of 4 | 99.1 / 99.7 |
| | 100 | 14987 → 3212 / 5351 | 0 / 0 of 3 | 99.5 / 100 |
| Spinareoviridae | 500 | 7662 → 1951 / 4110 | 5 / **6** of 8 | 92.8 / 98.5 |
| | 100 | 8222 → 5371 / 3997 | 5 / 5 of 8 | 99.7 / 99.9 |
| Xinmoviridae | 500 | 23428 → 3083 / 11636 | 3 / 3 of 4 | 99.5 / 100 |
| | 100 | 26880 → 3581 / 11831 | 3 / 3 of 4 | 99.8 / 100 |

- `-gt 0.5` was never worse on genera and better in 6 of 42. Median
  SH-aLRT was higher in 19 of 22 protein trees (lower: Geminiviridae
  tree_100 by 5, Tombusviridae tree_100 by 2.4, Hantaviridae tree_100 by
  0.6) and higher or equal in every nucleotide tree, where it sits near 100
  in either mode, so the nucleotide half is a weak test.
- Tree time over all alignments: protein ×1.22 (RdRp ×1.17, RT/Rep ×1.26),
  nucleotide ×1.40; the most, Retroviridae tree_500, ×2.6 (21 → 55 s).
  `-gt 0.5` is not always the larger: Tombusviridae tree_100 and
  Spinareoviridae nucleotide tree_100 keep fewer columns than with
  `-automated1`, whose strict/gappyout choice flips between alignments.
- Signal verdicts on these trees, `-automated1` → `-gt 0.5`: protein
  Partitiviridae tree_500 and tree_100 and Retroviridae tree_100 weak →
  saturated, Tombusviridae tree_100 weak → ok and tree_500 ok → weak;
  nucleotide Hepadnaviridae tree_500 and tree_100, Alphasatellitidae and
  Solemoviridae tree_100 weak → saturated; the other 32 unchanged. The
  kept columns are more variable, so the median distance rises; the
  core-fraction reason (Retroviridae tree_500: 18 % of a typical sequence
  kept) disappears. These are `--fast` trees, so the pipeline's own
  verdicts can differ. Verdicts from runs with different trimming are not
  comparable.
- trimAl drops sequences that are all gaps after trimming. In both modes
  that happened once in the 42 alignments (Tombusviridae tree_500, the
  same sequence).

### tree_500 search effort (2026-10-21)
The five-family rerun's trimmed tree_500 alignments, Tombusviridae (498 ×
236 aa) and Hantaviridae (345 × 922 aa), LG+G with `-alrt 1000`, 4 threads
while other runs used the rest of the machine. Reference: the best full
search. "Contradicted ≥ 80 / ≥ 95": splits with that SH-aLRT support which
the reference tree contradicts (incompatible splits). RF: share of splits
that differ from the reference.

| Tombusviridae | log-likelihood | CPU time | RF | contradicted ≥ 80 / ≥ 95 |
|---|---|---|---|---|
| `--fast`, seeds 1–5 | −104,656 … −104,863 | 1.5 min | 26–32 % | 49–57 / 6–13 |
| `--fast --runs 5` | −104,656 | 10 min | 26 % | 57 / 11 |
| `-n 20` | −104,584 | 4.4 min | 21 % | 34 / 5 |
| `-n 50` | −104,483 | 8.3 min | 16 % | 22 / 4 |
| `-n 100` | −104,474 | 14 min | 12 % | 18 / 3 |
| full search, seed 2 | −104,451 | 66 min | 13 % | 22 / 0 |
| full search, seed 1 (reference) | −104,429 | 53 min | — | — |

| Hantaviridae | log-likelihood | CPU time | RF | contradicted ≥ 80 / ≥ 95 |
|---|---|---|---|---|
| `--fast`, seeds 1–5 | −53,824 … −53,848 | 1.4–1.8 min | 12–18 % | 6–13 / 0–4 |
| `-n 20` | −53,805 | 8.5 min | 7 % | 2 / 0 |
| `-n 50` | −53,805 | 19 min | 8 % | 2 / 0 |
| `-n 100` | −53,801 | 36 min | 5 % | 1 / 0 |
| full search, seed 2 | −53,801 | 88 min | 8.5 % | 2 / 0 |
| full search, seed 1 (reference) | −53,800 | 125 min | — | — |

Genus monophyly was the same in every tree (Tombusviridae 9 of 15 genera;
Hantaviridae 0 of 5, because of NCBI labels), so `--fast` does not change
the genus-level picture; the differences are deeper splits. The full
search's first ~10 iterations already beat every `--fast` run, which is why
a capped search gets most of the gain. Per family, the eight validation
families needed 1–3 tree_500 builds (8–75 s each with `--fast`) and 1–4
tree_100 builds (9–26 min each). `-n 50` alone would have added ~15 min per
protein family; with the quick outlier check (Changed, above) Togaviridae
went from 92 to 52 min, `-n 50` included.

## [1.5.0] — 2026-10-20

Inference defaults from the 2026-10-20 accuracy review (benchmarks below).
Only templates changed (`DEFAULT_FAMILY_CONFIG`, `GLOBAL_CONFIG_TEMPLATE`);
existing YAML files keep their values, so edit `config/global.yaml` and use a
fresh per-family configs directory to pick these up. A global.yaml that still
says `tree_500: tool: fasttree` now fails at load time (see *Added*).

### Changed
- **tree_500 is inferred with IQ-TREE `--fast`** instead of FastTree: same wall
  time with a clearly better likelihood (benchmarks below). Support becomes
  SH-aLRT (`tree500_support_type` = `SH_aLRT`; was `SH_like`). FastTree remains
  available: `tool: fasttree` together with `options_nuc: ""` and
  `options_aa: ""`. FastTree is optional in `vfam_trees test` now.
- **tree_100 runs a full IQ-TREE search with `-alrt 1000 -B 1000 --bnni`** for
  both nucleotide and protein trees (was `--fast` for nucleotide, `-B 1000` for
  protein). The outputs carry the UFBoot value (`tree100_support_type` =
  `UFBoot`); the `.treefile` under `_work/100/tree/` keeps all values as
  `alrt/ufboot` labels.
- **Nucleotide MAFFT runs use `--adjustdirection`** (msa_500 and msa_100).
  GenBank records deposited on the opposite strand are reverse-complemented
  before alignment instead of entering the tree as nonsense long branches
  (2 of 202 hantavirus L segments and 1 of 500 alphasatellites in the
  2026-10-20 run). MAFFT orients against the first record; when that record
  is on the minority strand, the finished alignment is reverse-complemented
  as a whole so it sits on the majority strand (ties go to the strand of the
  smallest id), so the result never depends on record order and MAFFT runs
  only once. The flipped accessions are logged at WARNING level, also when a
  cached alignment is reused (recorded in `<alignment>.flipped`). The
  `*_sequences_raw_*.fasta` files keep the deposited orientation.
- msa_500 nucleotide options are `--retree 2 --adjustdirection` (the redundant
  `--6merpair` is gone; it is MAFFT's default distance anyway).
- **One support value per node, for any combination of IQ-TREE measures.**
  `-alrt`/`--alrt`, `-lbp`/`--lbp`, `-abayes`/`--abayes`, `-B`/`--ufboot`/`-bb`
  and `-b`/`--boot` are all recognised; IQ-TREE's slash-joined labels
  (`alrt/lbp/abayes/ufboot`) are reduced to the bootstrap-type value, else
  SH-aLRT, and `support_type_for` names that value. The `SH_aLRT_UFBoot`
  label is gone.
- The tree cache key includes the tree tool's version, so trees are rebuilt
  after a FastTree / IQ-TREE upgrade.

### Added
- **Tree options are checked at config load.** `tree_500` / `tree_100` must
  name `iqtree` or `fasttree`, and the options key the run uses must not hold
  the other tool's options (e.g. `--fast` under FastTree, `-gtr` under
  IQ-TREE). Previously such a config failed per family after fetch and
  alignment.
- `<prefix>.tree_info.json` beside every tree in `_work/<label>/tree/`: tool,
  options, the support measure the tree carries, and a note when IQ-TREE's
  UFBoot fallback was used. Written for every tool and removed before every
  run, so it always describes the tree next to it.

### Fixed
- **A UFBoot run that aborts no longer fails the family.** IQ-TREE 2.3.2 can
  abort in its post-bootstrap consensus step (`Assertion '!findSplit(sp)'
  failed` in `hashsplitset.cpp`) on particular taxon sets; the tree file it
  leaves holds the final ML tree but no support values. The wrapper now
  scores that fixed tree with SH-aLRT (`-t … --tree-fix`; ModelFinder's pick
  is read from the run's log, a partitioned run's from its
  `.best_scheme.nex`), which takes seconds, and repeats the search without
  UFBoot only if that is impossible. `summary.tsv`, the PDF and PhyloXML then
  report `SH_aLRT` and the options it came from (`-alrt 1000`), not UFBoot.
- **Branch-length outliers with IQ-TREE trees.** Zero-length tips are left out
  of the median / MAD, but the filter tested `> 0` and IQ-TREE writes them as
  1e-6 (FastTree writes 0.0). On a tree_500 with many identical tips (e.g.
  outbreak genomes) the threshold collapsed and ordinary tips were removed;
  in a synthetic tree of 30 identical, 30 near-identical and 40 ordinary tips
  all 40 ordinary tips were flagged. Branches ≤ 1e-5 now count as zero.
- **The tree model now reaches `summary.tsv`.** `tree{500,100}_tree_model` and
  `tree{500,100}_seq_type` were written by `build_summary_row` but missing from
  `summary.COLUMNS`, and `csv.DictWriter(extrasaction="ignore")` dropped them
  without a word (the per-family PDF was unaffected — it reads the in-memory
  row). All four are declared now; a test asserts every key the builder
  writes is declared.
- **Appending to a `summary.tsv` / status TSV from another version** no longer
  shifts columns: a file whose header differs is rewritten under the current
  columns (plus any extra old ones) before the row is added, rows matched by
  column name.
- `vfam_trees test` no longer fails on a machine without FastTree, and the
  FastTree version probe can no longer wait on standard input.

### Benchmarks (2026-10-20)
Hantavirus L segment from the 2026-10-20 run; 3 threads under load, so times
are conservative. IQ-TREE 2.3.2, FastTree 2.1.11, MAFFT 7.520, trimAl 1.5.
Likelihoods are for each topology re-optimised under GTR+G on the same
alignment.

| tree_500 set (202 × 6506 nt) | Time | lnL |
|---|---|---|
| FastTree `-nt -gtr -gamma` (old default) | 15 s | −418,410 |
| FastTree `-spr 4 -mlacc 2 -slownni` | 30 s | −418,385 |
| IQ-TREE `--fast -alrt 1000` (new default) | 17 s | −418,173 |
| IQ-TREE full search | 451 s | −418,168 |

| tree_100 set (99 × 6370 nt) | Time | Note |
|---|---|---|
| IQ-TREE `--fast` (old default) | 7 s | same topology as the full search |
| full search `-alrt 1000 -B 1000 --bnni` (new default) | 248 s | adds UFBoot |
| full search `-m MFP` | 1066 s | GTR+F+I+R10, lnL +680 |

MAFFT modes on the tree_100 raw set (`--retree 2` 63 s, FFT-NS-i 88 s,
E-INS-i 693 s, L-INS-i 1271 s) kept 74–80% of columns after trimAl at the same
median pairwise identity, so alignment mode is not where accuracy is lost;
whole-genome nucleotide alignment of divergent families is (e.g. Xinmoviridae
keeps 13% of columns at median identity 0.50).

## [1.4.0] — 2026-10-20

### Added
- **Archaeopteryx.js page per tree.** `<prefix>_tree_{500,100}.html` inlines the
  PhyloXML in a page that opens it in Archaeopteryx.js 3.28.0 (viewer, d3, sax,
  phyloxml and the PDF-export libraries loaded from jsDelivr at pinned versions;
  needs an internet connection). It opens coloured by Genus (by Species on a
  one-genus tree), as the PDF figures are, when the viewer can show that field in
  full (2–20 values); otherwise the viewer chooses. The URL hash carries the view,
  so a copied link reopens the tree as left. Toggle with `output.archaeopteryx_html`
  (default on).
- **Tip dates in PhyloXML.** Each tip with a parsable collection date gets
  `<date unit="year">` with the decimal year (midpoint of the stated day / month /
  year / range, leap-year aware — the convention both Archaeopteryx viewers use)
  and the raw GenBank string in `<desc>`. This is where both viewers read tip
  dates, e.g. for the root-to-tip clock plot (offered when every tip is dated).
- **`branch_length_unit="subs/site"`** on `<phylogeny>`: both viewers print it on
  the scale bar and use it as the clock-plot rate unit, and it keeps a tree with
  tip dates from being drawn on a calendar axis in Archaeopteryx.js.
- **`vipr:Country`** (PhyloXML) / **`country`** (Auspice): the country part of the
  location, for grouping and colouring.
- NCBI taxids on internal-node `<taxonomy>` (`<id provider="ncbi">`).

### Changed
- **`vipr:Species` is the ICTV species** (the species-rank entry of the NCBI
  lineage, e.g. "Orthoebolavirus zairense"), consistent with the internal LCA
  labels. It used to be the GenBank organism string ("Zaire ebolavirus", "Bundibugyo
  virus", "Sapovirus sp."), which both Archaeopteryx viewers reported as
  contradicting the internal nodes, and which split one species into several
  clades in "Clades by rank". A tip whose lineage has no species rank gets no
  `vipr:Species`; without a ranked lineage (pasted sequences) the organism is kept.
  The Auspice `species` trait follows. Leaf labels are unchanged.
- **Leaf `<taxonomy>` is the tip's own NCBI taxon**: its id, NCBI's current name
  for that id (GenBank organism strings can be stale), and `<rank>` when NCBI's
  rank is a phyloXML rank.
- **`vipr:Location`** is written only when it says more than the country
  ("USA: Texas, Austin"); a country-only location is `vipr:Country` alone. The
  Auspice `location` trait follows.
- **Accession sources.** Pathoplexus accessions are `source="pathoplexus"` with a
  `<uri>` to the record page; pasted / FASTA-file sequences are `source="user"`.
  Both were `source="ncbi"`, which gave dead "open database" links in the viewers.
- **`style:font_color` is off by default** (`output.phyloxml_font_colors`): both
  viewers open a tree with their own Color-by, which hides style colours, and these
  have no legend. PDF/PNG and Auspice keep the genus palette.

### Fixed
- **Location was empty for every NCBI record.** INSDC renamed the GenBank
  `/country` qualifier to `/geo_loc_name`; it is now read (with `/country` as the
  fallback), filling `vipr:Location`/`vipr:Country`, the Auspice traits, the
  metadata TSVs and the `{location}` label placeholder.
- **"unclassified X" internal nodes no longer get a Linnaean rank.** A rank was
  inferred from the name suffix ("unclassified Caliciviridae" → family) even when
  NCBI says "no rank"; NCBI's rank is now authoritative and the suffix is used
  only for names NCBI did not supply. Such nodes also no longer get a family
  label in the PDF tree images.
- **Placeholder values are no longer written as data.** "n/a", "none" and the
  INSDC missing-value terms ("not collected", "missing: control sample", …) are
  treated like "unknown" in PhyloXML properties and Auspice traits alike; before,
  PhyloXML skipped only "unknown", so the two files disagreed.
- **False "absent from the tree" warning for `additional_data_sources`.** Every
  injected record (e.g. Pathoplexus) was reported as a `manual.include` accession
  "not found … (will be absent from the tree)", though it was in the tree.
- A pair of sister tips with the same ICTV species no longer gets a redundant
  species-level internal label in the PhyloXML when their GenBank organism strings
  differ.

## [1.3.1] — 2026-06-11

### Changed
- **Poxviridae HMM capping-enzyme marker fixed.** A live 761-species Poxviridae
  HMM run scored **0/499** coverage on the mRNA capping enzyme (D1R): its bundled
  profile `mRNA_cap_enzyme` (PF01331, the generic catalytic domain) — and the
  alternative PF03919 — do not match the poxvirus D1 large subunit at all. The
  marker now uses **`mRNA_G-N7_MeTrfase`** (PF03291), the cap N7-methyltransferase
  domain carried by D1, which matches strongly (score ~58 vs GA ~27 on vaccinia).
  Paired with a lowered Poxviridae `hmm.relative_length_cutoff` of **0.40** (from
  0.50) so the divergent *Avipoxvirus* D1 — a strong hit (~53 vs GA ~27) that
  covers only 47 % of the model — is also admitted. Both Orthopox and Avipox D1
  now recover. The other 6 Poxviridae markers were already at 499/499 and are
  unaffected. Re-run `scripts/build_bundled_hmms.sh` to refresh the bundled set
  (the `mRNA_cap_enzyme.hmm` profile is replaced by `mRNA_G-N7_MeTrfase.hmm`).

## [1.3.0] — 2026-06-08

### Added
- **HMM-based marker-protein identification.** For large DNA virus families
  whose GenBank annotation is inconsistent, marker proteins can now be
  identified by **HMMER profile homology** instead of protein-name queries —
  far more robust against synonyms, gene-symbol vs descriptive names, and
  mis-/un-annotated or divergent orthologs. Enabled per family via a new
  top-level `hmm:` config block (`enabled`, `database`, `default_evalue`,
  `use_ga_when_available`, `relative_length_cutoff`,
  `multidomain_overlap_tolerance`, `threads`).
  - Markers carry an `hmms:` token list — a single profile (`"DNA_pol_B"`) or a
    **multidomain architecture** (`"Adeno_hexon--Adeno_hexon_C"`, HMMs in N-to-C
    order; tokens within a marker are alternatives, OR'd). The HMM tier is
    authoritative for markers that declare `hmms`; markers without fall back to
    name+alias matching, so mixed sets work.
  - Works on **both** the single-marker (`DNA_FAMILIES`) and concatenated-marker
    (`CONCATENATION_FAMILIES`) paths via the `markers.MarkerIdentifier` swap. In
    HMM mode the per-marker name queries are replaced by an all-proteins-per-
    species fetch; proteins are grouped by source genome and assigned to markers
    by HMM, with hits gated on the curated Pfam **GA** bit-score (or an E-value)
    plus an HMM-model coverage cutoff. `max_per_species` then bounds the number
    of genomes kept (RefSeq-first).
  - Ships curated **Pfam-A** profile sets (CC0) under
    `vfam_trees/data/hmms/<Family>/`: **Poxviridae** (8 markers) and
    **Adenoviridae** (hexon, as the 2-domain `Adeno_hexon--Adeno_hexon_C`
    architecture OR'd with a single-domain `Adeno_hexon` fallback), both
    enabled by default in their presets. Rebuildable via `scripts/build_bundled_hmms.sh`.
  - Requires **HMMER** (`hmmscan`/`hmmpress`) on PATH when `hmm.enabled`;
    `vfam_trees test` now reports their availability. The database is
    auto-`hmmpress`-ed on first use. The global sequence cache is bypassed in
    HMM mode (its key does not encode the HMM database).

## [1.2.51] — 2026-06-08

### Changed
- **Docs: Auspice JSON in `OUTPUT_FILES.md`.** The detailed per-family output-files
  reference never described the `*_auspice.json` output; it now has its own entry
  under §1 (divergence-only Auspice v2 dataset, viewable at auspice.us / Nextclade /
  `auspice view`, colorings/filters, and the **Outbreak** trait from
  `additional_data_sources` tips). Completes the v1.2.50 documentation pass.
  Documentation only; no behaviour change.

## [1.2.50] — 2026-06-08

### Changed
- **Docs: Auspice JSON viewer compatibility.** README and `run --help` now state
  that the exported `*_auspice.json` is a schema-valid Auspice v2 dataset and name
  the viewers it loads in — [auspice.us](https://auspice.us), Nextclade, and a
  local `auspice view` — and note the **Outbreak** colouring/filter contributed by
  `additional_data_sources` tips. Documentation only; no behaviour change.

## [1.2.49] — 2026-06-06

### Fixed
- **Auspice JSON now displays in Auspice / auspice.us / Nextclade.** The exported
  `*_auspice.json` violated the Auspice v2 dataset schema in two ways, and strict
  viewers reject the *entire* dataset on any violation (rendering nothing):
  - `node_attrs.accession` is a *reserved* schema field that must be a plain
    string; it was emitted as a `{"value": …}` object, which failed validation
    for the whole tree. It is now a plain string, guarded by the schema's
    accession pattern (a malformed value is dropped rather than sinking the file).
  - `display_defaults.layout` was `"rectangular"`, which is not an allowed value;
    it is now `"rect"` (one of `rect`/`radial`/`unrooted`/`clock`).

  Both the single-protein and concatenated-marker pipelines share the one writer,
  so every `*_auspice.json` is fixed. Re-run a family to regenerate valid output.
  A new regression test validates the export against the vendored official
  Auspice v2 schema (offline) to guard this bug class going forward.

## [1.2.48] — 2026-06-04

### Added
- **`additional_data_sources` gains two optional presentation fields.**
  - `name_prefix:` — a string prepended *verbatim* to every leaf label produced
    from that source (e.g. `name_prefix: PATHOPLEXUS_`), so externally-sourced
    tips stand out at a glance. The prefix flows into the leaf label, the Newick,
    the `*_id_map.tsv` display column, the PhyloXML `<name>`, and the Auspice
    `display_name`. You supply your own trailing separator.
  - `outbreak_name:` — a string (e.g. `outbreak_name: Bdbv-2026`) that tags each
    injected tip. It is emitted in PhyloXML as a `vipr:outbreak` property and, in
    the Auspice JSON, as an **Outbreak** colouring + filter — so the outbreak
    clade can be isolated or coloured in Nextstrain. Both are surfaced only when
    at least one tip carries the value, so ordinary trees are unchanged.

  Both default to empty (off) and apply only to `additional_data_sources` records.

## [1.2.47] — 2026-06-04

### Changed
- **Pathoplexus injection: outbreak tips now get a species label and genus
  colour even when the record has no NCBI taxon.**  Fresh outbreak genomes
  submitted directly to Pathoplexus often lack `ncbiVirusName` / `ncbiVirusTaxId`
  (NCBI hasn't catalogued them yet) — exactly the records this feature targets —
  so they previously rendered gray and unlabelled.  Each of the 14 Pathoplexus
  organisms now maps to its virus taxon (`PATHOPLEXUS_ORGANISM_TAXON` in
  `datasources.py`, e.g. `ebola-bdbv → Bundibugyo ebolavirus / taxid 565995`);
  when a record lacks an NCBI taxid, the organism's name and taxid are used so
  the tip is labelled and resolves to a genus colour.  Records that *do* carry
  NCBI metadata are unaffected.

## [1.2.46] — 2026-06-04

### Added
- **Pathoplexus as an additional sequence source (`additional_data_sources:`).**
  A new per-family config block pulls *forced-include* sequences from
  [Pathoplexus](https://pathoplexus.org) (via its GenSpectrum LAPIS API) and
  injects them into the tree like `manual.include` records: they bypass QC and
  clustering/subsampling and are protected from length-/branch-outlier removal.
  The typical use is adding the latest **outbreak** genomes to an existing
  family tree.  Unlike pasted sequences, these carry real host / location /
  collection-date / taxon metadata, so they are labelled and coloured by genus.

  Each entry names a `source` (currently only `pathoplexus`) and an `organism`
  slug (Pathoplexus is keyed by curated organism, not taxid, so a genus tree
  lists several entries — e.g. `ebola-zaire`, `ebola-sudan`, `ebola-bdbv`), with
  optional `country`, `host`, `date_from`, `date_to`, `max_seqs` (default 200),
  and `dedup_vs_ncbi` (default true) limiters.  Only the latest version of each
  accession is fetched, and records already present in the NCBI download (by
  INSDC accession) are dropped.  **Nucleotide-only**; not supported in concat
  mode.  Single source of truth: `vfam_trees/datasources.py`.

## [1.2.45] — 2026-06-01

### Fixed
- **`vfam_trees --version` now reports the actual package version.**  The
  version option read installed distribution metadata, so it showed whatever
  version was current at the last `pip install -e .` (e.g. `1.2.42`) rather
  than the running code.  It is now wired directly to
  `vfam_trees.__version__`.

### Changed
- Documentation: README now covers the Auspice v2 JSON export and the
  `manual.name`-driven output prefix / directory naming; the `run` command
  help text lists Auspice JSON among the per-family outputs.

## [1.2.44] — 2026-06-01

### Changed
- **`manual.name` now drives the output prefix and directory name.**  When set,
  every output file and the output directory are named from a filesystem-safe
  form of `manual.name` instead of the biological family name — e.g.
  `results/Ebola_3052317/Ebola_tree_100.nwk` instead of
  `results/Orthoebolavirus_3052317/Orthoebolavirus_tree_100.nwk`.  The taxid
  suffix on the directory is retained.  Families **without** `manual.name` are
  completely unaffected (the sanitized prefix equals the family name, so output
  paths are byte-identical to before).  NCBI queries, taxonomy resolution,
  sequence caching, the `.done_<family>` sentinel, and the `summary.tsv`
  `family` column continue to use the biological family name.

### Added
- Each output directory now carries a `.family` marker recording the biological
  family that owns it.  This lets `status` (and other CLI subcommands) locate a
  `manual.name`-renamed directory without loading config, and guards against two
  families resolving to the same directory: such a collision (only possible when
  names match *and* the taxid is unresolved) aborts the run with a clear error.

## [1.2.43] — 2026-06-01

### Added
- **Auspice v2 JSON export (Nextstrain interactive trees).**  Each family now
  emits `{Family}_tree_{500,100}_auspice.json` alongside the PhyloXML, loadable
  directly in [auspice.us](https://auspice.us) or a local `auspice view`.
  Divergence trees only — branch lengths are substitutions per site, so the
  tree is not time-resolved (no temporal layout).  The "color by genus" view
  reuses the same palette as the PDF/PNG/PhyloXML output (emitted as an explicit
  coloring `scale`), and genus/subfamily/species/host/location/year are
  available as colorings and filters.  Branch support is shown as a toggleable
  branch label; internal-node LCA names appear as `clade` branch labels.
  Controlled by the new `output.auspice_json` family-config flag (default
  `true`); supported in both the single-protein and concatenated-marker paths.

### Added
- **Numeric-field validation on family config load.**  Catches typos like
  `max_100: 4  00` (which YAML parses as the string `"4  00"`) at config-load
  time instead of letting them propagate into a downstream `%d` log format or
  arithmetic op many minutes into a run.  Validated fields cover `download`,
  `quality`, `clustering`, `targets`, `length_outlier`, `outlier_removal`, and
  `refseq_absorption` blocks; rejects non-numeric, out-of-range, and (where
  applicable) zero/negative values.  Error messages include the family name,
  the dotted field path, the offending value, and a hint when the value looks
  like a YAML typo.

## [1.2.41] — 2026-05-15

### Fixed
- **PhyloXML write crashed with `UnboundLocalError: display_name`** when
  `manual.name` was set.  The `display_name` value is computed in
  `run_family` (and `run_family_concat`) but was not threaded through to
  `_run_target` / `_run_target_concat`, where it is consumed when building
  the PhyloXML `<name>` element.  In `_run_target` the situation was made
  worse by local rebindings of `display_name` inside the branch-outlier
  logging loops, which made the parameter look local-from-entry to Python
  and raised `UnboundLocalError` at the first read.  Now `display_name` is
  passed explicitly to both helpers and the loop-local variable is renamed
  to `leaf_display` so it can no longer shadow the family-level value.

## [1.2.40] — 2026-05-15

### Changed
- **`manual.limit_lineages` renamed to `manual.restrict_to_lineages`** — the
  previous name suggested the wrong polarity (it sounded like the listed
  lineages would be *limited / removed*, when in fact the analysis is
  *restricted to* those lineages and everything else is dropped).  Behaviour
  is unchanged; only the key name moves.  The old `limit_lineages` key is now
  a **hard error** — configs must be renamed.  Pre-production, no migration
  path.

## [1.2.39] — 2026-05-15

### Changed
- **`manual.include_species` renamed to `manual.limit_lineages`** and
  generalized to match at any taxonomic rank.  Each entry is a scientific
  name or NCBI taxid for a taxon at any rank (species, genus, subfamily,
  ...); the pipeline queries NCBI Taxonomy for all species-rank descendants
  under each listed taxon and restricts the discovered species list to
  those falling under at least one entry.  The old `include_species` key
  is now a **hard error** — configs must be renamed.  Pre-production, no
  migration path.

## [1.2.38] — 2026-05-14

### Added
- **`manual.name`** — new per-family config key (default `""`) that overrides
  the display name used in all PDF/PNG titles, the PhyloXML `<name>` element,
  and the overview grid thumbnail.  When empty, the biological family name is
  used as before.  Output file names are unaffected.

## [1.2.37] — 2026-05-13

### Fixed
- **`manual.include_fasta_files` no longer duplicates single-token headers** —
  a FASTA header with no whitespace (e.g. a pre-formatted pipe-delimited
  label) was used as both the sequence id *and* the organism, so the
  leaf-label formatter rendered it twice (`HEADER|HEADER`).  Single-token
  headers now yield an empty organism, so the leaf label is just the header
  itself.  Such entries are bucketed by id rather than sharing one empty
  "species" bucket.

## [1.2.36] — 2026-05-13

### Changed
- **`manual.include_fasta` renamed to `manual.include_seq`** — the per-family
  YAML key for pasting inline sequences is now `include_seq`.  Configs that
  still use `include_fasta` emit a WARNING at load time and are otherwise
  unchanged (the stale key is ignored; `include_seq` defaults to empty).

### Added
- **`manual.include_fasta_files`** — new per-family YAML key that accepts a
  list of paths to FASTA-formatted files.  Sequences are injected identically
  to `include_seq` entries: the FASTA id field (first whitespace-delimited
  token of the header) becomes the sequence id; the remainder of the header
  line becomes the organism/name (a single-token header yields an empty
  organism, so the leaf label is just the header itself).
  Subject to the same collision rules, QC bypass, downstream protection, and
  concat-mode restriction as `include_seq`.  Paths are resolved relative to
  the working directory where `vfam_trees run` is invoked.

## [1.2.34] — 2026-05-12

### Changed
- **Format-string leaf labels** — the `labeling:` config block now accepts a
  `format:` string with named placeholders `{species}`, `{id}`, `{host}`,
  `{strain}`, `{location}`, `{year}`, `{genus}`.  The user controls which
  fields appear, in what order, and with what separators.  Default:
  `{species}|{id}|{host}` (preserves previous behaviour).
  Two accompanying boolean options:
  - `replace_whitespace` (default `true`) — spaces in field values become
    underscores so Newick labels remain parseable.
  - `keep_separator_on_empty` (default `false`) — when a field value is
    absent, `null`, `"unknown"`, or `"n/a"` (case-insensitive), its adjacent
    separator literal is also dropped so no leading or consecutive separators
    appear.  Set `true` to preserve separators regardless.
  Applies to PhyloXML `<name>` elements and all graphical output labels (PDF,
  PNG) in both single-protein and concat pipeline modes.

### Removed
- **`labeling.show_strain`** (added in v1.2.32) — superseded by the new
  format-string approach.  Use `format: "{species}|{strain}|{id}|{host}"`
  (or any ordering) instead.

## [1.2.33] — 2026-05-12

### Changed
- **`tree_500` config now mirrors `tree_100`** — the single `options:` key
  under `tree_500:` has been replaced with `options_nuc:` and `options_aa:`
  so nucleotide and protein runs can be tuned independently, exactly as
  `tree_100:` already allowed.  Existing per-family YAMLs that still carry
  the bare `options:` key continue to work via a backwards-compatible
  fallback in `_resolve_tree_options`.  The concat pipeline's tree_500 path
  now reads `options_aa` (concat is always protein) with the same fallback.

## [1.2.32] — 2026-05-12

### Added
- **`manual.include_species`** — new optional key in the per-family `manual:`
  config block. Accepts a list of species names (matched case-insensitively
  against NCBI species names) and/or numeric NCBI taxids (integer or digit
  string). When set, only the listed species proceed to download; all others
  are skipped before any network calls are made. Unmatched entries produce a
  WARNING; if the list matches nothing the family is skipped with an
  explanatory status entry. Duplicates (e.g. the same taxid supplied as both
  an integer and a digit string) are silently deduped. Empty or absent key
  preserves the existing behaviour (full discovered species list used).

### Fixed
- **`manual.include` sequences missing from `tree_100`** — sequences listed in
  `manual.include` were correctly included in `tree_500` but could be silently
  dropped from `tree_100` when clustering at a tighter threshold elected a
  different cluster representative. `_cluster_at` now force-appends any
  protected sequence absent from the clustering tool's output, so manual
  includes (and RefSeqs) are guaranteed present in both trees.
- **Spurious "Auto-configured" log messages** — `_apply_smart_defaults`
  previously logged "Auto-configured segment / concatenation / DNA marker for
  …" even when the family config file already specified those values. The
  messages are now suppressed whenever the file config contains the
  overriding key.

### Changed
- **Unclassified taxa coloring** — taxa whose NCBI `species` field starts with
  `"unclassified "` (case-insensitive, requires trailing space) are now always
  rendered as medium grey (`#999999`) in PDF/PNG tree images regardless of any
  inferred genus. This prevents spurious hue-colored leaves for sequences
  whose taxonomic classification is not yet resolved.

## [1.2.31] — 2026-05-11

### Added
- **Per-family `manual.include_fasta`** (renamed to `manual.include_seq` in
  a later release) — new optional block in per-family YAML configs that lets
  curators inject sequences not yet in GenBank directly into the pipeline.
  Each entry is a mapping with required keys `id`, `organism`, and `sequence`
  (whitespace is stripped and the sequence is uppercased).
  - Pasted entries are injected after the per-species fetch loop, fully
    bypass QC (length, ambiguity, organism exclusion), and their ids are
    added to the same protected set as `manual.include` accessions —
    they survive clustering, proportional merge, and length-outlier
    filtering.
  - Bucketed by organism: an entry whose `organism` matches a fetched
    species joins that species' bucket; otherwise a new species bucket
    is created.
  - Validation rejects non-mapping entries, missing / empty
    `id` / `organism` / `sequence`, duplicate ids within `include_seq`,
    and id collisions with `manual.include` / `manual.exclude`. At
    pipeline time, id collisions with any accession returned by NCBI
    for the family raise a hard error so the run aborts before MSA.
  - Pasted leaves render gray (no genus) and do not participate in LCA
    voting (no taxid / no species-rank lineage) — by design, since only
    `id` and `organism` are required.
  - Not supported when `sequence.region == "concatenated"` — a non-empty
    `include_seq` in concat mode is rejected at config-load time.
  - Existing per-family YAMLs without `include_seq` continue to load
    unchanged (defaults to an empty list).

## [1.2.28] — 2026-05-08

### Fixed
- **Final source-nuc length filter cleanup.** Two complementary changes:
  - `_source_nuc_accession`'s `db_source` fallback now only accepts known
    RefSeq nucleotide prefixes (`NC_`, `NG_`, `NM_`, `NR_`, `NS_`, `NT_`,
    `NW_`, `NZ_`, `AC_`, `XM_`, `XR_`).  In real data, `db_source` for a
    protein record describes the protein record's own home — RefSeq protein
    (`YP_`/`NP_`), GenBank protein (`ABO61246`), or UniProt (`P12345`) —
    never the source nucleotide.  An unrestricted shape match was leaking
    these as if they were source-nucs.  UniProt 1-letter+5-digit accessions
    in particular pass the shape filter but resolve to the protein db at
    NCBI, which is the residual source of the `Otherdb db=protein` errors
    in v1.2.27.
  - `fetch_nuc_lengths` now recognises deterministic NCBI errors (`Otherdb`,
    `Invalid uid`) and recovers via binary split instead of hammering the
    same poison batch five times.  A single bad UID in a 200-accession batch
    now costs ~16 esummary calls (≈2·log₂N) and still returns lengths for
    the 199 good UIDs, instead of erasing the whole batch's length map.

## [1.2.27] — 2026-05-08

### Fixed
- **Remaining "Otherdb db=protein" warnings from the v1.2.24 source-nuc length
  filter.** v1.2.26 caught RefSeq protein accessions (`YP_…`, `NP_…`, …) but
  missed the more common case: bare 3-letter GenBank protein accessions like
  `ABO61246.1` and `AAC54321.1`. These have nuccore-shape (letters + digits)
  but always resolve in the protein database, and are exactly the form that
  the partial single-gene submissions surface in `db_source` for non-RefSeq
  records. The shape regex now requires either an underscore separator
  (any 1-4 letter prefix) or a no-underscore prefix that's specifically 1-2
  letters (GenBank nucleotide) or 4 letters (WGS nucleotide) — 3-letter
  no-underscore (GenBank protein) is rejected. Combined with the v1.2.26
  RefSeq protein-prefix deny list, this should make Asfarviridae and similar
  fragmentation-prone families clean of nuc-length esummary warnings.

## [1.2.26] — 2026-05-08

### Fixed
- **Spurious "Otherdb db=protein" warnings from the v1.2.24 source-nuc length
  filter.** The `_source_nuc_accession` `db_source` fallback was extracting
  the protein record's *own* RefSeq accession (`YP_…`, `NP_…`, etc.) — RefSeq
  protein records have a `db_source` like `"REFSEQ: accession YP_009047263.1"`
  describing the protein itself, not its source nucleotide. Sending a protein
  accession to nuccore esummary makes NCBI return an `Otherdb uid=… db=protein`
  error that aborts the entire batch (echoing the resolved numeric UID back
  in `term=`, which is how the warning ended up looking like a bare GI). Both
  `_source_nuc_accession` and `fetch_nuc_lengths` now reject the seven known
  RefSeq protein prefixes (`NP_`, `XP_`, `YP_`, `AP_`, `WP_`, `ZP_`, `ELP_`).
  The extractor also keeps scanning past a rejected candidate so a real
  nucleotide accession appearing later in the same `coded_by` qualifier
  still gets picked up.

## [1.2.25] — 2026-05-08

### Fixed
- **Spurious "Invalid uid" warnings from the v1.2.24 source-nuc length filter.**
  Free-text fragments such as `Q8` were occasionally extracted from
  `coded_by` qualifiers as if they were nuccore accessions, then forwarded
  to `Entrez.esummary`, where one bad UID makes NCBI reject the entire
  batch and erase the length map for every well-formed accession alongside
  it. Two fixes: (1) the accession regex in `_source_nuc_accession` now
  requires at least 3 digits, so 1- or 2-digit fragments are no longer
  picked up; (2) `fetch_nuc_lengths` defensively re-validates each UID
  against the same shape and silently skips malformed entries before the
  esummary call, so a leaked fragment can't poison a batch even if it
  somehow survives the upstream check. No config changes required.

## [1.2.24] — 2026-05-08

### Added
- **Concat-mode source-nuc length filter** — new
  `concatenation.source_nuc_min_length_frac` knob (default `0.3`). After the
  per-marker fetch, source-nucleotide accessions whose parent record is
  shorter than this fraction of the longest parent in the species' fetch
  are dropped before the `min_fraction` marker-coverage check. This filters
  out partial single-gene submissions (very common for ASFV `p72`,
  papillomavirus `L1`, etc.) that would otherwise show up as 1-marker
  pseudo-genomes and prevent any complete genome from clearing
  `min_fraction`. Adds `n_dropped_short_source_nuc` to the per-species and
  family-level fetch logs. Set to `0` to disable.

### Changed
- **`max_per_species` default for concat families bumped from 300 → 3000**
  via the smart-default path. Concat mode fetches one query per marker, so
  the cap has to cover all markers' worth of partial submissions before the
  complete-genome RefSeq proteins surface in the result set; 300 was too
  tight for popular markers like ASFV `p72`. Existing per-family YAMLs with
  an explicit `max_per_species` still win — only auto-generated configs and
  files that omit the key pick up the new default. Regenerate with
  `init-configs --force` to refresh on-disk values.

## [1.2.23] — 2026-05-08

### Added
- **Per-family `manual.include` / `manual.exclude` lists** — new optional
  block in per-family YAML configs for curator overrides on record
  selection. Both lists hold exact accessions (with version, e.g.
  `NC_002617.1`).
  - `manual.exclude` accessions are dropped immediately after fetch,
    before QC.
  - `manual.include` accessions bypass all QC (length, ambiguity,
    organism exclusion) and are protected through clustering, proportional
    merge, and length/branch-length outlier removal — i.e. stronger than
    RefSeq at QC, equal to RefSeq downstream.
  - Validation rejects overlap between the two lists, non-string entries,
    and empty strings; whitespace is stripped and duplicates deduped.
  - Two new `summary.tsv` columns: `qc_manual_include_bypassed` and
    `qc_manual_exclude_dropped`.
  - Existing per-family YAMLs without a `manual:` block continue to load
    unchanged (defaults to two empty lists).
  - In concat mode, manual.include genomes still need to clear
    `min_fraction` at fetch time to be visible to the override.

### Changed
- **Concat-mode `status.tsv`** now lists the markers actually used in the
  alignment (those retained in ≥1 genome after per-marker length-outlier
  filtering) for the `OK` row, taken from the tree_100 set when present
  and falling back to tree_500. Skip-row markers stay as the target
  preset since no analysis ran.

## [1.2.22] — 2026-05-08

### Changed
- **Leaf-label naming convention** — single canonical format applied
  consistently across sequence FASTA outputs, PhyloXML `<name>` elements,
  Newick leaves, and rendered tree images (rooted/unrooted, PDF/PNG):
  `<species>|<accession>` when host is unknown / empty, or
  `<species>|<accession>|<host>` when host is known. The strain field is
  no longer part of the label. Both single-protein and concat modes share
  one builder (`canonical_leaf_label`) so byte-identical labels are
  produced for the same `(species, accession, host)` inputs.
- **Logging hygiene** — all path-bearing INFO log lines demoted to DEBUG
  except the global sequence cache message (`Global sequence cache: ...`),
  which is the only INFO line that may print a full filesystem path.
  Affected lines include id-map written, cache invalidation, sequence-
  length plot written, overview PNG written, tree image written, and
  Newick written.

### Fixed
- **Concat PhyloXML `<name>` and `<description>`** now report the number
  and list of markers actually present in the alignment (i.e. retained in
  ≥ 1 genome after per-marker length-outlier filtering), not the full
  target marker preset. `summary.tsv` `concat_n_markers_used` is
  numerically unchanged but is now derived from the same list, so all
  three outputs (PhyloXML, summary, report) agree.

## [1.2.21] — 2026-05-08

### Added
- New `summary.tsv` columns per tree: `length_filter_median`,
  `length_filter_lo_cutoff`, `length_filter_hi_cutoff` — the resolved
  median and final keep window from the pre-MSA length-outlier filter.
  Empty in concat-mode rows (per-marker breakdown stays in the log).

### Changed
- Pre-MSA length-outlier filter logs one user-readable INFO line per call
  (kept / dropped counts, median, resolved keep window, `k`, floor),
  always emitted when the filter is enabled.

## [1.2.20] — 2026-05-07

### Changed
- Length-outlier filter adds a hard-floor safety net unioned with the MAD
  window. Defaults: `min_lo_mult: 0.20`, `max_hi_mult: 5.0`. In tight
  families where MAD on log-lengths gives a very narrow window,
  moderately-truncated legitimate variants are kept; gross outliers are
  still caught. Either knob set to `0` disables that side of the floor.

## [1.2.19] — 2026-05-07

### Changed
- Pre-MSA length-outlier filter switched from a flat
  `lo_mult: 0.333 / hi_mult: 3.0` rule to a robust **MAD-on-log-lengths**
  window: `exp(median(log L) ± k · σ_log)` with `σ_log = 1.4826 · MAD(log L)`
  (default `k: 5.0`). Adapts to each family's natural length spread.
  Old config keys (`lo_mult` / `hi_mult`) are silently ignored — re-run
  `vfam_trees init-configs --force` to regenerate per-family yaml with
  the new schema.

## [1.2.18] — 2026-05-07

### Fixed
- Concat-mode `_safe_charset_name` `UnboundLocalError` from a local
  re-import.
- Rooted-tree PDF/PNG no longer overlays the family name on the figure
  caption (clears the auto-set matplotlib title from `Phylo.draw`).
- PhyloXML `<sequence type="…">` uses the schema-valid value `"protein"`
  (was `"aa"`).

## [1.2.16] — 2026-05-07

### Added
- Concat mode publishes the same per-tree outputs as single-protein mode:
  per-tree `id_map_*.tsv`, `metadata_*.tsv`, and per-marker FASTAs
  (`<safe_marker>_raw.fasta`, `<safe_marker>_alignment.fasta`) under
  `<Family>_markers_500/` and `<Family>_markers_100/`.

## [1.2.15] — 2026-05-07

### Changed
- Concat-mode iterative outlier-removal log lines match the single-
  protein detail: include branch length, MAD ratio, and threshold per
  removed leaf.

## [1.2.14] — 2026-05-07

### Changed
- Proportional cross-species merge keeps RefSeq-bearing species first
  when the species count exceeds the target; non-RefSeq species fill any
  remaining slots. The `summary.tsv` columns
  `tree{500,100}_n_species_dropped_at_cap` and
  `tree{500,100}_n_refseq_species_dropped_at_cap` quantify both classes.

## [1.2.13] — 2026-05-07

### Changed
- Tighter MAFFT default options for both tree_500 and tree_100.

## [1.2.12] — 2026-05-07

### Fixed
- Concat-mode `support_type` correctly reflects the actual measure
  computed (was sometimes mis-labelled).
- `vfam_trees run --force` truncates `summary.tsv` / `status.tsv` for
  re-run families instead of leaving stale rows.

## [1.2.11] — 2026-05-07

### Changed
- MMseqs2 / CD-HIT cluster wrappers honor the `-t` thread setting end to
  end.

## [1.2.10] — 2026-05-07

### Fixed
- `vfam_trees run -j N` clamps `--cores` correctly when N > available
  cores; tree-image captions render reliably; per-family default tweaks.

## [1.2.7] — 2026-05-06

### Changed
- Broad log-level cleanup; concat PhyloXML metadata reaches parity with
  single-protein output (`vipr:Host`, `vipr:Year`, etc.).

## [1.2.0] — 2026-05-02

### Added
- **Multi-marker protein concatenation** for large DNA virus families
  (CONCAT_DESIGN.md). Per-marker MAFFT + trimAl, gap-padded
  concatenation, partitioned IQ-TREE on tree_100 (`-p partitions.nex
  -m MFP` so each marker gets its own ModelFinder pick) with per-
  partition models recorded in `tree100_marker_models`. Shipped marker
  presets: Poxviridae (9), Herpesviridae and 3 other herpesvirus families
  (7), Asfarviridae (6), Iridoviridae (7), Baculoviridae + Nudi/Asco (7),
  6 NCLDV families (8-marker hallmark fallback). RefSeq genomes are
  protected at every step.

## [1.1.7] — 2026-04-30

### Added
- **RefSeq absorption**: non-RefSeq sequences ≥ 0.99 identical to a
  RefSeq within the same species are absorbed into the RefSeq before
  clustering, suppressing redundant near-zero-branch cherries. RefSeqs
  themselves are never removed. Configurable per family via the
  `refseq_absorption:` block; per-tree counts (`n_refseq_absorbed`)
  recorded in `summary.tsv`.

## [1.1.6] — 2026-04-24

### Changed
- RefSeqs are protected from both pre-MSA length-outlier removal and
  post-tree branch-length outlier removal. Flagged RefSeqs stay in the
  set and a warning is logged instead.

## [1.1.5] — 2026-04-23

### Added
- Lightweight `status.tsv` (one row per family analyzed) with `family`,
  `ncbi_taxid`, `molecule_region`, `status` (`OK` or skip reason),
  `lineage`, `baltimore_class`.
- Optional family-annotation TSV (`annotation_tsv` in `global.yaml`)
  joins extra per-family columns into `summary.tsv` and `status.tsv`;
  currently surfaces `baltimore_class` (Roman numeral I–VII).
- Additional summary columns: clustering thresholds, MSA tool/options,
  trim tool/options.

## [1.1.4] — 2026-04-23

### Added
- PhyloXML metadata expansion: `vipr:Host`, `vipr:Collection_Date`,
  `vipr:Year` (parsed from collection date), `vipr:Location`,
  `vipr:Strain`, plus rank properties `vipr:Species`, `vipr:Genus`,
  `vipr:Subgenus`, `vipr:Subfamily`. Absent values omitted.
- Unrooted radial PDF/PNG tree images (`<Family>_tree_{500,100}_ur.*`)
  alongside the rooted rectangular layouts.
- Negative result caching: species with zero GenBank hits are cached
  via a `_no_results` sentinel so they aren't re-queried every run.

## [1.1.0] — 2026-04-20

### Added
- Curated DNA virus family markers expanded to ~26 families: Papillomaviridae
  L1, Parvoviridae NS1, Polyomaviridae large T antigen, Circoviridae and
  Smacoviridae Rep, Anelloviridae ORF1, Poxviridae rpo147, Baculoviridae
  / Nudiviridae / Ascoviridae lef-8, plus 8 NCLDV families (Nimaviridae,
  Hytrosaviridae, Phycodnaviridae, Mimiviridae, Marseilleviridae,
  Pandoraviridae, Pithoviridae, Medusaviridae) using DNA polymerase as
  the universal marker.

## [1.0.15] — 2026-04-20

### Added
- **trimAl column trimming** between MAFFT and tree inference (on by
  default, `-automated1`). Pre-trim length, trim tool, and trim options
  recorded per tree in `summary.tsv`.
- **UFBoot ultrafast bootstrap** (`-B 1000`) for protein tree_100 (more
  robust on divergent protein families than `--fast` + SH-aLRT).
- Generic `tree{500,100}_support_{type,min,q1,median,q3,max,iqr}`
  columns; the PhyloXML `<confidence type="…">` mirrors the same label.

### Changed
- `download.max_per_species` default raised 200 → 300.

## [1.0.14] — 2026-04-19

### Added
- Hash-based MSA / tree checkpointing: changing the sequence set, MSA
  output, or relevant config (tool / model / options) automatically
  invalidates the cache.
- Segment query accepts any of `complete sequence`, `complete genome`, or
  `complete cds` so per-segment CDS records are not missed.
- RefSeq priority in the proportional cross-species merge.

### Changed
- Clustering / MSA / tree errors raise instead of silently falling back.
- FastTree honors `model_nuc` / `model_aa` (with warning fallback).

## [1.0.13] — 2026-04-18

### Added
- `vfam_trees status` reports the current processing stage for in-
  progress families (downloading/QC, MSA, tree inference, annotating).
- `vfam_trees run --dry-run` previews per-family parameters without
  executing anything.

## [1.0.11] — 2026-04-18

### Added
- Genus / subfamily HLS leaf coloring across PDF/PNG and PhyloXML.
- Cross-family `overview_tree_100.png` thumbnail grid shaded by viral
  realm.

## [1.0.10] — 2026-04-18

### Changed
- Branch-length outlier detection switched from absolute thresholds to
  **MAD-based** (`median + factor × MAD`, default `factor: 5.0`,
  iterative).

## [1.0.7] — 2026-04-17

### Added
- DNA virus protein-marker auto-config (initial set).
- Output directory naming `<Family>_<taxid>` (e.g.
  `Asfarviridae_137992`).
- `vfam_trees cache` CLI subcommands: `clear`, `clear --all`, `stats`.

## [1.0.5] — 2026-04-16

### Added
- Optional shared sequence download cache with per-entry TTL and lock
  files for parallel-safe use across `-j N` jobs and shared filesystems.

## [1.0.1] — 2026-04-15

### Added
- Initial public release of the vfam_trees pipeline: per-family species
  discovery, GenBank fetch, quality filtering, MMseqs2 clustering,
  proportional cross-species merge, MAFFT alignment, FastTree (tree_500)
  + IQ-TREE (tree_100), LCA-based internal-node taxonomy, taxonomy-
  guided rooting, Newick + PhyloXML output, per-family PDF report,
  cross-family `summary.tsv`. SH-aLRT support, adaptive `min_length`,
  per-target MSA options.
