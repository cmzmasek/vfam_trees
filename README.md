# vfam_trees

**vfam_trees** is a bioinformatics pipeline for building maximum-likelihood phylogenetic trees for viral families. For each family it discovers species via NCBI Taxonomy, downloads sequences from GenBank, applies quality filtering, clusters per species, aligns with MAFFT, infers trees with IQ-TREE (FastTree optional), annotates internal nodes with LCA-based taxonomy, roots taxonomy-aware, and writes Newick, PhyloXML, and Auspice v2 (Nextstrain) JSON output along with PDF/PNG figures and per-run TSV summaries.

Two trees are produced per family:

- **tree_500** — broad diversity tree (up to 500 sequences, IQ-TREE `--fast` with GTR+G for nucleotide, a 50-iteration search `-n 50` with LG+G for protein; SH-aLRT support)
- **tree_100** — collapsed representative tree (up to 100 sequences, full IQ-TREE search / GTR+G or `TEST`-selected best-fit AA model; UFBoot support, SH-aLRT also computed)

## How it works

The per-family pipeline runs in seven stages. Each stage is configurable per family (see [Configuration](#configuration)); RefSeqs are protected at every step where sequences can be dropped.

### 1. Species discovery

The family TaxID is resolved from NCBI Taxonomy and every descendant species is enumerated. The species list defines the universe of organisms eligible for download. An optional `manual.restrict_to_lineages` list in the per-family config restricts this universe to one or more taxonomic lineages — each entry is a taxon at any rank (species, genus, subfamily, ...) given as a scientific name or a taxid; only species whose taxid is a descendant of at least one listed lineage proceed to download; others are skipped entirely.

### 2. Sequence download

GenBank is queried per species with RefSeq priority. The query targets the configured molecule (`whole_genome` or a marker name) and, for segmented viruses, the configured segment keyword (records lacking it in their title are dropped). When `hmm.enabled`, the marker-name query is instead replaced by a single all-proteins-per-species fetch, and markers are identified downstream by HMMER profile homology rather than by name (see [HMM-based marker identification](#hmm-based-marker-identification)); this is how most RNA-virus families get their RdRp (see [RNA virus families](#rna-virus-families-rdrp)). RefSeqs are uncapped; non-RefSeq records are limited by `download.max_per_species`.

An optional shared cache (configurable TTL, per-entry locking, negative-result caching) avoids re-downloading the same species across runs and across parallel family jobs.

Every NCBI request is retried up to 5 times, and a connection that sends nothing for 180 s counts as failed. Without that limit, a connection lost while the computer slept or changed networks left the run waiting forever. A computer that sleeps still pauses a run, so keep it awake for long runs (macOS: `caffeinate -is vfam_trees run …`, on power).

### 3. Quality filtering

Each species' sequences are filtered in this order:

0. **Manual overrides** (per-family `manual.include` / `manual.include_seq` / `manual.include_fasta_files` / `manual.exclude` / `manual.restrict_to_lineages`) — accessions in `manual.exclude` are dropped immediately; accessions in `manual.include` bypass every QC step below and are protected from removal during clustering, proportional merge, and length-outlier filtering. Match is exact, version included (e.g. `NC_002617.1`). `manual.include_seq` accepts pasted sequences (records not yet in GenBank) as a list of `{id, organism, sequence}` mappings; `manual.include_fasta_files` accepts paths to one or more FASTA files whose sequences are injected the same way. Both are injected after fetch, receive the same bypass + downstream protection as `manual.include`, and their ids must not collide with any fetched accession (not supported in concat mode). `manual.restrict_to_lineages` restricts the set of species that proceed to download to those under one or more named taxonomic lineages (see [Species discovery](#1-species-discovery)).
1. **Organism exclusion** — case-insensitive substring match against `ORGANISM`, `SOURCE`, and `DEFINITION` (joined with newlines so terms cannot straddle field boundaries). Defaults: `synthetic construct`, `metagenome`, `MAG:`, `uncultured`, `unverified`, `vector`, `recombinant`, `patent`.
2. **Ambiguity** — `max_ambiguous` cap on the fraction of `N` / `X` / IUPAC degenerate characters.
3. **Minimum length** — `min_length: null` auto-sets the threshold to 50% of the per-species median, with relaxation fallback to 40% then 30% if too few sequences pass; in all cases a hard floor of 200 bp / 100 aa applies.
4. **Length-outlier filter** (post-merge, pre-MSA) — two-sided keep window that is the **union** of (a) a MAD-on-log-lengths window and (b) a hard floor `[min_lo_mult, max_hi_mult] × median`. See [Length-outlier behaviour](#length-outlier-behaviour) for the exact form. RefSeqs and `manual.include` records flagged by the filter are kept with a warning.

### 4. Selection and clustering

- **RefSeq absorption** — non-RefSeq sequences ≥ `refseq_absorption.threshold` (default 0.99) identical to a RefSeq within the same species are absorbed into the RefSeq, suppressing redundant near-zero-branch cherries. RefSeqs themselves are never removed.
- **Adaptive per-species clustering** — MMseqs2 with binary search for an identity threshold within `[clustering.threshold_min, threshold_max]` that yields ≤ `max_reps_500` (or `max_reps_100`) representatives.
- **Proportional cross-species merge** — when the species count exceeds the target tree size, species with at least one RefSeq are kept first; remaining slots fill by rep count. Species dropped at the cap (and how many of them carried a RefSeq) are recorded in `summary.tsv` and a WARNING is logged.

### 5. Alignment and tree inference

- **MAFFT** with separate options for nucleotide vs. protein, and for tree_500 vs. tree_100 (the latter typically slower / more accurate). Nucleotide alignments run with `--adjustdirection`: GenBank records deposited on the opposite strand are reverse-complemented before alignment (otherwise they land in the tree as nonsense long branches). MAFFT orients against the first record, so when that record is on the minority strand the finished alignment is reverse-complemented as a whole onto the majority strand (ties go to the strand of the smallest id); the strand never depends on record order and MAFFT runs once. The flipped accessions are logged, also when a cached alignment is reused, and the `*_sequences_raw_*.fasta` files keep the deposited orientation.
- **trimAl** column trimming between MAFFT and tree inference: by default (`-gt 0.5`) columns with gaps in more than half of the sequences are removed; pre-trim length, tool, and options are recorded per tree. Until 1.5.x the default was `-automated1`, which on alignments with many short fragments switches to a strict mode and keeps few columns (Retroviridae RT: 58 of 657), splitting genera; `-gt 0.5` was never worse on 42 test alignments and better on 6 (CHANGELOG 1.6.0).
- **Repeatable runs** — the same input gives the same alignment and tree on every run. MAFFT runs with `--threadit 0`: its threaded iterative refinement (FFT-NS-i, which `--auto` picks for a few hundred protein sequences, and L-INS-i, the tree_100 protein mode) otherwise returns a different alignment each time, which trimAl's adaptive `-automated1` (the trimming default before 1.6.0) turned into very different column counts (939 vs 441 for the same 100 Tombusviridae RdRps). With it, only that stage runs on one thread, and the alignment is the same for any `--threads`. It costs a few minutes per protein family; nucleotide modes (`--retree 2/3`) have no refinement stage and are unaffected. IQ-TREE runs with `-seed 42`, so a tree is repeatable at the same `--threads`. Options that set `--threadit` or `-seed` themselves win.
- **IQ-TREE for tree_500** (SH-aLRT support, auto-added by the wrapper). Nucleotide alignments use `--fast`: the same wall time as FastTree, a clearly better likelihood, and a full search gained nothing on them. Protein alignments use `-n 50`, IQ-TREE's normal search stopped after 50 iterations: on weak-signal protein sets (e.g. 500 RdRps) `--fast` lands in a different local optimum on each seed, with several strongly supported (SH-aLRT ≥ 95) splits that a full search contradicts, while `-n 50` comes about as close to the full search as a second full search does, at ~12× the `--fast` time (benchmark in CHANGELOG 1.6.0). FastTree remains available: set `tree_500.tool: fasttree` together with `options_nuc: ""` and `options_aa: ""` (config validation rejects IQ-TREE options such as `--fast` and `-n` under FastTree, and FastTree options under IQ-TREE). `+G` / `+GAMMA` enables discrete-gamma rate variation.
- **IQ-TREE full search** for tree_100 with `-alrt 1000 -B 1000 --bnni`: SH-aLRT and UFBoot (ultrafast bootstrap, with `--bnni` to curb overestimation under model violation). The outputs carry the UFBoot value; the IQ-TREE `.treefile` under `_work/100/tree/` keeps both as `alrt/ufboot` labels. IQ-TREE 2.3.2 can abort in its post-bootstrap consensus step on particular taxon sets; the wrapper then scores the ML tree that run already found with SH-aLRT (a fixed-topology run that takes seconds; ModelFinder's pick is read from the aborted run's log) instead of failing the family, and the outputs report `SH_aLRT`. A `WARNING` in the family log says when this happened. Setting `model_aa: TEST` triggers ModelFinder; the chosen model (e.g. `LG+I+G4`) is parsed from the IQ-TREE log and surfaces in `summary.tsv`, the per-family PDF, and the PhyloXML provenance.
- **Branch-support measure** is picked per tree (`SH_like` / `SH_aLRT` / `UFBoot` / `Bootstrap`) and recorded uniformly: `tree{500,100}_support_type` plus generic `support_{min,q1,median,q3,max,iqr}` columns; the PhyloXML `<confidence type="…">` attribute mirrors the same label. When several IQ-TREE measures are requested (any of `-alrt`/`--alrt`, `-lbp`, `-abayes`, `-B`/`--ufboot`, `-b`/`--boot`), each node keeps one value: the bootstrap-type one if present, else SH-aLRT. `_work/<label>/tree/tree_<label>.tree_info.json` records the tool, options and the measure the tree actually carries.
- **Iterative branch-length outlier removal** — between iterations, leaves with terminal branch length exceeding `median + factor × MAD` are removed and MSA + tree are re-run (up to `max_iterations`; only when ≥ `min_seqs` remain). With `quick_check: true` (default) each round is first checked on a quick IQ-TREE `--fast` tree (GTR+G / LG+G instead of ModelFinder); the configured search — tree_100's UFBoot search, tree_500's `-n 50` — runs once that finds nothing more to remove, and its tree is checked as before, so the final tree is always the configured one. Long-branch outliers are the same on both trees; the full search no longer repeats in every round (Togaviridae: 4 tree_100 searches of ~14 min each before). Zero-length tips (identical sequences, ≤ 1e-5 subs/site; IQ-TREE writes them as 1e-6) are left out of the median and MAD, so families with many identical genomes don't get an artificially low threshold. RefSeqs are protected: a flagged RefSeq stays with a warning. Detailed per-outlier log lines include length, ratio to median, and threshold.
- **Multi-marker protein concatenation** — large DNA virus families with curated marker presets concatenate per-marker MAFFT + trimAl alignments into a single matrix; tree_100 uses partitioned IQ-TREE (`-p partitions.nex -m MFP`) so each marker gets its own ModelFinder pick. Markers are found either by GenBank name (per-marker query) or, when `hmm.enabled`, by HMMER profile homology over one all-proteins-per-species fetch (see [HMM-based marker identification](#hmm-based-marker-identification)); either way the returned proteins are bucketed by source nucleotide accession (Policy A). A length filter (`concatenation.source_nuc_min_length_frac`, default 0.3) drops source-nuc accessions shorter than that fraction of the longest parent so partial single-gene submissions don't crowd out genome-scale proteins under `download.max_per_species` (default 3000 for concat families). See [CONCAT_DESIGN.md](CONCAT_DESIGN.md) for design details.

- **Phylogenetic-signal check** — every tree is measured once it is final: the median tip-to-tip distance on the tree (substitutions per site, the classic saturation measure), the median pairwise identity of the trimmed alignment (over shared A/C/G/T or standard amino-acid positions; ambiguity codes and X ignored), the core fraction (residues a typical sequence keeps after trimming) and the number of branches at IQ-TREE's 10 subs/site cap. A verdict of `ok`, `weak` or `saturated` goes to `summary.tsv` (`tree{500,100}_signal`, `_patristic_median`, `_identity_median`, `_core_frac`, `_branches_at_cap`), the report PDF and the log, with a `WARNING` naming the reasons when a tree is not `ok`. Nucleotide trees are `saturated` above a median of 1.5 subs/site (≈ 35% identity under JC; unrelated aligned sequence reaches about that) or with a branch at the cap and a median above 0.75, and `weak` above 0.75 or below 60% identity or 25% core; protein thresholds are 4.0 / 3.0 and 30%. A saturated nucleotide tree's deep structure is noise; a conserved protein marker carries it far better. Distance and identity are measured on the columns the tree was built from, so they depend on trimming: a stricter trim keeps the conserved columns and reads less saturated. Compare verdicts only between runs with the same `msa_trim` settings.

### 6. Annotation and rooting

- **LCA-based internal-node annotation** using NCBI ranked lineages, with a configurable lineage-depth filter (`taxonomy.lca_min_rank`, default `species`) so leaves whose lineage is too shallow are excluded from the LCA vote.
- **Genus inference** for taxa lacking a formal `genus` rank in NCBI lineage: `none`, `suffix` (single-word taxa ending in `virus` treated as genus, per ICTV convention), or `deepest` (default — suffix first, then deepest rank above species).
- **Taxonomy-guided rooting** via LCA specificity scoring, with MAD and midpoint fallbacks. The PhyloXML `<phylogeny>` is emitted with `rooted="true" rerootable="false"` so downstream viewers cannot re-root the carefully rooted output.
- **Genus / subfamily coloring** in HLS color space — one hue band per subfamily, lightness varies across genera within a subfamily; colors are applied to PDF/PNG figures and the Auspice JSON, and optionally to PhyloXML (`style:font_color`, `output.phyloxml_font_colors`).

### 7. Output

Per-family: Newick, PhyloXML (with `<taxonomy>`, `<sequence>`, tip `<date>`, `vipr:` metadata, and rank properties), a ready-to-open Archaeopteryx.js HTML page and an Auspice v2 (Nextstrain) JSON for interactive exploration, rooted-rectangular and unrooted-radial PDF + PNG tree images, a topology-only icon PNG, sequence and metadata FASTAs/TSVs, a per-family PDF report, a per-family log, the merged settings the family ran with (`<Family>.yaml`), and a list of species that have NCBI records but not the marker (`<Family>_species_without_marker.tsv`).

Cross-family: a row-per-family `summary.tsv` with full statistics, a lightweight `status.tsv` (OK, no tree, or skip reason), an `overview_tree_100.png` thumbnail grid shaded by viral realm, and a record of each `vfam_trees run` (`run_records/`: version, command, tool versions, effective global settings). See [Output](#output) for the full file inventory.

Per-family directories are named `<Family>_<taxid>` (e.g. `Asfarviridae_137992`) — or, when `manual.name` is set, `<manual.name>_<taxid>` (e.g. `Orthoebolavirus_11266`), with every output file prefixed by the same display name (see [Output prefix and directory naming](#output-prefix-and-directory-naming)). Failures and skips at any stage produce a row in `status.tsv` with the skip reason; the per-family work directory is preserved for inspection.

## Other capabilities

- **Interactive trees (Auspice / Nextstrain / Nextclade)** — every family also emits a schema-valid Auspice v2 JSON per tree (`<prefix>_tree_{500,100}_auspice.json`) loadable in the Nextstrain viewers: [auspice.us](https://auspice.us) (drag-and-drop), Nextclade, or a local `auspice view`. These are *divergence* trees (branch lengths are substitutions per site, not time-resolved); the "color by genus" view reuses the same palette as the PDF/PhyloXML output, and genus / subfamily / species / host / country / location / year (plus an **Outbreak** trait for `additional_data_sources` tips) are available as colorings and filters. Toggle with the `output.auspice_json` per-family flag (default on).
- **Interactive trees (Archaeopteryx / Archaeopteryx.js)** — the PhyloXML is written for the [Archaeopteryx](https://github.com/cmzmasek/forester) desktop viewer and [Archaeopteryx.js](https://github.com/cmzmasek/archaeopteryx-js): `branch_length_unit="subs/site"` (scale bar and clock-plot rate units), each tip's collection date as a decimal-year `<date unit="year">` (the viewers' root-to-tip clock plot needs every tip dated), ICTV species in `vipr:Species` consistent with the internal labels (so "Clades by rank" groups a species once), the leaf's own NCBI taxon (id, current name, rank) in `<taxonomy>`, NCBI ids on internal taxonomies, `vipr:Country` for grouping, and accession sources (`ncbi` / `pathoplexus` with a record `<uri>` / `user`). Each tree also gets `<prefix>_tree_{500,100}.html`: the tree inlined in a page that opens it in Archaeopteryx.js 3.28.0 loaded from jsDelivr — double-click to explore (needs an internet connection). Toggle with `output.archaeopteryx_html` (default on).
- **Custom output naming** — set `manual.name` to rename a family's outputs (e.g. run *Filoviridae* restricted to the *Orthoebolavirus* genus and label everything `Orthoebolavirus`); the directory and all files take the sanitized display name while NCBI/taxonomy/cache lookups keep using the biological family. See [Output prefix and directory naming](#output-prefix-and-directory-naming).
- **Additional data sources (outbreak sequences)** — supplement a family tree with forced-include sequences from [Pathoplexus](https://pathoplexus.org) via the `additional_data_sources` config block. These bypass QC and clustering (like `manual.include`) but carry real host/location/date/taxon metadata, so they are coloured by genus. Useful for grafting the latest outbreak genomes onto an existing tree. See [Additional data sources (Pathoplexus)](#additional-data-sources-pathoplexus).
- **HMM-based marker identification** — identify marker proteins by **HMMER profile homology** instead of GenBank protein names, via the per-family `hmm:` config block: by default the RdRp, RT and Rep of the RNA-virus, reverse-transcribing and CRESS-DNA families (see [RNA virus families](#rna-virus-families-rdrp)), and markers of large DNA virus families with inconsistent annotation. Each marker carries an `hmms:` token — a single curated Pfam-A profile (`"DNA_pol_B"`) or a multidomain architecture (`"Adeno_hexon--Adeno_hexon_C"`, HMMs in N-to-C order). Hits are gated on the curated Pfam GA bit-score (or an E-value; the RdRp, RT and Rep presets use E ≤ 1e-3) plus HMM-model coverage; the all-proteins-per-species pool is grouped by source genome and assigned to markers by HMM. Works on both the single-marker and concatenated-marker paths. Bundled CC0 Pfam-A sets ship for the **RdRp** (23 profiles), **RT** (4), **Rep** (4), **Poxviridae** (8 markers), **Adenoviridae** (hexon), **Mamonoviridae** (DNA polymerase B) and **NCLDV** (9 core-gene profiles for Asfarviridae, Pithoviridae, Mimiviridae and Marseilleviridae); requires HMMER on `$PATH`. See [HMM-based marker identification](#hmm-based-marker-identification).
- **Pre-configured family presets** — 134 families on an HMM-identified marker protein (107 RNA-virus families on the RdRp, 4 reverse-transcribing families on the RT, 23 CRESS-DNA and satellite families on Rep), 24 DNA virus families (sequence type, region/marker, and concat marker sets for the 15 families that use concatenation) and 28 segmented RNA virus families (segment keywords, used when such a family runs as nucleotide; 21 of them are on the RdRp by default) ship with curated defaults.
- **Checkpointing** — MSA and tree steps store a content-hashed sidecar of inputs + tool / model / options (and, for trees, the tool version), so any change, including a tool upgrade, auto-invalidates the cache. Each iterative outlier-removal pass gets its own hash; partially-completed runs resume in the right place.
- **Family-annotation TSV** (optional) — joins extra per-family columns (e.g. `baltimore_class`) into `summary.tsv` and `status.tsv`. Missing file or family is silent.
- **Stage tracking** — `vfam_trees status` reports the current processing stage for in-progress families (downloading/QC, MSA, tree inference, annotating), and `done`, `no tree` or `skipped` (with the reason) for finished ones.
- **Dry-run mode** — `vfam_trees run --dry-run` previews per-family parameters without executing.
- **Parallelism** — `-j N` runs N families concurrently; `-t T` sets threads per family job. The download cache and external tools are concurrency-safe.
- **Validation and warnings** — partial NCBI batches, config keys this version does not read (in `global.yaml` and family files, at any depth), configs that override a family preset, weak or saturated trees (signal check), species without the marker, and MAFFT / tree output sanity checks all surface as log warnings.

## Dependencies

### Python packages

```
biopython  >= 1.81
click      >= 8.1
pyyaml     >= 6.0
matplotlib >= 3.9    # PDF report and tree images; requires NumPy 2.x compatible build
```

### External tools

| Tool | Purpose |
|------|---------|
| `mafft` | Multiple sequence alignment |
| `trimal` | Alignment column trimming |
| `FastTree` | Optional alternative tree tool (`tool: fasttree`); not needed by the defaults, and `vfam_trees test` reports it as optional |
| `iqtree2` | ML tree inference (tree_500 and tree_100) |
| `mmseqs` | Sequence clustering |
| `hmmscan`, `hmmpress` (HMMER) | HMM-based marker identification: the RdRp of most RNA-virus families, the RT of the reverse-transcribing families, the Rep of the CRESS-DNA families, Poxviridae, Adenoviridae, Mamonoviridae, Asfarviridae, Pithoviridae, Mimiviridae and Marseilleviridae markers. Required for the default family list |

All required tools must be available on `$PATH`. Installation via conda is recommended:

```bash
conda install -c bioconda mafft iqtree mmseqs2 trimal hmmer
conda install -c bioconda fasttree     # only needed for tool: fasttree
```

## Installation

```bash
git clone https://github.com/cmzmasek/vfam_trees.git
cd vfam_trees
pip install -e .
```

## Quick start

```bash
# 1. Create and edit global config (sets NCBI email and API key)
vfam_trees init          # or: vfam_trees init --full (every setting listed)
# edit config/global.yaml

# 2. Check all dependencies
vfam_trees test

# 3. Generate per-family configs (review before running)
vfam_trees init-configs -f families.txt

# 4. Preview what will run without executing anything
vfam_trees run -f families.txt --dry-run

# 5. Run the pipeline
vfam_trees run -f families.txt -j 4 -t 4

# 6. Check progress
vfam_trees status -f families.txt
```

## Configuration

### Global config (`config/global.yaml`)

Generate a template with `vfam_trees init` (common settings) or `vfam_trees init --full` (every setting this version has, each with its default and a comment). Then edit it to set your NCBI credentials:

```yaml
ncbi:
  email: your.email@example.com     # REQUIRED
  api_key: your_ncbi_api_key        # optional but recommended (10 req/s vs 3 req/s)
```

An NCBI API key can be obtained for free at https://www.ncbi.nlm.nih.gov/account/.

The `defaults:` section of `global.yaml` overrides the built-in defaults for all families. Per-family configs further override individual parameters.

Settings are layered, lowest to highest priority: the built-in defaults, the `defaults:` section of `global.yaml`, the family preset (RdRp marker, DNA marker proteins, segment), and the family's own `configs/<Family>.yaml`. Nothing is migrated when a new version changes a default: a value written in `global.yaml` or a family file stays in force. A file written by `init --full` therefore pins every default of the version that wrote it; delete a block to follow the built-in default again.

Keys this version does not read (misspelled, or left over from an older version, such as `length_outlier.hi_mult`) are reported with a WARNING, at any depth, both for `global.yaml` and for family files. The output folder and log level are command-line options of `vfam_trees run` (`-o`, `--log-level`), not settings in `global.yaml`; older templates wrote `output_dir:` and `log_level:` there, which never had an effect.

Every run records what it used; see [Run records](#run-records).

#### Sequence download cache

```yaml
cache:
  dir: ~/.vfam_cache    # shared across all runs on this machine (or a lab filesystem)
  ttl_days: 90          # re-download after 90 days; null = never expire
```

Cache entries are keyed by `(taxid, db, region, segment, max_per_species)` so changing any query parameter triggers a fresh download. Parallel family jobs (`-j N`) coordinate via per-entry lock files. Negative results (species with zero hits) are cached too, so they are not re-queried on every run.

Management:

```bash
vfam_trees cache clear Asfarviridae
vfam_trees cache clear --all          # wipe entire cache
vfam_trees cache stats                # show entry count and size
```

#### Family-annotation TSV

```yaml
annotation_tsv: virus_families_annotation.tsv   # null → this name, next to global.yaml
```

Relative paths are resolved against the `global.yaml` directory. The file must have a `family` column (case-insensitive match) plus any extra columns to be picked up — currently `baltimore_class` (Roman numeral I–VII per Baltimore 1971) is the only one read by the pipeline. Missing file, missing key, or missing family → column simply left empty.

```tsv
family	baltimore_class	host_range	segmented	genome_size
Flaviviridae	IV	Vertebrates	No	~11 kb
Poxviridae	I	Vertebrates	No	~130–375 kb
```

### Per-family configs (`configs/<Family>.yaml`)

Per-family configs are auto-generated if missing. Generate them in advance to review and tune parameters before running:

```bash
vfam_trees init-configs -f families.txt

# Regenerate with current defaults (overwrites any manual edits):
vfam_trees init-configs -f families.txt --force
```

Key parameters:

```yaml
download:
  max_per_species: 300          # cap on non-RefSeq sequences per species

sequence:
  type: nucleotide              # nucleotide or protein
                                # auto-set to protein for large DNA virus families
  region: whole_genome          # whole_genome, or a marker name (e.g. "DNA polymerase", "hexon")
                                # auto-set for known large DNA virus families
  segment: null                 # segment keyword for segmented viruses (e.g. "segment L")
                                # auto-set for known segmented RNA families

quality:
  min_length: null              # null = auto (50% of per-species median, floor 200 bp/100 aa)
  max_ambiguous: 0.01           # maximum fraction of ambiguous bases/residues
  exclude_organisms:            # case-insensitive substring match against
    - synthetic construct       # ORGANISM + SOURCE + DEFINITION (joined with newline
    - metagenome                # so terms cannot straddle field boundaries)
    - "MAG:"                    # metagenome-assembled genomes (NCBI DEFINITION prefix)
    - uncultured
    - unverified
    - vector
    - recombinant
    - patent

clustering:
  tool: mmseqs2
  threshold_min: 0.70           # minimum clustering identity
  threshold_max: 0.99           # maximum clustering identity
  max_reps_500: 20              # max representatives per species for tree_500
  max_reps_100: 5               # max representatives per species for tree_100

targets:
  max_500: 500                  # target sequences for tree_500
  max_100: 100                  # target sequences for tree_100

msa_500:
  tool: mafft
  options_nuc: "--retree 2 --adjustdirection"   # nucleotide; --adjustdirection flips opposite-strand records
  options_aa: "--auto"                          # amino acid (MAFFT auto-selects strategy)

msa_100:
  tool: mafft
  options_nuc: "--retree 3 --adjustdirection"     # nucleotide
  options_aa: "--maxiterate 1000 --localpair"     # amino acid (MAFFT L-INS-i; high accuracy, slower)

msa_trim:
  enabled: true                  # drop poorly-aligned columns before tree inference
  tool: trimal
  options: "-gt 0.5"             # keep columns with at most 50% gaps; nucleotide and protein

tree_500:
  tool: iqtree                # FastTree: tool: fasttree plus options_nuc/options_aa: ""
  options_nuc: "--fast"       # nucleotide; SH-aLRT support auto-added by the wrapper
  options_aa: "-n 50"         # amino acid (marker proteins): search stopped after 50 iterations
  model_nuc: GTR+G
  model_aa: LG+G              # LG+G used for amino acid sequences

tree_100:
  tool: iqtree
  options_nuc: "-alrt 1000 -B 1000 --bnni"   # full search; outputs carry UFBoot, SH-aLRT stays in the .treefile
  options_aa: "-alrt 1000 -B 1000 --bnni"
  model_nuc: GTR+G
  model_aa: TEST               # TEST = IQ-TREE ModelFinder; the chosen best-fit model
                               # is recorded in summary.tsv / PDF / PhyloXML instead of "TEST"

length_outlier:
  enabled: true                 # pre-MSA length-based outlier removal
  k: 5.0                        # MAD-on-log keep window: exp(median(log L) ± k · σ_log)
                                # with σ_log = 1.4826 × MAD(log L); 0 disables MAD
  min_lo_mult: 0.20             # never drop seqs ≥ this × median (0 disables lower floor)
  max_hi_mult: 5.0              # never drop seqs ≤ this × median (0 disables upper floor)

outlier_removal:
  enabled: true                 # iterative post-tree branch-length outlier removal
  factor: 20.0                  # threshold = median + factor × MAD (Median Absolute Deviation)
  max_iterations: 3             # maximum MSA+tree iterations
  min_seqs: 40                  # only remove outliers when ≥ min_seqs sequences remain after removal
  quick_check: true             # check each round on a quick --fast tree; the configured
                                # search runs once no more outliers turn up

labeling:
  format: "{species}|{id}|{host}"  # format string for leaf labels (PhyloXML <name>,
                                    # PDF/PNG tree images, FASTA display names).
                                    # Placeholders: {species}, {id} (accession), {host},
                                    # {strain}, {location}, {year}, {genus}.
                                    # Literal text (separators etc.) is kept verbatim.
  replace_whitespace: true          # replace spaces in field values with underscores
  keep_separator_on_empty: false    # false: drop empty fields and their preceding
                                    #   separator so no leading/consecutive separators appear
                                    # true: keep separators regardless of field content

coloring:
  genus_inference: deepest      # none: only formal NCBI genus-rank entries are colored
                                # suffix: single-word taxa ending in "virus" treated as genus
                                #         (recovers ICTV genus names that NCBI still has at "no rank")
                                # deepest (default): suffix first, then deepest lineage entry
                                #         above species rank

taxonomy:
  lca_min_rank: species         # none: every leaf contributes to internal-node LCA voting
                                # subfamily / genus / species (default): exclude leaves whose
                                #   lineage does not reach this rank, so shallow lineages cannot
                                #   drag ancestor labels back toward the root

output:
  auspice_json: true            # also emit an Auspice v2 (Nextstrain) JSON per tree
                                # (<prefix>_tree_{500,100}_auspice.json) for interactive
                                # viewing at auspice.us; divergence tree only (branch lengths
                                # are substitutions/site, so it is not time-resolved)
  phyloxml_font_colors: false   # write the genus palette as style:font_color on PhyloXML
                                # tips (off: the Archaeopteryx viewers open a tree with their
                                # own Color-by, which hides style colours, and these have no legend)
  archaeopteryx_html: true      # also write <prefix>_tree_{500,100}.html, which opens the tree
                                # in Archaeopteryx.js (viewer loaded from cdn.jsdelivr.net)

manual:                         # curator overrides on per-family record selection
  include: []                   # exact accessions (with version, e.g. "NC_002617.1") that
                                # bypass all QC (length, ambiguity, organism exclusion) and
                                # are protected through clustering, proportional merge, and
                                # length-outlier filtering — stronger than RefSeq at QC,
                                # equal to RefSeq downstream
  include_seq: []               # pasted sequences not yet in GenBank — each entry is a
                                # {id, organism, sequence} mapping; injected after fetch,
                                # bypass QC, and receive the same downstream protection as
                                # manual.include records (not supported in concat mode)
  include_fasta_files: []       # paths to FASTA files whose sequences are injected identically
                                # to include_seq entries; FASTA id field → sequence id,
                                # remainder of header → organism (not supported in concat mode)
  exclude: []                   # exact accessions to drop immediately after fetch, before QC
  restrict_to_lineages: []      # restrict the pipeline to species under one or more named
                                # lineages — entries are taxon names at any rank (species,
                                # genus, subfamily, ...) or NCBI taxids (integer or digit
                                # string); only species descendant from at least one listed
                                # lineage proceed to download; empty or absent means use the
                                # full discovered species list
  name: ""                      # override the display name in PDF/PNG titles and the PhyloXML
                                # <name> element; e.g. "Hantaviridae 2026". When set it ALSO
                                # becomes the output prefix: every output file and the output
                                # directory are renamed to a filesystem-safe form of it (the
                                # directory keeps its _<taxid> suffix). NCBI queries, taxonomy,
                                # caching, and the summary.tsv 'family' column still use the
                                # biological family. Empty (default) uses the family name.
                                # See "Output prefix and directory naming" below.

additional_data_sources: []     # supplement the tree with forced-include sequences from an
                                # external database (currently Pathoplexus). Each entry is a
                                # mapping; see "Additional data sources (Pathoplexus)" below.
                                # Bypass QC + clustering like manual.include; nucleotide-only;
                                # not supported in concat mode.
```

The `coloring`, `taxonomy`, and `output` keys can also be set globally in the `defaults:` section of `global.yaml` — per-family configs inherit them automatically.

#### Output prefix and directory naming

By default a family's outputs live in `results/<Family>_<taxid>/` with every file prefixed `<Family>_…`. Setting `manual.name` overrides this prefix: the value is sanitized to a filesystem-safe form (runs of any character outside `[A-Za-z0-9._-]` collapse to `_`) and used for both the directory and every output file. The directory keeps its `_<taxid>` suffix for traceability.

For example, running `Filoviridae` restricted to the *Orthoebolavirus* genus with `manual.name: Orthoebolavirus`:

```
results/Orthoebolavirus_11266/
  Orthoebolavirus_tree_100.nwk
  Orthoebolavirus_tree_100.xml
  Orthoebolavirus_tree_100.html
  Orthoebolavirus_tree_100_auspice.json
  Orthoebolavirus_report.pdf
  ...
```

The biological family name is still used everywhere identity matters — NCBI species discovery and downloads, taxonomy resolution, the sequence cache, the `summary.tsv` / `status.tsv` `family` column, and `vfam_trees status` reporting. Each output directory carries a hidden `.family` marker recording the biological family, so `status` (and the other CLI subcommands) locate a renamed directory without re-reading config. If two families would resolve to the same directory (only possible when their `manual.name` values are identical *and* the NCBI taxid cannot be resolved), the run aborts with a clear error rather than overwriting the other family's outputs.

Families **without** `manual.name` are completely unaffected: the sanitized prefix equals the family name, so output paths are identical to earlier versions.

Example — adding a pasted reference and forcing a curator-picked accession in:

```yaml
manual:
  include:
    - NC_002617.1               # force-keep this RefSeq even if QC would drop it
  include_seq:
    - id: MyLab_Isolate_2026
      organism: Foo virus
      sequence: |
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    - id: PreprintSeq_42
      organism: Bar virus
      sequence: ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  include_fasta_files:
    - /data/unpublished/batch_2026_isolates.fasta
    - /data/unpublished/outgroup.fa
  exclude:
    - KF234567.1                # drop a known-bad isolate before QC
```

Notes on `include_seq` and `include_fasta_files`:

- `id` (for `include_seq`) or the first whitespace-delimited token of the FASTA header (for `include_fasta_files`) is used throughout the pipeline (FASTA headers, tree leaves, PhyloXML) and must not collide with any accession returned by NCBI for the family, with anything in `include` / `exclude`, or across the two options — collisions raise a hard error.
- `organism` (for `include_seq`) or the remainder of the FASTA header after the id (for `include_fasta_files`) controls the species bucket the entry joins (matching a fetched species name appends to that bucket; a novel name creates a new one) and shows up in leaf labels. If a FASTA header is a single token (no remainder), the organism is left empty and the leaf label is just the header itself.
- For `include_seq`, the sequence may be wrapped across lines — whitespace is stripped and the sequence is uppercased.
- Injected leaves render gray (no genus) and do not participate in LCA voting, since only id and organism are required.
- Not supported when `sequence.region == "concatenated"` — injected sequences cannot be split into per-marker proteins.
- The sequences must be of the family's type: a nucleotide sequence in a protein family, or a protein in a nucleotide family, is a hard error. Likewise a `manual.include` accession whose shape clearly belongs to the other NCBI database (a nuccore accession for a protein family, a protein accession for a nucleotide one; case does not matter) is left out with a warning instead of being fetched. Ids of no NCBI shape, such as Pathoplexus ids, are not judged by shape.

#### Additional data sources (Pathoplexus)

`additional_data_sources` pulls *forced-include* sequences from an external pathogen database and injects them into the tree exactly like `manual.include` records — they **bypass QC and clustering/subsampling and are protected** from length-/branch-outlier removal. The typical use is adding the latest **outbreak** genomes to an existing family tree. Unlike pasted sequences, these carry real host / location / collection-date / taxon metadata, so they are labelled and **coloured by genus** like fetched records.

Currently the only `source` is **Pathoplexus** ([pathoplexus.org](https://pathoplexus.org)), queried through its GenSpectrum LAPIS API. Pathoplexus organises data into a fixed set of curated *organisms* addressed by slug (not by taxid), so a genus-level tree lists one entry per organism — e.g. the *Orthoebolavirus* genus spans `ebola-zaire`, `ebola-sudan`, and `ebola-bdbv`:

```yaml
additional_data_sources:
  - source: pathoplexus
    organism: ebola-zaire         # Pathoplexus organism slug (required)
    country: Uganda               # optional, exact (LAPIS geoLocCountry)
    host: Homo sapiens            # optional, exact (LAPIS hostNameScientific)
    date_from: "2022-01-01"       # optional ISO date (collection-date lower bound)
    date_to: "2022-12-31"         # optional ISO date (collection-date upper bound)
    max_seqs: 200                 # optional cap (default 200) — these skip subsampling
    dedup_vs_ncbi: true           # optional (default true) — drop records already pulled
                                  #   from NCBI (matched by INSDC accession)
    name_prefix: "PATHOPLEXUS_"   # optional — prepended verbatim to each leaf label
    outbreak_name: "Bdbv-2026"    # optional — tags tips (vipr:outbreak + Auspice trait)
  - source: pathoplexus
    organism: ebola-sudan
```

Valid organism slugs: `andv`, `cchf`, `dengue`, `ebola-bdbv`, `ebola-sudan`, `ebola-zaire`, `hmpv`, `marburg`, `measles`, `mpox`, `rsv-a`, `rsv-b`, `west-nile`, `yellow-fever`.

Notes:

- **Nucleotide-only**, and the family's tree must be nucleotide: a family analysed as protein (including the RNA-virus families on the RdRp preset, e.g. Filoviridae) is refused at config load. Set `sequence.type: nucleotide` in its config first (see [RNA virus families](#rna-virus-families-rdrp)). Not supported when `sequence.region == "concatenated"`.
- `name_prefix` is prepended **verbatim** to every leaf label from that source (supply your own trailing separator, e.g. `PATHOPLEXUS_`), so external tips are easy to spot in the tree, FASTA, `*_id_map.tsv`, PhyloXML, and Auspice.
- `outbreak_name` tags each injected tip with the given label, emitted as a PhyloXML `vipr:outbreak` property and an **Outbreak** colouring + filter in the Auspice JSON — handy for isolating the outbreak clade in Nextstrain. Both fields appear only when set, so ordinary trees are unaffected.
- Only the **latest version** of each Pathoplexus accession is fetched; revoked records are excluded.
- The Pathoplexus accession (e.g. `PP_000PHDC.3`) is used as the sequence id and must not collide with a fetched NCBI accession.
- The fetch happens at pipeline runtime; the config is validated at load (known source, valid organism slug, ISO dates, positive `max_seqs`). A network failure aborts that family with a clear error rather than silently dropping the requested sequences.

Example — restricting a family run to specific lineages (mixing ranks, names, and taxids):

```yaml
manual:
  restrict_to_lineages:
    - Orthohantavirus               # genus — all species under it are kept
    - Mobatvirus                    # another genus in the same family
    - Hantaan orthohantavirus       # individual species also accepted
    - 11595                          # NCBI taxid (integer) — equivalent to the species name above
    - "1980519"                      # taxid as a quoted digit string — also accepted
```

For each entry, NCBI Taxonomy is queried for the species-rank descendants under that taxon; the family's discovered species list is then restricted to those that fall under at least one of the listed lineages. Entries whose subtree does not overlap the discovered species produce a WARNING in the log (typo, or a taxon outside this family). If `restrict_to_lineages` matches nothing in the discovered species list the family is skipped with an explanatory status entry.

#### Length-outlier behaviour

Just before alignment, vfam_trees removes sequences whose **length** is wildly out of line with the rest of the family — tiny fragments and unusually long records (often mis-annotated or concatenated entries). A handful of these distort the whole multiple-sequence alignment far out of proportion to their number, so they are dropped first. (RefSeq and `manual.include` records are never dropped here; if flagged they are kept with a warning.)

The filter keeps a **window of acceptable lengths** and removes anything outside it. The window is designed to be *generous by default but adaptive when the data call for it*. It is the **wider** of two ranges on each side:

1. **A fixed safety range** — by default **`0.20× to 5×` the median length**. Anything in this range is always kept; this is what stops the filter from over-trimming.
2. **An adaptive range** — from how spread out the family's lengths actually are (a robust, outlier-resistant estimate on log-lengths). When a family genuinely spans a broad range of lengths, this **widens** the window so real long/short members survive. It can only widen the window, never shrink it below the fixed range.

A sequence is removed **only if it falls outside both** — that is, only genuine extremes are cut.

*Why combine them?* If you used only the adaptive range, a very uniform family (almost every sequence the same length) would have a near-zero spread, collapsing the window to a sliver and wrongly flagging anything slightly off. The fixed floor prevents that. If you used only the fixed range, a genuinely diverse family would lose real members; the adaptive range widens it when warranted. Combining them gives both.

**Exact form.** Keep a sequence of length *L* when `lo ≤ L ≤ hi`:

```
lo = min( 0.20 × median ,  exp( median(log L) − k·σ ) )      # k = 5 by default
hi = max( 5.0  × median ,  exp( median(log L) + k·σ ) )
     └─── fixed floor ───┘  └──── adaptive (robust, on log-lengths) ────┘
σ = 1.4826 × MAD(log L)      # a robust stand-in for the std. dev. of log-lengths
```

Per-family knobs: `length_outlier.k` (default 5; set `0` to disable the adaptive range), `min_lo_mult` / `max_hi_mult` (default 0.20 / 5.0 — set either to `0` to disable that side of the floor).

**Worked examples** (median length = 1000, defaults `k=5`, floor `0.20×–5×`):

| Family's length distribution | Keep window | What happens |
|---|---|---|
| Very uniform (all ≈ 1000) | `200 – 5000` | Adaptive range ~vanishes; the floor governs. A 400-nt record (0.40×) is **kept** |
| Broadly spread (~2× range) | `≈190 – ≈5230` | Adaptive range widens past the floor to admit the family's real long/short members |
| Half the sequences tied at the median | `200 – 5000` | Robust spread = 0, so the floor alone applies |
| Uniform bulk + one 30,000-nt record (30×) | `200 – 5000` | 30,000 is outside the window → **dropped** |
| Uniform bulk + one 30-nt fragment (0.03×) | `200 – 5000` | 30 is outside the window → **dropped** |

The outcome — how many long and short outliers were dropped, the median, and the resolved `lo`/`hi` — is logged once per tree and written to `summary.tsv` (`tree{500,100}_n_length_outliers_{long,short}`, `length_filter_median`, `length_filter_lo_cutoff`, `length_filter_hi_cutoff`).

#### DNA virus families

Known DNA virus families are automatically configured with curated ICTV-aligned markers (the CRESS-DNA families and Hepadnaviridae are on the HMM-identified Rep and RT markers, [below](#reverse-transcribing-and-cress-dna-families-rt-rep)). Small ssDNA / dsDNA families use a single diagnostic protein (when one is well established), medium–large dsDNA families use a family-specific structural or replication protein, and nucleocytoplasmic large DNA viruses (NCLDVs) share DNA polymerase (family B) as a universal marker:

| Family group | Marker gene | Sequence type |
|---|---|---|
| Smacoviridae | Rep (by GenBank name) | protein |
| Anelloviridae | ORF1 | protein |
| Parvoviridae | NS1 | protein |
| Polyomaviridae | large T antigen | protein |
| Papillomaviridae | L1 | protein |
| Adenoviridae | hexon | protein |
| Orthoherpesviridae, Alloherpesviridae, Herpesviridae | DNA polymerase | protein |
| Malacoherpesviridae | major capsid protein | protein |
| Iridoviridae | major capsid protein | protein |
| Asfarviridae | B646L (p72) | protein |
| Baculoviridae, Nudiviridae | lef-8 | protein |
| Ascoviridae | DNA polymerase | protein |
| Poxviridae | rpo147 | protein |
| Nimaviridae, Hytrosaviridae, Phycodnaviridae, Mimiviridae, Marseilleviridae, Pandoraviridae, Pithoviridae | DNA polymerase | protein |
| Mamonoviridae (medusaviruses) | DNA polymerase B (HMM `DNA_pol_B`, cut to the domain ± 100 aa) | protein |

Families using **multi-marker protein concatenation** (`region: concatenated`) are documented in [CONCAT_DESIGN.md](CONCAT_DESIGN.md); when a family appears there, the concat preset takes precedence over the single-marker fallback in the table above. The presets were re-tuned against actual NCBI annotation coverage in the 2026-05 cache audit (v1.2.30 also added `[Title]` to the protein fetch query — `[Protein Name]` is sparsely populated for several viral families, so the actual product name had to be matched against the FASTA defline instead):

- **Ascoviridae** — split off from the 7-gene baculovirus core (LEF / PIF genes 0 % coverage in Ascoviridae) to its own 3-marker preset (DNA polymerase, major capsid protein, DNA helicase P143).
- **Iridoviridae** — trimmed from 7 to 6 markers (VLTF-3 dropped — genuinely absent from NCBI Iridoviridae annotations even after the `[Title]` query fix in v1.2.30).
- **Poxviridae** — trimmed from 9 to 8 markers (single-stranded DNA-binding protein / I3L at 1 %).
- **Alloherpesviridae** — split off from the Orthoherpesviridae 7-gene set to a 6-marker subset (ssDNA-binding / ICP8 at 0 % in fish herpesviruses).
- **Malacoherpesviridae** — removed from concatenation entirely (5 species, only MCP and DNA pol annotated); falls back to single-marker MCP.
- **Phycodnaviridae** — split off from the 8-marker NCLDV hallmark fallback to a custom 4-marker preset (DNA polymerase, major capsid protein, packaging ATPase, mRNA capping enzyme); the other four hallmarks remain at 0–1 % coverage even with the `[Title]` query fix.
- **Mimiviridae, Marseilleviridae** — markers identified by HMM since 2026-10-27 (bundled set `NCLDV`, six of its markers: DNA polymerase B, RNA polymerase largest and second subunits, D5-like primase-helicase, mRNA capping enzyme, topoisomerase II), `relative_length_cutoff` 0.30. Mimiviridae genomes carry several major capsid proteins and often two VLTF-3-like proteins, so those two markers are left out (one protein per marker is chosen, the longest); Marseilleviridae uses the same six. By GenBank name (until 1.6.2) Noumeavirus, Golden Marseillevirus and many Mimiviridae species were missing. Mimiviridae runs without branch-length outlier removal (`outlier_removal.enabled: false`): most of its genomes are near-identical Megamimivirinae, so the median + 20 × MAD threshold falls to ~0.10 subs/site and real members of the divergent Klosneuvirinae were removed. Its preset leaves out Megavirus haemonchus182742 (`manual.exclude: [BK066332.1]`), a record assembled from a nematode genome project that sits alone on a long branch at the base of Megamimivirinae.
- **Asfarviridae, Pithoviridae** — markers identified by HMM since 1.6.2 (bundled set `NCLDV`: DNA polymerase B, major capsid protein, RNA polymerase largest and second subunits, D5-like primase-helicase, VLTF-3, mRNA capping enzyme, topoisomerase II; Pithoviridae without the capsid protein, which pithoviruses lack in a Pfam-detectable form), with `relative_length_cutoff` 0.30 because HMMER splits their polymerase domain at an insertion. By GenBank name, African swine fever virus matched only its DNA polymerase and p72 and pithoviruses 2 genomes, so neither family had a tree.

- **Mamonoviridae** — single-marker family-B DNA polymerase identified by HMM (bundled `DNA_pol_B`, PF00136) and cut to the domain ± 100 aa. Two of the four NCBI genomes (2026-10) split the polymerase gene into two proteins, so a name-based pick would compare different pieces. The presets listed it under the proposed name *Medusaviridae*, which NCBI never used, so until 1.6.1 the family ran as whole-genome nucleotide (four 360–560 kb genomes; MAFFT ran for hours).

The 8-marker NCLDV-hallmark fallback now applies only to Pandoraviridae (no cache audit data yet). If a stale auto-generated config file exists with incorrect settings for any of these families, the program logs a warning and suggests deleting the file to regenerate it.

#### RNA virus families (RdRp)

Whole-genome nucleotide alignment is saturated at family scale for most RNA viruses: in a survey of 117 earlier nucleotide family trees, 84 had a median tip-to-tip distance above 1.5 substitutions per site and 24 more above 0.75 (see the [phylogenetic-signal check](#5-alignment-and-tree-inference)). Their deep structure was noise. So **107 RNA-virus families** (kingdom *Orthornavirae*) default to the **RNA-dependent RNA polymerase** — the one protein all of them encode and the standard family-level marker — as a protein tree (`RDRP_FAMILIES` in `config.py`):

- **Identified by HMM**, not by name, against the bundled Pfam set `RdRp` (23 profiles: `RdRP_1`–`RdRP_5`, `Mononeg_RNA_pol`, `Bunya_RdRp`, `Arena_RNA_pol`, `Flu_PB1`, `Flavi_NS5`, `Mitovir_RNA_pol`, `RNA_replicase_B`, the reovirus and birnavirus polymerases, …). GenBank names the protein in dozens of ways (NIb, nsP4, 2a, lambda-3, L, PB1, ORF1ab, polyprotein, "183 kDa protein") and HMM identification makes that irrelevant.
- **Only the family's own kind of RdRp.** Each family's marker tokens are the profiles its own RefSeq RdRps hit (`RDRP_FAMILY_PROFILES`: e.g. `Mononeg_RNA_pol` for Rhabdoviridae; `Bunya_RdRp` and `Flu_PB1` for Phenuiviridae; the four broad palm-domain profiles `RdRP_1`–`RdRP_4` together for +ssRNA families such as Togaviridae), and `hmm.best_hit_only` makes a protein count only when one of them fits it better than any other profile in the set (profiles for other domains of the same polymerase, such as the bunyavirus L thumb-ring and PA-C-like domains, count as the family's own: `hmm.best_hit_also`). Those two L domains outside the RdRp core never select or cut a protein, even where a family's RefSeq L proteins hit them (Peribunyaviridae, Phenuiviridae, Cruliviridae): a partial L covering only such a domain would otherwise enter the alignment with no RdRp columns. Every RefSeq L protein of those families also hits `Bunya_RdRp`. This keeps out other viruses filed under the family's taxa: in October 2026 NCBI had about 20 orbivirus and bunyavirus species under *Alphavirus*, whose VP1 and L proteins fit `Orbi_VP1` / `Bunya_RdRp` best and would otherwise have joined the Togaviridae tree (whole-genome nucleotide mode took them too). Each such species is named in a WARNING in the family log, with its accessions.
- **Cut to the RdRp domain** ± 100 residues (`sequence.domain_padding`): in polyprotein viruses the longest RdRp-carrying protein is the whole polyprotein, whose structural and protease regions would otherwise dominate the alignment. The span is the run of profile hits around the best one; a hit more than 500 residues from its neighbours (a stray hit, or a second RdRp copy in a chimeric record) is left out. Metadata is kept and a GenBank COMMENT records the residues used. Proteins given by hand (`manual.include` protein accessions, `include_seq`, `include_fasta_files`) are cut the same way.
- **Proteins under 100 aa are never candidates** (`hmm.min_protein_length`): quality control drops them anyway. They are removed after download, because a length clause makes each NCBI protein search several seconds slower.
- **Hit rule:** domain E-value ≤ 1e-3 and ≥ 20 % profile coverage, no Pfam GA (`hmm.use_ga_when_available: false`). Every profile in the set is an RdRp, so the looser rule costs little; on RefSeq proteins it raised the families with ≥ 70 % of viruses covered from 78 to 107 of 138, and an audit of the extra hits found RdRps under other names except ~4 in several thousand.
- **Families stay nucleotide** when fewer than 70 % of their RefSeq viruses had an RdRp hit (divergent RdRps Pfam does not model: ormycoviruses, ambiviruses, most botourmia-, narna- and hypo-like families, Phasmaviridae, Qinviridae, …; listed with their coverage in `config.py`), and by design for deltaviruses (Kolmioviridae), satellites and viroids. Reverse-transcribing families have no RdRp; they are on the RT marker ([below](#reverse-transcribing-and-cress-dna-families-rt-rep)).
- The segment keyword of segmented families is cleared: the HMM finds the RdRp on whichever segment carries it.

To keep nucleotide analysis for a family, set in its config:

```yaml
sequence:
  type: nucleotide
  region: whole_genome
  segment: segment L      # segmented families only: the segment to use
hmm:
  enabled: false
```

A file that sets `sequence.type` other than protein, `sequence.region` other than `RdRp`, or `hmm.enabled: false` gets no part of the RdRp preset: the family runs on the template defaults under the file's own values. A file written without any of these keys (as in v1.5.0) therefore gets whole-genome nucleotide with the usual segment. A file generated by v1.6.0 holds every key, so edit all of the lines above: its `segment: null` would otherwise put all segments into one alignment. The load step logs a warning naming the type, region and segment that will run. Two half-edited combinations stop the run at load, with the lines to fix: `hmm.enabled: true` with `sequence.type: nucleotide` (the run would fetch RdRp proteins and reject them all as nucleotide), and a `sequence.segment` with the HMM protein marker (protein titles name no segment, so every record would be dropped). Outbreak trees (`additional_data_sources`) need nucleotide analysis, as Pathoplexus supplies nucleotide genomes.

#### Reverse-transcribing and CRESS-DNA families (RT, Rep)

The same approach covers the two groups without an RdRp whose whole-genome trees were saturated or weak in the survey: the reverse-transcribing viruses (*Pararnavirae*) on the **reverse transcriptase**, and the circular Rep-encoding ssDNA viruses (*Cressdnaviricota*, plus Alphasatellitidae, whose satellites encode a circovirus-like Rep) on the **Rep** protein. Hit rule, 100-aa floor, opt-out (`sequence.region: RT` / `Rep` instead of `RdRp`) and the load-time checks are as for the RdRp. Coverage of RefSeq viruses was tested the same way (2026-10-21):

- **RT** (bundled set `RT`: `RVT_1`, `RVT_2`, `RVT_thumb`, `RVT_connect`): Retroviridae, Caulimoviridae, Hepadnaviridae, Metaviridae. `RVT_1` hits every RefSeq caulimo-, hepadna- and metavirus and 64 of 89 retroviruses; the other 25 have no RT gene at all (transforming viruses that carry an onc gene in place of *pol*, endogenous-retrovirus fragments), so they are not in the Retroviridae tree. The RT is cut to its domain ± 100 residues (~310–370 aa; Gag-Pol, the hepadnavirus P protein and caulimovirus polyproteins are 800–1850 aa). The thumb and connection domains count only for the best-hit rule; `RVT_2` (Ty1/copia RT) is there so such a protein is recognised as foreign. Pseudoviridae (no RefSeq proteins) stays as it was. Hepadnaviridae was whole-genome nucleotide before.
- **Rep** (bundled set `Rep`: `Viral_Rep`, `RNA_helicase` (the SF3 helicase of Rep), `Gemini_AL1`, `Gemini_AL1_M`): 23 families — Adamaviridae, Alphasatellitidae, Amesuviridae, Anicreviridae, Bacilladnaviridae, Circoviridae, Draupnirviridae, Endolinaviridae, Gandrviridae, Geminiviridae, Genomoviridae, Geplanaviridae, Kanorauviridae, Kirkoviridae, Mahapunaviridae, Nanoviridae, Naryaviridae, Nenyaviridae, Oomyviridae, Ouroboviridae, Pecoviridae, Redondoviridae, Vilyaviridae; ≥ 94 % of RefSeq viruses in each. Rep proteins are ~270–470 aa and are kept whole. Geminivirus-like and circovirus-like Reps hit different profiles, but the profiles do not tell families apart, so every Rep profile selects and there is no best-hit rule. Smacoviridae keeps its Rep chosen by GenBank name (Pfam finds only 38 of 88 smacovirus Reps); Metaxyviridae (no RefSeq proteins) and Tolecusatellitidae (betasatellites encode no Rep) stay whole-genome nucleotide. Circoviridae used the name-based Rep before.

#### Species without the marker

A tree built on one marker (RdRp, RT, Rep, a DNA marker protein, a segment, or complete genomes) can only show species NCBI has that sequence for. Two things keep that loss small and visible:

- **Translation fallback** (HMM markers: RdRp, RT, Rep, Adenoviridae hexon, Mamonoviridae DNA polymerase): when no annotated protein of a species carries the marker, the species' nucleotide records (up to 100, complete genomes first) are translated in all six frames and searched with the same profiles and rules. This finds markers NCBI never annotated as a protein: foamy-virus and endogenous-retrovirus genomes without a *pol* CDS, metagenomic records with no protein at all. The record is named after its genome accession, keeps its metadata, and a GenBank COMMENT says it was translated. On Retroviridae (2026-10-21) it recovered the galago prosimian foamy virus, bat spumavirus 1, the *Xenopus* endogenous retrovirus Xen1 and SHIV (whose 70,000 GenBank proteins are mostly env clones). It uses the standard genetic code.
- **Report:** a species with NCBI records but none of the marker is listed in a WARNING in the family log and in `<Family>_species_without_marker.tsv` (with the species found only by translation), and counted in `summary.tsv` (`species_without_marker`, `species_marker_translated`) and in the report PDF. For Retroviridae these are mostly the acutely transforming viruses, whose genomes carry a host oncogene in place of *pol*. In concatenated mode a species whose genomes carry some of the markers, but none as many as `concatenation.min_fraction` asks for, is listed as `too few markers` (Marseilleviridae with the name-based preset of 1.6.2: Noumeavirus and Golden Marseillevirus, whose packaging ATPase and RNA polymerase subunits NCBI names differently).

#### HMM-based marker identification

For large DNA virus families whose GenBank annotation is inconsistent (synonyms, gene-symbol vs descriptive names, mis-/un-annotated or divergent orthologs), marker proteins can be identified by **HMMER profile homology** instead of protein-name queries. Enable it per family with the `hmm:` block and give each marker an `hmms:` token list:

```yaml
hmm:
  enabled: true
  database: Poxviridae          # bundled set name, a dir of .hmm files, or a single .hmm path
  default_evalue: 1.0e-5        # used only for profiles lacking a curated Pfam GA
  use_ga_when_available: true   # prefer each profile's curated GA bit-score
  relative_length_cutoff: 0.5   # min fraction of the HMM model a hit must cover
  multidomain_overlap_tolerance: 30   # residue slack between adjacent domains of an "A--B" token
  best_hit_only: false          # true: a protein counts only if one of the marker's profiles fits it best
  best_hit_also: []             # profiles that also count as the marker's own for best_hit_only
                                #   (other domains of the same protein); not tokens
  min_protein_length: null      # proteins shorter than this (aa) are never candidates
  threads: null                 # hmmscan --cpu; null → the run's thread count

# Single-marker family (DNA_FAMILIES, RDRP_FAMILIES) — token on the sequence block.
# domain_padding: null keeps the whole protein; an integer keeps only the
# span of the token's hits ± that many residues (the RdRp and RT presets use 100).
# OR'd tokens: try the 2-domain architecture first, fall back to the single domain.
sequence:
  region: hexon
  type: protein
  hmms: ["Adeno_hexon--Adeno_hexon_C", "Adeno_hexon"]

# Concatenation family — token per marker:
concatenation:
  proteins:
    - {name: "DNA polymerase", hmms: ["DNA_pol_B"], locus_tag_hint: "E9L|polB"}
    # ...
```

**Tokens.** A token is either a single profile `NAME` (`"DNA_pol_B"`) or a **multidomain architecture** `"A--B--C"` whose HMMs must hit in N-to-C order along the protein (non-overlapping, within `multidomain_overlap_tolerance`). Multiple tokens on one marker are **alternatives** (OR'd) — useful when a protein has divergent architectures. A marker that declares `hmms` is HMM-authoritative (no name fallback); markers without `hmms` fall back to name+alias matching, so mixed sets work.

**How it works.** When `hmm.enabled`, the per-marker name queries are replaced by an all-proteins-per-species fetch. Each protein is scanned once (`hmmscan`, batched and sequence-deduplicated); a hit passes when it clears the curated Pfam **GA** bit-score (or `default_evalue`) **and** covers ≥ `relative_length_cutoff` of the HMM model. Proteins are grouped by source genome and the longest token-satisfying protein is chosen per marker; `max_per_species` then bounds the number of genomes kept (RefSeq-first). Works on both the single-marker and concatenated-marker paths; the global sequence cache is bypassed in HMM mode.

**Bundled profiles.** Curated CC0 **Pfam-A** sets ship under `vfam_trees/data/hmms/<Family>/` and are selected by bare family name (`database: Poxviridae`). **Poxviridae** (8 single-domain markers), **Adenoviridae** (hexon, as the 2-domain `Adeno_hexon--Adeno_hexon_C` architecture OR'd with the single `Adeno_hexon` domain) **RdRp** (23 RNA-polymerase profiles, shared by the 107 RNA-virus families above), **RT**, **Rep**, **Mamonoviridae** (`DNA_pol_B`, cut to the domain ± 100 aa) and **NCLDV** (core genes of large DNA viruses, used by the Asfarviridae, Pithoviridae, Mimiviridae and Marseilleviridae concatenations) are bundled and enabled in their presets. (A live 761-species Adenoviridae run showed the strict 2-domain gate is highly specific but under-sensitive — `Adeno_hexon_C` misses its GA cutoff on divergent hexons and partial deposits split the domains across records — hence the single-domain fallback. Likewise a live Poxviridae run showed the generic capping-enzyme Pfam families don't match the poxvirus D1 large subunit at all, so that marker uses the N7-methyltransferase domain `mRNA_G-N7_MeTrfase` with a per-family `relative_length_cutoff` of 0.40 — empirical validation per family matters.) Point `database` at your own directory of `.hmm` files (or a single `.hmm`) for other families; the set is auto-combined and `hmmpress`-ed on first use. Rebuild the bundled sets with `scripts/build_bundled_hmms.sh`. Requires **HMMER** (`hmmscan`/`hmmpress`) on `$PATH`.

## Usage

```bash
# Generate global config template
vfam_trees init

# ... listing every setting with its default
vfam_trees init --full

# Overwrite an existing global config
vfam_trees init --force

# Check all dependencies and NCBI connectivity
vfam_trees test

# Generate per-family configs without running
vfam_trees init-configs -f families.txt

# Regenerate configs with current defaults (overwrites manual edits)
vfam_trees init-configs -f families.txt --force

# Preview per-family parameters without running the pipeline
vfam_trees run -f families.txt --dry-run

# Run pipeline (1 family at a time).  Before Snakemake starts, `run` checks
# global.yaml and that the external tools the listed families need are on
# PATH (HMMER for the RdRp and other HMM families, FastTree if a tree block
# uses it), and stops at once if not.
vfam_trees run -f families.txt

# Run with 4 parallel families, 4 threads each
vfam_trees run -f families.txt -j 4 -t 4

# Force re-run families already marked as done.  Each species' download in
# the output folder is reused only when the sequence settings (type, region,
# segment, HMM marker settings and profile content, max_per_species, excluded
# organisms) are unchanged; otherwise it is downloaded again.  Lineages of
# taxa new to the folder are fetched and added to its taxonomy cache.
vfam_trees run -f families.txt --force

# Check progress (shows current stage for in-progress families)
vfam_trees status -f families.txt

# (Re-)generate the overview PNG of all tree_100 trees
vfam_trees overview

# Cache management
vfam_trees cache clear Asfarviridae
vfam_trees cache clear --all --yes
vfam_trees cache stats
```

The `families.txt` file should contain one ICTV family name per line. Lines beginning with `#` are treated as comments and ignored:

```
# Positive-sense RNA viruses
Flaviviridae
Coronaviridae
# Negative-sense RNA viruses
Filoviridae
```

## Output

### Per-family files (`results/<Family>_<taxid>/`)

> When `manual.name` is set, `<Family>` in the directory name and in every filename below is replaced by the sanitized display name (the directory still keeps its `_<taxid>` suffix). See [Output prefix and directory naming](#output-prefix-and-directory-naming).

**Trees**

| File | Description |
|------|-------------|
| `<Family>_tree_{500,100}.nwk` | Newick tree |
| `<Family>_tree_{500,100}.xml` | PhyloXML tree (rooted, non-rerootable, `branch_length_unit="subs/site"`; `<taxonomy>`, `<sequence>`, tip `<date>`, `vipr:` metadata, rank properties; `style:font_color` when `output.phyloxml_font_colors: true`) |
| `<Family>_tree_{500,100}.html` | The PhyloXML inlined in a page that opens it in Archaeopteryx.js (loaded from jsDelivr; needs internet). Emitted unless `output.archaeopteryx_html: false` |
| `<Family>_tree_{500,100}_auspice.json` | Auspice v2 (Nextstrain) JSON — interactive divergence tree, schema-valid for [auspice.us](https://auspice.us) / Nextclade / `auspice view`; genus / subfamily / species / host / country / location / year (and **Outbreak** for `additional_data_sources` tips) colorings and filters, genus palette matching the other outputs, branch support as toggleable branch labels. Emitted unless `output.auspice_json: false` |

**Tree images** — each tree is exported in two layouts at PDF and 150 dpi PNG:

| File | Layout |
|------|--------|
| `<Family>_tree_{500,100}.{pdf,png}` | Rooted rectangular: taxonomy-annotated internal labels (genus / subgenus / subfamily / family); support values < 50% suppressed; genus/subfamily color legend; two-line figure caption (method · alignment · stats) |
| `<Family>_tree_{500,100}_ur.{pdf,png}` | Unrooted radial: equal-angle layout with leaf labels rotated outward; no internal labels, no support values; same coloring as rooted |
| `<Family>_tree_icon.png` | Square topology-only thumbnail of tree_100 (no labels, uniform branch color; size and colors configurable) |

**Sequences and metadata**

| File | Description |
|------|-------------|
| `<Family>_sequences_raw_{500,100}.fasta` | Sequences entering the MSA (post-QC, post-clustering, post-merge, post-length-outlier-filter; before post-tree branch-length outlier removal) |
| `<Family>_alignment_{500,100}.fasta` | Final alignment fed to the tree (MAFFT + trimAl when `msa_trim.enabled: true`; reflects sequences after iterative branch-length outlier removal) |
| `<Family>_metadata_{500,100}.tsv` | Per-leaf metadata |
| `<Family>_id_map.tsv` | Short ID → display name mapping (single-protein / whole-genome runs) |

**Concat-mode only** (large DNA virus families; see [CONCAT_DESIGN.md](CONCAT_DESIGN.md))

| File | Description |
|------|-------------|
| `<Family>_id_map_{500,100}.tsv` | Per-tree id-map (concat keeps source-nuc accessions as leaf IDs) |
| `<Family>_partitions_{500,100}.nex` | NEXUS charset coordinates of each marker block in the concatenated alignment |
| `<Family>_markers_{500,100}/` | Per-marker FASTAs: `<safe_marker>_raw.fasta` (unaligned) and `<safe_marker>_alignment.fasta` (post-trim when enabled) |

**Reports and logs**

| File | Description |
|------|-------------|
| `<Family>_report.pdf` | Stats table (taxid, lineage, molecule/region, species/QC counts, sequence-length stats, per-tree MSA / tree / clustering / support stats), post-QC and per-tree length histograms, support histograms, tree_100 visualization with color legend |
| `<Family>.log` | Per-family log |
| `<Family>_species_without_marker.tsv` | Species NCBI has records for but none of the analysed marker / region, or too few of the concatenated markers (not in the trees), and species whose marker was found only by translating nucleotide records ([details](#species-without-the-marker)); written only when there are any |
| `<Family>.yaml` | The settings this family ran with: built-in defaults, `global.yaml` `defaults:`, family preset and family file, merged. A header names the vfam_trees version, start time, thread count, the files merged and any keys the version did not read. Usable as the family's config to repeat the analysis |

### Cross-family files (`results/`)

| File | Description |
|------|-------------|
| `summary.tsv` | One row per family analyzed (success or skip): species counts (discovered, with sequences, without the marker, marker found only by translation), QC breakdown, sequence-length stats, clustering thresholds, MSA / trim / tree program-model-options, leaf counts, support stats, phylogenetic-signal verdict and measures, length- and branch-outlier counts and cutoffs, genus/subfamily counts, plus joined columns from `annotation_tsv` (e.g. `baltimore_class`). Rows are written by column name, so appending to a file from an earlier version first rewrites its header to the current columns |
| `status.tsv` | Lightweight one-row-per-family table: `family`, `ncbi_taxid`, `molecule_region`, `status` (`OK`; `no tree: <reason>` when the family ran to the end but neither tree had the 4 sequences it needs, e.g. `no tree: too few sequences after per-species clustering (2 < 4)`; or the reason the family was skipped), `lineage`, `baltimore_class` |
| `overview_tree_100.png` | Thumbnail grid of all tree_100 trees, shaded by viral realm (ssDNA, dsDNA, –ssRNA, +ssRNA/dsRNA, RT viruses) |
| `run_records/<start time>.yaml` | One per `vfam_trees run`; see [Run records](#run-records) |

### Run records

Each `vfam_trees run` writes `run_records/<start time>.yaml` into the output folder before Snakemake starts, and completes it when the run ends. It holds:

- the vfam_trees version, command line, working directory, family list and options (`-j`, `-t`, `--force`, log level), and the paths of `global.yaml`, the configs folder and the output folder;
- `status` (`running`, `completed`, `failed` or `interrupted`), start and finish times and Snakemake's exit code; a killed run stays `running`;
- the versions of Python, the main Python packages and the external tools (MAFFT, trimAl, IQ-TREE, FastTree, mmseqs, CD-HIT, HMMER);
- `global_config`: every global setting the run used, i.e. the built-in defaults overlaid with `global.yaml`, plus the keys it did not read. The NCBI API key is not recorded.

Family presets and family files are applied per family on top of `global_config.defaults`; what each family ran with is in its `<Family>.yaml` (above). Together the two say what a results folder was made with, without having to reconstruct an older version's defaults. A family skipped as already done keeps the `<Family>.yaml` of the run that made it.

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.
