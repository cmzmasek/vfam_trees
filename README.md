# vfam_trees

**vfam_trees** is a bioinformatics pipeline for building maximum-likelihood phylogenetic trees for viral families. For each family it discovers species via NCBI Taxonomy, downloads sequences from GenBank, applies quality filtering, clusters per species, aligns with MAFFT, infers trees with FastTree and IQ-TREE, annotates internal nodes with LCA-based taxonomy, roots taxonomy-aware, and writes Newick and PhyloXML output along with PDF/PNG figures and per-run TSV summaries.

Two trees are produced per family:

- **tree_500** — broad diversity tree (up to 500 sequences, FastTree / GTR+G or LG+G, SH-like support)
- **tree_100** — collapsed representative tree (up to 100 sequences, IQ-TREE / GTR+G or `TEST`-selected best-fit AA model; SH-aLRT support for nucleotide trees, UFBoot for protein trees)

## How it works

The per-family pipeline runs in seven stages. Each stage is configurable per family (see [Configuration](#configuration)); RefSeqs are protected at every step where sequences can be dropped.

### 1. Species discovery

The family TaxID is resolved from NCBI Taxonomy and every descendant species is enumerated. The species list defines the universe of organisms eligible for download. An optional `manual.include_species` list in the per-family config restricts this universe to a named subset — only species whose name or NCBI taxid appears in the list proceed to download; others are skipped entirely.

### 2. Sequence download

GenBank is queried per species with RefSeq priority. The query targets the configured molecule (`whole_genome` or a marker name) and, for segmented viruses, the configured segment keyword (records lacking it in their title are dropped). RefSeqs are uncapped; non-RefSeq records are limited by `download.max_per_species`.

An optional shared cache (configurable TTL, per-entry locking, negative-result caching) avoids re-downloading the same species across runs and across parallel family jobs.

### 3. Quality filtering

Each species' sequences are filtered in this order:

0. **Manual overrides** (per-family `manual.include` / `manual.include_fasta` / `manual.exclude` / `manual.include_species`) — accessions in `manual.exclude` are dropped immediately; accessions in `manual.include` bypass every QC step below and are protected from removal during clustering, proportional merge, and length-outlier filtering. Match is exact, version included (e.g. `NC_002617.1`). `manual.include_fasta` accepts pasted sequences (records not yet in GenBank) as a list of `{id, organism, sequence}` mappings; entries are injected after fetch, receive the same bypass + downstream protection as `manual.include`, and their ids must not collide with any fetched accession (not supported in concat mode). `manual.include_species` restricts the set of species that proceed to download (see [Species discovery](#1-species-discovery)).
1. **Organism exclusion** — case-insensitive substring match against `ORGANISM`, `SOURCE`, and `DEFINITION` (joined with newlines so terms cannot straddle field boundaries). Defaults: `synthetic construct`, `metagenome`, `MAG:`, `uncultured`, `unverified`, `vector`, `recombinant`, `patent`.
2. **Ambiguity** — `max_ambiguous` cap on the fraction of `N` / `X` / IUPAC degenerate characters.
3. **Minimum length** — `min_length: null` auto-sets the threshold to 50% of the per-species median, with relaxation fallback to 40% then 30% if too few sequences pass; in all cases a hard floor of 200 bp / 100 aa applies.
4. **Length-outlier filter** (post-merge, pre-MSA) — two-sided keep window that is the **union** of (a) a MAD-on-log-lengths window and (b) a hard floor `[min_lo_mult, max_hi_mult] × median`. See [Length-outlier behaviour](#length-outlier-behaviour) for the exact form. RefSeqs and `manual.include` records flagged by the filter are kept with a warning.

### 4. Selection and clustering

- **RefSeq absorption** — non-RefSeq sequences ≥ `refseq_absorption.threshold` (default 0.99) identical to a RefSeq within the same species are absorbed into the RefSeq, suppressing redundant near-zero-branch cherries. RefSeqs themselves are never removed.
- **Adaptive per-species clustering** — MMseqs2 with binary search for an identity threshold within `[clustering.threshold_min, threshold_max]` that yields ≤ `max_reps_500` (or `max_reps_100`) representatives.
- **Proportional cross-species merge** — when the species count exceeds the target tree size, species with at least one RefSeq are kept first; remaining slots fill by rep count. Species dropped at the cap (and how many of them carried a RefSeq) are recorded in `summary.tsv` and a WARNING is logged.

### 5. Alignment and tree inference

- **MAFFT** with separate options for nucleotide vs. protein, and for tree_500 vs. tree_100 (the latter typically slower / more accurate).
- **trimAl** column trimming (`-automated1` by default) between MAFFT and tree inference; pre-trim length, tool, and options are recorded per tree.
- **FastTree** for tree_500 (configurable `model_nuc` / `model_aa`; `+G` / `+GAMMA` enables discrete-gamma rate variation).
- **IQ-TREE** for tree_100, with per-sequence-type options: `--fast` for nucleotides (SH-aLRT support, auto-added by the wrapper) and `-B 1000` for protein (UFBoot ultrafast bootstrap, more robust on divergent protein families). Setting `model_aa: TEST` triggers ModelFinder; the chosen model (e.g. `LG+I+G4`) is parsed from the IQ-TREE log and surfaces in `summary.tsv`, the per-family PDF, and the PhyloXML provenance.
- **Branch-support measure** is picked per tree (`SH_like` / `SH_aLRT` / `UFBoot`) and recorded uniformly: `tree{500,100}_support_type` plus generic `support_{min,q1,median,q3,max,iqr}` columns; the PhyloXML `<confidence type="…">` attribute mirrors the same label.
- **Iterative branch-length outlier removal** — between iterations, leaves with terminal branch length exceeding `median + factor × MAD` are removed and MSA + tree are re-run (up to `max_iterations`; only when ≥ `min_seqs` remain). RefSeqs are protected: a flagged RefSeq stays with a warning. Detailed per-outlier log lines include length, ratio to median, and threshold.
- **Multi-marker protein concatenation** — large DNA virus families with curated marker presets concatenate per-marker MAFFT + trimAl alignments into a single matrix; tree_100 uses partitioned IQ-TREE (`-p partitions.nex -m MFP`) so each marker gets its own ModelFinder pick. Per-marker fetches are bucketed by source nucleotide accession (Policy A); a length filter (`concatenation.source_nuc_min_length_frac`, default 0.3) drops source-nuc accessions shorter than that fraction of the longest parent so partial single-gene submissions don't crowd out genome-scale proteins under `download.max_per_species` (default 3000 for concat families). See [CONCAT_DESIGN.md](CONCAT_DESIGN.md) for design details.

### 6. Annotation and rooting

- **LCA-based internal-node annotation** using NCBI ranked lineages, with a configurable lineage-depth filter (`taxonomy.lca_min_rank`, default `species`) so leaves whose lineage is too shallow are excluded from the LCA vote.
- **Genus inference** for taxa lacking a formal `genus` rank in NCBI lineage: `none`, `suffix` (single-word taxa ending in `virus` treated as genus, per ICTV convention), or `deepest` (default — suffix first, then deepest rank above species).
- **Taxonomy-guided rooting** via LCA specificity scoring, with MAD and midpoint fallbacks. The PhyloXML `<phylogeny>` is emitted with `rooted="true" rerootable="false"` so downstream viewers cannot re-root the carefully rooted output.
- **Genus / subfamily coloring** in HLS color space — one hue band per subfamily, lightness varies across genera within a subfamily; colors are applied to PDF/PNG figures and to PhyloXML (`style:font_color`).

### 7. Output

Per-family: Newick, PhyloXML (with `<taxonomy>`, `<sequence>`, `vipr:` metadata, and rank properties), rooted-rectangular and unrooted-radial PDF + PNG tree images, a topology-only icon PNG, sequence and metadata FASTAs/TSVs, a per-family PDF report, and a per-family log.

Cross-family: a row-per-family `summary.tsv` with full statistics, a lightweight `status.tsv` (success/skip), and an `overview_tree_100.png` thumbnail grid shaded by viral realm. See [Output](#output) for the full file inventory.

Per-family directories are named `<Family>_<taxid>` (e.g. `Asfarviridae_137992`). Failures and skips at any stage produce a row in `status.tsv` with the skip reason; the per-family work directory is preserved for inspection.

## Other capabilities

- **Pre-configured family presets** — 28 segmented RNA virus families (segment keywords) and 26 DNA virus families (sequence type, region/marker, and concat marker sets for the 16 families that use concatenation) ship with curated defaults.
- **Checkpointing** — MSA and tree steps store a content-hashed sidecar of inputs + tool / model / options, so any change auto-invalidates the cache. Each iterative outlier-removal pass gets its own hash; partially-completed runs resume in the right place.
- **Family-annotation TSV** (optional) — joins extra per-family columns (e.g. `baltimore_class`) into `summary.tsv` and `status.tsv`. Missing file or family is silent.
- **Stage tracking** — `vfam_trees status` reports the current processing stage for in-progress families (downloading/QC, MSA, tree inference, annotating).
- **Dry-run mode** — `vfam_trees run --dry-run` previews per-family parameters without executing.
- **Parallelism** — `-j N` runs N families concurrently; `-t T` sets threads per family job. The download cache and external tools are concurrency-safe.
- **Validation and warnings** — partial NCBI batches, unrecognized per-family YAML keys, configs that override a recommended DNA-family setting, and MAFFT / tree output sanity checks all surface as log warnings.

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
| `FastTree` | Rapid ML tree inference (tree_500) |
| `iqtree2` | ML tree inference (tree_100) |
| `mmseqs` | Sequence clustering |

All tools must be available on `$PATH`. Installation via conda is recommended:

```bash
conda install -c bioconda mafft fasttree iqtree mmseqs2 trimal
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
vfam_trees init
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

Generate a template with `vfam_trees init`. Then edit it to set your NCBI credentials:

```yaml
ncbi:
  email: your.email@example.com     # REQUIRED
  api_key: your_ncbi_api_key        # optional but recommended (10 req/s vs 3 req/s)
```

An NCBI API key can be obtained for free at https://www.ncbi.nlm.nih.gov/account/.

The `defaults:` section of `global.yaml` overrides the built-in defaults for all families. Per-family configs further override individual parameters.

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
annotation_tsv: virus_families_annotation.tsv
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
    - uncultured
    - "MAG:"                    # metagenome-assembled genomes (NCBI DEFINITION prefix)
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
  options_nuc: "--6merpair --retree 2"   # used for nucleotide sequences
  options_aa: "--auto"                   # used for amino acid sequences (MAFFT auto-selects strategy)

msa_100:
  tool: mafft
  options_nuc: "--retree 3"                       # used for nucleotide sequences
  options_aa: "--maxiterate 1000 --localpair"     # used for amino acid sequences (MAFFT L-INS-i; high accuracy, slower)

msa_trim:
  enabled: true                  # drop poorly-aligned columns before tree inference
  tool: trimal
  options: "-automated1"         # adaptive; works for both nucleotide and protein

tree_500:
  tool: fasttree
  options_nuc: ""             # used for nucleotide sequences
  options_aa: ""              # used for amino acid sequences (e.g. marker-gene families)
  model_nuc: GTR+G
  model_aa: LG+G              # LG+G used for amino acid sequences

tree_100:
  tool: iqtree
  options_nuc: "--fast"        # nucleotide: SH-aLRT support (auto-added by wrapper)
  options_aa: "-B 1000"        # protein: UFBoot ultrafast bootstrap (stronger support
                               # for divergent protein families; --fast is incompatible with -B)
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

manual:                         # curator overrides on per-family record selection
  include: []                   # exact accessions (with version, e.g. "NC_002617.1") that
                                # bypass all QC (length, ambiguity, organism exclusion) and
                                # are protected through clustering, proportional merge, and
                                # length-outlier filtering — stronger than RefSeq at QC,
                                # equal to RefSeq downstream
  include_fasta: []             # pasted sequences not yet in GenBank — each entry is a
                                # {id, organism, sequence} mapping; injected after fetch,
                                # bypass QC, and receive the same downstream protection as
                                # manual.include records (not supported in concat mode)
  exclude: []                   # exact accessions to drop immediately after fetch, before QC
  include_species: []           # restrict the pipeline to a subset of species — entries are
                                # species names (case-insensitive) or NCBI taxids (integer or
                                # digit string); species not listed are skipped before download;
                                # empty or absent means use the full discovered species list
```

The `coloring` and `taxonomy` keys can also be set globally in the `defaults:` section of `global.yaml` — per-family configs inherit them automatically.

Example — adding a pasted reference and forcing a curator-picked accession in:

```yaml
manual:
  include:
    - NC_002617.1               # force-keep this RefSeq even if QC would drop it
  include_fasta:
    - id: MyLab_Isolate_2026
      organism: Foo virus
      sequence: |
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    - id: PreprintSeq_42
      organism: Bar virus
      sequence: ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  exclude:
    - KF234567.1                # drop a known-bad isolate before QC
```

Notes on `include_fasta`:

- `id` is used throughout the pipeline (FASTA headers, tree leaves, PhyloXML) and must not collide with any accession returned by NCBI for the family or with anything in `include` / `exclude` — collisions raise a hard error.
- `organism` controls the species bucket the entry joins (matching a fetched species name appends to that bucket; a novel name creates a new one) and shows up in leaf labels.
- The sequence may be wrapped across lines — whitespace is stripped and the sequence is uppercased.
- Pasted leaves render gray (no genus) and do not participate in LCA voting, since only `id` and `organism` are required.
- Not supported when `sequence.region == "concatenated"` — pasted sequences cannot be split into per-marker proteins.

Example — restricting a family run to a specific subset of species (mixing names and taxids):

```yaml
manual:
  include_species:
    - Hantaan orthohantavirus       # matched case-insensitively against NCBI species name
    - Seoul orthohantavirus
    - 11595                          # NCBI taxid (integer) — equivalent to the species name above
    - "1980519"                      # taxid as a quoted digit string — also accepted
```

Unmatched entries produce a WARNING in the log (no hard error — it may indicate a typo). If `include_species` matches nothing in the discovered species list the family is skipped with an explanatory status entry.

#### Length-outlier behaviour

The pre-MSA length-outlier filter takes the **union** of two windows around the median sequence length:

- **MAD-on-log-lengths**: `exp(median(log L) ± k · σ_log)` with `σ_log = 1.4826 · MAD(log L)`. Adapts to each family's natural spread; `k=0` disables.
- **Hard floor**: `[min_lo_mult, max_hi_mult] × median`. Guarantees a minimum keep window even when MAD is degenerate or the family is very tight; either knob set to `0` disables that side.

Behaviour across distribution shapes (median = 1000, defaults `k=5.0`, floor `[0.20×, 5.0×]`):

| Distribution | Effective keep window | Notes |
|---|---|---|
| Tight family + 0.40× truncation | `[0.20×, 5.0×]` | Floor kicks in; moderate truncations kept |
| Variable family (~2× spread) | `[0.19×, 5.23×]` | MAD widens beyond the floor |
| Bimodal (≥50% tied at median) | `[0.20×, 5.0×]` | MAD = 0; floor alone applies |
| Tight bulk + 30× extreme outlier | `[0.20×, 5.0×]` | Floor catches the outlier |
| Tight bulk + 0.03× extreme outlier | `[0.20×, 5.0×]` | Floor catches the outlier |

Filter outcome (kept / dropped, the median, and the resolved keep window) is logged at INFO once per tree and written to `summary.tsv` (`tree{500,100}_n_length_outliers_{long,short}`, `length_filter_median`, `length_filter_lo_cutoff`, `length_filter_hi_cutoff`).

#### DNA virus families

Known DNA virus families are automatically configured with curated ICTV-aligned markers. Small ssDNA / dsDNA families use a single diagnostic protein (when one is well established), medium–large dsDNA families use a family-specific structural or replication protein, and nucleocytoplasmic large DNA viruses (NCLDVs) share DNA polymerase (family B) as a universal marker:

| Family group | Marker gene | Sequence type |
|---|---|---|
| Circoviridae, Smacoviridae | Rep | protein |
| Anelloviridae | ORF1 | protein |
| Parvoviridae | NS1 | protein |
| Polyomaviridae | large T antigen | protein |
| Papillomaviridae | L1 | protein |
| Hepadnaviridae | whole genome | nucleotide |
| Adenoviridae | hexon | protein |
| Orthoherpesviridae, Alloherpesviridae, Herpesviridae | DNA polymerase | protein |
| Malacoherpesviridae | major capsid protein | protein |
| Iridoviridae | major capsid protein | protein |
| Asfarviridae | B646L (p72) | protein |
| Baculoviridae, Nudiviridae | lef-8 | protein |
| Ascoviridae | DNA polymerase | protein |
| Poxviridae | rpo147 | protein |
| Nimaviridae, Hytrosaviridae, Phycodnaviridae, Mimiviridae, Marseilleviridae, Pandoraviridae, Pithoviridae, Medusaviridae | DNA polymerase | protein |

Families using **multi-marker protein concatenation** (`region: concatenated`) are documented in [CONCAT_DESIGN.md](CONCAT_DESIGN.md); when a family appears there, the concat preset takes precedence over the single-marker fallback in the table above. The presets were re-tuned against actual NCBI annotation coverage in the 2026-05 cache audit (v1.2.30 also added `[Title]` to the protein fetch query — `[Protein Name]` is sparsely populated for several viral families, so the actual product name had to be matched against the FASTA defline instead):

- **Ascoviridae** — split off from the 7-gene baculovirus core (LEF / PIF genes 0 % coverage in Ascoviridae) to its own 3-marker preset (DNA polymerase, major capsid protein, DNA helicase P143).
- **Iridoviridae** — trimmed from 7 to 6 markers (VLTF-3 dropped — genuinely absent from NCBI Iridoviridae annotations even after the `[Title]` query fix in v1.2.30).
- **Poxviridae** — trimmed from 9 to 8 markers (single-stranded DNA-binding protein / I3L at 1 %).
- **Alloherpesviridae** — split off from the Orthoherpesviridae 7-gene set to a 6-marker subset (ssDNA-binding / ICP8 at 0 % in fish herpesviruses).
- **Malacoherpesviridae** — removed from concatenation entirely (5 species, only MCP and DNA pol annotated); falls back to single-marker MCP.
- **Phycodnaviridae** — split off from the 8-marker NCLDV hallmark fallback to a custom 4-marker preset (DNA polymerase, major capsid protein, packaging ATPase, mRNA capping enzyme); the other four hallmarks remain at 0–1 % coverage even with the `[Title]` query fix.
- **Mimiviridae, Marseilleviridae, Pithoviridae** — split off from the NCLDV fallback to 6-marker presets (per-family drop list).

The 8-marker NCLDV-hallmark fallback now applies only to Pandoraviridae and Medusaviridae (no cache audit data yet). If a stale auto-generated config file exists with incorrect settings for any of these families, the program logs a warning and suggests deleting the file to regenerate it.

## Usage

```bash
# Generate global config template
vfam_trees init

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

# Run pipeline (1 family at a time)
vfam_trees run -f families.txt

# Run with 4 parallel families, 4 threads each
vfam_trees run -f families.txt -j 4 -t 4

# Force re-run families already marked as done
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

**Trees**

| File | Description |
|------|-------------|
| `<Family>_tree_{500,100}.nwk` | Newick tree |
| `<Family>_tree_{500,100}.xml` | PhyloXML tree (rooted, non-rerootable; `<taxonomy>`, `<sequence>`, `vipr:` metadata, rank properties, `style:font_color`) |

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

### Cross-family files (`results/`)

| File | Description |
|------|-------------|
| `summary.tsv` | One row per family analyzed (success or skip): species counts, QC breakdown, sequence-length stats, clustering thresholds, MSA / trim / tree program-model-options, leaf counts, support stats, length- and branch-outlier counts and cutoffs, genus/subfamily counts, plus joined columns from `annotation_tsv` (e.g. `baltimore_class`) |
| `status.tsv` | Lightweight one-row-per-family table: `family`, `ncbi_taxid`, `molecule_region`, `status` (`OK` or skip reason), `lineage`, `baltimore_class` |
| `overview_tree_100.png` | Thumbnail grid of all tree_100 trees, shaded by viral realm (ssDNA, dsDNA, –ssRNA, +ssRNA/dsRNA, RT viruses) |

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.
