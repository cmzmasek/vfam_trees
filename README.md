# vfam_trees

**vfam_trees** is a bioinformatics pipeline for building maximum-likelihood phylogenetic trees for viral families. For each family it automatically discovers species via NCBI Taxonomy, downloads sequences from GenBank, applies quality filtering, clusters sequences per species, aligns them with MAFFT, infers trees with FastTree and IQ-TREE, annotates internal nodes with LCA-based taxonomy, and writes output in Newick and PhyloXML formats.

Two trees are produced per family:

- **tree_500** — broad diversity tree (up to 500 sequences, FastTree / GTR+G or LG+G, SH-like support)
- **tree_100** — collapsed representative tree (up to 100 sequences, IQ-TREE / GTR+G or TEST; SH-aLRT support for nucleotide trees, UFBoot (`-B 1000`) support for protein trees)

## Features

- Automatic species discovery from NCBI Taxonomy
- Per-species sequence download with RefSeq priority; RefSeq records are also preferred during cross-species proportional subsampling so reference sequences remain in the final tree set
- **Smart sequence type selection**: large DNA virus families automatically use protein marker genes (DNA polymerase, major capsid protein, hexon, etc.) instead of whole-genome nucleotide sequences; small DNA virus families use whole-genome nucleotide sequences
- Adaptive quality filtering: `min_length = null` auto-sets the threshold to 50% of the per-species median (with fallback to 40% and 30% if too few sequences pass), plus a hard floor of 200 bp / 100 aa
- Adaptive per-species clustering (MMseqs2) with binary search for optimal identity threshold
- Proportional cross-species sampling to fill target tree sizes
- Minimum sequence checks at multiple stages (post-QC, post-merge, post-outlier-removal); families and individual tree targets are skipped gracefully when too few sequences remain
- Length outlier removal before alignment: two-sided, drops sequences that are both much longer (`> hi_mult × median`, default 3×) **and** much shorter (`< lo_mult × median`, default ~1/3) than the median of the selected set; configurable per family via the `length_outlier:` block; counts (`n_length_outliers_long` / `n_length_outliers_short`) are recorded in `summary.tsv` per tree
- Iterative post-tree branch-length outlier removal: after each tree, leaves with branch length exceeding `median + factor × MAD` (Median Absolute Deviation — robust to skewed distributions) are removed and MSA+tree is re-run (up to max_iterations; removal only proceeds when at least min_seqs sequences remain; configurable per family, on by default); detailed per-outlier log messages include branch length, ratio to median, and threshold
- Separate MSA options for nucleotide vs. amino acid sequences (`options_nuc` / `options_aa`) in both msa_500 and msa_100 sections; IQ-TREE `TEST` model selects the best-fit substitution model automatically for amino acid tree_100 runs — the chosen model (e.g. `LG+I+G4`, `Q.yeast+F+I+G4`) is parsed from the IQ-TREE log and reported in `summary.tsv`, the per-family PDF, and the PhyloXML description in place of `TEST`
- FastTree respects `model_nuc` / `model_aa` from the config: `GTR` / `JC` for nucleotides, `LG` / `WAG` / `JTT` for amino acids (unsupported models fall back with a warning); the `+G` / `+GAMMA` suffix enables discrete-gamma rate variation
- MAFFT multiple sequence alignment (separate options for tree_500 and tree_100, and for nucleotide vs. protein)
- Alignment column trimming with **trimAl** (`-automated1` by default, on by default) between MAFFT and tree inference — drops poorly-aligned / ambiguous columns to improve signal-to-noise, applied uniformly to both nucleotide and protein alignments; pre-trim length, trim tool, and trim options are recorded per tree in `summary.tsv`; disable or retune per family via the `msa_trim:` config block
- FastTree (tree_500) and IQ-TREE (tree_100) tree inference; tree_100 has per-sequence-type options — `options_nuc: "--fast"` (SH-aLRT support, auto-added by the wrapper) and `options_aa: "-B 1000"` (UFBoot ultrafast bootstrap, more robust on divergent protein families)
- Branch-support measure is picked automatically per tree and recorded as `tree{500,100}_support_type` in `summary.tsv` (`SH_like` / `SH_aLRT` / `UFBoot`); support-value stats are reported under generic `tree{500,100}_support_{min,q1,median,q3,max,iqr}` columns so a single schema covers all three measures; the PhyloXML `<confidence type="…">` attribute mirrors the same label
- Taxonomy-guided tree rooting using LCA specificity scoring, with MAD and midpoint fallbacks
- LCA-based internal node annotation using NCBI ranked lineages
- PhyloXML output with `<confidence type="SH_like|SH_aLRT|UFBoot">` (type chosen per tree, matches the actual support measure computed), `<taxonomy>`, `<sequence>`, and `vipr:` metadata properties; leaf node labels are colored by genus using a `style:font_color` property
- **Genus/subfamily leaf coloring**: leaves are colored in HLS color space — one hue band per subfamily, lightness varies across genera within a subfamily; colors are applied in PDF/PNG tree images, standalone tree images, and PhyloXML output; a structured legend (grouped by subfamily) is included in all tree figures
- Segment keyword validation: for segmented RNA families, records not containing the expected segment keyword in their title are excluded. The segment query accepts any of "complete sequence", "complete genome", or "complete cds" in the record title, so per-segment CDS records are not missed
- Checkpointing: MSA and tree steps are skipped only if their inputs still match — checkpoint sidecars store a content hash of the sequence set, MSA output, and relevant config (tool / model / options), so changing any of them automatically invalidates the cache and forces a rerun. Each iterative outlier-removal pass gets its own hash, so resuming a partially-completed run picks up in the right place
- Validation of MAFFT and tree output files before continuing
- Warning when NCBI returns a partial batch
- Warning when a per-family YAML config contains unrecognized keys
- Warning when a config file overrides a recommended DNA-family setting (e.g. stale auto-generated config with `region: whole_genome` for a large DNA virus family)
- Per-family PDF report containing:
  - Statistics table (NCBI taxid, lineage, molecule/region, species counts, QC breakdown, post-QC sequence length stats, and per-tree: sequence type, MSA tool/options, tree program/model/options, leaf count, sequence length stats, MSA length/gap%, clustering thresholds, SH support stats)
  - Post-QC sequence length histogram
  - Per-tree sequence length histograms for tree_500 and tree_100 (sequences actually used to build each tree)
  - SH support value histograms for both trees
  - tree_100 visualization with genus/subfamily color legend
- Standalone PDF and PNG tree images for both tree_100 and tree_500, with genus/subfamily color legend
- **Overview PNG** (`overview_tree_100.png`): thumbnail grid of all tree_100 trees across all processed families, automatically generated at the end of `vfam_trees run`; thumbnails are shaded by viral realm (ssDNA, dsDNA, –ssRNA, +ssRNA/dsRNA, RT viruses) using NCBI lineage data; regenerate at any time with `vfam_trees overview`
- Output directories named `<Family>_<taxid>` (e.g. `Asfarviridae_137992`)
- Pre-configured support for 35+ segmented RNA virus families and 27 DNA virus families
- Per-run summary TSV with SH support statistics, MSA statistics, QC breakdown, clustering thresholds, outlier removal counts, and genus/subfamily diversity counts; skipped families are always included
- Optional shared sequence download cache keyed by query parameters, with configurable TTL and per-entry lock files for safe parallel use
- **Pipeline stage tracking**: `vfam_trees status` reports the current processing stage for in-progress families (downloading/QC, MSA, tree inference, annotating) in addition to done/pending/skipped
- **Dry-run mode**: `vfam_trees run --dry-run` previews per-family configuration parameters (sequence type, region, tree tools and models) without executing the pipeline

## Dependencies

### Python packages

```
biopython  >= 1.81
click      >= 8.1
pyyaml     >= 6.0
snakemake  >= 7.0
requests   >= 2.31
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

### 1. Global config (`config/global.yaml`)

Generate a template with:

```bash
vfam_trees init
```

Then edit it to set your NCBI credentials:

```yaml
ncbi:
  email: your.email@example.com     # REQUIRED
  api_key: your_ncbi_api_key        # optional but recommended (10 req/s vs 3 req/s)
```

An NCBI API key can be obtained for free at https://www.ncbi.nlm.nih.gov/account/

The `defaults:` section in `global.yaml` overrides the built-in defaults for all families. Per-family configs can further override individual parameters.

#### Sequence download cache

To avoid re-downloading the same species across runs or across families, enable the global cache in `global.yaml`:

```yaml
cache:
  dir: ~/.vfam_cache    # shared across all runs on this machine (or a lab filesystem)
  ttl_days: 90          # re-download after 90 days; null = never expire
```

Cache entries are keyed by `(taxid, db, region, segment, max_per_species)` so changing any query parameter automatically triggers a fresh download. Parallel family jobs (`-j N`) coordinate via per-entry lock files so the same species is never downloaded twice concurrently.

To clear the cache for a specific family (e.g. after a query fix):

```bash
vfam_trees cache clear Asfarviridae
vfam_trees cache clear --all          # wipe entire cache
vfam_trees cache stats                # show entry count and size
```

### 2. Per-family configs (`configs/<Family>.yaml`)

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
  exclude_organisms:
    - synthetic construct
    - metagenome
    - uncultured

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
  options_nuc: "--6merpair --retree 1"   # used for nucleotide sequences
  options_aa: "--6merpair --retree 1"   # used for amino acid sequences

msa_100:
  tool: mafft
  options_nuc: "--retree 2"             # used for nucleotide sequences
  options_aa: "--auto"                  # used for amino acid sequences (MAFFT auto-selects strategy)

msa_trim:
  enabled: true                  # drop poorly-aligned columns before tree inference
  tool: trimal
  options: "-automated1"         # adaptive; works for both nucleotide and protein

tree_500:
  tool: fasttree
  options: ""
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
  hi_mult: 3.0                  # drop seqs longer than hi_mult × median (0 disables)
  lo_mult: 0.333                # drop seqs shorter than lo_mult × median (0 disables)

outlier_removal:
  enabled: true                 # iterative post-tree branch-length outlier removal
  factor: 20.0                  # threshold = median + factor × MAD (Median Absolute Deviation)
  max_iterations: 3             # maximum MSA+tree iterations
  min_seqs: 40                  # only remove outliers when ≥ min_seqs sequences remain after removal
```

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
| Orthoherpesviridae, Alloherpesviridae, Malacoherpesviridae, Herpesviridae | DNA polymerase | protein |
| Iridoviridae | major capsid protein | protein |
| Asfarviridae | B646L (p72) | protein |
| Baculoviridae, Nudiviridae, Ascoviridae | lef-8 | protein |
| Poxviridae | rpo147 | protein |
| Nimaviridae, Hytrosaviridae, Phycodnaviridae, Mimiviridae, Marseilleviridae, Pandoraviridae, Pithoviridae, Medusaviridae | DNA polymerase | protein |

If a stale auto-generated config file exists with incorrect settings for these families, the program will log a warning and suggest deleting the file to regenerate it.

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

For each family, results are written to `results/<Family>_<taxid>/` (e.g. `results/Asfarviridae_137992/`):

| File | Description |
|------|-------------|
| `<Family>_tree_500.xml` | PhyloXML tree (broad, up to 500 sequences) |
| `<Family>_tree_100.xml` | PhyloXML tree (collapsed, up to 100 sequences) |
| `<Family>_tree_500.nwk` | Newick tree (broad) |
| `<Family>_tree_100.nwk` | Newick tree (collapsed) |
| `<Family>_tree_500.pdf` | Standalone PDF tree image with genus/subfamily color legend |
| `<Family>_tree_500.png` | Standalone PNG tree image (150 dpi) with genus/subfamily color legend |
| `<Family>_tree_100.pdf` | Standalone PDF tree image with genus/subfamily color legend |
| `<Family>_tree_100.png` | Standalone PNG tree image (150 dpi) with genus/subfamily color legend |
| `<Family>_alignment_500.fasta` | Final alignment fed to tree_500 (MAFFT + trimAl when `msa_trim.enabled: true`; reflects sequences after iterative outlier removal) |
| `<Family>_alignment_100.fasta` | Final alignment fed to tree_100 (MAFFT + trimAl when `msa_trim.enabled: true`; reflects sequences after iterative outlier removal) |
| `<Family>_sequences_raw_500.fasta` | Sequences entering the MSA (after QC, clustering, and proportional merge; before post-tree outlier removal) |
| `<Family>_sequences_raw_100.fasta` | Sequences entering the MSA (after QC, clustering, and proportional merge; before post-tree outlier removal) |
| `<Family>_metadata_500.tsv` | Sequence metadata (broad) |
| `<Family>_metadata_100.tsv` | Sequence metadata (collapsed) |
| `<Family>_id_map.tsv` | Short ID → display name mapping |
| `<Family>_report.pdf` | Per-family PDF report: stats table, post-QC length histogram, per-tree length histograms (tree_500 and tree_100), SH support histograms, tree_100 visualization with genus/subfamily color legend |
| `<Family>.log` | Per-family log |

At the cross-family level, in `results/`:

| File | Description |
|------|-------------|
| `summary.tsv` | One row per family: species counts, QC exclusion breakdown, post-QC sequence length stats, clustering thresholds, MSA tool/options, tree program/model/options, per-tree leaf counts, SH support stats, outlier removal counts, genus/subfamily diversity counts |
| `overview_tree_100.png` | Thumbnail grid of all tree_100 trees; thumbnails shaded by viral realm (ssDNA, dsDNA, –ssRNA, +ssRNA/dsRNA, RT viruses) |

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.
