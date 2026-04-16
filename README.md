# vfam_trees

**vfam_trees** is a bioinformatics pipeline for building maximum-likelihood phylogenetic trees for viral families. For each family it automatically discovers species via NCBI Taxonomy, downloads sequences from GenBank, applies quality filtering, clusters sequences per species, aligns them with MAFFT, infers trees with FastTree and IQ-TREE, annotates internal nodes with LCA-based taxonomy, and writes output in Newick and PhyloXML formats.

Two trees are produced per family:

- **tree_500** — broad diversity tree (up to 500 sequences, FastTree / GTR+G, SH-like support)
- **tree_100** — collapsed representative tree (up to 100 sequences, IQ-TREE / GTR+G, SH-aLRT support)

## Features

- Automatic species discovery from NCBI Taxonomy
- Per-species sequence download with RefSeq priority
- Adaptive quality filtering: `min_length = null` auto-sets the threshold to 50% of the per-species median (with fallback to 40% and 30% if too few sequences pass), plus a hard floor of 200 bp / 100 aa
- Adaptive per-species clustering (MMseqs2) with binary search for optimal identity threshold
- Proportional cross-species sampling to fill target tree sizes
- Length outlier removal before alignment (sequences >3× median length excluded)
- MAFFT multiple sequence alignment (separate options for tree_500 and tree_100)
- FastTree (tree_500) and IQ-TREE `--fast` (tree_100) tree inference
- SH-like support values (FastTree) and SH-aLRT support values (IQ-TREE); stored in PhyloXML `<confidence>` elements and reported in the summary TSV
- Taxonomy-guided tree rooting using LCA specificity scoring, with MAD and midpoint fallbacks
- LCA-based internal node annotation using NCBI ranked lineages
- PhyloXML output with `<confidence type="SH_like|SH_aLRT">`, `<taxonomy>`, and `vipr:` metadata properties
- Chimera / misannotation detection: warns when a terminal branch length exceeds 10× the median
- Segment keyword validation: for segmented RNA families, records not containing the expected segment keyword in their title are excluded
- Checkpointing: MSA and tree steps are skipped if their outputs already exist (resume after interruption)
- Validation of MAFFT and tree output files before continuing
- Warning when NCBI returns a partial batch
- Warning when a per-family YAML config contains unrecognized keys
- Per-family PDF report (statistics table, sequence length histogram, SH support histograms, tree_100 visualization)
- Pre-configured support for 35+ segmented RNA virus families and 19 DNA virus families
- Per-run summary TSV with SH support statistics, MSA statistics, QC breakdown, and clustering thresholds; skipped families are always included

## Dependencies

### Python packages

```
biopython  >= 1.81
click      >= 8.1
pyyaml     >= 6.0
snakemake  >= 7.0
requests   >= 2.31
matplotlib >= 3.9    # PDF report; requires NumPy 2.x compatible build
```

### External tools

| Tool | Purpose |
|------|---------|
| `mafft` | Multiple sequence alignment |
| `FastTree` | Rapid ML tree inference (tree_500) |
| `iqtree2` | ML tree inference (tree_100) |
| `mmseqs` | Sequence clustering |

All tools must be available on `$PATH`. Installation via conda is recommended:

```bash
conda install -c bioconda mafft fasttree iqtree mmseqs2
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

# 4. Run the pipeline
vfam_trees run -f families.txt -j 4 -t 4

# 5. Check progress
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

### 2. Per-family configs (`configs/<Family>.yaml`)

Per-family configs are auto-generated if missing. Generate them in advance to review and tune parameters before running:

```bash
vfam_trees init-configs -f families.txt
```

Key parameters:

```yaml
download:
  max_per_species: 200          # cap on non-RefSeq sequences per species

sequence:
  type: nucleotide              # nucleotide or protein
  region: whole_genome          # whole_genome, or a gene name (e.g. "hexon", "DNA polymerase")
  segment: null                 # segment keyword for segmented viruses (e.g. "segment L")
                                # auto-set for known segmented families

quality:
  min_length: null              # null = auto (50% of per-species median, with floor of 200 bp/100 aa)
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
  options: "--6merpair --retree 1"   # fast, good for large diverse sets

msa_100:
  tool: mafft
  options: "--retree 1"

tree_500:
  tool: fasttree
  options: ""
  model_nuc: GTR+G
  model_aa: WAG+G

tree_100:
  tool: iqtree
  options: "--fast"             # SH-aLRT support; compatible with --fast
  model_nuc: GTR+G
  model_aa: WAG+G
```

Known segmented and large-genome DNA families are auto-configured at config generation time. For segmented families the correct segment keyword is set automatically; for large DNA virus families the appropriate marker gene is selected.

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

# Run pipeline (1 family at a time)
vfam_trees run -f families.txt

# Run with 4 parallel families, 4 threads each
vfam_trees run -f families.txt -j 4 -t 4

# Force re-run families already marked as done
vfam_trees run -f families.txt --force

# Check progress
vfam_trees status -f families.txt
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

For each family, results are written to `results/<Family>/`:

| File | Description |
|------|-------------|
| `<Family>_tree_500.xml` | PhyloXML tree (broad, up to 500 sequences) |
| `<Family>_tree_100.xml` | PhyloXML tree (collapsed, up to 100 sequences) |
| `<Family>_tree_500.nwk` | Newick tree (broad) |
| `<Family>_tree_100.nwk` | Newick tree (collapsed) |
| `<Family>_alignment_500.fasta` | MAFFT alignment (broad) |
| `<Family>_alignment_100.fasta` | MAFFT alignment (collapsed) |
| `<Family>_sequences_raw_500.fasta` | QC-filtered sequences before alignment (broad) |
| `<Family>_sequences_raw_100.fasta` | QC-filtered sequences before alignment (collapsed) |
| `<Family>_metadata_500.tsv` | Sequence metadata (broad) |
| `<Family>_metadata_100.tsv` | Sequence metadata (collapsed) |
| `<Family>_id_map.tsv` | Short ID → display name mapping |
| `<Family>_report.pdf` | Per-family PDF report (stats, length histogram, SH support, tree plot) |
| `<Family>.log` | Per-family log |

A cross-family summary TSV is written to `results/summary.tsv` and updated after each family completes (including families that were skipped due to no species or too few sequences). Key columns include species counts, QC exclusion breakdown, sequence length statistics, clustering thresholds, MSA statistics, and SH support statistics for both trees.

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.
