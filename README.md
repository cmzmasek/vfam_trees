# vfam_trees

**vfam_trees** is a bioinformatics pipeline for building maximum-likelihood phylogenetic trees for viral families. For each family it automatically discovers species via NCBI Taxonomy, downloads sequences from GenBank, applies quality filtering, clusters sequences per species, aligns them with MAFFT, infers trees with FastTree and IQ-TREE, annotates internal nodes with LCA-based taxonomy, and writes output in Newick and PhyloXML formats.

Two trees are produced per family:

- **tree_500** — broad diversity tree (up to 500 sequences, FastTree / GTR+G)
- **tree_100** — collapsed representative tree (up to 100 sequences, IQ-TREE / GTR+G)

## Features

- Automatic species discovery from NCBI Taxonomy
- Per-species sequence download with RefSeq priority
- Adaptive per-species clustering (MMseqs2) with binary search for optimal identity threshold
- Proportional cross-species sampling to fill target tree sizes
- MAFFT multiple sequence alignment
- FastTree (tree_500) and IQ-TREE (tree_100) tree inference
- Taxonomy-guided tree rooting using LCA specificity scoring (with MAD and midpoint fallbacks)
- LCA-based internal node annotation using NCBI ranked lineages
- PhyloXML output with `<taxonomy>` elements and `vipr:` metadata properties
- Pre-configured support for 30+ segmented RNA virus families and 19 DNA virus families
- Per-run summary TSV with bootstrap statistics, MSA statistics, and sequence length distributions

## Dependencies

### Python packages

```
biopython >= 1.81
click >= 8.1
pyyaml >= 6.0
snakemake >= 7.0
```

### External tools

| Tool | Purpose |
|------|---------|
| `mafft` | Multiple sequence alignment |
| `FastTree` | Rapid ML tree inference (tree_500) |
| `iqtree2` | ML tree inference with model selection (tree_100) |
| `mmseqs` | Sequence clustering |

All tools must be available on `$PATH`. Installation via conda is recommended:

```bash
conda install -c bioconda mafft fasttree iqtree mmseqs2
```

## Installation

```bash
git clone https://github.com/<your-username>/vfam_trees.git
cd vfam_trees
pip install -e .
```

## Configuration

### 1. Global config (`config/global.yaml`)

Create this file before running (it is excluded from the repository as it may contain your NCBI API key):

```yaml
ncbi:
  email: your.email@example.com
  api_key: your_ncbi_api_key   # optional but recommended (10 req/s vs 3 req/s)
```

An NCBI API key can be obtained for free at https://www.ncbi.nlm.nih.gov/account/

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
  region: whole_genome          # whole_genome, or a gene name (e.g. hexon)
  segment: null                 # segment keyword for segmented viruses (e.g. "segment L")

quality:
  min_length: null              # minimum sequence length (null = no filter)
  max_ambiguous: 0.01           # maximum fraction of ambiguous bases

clustering:
  tool: mmseqs2
  threshold_min: 0.70           # minimum clustering identity
  threshold_max: 0.99           # maximum clustering identity
  max_reps_500: 20              # max representatives per species for tree_500
  max_reps_100: 5               # max representatives per species for tree_100

targets:
  max_500: 500                  # target sequences for tree_500
  max_100: 100                  # target sequences for tree_100

msa:
  tool: mafft
  options: "--retree 1"

tree_500:
  tool: fasttree
  model_nuc: GTR+G
  model_aa: WAG+G

tree_100:
  tool: iqtree
  options: "--fast"
  model_nuc: GTR+G
  model_aa: WAG+G
```

## Usage

```bash
# Check all dependencies
vfam_trees test

# Generate per-family configs (review before running)
vfam_trees init-configs -f families.txt

# Run the pipeline (1 family at a time)
vfam_trees run -f families.txt

# Run with 4 parallel families, 4 threads each
vfam_trees run -f families.txt -j 4 -t 4

# Check progress
vfam_trees status -f families.txt
```

The `families.txt` file should contain one ICTV family name per line:

```
Flaviviridae
Coronaviridae
Filoviridae
```

## Output

For each family, results are written to `results/<Family>/`:

| File | Description |
|------|-------------|
| `<Family>_tree_500.xml` | PhyloXML tree (broad, 500 sequences) |
| `<Family>_tree_100.xml` | PhyloXML tree (collapsed, 100 sequences) |
| `<Family>_tree_500.nwk` | Newick tree (broad) |
| `<Family>_tree_100.nwk` | Newick tree (collapsed) |
| `<Family>_alignment_500.fasta` | MAFFT alignment (broad) |
| `<Family>_alignment_100.fasta` | MAFFT alignment (collapsed) |
| `<Family>_sequences_raw_500.fasta` | QC-filtered sequences (broad) |
| `<Family>_sequences_raw_100.fasta` | QC-filtered sequences (collapsed) |
| `<Family>_metadata_500.tsv` | Sequence metadata (broad) |
| `<Family>_metadata_100.tsv` | Sequence metadata (collapsed) |
| `<Family>_id_map.tsv` | Short ID → display name mapping |
| `<Family>.log` | Per-family log |

A cross-family summary TSV is written to `results/summary.tsv` and updated after each family completes.

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.
