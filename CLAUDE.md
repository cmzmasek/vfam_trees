# vfam_trees — CLAUDE.md

## Project overview

`vfam_trees` (v1.2.37) is a Python/Snakemake pipeline that builds maximum-likelihood phylogenetic trees for viral families using sequences downloaded from NCBI.  For each family it produces two trees: a broad **tree_500** (≤500 seqs, FastTree GTR+G or WAG+G, SH-like support) and a tighter **tree_100** (≤100 seqs, IQ-TREE `--fast`, SH-aLRT support).

Pre-production — there is no deployed user base.  Changes that affect config schema only need to update templates (`DEFAULT_FAMILY_CONFIG` in `config.py`, `GLOBAL_CONFIG_TEMPLATE` in `cli.py`).  Do **not** add migration/backfill logic to load paths.

---

## Running

```bash
# Install
pip install -e .

# Typical workflow
vfam_trees test                              # verify all external tools
vfam_trees init-configs -f families.txt     # generate per-family YAMLs
vfam_trees run -f families.txt -j 4 -t 4   # run pipeline (4 families, 4 threads each)
vfam_trees status -f families.txt           # progress report
vfam_trees cache --info                     # inspect sequence cache

# Run tests
python -m pytest tests/ -v
```

`--force` on `run` deletes `.done_<family>` sentinels and reruns; without it the Snakefile skips already-finished families.

---

## Environment

- Python 3.11 (Anaconda), macOS (zsh)
- NumPy 2.x is installed — **avoid** `ete3`, `scipy`, `pandas+pyarrow` (compiled against NumPy 1.x)
- External tools that must be in PATH: `mafft`, `FastTree`, `iqtree2`, `mmseqs`

---

## Module pipeline (per family)

```
fetch.py → quality.py → rename.py → subsample.py → msa.py → trim.py
→ tree.py → branch_outliers.py → taxonomy.py → phyloxml_writer.py
→ colors.py → report.py → summary.py
```

Entry point: `vfam_trees/pipeline.py:run_family()` (called by Snakemake rules in `workflow/Snakefile`).  Concatenated-marker path: `pipeline_concat.py` / `concat.py`.

---

## Config system

Priority (lowest → highest): `DEFAULT_FAMILY_CONFIG` → global defaults → `DNA_FAMILIES`/`SEGMENTED_FAMILIES` overrides → user YAML file.  Missing keys are filled at runtime by `_merge_with_defaults`; no migration ever rewrites on-disk files.

- **Global**: `config/global.yaml` — NCBI email/API key, cache settings, default thresholds
- **Per-family**: `configs/{Family}.yaml` — auto-generated from `DEFAULT_FAMILY_CONFIG` if absent

Key config sections: `download`, `sequence`, `quality`, `clustering`, `targets`, `msa_500`, `msa_100`, `tree_500` (`options_nuc` / `options_aa`), `tree_100` (`options_nuc` / `options_aa`), `outlier_removal`, `length_outlier`, `labeling`, `manual`

---

## Key design decisions

### Short IDs
Short IDs (e.g. `S000001`) are assigned in `rename.py` and used throughout the pipeline; display names are restored at output time.

### Length-outlier filter (v1.2.20+)
`quality.remove_length_outliers` uses **MAD on log-lengths + hard floor** (union of two keep windows).
- MAD window: `exp(median(log L) ± k · σ_log)` where `σ_log = 1.4826 · MAD(log L)`
- Hard floor: `[min_lo_mult, max_hi_mult] × median`
- Defaults: `k=5.0`, `min_lo_mult=0.20`, `max_hi_mult=5.0`
- Single source of truth: `quality.log_mad_cutoffs`; both single-protein and concat paths import it
- `k=0` → floor-only; setting both multipliers to 0 → pure MAD

### Branch-length outlier removal
Iterative post-tree removal: threshold = `median + factor × MAD`.  Controlled by `outlier_removal:` YAML block (`enabled=True`, `factor=20.0`, `max_iterations=3`, `min_seqs=40`).  Checkpoints `.msa_done` / `.tree_done` are cleared when a removal triggers rerun.

### Inline sequence injection (`include_seq` / `include_fasta_files`)
`manual.include_seq` accepts inline `{id, organism, sequence}` mappings. `manual.include_fasta_files` accepts paths to FASTA files; `rec.id` (BioPython first token) is the sequence id, `rec.description[len(id):].strip()` is the organism. A single-token header yields an empty organism (the id is *not* reused — doing so duplicated it in the leaf label); such entries are bucketed by id. Both bypass QC and receive downstream protection identical to `manual.include`. Collision detection covers NCBI-fetched accessions and cross-source duplicates. Not supported in concat mode. Config validation checks path strings are non-empty; file I/O happens at pipeline runtime. Stale configs with `include_fasta` key emit a WARNING and are otherwise ignored. Single source of truth for injection: `pipeline._inject_pasted_sequences`; file loading: `pipeline._load_fasta_file_entries`.

### Leaf labeling (v1.2.34+)
`labeling.format` is a Python format string with placeholders `{species}`, `{id}`, `{host}`, `{strain}`, `{location}`, `{year}`, `{genus}`.  Default: `{species}|{id}|{host}`.  `replace_whitespace` (default `true`) and `keep_separator_on_empty` (default `false`) are companion options.

### Leaf coloring
`colors.py` uses HLS colorspace (stdlib `colorsys`).  Same base hue per subfamily, varying lightness per genus.  Gray (`#888888`) for leaves with no genus info.  Applied to PDF/PNG tree images and PhyloXML as `style:font_color` vipr property.  `Phylo.draw` must be called with `label_func=lambda c: c.name or ""` — BioPython's default truncates names >40 chars, breaking the name→color lookup.

### DNA families
`DNA_FAMILIES` in `config.py` maps large-genome families to marker protein queries; small/overlapping-ORF families use whole-genome nucleotide.  `SEGMENTED_FAMILIES` maps multi-segment families to their default phylogenetic segment (title keyword used in NCBI query).

### Sequence cache
`cache.py` — keyed by `(taxid, db, region, segment, max_per_species)`, TTL-based, per-entry lock files.

---

## Output per family (`results/{Family}_{taxid}/`)

| File | Contents |
|------|----------|
| `{Family}_sequences_raw_{500,100}.fasta` | Sequences entering MSA |
| `{Family}_alignment_{500,100}.fasta` | Final MAFFT alignment |
| `{Family}_tree_{500,100}.nwk` / `.xml` | Newick and PhyloXML trees |
| `{Family}_tree_{500,100}.pdf` / `.png` | Tree images with genus color legend |
| `{Family}_metadata_{500,100}.tsv` | Per-sequence metadata |
| `{Family}_id_map.tsv` | Short-ID ↔ accession map |
| `{Family}_report.pdf` | Stats, histograms, colored tree_100 |
| `{Family}.log` | Run log |

Cross-family: `results/summary.tsv`

---

## Tests

```bash
python -m pytest tests/ -v
```

Test files in `tests/`: `test_branch_outliers`, `test_cache`, `test_cli`, `test_colors`, `test_concat`, `test_config`, `test_fetch`, `test_logger`, `test_markers`, `test_phyloxml_writer`, `test_pipeline_concat`, `test_pipeline_helpers`, `test_quality`, `test_rename`, `test_report`, `test_subsample`, `test_summary`, `test_taxonomy`, `test_taxonomy_helpers`, `test_tree_helpers`.

Synthetic quality tests should use varying lengths (see `_NORMAL_LENGTHS` in `test_quality.py`) so MAD has signal — `[L]*N` with all identical lengths hits MAD=0 and exercises only the floor path.

---

## Versioning

Single source of truth: `vfam_trees/__init__.py:__version__`.  `setup.py` reads it via regex.  `CHANGELOG.md` documents user-visible changes.
