# Changelog

All notable user-facing changes are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The full per-commit history is in `git log`; this file summarizes only
changes that affect outputs, configuration, CLI behaviour, or other
user-visible surface area.

## [Unreleased]

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
