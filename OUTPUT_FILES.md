# Per-family output files — detailed description

For each viral family the pipeline creates a directory `results/<Family>_<taxid>/` (e.g. `Hantaviridae_1980413/`). Two trees are produced per family — a **broad tree (tree_500)** with up to 500 sequences and a **collapsed tree (tree_100)** with up to 100 representative sequences — so most outputs come in matched pairs.

## 1. Tree files

### Newick (`<Family>_tree_500.nwk`, `<Family>_tree_100.nwk`)
Plain text trees in standard Newick format. Internal nodes are annotated with the **last common ancestor (LCA)** taxon name (e.g. genus, subfamily) wherever the taxonomy of the descendant leaves agrees. Branch lengths are in expected substitutions per site. Support values are attached to internal nodes (range and meaning depend on the tree tool — see PhyloXML section). These are intended for downstream programmatic use (Bio.Phylo, ETE, dendropy, R `ape`, etc.).

### PhyloXML (`<Family>_tree_500.xml`, `<Family>_tree_100.xml`)
Richly-annotated phylogenetic XML. **This is the canonical output format** — it carries the tree topology *plus* per-leaf metadata, taxonomy, and visual styling that a Newick file cannot. Detailed structure is given in §6 below.

### Tree images — rooted rectangular (`*_tree_500.pdf` / `.png`, `*_tree_100.pdf` / `.png`)
Conventional left-to-right rectangular cladograms, **rooted** by the pipeline's taxonomy-guided algorithm (with MAD and midpoint fallbacks). Support values below 50% are suppressed for visual clarity. Leaf labels are colored by genus (one hue band per subfamily; lightness varies across genera within a subfamily). A structured legend grouped by subfamily appears alongside the tree.

### Tree images — unrooted radial (`*_tree_500_ur.pdf` / `.png`, `*_tree_100_ur.pdf` / `.png`)
Radial (fan) layout with no root and no internal-node or support labels — useful for inspecting **overall topology** and **divergence** without privileging any one root hypothesis.

### Tree icon (`<Family>_tree_icon.png`)
Square topology-only thumbnail of `tree_100` (no labels, uniform line color). Used to build the cross-family `overview_tree_100.png` mosaic at the top level.

## 2. Sequence and alignment files

### Raw QC'd sequences (`<Family>_sequences_raw_500.fasta`, `<Family>_sequences_raw_100.fasta`)
Sequences that *entered* the MSA, after: organism-name exclusion (e.g. synthetic constructs, patents), length filtering, ambiguity filtering, deduplication, per-species clustering with MMseqs2, and proportional cross-species subsampling. RefSeq records are preferred during subsampling and **never removed** by the length-outlier or branch-length-outlier filters. Headers carry pipe-delimited display names: `<species>|<strain>|<accession>|<host>`.

### Final alignment (`<Family>_alignment_500.fasta`, `<Family>_alignment_100.fasta`)
The MSA actually fed to the tree-building step — MAFFT alignment, optionally column-trimmed by **trimAl** (`-automated1` by default). Reflects the final iteration after any post-tree branch-length outlier removal. Length and gap percentage are reported in `summary.tsv`.

## 3. Metadata files

### Sequence metadata TSV (`<Family>_metadata_500.tsv`, `<Family>_metadata_100.tsv`)
One row per sequence in the corresponding tree. Columns include `short_id`, `display_name`, `accession`, `species`, `strain`, `host`, `collection_date`, `location`, `taxon_id`, `length`. This is a flat, spreadsheet-friendly view of what's encoded richly in the PhyloXML.

### ID map (`<Family>_id_map.tsv`)
Three columns: `short_id` (used internally during MSA/tree inference because phylogenetics tools choke on long names), `accession`, `display_name`. Restores the human-readable label after tree inference.

### Color map (`<Family>_colors_500.json`, `<Family>_colors_100.json`)
JSON mapping each leaf's `display_name` to a hex color. The same colors are reused in the PDF/PNG figures and the PhyloXML `style:font_color` property.

## 4. Configuration & log

### Resolved per-family config (`<Family>.yaml`)
The fully merged configuration actually used for this family — global config + family overrides + defaults. Useful for reproducibility and for re-running a single family with the same settings.

### Per-family log (`<Family>.log`)
Verbose, time-stamped log of every pipeline step for this family: NCBI fetches, QC counts, clustering thresholds tried, MSA + tree commands, outlier removal decisions, RefSeq-protection warnings, rooting choice, and timing.

## 5. Reports

### Per-family report (`<Family>_report.pdf`)
Multi-page PDF with:
1. **Stats table** — species discovered, sequences kept after each QC stage, clustering identity, tree leaf counts, support distribution.
2. **Post-QC sequence length histogram** — pre-MSA distribution.
3. **Per-tree sequence-length histograms** for tree_500 and tree_100.
4. **SH support histograms** for tree_500 and tree_100.
5. **tree_100 visualization** with the genus/subfamily color legend.

### Sequence-length plot (`<Family>_seqlen_plot.pdf`)
Standalone histogram of post-QC sequence lengths.

---

## 6. PhyloXML files — deep dive

PhyloXML is an XML schema for phylogenetic trees. Standard viewers include **Archaeopteryx** (and its web port phyloT/phylo.io readers), **Forester**, and **Dendroscope**. Below is a description of every element the pipeline emits, what biological information it carries, and how a viewer renders it.

### Document root

```xml
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="true" rerootable="false">
    ...
  </phylogeny>
</phyloxml>
```

- **`phyloxml`** — outer container with the PhyloXML namespace.
- **`phylogeny rooted="true" rerootable="false"`** — declares the tree is rooted *and* asks viewers not to reroot it. The pipeline's root is chosen by a taxonomy-guided LCA-specificity algorithm (with MAD and midpoint fallbacks); this attribute prevents the viewer from silently overriding that choice.

### Phylogeny header

```xml
<name>Hantaviridae [nucleotide|segment L|MAFFT|FastTree GTR+G]</name>
<description>vfam_trees v1.1.4 phylogenetic tree for Hantaviridae —
  generated 2026-04-23 20:45 UTC | nucleotide, segment L
  (target=500 seqs, cluster id 0.70–0.99) |
  MAFFT 7.520 --6merpair --retree 1 |
  FastTree 2.1.11 GTR+G</description>
```

- **`<name>`** — short tree title encoding family + molecule + region/segment + MSA tool + tree tool/model in pipe-delimited tags. Designed to be readable at a glance in tree viewers.
- **`<description>`** — full provenance: pipeline version, family name, UTC timestamp, target tree size, clustering identity range actually used, exact MSA tool + version + options, exact tree tool + version + substitution model. **This is your reproducibility breadcrumb** — if you publish a tree, cite this string verbatim.

### Clade element — the tree itself

The tree topology is a recursive nesting of `<clade>` elements. Each `<clade>` represents either a **leaf** (terminal taxon) or an **internal node** (LCA of its descendants). The PhyloXML schema fixes child-element order: `name`, `branch_length`, `confidence`, `taxonomy`, `sequence`, `property`, then nested `<clade>` children.

#### Common to leaves and internals

- **`<branch_length>`** — length of the branch *leading to this node*, in expected substitutions per site (the standard ML units). Visible as horizontal distance in rectangular layouts.
- **`<confidence type="...">`** — branch support value, integer 0–100. The `type` attribute tells you *how* to interpret it:
  - `type="SH_aLRT"` (IQ-TREE) — Shimodaira-Hasegawa-like approximate Likelihood Ratio Test. Conservative: ≥80 is well-supported.
  - `type="SH_like"` (FastTree) — FastTree's local Shimodaira-Hasegawa support. Calibrated similarly: ≥80 is well-supported, 50–80 is uncertain, <50 typically suppressed in figures.
  - `type="UFBoot"` (IQ-TREE with `-bb`) — Ultrafast bootstrap. Higher threshold: ≥95 is well-supported.

#### Leaves

```xml
<clade>
  <name>Orthohantavirus_andesense|CHI-7913|PV808472.1|Homo_sapiens</name>
  <branch_length>0.013183853</branch_length>
  <taxonomy>
    <id provider="ncbi">1980456</id>
    <scientific_name>Orthohantavirus andesense</scientific_name>
  </taxonomy>
  <sequence>
    <accession source="ncbi">PV808472.1</accession>
    <name>Orthohantavirus andesense isolate CHI-7913 segment L, complete sequence</name>
  </sequence>
  <property ref="vipr:Host"            datatype="xsd:string" applies_to="node">Homo sapiens</property>
  <property ref="vipr:Collection_Date" datatype="xsd:string" applies_to="node">1999</property>
  <property ref="vipr:Strain"          datatype="xsd:string" applies_to="node">CHI-7913</property>
  <property ref="vipr:Year"            datatype="xsd:string" applies_to="node">1999</property>
  <property ref="vipr:Species"         datatype="xsd:string" applies_to="node">Orthohantavirus andesense</property>
  <property ref="vipr:Genus"           datatype="xsd:string" applies_to="node">Orthohantavirus</property>
  <property ref="vipr:Subgenus"        datatype="xsd:string" applies_to="node">...</property>
  <property ref="vipr:Subfamily"       datatype="xsd:string" applies_to="node">...</property>
  <property ref="style:font_color"     datatype="xsd:token"  applies_to="node">#21e6e6</property>
</clade>
```

- **`<name>`** — pipe-delimited display name `species|strain|accession|host`. This is what shows up as the leaf label in tree viewers.
- **`<taxonomy>`**:
  - **`<id provider="ncbi">`** — NCBI Taxonomy ID of the species (clickable in some viewers; useful for joining to other datasets).
  - **`<scientific_name>`** — official species binomial.
- **`<sequence>`**:
  - **`<accession source="ncbi">`** — GenBank accession with version (e.g. `PV808472.1`).
  - **`<name>`** — full GenBank record title (the DEFINITION line — handy when scanning for "complete genome" vs. partial).
- **`<property>` elements** — every leaf carries a uniform set of typed key-value pairs. The `ref` namespace prefix (`vipr:`, `style:`) is a viewer-recognized convention; `datatype="xsd:string"` is the XSD string type; `applies_to="node"` means the property is associated with this leaf node (rather than the branch leading to it).
  - **`vipr:Host`** — host organism extracted from the GenBank `/host` feature qualifier. Often the most useful leaf annotation for ecologists / virologists.
  - **`vipr:Collection_Date`** — full collection date as recorded in GenBank (`/collection_date`, e.g. `21-Dec-2009` or `1999`).
  - **`vipr:Year`** — 4-digit collection year, parsed from `Collection_Date` for easy temporal coloring/filtering.
  - **`vipr:Location`** — geographic collection location from `/country` or `/geo_loc_name` (when present).
  - **`vipr:Strain`** — isolate/strain designation from `/strain` or `/isolate`.
  - **`vipr:Species`**, **`vipr:Genus`**, **`vipr:Subgenus`**, **`vipr:Subfamily`** — taxonomic ranks resolved against the NCBI Taxonomy lineage. These let viewers color or group leaves by any rank, independent of tree topology. (Empty ranks are simply omitted — don't be surprised if some leaves lack a `Subgenus` property.)
  - **`style:font_color`** — hex color (e.g. `#21e6e6`) applied to the leaf label by viewers that honor the `style:` namespace (Archaeopteryx does). Colors come from the same HLS scheme used in the PDF/PNG figures, so the family-level color story is consistent across all output formats.

> **Note on the `unknown` placeholder**: the pipeline's metadata extractor sometimes emits `"unknown"` as a fallback. The PhyloXML writer treats `unknown` as absent — it will *not* emit `<scientific_name>unknown</scientific_name>` or a `vipr:Host` of `unknown`. If a property is missing, the underlying GenBank record didn't carry that information.

#### Internal nodes (LCA-annotated)

```xml
<clade>
  <name>Orthohantavirus</name>
  <branch_length>0.114886822</branch_length>
  <confidence type="SH_like">50</confidence>
  <taxonomy>
    <scientific_name>Orthohantavirus</scientific_name>
    <rank>genus</rank>
  </taxonomy>
  <clade>...</clade>
  <clade>...</clade>
</clade>
```

- **`<name>`** — taxon name where the LCA of all descendants agrees on a *named* rank (e.g. genus `Orthohantavirus`, family `Hantaviridae`, subfamily `Mammantavirinae`). This is what lets you visually trace where each clade sits in the standard taxonomy without consulting an external resource.
- **`<taxonomy>`** with **`<scientific_name>`** + **`<rank>`** — the same name re-emitted in PhyloXML's structured taxonomy element, with `<rank>` drawn from the schema's controlled vocabulary (`family`, `subfamily`, `genus`, `subgenus`, `species`, …). Ranks not in the schema collapse to `other`.
- **Suppression of trivial labels**: when an internal node has exactly two leaf children of the same species, the LCA label is redundant and is **omitted** (no `<name>` element on that node). Keeps the tree visually uncluttered. Species-level internal labels are also filtered from PDF/PNG figures.

### What the PhyloXML lets you do that the Newick doesn't

1. **Color leaves by host, year, country, or genus** without parsing the leaf label string — point your viewer at the relevant `vipr:` property.
2. **Click through to NCBI** — Archaeopteryx and other viewers can hyperlink `<accession source="ncbi">` and `<id provider="ncbi">` directly to the appropriate NCBI page.
3. **Filter or hide clades by rank** — tools that read `<taxonomy><rank>` can collapse all subtrees below a chosen rank (e.g. show only genus-level internal nodes).
4. **Reproducibility** — the `<description>` string captures pipeline version, MSA tool + options, tree tool + model, and clustering identity used.
5. **Round-trip safety** — taxonomy IDs and full-precision branch lengths survive in PhyloXML; Newick label-encoded annotations often don't.
