# CONCAT_DESIGN.md — Multi-marker protein concatenation mode (MVP)

> Design specification for adding concatenated multi-protein phylogenies to *vfam_trees*. Targets large DNA virus families (Poxviridae, Herpesviridae, Asfarviridae, Iridoviridae, Baculoviridae, NCLDVs broadly) where single-protein analysis lacks sufficient phylogenetic signal. v2 will swap in HMMER-based marker identification; this design keeps that swap clean.

## 1. Goals and scope

- Build family-scale phylogenies from concatenated alignments of N curated marker proteins per genome (typically 6–9), where each leaf is one **genome** rather than one protein record.
- Provide reasonable per-family default marker sets so auto-generated configs work out of the box for the targeted families.
- Preserve all existing behavior for single-protein and whole-genome families.
- Use **partitioned tree inference** for tree_100 so each marker gets its own substitution model.
- Architect the marker-identification step behind a small interface so v2 (HMMER) drops in without touching the fetcher / concat / tree code.

## 2. Decisions locked in (from design discussion)

| # | Topic | Decision |
|---|---|---|
| Q1 | Marker sets | Curated per family; user can override or revert any family to single-protein mode |
| Q1a | Per-marker schema | `name`, `aliases[]`, optional `length_range`, optional `locus_tag_hint` |
| Q1b | Herpesviridae | Single 7-marker family-wide set (no subfamily split) |
| Q1c | Poxviridae | 9-marker set covering both subfamilies; entomopoxviruses retained, low-coverage handled by `min_fraction` |
| Q2 | Per-species policy | **Policy B** — adaptive clustering operates on the per-genome concatenated sequence; max_reps caps how many genomes any one species contributes |
| Q3 | Genome grouping | **Policy A** — strict source-nucleotide accession grouping; legacy split-submission specimens are dropped, count logged |
| Q4 | RefSeq policy in concat mode | Genome-level: a "RefSeq genome" is one whose source nuc accession matches the RefSeq prefix; uncapped fetch; absorption on concats ≥0.99; branch-length-outlier exempt |
| Q4a | Length-outlier filter | **Per-marker** before concatenation; RefSeq markers exempt |
| Q5 | Tree inference | Partitioned IQ-TREE for tree_100 (`-p partitions.nex -m MFP`); see §10 for tree_500 open item |

## 3. User-facing config schema

### 3.1 Family config

```yaml
# Existing single-protein mode (unchanged)
sequence:
  type: protein
  region: "DNA polymerase"   # single marker by name (existing behavior)
  segment: null

# NEW: concatenated mode
sequence:
  type: protein
  region: concatenated       # NEW value; triggers concat path
  segment: null
concatenation:
  proteins:
    - name: "DNA polymerase"
      aliases:
        - "DNA-directed DNA polymerase"
        - "DNA pol"
      aliases_Entomopoxvirinae:        # optional subfamily-aware aliases
        - "DNA polymerase B"
      length_range: [800, 1500]        # optional aa range; outside → reject
      locus_tag_hint: "UL30|polB"      # optional regex; tiebreaker for paralogs
    - name: "major capsid protein"
      aliases: ["MCP", "VP39"]
    # ... more markers
  min_fraction: 0.7          # genomes with < min_fraction × N markers are dropped
  partition_tree_100: true   # partitioned IQ-TREE on tree_100 (default)
  partition_tree_500: false  # FastTree cannot partition; tree_500 stays single-model
```

### 3.2 Auto-generated default

For families in the curated `CONCATENATION_FAMILIES` table (§4), auto-generated family configs ship `region: concatenated` with the family-specific marker list. For all other families, defaults remain `region: whole_genome` or single-protein per existing logic.

### 3.3 User overrides

Either direction works:
- **Concat → single-protein**: user edits family yaml, sets `region: <single protein name>`, deletes `concatenation:` block.
- **Single-protein → concat**: user adds `region: concatenated` plus their own `concatenation: { proteins: [...] }` list. No requirement for the family to be in `CONCATENATION_FAMILIES`.

## 4. Curated marker presets (`CONCATENATION_FAMILIES` in `config.py`)

The curated content is reproduced here for review. Locus tags use the canonical reference genome (HSV-1 for Herpesviridae, VACV-Cop for Poxviridae, ASFV-BA71V for Asfarviridae, AcMNPV for Baculoviridae, FV3 for Iridoviridae).

> **Status**: first-pass curation from training-data familiarity with the literature. Validate before merging. ASFV gene IDs are flagged TODO for explicit verification.

### 4.1 Poxviridae (9 markers)

References: Upton et al. 2003; Hughes & Friedman 2005; ICTV Poxviridae chapter.

| Marker | VACV gene | Aliases | Locus tag hint |
|---|---|---|---|
| DNA polymerase | E9L | "DNA polymerase", "DNA-directed DNA polymerase" | `E9L\|polB` |
| RNA polymerase RPO147 | A24R | "DNA-directed RNA polymerase 147 kDa subunit", "RPO147", "RNA polymerase subunit RPO147" | `A24R\|RPO147` |
| RNA polymerase RPO132 | J6R | "DNA-directed RNA polymerase 132 kDa subunit", "RPO132", "RNA polymerase subunit RPO132" | `J6R\|RPO132` |
| mRNA capping enzyme large subunit | D1R | "mRNA capping enzyme large subunit", "capping enzyme large subunit" | `D1R` |
| DNA helicase / NPH-II | A18R | "DNA helicase", "NPH-II", "transcript release factor" | `A18R\|NPH2` |
| Poly(A) polymerase catalytic subunit | E1L | "poly(A) polymerase catalytic subunit", "poly(A) polymerase large subunit" | `E1L\|VP55` |
| Late transcription factor VLTF-3 | A1L | "late transcription factor VLTF-3", "VLTF3", "late transcription factor 3" | `A1L\|VLTF` |
| Uracil-DNA glycosylase | D4R | "uracil-DNA glycosylase", "UNG" | `D4R\|UNG` |
| Single-stranded DNA-binding protein | I3L | "single-stranded DNA-binding protein", "ssDNA binding protein" | `I3L\|ssb` |

### 4.2 Herpesviridae (7 markers)

References: McGeoch et al. 1995, 2006; Davison 2010.

| Marker | HSV-1 ORF | Aliases | Locus tag hint |
|---|---|---|---|
| DNA polymerase catalytic subunit | UL30 | "DNA polymerase catalytic subunit", "DNA-directed DNA polymerase catalytic subunit" | `UL30` |
| Helicase-primase helicase subunit | UL5 | "helicase-primase helicase subunit", "DNA helicase" | `UL5\b` |
| Helicase-primase primase subunit | UL52 | "helicase-primase primase subunit", "primase" | `UL52` |
| Major capsid protein | UL19 | "major capsid protein", "MCP", "capsid protein VP5" | `UL19\|VP5\b` |
| Capsid triplex subunit 2 (VP23) | UL18 | "capsid triplex subunit 2", "VP23", "minor capsid protein" | `UL18\|VP23` |
| DNA packaging terminase ATPase | UL15 | "DNA packaging terminase subunit 1", "terminase ATPase subunit" | `UL15` |
| Single-stranded DNA-binding protein | UL29 | "single-stranded DNA-binding protein", "major DNA-binding protein", "ICP8" | `UL29\|ICP8` |

Glycoprotein B (UL27) deliberately excluded — too divergent across alpha/beta/gamma subfamilies for reliable family-wide alignment.

### 4.3 Asfarviridae (6 markers — ASFV gene IDs TBV)

References: Yutin & Koonin 2012; Iyer et al. 2006; ICTV Asfarviridae chapter.

| Marker | ASFV gene | Aliases | Locus tag hint |
|---|---|---|---|
| DNA polymerase B | (TBV) | "DNA polymerase B", "DNA-directed DNA polymerase" | `polB` |
| Major capsid protein p72 | B646L | "major capsid protein", "p72", "capsid protein p72" | `B646L\|p72` |
| Packaging ATPase | (TBV) | "packaging ATPase", "A32-like ATPase", "FtsK-like ATPase" | `A32` |
| D5/A18 superfamily-3 helicase | (TBV) | "primase-helicase", "D5-like helicase", "superfamily 3 helicase" | `D5\|A18` |
| Late transcription factor VLTF-3 | (TBV) | "late transcription factor 3", "VLTF-3" | `VLTF` |
| RNA polymerase largest subunit | NP1450L | "DNA-directed RNA polymerase subunit", "RPB1-like subunit", "RNA polymerase largest subunit" | `NP1450L\|RPB1` |

Gene IDs marked `(TBV)` to be verified against NC_001659 (ASFV BA71V RefSeq).

### 4.4 Iridoviridae (7 markers)

References: Tidona & Darai 1997; Eaton et al. 2007; ICTV Iridoviridae chapter.

| Marker | Aliases | Locus tag hint |
|---|---|---|
| Major capsid protein | "major capsid protein", "MCP" | `MCP` |
| DNA polymerase | "DNA polymerase", "DNA-directed DNA polymerase" | `polB` |
| Packaging ATPase | "packaging ATPase", "A32-like ATPase" | `A32` |
| RNase III-like | "ribonuclease III", "RNase III" | `rnc` |
| DNA helicase | "DNA helicase", "D5-like helicase" | `D5\|helicase` |
| Late transcription factor VLTF-3 | "late transcription factor 3", "VLTF-3" | `VLTF` |
| Immediate-early protein ICP-46 | "immediate-early protein ICP-46", "ICP46" | `ICP46` |

### 4.5 Baculoviridae (7 markers)

References: Herniman et al. 2003; Jehle et al. 2006; Miele et al. 2011.

| Marker | AcMNPV ORF | Aliases | Locus tag hint |
|---|---|---|---|
| DNA polymerase | ac65 | "DNA polymerase", "DNA-directed DNA polymerase" | `polB` |
| Late expression factor 8 | ac50 (lef-8) | "late expression factor 8", "LEF-8", "RNA polymerase subunit LEF-8" | `lef-?8` |
| Late expression factor 9 | ac62 (lef-9) | "late expression factor 9", "LEF-9", "RNA polymerase subunit LEF-9" | `lef-?9` |
| DNA helicase P143 | ac95 | "DNA helicase P143", "p143", "helicase" | `p143\|helicase` |
| Per os infectivity factor 1 (P74) | ac138 | "per os infectivity factor 1", "PIF-1", "P74" | `pif-?1\|p74` |
| Per os infectivity factor 2 | ac22 | "per os infectivity factor 2", "PIF-2" | `pif-?2` |
| VP39 / major capsid protein | ac89 | "major capsid protein", "VP39", "capsid protein VP39" | `vp39\|MCP` |

### 4.6 NCLDV hallmark fallback (8 markers)

References: Yutin & Koonin 2009, 2012; Koonin & Yutin 2019.

Used as the auto-generated default for any large-DNA-virus family not in §4.1–4.5: Mimiviridae, Phycodnaviridae, Marseilleviridae, Pithoviridae, Pandoraviridae, Ascoviridae, etc.

| Marker | Aliases |
|---|---|
| DNA polymerase B family | "DNA polymerase", "DNA-directed DNA polymerase", "polymerase B" |
| Major capsid protein (double jelly-roll) | "major capsid protein", "MCP", "capsid protein" |
| A32-like packaging ATPase | "packaging ATPase", "A32-like ATPase", "FtsK-like ATPase" |
| D5-like primase-helicase | "primase-helicase", "D5-like helicase", "superfamily 3 helicase" |
| Late transcription factor VLTF-3 | "late transcription factor 3", "VLTF-3" |
| mRNA capping enzyme | "mRNA capping enzyme", "capping enzyme large subunit" |
| RNA polymerase largest subunit (RPB1-like) | "DNA-directed RNA polymerase subunit alpha", "RNA polymerase RPB1", "largest subunit RNA polymerase" |
| RNA polymerase second-largest subunit (RPB2-like) | "DNA-directed RNA polymerase subunit beta", "RNA polymerase RPB2", "second-largest subunit RNA polymerase" |

### 4.7 Families that stay single-protein

Adenoviridae, Papillomaviridae, Polyomaviridae, all ssDNA families, all RNA virus families.

## 5. Module-level architecture

### 5.1 Marker-identification interface (the v2 swap point)

A small Protocol for identifying which proteins in a genome correspond to which markers. v1 implements name+alias matching; v2 will implement HMMER-based identification behind the same interface.

```python
# vfam_trees/markers.py
from typing import Protocol

class MarkerIdentifier(Protocol):
    def identify(
        self,
        proteins: list[SeqRecord],          # all proteins from one genome
        marker_set: list[MarkerSpec],       # the family's curated marker list
    ) -> dict[str, SeqRecord]:              # marker_name → chosen protein record
        ...

class NameMatchIdentifier:
    """v1 — matches on protein name + aliases, with locus_tag_hint and length_range tiebreakers."""

class HmmerIdentifier:
    """v2 — runs hmmsearch against profile HMMs per marker; not implemented in MVP."""
```

The pipeline picks the identifier from `family_cfg["concatenation"].get("identifier", "name_match")`. Default `"name_match"` for v1.

### 5.2 Fetcher contract

```python
# vfam_trees/fetch.py — new function alongside fetch_species_sequences
def fetch_species_genomes(
    taxid: int,
    species_name: str,
    marker_set: list[MarkerSpec],
    output_dir: Path,
    max_per_species: int,
    exclude_organisms: list[str] | None,
) -> dict[str, dict[str, SeqRecord]]:
    """Fetch protein records for each marker; group by source nucleotide accession.

    Returns:
        {source_nuc_accession: {marker_name: protein SeqRecord}}
    """
```

Implementation outline:
1. **Resolve subfamily-aware aliases.** Look up the species' NCBI ranked lineage (already fetched per family). For each marker spec, the effective alias list = `aliases` ∪ `aliases_<subfamily>` (where `<subfamily>` is the species' subfamily-rank entry, if any). This lets one Poxviridae marker spec serve both Chordopoxvirinae and Entomopoxvirinae with subfamily-specific annotation conventions.
2. For each marker, build an Entrez query of the form
   `txid{taxid}[Organism:exp] AND ("{name}"[Protein Name] OR "{alias1}"[Protein Name] OR ...)`
   using the resolved alias list.
3. Fetch RefSeq matches uncapped, non-RefSeq matches capped at `max_per_species`.
4. Parse `coded_by` and `GBSeq_source-db` from each protein record to identify the source nucleotide accession.
5. Group records into `{source_acc: {marker: [candidate proteins]}}`.
6. For each genome × marker with multiple candidates, apply the within-genome tiebreaker (§5.3).
7. Drop genomes whose markers span multiple source nuc accessions (Policy A) — increment `n_genomes_dropped_split_submission` counter.
8. Drop genomes with fewer than `min_fraction × N` markers — increment `n_genomes_dropped_min_fraction` counter.

### 5.3 Within-genome paralog tiebreaker

When a genome has multiple protein records matching one marker's name+aliases, choose by priority order:
1. **Locus-tag hint regex match** (if `locus_tag_hint` is set on the marker spec).
2. **Length closest to `length_range` midpoint** (if `length_range` set).
3. **Longest sequence**.
4. **Lowest accession** (deterministic tiebreaker).

### 5.4 Per-marker QC

Each marker's protein records pass through the existing `filter_sequences` (organism exclusion, ambiguity, deduplication) **before grouping into genomes**. Length filtering uses `length_range` from the marker spec when present, else falls back to the existing adaptive 50%/40%/30%-of-median floor with the 100 aa hard floor.

### 5.5 Per-marker length-outlier (Q4a)

After grouping genomes, compute per-marker median lengths across all surviving genomes. Drop a `(genome, marker)` cell where the marker length is a 2-sided outlier (`>3×` or `<1/3×` median, configurable). RefSeq protein records are exempt — flagged but kept, warning logged. The genome remains in the dataset; the dropped marker becomes a gap-padded block in the concat.

### 5.6 Concat alignment

```
For each marker:
    Run MAFFT on the surviving (genome → marker protein) records — short IDs are genome accessions
    Run trimAl on the per-marker MSA
    Record the trimmed alignment's column count = block_length
For each genome:
    Build the concatenated sequence by emitting each marker block in order:
        - if genome has the marker: the trimmed-aligned sequence for that genome
        - if genome lacks the marker: '-' × block_length
    The concatenated sequence is one row in the family-level MSA
Record partition map: marker_name → (start_col, end_col) in concat
```

Marker order is the order declared in `concatenation.proteins` — stable across runs.

### 5.7 Genome-level adaptive clustering (Policy B)

Run the existing `adaptive_cluster_species` machinery, but on the per-genome concatenated sequences (not on individual proteins). Each leaf coming out is one genome. RefSeq genomes are prioritized as cluster representatives via the existing `priority_ids` mechanism extended to the genome layer.

### 5.8 Genome-level RefSeq absorption

Extend `absorb_into_refseqs` to the genome layer: cluster concats at `refseq_absorption.threshold` (default 0.99); for each cluster containing one or more RefSeq genomes, drop non-RefSeq genomes in the cluster.

### 5.9 Cross-species proportional merge

Unchanged; operates on genome records keyed by species. RefSeq genomes prioritized.

### 5.10 Tree inference (partitioned tree_100; single-model tree_500)

**tree_500** — FastTree on the concatenated alignment, single model (GTR+G nuc / LG+G aa). FastTree does not support partitioned analysis; concat mode keeps the existing rigor-vs-speed split (broad/fast tree_500, representative/rigorous tree_100).

**tree_100** — IQ-TREE on the concatenated alignment, **partitioned**:
```
iqtree2 -s concat.fasta -p partitions.nex -m MFP -B 1000 [-alrt 1000]
```
- `-p` is the **edge-linked proportional** partition model: shared topology, per-partition rate scalars, per-partition substitution model.
- `-m MFP` runs ModelFinder per partition (BIC by default).
- `-B 1000` UFBoot for protein data; `-alrt 1000` SH-aLRT for nucleotide data — same as current logic.

Partition file (NEXUS):
```
#nexus
begin sets;
    charset DNApolymerase = 1-1234;
    charset MCP = 1235-2456;
    charset packagingATPase = 2457-3100;
    ...
end;
```

The partition file is generated automatically from §5.6's recorded marker boundaries.

### 5.11 IQ-TREE log parsing for partitioned models

Existing `parse_iqtree_best_model` returns a single best model. For partitioned analysis, IQ-TREE writes a per-partition model in the `.iqtree` file under a section like:

```
ID  Model            Parameters
1   LG+I+G4          ...
2   WAG+G4           ...
```

New helper `parse_iqtree_partition_models(log_path) -> dict[partition_name, model_string]` will read this section and return per-partition models. Recorded in `summary.tsv` and PhyloXML provenance.

### 5.12 Branch-length outlier removal

Unchanged in logic; leaves are now genomes. RefSeq genomes exempt.

### 5.13 Tree rooting + LCA annotation

Unchanged. LCA uses the source nucleotide record's NCBI lineage, which is identical to what the species-level lineage gives in single-protein mode (proteins inherit lineage from their source nuc). Internal-node labels remain at family/subfamily/genus/species rank.

## 6. Output additions

### 6.1 `summary.tsv` new columns

Per tree (prefix `tree500_` and `tree100_`):

| Column | Type | Description |
|---|---|---|
| `concat_n_markers_target` | int | Markers in the family's curated set |
| `concat_n_markers_used` | int | Markers with ≥1 leaf coverage |
| `concat_length` | int | Concatenated alignment column count post-trim |
| `concat_pct_gaps` | float | Average gap percentage in concat |
| `n_genomes_dropped_split_submission` | int | Specimens dropped by Policy A (multi-accession) |
| `n_genomes_dropped_min_fraction` | int | Specimens dropped for low marker coverage |
| `marker_models` | string | Comma-separated `<marker>:<model>` (tree_100 only, when partitioned) |

Per-marker counts go into a single denormalized column to avoid per-marker schema explosion:

| Column | Type | Description |
|---|---|---|
| `marker_coverage` | string | Comma-separated `<marker>:<n_genomes>` |

### 6.2 PhyloXML additions

- Leaf `<accession source="ncbi">` is the **source nucleotide accession** (the genome ID).
- New leaf property `vipr:Marker_Accessions` — comma-separated list of the protein accessions used in this genome's concat (e.g. `YP_000001.1,YP_000007.1,YP_000019.1`).
- New leaf property `vipr:N_Markers_Present` — integer count.
- Phylogeny `<description>` extended to include `concat: <N> markers (<family preset name or "custom">); partitioned IQ-TREE` etc.

### 6.3 Per-family report PDF

New page: per-marker coverage bar chart (genomes-with-marker, one bar per marker, plus a horizontal line at `min_fraction × n_genomes`).

## 7. Backward compatibility

- All existing single-protein and whole-genome configs work unchanged.
- The `region: concatenated` value triggers the new path; no default behavior changes for families not in `CONCATENATION_FAMILIES`.
- For families newly added to `CONCATENATION_FAMILIES`, users with a pre-existing custom family yaml (which already specifies `region`) keep their setting — auto-defaults only apply to auto-generated configs (per the existing project policy: no migration of existing yamls).

## 8. v2 plan (HMMER swap-in, separate effort)

- Add `HmmerIdentifier` implementing the `MarkerIdentifier` protocol.
- Add HMMER as a runtime dependency (`hmmsearch`).
- Bundle or fetch HMM profiles per family (proposed: a per-family asset directory `vfam_trees/markers/<family>/<marker>.hmm`, or pull from a curated repository).
- Fetcher path changes: instead of per-marker name queries, fetch the full proteome of each genome (one query per source nuc accession), run `hmmsearch` with each marker's HMM, take the best hit.
- Keep `NameMatchIdentifier` as a fallback path or for families without curated HMMs.

## 9. Implementation phases

| Phase | Work | Tests |
|---|---|---|
| 1 | Config schema + curated marker presets in `config.py`; auto-generation logic | Schema validation; default-marker lookup |
| 2 | Marker-identification interface + `NameMatchIdentifier` | Tiebreaker logic, alias matching, length_range / locus_tag_hint |
| 3 | New fetcher path (per-marker queries, source-nuc grouping, paralog tiebreaker, drop counters) | Mock Entrez responses; multi-accession-rejection; min_fraction-rejection |
| 4 | Per-marker QC + per-marker length-outlier (RefSeq exempt) + concat alignment + partition file generation | Per-marker outlier counts; gap-padding correctness; partition-file format |
| 5 | Genome-level adaptive clustering + RefSeq absorption + cross-species merge | Cluster on concats; genome-level RefSeq prioritization |
| 6 | Tree inference: partitioned IQ-TREE for tree_100; partition-aware model parsing; branch-outlier at genome level | Partition file → IQ-TREE invocation; multi-model parsing; RefSeq genome exemption |
| 7 | Outputs: PhyloXML leaf adjustments, summary columns, per-family report page | PhyloXML schema validation; column population |
| 10 | Wire concat mode into the global sequence cache (`~/.vfam_cache`).  Extend the cache key with a stable hash of `concatenation.proteins` (names, aliases, length_range, locus_tag_hint) when `region == "concatenated"`.  Verify the per-marker `.gb` files round-trip through cache → fetch → reparse.  Remove the one-time "concat bypasses cache" WARNING from `_fetch_all_species` once landed. | Cache key includes marker-set hash; changing the marker list invalidates the cache for that family; existing single-protein cache entries unaffected |

Phases 8–9 (iterative branch-outlier loop in concat mode + per-family report PDF coverage page) shipped as v1.2.0; phase 10 is deferred.

Each phase ships as its own commit and stays runnable end-to-end.

## 10. Open items

Items requiring user decision before/during implementation:

1. ~~**tree_500 partitioning policy.**~~ **DECIDED**: option (a) — tree_500 stays FastTree, single-model on concat. Partitioned analysis only on tree_100 (IQ-TREE). Preserves the existing broad/fast vs. representative/rigorous split. See §5.10.

2. ~~**`min_fraction` default value.**~~ **DECIDED**: `0.7` (genome must have ≥70% of markers).

3. ~~**Per-partition model output format.**~~ **DECIDED**: single denormalized `tree100_marker_models` column, comma-separated `name:model` pairs. Per-marker columns are infeasible because different families have different marker sets within a single TSV schema.

4. **ASFV gene IDs (§4.3).** Deferred — user will verify against NC_001659 before the curation lands. Until then the four `(TBV)` entries in §4.3 stay as `null` `locus_tag_hint` values; absorption falls back to name+alias matching only.

5. ~~**Marker-fetch retry policy.**~~ **DECIDED**: no retry in MVP. Markers not matched by name+alias are treated as missing for that genome; the gap is logged and `min_fraction` decides whether the genome stays. BLAST/HMM-based recovery is left for v2 (the HMMER swap-in in §8 is the right tool for this job).

6. **Concat-mode cache keying.** Decision: when `region == "concatenated"`, the cache key auto-extends to include a stable hash of the full `concatenation.proteins` spec (names, aliases, length_range, locus_tag_hint). Same flat cache layout, no user-visible change. **STATUS: deferred — not yet implemented.** Phases 1–9 shipped concat mode without integrating the global sequence cache (`~/.vfam_cache`); every concat run currently hits NCBI for every species. `pipeline_concat._fetch_all_species` emits a one-time WARNING noting this. Tracked as **phase 10** below.

7. ~~**Entomopoxvirinae handling.**~~ **DECIDED**: B1 — single shared marker set, with optional **subfamily-aware aliases**. The marker spec gains an optional `aliases_<subfamily>` override (e.g. `aliases_Entomopoxvirinae`); at fetch time the pipeline determines each species' subfamily from its NCBI ranked lineage and unions the base `aliases` with any applicable subfamily-specific aliases for that species. One concat alignment, one tree pair per family — Chordopox and Entomopox sit in the same Poxviridae tree. The mechanism is generic and can be applied to any family later if needed; currently only Poxviridae uses it.

## 11. Non-goals (MVP)

- HMMER-based marker identification (v2).
- Coalescent / multi-tree concordance methods (ASTRAL).
- Per-marker single-protein trees as a side output.
- Recombination detection across markers within concat.
- User-supplied per-family HMM profiles (v2 territory).
- Migration of existing user yamls (preserved per project policy).
