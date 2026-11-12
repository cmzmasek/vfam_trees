#!/usr/bin/env python3
"""Extract the NCBI identifiers of the external nodes (tips) of the trees in a
results directory — one TSV per family.

Every tip the pipeline writes carries its identifier in the phyloXML
``<sequence><accession source="...">`` element — the protein accession in
marker mode, the nucleotide accession in whole-genome mode — so that is what
this reads, not the leaf label (which is a ``labeling.format`` string and can
be changed per run).  ``source`` is ``ncbi`` for fetched records; sequences
grafted in from elsewhere (``additional_data_sources``, ``manual.include_seq``)
carry ``pathoplexus`` or ``user`` and are reported but left out of
``--ids-only``.

``--source-accession`` adds the accession of the record the tip's sequence was
read off — the isolate's nucleotide record.  That one is **not** in the tree
file, and not in the metadata TSV either: in marker mode the tip is a protein
accession and only the protein record itself says which genome it is annotated
on.  So it is resolved against NCBI with the pipeline's own
``fetch._source_nuc_accession`` (the GenPept ``coded_by`` qualifier, then the
RefSeq ``DBSOURCE`` fallback), and cached in a TSV so a second run is free.
Whole-genome-nucleotide tips are already the isolate's record, so they get
their own accession back.

Output: ``<out_dir>/<family>_tree_<label>_tips.tsv`` per family, or
``..._ids.txt`` with ``--ids-only``.  Nothing is written into the results
directory itself unless you point ``-o`` at it.

Usage (from anywhere):
    # one TSV per family, all the tree_500 trees of the run
    python scripts/extract_tree_accessions.py <results_dir> -o tips_500/
    # ... with the source isolate's accession as well (queries NCBI)
    python scripts/extract_tree_accessions.py <results_dir> -o tips_500/ \
        --source-accession
    # just the accessions, one file per family, ready for efetch
    python scripts/extract_tree_accessions.py <results_dir> -o ids_500/ --ids-only
    # the tree_100 trees instead, or both (two files per family)
    python scripts/extract_tree_accessions.py <results_dir> -o tips/ --tree 100
    python scripts/extract_tree_accessions.py <results_dir> -o tips/ --tree both

Writes a per-family count and the totals to stderr, so a family whose tree is
missing or unreadable cannot pass unnoticed.  Standard library only, unless
``--source-accession`` is given: that path needs Biopython and an installed
``vfam_trees``.
"""
from __future__ import annotations

import argparse
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

#: Columns of each per-family TSV, in order.  ``family`` and ``tree`` are
#: repeated in every row although the file name carries them, so the files can
#: be concatenated without losing track of where a row came from.
COLUMNS = ("family", "tree", "accession", "source", "seq_type", "taxid",
           "species", "label")

#: Extra column added by ``--source-accession``.
SOURCE_COLUMN = "source_accession"

#: Protein accessions per efetch call.  NCBI accepts a few hundred; GenPept
#: records are a few kB each, so this keeps a batch around half a megabyte.
BATCH_SIZE = 200

#: Retries per batch, and the pause between them, matching ``fetch.py``.
MAX_TRIES = 3
RETRY_DELAY = 10

#: Leaf rank/annotation properties are written as ``vipr:<Field>`` in any
#: namespace prefix, so match on the part after the colon.
SPECIES_PROPERTY = "Species"


def local(tag: str) -> str:
    """Local name of an element tag, namespace or not."""
    return tag.rsplit("}", 1)[-1]


def children(el: ET.Element, name: str) -> list[ET.Element]:
    return [c for c in el if local(c.tag) == name]


def child(el: ET.Element, name: str) -> ET.Element | None:
    found = children(el, name)
    return found[0] if found else None


def text(el: ET.Element | None) -> str:
    return (el.text or "").strip() if el is not None else ""


def species_of(clade: ET.Element) -> str:
    """ICTV species of a tip: its ``vipr:Species`` property."""
    for prop in children(clade, "property"):
        if prop.get("ref", "").rsplit(":", 1)[-1] == SPECIES_PROPERTY:
            return text(prop)
    return ""


def tips(path: Path) -> list[dict[str, str]]:
    """One record per external node of a phyloXML file."""
    root = ET.parse(path).getroot()
    rows = []
    for clade in root.iter():
        if local(clade.tag) != "clade" or children(clade, "clade"):
            continue  # internal node, or not a clade at all
        seq = child(clade, "sequence")
        acc = child(seq, "accession") if seq is not None else None
        taxonomy = child(clade, "taxonomy")
        rows.append({
            "accession": text(acc),
            "source": acc.get("source", "") if acc is not None else "",
            "seq_type": seq.get("type", "") if seq is not None else "",
            "taxid": text(child(taxonomy, "id")) if taxonomy is not None else "",
            "species": species_of(clade),
            "label": text(child(clade, "name")),
        })
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]],
              columns: tuple[str, ...] = COLUMNS) -> int:
    with path.open("w") as out:
        print("\t".join(columns), file=out)
        for row in rows:
            print("\t".join(row.get(c, "") for c in columns), file=out)
    return len(rows)


def write_ids(path: Path, rows: list[dict[str, str]]) -> int:
    seen: set[str] = set()
    with path.open("w") as out:
        for row in rows:
            acc = row["accession"]
            if row["source"] == "ncbi" and acc and acc not in seen:
                seen.add(acc)
                print(acc, file=out)
    return len(seen)


def load_source_cache(path: Path) -> dict[str, str]:
    """``accession -> source accession`` written by an earlier run."""
    known: dict[str, str] = {}
    if not path.exists():
        return known
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0] != "accession":
                known[parts[0]] = parts[1]
    return known


def configure_ncbi(global_config: Path | None, email: str | None,
                   api_key: str | None) -> None:
    """Set the NCBI contact details, from the flags or from a global.yaml."""
    from vfam_trees import fetch

    if not email and global_config and global_config.exists():
        import yaml
        ncbi = (yaml.safe_load(global_config.read_text()) or {}).get("ncbi", {})
        email, api_key = ncbi.get("email"), api_key or ncbi.get("api_key")
    if not email:
        raise SystemExit(
            "--source-accession needs an NCBI email: pass --email, or "
            "--global-config pointing at a global.yaml that has one")
    fetch.configure_entrez(email, api_key)


def in_batches(ids: list[str], size: int, query, what: str) -> dict[str, str]:
    """Run ``query`` over ``ids`` in batches, mapping id -> answer.

    A batch that keeps failing is retried id by id: NCBI answers a whole
    batch with HTTP 400 if a single id is not in the database being asked,
    so one odd accession would otherwise cost its 200 neighbours.
    """
    found: dict[str, str] = {}
    for start in range(0, len(ids), size):
        batch = ids[start:start + size]
        for attempt in range(1, MAX_TRIES + 1):
            try:
                found.update(query(batch))
                break
            except Exception as exc:                      # network, HTTP, parse
                if attempt < MAX_TRIES:
                    print(f"    {what} at {start}: {exc} — retrying in "
                          f"{RETRY_DELAY}s", file=sys.stderr)
                    time.sleep(RETRY_DELAY)
                    continue
                print(f"    {what} at {start}: {exc} — falling back to single "
                      f"queries", file=sys.stderr)
                for acc in batch:
                    try:
                        found.update(query([acc]))
                    except Exception as exc2:
                        print(f"      {acc}: {exc2}", file=sys.stderr)
        print(f"    {what}: {len(found)}/{len(ids)} answered", file=sys.stderr)
    return found


def fetch_source_accessions(accessions: list[str], batch_size: int) -> dict[str, str]:
    """Read the source nucleotide accession off each protein record at NCBI.

    GenPept records, parsed with the pipeline's own
    ``fetch._source_nuc_accession``, so the answer is the one the pipeline
    itself would have recorded.
    """
    from vfam_trees import fetch
    from Bio import Entrez, SeqIO

    def query(ids: list[str]) -> dict[str, str]:
        handle = Entrez.efetch(db="protein", id=",".join(ids),
                               rettype="gp", retmode="text")
        try:
            return {rec.id: fetch._source_nuc_accession(rec)
                    for rec in SeqIO.parse(handle, "genbank")}
        finally:
            handle.close()

    return in_batches(accessions, batch_size, query, "protein records")


def confirm_nucleotide(accessions: list[str], batch_size: int) -> dict[str, str]:
    """Keep the accessions NCBI really serves from nuccore, id -> accession.

    Accession shape alone cannot tell a two-letter nucleotide accession from
    a two-letter protein one, and a tip that is its own source has to be the
    former, so ask instead of guessing.  ``rettype="acc"`` answers with one
    accession per line and nothing else.
    """
    from Bio import Entrez

    def query(ids: list[str]) -> dict[str, str]:
        handle = Entrez.efetch(db="nuccore", id=",".join(ids),
                               rettype="acc", retmode="text")
        try:
            served = [line.strip() for line in handle if line.strip()]
        finally:
            handle.close()
        if len(served) != len(ids):      # cannot be matched up id by id
            raise RuntimeError(f"asked for {len(ids)} accessions, got {len(served)}")
        return dict(zip(ids, served))

    return in_batches(accessions, batch_size, query, "nuccore check")


def resolve_sources(rows: list[dict[str, str]], cache_path: Path,
                    batch_size: int) -> dict[str, str]:
    """Map every tip accession to the accession of the record it came from.

    A whole-genome-nucleotide tip already *is* that record.  So is a tip the
    translation fallback renamed to its genome accession — those look like
    protein tips but carry a nucleotide accession, which is why anything of
    nucleotide shape is confirmed against nuccore before being mapped to
    itself, and kept out of the protein batches either way.  Everything else
    is a protein accession and is looked up, unless the cache has it already.
    """
    from vfam_trees import fetch

    known = load_source_cache(cache_path)
    proteins: list[str] = []
    maybe_own: list[str] = []
    for row in rows:
        acc = row["accession"]
        if not acc or acc in known:
            continue
        if row["seq_type"] == "dna":
            known[acc] = acc                       # the tip is the isolate's record
        elif fetch.accession_db(acc) == "nuccore":
            if acc not in maybe_own:
                maybe_own.append(acc)
        elif acc not in proteins:
            proteins.append(acc)

    resolved: dict[str, str] = {}
    if maybe_own:
        print(f"  checking {len(maybe_own)} nucleotide-shaped tip accession(s) "
              f"against nuccore", file=sys.stderr)
        served = confirm_nucleotide(maybe_own, batch_size)
        resolved.update({acc: acc for acc in served})
        for acc in maybe_own:
            if acc not in served:
                print(f"    {acc}: not served by nuccore — left empty",
                      file=sys.stderr)
    if proteins:
        print(f"  resolving {len(proteins)} protein accession(s) at NCBI "
              f"({len(known)} known from the cache or the tree)", file=sys.stderr)
        resolved.update(fetch_source_accessions(proteins, batch_size))

    if resolved:
        is_new_file = not cache_path.exists()
        with cache_path.open("a") as out:
            if is_new_file:
                print("accession\tsource_accession", file=out)
            for acc, src in resolved.items():
                print(f"{acc}\t{src}", file=out)
        known.update(resolved)
    return known


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir", type=Path,
                    help="a results directory holding one subdirectory per family")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("."),
                    help="directory for the per-family files (created if needed; "
                         "default: the current directory)")
    ap.add_argument("--tree", choices=("500", "100", "both"), default="500",
                    help="which tree of each family to read (default: 500)")
    ap.add_argument("--ids-only", action="store_true",
                    help="write one NCBI accession per line, deduplicated within "
                         "the family, instead of the TSV")
    ap.add_argument("--strip-version", action="store_true",
                    help="drop the trailing .N version from each accession")
    ap.add_argument("--source-accession", action="store_true",
                    help=f"add a {SOURCE_COLUMN} column: the accession of the "
                         f"record the tip's sequence was read off. Queries NCBI "
                         f"for the protein tips and caches the answers")
    ap.add_argument("--source-cache", type=Path,
                    help="where to keep the resolved source accessions "
                         "(default: <out_dir>/_source_accessions.tsv)")
    ap.add_argument("--global-config", type=Path,
                    help="a global.yaml to take the NCBI email and API key from "
                         "(default: <results_dir>/../config/global.yaml)")
    ap.add_argument("--email", help="NCBI contact email, overriding the global.yaml")
    ap.add_argument("--api-key", help="NCBI API key, overriding the global.yaml")
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE,
                    help=f"accessions per NCBI query (default: {BATCH_SIZE})")
    args = ap.parse_args(argv)
    if args.ids_only and args.source_accession:
        ap.error("--ids-only writes bare accessions; it has no columns to add")

    labels = ("500", "100") if args.tree == "both" else (args.tree,)
    files = sorted(f for label in labels
                   for f in args.results_dir.glob(f"*/*_tree_{label}.xml"))
    if not files:
        ap.error(f"no *_tree_{'/'.join(labels)}.xml under {args.results_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Parse every tree first, so --source-accession can batch its NCBI
    # queries over the whole run instead of per family.
    sources: Counter[str] = Counter()
    per_file: list[tuple[Path, list[dict[str, str]]]] = []
    n_tips = n_written = n_bad = n_missing = 0
    for path in files:
        label = path.stem.rsplit("_tree_", 1)[1]
        family = path.name[: path.name.index("_tree_")]
        try:
            rows = tips(path)
        except ET.ParseError as exc:
            print(f"  {path}: UNREADABLE ({exc})", file=sys.stderr)
            n_bad += 1
            continue
        for row in rows:
            sources[row["source"] or "(none)"] += 1
            if not row["accession"]:
                n_missing += 1
            row["family"], row["tree"] = family, label
        suffix = "_ids.txt" if args.ids_only else "_tips.tsv"
        per_file.append((args.out_dir / f"{family}_tree_{label}{suffix}", rows))
        n_tips += len(rows)

    columns = COLUMNS
    if args.source_accession:
        columns = COLUMNS + (SOURCE_COLUMN,)
        default_config = args.results_dir.parent / "config" / "global.yaml"
        configure_ncbi(args.global_config or default_config,
                       args.email, args.api_key)
        cache = args.source_cache or args.out_dir / "_source_accessions.tsv"
        all_rows = [row for _, rows in per_file for row in rows]
        known = resolve_sources(all_rows, cache, args.batch_size)
        unresolved = 0
        for row in all_rows:
            row[SOURCE_COLUMN] = known.get(row["accession"], "")
            if not row[SOURCE_COLUMN]:
                unresolved += 1
        print(f"  source accessions: {len(all_rows) - unresolved} of "
              f"{len(all_rows)} tips, cached in {cache}", file=sys.stderr)
        if unresolved:
            print(f"  WARNING: {unresolved} tip(s) have no source accession",
                  file=sys.stderr)

    # Version stripping happens last, so the NCBI lookups above use the
    # versioned accessions NCBI expects.
    for dest, rows in per_file:
        if args.strip_version:
            for row in rows:
                for key in ("accession", SOURCE_COLUMN):
                    if row.get(key):
                        row[key] = row[key].rsplit(".", 1)[0]
        written = (write_ids(dest, rows) if args.ids_only
                   else write_tsv(dest, rows, columns))
        n_written += written
        print(f"  {dest.name}: {written} of {len(rows)} tips", file=sys.stderr)

    print(f"{len(files) - n_bad} tree file(s), {n_tips} tips, "
          f"{n_written} line(s) in {len(per_file)} file(s) under {args.out_dir}",
          file=sys.stderr)
    print("  identifier sources: "
          + ", ".join(f"{k} {v}" for k, v in sorted(sources.items())), file=sys.stderr)
    if n_missing:
        print(f"  WARNING: {n_missing} tip(s) carry no accession", file=sys.stderr)
    if n_bad:
        print(f"  WARNING: {n_bad} file(s) could not be parsed", file=sys.stderr)
    return 1 if (n_bad or n_missing) else 0


if __name__ == "__main__":
    sys.exit(main())
