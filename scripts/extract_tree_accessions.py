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

Output: ``<out_dir>/<family>_tree_<label>_tips.tsv`` per family, or
``..._ids.txt`` with ``--ids-only``.  Nothing is written into the results
directory itself unless you point ``-o`` at it.

Usage (from anywhere):
    # one TSV per family, all the tree_500 trees of the run
    python scripts/extract_tree_accessions.py <results_dir> -o tips_500/
    # just the accessions, one file per family, ready for efetch
    python scripts/extract_tree_accessions.py <results_dir> -o ids_500/ --ids-only
    # the tree_100 trees instead, or both (two files per family)
    python scripts/extract_tree_accessions.py <results_dir> -o tips/ --tree 100
    python scripts/extract_tree_accessions.py <results_dir> -o tips/ --tree both

Writes a per-family count and the totals to stderr, so a family whose tree is
missing or unreadable cannot pass unnoticed.  Standard library only.
"""
from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

#: Columns of each per-family TSV, in order.  ``family`` and ``tree`` are
#: repeated in every row although the file name carries them, so the files can
#: be concatenated without losing track of where a row came from.
COLUMNS = ("family", "tree", "accession", "source", "seq_type", "taxid",
           "species", "label")

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


def write_tsv(path: Path, rows: list[dict[str, str]]) -> int:
    with path.open("w") as out:
        print("\t".join(COLUMNS), file=out)
        for row in rows:
            print("\t".join(row[c] for c in COLUMNS), file=out)
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
    args = ap.parse_args(argv)

    labels = ("500", "100") if args.tree == "both" else (args.tree,)
    files = sorted(f for label in labels
                   for f in args.results_dir.glob(f"*/*_tree_{label}.xml"))
    if not files:
        ap.error(f"no *_tree_{'/'.join(labels)}.xml under {args.results_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sources: Counter[str] = Counter()
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
            elif args.strip_version:
                row["accession"] = row["accession"].rsplit(".", 1)[0]
            row["family"], row["tree"] = family, label
        suffix = "_ids.txt" if args.ids_only else "_tips.tsv"
        dest = args.out_dir / f"{family}_tree_{label}{suffix}"
        written = write_ids(dest, rows) if args.ids_only else write_tsv(dest, rows)
        n_tips += len(rows)
        n_written += written
        print(f"  {dest.name}: {written} of {len(rows)} tips", file=sys.stderr)

    print(f"{len(files) - n_bad} tree file(s), {n_tips} tips, "
          f"{n_written} line(s) in {len(files) - n_bad} file(s) under {args.out_dir}",
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
