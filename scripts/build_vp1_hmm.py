#!/usr/bin/env python3
"""Build vfam_trees/data/hmms/VP1/VP1_entero.hmm — the enterovirus VP1 profile.

Why this one profile is not a Pfam download (the other bundled sets are, see
``build_bundled_hmms.sh``): Pfam cannot single VP1 out.  PF00073 (Rhv) matches
VP2, VP3 *and* VP1 as three separate domains of the polyprotein, and Pfam
assigns VP1 itself to PF00073 or to PF22663 (Rhv_5) depending on the virus
(InterPro matches of UniProt P03300, P03303, Q66478, checked 2026-11-05), so a
Pfam token would crop the whole capsid block instead of the typing region.

Recipe: take enterovirus genome deposits that annotate a VP1 ``mat_peptide``,
keep one per NCBI organism (RefSeq preferred, then lowest accession — so the
seed set spans the genus rather than the sequencing effort), align with MAFFT
L-INS-i, and ``hmmbuild`` the profile under the NAME ``VP1_entero`` (the token
in ``config.VP1_HMM_TOKENS``).

The committed profile was built on 2026-11-05 from the cached genomes of the
2026-10-31 Enterovirus run (``--from-gb``): 242 seeds of 242 organisms (19
RefSeq), model length 282.  Validation, on the 7,431 polyproteins of that run:
a passing VP1 hit in 7,425 (99.9 %), model coverage median 0.98, domain
E-value median 5e-108; VP2/VP3 never pass, so the crop stays on VP1.

Usage (from the repo root):
    # from local GenBank files (how the committed profile was made)
    python scripts/build_vp1_hmm.py --from-gb '<run>/_work/species/sp_*/sequences.gb'
    # or fetch the genomes from NCBI
    python scripts/build_vp1_hmm.py --email you@example.org [--api-key KEY]

Needs mafft and hmmbuild in PATH.  GenBank records are public domain.
"""
from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import Entrez, SeqIO

#: VP1 goes by its name or its L-numbering (1D) in mat_peptide products.
VP1_PRODUCT = re.compile(r"\bVP1\b|\b1D\b", re.IGNORECASE)
#: Plausible VP1 length (aa); anything else is a mis-annotation.
VP1_LENGTH = (150, 400)
ENTEROVIRUS_TAXID = 12059
PROFILE_NAME = "VP1_entero"
OUT_PATH = Path("vfam_trees/data/hmms/VP1") / f"{PROFILE_NAME}.hmm"


def vp1_peptides(records) -> list[dict]:
    """Every annotated VP1 mat_peptide translation in *records*."""
    out = []
    for rec in records:
        for ft in rec.features:
            if ft.type != "mat_peptide":
                continue
            if not any(VP1_PRODUCT.search(p) for p in ft.qualifiers.get("product", [])):
                continue
            seq = ft.qualifiers.get("translation", [None])[0]
            if seq is None:
                try:
                    seq = str(ft.extract(rec.seq).translate(to_stop=True))
                except Exception:          # sequence-less record (WGS master)
                    break
            seq = (seq or "").replace("*", "")
            if VP1_LENGTH[0] <= len(seq) <= VP1_LENGTH[1]:
                out.append({
                    "acc": rec.id,
                    "organism": rec.annotations.get("organism", "unknown"),
                    "refseq": rec.id.startswith("NC_"),
                    "seq": seq,
                })
            break
    return out


def from_genbank(patterns: list[str]) -> list[dict]:
    paths = sorted({p for pat in patterns for p in glob.glob(pat)})
    if not paths:
        sys.exit(f"No GenBank files matched: {patterns}")
    print(f"Reading {len(paths)} GenBank file(s) ...")
    peptides = []
    for path in paths:
        peptides += vp1_peptides(SeqIO.parse(path, "genbank"))
    return peptides


def from_ncbi(email: str, api_key: str | None, max_genomes: int) -> list[dict]:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    query = (
        f"txid{ENTEROVIRUS_TAXID}[Organism:exp]"
        ' AND ("complete genome"[Title] OR "complete sequence"[Title])'
        " NOT patent[filter]"
    )
    ids: list[str] = []
    for refseq_only in (True, False):     # RefSeq first, then top up
        term = query + (" AND refseq[filter]" if refseq_only else "")
        with Entrez.esearch(db="nuccore", term=term, retmax=max_genomes) as h:
            found = Entrez.read(h)["IdList"]
        ids += [i for i in found if i not in set(ids)]
        ids = ids[:max_genomes]
    print(f"Fetching {len(ids)} genome records from NCBI ...")
    peptides = []
    for start in range(0, len(ids), 200):
        batch = ids[start:start + 200]
        with Entrez.efetch(db="nuccore", id=",".join(batch), rettype="gb",
                           retmode="text") as h:
            peptides += vp1_peptides(SeqIO.parse(h, "genbank"))
        print(f"  {min(start + 200, len(ids))}/{len(ids)}")
    return peptides


def pick_seeds(peptides: list[dict]) -> list[dict]:
    """One VP1 per organism: RefSeq first, then the lowest accession."""
    by_organism: dict[str, list[dict]] = defaultdict(list)
    for p in peptides:
        by_organism[p["organism"]].append(p)
    seeds = []
    for organism in sorted(by_organism):
        candidates = sorted(by_organism[organism], key=lambda p: (not p["refseq"], p["acc"]))
        seeds.append(candidates[0])
    return seeds


def build(seeds: list[dict], out_path: Path, threads: int) -> None:
    lengths = [len(s["seq"]) for s in seeds]
    print(f"{len(seeds)} seed VP1s of {len(seeds)} organisms "
          f"({sum(s['refseq'] for s in seeds)} RefSeq), lengths "
          f"{min(lengths)}-{max(lengths)} aa")
    if len(seeds) < 20:
        sys.exit("Too few seeds for a genus-wide profile — check the input.")
    with tempfile.TemporaryDirectory() as tmp:
        faa, aln = Path(tmp) / "seed.faa", Path(tmp) / "seed.aln"
        with open(faa, "w") as f:
            for s in seeds:
                f.write(f">{s['acc']}|{s['organism'].replace(' ', '_')}\n{s['seq']}\n")
        # --threadit 0 keeps MAFFT's refinement deterministic (see CLAUDE.md
        # "Reproducibility"), so a rebuild gives the same alignment.
        subprocess.run(
            f"mafft --maxiterate 1000 --localpair --threadit 0 --thread {threads} "
            f"--quiet {faa} > {aln}", shell=True, check=True,
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["hmmbuild", "--amino", "-n", PROFILE_NAME, "--informat", "afa",
             str(out_path), str(aln)],
            check=True, capture_output=True,
        )
    head = [l for l in out_path.read_text().splitlines()[:20]
            if l.startswith(("NAME", "LENG", "NSEQ"))]
    print(f"Wrote {out_path}\n  " + "\n  ".join(head))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-gb", nargs="+", metavar="GLOB",
                    help="build from local GenBank file(s) instead of fetching")
    ap.add_argument("--email", help="NCBI contact e-mail (required to fetch)")
    ap.add_argument("--api-key", default=None, help="NCBI API key")
    ap.add_argument("--max-genomes", type=int, default=3000,
                    help="genome records to fetch (default 3000)")
    ap.add_argument("--threads", type=int, default=8, help="MAFFT threads")
    ap.add_argument("-o", "--out", type=Path, default=OUT_PATH)
    args = ap.parse_args()

    if args.from_gb:
        peptides = from_genbank(args.from_gb)
    elif args.email:
        peptides = from_ncbi(args.email, args.api_key, args.max_genomes)
    else:
        ap.error("pass --from-gb or --email")
    print(f"{len(peptides)} annotated VP1 peptides, "
          f"{len({p['organism'] for p in peptides})} organisms")
    build(pick_seeds(peptides), args.out, args.threads)


if __name__ == "__main__":
    main()
