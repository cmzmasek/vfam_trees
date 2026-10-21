"""MSA inference using MAFFT (or other configured tool)."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .logger import get_logger

log = get_logger(__name__)


def get_mafft_version() -> str:
    """Return the MAFFT version string, or 'unknown' if it cannot be determined."""
    try:
        result = subprocess.run(
            ["mafft", "--version"], capture_output=True, text=True
        )
        # MAFFT prints version to stderr: "v7.520 (2023/Apr/16)"
        output = (result.stderr or result.stdout).strip()
        m = re.search(r"v?(\d+\.\d+\S*)", output)
        return m.group(1) if m else output.split()[0] if output else "unknown"
    except Exception:
        return "unknown"


def validate_msa(fasta_path: Path) -> None:
    """Validate that a FASTA MSA file is non-empty and all sequences are the same length.

    Raises:
        FileNotFoundError: if the file does not exist or is empty
        ValueError: if fewer than 2 sequences, or sequences differ in length
    """
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        raise FileNotFoundError(f"MSA output is missing or empty: {fasta_path}")
    from Bio import SeqIO
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) < 2:
        raise ValueError(f"MSA contains fewer than 2 sequences: {fasta_path}")
    lengths = {len(r.seq) for r in records}
    if len(lengths) > 1:
        raise ValueError(
            f"MSA sequences have inconsistent lengths {lengths} — "
            f"alignment may be corrupted: {fasta_path}"
        )
    log.debug("MSA validated: %d sequences, alignment length %d", len(records), lengths.pop())


# Prefix MAFFT --adjustdirection puts on the name of every sequence it
# reverse-complemented.
MAFFT_REVERSED_PREFIX = "_R_"

# Sidecar beside an alignment listing the ids it holds reverse-complemented,
# one per line, so a run that reuses a cached alignment can still report them.
FLIPPED_IDS_SUFFIX = ".flipped"

# IUPAC nucleotide complement, both cases; gaps and unknown symbols pass through.
_COMPLEMENT = str.maketrans(
    "ACGTURYKMBVDHNSWacgturykmbvdhnsw",
    "TGCAAYRMKVBHDNSWtgcaayrmkvbhdnsw",
)


def reverse_complement(seq: str) -> str:
    """Reverse-complement a nucleotide (alignment) row; gaps are kept."""
    return seq.translate(_COMPLEMENT)[::-1]


def _read_records(path: Path) -> list:
    from Bio import SeqIO
    return list(SeqIO.parse(str(path), "fasta"))


def _write_records(records: list, path: Path) -> None:
    from Bio import SeqIO
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fasta")


def strip_mafft_direction_prefix(fasta_path: Path) -> list[str]:
    """Remove MAFFT's ``_R_`` name prefix from *fasta_path* in place.

    ``--adjustdirection`` reverse-complements sequences deposited on the
    opposite strand and marks them by prefixing the name.  The sequences are
    left as MAFFT wrote them (every row on the same strand); only the names
    are restored so downstream id lookups keep working.

    Returns the ids of the flipped sequences, in file order; empty when
    nothing was flipped (the file is then left untouched).
    """
    records = _read_records(fasta_path)
    flipped: list[str] = []
    for rec in records:
        if rec.id.startswith(MAFFT_REVERSED_PREFIX):
            n = len(MAFFT_REVERSED_PREFIX)
            rec.id = rec.name = rec.id[n:]
            rec.description = rec.description[n:]
            flipped.append(rec.id)
    if flipped:
        _write_records(records, fasta_path)
    return flipped


def _flipped_ids_path(alignment_fasta: Path) -> Path:
    return alignment_fasta.with_name(alignment_fasta.name + FLIPPED_IDS_SUFFIX)


def write_flipped_ids(alignment_fasta: Path, ids: list[str]) -> None:
    """Record *ids* beside *alignment_fasta*; remove the record when empty."""
    path = _flipped_ids_path(alignment_fasta)
    if ids:
        path.write_text("".join(f"{i}\n" for i in ids))
    else:
        path.unlink(missing_ok=True)


def read_flipped_ids(alignment_fasta: Path) -> list[str]:
    """Return the ids recorded by :func:`write_flipped_ids`; empty if none."""
    try:
        text = _flipped_ids_path(alignment_fasta).read_text()
    except OSError:
        return []
    return [line.strip() for line in text.splitlines() if line.strip()]


def orient_to_majority_strand(alignment_fasta: Path, flipped: list[str]) -> list[str]:
    """Put the alignment on the strand most records were deposited on.

    MAFFT's ``--adjustdirection`` orients every row against the *first*
    input record.  When that record is on the minority strand, MAFFT flips
    the majority instead.  Every row is then on one strand already, so
    reverse-complementing the whole alignment puts it on the majority strand
    without a second MAFFT run; the tree is unaffected.  At an exact tie the
    strand of the record with the smallest id wins, so the choice does not
    depend on record order either way.

    Returns the ids now held reverse-complemented (in file order), and
    rewrites *alignment_fasta* only when the strand changes.
    """
    if not flipped:
        return []
    records = _read_records(alignment_fasta)
    ids = [r.id for r in records]
    flipped_set = set(flipped)
    n_flipped, n = len(flipped_set), len(records)
    flip_all = 2 * n_flipped > n or (2 * n_flipped == n and min(ids) in flipped_set)
    if not flip_all:
        return flipped
    from Bio.Seq import Seq
    for rec in records:
        rec.seq = Seq(reverse_complement(str(rec.seq)))
    _write_records(records, alignment_fasta)
    now_flipped = [i for i in ids if i not in flipped_set]
    log.info(
        "MAFFT --adjustdirection oriented %d of %d sequences against a "
        "minority-strand leader; reverse-complemented the alignment so that "
        "%d stay as deposited.",
        n_flipped, n, n - len(now_flipped),
    )
    return now_flipped


def run_msa(
    input_fasta: Path,
    output_fasta: Path,
    tool: str = "mafft",
    options: str = "--auto",
    threads: int = 1,
) -> list[str]:
    """Run multiple sequence alignment.

    Args:
        input_fasta: unaligned sequences (short IDs)
        output_fasta: output aligned FASTA
        tool: alignment tool name ('mafft' supported)
        options: tool-specific options string
        threads: number of CPU threads

    Returns:
        ids of sequences the alignment holds reverse-complemented (MAFFT
        ``--adjustdirection``); empty otherwise.  The same list is recorded
        beside the alignment (:func:`read_flipped_ids`).

    With ``--adjustdirection`` the alignment is put on the majority strand
    (:func:`orient_to_majority_strand`), so the strand never depends on
    which record happens to come first.
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if tool.lower() != "mafft":
        raise ValueError(f"Unsupported MSA tool: {tool}. Supported: mafft")

    flipped = _run_mafft(input_fasta, output_fasta, options, threads)
    flipped = orient_to_majority_strand(output_fasta, flipped)
    write_flipped_ids(output_fasta, flipped)

    log.debug("MSA complete: %s", output_fasta)
    return flipped


def _run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    options: str,
    threads: int,
) -> list[str]:
    cmd = ["mafft"]
    cmd += options.split()
    adjusts_direction = any(o.startswith("--adjustdirection") for o in options.split())
    cmd += ["--thread", str(threads)]
    cmd += ["--out", str(output_fasta)]
    cmd += [str(input_fasta)]

    log.info("Running MAFFT (%s, thread=%d)", options, threads)
    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("MAFFT stderr:\n%s", (result.stderr or "")[-4000:])
        log.error("MAFFT stdout:\n%s", (result.stdout or "")[-1000:])
        raise RuntimeError(
            f"MAFFT failed with exit code {result.returncode}. If MAFFT was "
            "upgraded recently, its CLI flags may have changed — see stderr "
            "above."
        )

    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        log.error("MAFFT stderr:\n%s", (result.stderr or "")[-4000:])
        raise FileNotFoundError(
            f"MAFFT exited 0 but produced no alignment at {output_fasta}. "
            "This may indicate a change in MAFFT's output behavior — see "
            "stderr above."
        )

    log.debug("MAFFT stdout: %s", result.stdout[:500] if result.stdout else "(none)")
    log.info("MAFFT complete")
    flipped = strip_mafft_direction_prefix(output_fasta) if adjusts_direction else []
    if flipped:
        log.debug(
            "MAFFT --adjustdirection reversed %d sequence(s) relative to the first record",
            len(flipped),
        )
    return flipped
