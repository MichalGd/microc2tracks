#!/usr/bin/env python3
"""Generate a locked, content-keyed restriction-fragment cache for an exact FASTA."""
from __future__ import annotations

import argparse
import json
import re
from array import array
from pathlib import Path

from artifact_state import FileLock
from reference_registry import read_sizes
from workflow_common import json_write, registry_data, sha256, signature

GENERATOR_VERSION = 1


def scan_fasta(fasta: Path, definitions: dict, chunk_size: int = 1024 * 1024, check=None):
    """Scan with a motif-length overlap; do not join an entire genome in memory."""
    longest = max(len(x["motif"]) for x in definitions.values())
    patterns = [(re.compile(f"(?={value['motif']})"), value["cut_offset"]) for value in definitions.values()]
    chrom, length, pending, tail, cuts = None, 0, "", "", set()

    def scan(block, final=False):
        nonlocal length, tail
        if check:
            check()
        sequence = tail + block.upper()
        origin = length - len(tail)
        for pattern, offset in patterns:
            motif_length = len(pattern.pattern) - 4
            for match in pattern.finditer(sequence):
                if match.start() + motif_length > len(tail):
                    cuts.add(origin + match.start() + offset)
        length += len(block)
        tail = sequence[-(longest - 1):] if longest > 1 else ""

    with fasta.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom is not None:
                    scan(pending, True)
                    yield chrom, length, sorted({0, length, *cuts})
                chrom, length, pending, tail, cuts = line[1:].split()[0], 0, "", "", set()
            else:
                if chrom is None:
                    raise ValueError("FASTA sequence occurs before a header")
                sequence = line.strip()
                if not re.fullmatch("[ACGTRYMKSWBDHVNacgtrymkswbdhvn]*", sequence):
                    raise ValueError(f"Invalid FASTA sequence for {chrom}")
                pending += sequence
                if len(pending) >= chunk_size:
                    scan(pending)
                    pending = ""
    if chrom is not None:
        scan(pending, True)
        yield chrom, length, sorted({0, length, *cuts})


def prepare(fasta: Path, chrom_sizes: Path, definitions: dict, cache: Path, fasta_hash: str | None = None, check=None) -> tuple[Path, dict]:
    identity = {"fasta_sha256": fasta_hash or sha256(fasta), "chrom_sizes_sha256": sha256(chrom_sizes),
                "enzymes": definitions, "generator_version": GENERATOR_VERSION,
                "coordinates": "BED half-open; plus-strand cleavage offset; pairs pos-1"}
    key = signature(identity)
    cache.mkdir(parents=True, exist_ok=True)
    bed, sidecar = cache / f"{key}.fragments.bed", cache / f"{key}.json"
    with FileLock(cache / f"{key}.lock", blocking=True):
        if bed.is_file() and sidecar.is_file():
            saved = json.loads(sidecar.read_text())
            if saved.get("identity") == identity and saved.get("bed_sha256") == sha256(bed):
                return bed, saved
        expected = read_sizes(chrom_sizes, "chromosome sizes")
        seen = set()
        tmp = bed.with_suffix(".partial.bed")
        with tmp.open("w", encoding="utf-8", newline="\n") as out:
            for chrom, length, cuts in scan_fasta(fasta, definitions, check=check):
                if chrom in seen or expected.get(chrom) != length:
                    raise ValueError(f"FASTA/digest length or duplicate chromosome inconsistency: {chrom}")
                seen.add(chrom)
                for index, (start, end) in enumerate(zip(cuts, cuts[1:])):
                    if check and index % 100000 == 0:
                        check()
                    if start < end:
                        out.write(f"{chrom}\t{start}\t{end}\t{index}\n")
        if seen != set(expected):
            raise ValueError("FASTA and chromosome sizes contain different chromosome names")
        tmp.replace(bed)
        result = {"identity": identity, "digest_id": key, "bed_sha256": sha256(bed)}
        json_write(sidecar, result)
        return bed, result


def load_digest(path: Path) -> dict[str, array]:
    result: dict[str, array] = {}
    with path.open() as handle:
        for line in handle:
            chrom, start, end, ordinal = line.split()
            values = result.setdefault(chrom, array("q", [0]))
            if int(start) != values[-1] or int(ordinal) != len(values) - 1 or int(end) <= int(start):
                raise ValueError(f"Digest contains gaps, overlaps or invalid fragment ordinals: {path}")
            values.append(int(end))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--enzymes", required=True, help="Explicit '+'-separated enzyme names")
    parser.add_argument("--registry", type=Path, default=Path(__file__).resolve().parent.parent / "config/hic_protocols.json")
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    data = registry_data(args.registry)
    definitions = {name: data["enzymes"][name] for name in args.enzymes.split("+")}
    print(prepare(args.fasta, args.chrom_sizes, definitions, args.cache)[0])


if __name__ == "__main__":
    main()
