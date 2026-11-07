#!/usr/bin/env python3
"""Stream paired FASTQs and fail if read names or record counts differ."""

from __future__ import annotations

import argparse
import gzip
import sys
from contextlib import contextmanager, ExitStack
from pathlib import Path
from typing import Iterator, TextIO


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix.lower() in {".gz", ".gzip"}:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            yield handle
    else:
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            yield handle


def normalized_name(header: str) -> str:
    token = header.rstrip("\r\n").split(maxsplit=1)[0]
    if not token.startswith("@"):
        raise ValueError(f"invalid FASTQ header {token!r}")
    token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def records(handle: TextIO, label: str) -> Iterator[str]:
    record_number = 0
    while True:
        lines = [handle.readline() for _ in range(4)]
        if not lines[0]:
            if any(lines[1:]):
                raise ValueError(f"{label} has a truncated final FASTQ record")
            return
        record_number += 1
        if any(line == "" for line in lines[1:]):
            raise ValueError(f"{label} record {record_number} is truncated")
        if not lines[2].startswith("+"):
            raise ValueError(f"{label} record {record_number} has an invalid separator")
        yield normalized_name(lines[0])


def check_pairs(r1: Path, r2: Path) -> int:
    count = 0
    with open_text(r1) as first, open_text(r2) as second:
        names1 = records(first, "R1")
        names2 = records(second, "R2")
        while True:
            name1 = next(names1, None)
            name2 = next(names2, None)
            if name1 is None or name2 is None:
                if name1 != name2:
                    raise ValueError(f"paired FASTQs have different record counts after {count} synchronized pair(s)")
                return count
            count += 1
            if name1 != name2:
                raise ValueError(f"paired read-name mismatch at pair {count}: R1={name1!r}, R2={name2!r}")


def inspect_pairs(r1: Path, r2: Path, *, writers=None, check=lambda: None) -> dict:
    """Strict BL reader contract, optionally copying/compressing in the same pass.

    The unchanged upstream reader has a 1000-byte line buffer. Validate full
    gzip streams, four LF lines, ASCII identifiers, sequence/quality and pairing.
    Existing check_pairs behavior is deliberately unchanged.
    """
    def record(handle, label):
        lines = [handle.readline(1001) for _ in range(4)]
        if not lines[0]:
            if any(lines[1:]):
                raise ValueError(f"{label}: truncated FASTQ")
            return None
        if any(not line or len(line) > 999 or not line.endswith(b"\n") or b"\r" in line or b"\x00" in line for line in lines):
            raise ValueError(f"{label}: truncated/non-LF FASTQ or line exceeds trimLinker's 999-byte limit")
        if not lines[2].startswith(b"+"):
            raise ValueError(f"{label}: invalid FASTQ separator")
        seq, quality = lines[1][:-1], lines[3][:-1]
        if len(seq) != len(quality):
            raise ValueError(f"{label}: sequence/quality lengths differ")
        if any(base not in b"ACGTN" for base in seq):
            raise ValueError(f"{label}: BL sequence must use uppercase A/C/G/T/N; no silent canonicalization")
        if any(value < 33 or value > 126 for value in quality):
            raise ValueError(f"{label}: invalid FASTQ quality")
        try:
            name = normalized_name(lines[0].decode("ascii"))
        except (UnicodeDecodeError, IndexError) as exc:
            raise ValueError(f"{label}: invalid FASTQ identifier") from exc
        if not name:
            raise ValueError(f"{label}: empty FASTQ identifier")
        return name, lines, len(seq)

    result = {"pairs": 0, "uncompressed_bytes": 0, "min_length": None, "max_length": 0,
              "length_histogram": {}, "validation": "strict-blhic-v1"}
    with ExitStack() as stack:
        handles = []
        for path in (r1, r2):
            raw = stack.enter_context(Path(path).open("rb"))
            magic = raw.read(2)
            raw.seek(0)
            handles.append(stack.enter_context(gzip.GzipFile(fileobj=raw)) if magic == b"\x1f\x8b" else raw)
        check()
        while True:
            first, second = record(handles[0], "R1"), record(handles[1], "R2")
            if first is None or second is None:
                if first != second:
                    raise ValueError("paired FASTQs have different record counts")
                break
            if first[0] != second[0]:
                raise ValueError(f"paired read-name mismatch at pair {result['pairs'] + 1}")
            result["pairs"] += 1
            for mate, item in enumerate((first, second)):
                _, lines, length = item
                result["uncompressed_bytes"] += sum(map(len, lines))
                result["min_length"] = length if result["min_length"] is None else min(result["min_length"], length)
                result["max_length"] = max(result["max_length"], length)
                histogram = result["length_histogram"]
                histogram[str(length)] = histogram.get(str(length), 0) + 1
                if writers:
                    writers[mate].writelines(lines)
            if result["pairs"] % 10000 == 0:
                check()
        check()
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("r1", type=Path)
    parser.add_argument("r2", type=Path)
    args = parser.parse_args()
    try:
        count = check_pairs(args.r1, args.r2)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(f"OK: {count} synchronized FASTQ pair(s)")


if __name__ == "__main__":
    main()
