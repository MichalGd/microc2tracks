#!/usr/bin/env python3
"""Stream paired FASTQs and fail if read names or record counts differ."""

from __future__ import annotations

import argparse
import gzip
import sys
from contextlib import contextmanager
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
