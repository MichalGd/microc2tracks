#!/usr/bin/env python3
"""Verify that pairs headers all declare one expected assembly."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


ASSEMBLY_PATTERN = re.compile(r"^#{1,2}\s*(?:genome_assembly|assembly)\s*:\s*(\S+)", re.I)


def declared_assembly(path: Path) -> str:
    opener = gzip.open if path.suffix.lower() in {".gz", ".gzip"} else Path.open
    kwargs = {"mode": "rt", "encoding": "utf-8", "errors": "replace"}
    with opener(path, **kwargs) as handle:  # type: ignore[arg-type]
        for line in handle:
            if not line.startswith("#"):
                break
            match = ASSEMBLY_PATTERN.match(line.strip())
            if match:
                return match.group(1)
    raise ValueError(f"pairs header does not declare genome_assembly/assembly: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("pairs", nargs="+", type=Path)
    args = parser.parse_args()
    errors: list[str] = []
    for path in args.pairs:
        try:
            assembly = declared_assembly(path)
            if assembly != args.assembly:
                errors.append(f"{path}: declared assembly {assembly!r}, expected {args.assembly!r}")
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    if errors:
        for error in errors:
            print(f"ERROR: prohibited cross-reference merge: {error}", file=sys.stderr)
        raise SystemExit(1)
    print(f"OK: {len(args.pairs)} pairs file(s) declare assembly {args.assembly}")


if __name__ == "__main__":
    main()
