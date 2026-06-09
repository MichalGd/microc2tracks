#!/usr/bin/env python3
"""Call TAD-like intervals from a cooltools insulation table."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


TRUE_VALUES = {"true", "1", "yes", "y", "t"}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_int(value: str, label: str) -> int:
    try:
        return int(value)
    except ValueError:
        fail(f"Could not parse {label} as integer: {value!r}")


def parse_float(value: str) -> float | None:
    if value == "" or value.lower() in {"nan", "none", "null"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def is_true(value: str) -> bool:
    return value.strip().lower() in TRUE_VALUES


def choose_column(fieldnames: list[str], prefix: str, window_bp: str | None) -> str:
    if window_bp:
        column = f"{prefix}_{window_bp}"
        if column in fieldnames:
            return column
        fail(f"Column {column!r} not found in insulation table")

    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)$")
    matches = []
    for field in fieldnames:
        match = pattern.match(field)
        if match:
            matches.append((int(match.group(1)), field))
    if not matches:
        fail(f"No {prefix}_<window> column found in insulation table")
    matches.sort()
    return matches[-1][1]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create BED-like TAD intervals between adjacent strong insulation "
            "boundaries reported by cooltools insulation."
        )
    )
    parser.add_argument("--insulation", required=True, help="cooltools insulation TSV")
    parser.add_argument("--out", required=True, help="output BED file")
    parser.add_argument(
        "--window-bp",
        default="",
        help="boundary window in bp; must match is_boundary_<window> columns",
    )
    parser.add_argument(
        "--min-size-bp",
        type=int,
        default=0,
        help="minimum TAD interval size to report",
    )
    parser.add_argument(
        "--min-boundary-strength",
        type=float,
        default=None,
        help="optional additional boundary-strength cutoff",
    )
    args = parser.parse_args()

    insulation_path = Path(args.insulation)
    out_path = Path(args.out)

    if not insulation_path.exists():
        fail(f"Insulation table not found: {insulation_path}")

    with insulation_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("Insulation table is empty")

        required = {"chrom", "start", "end"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            fail(f"Insulation table missing required column(s): {', '.join(missing)}")

        window = args.window_bp or None
        boundary_col = choose_column(reader.fieldnames, "is_boundary", window)
        strength_col = choose_column(
            reader.fieldnames,
            "boundary_strength",
            boundary_col.rsplit("_", 1)[-1],
        )
        selected_window = boundary_col.rsplit("_", 1)[-1]

        chrom_rows: dict[str, list[dict[str, str]]] = {}
        for row in reader:
            chrom = row["chrom"]
            chrom_rows.setdefault(chrom, []).append(row)

    tad_count = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")

        for chrom in sorted(chrom_rows):
            rows = sorted(chrom_rows[chrom], key=lambda item: parse_int(item["start"], "start"))
            boundaries: list[tuple[int, float | None]] = []

            for row in rows:
                if not is_true(row.get(boundary_col, "")):
                    continue

                strength = parse_float(row.get(strength_col, ""))
                if args.min_boundary_strength is not None:
                    if strength is None or strength < args.min_boundary_strength:
                        continue

                start = parse_int(row["start"], "start")
                end = parse_int(row["end"], "end")
                boundaries.append(((start + end) // 2, strength))

            for index in range(len(boundaries) - 1):
                left_pos, left_strength = boundaries[index]
                right_pos, right_strength = boundaries[index + 1]
                size = right_pos - left_pos
                if size <= 0 or size < args.min_size_bp:
                    continue

                tad_count += 1
                strengths = [value for value in (left_strength, right_strength) if value is not None]
                score = 0
                if strengths:
                    score = min(1000, max(0, int(round(sum(strengths) / len(strengths) * 1000))))

                writer.writerow(
                    [
                        chrom,
                        left_pos,
                        right_pos,
                        f"TAD_{tad_count}",
                        score,
                        ".",
                        left_pos,
                        right_pos,
                        f"{left_strength:.6g}" if left_strength is not None else "nan",
                        f"{right_strength:.6g}" if right_strength is not None else "nan",
                        selected_window,
                    ]
                )

    print(f"Wrote {tad_count} TAD-like interval(s) to {out_path}")


if __name__ == "__main__":
    main()
