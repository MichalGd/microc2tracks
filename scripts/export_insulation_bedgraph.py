#!/usr/bin/env python3
"""Export bedGraph tracks from a cooltools insulation table."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_windows(value: str) -> list[str]:
    return [item for item in re.split(r"[,\s]+", value.strip()) if item]


def parse_float(value: str) -> float | None:
    if value == "" or value.lower() in {"nan", "none", "null"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def discover_windows(fieldnames: list[str]) -> list[str]:
    windows = []
    pattern = re.compile(r"^log2_insulation_score_(\d+)$")
    for field in fieldnames:
        match = pattern.match(field)
        if match:
            windows.append(match.group(1))
    if not windows:
        fail("No log2_insulation_score_<window> columns found")
    return sorted(windows, key=int)


def sanitize_track_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export log2 insulation score bedGraph tracks from cooltools insulation output."
    )
    parser.add_argument("--insulation", required=True, help="cooltools insulation TSV")
    parser.add_argument("--out-dir", required=True, help="output directory for bedGraph files")
    parser.add_argument("--sample", required=True, help="sample or merge name for output files")
    parser.add_argument(
        "--windows",
        default="",
        help="space- or comma-separated insulation windows to export; default: all detected windows",
    )
    parser.add_argument(
        "--track-line",
        action="store_true",
        help="prepend a UCSC bedGraph track line",
    )
    args = parser.parse_args()

    insulation_path = Path(args.insulation)
    out_dir = Path(args.out_dir)
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

        windows = parse_windows(args.windows) if args.windows.strip() else discover_windows(reader.fieldnames)
        columns = {window: f"log2_insulation_score_{window}" for window in windows}
        missing_columns = [column for column in columns.values() if column not in reader.fieldnames]
        if missing_columns:
            fail(f"Insulation table missing column(s): {', '.join(missing_columns)}")

        rows = list(reader)

    out_dir.mkdir(parents=True, exist_ok=True)
    sample = args.sample
    track_sample = sanitize_track_name(args.sample)

    for window in windows:
        column = columns[window]
        out_path = out_dir / f"{sample}.insulation_score.{window}.bedGraph"
        written = 0
        with out_path.open("w", newline="") as out_handle:
            if args.track_line:
                out_handle.write(
                    'track type=bedGraph '
                    f'name="{track_sample}_insulation_{window}" '
                    f'description="{track_sample} log2 insulation score {window} bp" '
                    "visibility=full autoScale=on color=180,40,40\n"
                )

            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            for row in rows:
                value = parse_float(row.get(column, ""))
                if value is None:
                    continue
                writer.writerow([row["chrom"], row["start"], row["end"], f"{value:.8g}"])
                written += 1

        print(f"Wrote {written} bedGraph row(s): {out_path}")


if __name__ == "__main__":
    main()
