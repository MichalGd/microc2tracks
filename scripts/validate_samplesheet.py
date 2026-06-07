#!/usr/bin/env python3
"""Validate the microc2tracks sample sheet."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample",
    "assay",
    "condition",
    "biological_replicate",
    "technical_replicate",
    "fastq_r1",
    "fastq_r2",
]
VALID_ASSAYS = {"microc", "hic"}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    require_files = False
    args = sys.argv[1:]

    if "--require-files" in args:
        require_files = True
        args.remove("--require-files")

    if len(args) != 1:
        fail("Usage: validate_samplesheet.py [--require-files] config/samplesheet.csv")

    sheet = Path(args[0])
    if not sheet.exists():
        fail(f"Sample sheet does not exist: {sheet}")

    with sheet.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail("Sample sheet is empty")

        reader.fieldnames = [column.strip().lstrip("\ufeff") for column in reader.fieldnames]
        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            fail(f"Missing required column(s): {', '.join(missing)}")

        seen_samples: set[str] = set()
        seen_technical_replicates: set[tuple[str, str, str, str]] = set()
        rows = 0
        errors: list[str] = []

        for line_number, row in enumerate(reader, start=2):
            rows += 1
            sample = (row.get("sample") or "").strip()
            assay = (row.get("assay") or "").strip().lower()
            condition = (row.get("condition") or "").strip()
            biological_replicate = (row.get("biological_replicate") or "").strip()
            technical_replicate = (row.get("technical_replicate") or "").strip()
            r1 = (row.get("fastq_r1") or "").strip()
            r2 = (row.get("fastq_r2") or "").strip()

            if not sample:
                errors.append(f"line {line_number}: sample is empty")
            elif sample in seen_samples:
                errors.append(f"line {line_number}: duplicate sample '{sample}'")
            else:
                seen_samples.add(sample)

            if assay not in VALID_ASSAYS:
                errors.append(f"line {line_number}: assay must be one of {sorted(VALID_ASSAYS)}")

            if not condition:
                errors.append(f"line {line_number}: condition is empty")

            if not biological_replicate:
                errors.append(f"line {line_number}: biological_replicate is empty")
            elif not biological_replicate.isdigit():
                errors.append(f"line {line_number}: biological_replicate must be an integer-like value")

            if not technical_replicate:
                errors.append(f"line {line_number}: technical_replicate is empty")
            elif not technical_replicate.isdigit():
                errors.append(f"line {line_number}: technical_replicate must be an integer-like value")

            tech_key = (assay, condition, biological_replicate, technical_replicate)
            if all(tech_key) and tech_key in seen_technical_replicates:
                errors.append(
                    f"line {line_number}: duplicate technical_replicate '{technical_replicate}' "
                    f"for assay={assay}, condition={condition}, biological_replicate={biological_replicate}"
                )
            else:
                seen_technical_replicates.add(tech_key)

            if not r1:
                errors.append(f"line {line_number}: fastq_r1 is empty")
            if not r2:
                errors.append(f"line {line_number}: fastq_r2 is empty")

            for label, value in [("fastq_r1", r1), ("fastq_r2", r2)]:
                if not value:
                    continue

                path = Path(value)
                if not path.exists():
                    message = f"line {line_number}: {label} path not found: {value}"
                    if require_files:
                        errors.append(message)
                    else:
                        print(f"WARNING: {message}", file=sys.stderr)
                elif not path.is_file():
                    errors.append(f"line {line_number}: {label} is not a regular file: {value}")
                elif not os.access(path, os.R_OK):
                    errors.append(f"line {line_number}: {label} is not readable: {value}")

        if rows == 0:
            fail("Sample sheet has no data rows")

        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            raise SystemExit(1)

    print(f"OK: {sheet} has {rows} sample row(s)")


if __name__ == "__main__":
    main()
