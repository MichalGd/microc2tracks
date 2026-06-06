#!/usr/bin/env python3
"""Validate the microc2tracks sample sheet."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = ["sample", "assay", "replicate_group", "condition", "fastq_r1", "fastq_r2"]
VALID_ASSAYS = {"microc", "hic"}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    if len(sys.argv) != 2:
        fail("Usage: validate_samplesheet.py config/samplesheet.csv")

    sheet = Path(sys.argv[1])
    if not sheet.exists():
        fail(f"Sample sheet does not exist: {sheet}")

    with sheet.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail("Sample sheet is empty")

        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            fail(f"Missing required column(s): {', '.join(missing)}")

        seen_samples: set[str] = set()
        rows = 0
        errors: list[str] = []

        for line_number, row in enumerate(reader, start=2):
            rows += 1
            sample = (row.get("sample") or "").strip()
            assay = (row.get("assay") or "").strip().lower()
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

            if not r1:
                errors.append(f"line {line_number}: fastq_r1 is empty")
            if not r2:
                errors.append(f"line {line_number}: fastq_r2 is empty")

            for label, value in [("fastq_r1", r1), ("fastq_r2", r2)]:
                if value and not Path(value).exists():
                    print(f"WARNING: line {line_number}: {label} path not found now: {value}", file=sys.stderr)

        if rows == 0:
            fail("Sample sheet has no data rows")

        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            raise SystemExit(1)

    print(f"OK: {sheet} has {rows} sample row(s)")


if __name__ == "__main__":
    main()
