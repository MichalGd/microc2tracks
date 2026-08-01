#!/usr/bin/env python3
"""Validate and normalize a microc2tracks sample sheet."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path

from reference_registry import Registry, RegistryError


REQUIRED_COLUMNS = [
    "sample",
    "assay",
    "condition",
    "biological_replicate",
    "technical_replicate",
    "fastq_r1",
    "fastq_r2",
]
NORMALIZED_COLUMNS = [
    "sample",
    "assay",
    "reference_genome",
    "condition",
    "biological_replicate",
    "technical_replicate",
    "fastq_r1",
    "fastq_r2",
    "merge_group",
]
VALID_ASSAYS = {"microc", "hic"}
GZIP_SUFFIXES = {".gz", ".gzip"}
SAFE_MERGE_GROUP = re.compile(r"^[A-Za-z0-9_.-]+$")
SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def format_bytes(size: int) -> str:
    value = float(size)
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if value < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def check_gzip(path: Path, label: str, line_number: int, errors: list[str]) -> None:
    if path.suffix.lower() not in GZIP_SUFFIXES:
        return
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
    except (EOFError, OSError) as exc:
        errors.append(f"line {line_number}: {label} gzip integrity failed for {path}: {exc}")


def validate_sheet(
    sheet: Path,
    registry: Registry,
    default_reference: str,
    require_files: bool = False,
    check_gzip_files: bool = False,
) -> tuple[list[dict[str, str]], list[tuple[str, str, Path, int]], list[str]]:
    try:
        default_id = registry.resolve(default_reference).reference_id
    except RegistryError as exc:
        raise RegistryError(f"invalid default reference: {exc}") from exc

    warnings: list[str] = []
    file_summaries: list[tuple[str, str, Path, int]] = []
    normalized_rows: list[dict[str, str]] = []
    errors: list[str] = []

    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RegistryError("sample sheet is empty")
        reader.fieldnames = [column.strip().lstrip("\ufeff") for column in reader.fieldnames]
        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise RegistryError(f"missing required column(s): {', '.join(missing)}")
        has_reference_column = "reference_genome" in reader.fieldnames
        if not has_reference_column:
            warnings.append(
                f"reference_genome column is absent; defaulting every row to {default_id} for backward compatibility"
            )

        seen_samples: set[str] = set()
        seen_technical_replicates: set[tuple[str, str, str, str, str, str]] = set()
        merge_references: dict[tuple[str, str, str, str], str] = {}
        explicit_merge_specs: dict[str, tuple[str, str, str, str]] = {}

        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            assay = (row.get("assay") or "").strip().lower()
            condition = (row.get("condition") or "").strip()
            biological_replicate = (row.get("biological_replicate") or "").strip()
            technical_replicate = (row.get("technical_replicate") or "").strip()
            r1 = (row.get("fastq_r1") or "").strip()
            r2 = (row.get("fastq_r2") or "").strip()
            merge_group = (row.get("merge_group") or "").strip()
            requested_reference = (row.get("reference_genome") or "").strip()
            if not requested_reference:
                requested_reference = default_id
                if has_reference_column:
                    warnings.append(
                        f"line {line_number}: reference_genome is blank; defaulting to {default_id}"
                    )
            try:
                reference_id = registry.resolve(requested_reference).reference_id
            except RegistryError as exc:
                errors.append(f"line {line_number}: {exc}")
                reference_id = ""

            if not sample:
                errors.append(f"line {line_number}: sample is empty")
            elif not SAFE_SAMPLE.fullmatch(sample):
                errors.append(
                    f"line {line_number}: sample may contain only letters, digits, dot, underscore, and hyphen"
                )
            elif sample in seen_samples:
                errors.append(f"line {line_number}: duplicate sample {sample!r}")
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
            if merge_group and not SAFE_MERGE_GROUP.fullmatch(merge_group):
                errors.append(
                    f"line {line_number}: merge_group may contain only letters, digits, dot, underscore, and hyphen"
                )
            for label, value in (
                ("condition", condition),
                ("fastq_r1", r1),
                ("fastq_r2", r2),
            ):
                if any(character in value for character in "\t\r\n"):
                    errors.append(f"line {line_number}: {label} contains a prohibited tab or newline")

            tech_key = (reference_id, assay, condition, biological_replicate, technical_replicate, merge_group)
            if all(tech_key[:5]) and tech_key in seen_technical_replicates:
                errors.append(
                    f"line {line_number}: duplicate technical_replicate {technical_replicate!r} "
                    f"for reference={reference_id}, assay={assay}, condition={condition}, "
                    f"biological_replicate={biological_replicate}"
                )
            else:
                seen_technical_replicates.add(tech_key)

            merge_key = (merge_group, assay, condition, biological_replicate)
            previous_reference = merge_references.get(merge_key)
            if reference_id and previous_reference and previous_reference != reference_id:
                label = merge_group or f"assay={assay}, condition={condition}, biological_replicate={biological_replicate}"
                errors.append(
                    f"line {line_number}: prohibited cross-reference replicate merge for {label}: "
                    f"{previous_reference} versus {reference_id}"
                )
            elif reference_id:
                merge_references[merge_key] = reference_id

            if merge_group and reference_id:
                merge_spec = (reference_id, assay, condition, biological_replicate)
                previous_spec = explicit_merge_specs.get(merge_group)
                if previous_spec and previous_spec != merge_spec:
                    if previous_spec[0] != reference_id:
                        errors.append(
                            f"line {line_number}: prohibited cross-reference replicate merge for "
                            f"merge_group={merge_group}: {previous_spec[0]} versus {reference_id}"
                        )
                    else:
                        errors.append(
                            f"line {line_number}: incompatible rows in merge_group={merge_group}; "
                            "assay, condition, and biological_replicate must match"
                        )
                else:
                    explicit_merge_specs[merge_group] = merge_spec

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
                        warnings.append(message)
                elif not path.is_file():
                    errors.append(f"line {line_number}: {label} is not a regular file: {value}")
                elif not os.access(path, os.R_OK):
                    errors.append(f"line {line_number}: {label} is not readable: {value}")
                else:
                    file_summaries.append((sample, label, path, path.stat().st_size))
                    if check_gzip_files:
                        check_gzip(path, label, line_number, errors)

            normalized_rows.append(
                {
                    "sample": sample,
                    "assay": assay,
                    "reference_genome": reference_id,
                    "condition": condition,
                    "biological_replicate": biological_replicate,
                    "technical_replicate": technical_replicate,
                    "fastq_r1": r1,
                    "fastq_r2": r2,
                    "merge_group": merge_group,
                }
            )

    if not normalized_rows:
        raise RegistryError("sample sheet has no data rows")
    if errors:
        raise RegistryError("\n".join(errors))
    return normalized_rows, file_summaries, warnings


def write_normalized(path: Path, rows: list[dict[str, str]], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_COLUMNS, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-files", action="store_true")
    parser.add_argument("--summarize-files", action="store_true")
    parser.add_argument("--check-gzip", action="store_true")
    parser.add_argument(
        "--reference-registry",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "config" / "references.tsv",
    )
    parser.add_argument("--default-reference", default="mm39")
    parser.add_argument("--normalized-output", type=Path)
    parser.add_argument("--output-format", choices=["csv", "tsv"], default="csv")
    parser.add_argument("samplesheet", type=Path)
    args = parser.parse_args()

    if not args.samplesheet.exists():
        fail(f"sample sheet does not exist: {args.samplesheet}")
    try:
        registry = Registry(args.reference_registry)
        rows, file_summaries, warnings = validate_sheet(
            args.samplesheet,
            registry,
            args.default_reference,
            require_files=args.require_files,
            check_gzip_files=args.check_gzip,
        )
    except (OSError, RegistryError) as exc:
        for line in str(exc).splitlines():
            print(f"ERROR: {line}", file=sys.stderr)
        raise SystemExit(1) from exc

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if args.normalized_output:
        write_normalized(args.normalized_output, rows, "\t" if args.output_format == "tsv" else ",")
    print(f"OK: {args.samplesheet} has {len(rows)} sample row(s)")
    print("Resolved references: " + ", ".join(sorted({row["reference_genome"] for row in rows})))
    if args.summarize_files and file_summaries:
        print("FASTQ input summary:")
        for sample, label, path, size in file_summaries:
            print(f"  {sample} {label}: {format_bytes(size)} ({size} bytes) {path}")
    if args.check_gzip:
        print("FASTQ gzip integrity check finished")


if __name__ == "__main__":
    main()
