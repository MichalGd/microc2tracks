#!/usr/bin/env python3
"""Read, validate, and query the microc2tracks reference registry."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path


REGISTRY_COLUMNS = [
    "reference_id",
    "aliases",
    "species",
    "assembly",
    "fasta",
    "bwa_index_prefix",
    "chrom_sizes",
    "canonical_regex",
    "browser_preset",
    "phasing_track",
    "phasing_assembly",
    "annotation_metadata",
    "blacklist_metadata",
]
SAFE_ID = re.compile(r"^[a-z0-9][a-z0-9_.-]*$")
EXPECTED_CANONICAL = {
    "mm39": [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"],
    "hg38": [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"],
}
EXPECTED_BROWSER_PRESET = {"mm39": "mouse", "hg38": "human"}
BWA_MEM2_INDEX_SUFFIXES = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]


class RegistryError(ValueError):
    """A user-facing registry error."""


@dataclass(frozen=True)
class Reference:
    reference_id: str
    aliases: tuple[str, ...]
    species: str
    assembly: str
    fasta: str
    bwa_index_prefix: str
    chrom_sizes: str
    canonical_regex: str
    browser_preset: str
    phasing_track: str
    phasing_assembly: str
    annotation_metadata: str
    blacklist_metadata: str

    def output_lines(self) -> list[str]:
        return [
            self.reference_id,
            self.species,
            self.assembly,
            self.fasta,
            self.bwa_index_prefix,
            self.chrom_sizes,
            self.canonical_regex,
            self.browser_preset,
            self.phasing_track,
            self.annotation_metadata,
            self.blacklist_metadata,
        ]


class Registry:
    def __init__(self, path: Path):
        self.path = path
        self.references: dict[str, Reference] = {}
        self.aliases: dict[str, str] = {}
        self._load()

    def _load(self) -> None:
        if not self.path.is_file():
            raise RegistryError(f"reference registry not found: {self.path}")
        with self.path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise RegistryError("reference registry is empty")
            reader.fieldnames = [name.strip().lstrip("\ufeff") for name in reader.fieldnames]
            missing = [name for name in REGISTRY_COLUMNS if name not in reader.fieldnames]
            if missing:
                raise RegistryError(f"reference registry missing column(s): {', '.join(missing)}")
            for line_number, row in enumerate(reader, start=2):
                values = {name: (row.get(name) or "").strip() for name in REGISTRY_COLUMNS}
                if not any(values.values()):
                    continue
                reference_id = values["reference_id"].lower()
                if not SAFE_ID.fullmatch(reference_id):
                    raise RegistryError(f"line {line_number}: unsafe reference_id {reference_id!r}")
                if reference_id in self.references:
                    raise RegistryError(f"line {line_number}: duplicate reference_id {reference_id!r}")
                try:
                    re.compile(values["canonical_regex"])
                except re.error as exc:
                    raise RegistryError(
                        f"line {line_number}: invalid canonical_regex for {reference_id}: {exc}"
                    ) from exc
                aliases = tuple(
                    alias.strip().lower()
                    for alias in values["aliases"].split(",")
                    if alias.strip()
                )
                ref = Reference(
                    reference_id=reference_id,
                    aliases=aliases,
                    species=values["species"],
                    assembly=values["assembly"],
                    fasta=values["fasta"],
                    bwa_index_prefix=values["bwa_index_prefix"],
                    chrom_sizes=values["chrom_sizes"],
                    canonical_regex=values["canonical_regex"],
                    browser_preset=values["browser_preset"],
                    phasing_track=values["phasing_track"],
                    phasing_assembly=values["phasing_assembly"],
                    annotation_metadata=values["annotation_metadata"],
                    blacklist_metadata=values["blacklist_metadata"],
                )
                for required in ("species", "assembly", "fasta", "bwa_index_prefix", "chrom_sizes", "canonical_regex", "browser_preset"):
                    if not getattr(ref, required):
                        raise RegistryError(f"line {line_number}: {required} is empty for {reference_id}")
                if reference_id in EXPECTED_BROWSER_PRESET and ref.browser_preset != EXPECTED_BROWSER_PRESET[reference_id]:
                    raise RegistryError(
                        f"line {line_number}: {reference_id} browser_preset must be "
                        f"{EXPECTED_BROWSER_PRESET[reference_id]!r}, got {ref.browser_preset!r}"
                    )
                if ref.phasing_track and ref.phasing_assembly != ref.assembly:
                    raise RegistryError(
                        f"line {line_number}: phasing_assembly must equal assembly for {reference_id}"
                    )
                self.references[reference_id] = ref
                for alias in (reference_id, *aliases):
                    if not SAFE_ID.fullmatch(alias):
                        raise RegistryError(f"line {line_number}: unsafe alias {alias!r}")
                    if alias in self.aliases:
                        raise RegistryError(f"line {line_number}: duplicate reference alias {alias!r}")
                    self.aliases[alias] = reference_id
        if not self.references:
            raise RegistryError("reference registry has no entries")

    def resolve(self, value: str) -> Reference:
        alias = value.strip().lower()
        reference_id = self.aliases.get(alias)
        if reference_id is None:
            accepted = ", ".join(sorted(self.aliases))
            raise RegistryError(f"unknown reference {value!r}; accepted values: {accepted}")
        return self.references[reference_id]


def read_sizes(path: Path, label: str) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                raise RegistryError(f"{label} line {line_number} has fewer than two tab-separated fields")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise RegistryError(f"{label} line {line_number} has invalid length {fields[1]!r}") from exc
            if name in result:
                raise RegistryError(f"{label} contains duplicate chromosome {name!r}")
            result[name] = length
    return result


def require_readable(path: Path, label: str, errors: list[str]) -> bool:
    if not path.is_file():
        errors.append(f"{label} not found: {path}")
        return False
    if not os.access(path, os.R_OK):
        errors.append(f"{label} is not readable: {path}")
        return False
    return True


def validate_reference_files(ref: Reference) -> list[str]:
    errors: list[str] = []
    fasta = Path(ref.fasta)
    fai = Path(f"{ref.fasta}.fai")
    chrom_sizes = Path(ref.chrom_sizes)
    fasta_ok = require_readable(fasta, f"{ref.reference_id} FASTA", errors)
    fai_ok = require_readable(fai, f"{ref.reference_id} FASTA index (.fai)", errors)
    chrom_ok = require_readable(chrom_sizes, f"{ref.reference_id} chromosome sizes", errors)
    for suffix in BWA_MEM2_INDEX_SUFFIXES:
        require_readable(Path(f"{ref.bwa_index_prefix}{suffix}"), f"{ref.reference_id} BWA-MEM2 index component", errors)

    if fai_ok and chrom_ok:
        try:
            fai_sizes = read_sizes(fai, f"{ref.reference_id} FASTA index")
            chrom_values = read_sizes(chrom_sizes, f"{ref.reference_id} chromosome sizes")
            if fai_sizes != chrom_values:
                missing = sorted(set(fai_sizes) - set(chrom_values))
                extra = sorted(set(chrom_values) - set(fai_sizes))
                length_mismatch = sorted(
                    name for name in set(fai_sizes) & set(chrom_values) if fai_sizes[name] != chrom_values[name]
                )
                errors.append(
                    f"{ref.reference_id} FASTA index/chromosome-size mismatch "
                    f"(missing={missing}, extra={extra}, length_mismatch={length_mismatch})"
                )
            canonical = [name for name in chrom_values if re.fullmatch(ref.canonical_regex, name)]
            expected = EXPECTED_CANONICAL.get(ref.reference_id)
            if not canonical:
                errors.append(f"{ref.reference_id} canonical_regex retains no chromosomes")
            if expected is not None and canonical != expected:
                errors.append(
                    f"{ref.reference_id} canonical chromosomes must be {','.join(expected)}; "
                    f"retained {','.join(canonical) or 'none'}"
                )
            mito_variants = [name for name in chrom_values if name in {"chrM", "MT", "chrMT"}]
            if mito_variants != ["chrM"]:
                errors.append(
                    f"{ref.reference_id} must use chrM deliberately for UCSC-compatible matrices; "
                    f"found mitochondrial names {mito_variants or 'none'}"
                )
        except (OSError, RegistryError) as exc:
            errors.append(str(exc))

    if fasta_ok and fasta.stat().st_size == 0:
        errors.append(f"{ref.reference_id} FASTA is empty: {fasta}")
    if ref.phasing_track:
        require_readable(Path(ref.phasing_track), f"{ref.reference_id} phasing track", errors)
    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    resolve_parser = subparsers.add_parser("resolve", help="resolve an alias and print one field per line")
    resolve_parser.add_argument("--registry", required=True, type=Path)
    resolve_parser.add_argument("reference")

    validate_parser = subparsers.add_parser("validate", help="validate registry structure and optional files")
    validate_parser.add_argument("--registry", required=True, type=Path)
    validate_parser.add_argument("--reference", action="append", default=[])
    validate_parser.add_argument("--check-files", action="store_true")
    validate_parser.add_argument("--override-reference")
    validate_parser.add_argument("--override-fasta")
    validate_parser.add_argument("--override-bwa-index-prefix")
    validate_parser.add_argument("--override-chrom-sizes")
    validate_parser.add_argument("--override-canonical-regex")
    validate_parser.add_argument("--override-phasing-track")

    args = parser.parse_args()
    try:
        registry = Registry(args.registry)
        if args.command == "resolve":
            print("\n".join(registry.resolve(args.reference).output_lines()))
            return
        refs = [registry.resolve(value) for value in args.reference] if args.reference else list(registry.references.values())
        if args.override_reference:
            override_id = registry.resolve(args.override_reference).reference_id
            override_values = {
                "fasta": args.override_fasta,
                "bwa_index_prefix": args.override_bwa_index_prefix,
                "chrom_sizes": args.override_chrom_sizes,
                "canonical_regex": args.override_canonical_regex,
                "phasing_track": args.override_phasing_track,
            }
            refs = [
                replace(ref, **{key: value for key, value in override_values.items() if value is not None})
                if ref.reference_id == override_id
                else ref
                for ref in refs
            ]
        errors: list[str] = []
        if args.check_files:
            for ref in refs:
                errors.extend(validate_reference_files(ref))
        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            raise SystemExit(1)
        print(f"OK: {registry.path} defines {len(registry.references)} reference(s)")
        if args.check_files:
            print("OK: reference files, chromosome naming, and BWA-MEM2 indexes validated")
    except RegistryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
