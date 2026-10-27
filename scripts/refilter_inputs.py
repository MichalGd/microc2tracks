#!/usr/bin/env python3
"""Validate and stream existing deduplicated pairs; never align or deduplicate."""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

from filter_contact_pairs import REQUIRED
from reference_registry import read_sizes
from workflow_common import absolute, identity, json_write, safe_id, text_open

NORMALIZER_VERSION = 1
FIRST_COLUMNS = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
STALE_FRAGMENT_COLUMNS = {"rfrag1", "rfrag2", "rfrag_start1", "rfrag_start2", "rfrag_end1", "rfrag_end2"}
RESERVED = {"merged", "dedup_groups", "run_metadata", "final_report", "tmp", "reference_cache"}


def source_identity(path: Path) -> dict:
    result = identity(path)
    stat = path.stat()
    result.update(device=stat.st_dev, inode=stat.st_ino)
    return result


def read_header(handle, label: Path) -> tuple[list[str], list[str], dict]:
    header, columns, assembly, sizes = [], None, None, {}
    first_row = None
    for line in handle:
        if not line.startswith("#"):
            if line.strip():
                first_row = line
                break
            continue
        line = line.rstrip("\r\n")
        header.append(line)
        if line.startswith("#columns:"):
            if columns is not None:
                raise ValueError(f"Duplicate #columns header: {label}")
            columns = line.split(":", 1)[1].split()
        elif line.startswith("#genome_assembly:"):
            if assembly is not None:
                raise ValueError(f"Duplicate assembly declaration: {label}")
            assembly = line.split(":", 1)[1].strip()
        elif line.startswith("#chromsize:"):
            fields = line.split(":", 1)[1].split()
            if len(fields) != 2 or fields[0] in sizes or not fields[1].isdigit() or int(fields[1]) <= 0:
                raise ValueError(f"Malformed/duplicate chromosome-size declaration: {label}")
            sizes[fields[0]] = int(fields[1])
    if columns is None or columns[:7] != FIRST_COLUMNS or len(set(columns)) != len(columns):
        raise ValueError(f"Expected standard pairs columns starting with {' '.join(FIRST_COLUMNS)}: {label}")
    missing = set(REQUIRED) - set(columns)
    if missing:
        raise ValueError(f"Pairs lack required mapping fields {sorted(missing)}: {label}; no MAPQ or pair-type values will be invented")
    if not assembly:
        raise ValueError(f"Pairs need an explicit #genome_assembly header matching the selected reference: {label}")
    return header, columns, {"assembly": assembly, "chrom_sizes": sizes, "first_row": first_row}


def validate_header(path: Path, reference: dict) -> dict:
    with text_open(path) as handle:
        _, columns, metadata = read_header(handle, path)
    if metadata["assembly"] != reference["assembly"]:
        raise ValueError(f"Input assembly {metadata['assembly']!r} does not match selected reference {reference['assembly']!r}: {path}")
    sizes = read_sizes(Path(reference["chrom_sizes"]), "reference")
    if any(sizes.get(chrom) != length for chrom, length in metadata["chrom_sizes"].items()):
        raise ValueError(f"Input chromosome-size declarations disagree with the reference: {path}")
    sidecar = Path(str(path) + ".provenance.json")
    previous = json.loads(sidecar.read_text()) if sidecar.is_file() else {}
    if not isinstance(previous, dict):
        raise ValueError(f"Expected a JSON object in upstream provenance: {sidecar}")
    if previous.get("pre_dedup") or previous.get("stage") == "sorted":
        raise ValueError(f"Input provenance identifies pre-dedup pairs, not deduplicated input: {path}")
    reference_hash = reference.get("identity", {}).get("fasta", {}).get("sha256")
    if reference_hash and previous.get("reference_identity", {}).get("fasta", {}).get("sha256") not in (None, reference_hash):
        raise ValueError(f"Input provenance identifies a different reference FASTA: {path}")
    members = previous.get("members") or None
    if members is not None:
        if not isinstance(members, list) or not members or any(not isinstance(member, str) for member in members):
            raise ValueError(f"Malformed upstream member metadata: {sidecar}")
        members = sorted(safe_id(member) for member in members)
        if len(set(members)) != len(members):
            raise ValueError(f"Repeated member in upstream provenance: {sidecar}")
    return {"columns": columns, "assembly": metadata["assembly"], "declared_chromosome_count": len(metadata["chrom_sizes"]),
            "previous_filter_policy": previous.get("policy") if previous.get("stage") == "accepted" else None, "previous_stage": previous.get("stage"),
            "known_members": members, "sidecar_identity": identity(sidecar, True) if sidecar.is_file() else None,
            "deduplication_status": "asserted by refilter-deduplicated mode; not independently verified",
            "alignment_reference_match": "assembly and declared chromosome lengths checked; original FASTA content unverified without provenance"}


def read_sheet(sheet: Path) -> tuple[list[dict], list[str]]:
    rows, warnings, seen_names, seen_files = [], [], set(), set()
    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Refilter sample sheet is empty")
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        required = {"sample", "assay", "reference_genome", "deduplicated_pairs"}
        if len(set(reader.fieldnames)) != len(reader.fieldnames) or not required.issubset(reader.fieldnames):
            raise ValueError("Refilter CSV requires unique columns: sample,assay,reference_genome,deduplicated_pairs")
        for number, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f"Line {number}: more values than CSV columns")
            raw = {name: (value or "").strip() for name, value in raw.items()}
            if not any(raw.values()):
                continue
            if any(any(ch in value for ch in "\r\n\t") for value in raw.values()):
                raise ValueError(f"Line {number}: prohibited control character in metadata")
            name = safe_id(raw["sample"])
            if name in RESERVED or name in seen_names:
                raise ValueError(f"Line {number}: reserved/duplicate sample name {name!r}")
            seen_names.add(name)
            if raw.get("dedup_group"):
                raise ValueError("refilter-deduplicated never runs deduplication; leave dedup_group empty and import each already-deduplicated unit once")
            if raw.get("fastq_r1") or raw.get("fastq_r2"):
                raise ValueError("Refilter CSV must not mix FASTQ and deduplicated-pairs inputs")
            if raw["assay"].lower() not in {"microc", "hic"} or not raw["reference_genome"]:
                raise ValueError(f"Line {number}: declare assay=microc/hic and an explicit registered reference_genome")
            path = Path(absolute(raw["deduplicated_pairs"]))
            if path.suffix.lower() in {".bam", ".sam", ".cram", ".cool", ".mcool", ".hic"}:
                raise ValueError("Refilter mode requires mapped deduplicated .pairs/.pairs.gz, not BAM/SAM or a matrix; no implicit alignment conversion is performed")
            if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
                raise ValueError(f"Missing, empty or unreadable deduplicated pairs: {path}")
            stat = path.stat()
            file_key = (stat.st_dev, stat.st_ino) if stat.st_ino else str(path)
            if file_key in seen_files:
                raise ValueError(f"Duplicate physical input file (including symlink/hardlink aliases): {path}")
            seen_files.add(file_key)
            row = {key: raw.get(key, "") for key in ("condition", "biological_replicate", "technical_replicate", "merge_group", "hic_protocol", "restriction_enzymes", "dedup_group")}
            for key in ("biological_replicate", "technical_replicate"):
                if row[key] and not row[key].isdigit():
                    raise ValueError(f"Line {number}: {key} must be integer-like when provided")
            if row["merge_group"]:
                safe_id(row["merge_group"])
            members = raw.get("input_members", "").split(";") if raw.get("input_members") else [name]
            members = [safe_id(value.strip()) for value in members]
            if len(set(members)) != len(members):
                raise ValueError(f"Line {number}: input_members repeats a member")
            rows.append({**row, "sample": name, "assay": raw["assay"].lower(), "reference_genome": raw["reference_genome"],
                         "fastq_r1": "", "fastq_r2": "", "deduplicated_pairs": str(path), "input_members": sorted(members),
                         "input_members_declared": bool(raw.get("input_members")),
                         "input_identity": source_identity(path)})
    if not rows:
        raise ValueError("Refilter sample sheet has no data rows")
    warnings.append("Refilter mode: input deduplication is assumed as declared; no alignment, parsing or deduplication is performed. Previously removed contacts cannot be recovered.")
    return rows, warnings


def validate_snapshot(row: dict) -> None:
    source = Path(row["deduplicated_pairs"])
    if source_identity(source) != row["input_identity"]:
        raise ValueError(f"Read-only input changed during reprocessing: {source}")
    sidecar = Path(str(source) + ".provenance.json")
    current = identity(sidecar, True) if sidecar.is_file() else None
    if current != row["input_header"]["sidecar_identity"]:
        raise ValueError(f"Upstream provenance changed during reprocessing: {sidecar}")


def normalize(source: Path, output, sizes: dict, expected: dict | None = None, check=None) -> dict:
    before = source_identity(source)
    if expected is not None and before != expected:
        raise ValueError(f"Input identity changed since planning: {source}")
    ranks = {chrom: rank for rank, chrom in enumerate(sizes)}
    ranks["!"] = -1
    count, flipped = 0, 0
    with text_open(source) as handle:
        header, columns, metadata = read_header(handle, source)
        kept = [column for column in columns if column not in STALE_FRAGMENT_COLUMNS]
        for line in header:
            if not line.startswith(("#columns:", "#sorted:", "#shape:", "#chromsize:")):
                output.write(line + "\n")
        for chrom, length in sizes.items():
            output.write(f"#chromsize: {chrom} {length}\n")
        output.write("#shape: upper triangle\n#columns: " + " ".join(kept) + "\n")

        def records():
            if metadata["first_row"] is not None:
                yield metadata["first_row"]
            yield from handle

        for line in records():
            if not line.strip():
                continue
            if line.startswith("#"):
                raise ValueError("Pairs header/comment appears after data")
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(columns):
                raise ValueError(f"Malformed input pairs row {count+1}: {source}")
            row = dict(zip(columns, values))
            count += 1
            for side in (1, 2):
                chrom, pos, strand, mapq = row[f"chrom{side}"], int(row[f"pos{side}"]), row[f"strand{side}"], int(row[f"mapq{side}"])
                if mapq < 0 or strand not in {"+", "-"} or (chrom == "!" and pos != 0):
                    raise ValueError(f"Invalid mapping fields at input row {count}: {source}")
                if chrom != "!" and (chrom not in sizes or not 1 <= pos <= sizes[chrom]):
                    raise ValueError(f"Mapped coordinate outside configured reference at input row {count}: {source}")
            if len(row["pair_type"]) != 2:
                raise ValueError(f"Invalid pair_type at input row {count}: {source}")
            if (ranks[row["chrom1"]], int(row["pos1"])) > (ranks[row["chrom2"]], int(row["pos2"])):
                for name in columns:
                    if name.endswith("1") and name[:-1] + "2" in row:
                        other = name[:-1] + "2"
                        row[name], row[other] = row[other], row[name]
                row["pair_type"] = row["pair_type"][::-1]
                flipped += 1
            output.write("\t".join(row[name] for name in kept) + "\n")
            if check and count % 100000 == 0:
                check()
    if source_identity(source) != before:
        raise ValueError(f"Input changed during reprocessing: {source}")
    return {"normalizer_version": NORMALIZER_VERSION, "input_identity": before, "input_records": count,
            "end_swapped_records": flipped, "removed_annotation_columns": sorted(set(columns) - set(kept)),
            "alignment_performed": False, "parsing_performed": False, "deduplication_performed": False,
            "deduplication_status": "assumed from user-selected input mode; not verified or repaired"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--expected-identity", type=Path, required=True)
    parser.add_argument("--qc", type=Path, required=True)
    args = parser.parse_args()
    qc = normalize(args.input, sys.stdout, read_sizes(args.chrom_sizes, "reference"), json.loads(args.expected_identity.read_text()))
    json_write(args.qc, qc)


if __name__ == "__main__":
    main()
