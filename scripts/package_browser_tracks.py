#!/usr/bin/env python3
"""Package UCSC and HiGlass-ready bedGraph tracks from downstream outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path


MOUSE_CANONICAL = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
HUMAN_CANONICAL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


@dataclass(frozen=True)
class SampleSpec:
    label: str
    merged_dir: Path


@dataclass(frozen=True)
class BedGraphRow:
    chrom: str
    start: int
    end: int
    value: float


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_sample(value: str) -> SampleSpec:
    if "=" not in value:
        fail(f"--sample must be LABEL=/path/to/merged_dir, got: {value}")
    label, path = value.split("=", 1)
    label = sanitize_token(label)
    if not label:
        fail(f"Empty sample label in --sample {value!r}")
    merged_dir = Path(path)
    if not merged_dir.exists():
        fail(f"Merged directory not found for {label}: {merged_dir}")
    return SampleSpec(label=label, merged_dir=merged_dir)


def sanitize_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value.strip()).strip("_")


def parse_float(value: str) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def parse_int(value: str) -> int | None:
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return None


def canonical_chroms(preset: str, custom: str) -> list[str]:
    if custom.strip():
        return [normalize_ucsc_chrom(item.strip()) for item in custom.split(",") if item.strip()]
    if preset == "mouse":
        return MOUSE_CANONICAL
    if preset == "human":
        return HUMAN_CANONICAL
    fail(f"Unsupported canonical preset: {preset}")


def normalize_ucsc_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        if chrom in {"chrMT", "chrMt"}:
            return "chrM"
        return chrom
    if chrom in {"M", "MT", "Mt", "mitochondria"}:
        return "chrM"
    return f"chr{chrom}"


def chrom_sort_key(chrom_order: dict[str, int], row: BedGraphRow) -> tuple[int, int, int]:
    return (chrom_order.get(row.chrom, 10**9), row.start, row.end)


def find_one(patterns: list[str], root: Path, label: str, kind: str) -> Path:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(root.glob(pattern)))
    matches = [path for path in matches if path.is_file()]
    if not matches:
        fail(f"Could not find {kind} for {label} under {root}")
    if len(matches) > 1:
        names = "\n  ".join(str(path) for path in matches)
        fail(f"Found multiple {kind} files for {label}; pass a more specific resolution:\n  {names}")
    return matches[0]


def find_insulation_table(sample: SampleSpec, tad_resolution: str) -> Path:
    root = sample.merged_dir / "05_downstream" / "insulation"
    patterns = [f"*.insulation.{tad_resolution}.tsv"] if tad_resolution else ["*.insulation.*.tsv"]
    return find_one(patterns, root, sample.label, "insulation table")


def find_eigs_table(sample: SampleSpec, compartment_resolution: str) -> Path:
    root = sample.merged_dir / "05_downstream" / "compartments"
    patterns = (
        [f"*.eigs.{compartment_resolution}.cis.vecs.tsv"]
        if compartment_resolution
        else ["*.eigs.*.cis.vecs.tsv", "*.cis.vecs.tsv"]
    )
    return find_one(patterns, root, sample.label, "compartment eigenvector table")


def choose_pc_column(fieldnames: list[str], requested: str) -> str:
    if requested:
        if requested in fieldnames:
            return requested
        fail(f"Requested PC column not found: {requested}")
    for column in ["E1", "PC1", "eig1", "Eigenvector1", "e1", "pc1"]:
        if column in fieldnames:
            return column
    matches: list[tuple[int, str]] = []
    for column in fieldnames:
        match = re.match(r"^[Ee](\d+)$", column)
        if match:
            matches.append((int(match.group(1)), column))
    if matches:
        matches.sort()
        return matches[0][1]
    fail("Could not detect PC1/E1 column in compartment table")


def read_track_rows(
    path: Path,
    value_column: str,
    canonical: set[str],
    flip_value: bool = False,
) -> list[BedGraphRow]:
    rows: list[BedGraphRow] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"Input table is empty: {path}")
        missing = {"chrom", "start", "end", value_column} - set(reader.fieldnames)
        if missing:
            fail(f"{path} missing required column(s): {', '.join(sorted(missing))}")
        for record in reader:
            chrom = normalize_ucsc_chrom(record["chrom"])
            if chrom not in canonical:
                continue
            start = parse_int(record["start"])
            end = parse_int(record["end"])
            value = parse_float(record.get(value_column, ""))
            if start is None or end is None or value is None:
                continue
            if end <= start:
                continue
            if flip_value:
                value *= -1.0
            rows.append(BedGraphRow(chrom=chrom, start=start, end=end, value=value))
    return rows


def write_bedgraph_gz(
    path: Path,
    rows: list[BedGraphRow],
    chrom_order: dict[str, int],
    track_name: str | None,
    description: str | None,
    color: str,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(rows, key=lambda row: chrom_sort_key(chrom_order, row))
    with gzip.open(path, "wt", newline="") as handle:
        if track_name is not None and description is not None:
            handle.write(
                "track type=bedGraph "
                f"name={sanitize_token(track_name)} "
                f"description={sanitize_token(description)} "
                f"visibility=full autoScale=on color={color}\n"
            )
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in sorted_rows:
            writer.writerow([row.chrom, row.start, row.end, f"{row.value:.8g}"])
    return len(sorted_rows)


def chrom_sizes_candidates(sample: SampleSpec) -> list[Path]:
    return [
        *sample.merged_dir.glob("03_pairs/*.chrom.sizes"),
        *sample.merged_dir.glob("03_pairs/*.matrix.chrom.sizes"),
        *sample.merged_dir.glob("**/*.chrom.sizes"),
    ]


def read_chrom_sizes(samples: list[SampleSpec], canonical: set[str]) -> dict[str, int]:
    sizes: dict[str, int] = {}
    for sample in samples:
        seen: set[Path] = set()
        for path in chrom_sizes_candidates(sample):
            if path in seen or not path.is_file():
                continue
            seen.add(path)
            with path.open(newline="") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        fields = line.split()
                    if len(fields) < 2:
                        continue
                    chrom = normalize_ucsc_chrom(fields[0])
                    size = parse_int(fields[1])
                    if chrom in canonical and size is not None:
                        sizes[chrom] = max(size, sizes.get(chrom, 0))
    return sizes


def fill_missing_chrom_sizes(
    sizes: dict[str, int],
    rows_by_sample: dict[str, dict[str, list[BedGraphRow]]],
) -> dict[str, int]:
    filled = dict(sizes)
    for tracks in rows_by_sample.values():
        for rows in tracks.values():
            for row in rows:
                filled[row.chrom] = max(row.end, filled.get(row.chrom, 0))
    return filled


def write_chrom_sizes(path: Path, sizes: dict[str, int], canonical_order: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom in canonical_order:
            size = sizes.get(chrom)
            if size:
                writer.writerow([chrom, size])


def write_combined_higlass_pc1(
    path: Path,
    rows_by_label: dict[str, list[BedGraphRow]],
    labels: list[str],
    chrom_order: dict[str, int],
) -> int:
    maps = {
        label: {(row.chrom, row.start, row.end): row.value for row in rows}
        for label, rows in rows_by_label.items()
    }
    common_keys: set[tuple[str, int, int]] | None = None
    for label in labels:
        keys = set(maps[label])
        common_keys = keys if common_keys is None else common_keys & keys
    common_keys = common_keys or set()

    def key_sort(item: tuple[str, int, int]) -> tuple[int, int, int]:
        chrom, start, end = item
        return (chrom_order.get(chrom, 10**9), start, end)

    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, start, end in sorted(common_keys, key=key_sort):
            writer.writerow(
                [chrom, start, end]
                + [f"{maps[label][(chrom, start, end)]:.8g}" for label in labels]
            )
    return len(common_keys)


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["label", "reference_id", "browser_preset", "track_type", "file", "value", "source"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_higlass_ingest_notes(
    path: Path,
    chrom_sizes_path: Path,
    higlass_files: list[tuple[str, Path]],
) -> None:
    lines = [
        "# HiGlass ingestion helper",
        "# These files have no UCSC track line and contain only canonical chromosomes.",
        "# If your server uses clodius, create BEDDB tilesets like this:",
        "",
    ]
    for label, file_path in higlass_files:
        beddb = file_path.with_suffix("").with_suffix(".beddb")
        tmp = file_path.with_suffix("")
        lines.extend(
            [
                f"gunzip -c {shlex.quote(str(file_path))} > {shlex.quote(str(tmp))}",
                (
                    "clodius aggregate bedfile "
                    f"--chromsizes-filename {shlex.quote(str(chrom_sizes_path))} "
                    f"--output-file {shlex.quote(str(beddb))} {shlex.quote(str(tmp))}"
                ),
                f"# Then ingest {shlex.quote(str(beddb))} as the {label} compartment_PC1 track in your HiGlass server.",
                "",
            ]
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create canonical, gzipped UCSC bedGraph tracks and HiGlass-clean "
            "compartment PC1 tracks from microc2tracks merged downstream outputs."
        )
    )
    parser.add_argument(
        "--sample",
        action="append",
        required=True,
        help="LABEL=/path/to/merged_dir; repeat for each merged matrix",
    )
    parser.add_argument("--out-dir", required=True, help="output package directory")
    parser.add_argument("--tad-resolution", default="10000", help="insulation table resolution")
    parser.add_argument("--insulation-window", default="100000", help="insulation score window column")
    parser.add_argument("--compartment-resolution", default="100000", help="compartment eigenvector resolution")
    parser.add_argument("--pc-column", default="", help="PC column to export; default: auto-detect E1/PC1")
    parser.add_argument(
        "--flip-pc1",
        action="append",
        default=[],
        help="sample label whose PC1 sign should be multiplied by -1; repeat if needed",
    )
    parser.add_argument(
        "--reference-id",
        choices=["mm39", "hg38"],
        required=True,
        help="assembly identity; prevents human tracks from silently using the mouse preset",
    )
    parser.add_argument("--canonical-preset", choices=["mouse", "human"], default="")
    parser.add_argument(
        "--canonical-chroms",
        default="",
        help="optional comma-separated canonical chromosome list overriding preset",
    )
    parser.add_argument("--color", default="64,64,64", help="UCSC track color, default dark gray")
    args = parser.parse_args()

    samples = [parse_sample(value) for value in args.sample]
    if len({sample.label for sample in samples}) != len(samples):
        fail("Sample labels must be unique")

    expected_preset = {"mm39": "mouse", "hg38": "human"}[args.reference_id]
    if args.canonical_preset and args.canonical_preset != expected_preset:
        fail(
            f"--reference-id {args.reference_id} requires canonical preset {expected_preset}, "
            f"not {args.canonical_preset}"
        )
    canonical_preset = args.canonical_preset or expected_preset
    canonical_order = canonical_chroms(canonical_preset, args.canonical_chroms)
    canonical = set(canonical_order)
    chrom_order = {chrom: index for index, chrom in enumerate(canonical_order)}
    out_dir = Path(args.out_dir)
    ucsc_dir = out_dir / "ucsc_bedgraph"
    higlass_dir = out_dir / "higlass"
    manifest_rows: list[dict[str, str]] = []
    rows_by_sample: dict[str, dict[str, list[BedGraphRow]]] = {}
    pc1_higlass_files: list[tuple[str, Path]] = []

    for sample in samples:
        insulation_table = find_insulation_table(sample, args.tad_resolution)
        insulation_column = f"log2_insulation_score_{args.insulation_window}"
        eigs_table = find_eigs_table(sample, args.compartment_resolution)
        with eigs_table.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                fail(f"Compartment eigenvector table is empty: {eigs_table}")
            pc_column = choose_pc_column(reader.fieldnames, args.pc_column)

        flip_pc1 = sample.label in set(args.flip_pc1)
        insulation_rows = read_track_rows(insulation_table, insulation_column, canonical)
        pc1_rows = read_track_rows(eigs_table, pc_column, canonical, flip_value=flip_pc1)
        rows_by_sample[sample.label] = {"insulation": insulation_rows, "pc1": pc1_rows}

        if not insulation_rows:
            fail(f"No valid canonical insulation rows written for {sample.label}")
        if not pc1_rows:
            fail(f"No valid canonical PC1 rows written for {sample.label}")

        ucsc_insulation = (
            ucsc_dir / f"{sample.label}.insulation.{args.insulation_window}.ucsc.bedGraph.gz"
        )
        ucsc_pc1 = ucsc_dir / f"{sample.label}.pc1.{args.compartment_resolution}.ucsc.bedGraph.gz"
        higlass_pc1 = (
            higlass_dir / f"{sample.label}.pc1.{args.compartment_resolution}.higlass.bedGraph.gz"
        )

        n_ins = write_bedgraph_gz(
            ucsc_insulation,
            insulation_rows,
            chrom_order,
            f"{sample.label}_insulation_{args.insulation_window}",
            f"{sample.label}_log2_insulation_score_{args.insulation_window}",
            args.color,
        )
        n_pc1_ucsc = write_bedgraph_gz(
            ucsc_pc1,
            pc1_rows,
            chrom_order,
            f"{sample.label}_compartment_PC1_{args.compartment_resolution}",
            f"{sample.label}_compartment_PC1_{args.compartment_resolution}",
            args.color,
        )
        n_pc1_higlass = write_bedgraph_gz(
            higlass_pc1,
            pc1_rows,
            chrom_order,
            None,
            None,
            args.color,
        )

        manifest_rows.extend(
            [
                {
                    "label": sample.label,
                    "track_type": "ucsc_bedgraph",
                    "file": str(ucsc_insulation),
                    "value": insulation_column,
                    "source": str(insulation_table),
                },
                {
                    "label": sample.label,
                    "track_type": "ucsc_bedgraph",
                    "file": str(ucsc_pc1),
                    "value": pc_column + ("_flipped" if flip_pc1 else ""),
                    "source": str(eigs_table),
                },
                {
                    "label": sample.label,
                    "track_type": "higlass_bedgraph_no_trackline",
                    "file": str(higlass_pc1),
                    "value": pc_column + ("_flipped" if flip_pc1 else ""),
                    "source": str(eigs_table),
                },
            ]
        )
        pc1_higlass_files.append((sample.label, higlass_pc1))
        print(f"{sample.label}: wrote {n_ins} UCSC insulation rows -> {ucsc_insulation}")
        print(f"{sample.label}: wrote {n_pc1_ucsc} UCSC PC1 rows -> {ucsc_pc1}")
        print(f"{sample.label}: wrote {n_pc1_higlass} HiGlass PC1 rows -> {higlass_pc1}")

    sizes = fill_missing_chrom_sizes(read_chrom_sizes(samples, canonical), rows_by_sample)
    chrom_sizes_path = higlass_dir / f"{canonical_preset}.canonical.chrom.sizes"
    write_chrom_sizes(chrom_sizes_path, sizes, canonical_order)

    labels = [sample.label for sample in samples]
    combined_path = (
        higlass_dir
        / f"{'_'.join(labels)}.pc1.{args.compartment_resolution}.higlass.multisample.bed.gz"
    )
    combined_rows = write_combined_higlass_pc1(
        combined_path,
        {label: rows_by_sample[label]["pc1"] for label in labels},
        labels,
        chrom_order,
    )
    manifest_rows.append(
        {
            "label": "_".join(labels),
            "track_type": "higlass_multisample_bed_no_trackline",
            "file": str(combined_path),
            "value": ",".join(f"column_{i + 4}:{label}_PC1" for i, label in enumerate(labels)),
            "source": "intersection_of_pc1_bins",
        }
    )

    for row in manifest_rows:
        row["reference_id"] = args.reference_id
        row["browser_preset"] = canonical_preset

    manifest_path = out_dir / "track_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)
    write_higlass_ingest_notes(higlass_dir / "ingest_commands.sh", chrom_sizes_path, pc1_higlass_files)

    print(f"Wrote canonical chrom sizes -> {chrom_sizes_path}")
    print(f"Wrote {combined_rows} combined HiGlass multisample PC1 rows -> {combined_path}")
    print(f"Wrote manifest -> {manifest_path}")


if __name__ == "__main__":
    main()
