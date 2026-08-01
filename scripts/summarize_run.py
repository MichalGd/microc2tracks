#!/usr/bin/env python3
"""Create a compact microc2tracks run report from pipeline outputs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from pathlib import Path
from typing import Iterable


def parse_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    pattern = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$")

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        match = pattern.match(line)
        if not match:
            continue

        key, value = match.groups()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        values[key] = value

    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_technical_manifest(path: Path) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    if not path.exists():
        return groups

    for row in read_tsv(path):
        group = row.get("merge_group", "")
        sample = row.get("sample", "")
        if group and sample:
            groups.setdefault(group, []).append(sample)

    return groups


def parse_pairtools_stats(path: Path) -> dict[str, str]:
    stats: dict[str, str] = {}
    if not path.exists():
        return stats

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if "\t" in line:
                key, value = line.split("\t", 1)
            elif ":" in line:
                key, value = line.split(":", 1)
            else:
                fields = line.split()
                if len(fields) < 2:
                    continue
                key, value = fields[0], fields[1]

            stats[key.strip()] = value.strip().split()[0]

    return stats


def parse_fastp_json(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}

    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    duplication = data.get("duplication", {})
    insert_size = data.get("insert_size", {})

    raw_reads = before.get("total_reads")
    passed_reads = after.get("total_reads")

    summary: dict[str, str] = {
        "raw_reads": as_text(raw_reads),
        "raw_read_pairs": as_text(divide(raw_reads, 2)),
        "passed_reads": as_text(passed_reads),
        "passed_read_pairs": as_text(divide(passed_reads, 2)),
        "fastp_pass_rate": as_text(percent(divide(passed_reads, raw_reads))),
        "duplication_rate": as_text(percent(duplication.get("rate"))),
        "insert_size_peak": as_text(insert_size.get("peak")),
    }

    for key in [
        "passed_filter_reads",
        "low_quality_reads",
        "too_many_N_reads",
        "too_short_reads",
        "too_long_reads",
    ]:
        if key in filtering:
            summary[key] = as_text(filtering[key])

    return summary


def parse_cooler_info(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}

    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    return {
        "cooler_nbins": as_text(data.get("nbins")),
        "cooler_nnz": as_text(data.get("nnz")),
        "cooler_sum": as_text(data.get("sum")),
    }


def parse_status(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    return rows[-1] if rows else {}


def divide(value: object, denominator: object) -> float | None:
    try:
        numerator = float(value)  # type: ignore[arg-type]
        denominator_float = float(denominator)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if denominator_float == 0:
        return None
    return numerator / denominator_float


def percent(value: object) -> float | None:
    try:
        return float(value) * 100  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def as_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.3f}"
    return str(value)


def file_size(path: Path) -> str:
    return str(path.stat().st_size) if path.exists() else ""


def exists_text(path: Path) -> str:
    return "yes" if path.exists() and path.stat().st_size > 0 else "no"


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def sample_summary_rows(
    manifest_rows: list[dict[str, str]],
    outdir: Path,
    mapq: str,
    base_resolution: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for row in manifest_rows:
        sample = row["sample"]
        sample_dir = Path(row.get("sample_dir") or outdir / sample)
        qc_dir = sample_dir / "01_qc"
        pairs_dir = sample_dir / "03_pairs"
        matrix_dir = sample_dir / "04_matrices"
        log_dir = sample_dir / "logs"

        fastp = parse_fastp_json(qc_dir / f"{sample}.fastp.json")
        dedup_stats = parse_pairtools_stats(pairs_dir / f"{sample}.dedup.stats.txt")
        valid_stats = parse_pairtools_stats(pairs_dir / f"{sample}.pairs.stats.txt")
        cooler_info = parse_cooler_info(matrix_dir / f"{sample}.raw.{base_resolution}.cool.info.json")
        status = parse_status(log_dir / f"{sample}.status.tsv")

        valid_pairs = pairs_dir / f"{sample}.valid.mapq{mapq}.pairs.gz"
        norm_mcool = matrix_dir / f"{sample}.norm.mcool"
        norm_hic = matrix_dir / f"{sample}.norm.hic"

        cis = valid_stats.get("cis", "")
        trans = valid_stats.get("trans", "")

        summary = {
            "sample": sample,
            "assay": row.get("assay", ""),
            "reference_id": row.get("reference_id", ""),
            "species": row.get("species", ""),
            "assembly": row.get("assembly", ""),
            "browser_preset": row.get("browser_preset", ""),
            "condition": row.get("condition", ""),
            "biological_replicate": row.get("biological_replicate", ""),
            "technical_replicate": row.get("technical_replicate", ""),
            "merge_group": row.get("merge_group", ""),
            "status": status.get("status", "unknown"),
            "runtime_seconds": status.get("runtime_seconds", ""),
            "raw_read_pairs": fastp.get("raw_read_pairs", ""),
            "passed_read_pairs": fastp.get("passed_read_pairs", ""),
            "fastp_pass_rate_percent": fastp.get("fastp_pass_rate", ""),
            "fastp_duplication_rate_percent": fastp.get("duplication_rate", ""),
            "insert_size_peak": fastp.get("insert_size_peak", ""),
            "dedup_total_nodups": dedup_stats.get("total_nodups", dedup_stats.get("total", "")),
            "valid_pairs_total": valid_stats.get("total", ""),
            "valid_pairs_cis": cis,
            "valid_pairs_trans": trans,
            "cis_trans_ratio": as_text(divide(cis, trans)),
            "cooler_sum": cooler_info.get("cooler_sum", ""),
            "cooler_nnz": cooler_info.get("cooler_nnz", ""),
            "norm_mcool": exists_text(norm_mcool),
            "norm_mcool_size_bytes": file_size(norm_mcool),
            "norm_hic": exists_text(norm_hic),
            "norm_hic_size_bytes": file_size(norm_hic),
            "valid_pairs_size_bytes": file_size(valid_pairs),
        }
        rows.append(summary)

    return rows


def merge_summary_rows(
    groups: dict[str, list[str]],
    merge_metadata: dict[str, dict[str, str]],
    outdir: Path,
    mapq: str,
    base_resolution: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for group, samples in sorted(groups.items()):
        metadata = merge_metadata.get(group, {})
        merge_dir = outdir / "merged" / group
        pairs_dir = merge_dir / "03_pairs"
        matrix_dir = merge_dir / "04_matrices"
        valid_pairs = pairs_dir / f"{group}.valid.mapq{mapq}.pairs.gz"
        norm_mcool = matrix_dir / f"{group}.norm.mcool"
        norm_hic = matrix_dir / f"{group}.norm.hic"
        stats = parse_pairtools_stats(pairs_dir / f"{group}.pairs.stats.txt")
        cooler_info = parse_cooler_info(matrix_dir / f"{group}.raw.{base_resolution}.cool.info.json")

        rows.append(
            {
                "merge_group": group,
                "reference_id": metadata.get("reference_id", ""),
                "assembly": metadata.get("assembly", ""),
                "browser_preset": metadata.get("browser_preset", ""),
                "technical_replicate_count": str(len(samples)),
                "samples": ",".join(samples),
                "status": "created" if norm_mcool.exists() or norm_hic.exists() else "not_created",
                "valid_pairs_total": stats.get("total", ""),
                "valid_pairs_cis": stats.get("cis", ""),
                "valid_pairs_trans": stats.get("trans", ""),
                "cooler_sum": cooler_info.get("cooler_sum", ""),
                "cooler_nnz": cooler_info.get("cooler_nnz", ""),
                "norm_mcool": exists_text(norm_mcool),
                "norm_mcool_size_bytes": file_size(norm_mcool),
                "norm_hic": exists_text(norm_hic),
                "norm_hic_size_bytes": file_size(norm_hic),
                "valid_pairs_size_bytes": file_size(valid_pairs),
            }
        )

    return rows


def file_manifest_rows(
    sample_rows: list[dict[str, str]],
    merge_rows: list[dict[str, str]],
    outdir: Path,
    mapq: str,
) -> list[dict[str, str]]:
    wanted: list[tuple[str, Path]] = []

    for row in sample_rows:
        sample = row["sample"]
        sample_dir = outdir / sample
        wanted.extend(
            [
                ("sample_valid_pairs", sample_dir / "03_pairs" / f"{sample}.valid.mapq{mapq}.pairs.gz"),
                ("sample_norm_mcool", sample_dir / "04_matrices" / f"{sample}.norm.mcool"),
                ("sample_norm_hic", sample_dir / "04_matrices" / f"{sample}.norm.hic"),
                ("sample_fastp_json", sample_dir / "01_qc" / f"{sample}.fastp.json"),
                ("sample_multiqc", sample_dir / "01_qc" / f"{sample}.multiqc.html"),
            ]
        )

    for row in merge_rows:
        group = row["merge_group"]
        merge_dir = outdir / "merged" / group
        wanted.extend(
            [
                ("merged_valid_pairs", merge_dir / "03_pairs" / f"{group}.valid.mapq{mapq}.pairs.gz"),
                ("merged_norm_mcool", merge_dir / "04_matrices" / f"{group}.norm.mcool"),
                ("merged_norm_hic", merge_dir / "04_matrices" / f"{group}.norm.hic"),
            ]
        )

    rows = []
    for category, path in wanted:
        if path.exists():
            rows.append({"category": category, "path": str(path), "size_bytes": str(path.stat().st_size)})

    global_multiqc = first_existing([outdir / "final_report" / "microc2tracks_multiqc.html"])
    if global_multiqc:
        rows.append(
            {
                "category": "global_multiqc",
                "path": str(global_multiqc),
                "size_bytes": str(global_multiqc.stat().st_size),
            }
        )

    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def html_table(title: str, rows: list[dict[str, str]]) -> str:
    if not rows:
        return f"<h2>{html.escape(title)}</h2><p>No rows.</p>"

    headers = list(rows[0].keys())
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(row.get(header, ''))}</td>" for header in headers) + "</tr>")

    return (
        f"<h2>{html.escape(title)}</h2>"
        "<div class=\"table-wrap\"><table>"
        f"<thead><tr>{header_html}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody>"
        "</table></div>"
    )


def write_html(path: Path, sample_rows: list[dict[str, str]], merge_rows: list[dict[str, str]], file_rows: list[dict[str, str]]) -> None:
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>microc2tracks final report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #1f2933; }}
    h1, h2 {{ color: #0f172a; }}
    .table-wrap {{ overflow-x: auto; margin-bottom: 28px; }}
    table {{ border-collapse: collapse; font-size: 13px; min-width: 100%; }}
    th, td {{ border: 1px solid #d9e2ec; padding: 6px 8px; text-align: left; white-space: nowrap; }}
    th {{ background: #f0f4f8; }}
    tr:nth-child(even) {{ background: #f8fafc; }}
  </style>
</head>
<body>
  <h1>microc2tracks final report</h1>
  {html_table("Sample Summary", sample_rows)}
  {html_table("Merged Technical Replicates", merge_rows)}
  {html_table("Important Files", file_rows)}
</body>
</html>
"""
    path.write_text(document, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize a completed microc2tracks run.")
    parser.add_argument("-c", "--config", required=True, type=Path)
    parser.add_argument("-s", "--samplesheet", required=True, type=Path)
    parser.add_argument("-o", "--outdir", required=True, type=Path)
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--technical-manifest", required=True, type=Path)
    parser.add_argument("--merge-manifest", required=True, type=Path)
    args = parser.parse_args()

    config = parse_config(args.config)
    run_outdir = Path(config.get("OUTDIR", "results"))
    mapq = config.get("MAPQ_THRESHOLD", "30")
    base_resolution = config.get("BASE_RESOLUTION", "1000")

    args.outdir.mkdir(parents=True, exist_ok=True)

    manifest_rows = read_tsv(args.sample_manifest)
    groups = read_technical_manifest(args.technical_manifest)
    merge_metadata = {
        row.get("merge_group", ""): row
        for row in read_tsv(args.merge_manifest)
        if row.get("merge_group")
    }
    sample_rows = sample_summary_rows(manifest_rows, run_outdir, mapq, base_resolution)
    merge_rows = merge_summary_rows(groups, merge_metadata, run_outdir, mapq, base_resolution)
    file_rows = file_manifest_rows(sample_rows, merge_rows, run_outdir, mapq)

    write_tsv(args.outdir / "microc2tracks_sample_summary.tsv", sample_rows)
    write_tsv(args.outdir / "microc2tracks_merge_summary.tsv", merge_rows)
    write_tsv(args.outdir / "microc2tracks_file_manifest.tsv", file_rows)
    write_html(args.outdir / "microc2tracks_summary.html", sample_rows, merge_rows, file_rows)

    print(f"Wrote final report to {args.outdir}")


if __name__ == "__main__":
    main()
