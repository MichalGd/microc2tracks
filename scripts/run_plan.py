#!/usr/bin/env python3
"""Resolve general references, profiles, disjoint dedup units, pools and storage paths."""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path

from reference_registry import Registry, read_sizes, validate_reference_files
from validate_samplesheet import validate_sheet
from workflow_common import VERSION, absolute, identity, json_write, load_config, resolve_policy, safe_id, signature, truth


def resolve_roots(c: dict) -> dict:
    out = absolute(c["OUTDIR"])
    roots = {"out": out, "pairs": absolute(c["PAIRS_OUTDIR"] or out),
            "matrices": absolute(c["MATRIX_OUTDIR"] or out), "trimmed": absolute(c["TRIMMED_OUTDIR"] or out),
            "cache": absolute(c["REFERENCE_CACHE_DIR"] or Path(out) / "reference_cache"),
            "scratch": absolute(c["TMPDIR"] or Path(out) / "tmp")}
    if c["PUBLIC_HIC_DIR"]:
        roots["public_hic"] = absolute(c["PUBLIC_HIC_DIR"])
    return roots


def artifact_paths(roots: dict, relative: str, name: str, c: dict) -> dict:
    meta, pairs, matrix = Path(roots["out"]) / relative, Path(roots["pairs"]) / relative / "03_pairs", Path(roots["matrices"]) / relative / "04_matrices"
    result = {"entity_dir": meta, "pairs_dir": pairs, "matrix_dir": matrix, "logs": meta / "logs",
              "done": meta / "logs/done", "qc_dir": meta / "01_qc",
              "sorted": pairs / f"{name}.sorted.pairs.gz", "dedup": pairs / f"{name}.dedup.pairs.gz",
              "dedup_stats": pairs / f"{name}.dedup.stats.txt",
              "accepted": pairs / f"{name}.valid.mapq{c['MAPQ_THRESHOLD']}.pairs.gz",
              "filter_qc": pairs / f"{name}.filter_qc.json", "stats": pairs / f"{name}.pairs.stats.txt",
              "chrom_sizes": pairs / f"{name}.matrix.chrom.sizes", "status": meta / "logs" / f"{name}.status.tsv",
              "policy": meta / "logs" / f"{name}.policy.json",
              "raw_cool": matrix / f"{name}.raw.{c['BASE_RESOLUTION']}.cool",
              "norm_cool": matrix / f"{name}.norm.{c['BASE_RESOLUTION']}.cool",
              "raw_mcool": matrix / f"{name}.raw.mcool", "norm_mcool": matrix / f"{name}.norm.mcool",
              "raw_hic": matrix / f"{name}.raw.hic", "norm_hic": matrix / f"{name}.norm.hic",
              "cool_info": matrix / f"{name}.raw.{c['BASE_RESOLUTION']}.cool.info.json",
              "juicer": pairs / f"{name}.valid.mapq{c['MAPQ_THRESHOLD']}.juicer.pairs.gz",
              "trim_r1": Path(roots["trimmed"]) / relative / "02_trimmed" / f"{name}.trim.R1.fastq.gz",
              "trim_r2": Path(roots["trimmed"]) / relative / "02_trimmed" / f"{name}.trim.R2.fastq.gz"}
    return {key: str(path) for key, path in result.items()}


def units_for_pool(members: list[str], units: dict) -> list[str]:
    selected, resolved = set(members), []
    for key, unit in units.items():
        overlap = selected.intersection(unit["members"])
        if overlap:
            if overlap != set(unit["members"]):
                raise ValueError(f"Pool selects only part of dedup unit {key}: selected={sorted(overlap)}, complete={unit['members']}")
            resolved.append(key)
    return sorted(resolved)


def build_plan(c: dict, sheet: Path, *, check_files=True) -> dict:
    registry = Registry(Path(c["REFERENCE_REGISTRY"]))
    default_id = registry.resolve(c["DEFAULT_REFERENCE_ID"]).reference_id
    rows, _, warnings = validate_sheet(sheet, registry, default_id, require_files=check_files,
                                       check_gzip_files=check_files and truth(c["CHECK_FASTQ_GZIP"]))
    roots = resolve_roots(c)
    references = {}
    for reference_id in sorted({r["reference_genome"] for r in rows}):
        ref = registry.resolve(reference_id)
        if reference_id == default_id and c["REFERENCE_FASTA"] and c["CHROM_SIZES"]:
            ref = replace(ref, fasta=c["REFERENCE_FASTA"], chrom_sizes=c["CHROM_SIZES"],
                          bwa_index_prefix=c["BWA_INDEX_PREFIX"] or c["REFERENCE_FASTA"],
                          canonical_regex=c["CANONICAL_CHROMS_REGEX"] or ref.canonical_regex,
                          phasing_track=c["PHASING_TRACK"] or ref.phasing_track)
        # Registry-relative reference paths; legacy overrides retain launch-directory semantics.
        base = Path.cwd() if reference_id == default_id and c["REFERENCE_FASTA"] else registry.path.resolve().parent
        ref = replace(ref, **{key: absolute(getattr(ref, key), base) for key in ("fasta", "bwa_index_prefix", "chrom_sizes")},
                      phasing_track=absolute(ref.phasing_track, base) if ref.phasing_track else "")
        if check_files:
            errors = validate_reference_files(ref)
            if errors:
                raise ValueError("\n".join(errors))
        values = asdict(ref)
        values["identity"] = {"fasta": identity(Path(ref.fasta), True), "chrom_sizes": identity(Path(ref.chrom_sizes), True)} if check_files else {"reference": values.copy()}
        if check_files:
            values["identity"]["bwa_index"] = [identity(Path(ref.bwa_index_prefix + ext)) for ext in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]
        references[reference_id] = values
    samples, entities, units, group_members, technical = {}, {}, {}, defaultdict(list), defaultdict(list)
    technical_specs = {}
    for row in rows:
        name = row["sample"]
        policy = resolve_policy(row, c)
        reference = references[row["reference_genome"]]
        policy["canonical_regex"] = reference["canonical_regex"] if policy["filter_canonical"] else ".*"
        row["fastq_r1"], row["fastq_r2"] = absolute(row["fastq_r1"]), absolute(row["fastq_r2"])
        group = row["dedup_group"]
        unit_id = f"dedup:{group}" if group else f"sample:{name}"
        unit = units.setdefault(unit_id, {"members": [], "policy": policy, "reference_id": row["reference_genome"], "entity": unit_id})
        if unit["policy"] != policy or unit["reference_id"] != row["reference_genome"]:
            raise ValueError(f"Incompatible reference/assay/profile/digest/filter policy in dedup_group {group}")
        unit["members"].append(name)
        samples[name] = {**row, "policy": policy, "unit": unit_id}
        entities[f"sample:{name}"] = {"id": f"sample:{name}", "name": name, "kind": "sample", "members": [name],
                                      "reference_id": row["reference_genome"], "policy": policy, "units": [f"sample:{name}"],
                                      "dedup_group": group, "paths": artifact_paths(roots, name, name, c)}
        if group:
            group_members[group].append(name)
        label = row["merge_group"] or re.sub(r"[^A-Za-z0-9_.-]", "_", f"{row['condition']}_{row['assay']}_B{row['biological_replicate']}_tech_merged")
        if not label[0].isalnum():
            label = "group_" + label
        safe_id(label)
        spec = (row["reference_genome"], row["assay"], row["condition"], row["biological_replicate"])
        if label in technical_specs and technical_specs[label] != spec:
            raise ValueError(f"Derived merge-name collision: {label}; choose distinct explicit merge_group values")
        technical_specs[label] = spec
        technical[label].append(name)
        samples[name]["resolved_merge_group"] = label
    for group, members in group_members.items():
        first = entities[f"sample:{members[0]}"]
        key = f"dedup:{group}"
        entities[key] = {**first, "id": key, "name": group, "kind": "dedup", "members": sorted(members),
                         "units": [key], "paths": artifact_paths(roots, f"dedup_groups/{group}", group, c)}
    pools = {}

    def add_pool(name, members, pool_type):
        if name in pools:
            raise ValueError(f"Pool name collision: {name}")
        resolved_units = units_for_pool(members, units)
        first = entities[resolved_units[0]]
        if any(entities[key]["policy"] != first["policy"] or entities[key]["reference_id"] != first["reference_id"] for key in resolved_units):
            raise ValueError(f"Pool {name} has incompatible reference, assay, profile, digest or filtering policy")
        key = f"pool:{safe_id(name)}"
        entities[key] = {**first, "id": key, "name": name, "kind": "pool", "pool_type": pool_type,
                         "members": sorted(members), "units": resolved_units, "dedup_group": "",
                         "paths": artifact_paths(roots, f"merged/{name}", name, c)}
        pools[name] = key

    if truth(c["RUN_MERGE_TECHNICAL_REPLICATES"]):
        for name, members in technical.items():
            if len(members) > 1:
                add_pool(name, members, "technical")
    if truth(c["RUN_MERGE_ALL_SAMPLES"]) and len(samples) > 1:
        add_pool(c["ALL_SAMPLES_MERGE_NAME"], list(samples), "all_samples")
    for key, unit in units.items():
        unit["members"].sort()
    if any(s["policy"]["profile"] == "legacy" for s in samples.values()):
        warnings.append("legacy Hi-C mode: no restriction-fragment filtering. Select restriction_hic with explicit enzymes to enable it.")
    return {"schema_version": 1, "workflow_version": VERSION, "roots": roots, "settings": c,
            "samplesheet": str(sheet.resolve()), "references": references, "samples": samples,
            "entities": entities, "dedup_units": units, "pools": pools, "warnings": warnings,
            "dependencies": {key: ([f"sample:{s}:sorted" for s in e["members"]] if e["kind"] == "dedup" else e["units"] if e["kind"] == "pool" else []) for key, e in entities.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", type=Path, required=True)
    parser.add_argument("-s", "--samplesheet", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()
    plan = build_plan(load_config(args.config), args.samplesheet)
    if args.output:
        json_write(args.output, plan)
    else:
        print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
