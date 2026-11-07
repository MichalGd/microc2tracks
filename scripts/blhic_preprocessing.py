"""Checked HaeIII/two-step BL preparation; mapping remains in PairProcessor."""
from __future__ import annotations

import gzip
import json
import shutil
from contextlib import ExitStack
from pathlib import Path

from check_fastq_pair_names import inspect_pairs
from workflow_common import ROOT, identity, json_write, sha256, signature, truth

SOURCE_COMMIT = "e120726d6440b24034f70bc3c51c17f351fef496"
PRESET = {"id": "haeiii_two_step_v1", "linker_a": "ACGCGATATCTTATC", "linker_b": "AGTCAGATAAGATAT",
          "mode": 1, "keep_empty": 2, "substitutions": 1, "min_length": 15,
          "matching": "case-sensitive-exact-first-substitution-only", "validation": "strict-blhic-v1"}


def raw_bl_rows(plan):
    if plan.get("workflow_mode") != "fastq":
        return []
    return [row for row in plan["samples"].values() if row.get("preprocessing", {}).get("input_state") == "raw"]


def scientific_tool_identity(tool):
    if not isinstance(tool, dict):
        raise ValueError("BL preprocessing needs a pinned trimLinker build receipt")
    return {key: tool[key] for key in ("source_commit", "source_hashes", "binary_sha256")}


def identify_tool(c):
    configured = c.get("BLHIC_TRIMLINKER")
    bundled = ROOT / "tools/blhic/trimLinker"
    candidate = configured or (str(bundled) if bundled.is_file() else "trimLinker")
    found = shutil.which(candidate)
    if not found:
        raise RuntimeError("BL-Hi-C raw inputs need trimLinker; build the bundled source once with scripts/build_trimlinker.py, or set BLHIC_TRIMLINKER")
    path = Path(found).resolve()
    receipt_path = Path(str(path) + ".build.json")
    if not receipt_path.is_file():
        raise ValueError(f"Missing pinned build receipt: {receipt_path}; use scripts/build_trimlinker.py")
    receipt = json.loads(receipt_path.read_text())
    expected = json.loads((ROOT / "vendor/trimLinker/source.json").read_text())
    if (receipt.get("source_commit") != SOURCE_COMMIT or receipt.get("source_hashes") != expected["sha256"]
            or receipt.get("binary_sha256") != sha256(path)):
        raise ValueError("trimLinker executable/source identity differs from the pinned build receipt")
    return {"path": str(path), **scientific_tool_identity(receipt)}


def parse_stats(path, inputs, outputs):
    values = {}
    for line in path.read_text().splitlines():
        name, value = line.split(":", 1)
        if name in values or not value.strip().isdigit():
            raise ValueError("Invalid/duplicate trimLinker statistics")
        values[name] = int(value)
    required = {"Total PETs", "Expect PETs", "Expect both PETs", "Chim PETs", "1Empty PETs", "2Empty PETs", "Valid PETs"}
    if set(values) != required:
        raise ValueError("Incomplete trimLinker statistics")
    values["ambiguous"] = values["Total PETs"] - values["Expect PETs"] - values["Chim PETs"]
    values["too_short"] = values["Expect PETs"] - values["Valid PETs"]
    if (min(values.values()) < 0 or values["Expect PETs"] != values["Expect both PETs"] + values["1Empty PETs"] + values["2Empty PETs"]
            or values["Total PETs"] != inputs["pairs"] or values["Valid PETs"] != outputs["pairs"]):
        raise ValueError("trimLinker statistics do not reconcile with validated FASTQ pairs")
    return values


def validate_inputs(rt, row):
    from artifact_state import FileLock
    before = [identity(Path(row[k])) for k in ("fastq_r1", "fastq_r2")]
    cache = Path(rt.plan.get("metadata_dir", str(Path(rt.plan["roots"]["out"]) / "run_metadata"))) / "blhic_input_checks" / (signature({"inputs": before, "validator": PRESET["validation"]}) + ".json")
    with FileLock(Path(str(cache) + ".lock"), blocking=True):
        if cache.is_file():
            saved = json.loads(cache.read_text())
            if saved.get("inputs") == before and saved.get("validation") == PRESET["validation"]:
                return saved
        result = inspect_pairs(Path(row["fastq_r1"]), Path(row["fastq_r2"]), check=rt.check)
        if not result["pairs"]:
            raise ValueError("BL-Hi-C supplied input pair is empty")
        if before != [identity(Path(row[k])) for k in ("fastq_r1", "fastq_r2")]:
            raise ValueError("BL FASTQs changed during validation")
        result["inputs"] = before
        json_write(cache, result)
        return result


def prepare(rt, entity, row, temporary, work, run_fastp):
    """Execute inside the existing owned trim transaction; caller publishes it."""
    c = rt.c
    inputs = validate_inputs(rt, row)
    contract = row["preprocessing"]
    if contract["preset"] != PRESET:
        raise ValueError("Unsupported BL preset contract; use the version supplied by this release")
    raw = contract["input_state"] == "raw"
    source_paths = [Path(row[k]) for k in ("fastq_r1", "fastq_r2")]
    counts, bridge = None, None
    if raw:
        # The plain output cannot exceed validated input bytes. Reserve a second
        # conservative copy for gzip/publication and respect concurrent tasks.
        amount = inputs["uncompressed_bytes"]
        requests = {rt.plan["roots"]["scratch"]: amount}
        target_root = rt.plan["roots"]["scratch"] if c["LARGE_OUTPUT_STAGING"] == "scratch" else rt.plan["roots"]["trimmed"]
        requests[target_root] = requests.get(target_root, 0) + int(amount * 1.02) + 65536
        if c["LARGE_OUTPUT_STAGING"] == "scratch":
            destination = rt.plan["roots"]["trimmed"]
            requests[destination] = requests.get(destination, 0) + int(amount * 1.02) + 65536
        if rt.guard.scratch_budget_bytes and requests[rt.plan["roots"]["scratch"]] > rt.guard.scratch_budget_bytes:
            raise RuntimeError("BL plain-output bound exceeds SCRATCH_BUDGET_GB; choose bulk TMPDIR or increase the explicit budget")
        if not rt.guard.grow(f"{entity['id']}/trim", requests):
            raise RuntimeError("Insufficient reserved disk for BL plain FASTQs and compressed outputs; choose bulk TMPDIR and retry")
        tool = rt.versions["trimLinker"]
        if sha256(Path(tool["path"])) != tool["binary_sha256"]:
            raise ValueError("trimLinker executable changed after preflight")
        attempt = work / "blhic"
        attempt.mkdir()
        command = [tool["path"], "-m", str(PRESET["mode"]), "-k", str(PRESET["keep_empty"]), "-e", str(PRESET["substitutions"]), "-l", str(PRESET["min_length"]), "-t", c["THREADS_FASTP"],
                   "-A", PRESET["linker_a"], "-B", PRESET["linker_b"], "-o", str(attempt), "-n", "reads", *map(str, source_paths)]
        rt.command(command, rt.log(entity, "blhic"))
        source_paths = [attempt / f"reads_{mate}.valid.fastq" for mate in (1, 2)]
        bridge = inspect_pairs(*source_paths, check=rt.check) if truth(c["RUN_FASTP"]) else None
        if bridge is not None:
            counts = parse_stats(attempt / "reads.trim.stat", inputs, bridge)
        shutil.copyfile(attempt / "reads.trim.stat", temporary["blhic_stats"])
    else:
        command = None
        temporary["blhic_stats"].write_text("Bridge processing not performed: user-declared pretrimmed inputs; historical bridge counts unknown.\n")
    if truth(c["RUN_FASTP"]):
        run_fastp(source_paths)
        final = inspect_pairs(temporary["r1"], temporary["r2"], check=rt.check)
        if final["pairs"] > (bridge or inputs)["pairs"]:
            raise ValueError("fastp output pair count exceeds its validated input")
    else:
        with ExitStack() as stack:
            writers = [stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=stack.enter_context(temporary[label].open("wb")), mtime=0, compresslevel=3)) for label in ("r1", "r2")]
            final = inspect_pairs(*source_paths, writers=writers, check=rt.check)
        if raw:
            bridge = final
            counts = parse_stats(work / "blhic/reads.trim.stat", inputs, bridge)
    if not final["pairs"]:
        raise ValueError("BL preprocessing retained no paired reads")
    if inputs["inputs"] != [identity(Path(row[k])) for k in ("fastq_r1", "fastq_r2")]:
        raise ValueError("BL FASTQs changed during preprocessing")
    receipt = {"schema": 1, "preprocessing": contract, "input_state": contract["input_state"],
               "supplied_input": inputs, "original_raw_pairs": inputs["pairs"] if raw else None,
               "bridge": bridge, "bridge_counts": counts, "command": command,
               "tool": rt.versions.get("trimLinker") if raw else None,
               "fastp_performed": truth(c["RUN_FASTP"]), "mapping_input": final,
               "outputs": {label: {**identity(temporary[label]), "path": entity["paths"][f"trim_{label}"]} for label in ("r1", "r2")},
               "provenance": "internally verified bridge processing" if raw else "pretrimming user-declared; original archive counts unknown"}
    json_write(temporary["blhic_qc"], receipt)
    if raw:
        # Only this successful attempt's explicitly named temporary payloads.
        for path in source_paths:
            path.unlink()
        (work / "blhic/reads.trim.stat").unlink()
        (work / "blhic").rmdir()
