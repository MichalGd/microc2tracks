#!/usr/bin/env python3
"""Integrated microc2tracks runner behind the existing Bash entry points."""
from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import importlib.metadata
import json
import os
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import asdict, replace
from pathlib import Path

from artifact_state import FileLock, StageStore, StorageGuard, publish
from matrix_products import MatrixProducts
from pairs_processing import PairProcessor, proof_path, read_proof, validate_pool_inputs, verify_proof_payload
from reference_registry import Registry
from run_plan import artifact_paths, build_plan, resolve_roots
from workflow_common import ROOT, VERSION, absolute, identity, json_write, load_config, resolve_policy, safe_id, sha256, signature, text_open, truth, tsv_write


def log(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def storage_guard(plan):
    c = plan["settings"]
    return StorageGuard(plan["roots"], shlex.split(c["REQUIRED_STORAGE_MOUNTS"]), float(c["MIN_FREE_OUTPUT_GB"]), float(c["MIN_FREE_SCRATCH_GB"]))


def validate_existing_layout(plan):
    """Fail before writing statuses/manifests if this would clobber legacy outputs."""
    metadata = Path(plan.get("metadata_dir", str(Path(plan["roots"]["out"]) / "run_metadata")))
    saved_plan = metadata / "run_plan.json"
    if metadata.exists() and any(p.is_file() and not p.name.startswith(".") for p in metadata.iterdir()):
        if not saved_plan.is_file() or json.loads(saved_plan.read_text()).get("schema_version") != 1:
            raise RuntimeError(f"Unowned/legacy run metadata at {metadata}; use a new OUTDIR and preserve originals")
    for entity in plan["entities"].values():
        if entity["kind"] in {"external", "external_sorted"}:
            continue
        paths = entity["paths"]
        owned = set()
        for record_path in Path(paths["done"]).glob("*.json"):
            record = json.loads(record_path.read_text())
            owned.update(p["path"] for p in record.get("outputs", {}).values())
        for label in ("accepted", "sorted", "dedup", "chrom_sizes", "raw_cool", "norm_cool", "raw_mcool", "norm_mcool", "raw_hic", "norm_hic", "trim_r1", "trim_r2"):
            path = Path(paths[label])
            if path.exists() and str(path.resolve()) not in owned:
                raise RuntimeError(f"Unowned/legacy output {path}; use a new output tree. Original files were not modified.")


def tool_versions(plan):
    c = plan["settings"]
    matrices = bool(plan["pools"]) or truth(c["BUILD_SAMPLE_MATRICES"]) or any(e.get("build_matrices") for e in plan["entities"].values())
    tools = ["pairtools", "bgzip", "pairix"]
    if plan["samples"] and plan.get("workflow_mode") != "refilter-deduplicated":
        tools += ["bwa-mem2"]
        if truth(c["RUN_FASTP"]):
            tools += ["fastp"]
    if truth(c["RUN_GLOBAL_MULTIQC"]) or (plan["samples"] and truth(c["RUN_MULTIQC"])):
        tools += ["multiqc"]
    if matrices:
        tools += ["cooler"]
        if truth(c["RUN_HIC"]):
            tools += ["java"]
            if not Path(c["JUICER_TOOLS_JAR"]).is_file():
                raise ValueError("RUN_HIC=true requires JUICER_TOOLS_JAR")
    missing = [name for name in tools if not shutil.which(name)]
    if missing:
        raise RuntimeError("Missing required tools: " + ", ".join(missing) + "; create the release environment first")
    versions = {"python": platform.python_version(), "workflow": VERSION}
    for library in ("numpy", "scipy", "pandas", "h5py", "pysam", "cooler"):
        try:
            versions[library] = importlib.metadata.version(library)
        except importlib.metadata.PackageNotFoundError:
            versions[library] = "not installed"
    for optional in ("cooltools", "multiqc"):
        if shutil.which(optional):
            tools.append(optional)
    for name in tools:
        flag = "-version" if name == "java" else "version" if name == "bwa-mem2" else "--version"
        try:
            result = subprocess.run([name, flag], capture_output=True, text=True, timeout=10)
            versions[name] = (result.stdout + result.stderr).strip()[:2000]
        except subprocess.TimeoutExpired:
            versions[name] = "version query timed out"
    if not re.search(r"\b1\.1\.3\b", versions["pairtools"]):
        raise RuntimeError(f"This release targets pairtools 1.1.3; found: {versions['pairtools']}. Do not update a running production environment in place.")
    return versions


class Runtime:
    def __init__(self, plan, config_path: Path, versions=None):
        self.plan, self.c = plan, plan["settings"]
        self.config_path = config_path.resolve()
        self.guard = storage_guard(plan)
        self.guard.check()
        self.versions = versions or tool_versions(plan)
        matrices = bool(plan["pools"]) or truth(self.c["BUILD_SAMPLE_MATRICES"]) or any(e.get("build_matrices") for e in plan["entities"].values())
        self.jar_identity = identity(Path(self.c["JUICER_TOOLS_JAR"]), True) if matrices and self.c["JUICER_TOOLS_JAR"] and truth(self.c["RUN_HIC"]) else {}
        scratch_root = Path(plan["roots"]["scratch"])
        scratch_root.mkdir(parents=True, exist_ok=True)
        self.scratch = Path(tempfile.mkdtemp(prefix="run_", dir=scratch_root))
        json_write(self.scratch / "owner.json", {"run_id": self.scratch.name, "pid": os.getpid(), "version": VERSION})
        plan["run_id"], plan["scratch"] = self.scratch.name, str(self.scratch)
        self.slots = {key: threading.BoundedSemaphore(int(self.c[f"MAX_PARALLEL_{key.upper()}"])) for key in ("fastp", "matrix", "hic")}
        self.digest_lock = threading.Lock()
        self.digests = {}
        self.code_hashes = {p.name: sha256(p) for p in (ROOT / "scripts").iterdir() if p.suffix in {".py", ".sh"}}
        self.stop = threading.Event()
        self.storage_check_lock = threading.Lock()
        self.last_storage_check = time.monotonic()
        self.warning_lock, self.warnings = threading.Lock(), []

    def key(self, stage, parameters, tools=(), code=()):
        return signature({"stage": stage, "schema": 1, "parameters": parameters,
                          "tools": {name: self.versions.get(name, "unknown") for name in tools},
                          "libraries": {name: self.versions.get(name) for name in (
                              ["numpy", "scipy", "pandas", "pysam"] if "pairtools" in tools else
                              ["numpy", "pandas", "h5py"] if "cooler" in tools else [])},
                          "code": {name: self.code_hashes[name] for name in ["workflow_common.py", *code]}})

    def store(self, entity):
        roots = [Path(p) for k, p in self.plan["roots"].items() if k != "scratch"]
        if self.c["PUBLIC_HIC_DIR"]:
            roots.append(Path(self.c["PUBLIC_HIC_DIR"]).resolve())
        return StageStore(Path(entity["paths"]["done"]), self.scratch / entity["id"].replace(":", "_"), roots,
                          self.c["LARGE_OUTPUT_STAGING"], self.guard)

    def log(self, entity, stage):
        path = Path(entity["paths"]["logs"]) / f"{entity['name']}.{stage}.log"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def check(self):
        if self.stop.is_set():
            raise InterruptedError("Run interrupted; retaining unfinished stage data")
        if time.monotonic() - self.last_storage_check >= float(self.c["STORAGE_CHECK_INTERVAL"]):
            with self.storage_check_lock:
                if time.monotonic() - self.last_storage_check >= float(self.c["STORAGE_CHECK_INTERVAL"]):
                    self.guard.check()
                    self.last_storage_check = time.monotonic()

    def _terminate(self, process):
        if process.poll() is None:
            if os.name != "nt":
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGTERM)
            else:
                process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                if os.name != "nt":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
                process.wait()

    def pipeline(self, commands, output, log_path, timeout=None):
        """Pipe real subprocesses, checking every exit code (not only the last writer)."""
        self.check()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        processes, previous = [], None
        with contextlib.ExitStack() as handles:
            errors = handles.enter_context(log_path.open("ab"))
            out = handles.enter_context(Path(output).open("wb")) if output is not None else errors
            try:
                for i, args in enumerate(commands):
                    errors.write(("\n$ " + shlex.join(map(str, args)) + "\n").encode())
                    errors.flush()
                    process = subprocess.Popen(list(map(str, args)), stdin=previous,
                        stdout=out if i == len(commands)-1 else subprocess.PIPE, stderr=errors,
                        start_new_session=os.name != "nt")
                    if previous is not None:
                        previous.close()
                    previous = process.stdout
                    processes.append(process)
                started = time.monotonic()
                while any(p.poll() is None for p in processes):
                    if self.stop.is_set():
                        raise InterruptedError("Run interrupted or another worker failed; stopping this run's child processes")
                    self.check()
                    failed = [p.returncode for p in processes if p.poll() not in (None, 0)]
                    if failed:
                        raise RuntimeError(f"Pipeline failed with exit code {failed[0]}; log: {log_path}")
                    if timeout and time.monotonic()-started > timeout:
                        raise TimeoutError(f"Command exceeded {timeout}s; log: {log_path}")
                    time.sleep(min(float(self.c["STORAGE_CHECK_INTERVAL"]), 0.1))
                if any(p.returncode for p in processes):
                    raise RuntimeError(f"Pipeline failed ({[p.returncode for p in processes]}); log: {log_path}")
            except BaseException:
                for p in processes:
                    self._terminate(p)
                raise

    def command(self, args, log_path, output=None, timeout=None):
        # Tool stdout/stderr stay live in the stage log (including long addNorm jobs).
        return self.pipeline([args], output, log_path, timeout)

    def warn(self, message):
        with self.warning_lock:
            self.warnings.append(message)
        log("WARNING: " + message)


def write_manifests(plan):
    out = Path(plan.get("metadata_dir", str(Path(plan["roots"]["out"]) / "run_metadata")))
    json_write(out / "run_plan.json", plan)
    rows, tech, merges, units = [], [], [], []
    for name, sample in plan["samples"].items():
        entity = plan["entities"][f"sample:{name}"]
        ref = plan["references"][entity["reference_id"]]
        rows.append({**sample, "reference_id": entity["reference_id"], "species": ref["species"], "assembly": ref["assembly"],
                     "browser_preset": ref["browser_preset"], "sample_dir": entity["paths"]["entity_dir"],
                     "pairs_dir": entity["paths"]["pairs_dir"], "matrix_dir": entity["paths"]["matrix_dir"],
                     "merge_group": sample["resolved_merge_group"], "workflow_mode": plan.get("workflow_mode", "fastq"),
                     "input_stage": entity.get("input_stage", "fastq"), "deduplicated_pairs": sample.get("deduplicated_pairs", ""),
                     "input_members": ";".join(sample.get("input_members", [name]))})
        tech.append({"merge_group": sample["resolved_merge_group"], "sample": name, "pairs": entity["paths"]["accepted"],
                     "reference_id": entity["reference_id"], "dedup_unit": sample["unit"]})
    for name, key in plan["pools"].items():
        entity = plan["entities"][key]
        ref = plan["references"][entity["reference_id"]]
        merges.append({"merge_group": name, "reference_id": entity["reference_id"], "assembly": ref["assembly"],
                       "browser_preset": ref["browser_preset"], "members": ",".join(entity["members"]),
                       "dedup_units": ",".join(entity["units"]), "merge_dir": entity["paths"]["entity_dir"],
                       "pairs_dir": entity["paths"]["pairs_dir"], "matrix_dir": entity["paths"]["matrix_dir"], "pool_type": entity["pool_type"]})
    tsv_write(out / "sample_manifest.tsv", rows, ["sample", "assay", "reference_id", "species", "assembly", "browser_preset", "condition",
        "biological_replicate", "technical_replicate", "merge_group", "dedup_group", "hic_protocol", "restriction_enzymes", "fastq_r1", "fastq_r2", "sample_dir", "pairs_dir", "matrix_dir",
        "workflow_mode", "input_stage", "deduplicated_pairs", "input_members"])
    tsv_write(out / "technical_replicates.tsv", tech, ["merge_group", "sample", "pairs", "reference_id", "dedup_unit"])
    tsv_write(out / "merge_manifest.tsv", merges, ["merge_group", "reference_id", "assembly", "browser_preset", "members", "dedup_units", "merge_dir", "pairs_dir", "matrix_dir", "pool_type"])
    for key, unit in plan["dedup_units"].items():
        units.append({"dedup_unit": key, "members": ",".join(unit["members"]), "reference_id": unit["reference_id"], "profile": unit["policy"]["profile"]})
    tsv_write(out / "dedup_units.tsv", units, ["dedup_unit", "members", "reference_id", "profile"])


def downstream(rt, entity, accepted_key):
    if not truth(rt.c["RUN_PRELIM_DOWNSTREAM"]):
        return True
    paths, store = entity["paths"], rt.store(entity)
    if not truth(rt.c["RUN_MCOOL"]):
        rt.warn(f"{entity['id']}: downstream skipped because RUN_MCOOL=false")
        return True
    if not shutil.which("cooltools"):
        rt.warn(f"{entity['id']}: downstream unavailable (cooltools missing)")
        return False
    parameters = {key: value for key, value in rt.c.items() if key.startswith(("TAD_", "INSULATION_", "COMPARTMENT_", "LOOP_", "MCOOL_"))}
    ref = rt.plan["references"][entity["reference_id"]]
    parameters["phasing_track"] = identity(Path(ref["phasing_track"]), True) if ref["phasing_track"] else None
    key = rt.key("downstream", {"accepted": accepted_key, "parameters": parameters}, ["cooler", "cooltools"], ["pipeline.py", "run_downstream.sh"])
    dest = Path(paths["logs"]) / f"{entity['name']}.downstream_manifest.json"
    outputs = {"manifest": dest}
    if store.complete("downstream", key, outputs):
        previous = json.loads(dest.read_text())
        if all(Path(item["path"]).is_file() and identity(Path(item["path"])) == item for item in previous["artifacts"]):
            return True
        # A lost downstream file must not be certified by the old manifest.
        key = rt.key("downstream_missing", key)
    def build(t, work):
        local = work / "products"
        # Per-entity overrides are passed as command arguments through env(1), not global thread state.
        command = ["env", f"MICROC2TRACKS_REFERENCE_ID={entity['reference_id']}", f"MICROC2TRACKS_PHASING_TRACK_OVERRIDE={ref['phasing_track']}",
                   "bash", str(ROOT / "scripts/run_downstream.sh"), "-l", "-c", str(rt.config_path), "-s", entity["name"],
                   "-m", paths["norm_mcool"], "-o", str(local)]
        rt.command(command, rt.log(entity, "prelim_downstream"))
        old = json.loads(dest.read_text()) if dest.is_file() else {"artifacts": []}
        owned = {a["path"] for a in old["artifacts"]}
        artifacts = []
        for source in local.rglob("*"):
            if not source.is_file():
                continue
            target = Path(paths["entity_dir"]) / source.relative_to(local)
            if target.exists() and str(target.resolve()) not in owned:
                raise RuntimeError(f"Refusing to overwrite unowned downstream output {target}")
            artifacts.append(publish(source, target))
        json_write(t["manifest"], {"status": "success", "artifacts": artifacts, "signature": key})
    try:
        store.run("downstream", key, outputs, build)
        return True
    except Exception as exc:
        rt.warn(f"{entity['id']}: downstream failed: {exc}")
        return False


def run_entity(rt, processor, matrices, key):
    entity = rt.plan["entities"][key]
    locks = sorted({str(Path(entity["paths"][p]) / ".microc2tracks.lock") for p in ("entity_dir", "pairs_dir", "matrix_dir")})
    started, status = time.time(), "running"
    with contextlib.ExitStack() as stack:
        for path in locks:
            stack.enter_context(FileLock(Path(path)))
        def write_status(value, error=""):
            tsv_write(Path(entity["paths"]["status"]), [{"sample": entity["name"], "merge_group": entity["name"] if entity["kind"] == "pool" else "",
                "reference_id": entity["reference_id"], "assembly": rt.plan["references"][entity["reference_id"]]["assembly"],
                "status": value, "start_epoch": int(started), "end_epoch": int(time.time()) if value != "running" else "",
                "runtime_seconds": int(time.time()-started), "error": error}],
                ["sample", "merge_group", "reference_id", "assembly", "status", "start_epoch", "end_epoch", "runtime_seconds", "error"])
        write_status(status)
        log(f"{key}: processing {len(entity['members'])} source member(s); profile={entity['policy']['profile']}")
        try:
            processor.ensure_accepted(key)
            processor.stats(key)
            build = entity.get("build_matrices") or entity["kind"] == "pool" or (entity["kind"] == "sample" and truth(rt.c["BUILD_SAMPLE_MATRICES"]))
            if build:
                matrices.build(entity, processor.keys[key]["accepted"])
                if not downstream(rt, entity, processor.keys[key]["accepted"]):
                    status = "partial_success"
            if truth(rt.c["RUN_MULTIQC"]) and entity["kind"] == "sample" and shutil.which("multiqc"):
                try:
                    rt.command(["multiqc", "-f", entity["paths"]["entity_dir"], entity["paths"]["pairs_dir"], "-n", f"{entity['name']}.multiqc.html", "-o", entity["paths"]["qc_dir"]], rt.log(entity, "multiqc"))
                except Exception as exc:
                    rt.warn(f"{key}: MultiQC failed: {exc}")
                    status = "partial_success"
            write_status("success" if status == "running" else status)
            log(f"{key}: finished")
        except BaseException as exc:
            print(f"ERROR: {key}: {exc}", file=sys.stderr, flush=True)
            with contextlib.suppress(OSError):
                write_status("interrupted" if isinstance(exc, (KeyboardInterrupt, InterruptedError)) else "failed", str(exc))
            raise


def batch(rt, keys, processor, matrices, limit):
    with concurrent.futures.ThreadPoolExecutor(max_workers=limit) as pool:
        futures = [pool.submit(run_entity, rt, processor, matrices, key) for key in keys]
        try:
            for future in concurrent.futures.as_completed(futures):
                future.result()
        except BaseException:
            rt.stop.set()
            for future in futures:
                future.cancel()
            raise


def cleanup(rt):
    for entity in rt.plan["entities"].values():
        if entity["kind"] not in {"sample", "dedup", "pool"}:
            continue
        store, paths = rt.store(entity), entity["paths"]
        remove = []
        for flag, labels in (("KEEP_SORTED_PAIRS", ["sorted"]), ("KEEP_DEDUP_PAIRS", ["dedup"]),
                            ("KEEP_TRIMMED_FASTQ", ["trim_r1", "trim_r2"]), ("KEEP_JUICER_PAIRS", ["juicer"]),
                            ("KEEP_RAW_MCOOL", ["raw_mcool"]), ("KEEP_RAW_HIC", ["raw_hic"])):
            if not truth(rt.c[flag]):
                remove += labels
        if not truth(rt.c["KEEP_SINGLE_RES_COOL"]):
            remove += ["raw_cool"] + (["norm_cool"] if truth(rt.c["RUN_MCOOL"]) else [])
        for label in remove:
            store.remove_owned(Path(paths[label]))
    # A shared temp root is never recursively removed. Remove only empty directories
    # inside this UUID run, and the explicitly generated owner marker.
    if truth(rt.c["CLEAN_TMP_ON_SUCCESS"]):
        (rt.scratch / "owner.json").unlink(missing_ok=True)
        for directory in sorted((p for p in rt.scratch.rglob("*") if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
            with contextlib.suppress(OSError):
                directory.rmdir()
        with contextlib.suppress(OSError):
            rt.scratch.rmdir()


def execute(plan, config):
    out = Path(plan["roots"]["out"])
    guard = storage_guard(plan)
    guard.check()  # before mkdir: an absent mount must not receive an output directory
    with FileLock(out / "run_metadata/.run.lock"):
        validate_existing_layout(plan)
        rt = Runtime(plan, config)
        plan["tool_versions"] = rt.versions
        plan["source_hashes"] = rt.code_hashes
        plan["source_config"] = str(config.resolve())
        write_manifests(plan)
        meta_dir = Path(plan.get("metadata_dir", str(out / "run_metadata")))
        json_write(meta_dir / "tool_versions.json", rt.versions)
        json_write(meta_dir / "source_hashes.json", rt.code_hashes)
        # Downstream receives the fully resolved defaults, not a partial/CRLF source
        # config. This file is generated metadata; the user's config is not edited.
        effective_config = meta_dir / "effective_config.conf"
        effective_config.write_text("".join(f"{key}={shlex.quote(value)}\n" for key, value in sorted(rt.c.items())), encoding="utf-8", newline="\n")
        rt.config_path = effective_config
        processor, matrices = PairProcessor(rt), MatrixProducts(rt)
        metadata = Path(plan.get("metadata_dir", str(out / "run_metadata")))
        report_dir = metadata / "final_report" if plan.get("manual") else out / "final_report"
        state = metadata / "run_status.json"
        json_write(state, {"status": "running", "run_id": plan["run_id"]})
        previous_handler = signal.getsignal(signal.SIGTERM)
        def interrupted(signum, frame):
            rt.stop.set()
            raise KeyboardInterrupt("Terminated; only this run's child processes will be stopped")
        signal.signal(signal.SIGTERM, interrupted)
        try:
            batch(rt, [k for k, e in plan["entities"].items() if e["kind"] == "sample"], processor, matrices, int(rt.c["MAX_PARALLEL_SAMPLES"]))
            batch(rt, [k for k, e in plan["entities"].items() if e["kind"] == "dedup"], processor, matrices, int(rt.c["MAX_PARALLEL_MERGES"]))
            batch(rt, list(plan["pools"].values()), processor, matrices, int(rt.c["MAX_PARALLEL_MERGES"]))
            if truth(rt.c["RUN_GLOBAL_MULTIQC"]) and shutil.which("multiqc"):
                try:
                    rt.command(["multiqc", "-f", str(out), "-n", "microc2tracks_multiqc.html", "-o", str(report_dir)], report_dir / "microc2tracks_multiqc.log")
                except Exception as exc:
                    rt.warn(f"Global MultiQC failed: {exc}")
            from summarize_run import write_plan_summary
            if truth(rt.c["RUN_FINAL_REPORT"]):
                write_plan_summary(plan, report_dir)
            cleanup(rt)
            export_artifacts(plan)
            json_write(state, {"status": "partial_success" if rt.warnings else "success", "run_id": plan["run_id"], "warnings": rt.warnings})
            return 2 if rt.warnings else 0
        except BaseException as exc:
            rt.stop.set()
            print(f"ERROR: preserving original failure: {exc}", file=sys.stderr, flush=True)
            with contextlib.suppress(Exception):
                export_artifacts(plan)
            if truth(rt.c["RUN_FINAL_REPORT"]):
                with contextlib.suppress(Exception):
                    from summarize_run import write_plan_summary
                    write_plan_summary(plan, report_dir)
            with contextlib.suppress(OSError):
                json_write(state, {"status": "interrupted" if isinstance(exc, (KeyboardInterrupt, InterruptedError)) else "failed", "error": str(exc), "scratch": str(rt.scratch)})
            raise
        finally:
            signal.signal(signal.SIGTERM, previous_handler)


def export_artifacts(plan):
    rows = []
    for key, entity in plan["entities"].items():
        if entity["kind"] in {"external", "external_sorted"}:
            continue
        for record_path in Path(entity["paths"]["done"]).glob("*.json"):
            record = json.loads(record_path.read_text())
            for label, item in record.get("outputs", {}).items():
                path = Path(item["path"])
                rows.append({"entity": key, "stage": record["stage"], "category": label, "path": item["path"],
                             "size_bytes": path.stat().st_size if path.exists() else "", "retained": path.exists(),
                             "stage_status": record["status"], "signature": record["signature"]})
    tsv_write(Path(plan.get("metadata_dir", str(Path(plan["roots"]["out"]) / "run_metadata"))) / "artifacts.tsv", rows,
              ["entity", "stage", "category", "path", "size_bytes", "retained", "stage_status", "signature"])


def manual_plan(c, args):
    """Standalone merges require declared input modes; never guess stage from filenames."""
    name = safe_id(args.name)
    sources = [Path(p).resolve() for p in args.pairs.split(",")]
    if not sources or len(set(sources)) != len(sources) or any(not p.is_file() for p in sources):
        raise ValueError("Provide distinct existing input pair files")
    if args.mode == "pool":
        proofs = validate_pool_inputs(sources)
    elif args.mode == "joint-dedup":
        proofs = [read_proof(p) for p in sources]
        for source, proof in zip(sources, proofs):
            verify_proof_payload(source, proof)
        if any(p.get("stage") != "sorted" or not p.get("pre_dedup") or len(p.get("members", [])) != 1 for p in proofs):
            raise ValueError("joint-dedup requires original per-row sorted pre-dedup streams with provenance, not deduplicated outputs")
        if len({signature({"reference": p["reference_identity"], "assay": p["policy"]["assay"]}) for p in proofs}) != 1:
            raise ValueError("Incompatible reference/assay in joint-dedup inputs")
        if len({signature(p["parse_settings"]) for p in proofs}) != 1:
            raise ValueError("Incompatible parse settings in joint-dedup inputs")
        if len({m for p in proofs for m in p["members"]}) != len(proofs):
            raise ValueError("Repeated source row in joint-dedup inputs")
    else:
        proofs = []
    registry = Registry(Path(c["REFERENCE_REGISTRY"]))
    requested = args.reference or (proofs[0]["reference_id"] if proofs else c["DEFAULT_REFERENCE_ID"])
    ref = registry.resolve(requested)
    legacy_override = ref.reference_id == registry.resolve(c["DEFAULT_REFERENCE_ID"]).reference_id and c["REFERENCE_FASTA"] and c["CHROM_SIZES"]
    if legacy_override:
        ref = replace(ref, fasta=c["REFERENCE_FASTA"], chrom_sizes=c["CHROM_SIZES"], canonical_regex=c["CANONICAL_CHROMS_REGEX"] or ref.canonical_regex)
    base = Path.cwd() if legacy_override else registry.path.resolve().parent
    ref = replace(ref, fasta=absolute(ref.fasta, base), chrom_sizes=absolute(ref.chrom_sizes, base),
                  phasing_track=absolute(ref.phasing_track, base) if ref.phasing_track else "")
    reference = asdict(ref)
    reference["identity"] = {"fasta": identity(Path(ref.fasta), True), "chrom_sizes": identity(Path(ref.chrom_sizes), True)}
    from reference_registry import read_sizes
    if read_sizes(Path(ref.fasta + ".fai"), "FASTA index") != read_sizes(Path(ref.chrom_sizes), "chromosome sizes"):
        raise ValueError("Configured FASTA index and chromosome sizes disagree")
    from validate_pairs_reference import declared_assembly
    if any(declared_assembly(p) != ref.assembly for p in sources):
        raise ValueError("Pair headers do not declare the selected reference assembly")
    if proofs and any(p["reference_id"] != ref.reference_id or any(p["reference_identity"][key]["sha256"] != reference["identity"][key]["sha256"] for key in ("fasta", "chrom_sizes")) for p in proofs):
        raise ValueError("Input provenance does not match configured reference")
    if args.mode == "joint-dedup":
        policy = resolve_policy({"assay": proofs[0]["policy"]["assay"]}, c)
    elif proofs:
        policy = proofs[0]["policy"].copy()  # pooling inherits verified upstream selection
    else:
        policy = {"assay": "unknown", "profile": "legacy_unverified", "enzymes": [], "enzyme_definitions": {},
                  "mapq": 0, "min_cis_dist": 0, "remove_same_fragment": False, "remove_adjacent_inward": False,
                  "filter_canonical": truth(c["FILTER_CANONICAL_CHROMS"]), "upstream_policy": "unknown"}
    if args.mode != "pool":
        policy["canonical_regex"] = ref.canonical_regex if truth(c["FILTER_CANONICAL_CHROMS"]) else ".*"
    roots, entities, units = resolve_roots(c), {}, {}
    members = [m for p in proofs for m in p["members"]] if proofs else [f"legacy_input_{i}" for i in range(len(sources))]
    source_keys = []
    for i, source in enumerate(sources):
        label = proofs[i]["members"][0] if args.mode == "joint-dedup" else f"source_{i}"
        key = f"sample:{label}" if args.mode == "joint-dedup" else f"external:{i}"
        source_keys.append(key)
        entities[key] = {"kind": "external_sorted" if args.mode == "joint-dedup" else "external", "id": key, "name": label,
            "source_signature": proofs[i]["signature"] if proofs else signature(identity(source)),
            "paths": {"sorted": str(source), "accepted": str(source), "filter_qc": proofs[i].get("qc", "") if proofs else ""}}
    key = f"dedup:{name}" if args.mode == "joint-dedup" else f"pool:{name}"
    entities[key] = {"id": key, "kind": "dedup" if args.mode == "joint-dedup" else "pool", "name": name, "members": sorted(members),
        "reference_id": ref.reference_id, "policy": policy, "units": [key] if args.mode == "joint-dedup" else source_keys,
        "dedup_group": name if args.mode == "joint-dedup" else "", "pool_type": args.mode,
        "legacy_pool": args.mode == "legacy-pool", "build_matrices": True,
        "paths": artifact_paths(roots, f"merged/{name}", name, c)}
    if args.mode == "pool":
        entities[key]["provenance_units"] = sorted({unit for proof in proofs for unit in proof["dedup_units"]})
        entities[key]["import_sources"] = [source for proof in proofs for source in proof.get("import_sources", [])]
    return {"schema_version": 1, "workflow_version": VERSION, "settings": c, "roots": roots,
            "metadata_dir": str(Path(roots["out"]) / "run_metadata/manual_merges" / name), "references": {ref.reference_id: reference},
            "samples": {}, "entities": entities, "dedup_units": units, "pools": {name: key} if args.mode != "joint-dedup" else {}, "warnings": [], "manual": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "preflight", "merge"):
        p = sub.add_parser(command)
        p.add_argument("--version", action="version", version=f"microc2tracks {VERSION}")
        p.add_argument("-c", "--config", type=Path, required=True)
        if command != "merge":
            p.add_argument("-s", "--samplesheet", type=Path, required=True)
            p.add_argument("--mode", choices=["fastq", "refilter-deduplicated"], help="Override WORKFLOW_MODE; refilter skips alignment and deduplication")
        else:
            p.add_argument("-n", "--name", required=True)
            p.add_argument("-r", "--reference", default="")
            p.add_argument("-p", "--pairs", required=True, help="Comma-separated paths")
            p.add_argument("--mode", choices=["pool", "joint-dedup", "legacy-pool"], default="pool")
    args = parser.parse_args()
    try:
        c = load_config(args.config)
        plan = manual_plan(c, args) if args.command == "merge" else build_plan(c, args.samplesheet, mode=args.mode)
        for warning in plan["warnings"]:
            log("WARNING: " + warning)
        if args.command == "preflight":
            validate_existing_layout(plan)
            print(json.dumps({"workflow_mode": plan.get("workflow_mode", "merge"), "storage": storage_guard(plan).check(), "samples": len(plan["samples"]),
                "dedup_units": plan["dedup_units"], "pools": plan["pools"], "tool_versions": tool_versions(plan)}, indent=2))
            log("Preflight successful; no analysis launched or output files adopted")
            return 0
        return execute(plan, args.config)
    except (Exception, KeyboardInterrupt) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 130 if isinstance(exc, KeyboardInterrupt) else 1


if __name__ == "__main__":
    raise SystemExit(main())
