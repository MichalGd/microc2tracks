"""Stage graph for ordinary runs, refilter imports and standalone merges."""
from __future__ import annotations

import contextlib
from pathlib import Path
import shutil
import time

from stage_scheduler import StageScheduler, Task, reservation, resource_profile, size_gib
from workflow_common import builds_hic, builds_matrices, json_write, truth, tsv_write


def write_status(rt, entity, status, error=""):
    started = entity.setdefault("start_epoch", int(time.time()))
    tsv_write(Path(entity["paths"]["status"]), [{"sample": entity["name"],
        "merge_group": entity["name"] if entity["kind"] == "pool" else "", "reference_id": entity["reference_id"],
        "assembly": rt.plan["references"][entity["reference_id"]]["assembly"], "status": status,
        "start_epoch": started, "end_epoch": "" if status == "running" else int(time.time()),
        "runtime_seconds": int(time.time())-started, "error": error}],
        ["sample", "merge_group", "reference_id", "assembly", "status", "start_epoch", "end_epoch", "runtime_seconds", "error"])


def graph(rt, processor, matrices):
    from pipeline import downstream
    plan, c, tasks = rt.plan, rt.c, {}
    entities = {key: e for key, e in plan["entities"].items() if e["kind"] not in {"external", "external_sorted"}}

    def add(key, stage, kind, action, dependencies=(), priority=10, cached=None, prepares="", consumes=""):
        cpus, memory = reservation(kind, c, plan)
        if kind == "filter" and plan.get("workflow_mode") == "refilter-deduplicated":
            cpus += int(c["THREADS_SORT"]) + 2
            memory += size_gib(c["PAIRTOOLS_SORT_MEMORY"])
        task = Task(f"{key}/{stage}", kind, key, action, set(dependencies), cpus, memory, priority, prepares, consumes, cached)
        e = plan["entities"][key]
        inputs = [Path(e["paths"][label]) for label in ("sorted", "dedup", "accepted") if label in e["paths"] and Path(e["paths"][label]).is_file()]
        if e["kind"] in {"pool", "dedup"}:
            source_keys = e["units"] if e["kind"] == "pool" else [f"sample:{name}" for name in e["members"]]
            for unit in source_keys:
                source = plan["entities"][unit]["paths"]
                inputs += [Path(source[label]) for label in ("sorted", "accepted") if label in source and Path(source[label]).is_file()]
        raw_bytes = sum(item.get("size", 0) for name in e.get("source_rows", e["members"])
                        for item in plan["samples"].get(name, {}).get("fastq_identities", {}).values())
        amount = max(raw_bytes, sum(p.stat().st_size for p in set(inputs)))
        task.prepared_bytes = int(raw_bytes * float(c["OUTPUT_ESTIMATE_FACTOR"]))
        scratch = int(float(c["STAGE_SCRATCH_RESERVATION_GB"]) * 1024**3)
        fan_in = len(e["units"]) if e["kind"] == "pool" else len(e["members"])
        if kind == "alignment" or (kind in {"merge", "joint_dedup"} and fan_in > int(c["MERGE_FAN_IN"])):
            scratch = max(scratch, int(amount * float(c["SCRATCH_ESTIMATE_FACTOR"])))
        output = max(int(float(c["STAGE_OUTPUT_RESERVATION_GB"]) * 1024**3), int(amount * float(c["OUTPUT_ESTIMATE_FACTOR"]))) if kind != "report" else 0
        scratch_root = plan["roots"]["scratch"]
        destination = plan["roots"]["trimmed" if kind == "fastp" else "matrices" if kind in {"matrix", "hic"} else "pairs"]
        task.disk = {scratch_root: scratch}
        task.disk[destination] = task.disk.get(destination, 0) + output
        if c["LARGE_OUTPUT_STAGING"] == "scratch":
            task.disk[scratch_root] += output  # publication source plus destination copy
        if kind == "hic" and c["PUBLIC_HIC_DIR"]:
            publication = str(Path(c["PUBLIC_HIC_DIR"]).resolve())
            task.disk[publication] = task.disk.get(publication, 0) + output
        budget = float(c["SCRATCH_BUDGET_GB"]) * 1024**3
        if budget and task.disk[scratch_root] > budget and not (cached and cached()):
            raise ValueError(f"{task.id} scratch estimate exceeds SCRATCH_BUDGET_GB; choose larger/explicit HDD scratch or revise the estimate")
        tasks[task.id] = task
        return task.id

    def known(key, stage):
        return processor.reusable(key, stage)

    # A retained downstream checkpoint can remove prerequisites from the work
    # graph. In particular, cleanup of disposable trims/sorted pairs must not
    # force a successful library to be mapped again on resume.
    needed_rows = set()
    for key, entity in entities.items():
        if entity["kind"] == "pool" or entity.get("input_stage") == "deduplicated":
            continue
        if entity["kind"] == "sample" and not entity.get("contact_processing", True):
            continue
        if processor.available(key, "accepted") or processor.available(key, "dedup"):
            continue
        if entity["kind"] == "sample":
            needed_rows.add(key)
        elif not processor.available(key, "sorted"):
            needed_rows.update(f"sample:{name}" for name in entity["members"])

    for key, entity in entities.items():
        if entity["kind"] != "sample" or entity.get("input_stage") == "deduplicated":
            continue
        needed = key in needed_rows
        trim = add(key, "trim", "fastp", lambda k=key, n=needed: processor.ensure_trim(k) if n else None,
            priority=1, cached=lambda k=key, n=needed: not n or known(k, "trim"), prepares=key)
        def alignment(k=key, n=needed):
            if n:
                processor.ensure_sorted(k)
                # Only newly generated trims are released early. Original RC4
                # outputs are never deleted by changing release/configuration.
                if not truth(c["KEEP_TRIMMED_FASTQ"]):
                    record = rt.store(plan["entities"][k]).read("trim")
                    if record.get("run_id") == plan["run_id"]:
                        for label in ("trim_r1", "trim_r2"):
                            rt.store(plan["entities"][k]).remove_owned(Path(plan["entities"][k]["paths"][label]))
        add(key, "sorted", "alignment", alignment, [trim], priority=0,
            cached=lambda k=key, n=needed: not n or known(k, "sorted"), consumes=key)

    # Construct contact tasks independently of declaration order; explicit
    # dependencies, not worker recursion, control admission.
    for key, entity in entities.items():
        if entity["kind"] == "pool" or not entity.get("contact_processing", True):
            continue
        imported = entity.get("input_stage") == "deduplicated"
        stream = not imported and processor.streaming(key)
        before = []
        if not imported:
            sorted_id = f"{key}/sorted"
            if entity["kind"] == "dedup":
                dependencies = [f"sample:{name}/sorted" for name in entity["members"] if f"sample:{name}/sorted" in tasks]
                if not stream:
                    add(key, "sorted", "merge", lambda k=key: None if processor.available(k, "accepted") or processor.available(k, "dedup") else processor.ensure_sorted(k),
                        dependencies, priority=4, cached=lambda k=key: processor.available(k, "accepted") or processor.available(k, "dedup") or known(k, "sorted"))
            before = dependencies if stream else [sorted_id] if sorted_id in tasks else []
            def deduplicate(k=key):
                if processor.available(k, "accepted"):
                    return
                processor.ensure_dedup(k)
                e = entities[k]
                if e["kind"] == "dedup" and not truth(c["KEEP_RAW_MERGED_PAIRS"]):
                    record = rt.store(e).read("sorted")
                    if record.get("run_id") == plan["run_id"]:
                        rt.store(e).remove_owned(Path(e["paths"]["sorted"]))
            dedup = add(key, "dedup", "joint_dedup" if stream else "dedup", deduplicate,
                before, priority=5, cached=lambda k=key: processor.available(k, "accepted") or known(k, "dedup"))
            before = [dedup]
        add(key, "accepted", "filter", lambda k=key: processor.ensure_accepted(k), before, priority=6,
            cached=lambda k=key: known(k, "accepted"))

    for key, entity in entities.items():
        if entity["kind"] == "pool":
            dependencies = [f"{unit}/accepted" for unit in entity["units"] if unit in entities]
            if entity.get("alias_of"):
                dependencies.append(f"{entity['alias_of']}/stats")
            dependencies += [f"{parent}/finished" for parent in entity.get("biological_dependencies", [])]
            add(key, "accepted", "merge", lambda k=key: processor.ensure_accepted(k), dependencies, priority=7,
                cached=(lambda: True) if entity.get("alias_of") else lambda k=key: known(k, "accepted"))

    for key, entity in entities.items():
        before = [f"{key}/accepted"] if entity.get("contact_processing", True) else [f"{key}/sorted"]
        if entity.get("contact_processing", True):
            stats_before = list(before)
            if entity["kind"] == "pool" and not entity.get("legacy_pool") and c["POOL_STATS_MODE"] == "merge":
                stats_before += [f"{unit}/stats" for unit in entity["units"] if unit in entities]
            stats = add(key, "stats", "report", lambda k=key: None if entities[k].get("alias_of") else processor.stats(k), stats_before, priority=8)
            # Full descriptive statistics may overlap product generation, but
            # remain required before this entity can be declared finished.
            required_stats = [stats]
        else:
            required_stats = []
        if builds_matrices(entity, c):
            matrix = add(key, "matrix", "matrix", lambda k=key: matrices.build(entities[k], processor.keys[k]["accepted"], kind="matrix"), before, priority=10)
            before = [matrix]
            if truth(c["RUN_PRELIM_DOWNSTREAM"]):
                def analyze(k=key):
                    if not downstream(rt, entities[k], processor.keys[k]["accepted"]):
                        entities[k]["partial_success"] = True
                before = [add(key, "downstream", "downstream", analyze, before, priority=12)]
            if builds_hic(entity, c):
                hic = add(key, "hic", "hic", lambda k=key: matrices.build(entities[k], processor.keys[k]["accepted"], kind="hic"),
                          [f"{key}/accepted"], priority=11)
                before.append(hic)
        def finish(k=key):
            e = entities[k]
            if truth(c["RUN_MULTIQC"]) and e["kind"] == "sample" and shutil.which("multiqc"):
                try:
                    rt.command(["multiqc", "-f", e["paths"]["entity_dir"], e["paths"]["pairs_dir"], "-n", f"{e['name']}.multiqc.html",
                                "-o", e["paths"]["qc_dir"]], rt.log(e, "multiqc"))
                except Exception as exc:
                    rt.warn(f"{k}: MultiQC failed: {exc}")
                    e["partial_success"] = True
            write_status(rt, e, "partial_success" if e.get("partial_success") else "success")
        add(key, "finished", "report", finish, [*before, *required_stats], priority=15)
    return list(tasks.values())


def execute_stages(rt, processor, matrices):
    metadata = Path(rt.plan.get("metadata_dir", str(Path(rt.plan["roots"]["out"]) / "run_metadata")))
    events = metadata / "stage_events.jsonl"
    states, started = {}, set()
    def notify(event):
        states[event["task"]] = event
        with events.open("a", encoding="utf-8") as handle:
            import json
            handle.write(json.dumps({"run_id": rt.plan["run_id"], **event}) + "\n")
        json_write(metadata / "stage_status.json", {"run_id": rt.plan["run_id"], "tasks": states})
        entity = rt.plan["entities"][event["entity"]]
        if event["state"] in {"running", "reused"} and event["entity"] not in started:
            started.add(event["entity"])
            write_status(rt, entity, "running")
        if event["state"] in {"failed", "interrupted"}:
            with contextlib.suppress(OSError):
                write_status(rt, entity, event["state"], f"Stage {event['task']} did not finish")
        if event["state"] in {"running", "complete", "reused", "failed"}:
            from pipeline import log
            log(f"{event['task']}: {event['state']} ({event['cpus']} CPU slots; {event['memory_gib']:.1f} GiB reservation)")
    profile = resource_profile(rt.c)
    profile["max_prepared_bytes"] = int(float(rt.c["MAX_PREPARED_GB"]) * 1024**3)
    rt.plan["resource_profile"] = profile
    scheduler = StageScheduler(profile, rt.stop, notify, rt.guard)
    tasks = graph(rt, processor, matrices)
    scheduler.validate(tasks)  # fail before any scientific stage starts
    json_write(metadata / "resource_plan.json", {**profile, "reservations": {
        kind: {"cpu_slots": reservation(kind, rt.c, rt.plan)[0], "memory_gib": reservation(kind, rt.c, rt.plan)[1]}
        for kind in profile["limits"]}, "rc4_compatibility_enabled": bool(rt.rc4_hashes),
        "tasks": [{"task": t.id, "kind": t.kind, "dependencies": sorted(t.dependencies),
                   "cpu_slots": t.cpus, "memory_gib": t.memory} for t in tasks]})
    return scheduler.run(tasks)
