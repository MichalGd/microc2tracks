"""Dependency-ready scheduling. Waiting tasks never occupy workers or credits."""
from __future__ import annotations

import concurrent.futures
from dataclasses import dataclass, field
import math
import os
from pathlib import Path
import re
import threading
import time
from typing import Callable


def size_gib(value: str) -> float:
    match = re.fullmatch(r"\s*([0-9.]+)\s*([KMGTPE]?)(?:i?B)?\s*", str(value), re.I)
    if not match:
        raise ValueError(f"Invalid memory size: {value!r}")
    return float(match[1]) * 1024 ** ("KMGTPE".index(match[2].upper()) - 2) if match[2] else float(match[1]) / 1024**3


def available_cpus() -> int:
    count = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    quota = Path("/sys/fs/cgroup/cpu.max")
    if quota.is_file():
        amount, period = quota.read_text().split()[:2]
        if amount != "max":
            count = min(count, max(1, math.floor(int(amount) / int(period))))
    return count


def available_memory_gib() -> float:
    if os.name == "nt":
        import ctypes
        class Memory(ctypes.Structure):
            _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong)] + [
                (name, ctypes.c_ulonglong) for name in ("total", "available", "page", "available_page", "virtual", "available_virtual", "extended")]
        value = Memory()
        value.length = ctypes.sizeof(value)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(value)):
            raise OSError("Cannot determine available RAM; set TOTAL_MEMORY_GB explicitly")
        total = value.available / 1024**3
    else:
        memory = {line.split(":")[0]: int(line.split()[1]) * 1024
                  for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith(("MemAvailable:", "MemTotal:"))}
        total = memory.get("MemAvailable", memory["MemTotal"]) / 1024**3
        limit, used = Path("/sys/fs/cgroup/memory.max"), Path("/sys/fs/cgroup/memory.current")
        if limit.is_file() and limit.read_text().strip() != "max":
            total = min(total, max(0, int(limit.read_text()) - int(used.read_text())) / 1024**3)
    return total


def resource_profile(c: dict) -> dict:
    cpus = int(c["TOTAL_CPU_BUDGET"]) or max(1, available_cpus() - 1)
    requested_memory = float(c["TOTAL_MEMORY_GB"])
    available_budget = max(0.1, available_memory_gib() * 0.8)
    memory = min(requested_memory, available_budget) if requested_memory else available_budget
    # Explicit budgets are submitting-user reservations, not a host-wide allocator.
    if int(c["TOTAL_CPU_BUDGET"]) and cpus > available_cpus():
        raise ValueError(f"TOTAL_CPU_BUDGET={cpus} exceeds available affinity/quota CPUs ({available_cpus()})")
    limits = {kind: int(c[key]) for kind, key in {
        "fastp": "MAX_PARALLEL_FASTP", "alignment": "MAX_PARALLEL_ALIGNMENTS", "merge": "MAX_PARALLEL_MERGES",
        "dedup": "MAX_PARALLEL_DEDUP", "filter": "MAX_PARALLEL_FILTER", "matrix": "MAX_PARALLEL_MATRIX",
        "joint_dedup": "MAX_PARALLEL_JOINT_DEDUP",
        "hic": "MAX_PARALLEL_HIC", "downstream": "MAX_PARALLEL_DOWNSTREAM", "report": "MAX_PARALLEL_REPORTS"}.items()}
    return {"cpu_budget": cpus, "memory_gib": memory, "limits": limits,
            "requested_memory_gib": requested_memory,
            "max_prepared": int(c["MAX_PREPARED_SAMPLES"]), "memory_headroom_gib": float(c["MIN_AVAILABLE_MEMORY_GB"])}


def reservation(kind: str, c: dict, plan: dict) -> tuple[int, float]:
    compression = int(c["THREADS_COMPRESS"])
    cpu = {"fastp": int(c["THREADS_FASTP"]) + 1,
           "alignment": int(c["THREADS_ALIGN"]) + int(c["THREADS_SORT"]) + compression + 4,
           "merge": int(c["THREADS_MERGE"]) + compression + 2,
           "dedup": compression + 3, "joint_dedup": int(c["THREADS_MERGE"]) + compression + 3,
           "filter": compression + 2,
           "matrix": max(int(c["THREADS_MATRIX"]), int(c["THREADS_BALANCE"])) + 2,
           "hic": int(c["THREADS_HIC_NORM"]) + 4,
           "downstream": int(c["COOLTOOLS_THREADS"]) + 2, "report": 3}[kind]
    memory = float(c[f"MEMORY_{'ALIGN' if kind == 'alignment' else kind.upper()}_GB"])
    if kind in {"merge", "dedup", "joint_dedup"}:
        readers = max((len(e.get("units", [])) if e.get("kind") == "pool" else len(e["members"]) for e in plan.get("entities", {}).values()
                       if e.get("kind") in {"dedup", "pool"}), default=1) if kind != "dedup" else 1
        cpu += min(readers, int(c["MERGE_READERS"]))
        if kind == "joint_dedup" and not memory:
            memory = (float(c["MEMORY_MERGE_GB"]) or size_gib(c["PAIRTOOLS_MERGE_MEMORY"]) + 1) + float(c["MEMORY_DEDUP_GB"]) + 4
    if not memory:
        if kind == "alignment":
            # Conservative index + sort allowance; site profiles should override
            # it using observed peaks. This is not a hard RSS limit.
            largest = max((sum(Path(ref["bwa_index_prefix"] + suffix).stat().st_size
                                  for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
                                  if Path(ref["bwa_index_prefix"] + suffix).is_file())
                           for ref in plan["references"].values()), default=0)
            memory = size_gib(c["PAIRTOOLS_SORT_MEMORY"]) + largest * 2 / 1024**3 + 1
        elif kind == "merge":
            memory = size_gib(c["PAIRTOOLS_MERGE_MEMORY"]) + 1
        elif kind == "hic":
            memory = size_gib(c["JAVA_HEAP"]) + 4
    return cpu, memory


@dataclass
class Task:
    id: str
    kind: str
    entity: str
    action: Callable
    dependencies: set[str] = field(default_factory=set)
    cpus: int = 1
    memory: float = 0
    priority: int = 10
    prepares: str = ""
    consumes: str = ""
    cached: Callable | None = None
    disk: dict = field(default_factory=dict)
    prepared_bytes: int = 0


class StageScheduler:
    def __init__(self, profile: dict, stop: threading.Event, notify: Callable | None = None, guard=None):
        self.profile, self.stop, self.notify = profile, stop, notify or (lambda event: None)
        self.used_cpu, self.used_memory = 0, 0.0
        self.active, self.prepared = {}, set()
        self.states = {}
        self.guard = guard
        self.prepared_sizes = {}

    def emit(self, task, state, reason="", **extra):
        previous = self.states.get(task.id, {})
        if previous.get("state") == state and previous.get("reason") == reason and not extra:
            return
        event = {"task": task.id, "entity": task.entity, "kind": task.kind, "state": state, "reason": reason,
                 "cpus": task.cpus, "memory_gib": task.memory, "time": time.time(),
                 "used_cpus": self.used_cpu, "used_memory_gib": self.used_memory, **extra}
        self.states[task.id] = event
        self.notify(event)

    def validate(self, tasks):
        by_id = {task.id: task for task in tasks}
        if len(by_id) != len(tasks):
            raise ValueError("Duplicate scheduler task ID")
        for task in tasks:
            if task.dependencies - by_id.keys():
                raise ValueError(f"Missing dependencies for {task.id}: {task.dependencies - by_id.keys()}")
            if task.cpus > self.profile["cpu_budget"] or task.memory > self.profile["memory_gib"]:
                raise ValueError(f"Task {task.id} requires {task.cpus} CPU slots / {task.memory:.1f} GiB, exceeding "
                                 f"budget {self.profile['cpu_budget']} / {self.profile['memory_gib']:.1f}. Adjust resources before starting.")
        remaining, complete = dict(by_id), set()
        while remaining:
            ready = [key for key, task in remaining.items() if task.dependencies <= complete]
            if not ready:
                raise ValueError("Cyclic scheduler dependencies")
            complete.update(ready)
            for key in ready:
                del remaining[key]
        return by_id

    def run(self, tasks):
        pending = self.validate(tasks)
        order = {key: i for i, key in enumerate(pending)}
        completed, running = set(), {}
        workers = max(1, min(len(tasks), sum(self.profile["limits"].values())))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            try:
                while pending or running:
                    if self.stop.is_set():
                        raise InterruptedError("Stage scheduler interrupted")
                    made_progress = False
                    waiting_for_ram = False
                    waiting_for_disk = False
                    for task in sorted(pending.values(), key=lambda t: (t.priority, order[t.id])):
                        if not task.dependencies <= completed:
                            self.emit(task, "queued", "dependencies")
                            continue
                        if task.cached and task.cached():
                            # A completed task needs neither a worker nor credits.
                            task.action()
                            completed.add(task.id)
                            made_progress = True
                            if task.prepares:
                                self.prepared.add(task.prepares)
                                self.prepared_sizes[task.prepares] = task.prepared_bytes
                            if task.consumes:
                                self.prepared.discard(task.consumes)
                                self.prepared_sizes.pop(task.consumes, None)
                            self.emit(task, "reused")
                            del pending[task.id]
                            continue
                        reason = ""
                        if self.active.get(task.kind, 0) >= self.profile["limits"][task.kind]:
                            reason = "stage concurrency"
                        elif self.used_cpu + task.cpus > self.profile["cpu_budget"]:
                            reason = "CPU budget"
                        elif self.used_memory + task.memory > self.profile["memory_gib"]:
                            reason = "RAM budget"
                        elif ("memory_headroom_gib" in self.profile and
                              task.memory + self.profile["memory_headroom_gib"] > available_memory_gib()):
                            reason = "available RAM"
                            waiting_for_ram = True
                        elif task.prepares and task.prepares not in self.prepared and len(self.prepared) >= self.profile["max_prepared"]:
                            reason = "prepared-input backlog"
                        elif (task.prepares and task.prepares not in self.prepared and self.profile.get("max_prepared_bytes", 0)
                              and sum(self.prepared_sizes.values()) + task.prepared_bytes > self.profile["max_prepared_bytes"]):
                            reason = "prepared-input byte budget"
                        elif self.guard and not self.guard.reserve(task.id, task.disk):
                            reason = "filesystem reservations"
                            waiting_for_disk = True
                        if reason:
                            self.emit(task, "queued", reason)
                            continue
                        self.used_cpu += task.cpus
                        self.used_memory += task.memory
                        self.active[task.kind] = self.active.get(task.kind, 0) + 1
                        if task.prepares:
                            self.prepared.add(task.prepares)
                            self.prepared_sizes[task.prepares] = task.prepared_bytes
                        self.emit(task, "running")
                        running[pool.submit(task.action)] = (task, time.monotonic())
                        del pending[task.id]
                    if not running:
                        if pending:
                            if made_progress:
                                continue
                            if waiting_for_ram or waiting_for_disk:
                                # Another user's process can temporarily consume
                                # RAM. Keep this distinct from a dependency deadlock.
                                self.stop.wait(0.2)
                                continue
                            raise RuntimeError("Ready tasks cannot make progress; inspect backlog and resource settings")
                        break
                    finished, _ = concurrent.futures.wait(running, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED)
                    for future in finished:
                        task, started = running.pop(future)
                        self.used_cpu -= task.cpus
                        self.used_memory = max(0, self.used_memory - task.memory)
                        self.active[task.kind] -= 1
                        if self.guard:
                            self.guard.release(task.id)
                        try:
                            future.result()
                        except BaseException:
                            self.emit(task, "failed", duration_seconds=time.monotonic()-started)
                            raise
                        completed.add(task.id)
                        if task.consumes:
                            self.prepared.discard(task.consumes)
                            self.prepared_sizes.pop(task.consumes, None)
                        self.emit(task, "complete", duration_seconds=time.monotonic()-started)
            except BaseException:
                self.stop.set()
                for future, (task, _) in running.items():
                    future.cancel()
                    self.emit(task, "interrupted")
                for task in pending.values():
                    self.emit(task, "cancelled")
                # Runtime subprocess loops observe stop and terminate their own
                # process groups. The executor waits before returning to caller.
                raise
            finally:
                self.used_cpu, self.used_memory = 0, 0.0
                self.active.clear()
                if self.guard:
                    for task in tasks:
                        self.guard.release(task.id)
        return self.states
