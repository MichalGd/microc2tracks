"""Fingerprinted stages, guarded publication and narrowly owned cleanup."""
from __future__ import annotations

import contextlib
import json
import os
import shutil
import time
import uuid
from pathlib import Path

from workflow_common import identity, json_write, sha256


class FileLock:
    """OS-released lock; an interrupted process cannot leave a stale busy flag."""
    def __init__(self, path: Path, blocking=False):
        self.path = path
        self.blocking = blocking

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("a+b")
        try:
            if os.name == "nt":
                import msvcrt
                if self.path.stat().st_size == 0:
                    self.handle.write(b" ")
                    self.handle.flush()
                self.handle.seek(0)
                msvcrt.locking(self.handle.fileno(), msvcrt.LK_LOCK if self.blocking else msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(self.handle, fcntl.LOCK_EX | (0 if self.blocking else fcntl.LOCK_NB))
        except OSError as exc:
            self.handle.close()
            raise RuntimeError(f"Target is already locked by another run: {self.path}") from exc
        return self

    def __exit__(self, *args):
        self.handle.close()


def nearest_existing(path: Path) -> Path:
    while not path.exists():
        if path.parent == path:
            raise ValueError(f"Cannot locate filesystem for {path}")
        path = path.parent
    return path


class StorageGuard:
    def __init__(self, roots: dict[str, str], mounts: list[str], output_gb: float, scratch_gb: float):
        self.roots, self.mounts = roots, [Path(p).resolve() for p in mounts]
        self.output_gb, self.scratch_gb = output_gb, scratch_gb
        self.devices: dict[str, int] = {}

    def check(self) -> dict:
        for mount in self.mounts:
            if not mount.is_mount():
                raise RuntimeError(f"Required storage mount is absent: {mount}; refusing to write to its underlying filesystem")
        report = {}
        for key, value in self.roots.items():
            target = Path(value)
            if target.exists() and not target.is_dir():
                raise RuntimeError(f"Storage root is not a directory: {value}")
            existing = nearest_existing(target)
            if not os.access(existing, os.W_OK | os.X_OK):
                raise RuntimeError(f"Storage root/parent is not writable: {value} (existing parent {existing})")
            device = existing.stat().st_dev
            if key in self.devices and self.devices[key] != device:
                raise RuntimeError(f"Filesystem identity changed for {key}: {value}")
            self.devices[key] = device
            free = shutil.disk_usage(existing).free
            reserve = self.scratch_gb if key == "scratch" else self.output_gb
            if free < reserve * 1024**3:
                raise RuntimeError(f"Insufficient free space for {key} at {value}: {free/1024**3:.2f} GiB; reserve={reserve} GiB")
            report[key] = {"path": value, "device": device, "free_bytes": free, "reserve_gib": reserve}
        return report


def same_filesystem(source: Path, directory: Path) -> bool:
    return source.stat().st_dev == directory.stat().st_dev


def publish(source: Path, destination: Path) -> dict:
    """The source must be stage-owned. Keep it until the destination is verified."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r+b" if os.name == "nt" else "rb") as handle:
        os.fsync(handle.fileno())
    if same_filesystem(source, destination.parent):
        os.replace(source, destination)
    else:
        partial = destination.with_name(f".partial.{uuid.uuid4().hex}.{destination.name}")
        expected = sha256(source)
        shutil.copyfile(source, partial)
        if partial.stat().st_size != source.stat().st_size or sha256(partial) != expected:
            raise RuntimeError(f"Publication checksum mismatch: {partial}; retained source: {source}")
        source_stat = source.stat()
        os.utime(partial, ns=(source_stat.st_atime_ns, source_stat.st_mtime_ns))
        with partial.open("r+b" if os.name == "nt" else "rb") as handle:
            os.fsync(handle.fileno())
        os.replace(partial, destination)
        source.unlink()
    if os.name != "nt":
        descriptor = os.open(destination.parent, os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    return identity(destination)


class StageStore:
    def __init__(self, metadata: Path, scratch: Path, roots: list[Path], staging: str, guard: StorageGuard):
        self.metadata, self.scratch = metadata, scratch
        self.roots = [p.resolve() for p in roots]
        self.staging, self.guard = staging, guard
        self.metadata.mkdir(parents=True, exist_ok=True)

    def record_path(self, stage: str) -> Path:
        return self.metadata / f"{stage}.json"

    def read(self, stage: str) -> dict:
        path = self.record_path(stage)
        return json.loads(path.read_text()) if path.is_file() else {}

    def complete(self, stage: str, key: str, outputs: dict[str, Path]) -> bool:
        record = self.read(stage)
        if record.get("signature") != key or record.get("status") != "success":
            return False
        for label, path in outputs.items():
            saved = record.get("outputs", {}).get(label, {})
            if not path.is_file() or identity(path) != saved:
                return False
        return True

    def ensure_reusable(self, stage: str, key: str, outputs: dict[str, Path]) -> bool:
        """A storage-only move can republish a verified managed stage before invoking parents."""
        if self.complete(stage, key, outputs):
            return True
        old = self.read(stage)
        if (old.get("status") == "success" and old.get("signature") == key
                and set(old.get("outputs", {})) == set(outputs)
                and all(Path(item["path"]).is_file() and identity(Path(item["path"])) == item
                        for item in old.get("outputs", {}).values())):
            def unexpected(t, work):
                raise RuntimeError("Reusable stage changed while acquiring its lock; retry the run")
            self.run(stage, key, outputs, unexpected)
            return True
        return False

    def run(self, stage: str, key: str, outputs: dict[str, Path], action, validate=None) -> dict:
        with contextlib.ExitStack() as locks:
            locks.enter_context(FileLock(self.metadata / f".{stage}.lock"))
            for path in sorted(set(outputs.values())):
                locks.enter_context(FileLock(path.parent / f".{path.name}.lock"))
            if self.complete(stage, key, outputs):
                return self.read(stage)
            old = self.read(stage)
            relocatable = (old.get("status") == "success" and old.get("signature") == key
                           and set(old.get("outputs", {})) == set(outputs)
                           and all(Path(item["path"]).is_file() and identity(Path(item["path"])) == item
                                   for item in old.get("outputs", {}).values()))
            owned = {item["path"] for item in old.get("outputs", {}).values()}
            for path in outputs.values():
                if path.exists() and str(path.resolve()) not in owned:
                    raise RuntimeError(f"Refusing to overwrite an unowned/legacy output: {path}. Use a new OUTDIR; originals are preserved.")
                path.parent.mkdir(parents=True, exist_ok=True)
            self.guard.check()
            token = uuid.uuid4().hex
            work = self.scratch / stage / token
            work.mkdir(parents=True)
            temporary = {
                label: (work / path.name if self.staging == "scratch" else path.with_name(f".partial.{token}.{path.name}"))
                for label, path in outputs.items()
            }
            # Explicit ownership list is the only cleanup authority.
            owner = {"stage": stage, "run_scratch": str(self.scratch), "files": [str(p) for p in temporary.values()]}
            json_write(work / "owner.json", owner)
            started = time.time()
            published = {}
            try:
                if relocatable:
                    mapping = {old["outputs"][label]["path"]: str(path.resolve()) for label, path in outputs.items()}
                    def relocate_metadata(value):
                        if isinstance(value, dict):
                            return {k: relocate_metadata(v) for k, v in value.items()}
                        if isinstance(value, list):
                            return [relocate_metadata(v) for v in value]
                        return mapping.get(value, value) if isinstance(value, str) else value
                    for label, target in temporary.items():
                        original = Path(old["outputs"][label]["path"])
                        shutil.copy2(original, target)
                        if sha256(original) != sha256(target):
                            raise RuntimeError(f"Relocation verification failed: {original} -> {target}")
                        if target.name.endswith(".provenance.json"):
                            json_write(target, relocate_metadata(json.loads(target.read_text())))
                else:
                    action(temporary, work)
                for path in temporary.values():
                    if not path.is_file() or path.stat().st_size == 0:
                        raise RuntimeError(f"Missing/empty stage output: {path}")
                if validate:
                    validate(temporary)
                self.guard.check()
                for label, path in outputs.items():
                    published[label] = publish(temporary[label], path)
                record = {"schema": 1, "signature": key, "stage": stage, "status": "success",
                          "started": started, "completed": time.time(), "outputs": published,
                          "storage_only_republication": relocatable}
                json_write(self.record_path(stage), record)
                # Only our now-empty stage directory and owner file are removed.
                (work / "owner.json").unlink()
                with contextlib.suppress(OSError):
                    work.rmdir()
                return record
            except BaseException as exc:
                # Preserve prior ownership even after partial publication or failure.
                record = {"schema": 1, "stage": stage, "signature": key,
                          "status": "interrupted" if isinstance(exc, (KeyboardInterrupt, InterruptedError)) else "failed",
                          "error": str(exc), "started": started, "outputs": old.get("outputs", {}),
                          "retained_work": str(work), "temporary": {k: str(v) for k, v in temporary.items()}}
                record["outputs"].update(published)
                # A published partial bundle remains owned but never marked successful.
                for label, path in outputs.items():
                    if path.is_file() and str(path.resolve()) in owned:
                        record["outputs"][label] = identity(path)
                with contextlib.suppress(OSError):
                    json_write(self.record_path(stage), record)
                raise

    def remove_owned(self, path: Path) -> bool:
        resolved = path.resolve()
        if path.is_symlink() or not any(resolved != root and resolved.is_relative_to(root) for root in self.roots):
            raise RuntimeError(f"Cleanup target is outside allowed roots or a symlink: {path}")
        for record_path in self.metadata.glob("*.json"):
            record = json.loads(record_path.read_text())
            for item in record.get("outputs", {}).values():
                if item.get("path") == str(resolved) and path.is_file() and identity(path) == item:
                    path.unlink()
                    return True
        return False
