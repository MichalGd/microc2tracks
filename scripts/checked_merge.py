#!/usr/bin/env python3
"""Bounded GNU merge with explicitly waited, CRC-checked input readers.

Retains pairtools 1.1.3's comparator and header merger. FIFOs carry data,
not checkpoints. Every reader must commit a receipt before the stage succeeds.
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from workflow_common import identity, json_write


def binary_open(path):
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


@contextlib.contextmanager
def decoder_slot(directory=None, slots=2):
    """Limit active decompression, not FIFO lifetime (which would deadlock merges)."""
    if directory is None:
        yield
        return
    import fcntl
    while True:
        for i in range(slots):
            lock = (directory / f"decoder_{i}.lock").open("a+b")
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                lock.close()
                continue
            try:
                yield
            finally:
                lock.close()
            return
        time.sleep(0.005)


def reader(source, fifo, receipt, lock_directory=None, slots=2):
    before = identity(source)
    count, digest, last = 0, hashlib.sha256(), b""
    with binary_open(source) as inp, (fifo.open("wb") if fifo else contextlib.nullcontext(sys.stdout.buffer)) as out:
        with decoder_slot(lock_directory, slots):
            line = inp.readline()
        while line.startswith(b"#"):
            with decoder_slot(lock_directory, slots):
                line = inp.readline()
        while line:
            count += line.count(b"\n")
            digest.update(line)
            out.write(line)
            last = line[-1:]
            with decoder_slot(lock_directory, slots):
                line = inp.read(1024 * 1024)
        if last and last != b"\n":
            raise ValueError(f"Truncated/non-terminated pairs record: {source}")
        if identity(source) != before:
            raise ValueError(f"Input changed during merge: {source}")
        # Receipt precedes FIFO close/EOF, so downstream success cannot race it.
        json_write(receipt, {"input": before, "records": count, "data_sha256": digest.hexdigest()})


def headers(sources):
    from pairtools.lib import headerops
    result = []
    for source in sources:
        with binary_open(source) as inp:
            header = []
            for line in inp:
                if not line.startswith(b"#"):
                    break
                header.append(line.decode().rstrip("\n"))
        result.append(header)
    if not headerops.all_same_columns(result):
        raise ValueError("Input pairs cannot contain different columns")
    expected = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type".split()
    for header in result:
        if headerops.extract_column_names(header)[:8] != expected:
            raise ValueError("Merge needs canonical sorted-pairs column order")
    if len(sources) == 1:
        return result[0]
    return headerops.append_new_pg(headerops.merge_headers(result), ID="pairtools_merge", PN="pairtools_merge")


def verify_digest(source, record):
    import json
    proof = Path(str(source) + ".provenance.json")
    if proof.is_file():
        expected = json.loads(proof.read_text()).get("logical_data_sha256")
        if expected and expected != record["data_sha256"]:
            raise ValueError(f"Sorted pairs logical checksum changed: {source}")


def merge(sources, work, report, memory="8G", threads=4, fan_in=8, readers=2, max_inputs=256):
    if os.name == "nt":
        raise RuntimeError("Checked pairs merging requires Linux/POSIX FIFOs")
    if len(sources) > max_inputs:
        raise ValueError("Input count exceeds MAX_MERGE_INPUTS; pre-merge sorted inputs in bounded batches, then perform ONE joint dedup pass")
    work.mkdir(parents=True, exist_ok=True)
    header = headers(sources)
    sys.stdout.buffer.write(("\n".join(header) + "\n").encode())
    sys.stdout.buffer.flush()
    if len(sources) == 1:
        receipt = work / "single.receipt.json"
        reader(sources[0], None, receipt)
        import json
        record = json.loads(receipt.read_text())
        verify_digest(sources[0], record)
        json_write(report, {"schema": 1, "reader_receipts": [record], "records": record["records"], "comparator": "single-input"})
        return
    children, receipts, fifos = [], [], []
    previous = signal.getsignal(signal.SIGTERM)
    def interrupted(*_):
        raise InterruptedError("Checked merge interrupted")
    signal.signal(signal.SIGTERM, interrupted)
    try:
        for i, source in enumerate(sources):
            fifo, receipt = work / f"input_{i}.fifo", work / f"input_{i}.receipt.json"
            os.mkfifo(fifo, 0o600)
            fifos.append(fifo)
            receipts.append(receipt)
            children.append(subprocess.Popen([sys.executable, "-B", __file__, "--reader", str(source),
                "--fifo", str(fifo), "--report", str(receipt), "--reader-lock-dir", str(work), "--readers", str(readers)]))
        environment = {**os.environ, "LC_ALL": "C", "LC_COLLATE": "C", "LANG": "C"}
        # Exact pairtools comparator including GNU sort's full-line fallback.
        command = ["sort", "-k2,2", "-k4,4", "-k3,3n", "-k5,5n", "-k8,8", "--merge", "-t", "\t",
                   f"--parallel={threads}", f"--batch-size={fan_in}", "-S", memory,
                   f"--temporary-directory={work}", *map(str, fifos)]
        children.append(subprocess.Popen(command, stdout=sys.stdout.buffer, env=environment))
        while any(p.poll() is None for p in children):
            if any(p.poll() not in (None, 0) for p in children):
                raise RuntimeError("Merge input reader or GNU sort failed; output is incomplete")
            time.sleep(0.05)
        if any(p.returncode for p in children) or any(not p.is_file() for p in receipts):
            raise RuntimeError("Checked merge failed: missing successful reader receipt")
        import json
        records = [json.loads(path.read_text()) for path in receipts]
        for source, record in zip(sources, records):
            verify_digest(source, record)
        json_write(report, {"schema": 1, "reader_receipts": records,
                            "records": sum(r["records"] for r in records), "comparator": "pairtools-1.1.3/GNU/C"})
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
        for child in children:
            try:
                child.wait(timeout=3)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
        for fifo in fifos:
            with contextlib.suppress(FileNotFoundError):
                fifo.unlink()
        signal.signal(signal.SIGTERM, previous)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reader", type=Path)
    p.add_argument("--fifo", type=Path)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--tmpdir", type=Path)
    p.add_argument("--memory", default="8G")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--fan-in", type=int, default=8)
    p.add_argument("--readers", type=int, default=2)
    p.add_argument("--reader-lock-dir", type=Path)
    p.add_argument("--max-inputs", type=int, default=256)
    p.add_argument("sources", type=Path, nargs="*")
    args = p.parse_args()
    if args.reader:
        reader(args.reader, args.fifo, args.report, args.reader_lock_dir, args.readers)
    else:
        if not args.sources or args.fan_in < 2:
            p.error("Provide inputs and fan-in >= 2")
        merge(args.sources, args.tmpdir, args.report, args.memory, args.threads, args.fan_in, args.readers, args.max_inputs)


if __name__ == "__main__":
    main()
