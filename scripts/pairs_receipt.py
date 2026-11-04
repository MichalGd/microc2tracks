#!/usr/bin/env python3
"""Bounded pass-through count/digest receipt collected during an existing write."""
import argparse
import hashlib
from pathlib import Path
import sys
from workflow_common import json_write


def copy_receipt(inp, out):
    digest, count, last = hashlib.sha256(), 0, b""
    line = inp.readline()
    while line.startswith(b"#"):
        out.write(line)
        line = inp.readline()
    while line:
        out.write(line)
        digest.update(line)
        count += line.count(b"\n")
        last = line[-1:]
        line = inp.read(1024 * 1024)
    if last and last != b"\n":
        raise ValueError("Incomplete final pairs record")
    out.flush()
    return {"records": count, "data_sha256": digest.hexdigest()}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report", type=Path, required=True)
    args = p.parse_args()
    json_write(args.report, copy_receipt(sys.stdin.buffer, sys.stdout.buffer))
