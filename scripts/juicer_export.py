#!/usr/bin/env python3
"""Bounded short-format projection; compression is a separately checked process."""
import argparse
from pathlib import Path
import sys

from workflow_common import bounded_lines, json_write, text_open


def project(source, output, batch_bytes=1048576, check=lambda: None):
    columns, selected, buffer, count, amount = None, None, [], 0, 0
    with text_open(source) as handle:
        for line in bounded_lines(handle, max(batch_bytes, 16777216)):
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].split()
                names = ["chrom1", "chrom2", "pos1", "pos2", "strand1", "strand2"]
                if len(set(columns)) != len(columns) or any(n not in columns for n in names):
                    raise ValueError("Invalid Juicer export schema")
                selected = [columns.index(n) for n in names]
            elif not line.startswith("#") and line.strip():
                if columns is None:
                    raise ValueError("Juicer export requires a #columns schema")
                fields = line.rstrip().split("\t")
                if len(fields) != len(columns):
                    raise ValueError("Malformed accepted pairs record")
                c1, c2, p1, p2, s1, s2 = (fields[i] for i in selected)
                p1, p2 = int(p1) - 1, int(p2) - 1
                if min(p1, p2) < 0:
                    raise ValueError("Accepted pairs must have positive one-based positions")
                record = f"{'0' if s1 == '+' else '16'}\t{c1}\t{p1}\t0\t{'0' if s2 == '+' else '16'}\t{c2}\t{p2}\t1\n"
                buffer.append(record)
                count += 1
                amount += len(record)
                if amount >= batch_bytes:
                    output.writelines(buffer)
                    buffer, amount = [], 0
                    check()
    output.writelines(buffer)
    check()
    if columns is None:
        raise ValueError("Missing #columns in Juicer export")
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--batch-bytes", type=int, default=1048576)
    args = p.parse_args()
    if args.batch_bytes < 1:
        p.error("Batch size must be positive")
    count = project(args.input, sys.stdout, args.batch_bytes)
    sys.stdout.flush()
    json_write(args.report, {"records": count, "coordinates": "pos-1", "format": "juicer-short"})


if __name__ == "__main__":
    main()
