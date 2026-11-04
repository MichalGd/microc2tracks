#!/usr/bin/env python3
"""Bounded synthetic filter comparison, not a production workflow speed claim."""
import argparse
from array import array
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from workflow_common import json_write


class DigestWriter:
    def __init__(self):
        self.digest = hashlib.sha256()
    def write(self, value):
        self.digest.update(value.encode())
    def writelines(self, values):
        for value in values:
            self.write(value)


def worker(engine, source):
    import numpy  # exclude one-time import cost equally for this kernel comparison
    from filter_contact_pairs import filter_contacts
    from batched_filter import filter_batched
    policy = {"enzymes": ["DpnII"], "mapq": 30, "min_cis_dist": 0,
              "remove_same_fragment": True, "remove_adjacent_inward": False}
    sizes = {"chr1": 100000000, "chr2": 100000000}
    digest = {c: array("q", [*range(0, length, 256), length]) for c, length in sizes.items()}
    output = DigestWriter()
    start = time.perf_counter()
    qc = (filter_contacts if engine == "scalar" else filter_batched)(source, output, policy, sizes, set(sizes), digest)
    seconds = time.perf_counter() - start
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    print(json.dumps({"engine": engine, "seconds": seconds, "peak_process_rss_bytes": peak,
                      "output_sha256": output.digest.hexdigest(), "qc": qc}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--records", type=int, default=200000)
    p.add_argument("--report", type=Path)
    p.add_argument("--engine", choices=["scalar", "batched"])
    p.add_argument("--input", type=Path)
    args = p.parse_args()
    if args.engine:
        worker(args.engine, args.input)
        return
    if not 1 <= args.records <= 2000000:
        p.error("Use 1..2000000 records; no production files are accepted")
    with tempfile.TemporaryDirectory(prefix="microc2tracks_rc6_benchmark_") as directory:
        source = Path(directory) / "synthetic.pairs"
        with source.open("w") as out:
            out.write("## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n")
            for i in range(args.records):
                p1 = (i * 7919) % 90000000 + 1
                p2 = (p1 + 500 + (i * 577) % 9000000) % 100000000 + 1
                out.write(f"r{i}\tchr1\t{p1}\t{'chr2' if i % 5 == 0 else 'chr1'}\t{p2}\t+\t-\tUU\t60\t60\n")
        results = []
        for engine in ("scalar", "batched"):
            output = subprocess.check_output([sys.executable, "-B", __file__, "--engine", engine, "--input", str(source)], text=True)
            results.append(json.loads(output))
        if results[0]["output_sha256"] != results[1]["output_sha256"] or results[0]["qc"] != results[1]["qc"]:
            raise ValueError("Synthetic benchmark detected non-equivalent filtering")
        report = {"scope": "Synthetic filtering only; not production throughput or whole-workflow acceleration",
                  "input_records": args.records, "input_bytes": source.stat().st_size, "results": results}
        print(json.dumps(report, indent=2))
        if args.report:
            json_write(args.report, report)


if __name__ == "__main__":
    main()
