#!/usr/bin/env python3
"""Pool raw integer counts only; caller certifies disjoint accepted membership."""
import argparse
from pathlib import Path


def merge_counts(sources, destination, buffer=1000000, expected_total=None):
    import cooler
    import numpy as np
    from cooler._reduce import CoolerMerger
    coolers = [cooler.Cooler(str(p)) for p in sources]
    if not coolers:
        raise ValueError("No count matrices")
    first, total = coolers[0], 0
    for item in coolers:
        if item.storage_mode != "symmetric-upper" or item.binsize != first.binsize or item.info.get("genome-assembly") != first.info.get("genome-assembly"):
            raise ValueError("Incompatible Cooler reference/resolution/storage mode")
        if item.info["nbins"] != first.info["nbins"]:
            raise ValueError("Incompatible Cooler bin axes")
        if not np.issubdtype(item.pixels().dtypes["count"], np.integer):
            raise ValueError("Pooling requires integer raw counts, not normalized values")
        for lo in range(0, first.info["nbins"], buffer):
            a = first.bins()[["chrom", "start", "end"]][lo:lo + buffer].reset_index(drop=True)
            b = item.bins()[["chrom", "start", "end"]][lo:lo + buffer].reset_index(drop=True)
            if not a.equals(b):
                raise ValueError("Incompatible Cooler bin axes")
        value = int(item.info["sum"])
        if value < 0:
            raise ValueError("Negative raw-count total")
        total += value
    if total > np.iinfo(np.int64).max:
        raise OverflowError("Raw-count total exceeds signed int64; pooling refused")
    if expected_total is not None and total != expected_total:
        raise ValueError("Constituent Cooler totals disagree with accepted-pair QC")
    seen = 0
    def chunks():
        nonlocal seen
        for chunk in CoolerMerger(coolers, mergebuf=buffer, columns=["count"]):
            values = chunk["count"]
            if (values < 0).any():
                raise OverflowError("Negative/overflowed pooled raw counts")
            seen += sum(map(int, values))
            yield chunk
    # Weight columns are deliberately excluded from both bins and pixels.
    cooler.create_cooler(str(destination), first.bins()[["chrom", "start", "end"]][:], chunks(),
        dtypes={"count": np.int64}, ordered=True, symmetric_upper=True,
        assembly=first.info.get("genome-assembly"))
    if seen != total or int(cooler.Cooler(str(destination)).info["sum"]) != total:
        raise ValueError("Pooled raw counts do not reconcile")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--buffer", type=int, default=1000000)
    p.add_argument("--expected-total", type=int)
    p.add_argument("sources", nargs="+")
    args = p.parse_args()
    if args.buffer < 1:
        p.error("Buffer must be positive")
    merge_counts(args.sources, args.output, args.buffer, args.expected_total)


if __name__ == "__main__":
    main()
