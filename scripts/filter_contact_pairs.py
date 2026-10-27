#!/usr/bin/env python3
"""Shared accepted-contact selection for Micro-C, legacy Hi-C and restriction Hi-C."""
from __future__ import annotations

import argparse
import bisect
import json
import math
import sys
from collections import Counter
from pathlib import Path

from prepare_restriction_fragments import load_digest
from reference_registry import read_sizes
from workflow_common import json_write, text_open

FILTER_VERSION = 1
REQUIRED = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "pair_type", "mapq1", "mapq2"]


def classify(chrom1, pos1, strand1, frag1, chrom2, pos2, strand2, frag2):
    if chrom1 != chrom2:
        return "trans"
    if pos1 > pos2:
        strand1, strand2, frag1, frag2 = strand2, strand1, frag2, frag1
    if frag1 == frag2:
        if pos1 == pos2:
            return "same_fragment_other"
        return {( "+", "-"): "dangling_end_candidate", ("-", "+"): "self_circle_candidate"}.get(
            (strand1, strand2), "same_fragment_other")
    if abs(frag1 - frag2) == 1 and (strand1, strand2) == ("+", "-"):
        return "adjacent_inward_candidate"
    return "different_fragment_cis"


def filter_contacts(source: Path, output, policy: dict, sizes: dict[str, int], retained: set[str], digest=None) -> dict:
    columns = None
    counts, flags, distances = Counter(), Counter(), Counter()
    accepted_cis, accepted_trans = 0, 0
    restriction = bool(policy["enzymes"])
    if restriction and (digest is None or set(digest) != set(sizes) or any(digest[c][-1] != sizes[c] for c in sizes)):
        raise ValueError("Digest/reference chromosome or length inconsistency")
    with text_open(source) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].split()
                if len(set(columns)) != len(columns) or any(name not in columns for name in REQUIRED):
                    raise ValueError(f"Unsupported pairs schema in {source}: {columns}")
                if restriction and any(name in columns for name in ("rfrag1", "rfrag2")):
                    raise ValueError("Restriction annotation must be built from unannotated deduplicated pairs")
                output.write("#columns: " + " ".join(columns + (["rfrag1", "rfrag2"] if restriction else [])) + "\n")
                continue
            if line.startswith("#"):
                output.write(line)
                continue
            if not line.strip():
                continue
            if columns is None:
                raise ValueError("Pairs data precedes #columns header")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError(f"Malformed pairs row {number} in {source}")
            row = dict(zip(columns, fields))
            counts["input"] += 1
            if row["pair_type"] not in {"UU", "UR", "RU"}:
                counts["rejected_pair_type"] += 1
                continue
            if min(int(row["mapq1"]), int(row["mapq2"])) < policy["mapq"]:
                counts["rejected_mapq"] += 1
                continue
            c1, c2, p1, p2 = row["chrom1"], row["chrom2"], int(row["pos1"]), int(row["pos2"])
            if any(c not in sizes or not 1 <= p <= sizes[c] for c, p in ((c1, p1), (c2, p2))):
                raise ValueError(f"Eligible mapped end outside reference at pairs row {number}")
            if row["strand1"] not in {"+", "-"} or row["strand2"] not in {"+", "-"}:
                raise ValueError(f"Invalid mapped strand at pairs row {number}")
            if c1 not in retained or c2 not in retained:
                counts["rejected_chromosome"] += 1
                continue
            distance = abs(p2-p1)
            if c1 == c2:
                orientations = row["strand1"] + row["strand2"] if p1 <= p2 else row["strand2"] + row["strand1"]
                distances[f"{orientations}:log2_{int(math.log2(max(1, distance)))}"] += 1
            f1 = f2 = -1
            if restriction:
                f1, f2 = bisect.bisect_right(digest[c1], p1-1)-1, bisect.bisect_right(digest[c2], p2-1)-1
                category = classify(c1, p1, row["strand1"], f1, c2, p2, row["strand2"], f2)
                flags[category] += 1
                if category in {"dangling_end_candidate", "self_circle_candidate", "same_fragment_other"} and policy["remove_same_fragment"]:
                    counts["rejected_" + category] += 1
                    continue
                if category == "adjacent_inward_candidate" and policy["remove_adjacent_inward"]:
                    counts["rejected_adjacent_inward_candidate"] += 1
                    continue
            if c1 == c2 and policy["min_cis_dist"] > 0 and distance < policy["min_cis_dist"]:
                counts["rejected_cis_distance"] += 1
                continue
            counts["accepted"] += 1
            accepted_cis += c1 == c2
            accepted_trans += c1 != c2
            output.write("\t".join(fields + ([str(f1), str(f2)] if restriction else [])) + "\n")
    if columns is None:
        raise ValueError("Missing #columns pairs header")
    counts["input"] += 0
    counts["accepted"] += 0
    if sum(v for k, v in counts.items() if k != "input") != counts["input"]:
        raise AssertionError("Exclusive filter counts do not reconcile")
    return {"filter_version": FILTER_VERSION, "policy": policy, "exclusive_counts": dict(counts),
            "diagnostic_flags": dict(flags), "eligible_cis_distance_orientation_log2": dict(distances),
            "accepted_cis": accepted_cis, "accepted_trans": accepted_trans,
            "interpretation": "Orientation/fragment classes are candidates, not proof of molecular origin."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--matrix-chrom-sizes", type=Path, required=True)
    parser.add_argument("--digest", type=Path)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--legacy-pool", action="store_true", help="Chromosome selection only; upstream filtering is unverified")
    args = parser.parse_args()
    if args.legacy_pool:
        retained = set(read_sizes(args.matrix_chrom_sizes, "retained"))
        counts, cis, trans, columns = Counter(input=0, accepted=0), 0, 0, None
        with text_open(args.input) as handle:
            for line in handle:
                if line.startswith("#columns:"):
                    columns = line.split(":", 1)[1].split()
                if line.startswith("#"):
                    sys.stdout.write(line)
                elif line.strip():
                    if columns is None:
                        raise ValueError("Legacy pairs require a #columns header")
                    row = dict(zip(columns, line.rstrip().split("\t")))
                    counts["input"] += 1
                    if row["chrom1"] in retained and row["chrom2"] in retained:
                        counts["accepted"] += 1
                        cis += row["chrom1"] == row["chrom2"]
                        trans += row["chrom1"] != row["chrom2"]
                        sys.stdout.write(line)
                    else:
                        counts["rejected_chromosome"] += 1
        qc = {"exclusive_counts": dict(counts), "accepted_cis": cis, "accepted_trans": trans,
              "limited_provenance": "Legacy pool: upstream quality, distance, restriction filtering and duplicate scope are unverified"}
    else:
        qc = filter_contacts(args.input, sys.stdout, json.loads(args.policy.read_text()),
                             read_sizes(args.chrom_sizes, "reference"), set(read_sizes(args.matrix_chrom_sizes, "retained")),
                             load_digest(args.digest) if args.digest else None)
    json_write(args.qc, qc)


if __name__ == "__main__":
    main()
