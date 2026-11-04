"""Bounded vectorized equivalent of the RC5 scalar contact predicate."""
from collections import Counter
import math

from workflow_common import bounded_lines, text_open


def filter_batched(source, output, policy, sizes, retained, digest=None, batch_records=50000, batch_bytes=16777216):
    import numpy as np
    from filter_contact_pairs import FILTER_VERSION, REQUIRED
    restriction = bool(policy["enzymes"])
    if restriction and (digest is None or set(digest) != set(sizes) or any(digest[c][-1] != sizes[c] for c in sizes)):
        raise ValueError("Digest/reference chromosome or length inconsistency")
    if restriction:
        digest = {chrom: np.asarray(cuts, dtype=np.int64) for chrom, cuts in digest.items()}
    counts, flags, distances = Counter(input=0, accepted=0), Counter(), Counter()
    cis_total = trans_total = 0
    columns, pending, pending_bytes = None, [], 0

    def process(lines):
        nonlocal cis_total, trans_total
        if not lines:
            return
        fields = [line.rstrip("\r\n").split("\t") for line in lines]
        if any(len(row) != len(columns) for row in fields):
            raise ValueError("Malformed pairs row")
        data = np.asarray(fields, dtype=object)
        index = {name: i for i, name in enumerate(columns)}
        n = len(fields)
        counts["input"] += n
        alive = np.ones(n, dtype=bool)
        def reject(mask, name):
            selected = alive & mask
            count = int(selected.sum())
            if count:
                counts[name] += count
            alive[selected] = False
        reject(~np.isin(data[:, index["pair_type"]], ["UU", "UR", "RU"]), "rejected_pair_type")
        q1, q2 = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
        q1[alive] = np.asarray(data[alive, index["mapq1"]], dtype=np.int64)
        q2[alive] = np.asarray(data[alive, index["mapq2"]], dtype=np.int64)
        reject(np.minimum(q1, q2) < policy["mapq"], "rejected_mapq")
        c1, c2 = data[:, index["chrom1"]], data[:, index["chrom2"]]
        s1, s2 = data[:, index["strand1"]], data[:, index["strand2"]]
        p1, p2 = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
        p1[alive] = np.asarray(data[alive, index["pos1"]], dtype=np.int64)
        p2[alive] = np.asarray(data[alive, index["pos2"]], dtype=np.int64)
        lengths1 = np.fromiter((sizes.get(c, 0) for c in c1[alive]), dtype=np.int64)
        lengths2 = np.fromiter((sizes.get(c, 0) for c in c2[alive]), dtype=np.int64)
        if ((p1[alive] < 1) | (p2[alive] < 1) | (p1[alive] > lengths1) | (p2[alive] > lengths2)).any():
            raise ValueError("Eligible mapped end outside reference")
        if not (np.isin(s1[alive], ["+", "-"]).all() and np.isin(s2[alive], ["+", "-"]).all()):
            raise ValueError("Invalid mapped strand")
        reject(~np.isin(c1, list(retained)) | ~np.isin(c2, list(retained)), "rejected_chromosome")
        cis, distance = c1 == c2, np.abs(p2 - p1)
        reverse = p1 > p2
        left, right = np.where(reverse, s2, s1), np.where(reverse, s1, s2)
        eligible_cis = alive & cis
        unique_distances, inverse = np.unique(distance[eligible_cis], return_inverse=True)
        # Use the SAME scalar math.log2 rule once per unique value; frexp/float
        # approximations can change boundary QC for very large coordinates.
        bins = np.fromiter((int(math.log2(max(1, int(d)))) for d in unique_distances), dtype=np.int64)
        orientation = (left[eligible_cis] == "-").astype(np.int64) * 2 + (right[eligible_cis] == "-").astype(np.int64)
        histogram = np.bincount(bins[inverse] * 4 + orientation)
        names = ["++", "+-", "-+", "--"]
        for value in np.flatnonzero(histogram):
            distances[f"{names[value % 4]}:log2_{value // 4}"] += int(histogram[value])
        f1, f2 = np.full(n, -1, dtype=np.int64), np.full(n, -1, dtype=np.int64)
        if restriction:
            for chrom in set(c1[alive]) | set(c2[alive]):
                a, b = alive & (c1 == chrom), alive & (c2 == chrom)
                f1[a] = np.searchsorted(digest[chrom], p1[a] - 1, side="right") - 1
                f2[b] = np.searchsorted(digest[chrom], p2[b] - 1, side="right") - 1
            category = np.full(n, "different_fragment_cis", dtype=object)
            category[~cis] = "trans"
            same = cis & (f1 == f2)
            category[same] = "same_fragment_other"
            category[same & (p1 != p2) & (left == "+") & (right == "-")] = "dangling_end_candidate"
            category[same & (p1 != p2) & (left == "-") & (right == "+")] = "self_circle_candidate"
            category[cis & ~same & (np.abs(f1 - f2) == 1) & (left == "+") & (right == "-")] = "adjacent_inward_candidate"
            categories, multiplicities = np.unique(category[alive], return_counts=True)
            flags.update({str(k): int(v) for k, v in zip(categories, multiplicities)})
            if policy["remove_same_fragment"]:
                for label in ("dangling_end_candidate", "self_circle_candidate", "same_fragment_other"):
                    reject(category == label, "rejected_" + label)
            if policy["remove_adjacent_inward"]:
                reject(category == "adjacent_inward_candidate", "rejected_adjacent_inward_candidate")
        if policy["min_cis_dist"] > 0:
            reject(cis & (distance < policy["min_cis_dist"]), "rejected_cis_distance")
        indices = np.flatnonzero(alive)
        counts["accepted"] += len(indices)
        cis_total += int((alive & cis).sum())
        trans_total += int((alive & ~cis).sum())
        if restriction:
            fragments1, fragments2 = f1[indices].astype(str), f2[indices].astype(str)
            output.writelines("\t".join(fields[i] + [a, b]) + "\n" for i, a, b in zip(indices, fragments1, fragments2))
        else:
            output.writelines("\t".join(fields[i]) + "\n" for i in indices)

    with text_open(source) as handle:
        for line in bounded_lines(handle, batch_bytes):
            if line.startswith("#"):
                process(pending)
                pending, pending_bytes = [], 0
                if line.startswith("#columns:"):
                    columns = line.split(":", 1)[1].split()
                    if len(set(columns)) != len(columns) or any(c not in columns for c in REQUIRED):
                        raise ValueError(f"Unsupported pairs schema: {columns}")
                    if restriction and any(c in columns for c in ("rfrag1", "rfrag2")):
                        raise ValueError("Restriction annotation must be built from unannotated deduplicated pairs")
                    output.write("#columns: " + " ".join(columns + (["rfrag1", "rfrag2"] if restriction else [])) + "\n")
                else:
                    output.write(line)
            elif line.strip():
                if columns is None:
                    raise ValueError("Pairs data precedes #columns header")
                amount = len(line.encode("utf-8"))
                if amount > batch_bytes:
                    raise ValueError("Single pairs record exceeds FILTER_BATCH_BYTES")
                if pending and (len(pending) >= batch_records or pending_bytes + amount > batch_bytes):
                    process(pending)
                    pending, pending_bytes = [], 0
                pending.append(line)
                pending_bytes += amount
    process(pending)
    if columns is None:
        raise ValueError("Missing #columns pairs header")
    if sum(v for k, v in counts.items() if k != "input") != counts["input"]:
        raise AssertionError("Exclusive filter counts do not reconcile")
    return {"filter_version": FILTER_VERSION, "policy": policy, "exclusive_counts": dict(counts),
            "diagnostic_flags": dict(flags), "eligible_cis_distance_orientation_log2": dict(distances),
            "accepted_cis": cis_total, "accepted_trans": trans_total,
            "interpretation": "Orientation/fragment classes are candidates, not proof of molecular origin."}
