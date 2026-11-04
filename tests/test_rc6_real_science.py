"""Small real pairtools/Cooler equivalence checks; no FASTQs, JAR or real genomes.

MICROC2TRACKS_PAIR_TOOLS_TESTS=1 python -m unittest discover -s tests -p test_rc6_real_science.py -v
"""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from pool_cooler import merge_counts


@unittest.skipUnless(os.environ.get("MICROC2TRACKS_PAIR_TOOLS_TESTS") == "1", "opt-in pinned pairtools/Cooler checks")
class RealScience(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    @staticmethod
    def data(path):
        return [line for line in path.read_text().splitlines() if not line.startswith("#")]

    def test_streamed_joint_dedup_matches_materialized_including_chunk_bridge(self):
        import pairtools
        self.assertEqual(pairtools.__version__, "1.1.3")
        header = ("## pairs format v1.0\n#genome_assembly: toy\n#chromsize: scaffoldA 200000\n"
                  "#chromsize: plasmidZ 200000\n#shape: upper triangle\n#sorted: chr1-chr2-pos1-pos2\n"
                  "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n")
        rows = [[f"f{i}\tscaffoldA\t{i*10+1}\tplasmidZ\t{i*10+5}\t+\t-\tUU\t60\t60\n" for i in range(9998)], []]
        rows[0] += ["bridgeA\tscaffoldA\t150000\tplasmidZ\t150000\t+\t-\tUU\t60\t60\n", "bridgeC\tscaffoldA\t150006\tplasmidZ\t150006\t+\t-\tUU\t60\t60\n"]
        rows[1] += ["bridgeB\tscaffoldA\t150003\tplasmidZ\t150003\t+\t-\tUU\t60\t60\n"]
        sources, compressed = [], []
        for i, records in enumerate(rows):
            path = self.root / f"input{i}.pairs"
            path.write_text(header + "".join(records))
            sources.append(path)
            gz = self.root / f"input{i}.pairs.gz"
            gz.write_bytes(gzip.compress(path.read_bytes()))
            compressed.append(gz)
        merged = self.root / "materialized.pairs"
        with merged.open("w") as out:
            subprocess.run(["pairtools", "merge", *map(str, sources)], stdout=out, check=True)
        options = ["pairtools", "dedup", "--backend", "scipy", "--method", "max", "--max-mismatch", "3", "--chunksize", "10000", "--carryover", "100"]
        baseline, streamed = self.root / "baseline", self.root / "streamed"
        with baseline.open("w") as out:
            subprocess.run([*options, str(merged)], stdout=out, check=True)
        checked = subprocess.Popen([sys.executable, str(ROOT / "scripts/checked_merge.py"), "--tmpdir", str(self.root / "work"),
            "--report", str(self.root / "receipt.json"), *map(str, compressed)], stdout=subprocess.PIPE)
        try:
            with streamed.open("w") as out:
                subprocess.run(options, stdin=checked.stdout, stdout=out, check=True)
            checked.stdout.close()
            self.assertEqual(checked.wait(timeout=20), 0)
        finally:
            if checked.poll() is None:
                checked.kill()
                checked.wait()
        self.assertEqual(self.data(baseline), self.data(streamed))
        self.assertEqual(sum(line.startswith("bridge") for line in self.data(streamed)), 2)
        self.assertEqual(json.loads((self.root / "receipt.json").read_text())["records"], 10001)

    def test_raw_count_pooling_and_normalization_match_direct_counts(self):
        import cooler
        import numpy as np
        import pandas as pd
        bins = cooler.binnify(pd.Series({"chr1": 10000, "chr2": 10000}), 1000)
        a, b = np.triu_indices(len(bins))
        left = pd.DataFrame({"bin1_id": a, "bin2_id": b, "count": (a + b) % 7 + 1})
        right = pd.DataFrame({"bin1_id": a, "bin2_id": b, "count": (a * 3 + b) % 5 + 1})
        files = [self.root / "a.cool", self.root / "b.cool"]
        for path, pixels in zip(files, [left, right]):
            weighted = bins.copy()
            weighted["weight"] = np.arange(len(bins)) + 1.0
            cooler.create_cooler(str(path), weighted, pixels, assembly="toy", ordered=True)
        expected = self.root / "direct.cool"
        direct = left.copy()
        direct["count"] += right["count"]
        cooler.create_cooler(str(expected), bins, direct, assembly="toy", ordered=True)
        pooled = self.root / "pooled.cool"
        merge_counts(files, pooled, buffer=13, expected_total=int(direct["count"].sum()))
        x, y = cooler.Cooler(str(expected)), cooler.Cooler(str(pooled))
        np.testing.assert_array_equal(x.pixels()[:].to_numpy(), y.pixels()[:].to_numpy())
        self.assertNotIn("weight", y.bins().columns)
        for value in (x, y):
            weights, stats = cooler.balance_cooler(value, min_nnz=1, min_count=0, ignore_diags=0, tol=1e-6, max_iters=1000)
            if value is x:
                first_weights, first_stats = weights, stats
            else:
                np.testing.assert_array_equal(np.isnan(first_weights), np.isnan(weights))
                np.testing.assert_allclose(first_weights, weights, rtol=1e-6, atol=1e-10, equal_nan=True)
                np.testing.assert_array_equal(first_stats["converged"], stats["converged"])
        for label, path in (("direct", expected), ("pooled", pooled)):
            cooler.zoomify_cooler(str(path), str(self.root / f"{label}.mcool"), [1000, 2000, 5000], chunksize=23,
                                  nproc=1, columns=["count"], dtypes={"count": np.int64})
        for resolution in (1000, 2000, 5000):
            a = cooler.Cooler(f"{self.root / 'direct.mcool'}::/resolutions/{resolution}")
            b = cooler.Cooler(f"{self.root / 'pooled.mcool'}::/resolutions/{resolution}")
            np.testing.assert_array_equal(a.pixels()[:].to_numpy(), b.pixels()[:].to_numpy())
            wa, sa = cooler.balance_cooler(a, min_nnz=1, min_count=0, ignore_diags=0, tol=1e-6, max_iters=1000)
            wb, sb = cooler.balance_cooler(b, min_nnz=1, min_count=0, ignore_diags=0, tol=1e-6, max_iters=1000)
            np.testing.assert_array_equal(np.isnan(wa), np.isnan(wb))
            np.testing.assert_array_equal(sa["converged"], sb["converged"])
            np.testing.assert_allclose(wa, wb, rtol=1e-6, atol=1e-10, equal_nan=True)
        with self.assertRaisesRegex(ValueError, "QC"):
            merge_counts(files, self.root / "wrong.cool", expected_total=0)
        huge = pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [2**62]})
        big_files = [self.root / "big_a.cool", self.root / "big_b.cool"]
        for path in big_files:
            cooler.create_cooler(str(path), bins, huge, assembly="toy", ordered=True, dtypes={"count": np.int64})
        with self.assertRaises(OverflowError):
            merge_counts(big_files, self.root / "overflow.cool")

    def test_merged_accepted_stats_match_direct_pool(self):
        header = "## pairs format v1.0\n#genome_assembly: toy\n#chromsize: chr1 100000\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
        inputs, stats, records = [], [], []
        for i in range(2):
            data = [f"r{i}_{j}\tchr1\t{j+1}\tchr1\t{100+j*500}\t+\t-\tUU\n" for j in range(40)]
            source, stat = self.root / f"pairs{i}", self.root / f"stats{i}"
            source.write_text(header + "".join(data))
            subprocess.run(["pairtools", "stats", "-o", str(stat), str(source)], check=True)
            records += data
            inputs.append(source)
            stats.append(stat)
        pooled, direct, merged = self.root / "pool", self.root / "direct.stats", self.root / "merged.stats"
        pooled.write_text(header + "".join(records))
        subprocess.run(["pairtools", "stats", "-o", str(direct), str(pooled)], check=True)
        subprocess.run(["pairtools", "stats", "--merge", "-o", str(merged), *map(str, stats)], check=True)
        def values(path):
            return dict(line.split("\t", 1) for line in path.read_text().splitlines())
        self.assertEqual(values(direct), values(merged))


if __name__ == "__main__":
    unittest.main()
