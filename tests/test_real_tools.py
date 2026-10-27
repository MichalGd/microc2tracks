"""Opt-in tiny real-tool checks. No production FASTQs or existing outputs are used.

MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR=/path/to/juicer_tools.2.18.00.jar \
  python -m unittest discover -s tests -p test_real_tools.py -v
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from reference_registry import REGISTRY_COLUMNS
from workflow_common import DEFAULTS, text_open, sha256


@unittest.skipUnless(os.environ.get("MICROC2TRACKS_REAL_TOOLS") == "1", "opt-in Linux scientific tools not requested")
class RealTools(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        missing = [p for p in ("bwa-mem2", "pairtools", "pairix", "bgzip", "cooler", "java") if not shutil.which(p)]
        cls.jar = Path(os.environ.get("JUICER_TOOLS_JAR", ""))
        if missing or not cls.jar.is_file():
            raise RuntimeError(f"Real-tool test requested but tools/JAR missing: {missing}; jar={cls.jar}")
        cls.temp = tempfile.TemporaryDirectory(prefix="microc2tracks_real_test_")
        cls.root = Path(cls.temp.name)
        rng = random.Random(1048)
        cls.sequences = {name: "".join(rng.choices("ACGT", k=200000)) for name in ("scaffoldA", "plasmidZ")}
        fasta = cls.root / "reference.fa"
        with fasta.open("w") as out:
            for chrom, seq in cls.sequences.items():
                out.write(f">{chrom}\n")
                for i in range(0, len(seq), 60): out.write(seq[i:i+60] + "\n")
        sizes = cls.root / "sizes.tsv"
        sizes.write_text("scaffoldA\t200000\nplasmidZ\t200000\n")
        # The registry validates names/lengths; BWA builds its actual index below.
        Path(str(fasta) + ".fai").write_text(sizes.read_text())
        subprocess.run(["bwa-mem2", "index", str(fasta)], check=True, stdout=subprocess.DEVNULL)
        cls.registry = cls.root / "references.tsv"
        with cls.registry.open("w", newline="") as out:
            writer = csv.DictWriter(out, REGISTRY_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"reference_id": "toy", "species": "Synthetic", "assembly": "toy", "fasta": str(fasta),
                "bwa_index_prefix": str(fasta), "chrom_sizes": str(sizes), "canonical_regex": ".*", "browser_preset": "custom"})

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_fixture(self, assay, profile):
        import cooler
        work = self.root / (assay + "_" + profile)
        work.mkdir()
        sheet = work / "samples.csv"
        rows = []
        for lane in (1, 2):
            r1, r2 = work / f"part{lane}.R1.fastq", work / f"part{lane}.R2.fastq"
            with r1.open("w") as one, r2.open("w") as two:
                for i in range(48):
                    c1, c2 = (("scaffoldA", "scaffoldA"), ("plasmidZ", "plasmidZ"), ("scaffoldA", "plasmidZ"))[i%3]
                    offset = 0 if lane == 1 or i<12 else 700
                    p1, p2 = 2000+i*1700+offset, 12001+i*1700+offset
                    seq1 = self.sequences[c1][p1:p1+100]
                    seq2 = self.sequences[c2][p2:p2+100].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                    one.write(f"@part{lane}_{i}/1\n{seq1}\n+\n{'I'*100}\n")
                    two.write(f"@part{lane}_{i}/2\n{seq2}\n+\n{'I'*100}\n")
            rows.append({"sample": f"part{lane}", "assay": assay, "reference_genome": "toy", "condition": "synthetic",
                "biological_replicate": "1", "technical_replicate": str(lane), "fastq_r1": str(r1), "fastq_r2": str(r2), "dedup_group": "joint"})
        with sheet.open("w", newline="") as out:
            writer = csv.DictWriter(out, list(rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        config = work / "config.conf"
        values = {**DEFAULTS, "OUTDIR": str(work / "out"), "REFERENCE_REGISTRY": str(self.registry), "DEFAULT_REFERENCE_ID": "toy",
            "TMPDIR": str(work / "ssd_scratch"), "PAIRS_OUTDIR": str(work / "hdd_pairs"), "MATRIX_OUTDIR": str(work / "hdd_matrices"),
            "TRIMMED_OUTDIR": str(work / "hdd_trimmed"), "RUN_FASTP": "false", "RUN_MULTIQC": "false", "RUN_GLOBAL_MULTIQC": "false",
            "RUN_PRELIM_DOWNSTREAM": "false", "RUN_MERGE_ALL_SAMPLES": "true", "HIC_PROTOCOL": profile,
            "HIC_RESTRICTION_ENZYMES": "DpnII" if profile == "restriction_hic" else "", "MIN_FREE_OUTPUT_GB": "0", "MIN_FREE_SCRATCH_GB": "0",
            "RESOLUTIONS": "1000", "MCOOL_RESOLUTIONS": "1000", "HIC_RESOLUTIONS": "1000", "JAVA_HEAP": "1g", "HIC_NORMALIZATIONS": "VC",
            "JUICER_TOOLS_JAR": str(self.jar.resolve()), "THREADS_ALIGN": "1", "THREADS_MATRIX": "1", "THREADS_SORT": "1", "THREADS_MERGE": "1", "THREADS_HIC_NORM": "1"}
        import shlex
        config.write_text("".join(f"{key}={shlex.quote(str(value))}\n" for key, value in values.items()))
        run_log = work / "run.log"
        with run_log.open("w") as out:
            result = subprocess.run([sys.executable, str(ROOT / "scripts/pipeline.py"), "run", "-c", str(config), "-s", str(sheet)], stdout=out, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, run_log.read_text())
        plan = json.loads((work / "out/run_metadata/run_plan.json").read_text())
        pool = plan["entities"]["pool:all_samples_merged"]
        paths = pool["paths"]
        qc = json.loads(Path(paths["filter_qc"]).read_text())
        expected_count = qc["exclusive_counts"]["accepted"]
        clr = cooler.Cooler(paths["norm_mcool"] + "::resolutions/1000")
        self.assertGreater(expected_count, 0)
        self.assertEqual(clr.info["sum"], expected_count)
        standalone = sum(json.loads(Path(plan["entities"][f"sample:part{lane}"]["paths"]["filter_qc"]).read_text())["exclusive_counts"]["accepted"] for lane in (1, 2))
        self.assertLess(expected_count, standalone)
        self.assertEqual(pool["units"], ["dedup:joint"])
        bins = clr.bins()[:]
        cooler_counts = {}
        for record in clr.pixels()[:].itertuples():
            a, b = bins.iloc[record.bin1_id], bins.iloc[record.bin2_id]
            cooler_counts[(str(a.chrom), int(a.start), str(b.chrom), int(b.start))] = int(record.count)
        hic_counts = {}
        for i, c1 in enumerate(self.sequences):
            for c2 in list(self.sequences)[i:]:
                dump = work / f"{c1}_{c2}.dump.tsv"
                subprocess.run(["java", "-Xmx1g", "-jar", str(self.jar), "dump", "observed", "NONE", paths["norm_hic"], c1, c2, "BP", "1000", str(dump)], check=True)
                for line in dump.read_text().splitlines():
                    x, y, count = line.split()
                    if float(count):
                        x, y = int(x), int(y)
                        if c1 == c2 and x > y: x, y = y, x
                        hic_counts[(c1, x, c2, y)] = float(count)
        self.assertEqual(cooler_counts, hic_counts)
        before = Path(paths["norm_mcool"]).stat().st_mtime_ns
        with run_log.open("a") as out:
            result = subprocess.run([sys.executable, str(ROOT / "scripts/pipeline.py"), "run", "-c", str(config), "-s", str(sheet)], stdout=out, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, run_log.read_text())
        self.assertEqual(Path(paths["norm_mcool"]).stat().st_mtime_ns, before)

        if assay == "hic" and profile == "restriction_hic":
            # Refilter the already jointly deduplicated stream, before cis/fragment
            # selection, without repeating alignment or deduplication.
            source = Path(plan["entities"]["dedup:joint"]["paths"]["dedup"])
            original_hash = sha256(source)
            refilter_sheet = work / "refilter.csv"
            with refilter_sheet.open("w", newline="") as out:
                writer = csv.DictWriter(out, ["sample", "assay", "reference_genome", "deduplicated_pairs", "input_members"])
                writer.writeheader()
                writer.writerow({"sample": "existing_joint", "assay": "hic", "reference_genome": "toy",
                                 "deduplicated_pairs": str(source), "input_members": "part1;part2"})
            imported_values = {**values, "WORKFLOW_MODE": "refilter-deduplicated", "OUTDIR": str(work / "refilter_out"),
                "PAIRS_OUTDIR": str(work / "refilter_pairs"), "MATRIX_OUTDIR": str(work / "refilter_matrices"),
                "TRIMMED_OUTDIR": "", "RUN_MERGE_ALL_SAMPLES": "false", "RUN_FASTP": "true", "HIC_MIN_CIS_DIST": "5000"}
            refilter_config = work / "refilter.conf"
            refilter_config.write_text("".join(f"{key}={shlex.quote(str(value))}\n" for key,value in imported_values.items()))
            with run_log.open("a") as out:
                result = subprocess.run([sys.executable, str(ROOT / "scripts/pipeline.py"), "run", "-c", str(refilter_config), "-s", str(refilter_sheet)], stdout=out, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, run_log.read_text())
            imported_plan = json.loads((work / "refilter_out/run_metadata/run_plan.json").read_text())
            imported_paths = imported_plan["entities"]["sample:existing_joint"]["paths"]
            imported_clr = cooler.Cooler(imported_paths["norm_mcool"] + "::resolutions/1000")
            self.assertEqual(imported_clr.info["sum"], expected_count)
            for skipped in ("trim", "sorted", "dedup"):
                self.assertFalse((Path(imported_paths["done"]) / f"{skipped}.json").exists())
            imported_hic_counts = {}
            for i, c1 in enumerate(self.sequences):
                for c2 in list(self.sequences)[i:]:
                    dump = work / f"refilter_{c1}_{c2}.dump.tsv"
                    subprocess.run(["java", "-Xmx1g", "-jar", str(self.jar), "dump", "observed", "NONE", imported_paths["norm_hic"], c1, c2, "BP", "1000", str(dump)], check=True)
                    for line in dump.read_text().splitlines():
                        x,y,count = line.split()
                        if float(count):
                            x,y = int(x),int(y)
                            if c1==c2 and x>y: x,y=y,x
                            imported_hic_counts[(c1,x,c2,y)] = float(count)
            self.assertEqual(imported_hic_counts, cooler_counts)
            self.assertEqual(sha256(source), original_hash)

    def test_tiny_microc_full_workflow(self):
        self.run_fixture("microc", "restriction_hic")

    def test_tiny_generic_hic_full_workflow(self):
        self.run_fixture("hic", "restriction_hic")

    def test_tiny_dovetail_hic_full_workflow(self):
        self.run_fixture("hic", "dovetail_hic")

    def test_single_joint_pass_bridge_chunk_boundary_and_independent_libraries(self):
        work = self.root / "dedup_bridge"
        work.mkdir()
        # Three near-matches form a bridge across the production chunk boundary.
        # The expected answer is always the pinned single joint pass, not a second
        # pass over representatives discarded by previous independent deduplication.
        header = "## pairs format v1.0\n#genome_assembly: toy\n#chromsize: scaffoldA 200000\n#chromsize: plasmidZ 200000\n#shape: upper triangle\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n"
        paths = [work / "a.pairs", work / "b.pairs"]
        records = [[f"f{i}\tscaffoldA\t{i*10+1}\tplasmidZ\t{i*10+5}\t+\t-\tUU\t60\t60\n" for i in range(9998)], []]
        records[0] += ["bridgeA\tscaffoldA\t150000\tplasmidZ\t150000\t+\t-\tUU\t60\t60\n", "bridgeC\tscaffoldA\t150006\tplasmidZ\t150006\t+\t-\tUU\t60\t60\n"]
        records[1] += ["bridgeB\tscaffoldA\t150003\tplasmidZ\t150003\t+\t-\tUU\t60\t60\n"]
        for path, data in zip(paths, records): path.write_text(header + "".join(data))
        combined = work / "combined.pairs"
        with combined.open("w") as out:
            subprocess.run(["pairtools", "merge", *map(str, paths)], check=True, stdout=out)
        def dedup(source, dest, chunksize="10000"):
            with dest.open("w") as out:
                subprocess.run(["pairtools", "dedup", "--backend", "scipy", "--method", "max", "--max-mismatch", "3", "--chunksize", chunksize, "--carryover", "100", str(source)], check=True, stdout=out)
            return [line.split("\t")[0] for line in dest.read_text().splitlines() if not line.startswith("#")]
        joint = dedup(combined, work / "joint.pairs")
        # pairtools 1.1.3 carries representatives, not discarded bridge records.
        # Pin its chunk-sensitive behavior explicitly; do not claim global closure.
        self.assertEqual(len([name for name in joint if name.startswith("bridge")]), 2)
        one_chunk = dedup(combined, work / "one_chunk.pairs", chunksize="20000")
        self.assertEqual(len([name for name in one_chunk if name.startswith("bridge")]), 1)
        independent = dedup(paths[0], work / "independent_a.pairs") + dedup(paths[1], work / "independent_b.pairs")
        self.assertEqual(len([name for name in independent if name.startswith("bridge")]), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
