"""Small real-linker and integration regressions; no sequencing downloads."""
import copy
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import test_integration_features as fixtures
from artifact_state import StorageGuard
from blhic_preprocessing import PRESET, identify_tool, parse_stats, scientific_tool_identity
from check_fastq_pair_names import inspect_pairs
from pairs_processing import PairProcessor
from pipeline import Runtime, tool_versions
from run_plan import build_plan
from stage_scheduler import reservation
from summarize_run import write_plan_summary
from workflow_common import identity, json_write, resolve_policy, resolve_preprocessing


def write_fastq(path, records, mate):
    path.write_bytes(b"".join(f"@{name}/{mate}\n{sequence}\n+\n{'I'*len(sequence)}\n".encode() for name, sequence in records))


def read_fastq(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        lines = handle.read().splitlines()
    return {lines[i].split('/')[0][1:]: (lines[i+1], lines[i+3]) for i in range(0, len(lines), 4)}


@unittest.skipUnless(sys.platform.startswith("linux"), "The optional BL executable is tested on Linux")
class BLHiCTests(unittest.TestCase):
    def setUp(self):
        if not (ROOT / "tools/blhic/trimLinker").is_file() and not shutil.which("trimLinker"):
            self.skipTest("Build the optional tool with python3 scripts/build_trimlinker.py before BL tests")
        self.fixture = fixtures.IntegrationFeatures("test_profiles_generic_dovetail_legacy_and_explicit_enzyme_requirement")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.addCleanup(self.fixture.tearDown)
        self.root = self.fixture.root
        self.c = {**self.fixture.c, "HIC_PROTOCOL": "blhic", "HIC_RESTRICTION_ENZYMES": "", "RUN_FASTP": "false",
                  "RUN_HIC": "false", "RUN_MULTIQC": "false", "RUN_GLOBAL_MULTIQC": "false",
                  "RUN_MERGE_TECHNICAL_REPLICATES": "false"}
        self.r1, self.r2 = self.root / "raw_R1.fastq", self.root / "raw_R2.fastq"
        a, b, g = PRESET["linker_a"], PRESET["linker_b"], "C"*20
        self.input1 = [("ab", g+a+"TT"), ("aa", g+a), ("one", g+a), ("empty", g), ("short", "C"*14+a), ("edge", "C"*15+a)]
        self.input2 = [("ab", g+b+"TT"), ("aa", g+a), ("one", g), ("empty", g), ("short", g+b), ("edge", g+b)]
        write_fastq(self.r1, self.input1, 1)
        write_fastq(self.r2, self.input2, 2)
        self.fixture.rows = [{**self.fixture.rows[0], "assay": "hic", "fastq_r1": str(self.r1), "fastq_r2": str(self.r2),
                              "dedup_group": "", "hic_protocol": "", "restriction_enzymes": ""}]
        self.fixture.write_sheet()
        self.calls = []

    def runtime(self, c=None):
        plan = build_plan(c or self.c, self.fixture.sheet)
        versions = {"pairtools": "1.1.3", "bwa-mem2": "fixture", "fastp": "fixture"}
        if any(r.get("preprocessing", {}).get("input_state") == "raw" for r in plan["samples"].values()):
            versions["trimLinker"] = identify_tool(c or self.c)
        return Runtime(plan, self.root / "config", versions=versions)

    def test_profile_defaults_precedence_and_conflicts(self):
        policy = resolve_policy({"assay": "hic"}, self.c)
        self.assertEqual(policy["enzymes"], ["HaeIII"])
        self.assertEqual(policy["enzyme_definitions"], {"HaeIII": {"motif": "GGCC", "cut_offset": 2}})
        from prepare_restriction_fragments import scan_fasta
        fasta = self.root / "haeiii.fa"
        fasta.write_text(">toy\nAAggccTT\n")
        self.assertEqual(list(scan_fasta(fasta, policy["enzyme_definitions"], chunk_size=3))[0][2], [0, 4, 8])
        with self.assertRaisesRegex(ValueError, "requires HaeIII"):
            resolve_policy({"assay": "hic"}, {**self.c, "HIC_RESTRICTION_ENZYMES": "DpnII"})
        self.assertEqual(resolve_policy({"assay": "hic", "restriction_enzymes": "HaeIII"}, {**self.c, "HIC_RESTRICTION_ENZYMES": "DpnII"})["enzymes"], ["HaeIII"])
        with self.assertRaisesRegex(ValueError, "only valid"):
            resolve_policy({"assay": "microc", "blhic_input_state": "raw"}, self.c)
        with self.assertRaisesRegex(ValueError, "paired reads"):
            resolve_preprocessing({}, {**self.c, "FASTP_EXTRA_ARGS": "--merge"}, policy)
        self.assertNotIn("preprocessing", build_plan(self.fixture.c, self.fixture.sheet)["samples"]["row1"])  # ordinary hic

    def test_strict_input_errors_and_corrupt_gzip(self):
        original = self.r2.read_bytes()
        variants = [original.replace(b"@ab/2", b"@wrong/2"), original[:-1], original.replace(b"\n", b"\r\n"),
                    original.replace(b"C"*20, b"c"*20, 1), original.replace(b"\n+\n", b"\n?\n", 1),
                    original.replace(b"I", b"", 1), b"@"+b"x"*1000+b"\nA\n+\nI\n", b""]
        for value in variants:
            self.r2.write_bytes(value)
            with self.assertRaises(ValueError):
                inspect_pairs(self.r1, self.r2)
        bad = self.root / "bad.gz"
        bad.write_bytes(gzip.compress(original)[:-5])
        with self.assertRaises((EOFError, OSError)):
            inspect_pairs(self.r1, bad)

    def test_real_linker_semantics_one_and_two_threads(self):
        tool = identify_tool(self.c)
        a, b, g = PRESET["linker_a"], PRESET["linker_b"], "C"*20
        mutate = lambda word, count: ''.join(('T' if char != 'T' else 'C') if i < count else char for i, char in enumerate(word))
        extra1 = [("ba", g+b), ("bb", g+b), ("two", g+mutate(a, 2)), ("fuzzy", g+mutate(a, 1)),
                  ("exactfirst", g+mutate(a, 1)+"CCCC"+a), ("ambiguous", g+a+"CCCC"+b)]
        extra2 = [("ba", g+a), ("bb", g+b), ("two", g+b), ("fuzzy", g+b), ("exactfirst", g+b), ("ambiguous", g+b)]
        write_fastq(self.r1, self.input1+extra1, 1)
        write_fastq(self.r2, self.input2+extra2, 2)
        expected = {"ab": g, "one": g, "empty": g, "edge": "C"*15, "ba": g,
                    "two": g+mutate(a, 2), "fuzzy": g, "exactfirst": g+mutate(a, 1)+"CCCC"}
        for threads in (1, 2):
            attempt = self.root / f"real_{threads}"
            attempt.mkdir()
            subprocess.run([tool["path"], "-m", "1", "-k", "2", "-e", "1", "-l", "15", "-t", str(threads),
                            "-A", a, "-B", b, "-o", str(attempt), "-n", "reads", str(self.r1), str(self.r2)],
                           check=True, stdout=subprocess.DEVNULL)
            out1, out2 = attempt / "reads_1.valid.fastq", attempt / "reads_2.valid.fastq"
            pairs = inspect_pairs(out1, out2)
            self.assertEqual(read_fastq(out1), {name: (seq, 'I'*len(seq)) for name, seq in expected.items()})
            counts = parse_stats(attempt / "reads.trim.stat", inspect_pairs(self.r1, self.r2), pairs)
            self.assertEqual(counts["Valid PETs"], len(expected))

    def test_automatic_stage_runs_with_fastp_disabled_and_reuses(self):
        rt = self.runtime()
        p = PairProcessor(rt)
        outputs = p.ensure_trim("sample:row1")
        self.assertEqual(inspect_pairs(*outputs)["pairs"], 4)
        receipt = json.loads(Path(rt.plan["entities"]["sample:row1"]["paths"]["blhic_qc"]).read_text())
        self.assertEqual(receipt["original_raw_pairs"], 6)
        self.assertEqual(receipt["bridge_counts"]["too_short"], 1)
        self.assertFalse(receipt["fastp_performed"])
        self.assertNotIn("html", p.trim_outputs(rt.plan["entities"]["sample:row1"]))
        self.assertEqual(self.r1.read_bytes().count(b"@"), 6)
        with patch.object(rt, "command", side_effect=AssertionError("Must reuse completed stage")):
            self.assertEqual(p.ensure_trim("sample:row1"), outputs)

    def test_fastp_receives_bridge_outputs_before_alignment(self):
        rt = self.runtime({**self.c, "RUN_FASTP": "true"})
        p = PairProcessor(rt)
        original = Runtime.pipeline
        def command(runtime, commands, output, log_path, timeout=None):
            args = commands[0]
            self.calls.append(Path(args[0]).name)
            if args[0] == "fastp":
                for src, dst in (("-i", "-o"), ("-I", "-O")):
                    content = Path(args[args.index(src)+1]).read_bytes()
                    Path(args[args.index(dst)+1]).write_bytes(gzip.compress(content))
                Path(args[args.index("--html")+1]).write_text("<html>fixture</html>")
                Path(args[args.index("--json")+1]).write_text('{"summary":{}}')
            elif args[0] == "bwa-mem2":
                self.assertEqual(inspect_pairs(Path(args[-2]), Path(args[-1]))["pairs"], 4)
                from test_stage_performance import OptimizedWorkflowTests
                OptimizedWorkflowTests.emit(output, [fixtures.contact("mapped", 10, 30)])
                receipt = next(c for c in commands if len(c) > 1 and str(c[1]).endswith("pairs_receipt.py"))
                json_write(Path(receipt[receipt.index("--report")+1]), {"data_sha256": "fixture", "records": 1})
            else:
                return original(runtime, commands, output, log_path, timeout)
        with patch.object(Runtime, "pipeline", command):
            p.ensure_sorted("sample:row1")
        self.assertEqual(self.calls, ["trimLinker", "fastp", "bwa-mem2"])

    def test_pretrimmed_no_dependency_and_unknown_raw_statistics(self):
        rt = self.runtime({**self.c, "BLHIC_INPUT_STATE": "pretrimmed"})
        with patch.object(rt, "command", side_effect=AssertionError("No trimmer or fastp expected")):
            outputs = PairProcessor(rt).ensure_trim("sample:row1")
        self.assertEqual(inspect_pairs(*outputs)["pairs"], 6)
        receipt = json.loads(Path(rt.plan["entities"]["sample:row1"]["paths"]["blhic_qc"]).read_text())
        self.assertIsNone(receipt["original_raw_pairs"])
        self.assertIsNone(receipt["bridge_counts"])
        self.assertIsNone(receipt["tool"])

    def test_failure_does_not_publish_or_accept(self):
        rt = self.runtime()
        p = PairProcessor(rt)
        with patch.object(rt, "command", side_effect=RuntimeError("fixture interrupted")):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                p.ensure_trim("sample:row1")
        entity = rt.plan["entities"]["sample:row1"]
        self.assertFalse(Path(entity["paths"]["trim_r1"]).exists())
        self.assertFalse(p.reusable("sample:row1", "trim"))
        self.assertEqual(rt.store(entity).read("trim")["status"], "failed")

    def test_linker_science_changes_mapping_but_filter_changes_do_not(self):
        rt = self.runtime()
        old = PairProcessor(rt).keys["sample:row1"]
        filtered = self.runtime({**self.c, "MAPQ_THRESHOLD": "40"})
        new = PairProcessor(filtered).keys["sample:row1"]
        self.assertEqual(old["sorted"], new["sorted"])
        self.assertNotEqual(old["accepted"], new["accepted"])
        changed = self.runtime()
        changed.plan["samples"]["row1"]["preprocessing"]["preset"]["min_length"] = 16
        self.assertNotEqual(old["sorted"], PairProcessor(changed).keys["sample:row1"]["sorted"])
        moved = self.runtime()
        moved.versions["trimLinker"]["path"] = "/another/location/trimLinker"
        self.assertEqual(old["sorted"], PairProcessor(moved).keys["sample:row1"]["sorted"])

    def test_missing_raw_resume_checks_protocol_and_binary(self):
        rt = self.runtime()
        rt.plan["checkpoint_keys"] = PairProcessor(rt).keys
        rt.plan["tool_versions"] = rt.versions
        json_write(self.root / "out/run_metadata/run_plan.json", rt.plan)
        self.r1.unlink()
        self.r2.unlink()
        current = self.runtime({**self.c, "RESUME_FROM_CHECKPOINTS": "true"})
        self.assertEqual(PairProcessor(current).keys["sample:row1"]["sorted"], rt.plan["checkpoint_keys"]["sample:row1"]["sorted"])
        changed = self.runtime({**self.c, "RESUME_FROM_CHECKPOINTS": "true", "BLHIC_INPUT_STATE": "pretrimmed"})
        with self.assertRaisesRegex(ValueError, "preprocessing contract changed"):
            PairProcessor(changed)
        current.versions = copy.deepcopy(current.versions)
        current.versions["trimLinker"]["binary_sha256"] = "changed"
        with self.assertRaisesRegex(ValueError, "tool changed"):
            PairProcessor(current)

    def test_resources_storage_growth_and_mixed_summary(self):
        rt = self.runtime()
        self.assertGreaterEqual(reservation("fastp", self.c, rt.plan)[1], 16)
        self.assertEqual(reservation("fastp", self.fixture.c, build_plan(self.fixture.c, self.fixture.sheet))[1], 2)
        p = PairProcessor(rt)
        p.ensure_trim("sample:row1")
        # Mixed report keeps the original denominator and a consistent schema.
        self.fixture.rows.append({**self.fixture.rows[0], "sample": "ordinary", "condition": "ordinary", "hic_protocol": "restriction_hic", "restriction_enzymes": "DpnII"})
        self.fixture.write_sheet()
        mixed = build_plan(self.c, self.fixture.sheet)
        output = self.root / "report"
        write_plan_summary(mixed, output)
        import csv
        with (output / "microc2tracks_sample_summary.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["blhic_original_raw_pairs"], "6")
        self.assertEqual(rows[1]["blhic_original_raw_pairs"], "")
        from collections import namedtuple
        usage = namedtuple("Usage", "total used free")
        guard = StorageGuard({"out": str(self.root), "scratch": str(self.root)}, [], 0, 0)
        with patch("artifact_state.shutil.disk_usage", return_value=usage(1000, 0, 1000)):
            self.assertTrue(guard.reserve("a", {str(self.root): 600}))
            self.assertTrue(guard.grow("a", {str(self.root): 700}))
            self.assertFalse(guard.grow("b", {str(self.root): 400}))

    def test_pinned_tool_identity_and_refilter_dependency_gating(self):
        tool = identify_tool(self.c)
        self.assertEqual(scientific_tool_identity(tool)["source_commit"], "e120726d6440b24034f70bc3c51c17f351fef496")
        self.assertIn("binary_sha256", tool)
        with self.assertRaisesRegex(RuntimeError, "need trimLinker"):
            identify_tool({**self.c, "BLHIC_TRIMLINKER": "/missing/trimmer"})
        from blhic_preprocessing import raw_bl_rows
        plan = self.runtime().plan
        self.assertEqual(len(raw_bl_rows(plan)), 1)
        plan["workflow_mode"] = "refilter-deduplicated"
        self.assertEqual(raw_bl_rows(plan), [])

    def test_exit_zero_with_missing_mate_is_not_success(self):
        rt = self.runtime()
        p = PairProcessor(rt)
        original = rt.command
        def damaged(args, log_path, **kwargs):
            original(args, log_path, **kwargs)
            (Path(args[args.index("-o")+1]) / "reads_2.valid.fastq").unlink()
        with patch.object(rt, "command", damaged):
            with self.assertRaises(FileNotFoundError):
                p.ensure_trim("sample:row1")
        self.assertFalse(p.reusable("sample:row1", "trim"))
        self.assertFalse(Path(rt.plan["entities"]["sample:row1"]["paths"]["trim_r1"]).exists())

    def test_scheduled_bl_run_and_checkpoint_after_raw_cleanup(self):
        """Actual linker through the normal DAG; only mapping/contact tools mocked."""
        from mock_rc6_transport import adapt
        from pipeline import execute
        from test_stage_performance import OptimizedWorkflowTests
        helper = OptimizedWorkflowTests("test_one_joint_pass_per_library_outputs_and_final_barrier")
        helper.calls = []
        c = {**self.c, "BUILD_SAMPLE_MATRICES": "false", "KEEP_SORTED_PAIRS": "false"}
        # The mock contact at 10/30 must cross a real HaeIII fragment boundary.
        length = len(self.fixture.sequence)
        sequence = "C"*18 + "GGCC" + "A"*(length-22)
        self.fixture.fasta.write_text(f">scaffoldA\n{sequence}\n>plasmidZ\n{'A'*length}\n")
        versions = self.runtime(c).versions
        original = Runtime.pipeline
        def fake(runtime, commands, output, log_path, timeout=None):
            if Path(commands[0][0]).name == "trimLinker":
                helper.calls.append(commands[0])
                return original(runtime, commands, output, log_path, timeout)
            return helper.fake(runtime, commands, output, log_path, timeout)
        with patch("pipeline.tool_versions", return_value=versions), patch.object(Runtime, "pipeline", adapt(fake)), patch("pipeline.log"):
            plan = build_plan(c, self.fixture.sheet)
            self.assertEqual(execute(plan, self.root / "config"), 0)
            names = [Path(args[0]).name for args in helper.calls]
            self.assertEqual(names[:2], ["trimLinker", "bwa-mem2"])
            self.assertNotIn("fastp", names)
            paths = plan["entities"]["sample:row1"]["paths"]
            self.assertTrue(Path(paths["blhic_qc"]).is_file())
            self.assertIn("haeiii_two_step_v1", (self.root / "out/run_metadata/sample_manifest.tsv").read_text())
            self.assertFalse(Path(paths["trim_r1"]).exists())
            self.assertTrue(self.r1.is_file())
            self.r1.unlink()
            self.r2.unlink()
            helper.calls.clear()
            resumed = build_plan({**c, "RESUME_FROM_CHECKPOINTS": "true"}, self.fixture.sheet)
            self.assertEqual(execute(resumed, self.root / "config"), 0)
            self.assertEqual(helper.calls, [])

    def test_same_scratch_and_trimmed_root_reserves_all_copies(self):
        bulk = str(self.root / "bulk")
        rt = self.runtime({**self.c, "TMPDIR": bulk, "TRIMMED_OUTDIR": bulk, "LARGE_OUTPUT_STAGING": "scratch"})
        calls, original = [], rt.guard.grow
        def capture(token, requests):
            calls.append(requests.copy())
            return original(token, requests)
        with patch.object(rt.guard, "grow", capture):
            out = PairProcessor(rt).ensure_trim("sample:row1")
        self.assertEqual(inspect_pairs(*out)["pairs"], 4)
        amount = inspect_pairs(self.r1, self.r2)["uncompressed_bytes"]
        self.assertGreaterEqual(calls[0][bulk], amount + 2 * (int(amount * 1.02) + 65536))


if __name__ == "__main__":
    unittest.main()
