"""Small output-scope regressions; external scientific commands are mocked."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import test_refilter_mode as refilter_tests
from matrix_products import export_juicer
from pipeline import Runtime, execute, tool_versions
from workflow_common import DEFAULTS, builds_hic, hic_targets, sha256, validate_config


class HicOutputScopeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = refilter_tests.RefilterModeTests("test_minimal_sheet_needs_no_fastqs_or_bwa_index")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.fixture.c["HIC_OUTPUT_SCOPE"] = "merged"

    def add_technical_rows(self, count=2):
        f = self.fixture
        f.rows = []
        for i in range(count):
            source = self.root / f"input_{i}.pairs.gz"
            f.emit_source(source, [[f"input{i}_" + row[0], *row[1:]] for row in f.records])
            f.rows.append({"sample": f"input_{i}", "assay": "hic", "reference_genome": "toy", "deduplicated_pairs": str(source),
                "input_members": f"independent_library_{i}", "condition": "condition", "biological_replicate": str(i//2+1),
                "technical_replicate": str(i%2+1), "merge_group": f"technical_{i//2+1}"})
        f.write_sheet()

    def mock_execution(self, calls, *, fail_normalization=False):
        fake = self.fixture.fake_pipeline(calls)
        def pipeline(rt, commands, output, log_path, timeout=None):
            if fail_normalization and commands[0][0] == "java" and "addNorm" in commands[0]:
                raise RuntimeError("synthetic addNorm failure")
            return fake(rt, commands, output, log_path, timeout)
        import contextlib
        stack = contextlib.ExitStack()
        stack.enter_context(patch("pipeline.tool_versions", return_value={"pairtools": "1.1.3", "cooler": "mock", "java": "mock"}))
        stack.enter_context(patch.object(Runtime, "pipeline", pipeline))
        stack.enter_context(patch("pipeline.log"))
        return stack

    def test_default_scope_validation_and_explicit_merge_targets(self):
        self.assertEqual(DEFAULTS["HIC_OUTPUT_SCOPE"], "merged")
        valid = {**DEFAULTS, "MCOOL_RESOLUTIONS": DEFAULTS["RESOLUTIONS"], "HIC_RESOLUTIONS": DEFAULTS["RESOLUTIONS"]}
        validate_config(valid)
        with self.assertRaisesRegex(ValueError, "HIC_OUTPUT_SCOPE"):
            validate_config({**valid, "HIC_OUTPUT_SCOPE": "typo"})
        c = dict(DEFAULTS)
        self.assertFalse(builds_hic({"kind": "sample"}, c))
        self.assertTrue(builds_hic({"kind": "pool"}, c))
        self.assertFalse(builds_hic({"kind": "dedup"}, c))
        self.assertTrue(builds_hic({"kind": "dedup", "build_matrices": True}, c))
        self.assertTrue(builds_hic({"kind": "sample"}, {**c, "HIC_OUTPUT_SCOPE": "all"}))
        self.assertFalse(builds_hic({"kind": "sample"}, {**c, "HIC_OUTPUT_SCOPE": "all", "BUILD_SAMPLE_MATRICES": "false"}))
        self.assertFalse(builds_hic({"kind": "pool"}, {**c, "RUN_HIC": "false"}))
        self.assertFalse(builds_hic({"kind": "external"}, c))

    def test_only_merged_outputs_receive_java_exports_tracks_and_public_copies(self):
        self.add_technical_rows(4)
        public = self.root / "public_hic"
        plan = self.fixture.plan(RUN_MERGE_ALL_SAMPLES="true", PUBLIC_HIC_DIR=str(public), PUBLIC_HIC_BASE_URL="https://example.invalid/hic")
        self.assertEqual(set(hic_targets(plan)), {"pool:technical_1", "pool:technical_2", "pool:all_samples_merged"})
        self.assertEqual(plan["hic_targets"], hic_targets(plan))
        calls, exports = [], []
        def export(source, destination):
            # Technical/all-sample exports can only happen after the source pool exists.
            self.assertTrue(Path(source).is_file())
            exports.append(str(source))
            return export_juicer(source, destination)
        with self.mock_execution(calls), patch("matrix_products.export_juicer", side_effect=export):
            self.assertEqual(execute(plan, self.root / "config"), 0)
        for entity in plan["entities"].values():
            p = entity["paths"]
            self.assertTrue(Path(p["accepted"]).is_file())
            self.assertTrue(Path(p["norm_mcool"]).is_file())
            if entity["kind"] == "sample":
                for label in ("raw_hic", "norm_hic", "juicer"):
                    self.assertFalse(Path(p[label]).exists())
                self.assertFalse((Path(p["done"]) / "juicer.json").exists())
                self.assertFalse((Path(p["matrix_dir"]) / (entity["name"] + ".ucsc.hic.track.txt")).exists())
                self.assertFalse((public / Path(p["norm_hic"]).name).exists())
            else:
                self.assertTrue(Path(p["norm_hic"]).is_file())
                self.assertTrue((public / Path(p["norm_hic"]).name).is_file())
                self.assertTrue((Path(p["matrix_dir"]) / (entity["name"] + ".ucsc.hic.track.txt")).is_file())
        self.assertEqual(set(exports), {plan["entities"][key]["paths"]["accepted"] for key in hic_targets(plan)})
        self.assertEqual(sum(command[0] == "java" and "pre" in command for command in calls), 3)
        self.assertEqual(sum(command[0] == "java" and "addNorm" in command for command in calls), 3)

    def test_singleton_has_no_java_requirement_and_report_marks_intentional_skip(self):
        plan = self.fixture.plan(JUICER_TOOLS_JAR="")
        self.assertEqual(hic_targets(plan), [])
        self.assertTrue(any("No .hic outputs" in warning for warning in plan["warnings"]))
        queries = []
        def query(args, **kwargs):
            queries.append(args[0])
            self.assertNotEqual(args[0], "java")
            return SimpleNamespace(stdout="pairtools, version 1.1.3" if args[0] == "pairtools" else "mock", stderr="")
        def which(name):
            return None if name in {"java", "bwa-mem2", "fastp", "cooltools", "multiqc"} else name
        with patch("pipeline.shutil.which", side_effect=which), patch("pipeline.subprocess.run", side_effect=query):
            self.assertNotIn("java", tool_versions(plan))
        calls = []
        with self.mock_execution(calls):
            self.assertEqual(execute(plan, self.root / "config"), 0)
        self.assertFalse(any(command[0] == "java" for command in calls))
        with (self.root / "out/final_report/microc2tracks_sample_summary.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["status"], "success")
        self.assertEqual(row["hic_requested"], "false")
        self.assertEqual(row["hic_output_status"], "not_requested")
        self.assertEqual(row["norm_hic_path"], "")

    def test_scope_changes_reuse_pairs_coolers_and_preserve_existing_hic_products(self):
        calls = []
        f = self.fixture
        with self.mock_execution(calls):
            first = f.plan(KEEP_RAW_HIC="true", KEEP_JUICER_PAIRS="true")
            self.assertEqual(execute(first, self.root / "config"), 0)
            paths = first["entities"]["sample:existing_library"]["paths"]
            identities = {label: (sha256(Path(paths[label])), Path(paths[label]).stat().st_mtime_ns)
                for label in ("accepted", "norm_mcool")}
            calls.clear()
            all_outputs = f.plan(HIC_OUTPUT_SCOPE="all", KEEP_RAW_HIC="true", KEEP_JUICER_PAIRS="true")
            self.assertEqual(execute(all_outputs, self.root / "config"), 0)
            self.assertTrue(calls)
            self.assertTrue(all(command[0] == "java" for command in calls))
            self.assertEqual(identities, {label: (sha256(Path(paths[label])), Path(paths[label]).stat().st_mtime_ns) for label in identities})
            hic_identities = {label: (sha256(Path(paths[label])), Path(paths[label]).stat().st_mtime_ns)
                for label in ("norm_hic", "raw_hic", "juicer")}
            calls.clear()
            self.assertEqual(execute(f.plan(), self.root / "config"), 0)
            self.assertEqual(calls, [])
            self.assertEqual(hic_identities, {label: (sha256(Path(paths[label])), Path(paths[label]).stat().st_mtime_ns) for label in hic_identities})
        with (self.root / "out/final_report/microc2tracks_sample_summary.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["hic_output_status"], "existing_preserved_not_requested")

    def test_disabled_hic_keeps_merged_cooler_outputs_without_java(self):
        self.add_technical_rows()
        plan = self.fixture.plan(RUN_HIC="false", JUICER_TOOLS_JAR="")
        self.assertEqual(hic_targets(plan), [])
        calls = []
        with self.mock_execution(calls):
            self.assertEqual(execute(plan, self.root / "config"), 0)
        self.assertFalse(any(command[0] == "java" for command in calls))
        for entity in plan["entities"].values():
            self.assertTrue(Path(entity["paths"]["norm_mcool"]).is_file())
            self.assertFalse(Path(entity["paths"]["norm_hic"]).exists())

    def test_merged_hic_failure_preserves_sample_outputs_and_resumes_only_normalization(self):
        self.add_technical_rows()
        plan = self.fixture.plan()
        calls = []
        with self.mock_execution(calls, fail_normalization=True):
            with self.assertRaisesRegex(RuntimeError, "synthetic addNorm failure"):
                execute(plan, self.root / "config")
        self.assertEqual(json.loads((self.root / "out/run_metadata/run_status.json").read_text())["status"], "failed")
        saved = {key: {label: sha256(Path(entity["paths"][label])) for label in ("accepted", "norm_mcool")}
            for key, entity in plan["entities"].items() if entity["kind"] == "sample"}
        calls.clear()
        with self.mock_execution(calls):
            self.assertEqual(execute(self.fixture.plan(), self.root / "config"), 0)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][0], "java")
        self.assertIn("addNorm", calls[0])
        self.assertEqual(saved, {key: {label: sha256(Path(plan["entities"][key]["paths"][label])) for label in values} for key, values in saved.items()})


if __name__ == "__main__":
    unittest.main()
