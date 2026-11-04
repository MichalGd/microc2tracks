"""Bounded-processing equivalence, checkpoint recovery, and publication safety."""
import gzip
import io
import json
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import test_integration_features as fixtures
from artifact_state import StorageGuard
from batched_filter import filter_batched
from filter_contact_pairs import filter_contacts
from juicer_export import project
from matrix_products import export_juicer
from pairs_processing import PairProcessor
from pipeline import Runtime, execute
from run_plan import build_plan
from workflow_common import identity, json_write, signature


class SafetyTests(fixtures.IntegrationFeatures):
    # Reuse only the fixture setup, not the inherited test collection.
    def test_batched_filter_matches_scalar_for_protocols_cutoffs_and_reordered_columns(self):
        rng = random.Random(15)
        records = [fixtures.contact(f"r{i}", rng.randint(1, 120), rng.randint(1, 120),
            c2=rng.choice(["scaffoldA", "plasmidZ"]), s1=rng.choice(["+", "-"]), s2=rng.choice(["+", "-"]),
            kind=rng.choice(["UU", "UR", "RU", "DD"]), q1=rng.choice([0, 29, 30, 60])) for i in range(1200)]
        source = self.root / "pairs"
        fixtures.write_pairs(source, records)
        sizes = {"scaffoldA": len(self.sequence), "plasmidZ": len(self.sequence)}
        digest = {c: [0, 2, 6, 19, length] for c, length in sizes.items()}
        for assay in ("microc", "hic"):
            for cutoff in (0, 1, 10, 1000):
                for adjacent in (False, True):
                    policy = {**self.policy(assay), "min_cis_dist": cutoff, "remove_adjacent_inward": adjacent}
                    old, new = io.StringIO(), io.StringIO()
                    args = (policy, sizes, set(sizes), digest if assay == "hic" else None)
                    a = filter_contacts(source, old, *args)
                    b = filter_batched(source, new, *args, batch_records=37, batch_bytes=4096)
                    self.assertEqual(a, b)
                    self.assertEqual(old.getvalue(), new.getvalue())
        lines = source.read_text().splitlines()
        columns = fixtures.COLS[::-1]
        source.write_text("#columns: " + " ".join(columns) + "\n" + "\n".join("\t".join(row.split("\t")[::-1]) for row in lines if not row.startswith("#")) + "\n")
        a, b = io.StringIO(), io.StringIO()
        self.assertEqual(filter_contacts(source, a, self.policy(), sizes, set(sizes)), filter_batched(source, b, self.policy(), sizes, set(sizes), batch_records=11))
        self.assertEqual(a.getvalue(), b.getvalue())

    def test_optimized_export_matches_rc5_projection(self):
        source, old = self.root / "pairs", self.root / "old.gz"
        fixtures.write_pairs(source, [fixtures.contact("one", 1, 11), fixtures.contact("edge", 10, 20)])
        export_juicer(source, old)
        output, checks = io.StringIO(), []
        self.assertEqual(project(source, output, batch_bytes=1, check=lambda: checks.append(True)), 2)
        with gzip.open(old, "rt") as inp:
            self.assertEqual(output.getvalue(), inp.read())
        self.assertGreaterEqual(len(checks), 2)

    def test_publication_journal_recovers_sigkill_window_without_adopting_success(self):
        store = self.store()
        target = self.root / "outputs/data"
        target.parent.mkdir(parents=True)
        target.write_text("newly published")
        json_write(store.metadata / ".bundle.publication.json", {"generation": "dead", "signature": "key", "outputs": {"data": identity(target)}})
        self.assertFalse(store.complete("bundle", "key", {"data": target}))
        store.run("bundle", "key", {"data": target}, lambda t, w: t["data"].write_text("retry"))
        self.assertTrue(store.complete("bundle", "key", {"data": target}))
        self.assertFalse((store.metadata / ".bundle.publication.json").exists())

    def test_per_device_disk_reservations_share_capacity(self):
        guard = StorageGuard({"out": str(self.root), "scratch": str(self.root / "scratch")}, [], 0, 0)
        from collections import namedtuple
        Usage = namedtuple("Usage", "total used free")
        with patch("artifact_state.shutil.disk_usage", return_value=Usage(1000, 0, 1000)):
            self.assertTrue(guard.reserve("a", {str(self.root): 600}))
            self.assertFalse(guard.reserve("b", {str(self.root / "scratch"): 500}))
            guard.release("a")
            self.assertTrue(guard.reserve("b", {str(self.root / "scratch"): 500}))

    def test_resume_with_removed_fastqs_uses_recorded_lineage_and_rejects_science_changes(self):
        c = {**self.c, "RUN_FASTP": "false", "RUN_HIC": "false"}
        old = build_plan(c, self.sheet)
        rt = Runtime(old, self.root / "config", versions={"pairtools": "1.1.3"})
        keys = PairProcessor(rt).keys
        old["checkpoint_keys"] = keys
        json_write(self.root / "out/run_metadata/run_plan.json", old)
        self.fastq.unlink()
        resumed = build_plan({**c, "RESUME_FROM_CHECKPOINTS": "true"}, self.sheet)
        current = Runtime(resumed, self.root / "config", versions=rt.versions)
        self.assertEqual(PairProcessor(current).keys["sample:row1"]["sorted"], keys["sample:row1"]["sorted"])
        changed = Runtime(build_plan({**c, "RESUME_FROM_CHECKPOINTS": "true", "BWA_EXTRA_ARGS": "-changed"}, self.sheet), self.root / "config", versions=rt.versions)
        with self.assertRaisesRegex(ValueError, "contract changed"):
            PairProcessor(changed)

    def test_completed_library_resumes_after_raw_and_sorted_inputs_removed(self):
        from test_stage_performance import OptimizedWorkflowTests
        fixture = OptimizedWorkflowTests("test_one_joint_pass_per_library_outputs_and_final_barrier")
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        self.addCleanup(fixture.doCleanups)
        with fixture.execution():
            plan = build_plan(fixture.c, fixture.fixture.sheet)
            self.assertEqual(execute(plan, fixture.root / "config"), 0)
            fixture.fixture.fastq.unlink()
            for e in plan["entities"].values():
                if e["kind"] == "sample":
                    Path(e["paths"]["sorted"]).unlink(missing_ok=True)
            fixture.calls.clear()
            resumed = build_plan({**fixture.c, "RESUME_FROM_CHECKPOINTS": "true"}, fixture.fixture.sheet)
            self.assertEqual(execute(resumed, fixture.root / "config"), 0)
            self.assertEqual(fixture.calls, [])

    def test_exact_rc5_sorted_contract_migrates_without_rewriting_payload_or_proof(self):
        from rc5_compatibility import RC5_HASHES
        from mock_rc6_transport import adapt
        from test_stage_performance import OptimizedWorkflowTests
        from workflow_common import sha256
        c = {**self.c, "RUN_FASTP": "false", "RUN_HIC": "false"}
        plan = build_plan(c, self.sheet)
        old = Runtime(plan, self.root / "config", versions={"pairtools": "1.1.3", "bwa-mem2": "fixture"})
        old.code_hashes.update(RC5_HASHES)
        def rc5_key(stage, parameters, tools=(), code=()):
            return signature({"stage": stage, "schema": 1, "parameters": parameters,
                "tools": {name: old.versions.get(name, "unknown") for name in tools},
                "libraries": {name: old.versions.get(name) for name in (["numpy", "scipy", "pandas", "pysam"] if "pairtools" in tools else ["numpy", "pandas", "h5py"] if "cooler" in tools else [])},
                "code": {name: old.code_hashes[name] for name in ["workflow_common.py", *code]}})
        # This fixture deliberately constructs the original RC5 key formula.
        with patch.object(old, "key", side_effect=rc5_key):
            processor = PairProcessor(old)
            def fake(rt, commands, output, log_path, timeout=None):
                OptimizedWorkflowTests.emit(output, [fixtures.contact("one", 10, 30)])
            with patch.object(Runtime, "pipeline", adapt(fake)):
                processor.ensure_sorted("sample:row1")
        paths = plan["entities"]["sample:row1"]["paths"]
        protected = [Path(paths["sorted"]), Path(paths["sorted"] + ".provenance.json"), Path(paths["done"]) / "sorted.json"]
        before = [(sha256(path), path.stat().st_mtime_ns) for path in protected]
        json_write(self.root / "out/run_metadata/source_hashes.json", RC5_HASHES)
        current = Runtime(build_plan(c, self.sheet), self.root / "config", versions=old.versions)
        resumed = PairProcessor(current)
        self.assertTrue(resumed.reusable("sample:row1", "sorted"))
        with patch.object(Runtime, "pipeline", side_effect=AssertionError("Must not remap")):
            resumed.ensure_sorted("sample:row1")
        self.assertEqual(before, [(sha256(path), path.stat().st_mtime_ns) for path in protected])


# Prevent inherited fixture tests from being executed a second time.
for name in list(vars(fixtures.IntegrationFeatures)):
    if name.startswith("test_") and name not in SafetyTests.__dict__:
        setattr(SafetyTests, name, None)


if __name__ == "__main__":
    unittest.main()
