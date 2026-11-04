"""Focused scheduler/grouping/reuse regressions; mocks are not speed benchmarks."""
from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import test_integration_features as fixtures
from filter_contact_pairs import filter_contacts
from pairs_processing import PairProcessor, proof_path
from pipeline import Runtime, execute, subprocess_environment
from prepare_restriction_fragments import load_digest
from rc4_compatibility import RC4_HASHES
from reference_registry import read_sizes
from run_plan import build_plan
from stage_scheduler import StageScheduler, Task, resource_profile, size_gib
from workflow_common import hic_targets, identity, json_write, load_config, sha256, text_open


class SchedulerTests(unittest.TestCase):
    def test_nested_native_threads_are_bounded_without_overriding_bwa(self):
        with patch.dict("os.environ", {"OPENBLAS_NUM_THREADS": "64", "OMP_NUM_THREADS": "32"}):
            self.assertEqual(subprocess_environment(["cooler", "balance"])["OPENBLAS_NUM_THREADS"], "1")
            self.assertEqual(subprocess_environment(["pairtools", "dedup"])["OMP_NUM_THREADS"], "1")
            self.assertEqual(subprocess_environment(["bwa-mem2", "mem", "-t", "32"])["OMP_NUM_THREADS"], "32")

    def profile(self, cpu=4, memory=8, prepared=2):
        return {"cpu_budget": cpu, "memory_gib": memory, "max_prepared": prepared,
                "limits": {"alignment": 2, "dedup": 2, "fastp": 2, "report": 1}}

    def test_downstream_does_not_occupy_alignment_workers(self):
        release, started, lock, events = threading.Event(), [], threading.Lock(), []
        def dedup():
            self.assertTrue(release.wait(3), "Ready alignments were starved by downstream work")
        def align():
            with lock:
                started.append("alignment")
                if len(started) == 2:
                    release.set()
        tasks = [Task(f"dedup{i}", "dedup", str(i), dedup, priority=0) for i in range(2)]
        tasks += [Task(f"align{i}", "alignment", str(i), align, priority=1) for i in range(2)]
        scheduler = StageScheduler(self.profile(), threading.Event(), events.append)
        scheduler.run(tasks)
        self.assertEqual(len(started), 2)
        self.assertTrue(all(e["used_cpus"] <= 4 and e["used_memory_gib"] <= 8 for e in events))

    def test_cpu_memory_caps_backlog_and_dependency_readiness(self):
        events, tasks = [], []
        for i in range(5):
            tasks += [Task(f"trim{i}", "fastp", str(i), lambda: time.sleep(0.002), cpus=1, memory=2, prepares=str(i)),
                      Task(f"align{i}", "alignment", str(i), lambda: time.sleep(0.002), {f"trim{i}"}, cpus=1, memory=2,
                           priority=0, consumes=str(i))]
        scheduler = StageScheduler(self.profile(cpu=2, memory=4, prepared=1), threading.Event(), events.append)
        scheduler.run(tasks)
        self.assertTrue(all(e["used_cpus"] <= 2 and e["used_memory_gib"] <= 4 for e in events))
        self.assertTrue(any(e["reason"] == "prepared-input backlog" for e in events))
        self.assertEqual(scheduler.prepared, set())
        with self.assertRaisesRegex(ValueError, "exceeding"):
            scheduler.run([Task("too-big", "alignment", "x", lambda: None, cpus=3)])
        with self.assertRaisesRegex(ValueError, "Cyclic"):
            scheduler.run([Task("a", "alignment", "x", lambda: None, {"b"}), Task("b", "alignment", "x", lambda: None, {"a"})])

    def test_failure_cancels_owned_work_and_releases_credits(self):
        stop, active = threading.Event(), threading.Event()
        def sibling():
            active.set()
            self.assertTrue(stop.wait(3))
        def fail():
            self.assertTrue(active.wait(3))
            raise RuntimeError("synthetic stage failure")
        scheduler = StageScheduler(self.profile(), stop)
        with self.assertRaisesRegex(RuntimeError, "synthetic"):
            scheduler.run([Task("sibling", "alignment", "x", sibling), Task("failure", "dedup", "y", fail)])
        self.assertTrue(stop.is_set())
        self.assertEqual((scheduler.used_cpu, scheduler.used_memory, scheduler.active), (0, 0, {}))


class OptimizedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.IntegrationFeatures("test_technical_and_all_pools_use_disjoint_units_once")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.addCleanup(self.fixture.tearDown)
        self.root, self.c = self.fixture.root, dict(self.fixture.c)
        self.c.update(RUN_FASTP="true", BUILD_SAMPLE_MATRICES="false", HIC_OUTPUT_SCOPE="final", RUN_HIC="true",
                      RUN_MULTIQC="false", RUN_GLOBAL_MULTIQC="false", RUN_PRELIM_DOWNSTREAM="false",
                      KEEP_RAW_HIC="true", KEEP_JUICER_PAIRS="true", KEEP_TRIMMED_FASTQ="false")
        self.jar = self.root / "toy.jar"
        self.jar.write_bytes(b"mock, not an installed JAR")
        self.c["JUICER_TOOLS_JAR"] = str(self.jar)
        for i, row in enumerate(self.fixture.rows):
            row.update(dedup_group=f"library{i//2+1}", merge_group=f"biological{i//2+1}", final_pool="final_selected")
        self.fixture.write_sheet()
        self.calls = []

    @staticmethod
    def rows(path):
        with text_open(Path(path)) as handle:
            return [line.rstrip().split("\t") for line in handle if line.strip() and not line.startswith("#")]

    @staticmethod
    def emit(path, rows):
        with gzip.open(path, "wt") as handle:
            handle.write("## pairs format v1.0\n#genome_assembly: test\n#columns: " + " ".join(fixtures.COLS) + "\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")

    def fake(self, rt, commands, output, log_path, timeout=None):
        args = list(commands[0])
        self.calls.append(args)
        if args[0] == "fastp":
            for flag, source_flag in (("-o", "-i"), ("-O", "-I")):
                shutil.copyfile(args[args.index(source_flag)+1], args[args.index(flag)+1])
            Path(args[args.index("--html")+1]).write_text("<html>mock QC</html>")
            Path(args[args.index("--json")+1]).write_text('{"summary": {}}')
        elif args[0] == "bwa-mem2":
            i = int(log_path.name.split(".")[0][3:])
            self.emit(output, [fixtures.contact(f"common{i}", 10, 30), fixtures.contact(f"unique{i}", 30+i*10, 50+i*10)])
        elif args[:2] == ["pairtools", "merge"]:
            sources = [Path(p) for p in args[2:] if Path(p).is_file()]
            data = [row for source in sources for row in self.rows(source)]
            data.sort(key=lambda r: (r[1], r[3], int(r[2]), int(r[4])))
            self.emit(output, data)
        elif args[:2] == ["pairtools", "dedup"]:
            data, seen, kept = self.rows(args[-1]), set(), []
            for row in data:
                coordinate = tuple(row[1:7])
                if coordinate not in seen:
                    kept.append(row)
                    seen.add(coordinate)
            self.emit(output, kept)
            Path(args[args.index("--output-stats")+1]).write_text(f"total\t{len(data)}\ntotal_nodups\t{len(kept)}\n")
        elif len(args)>1 and str(args[1]).endswith("filter_contact_pairs.py"):
            def get(flag): return Path(args[args.index(flag)+1])
            data = io.StringIO()
            qc = filter_contacts(get("--input"), data, json.loads(get("--policy").read_text()),
                read_sizes(get("--chrom-sizes"), "reference"), set(read_sizes(get("--matrix-chrom-sizes"), "matrix")),
                load_digest(get("--digest")) if "--digest" in args else None)
            with gzip.open(output, "wt") as handle:
                handle.write(data.getvalue())
            json_write(get("--qc"), qc)
        elif args[0] == "pairix":
            Path(args[-1] + ".px2").write_bytes(b"mock index")
        elif args[:2] == ["pairtools", "stats"]:
            Path(output).write_text(f"total\t{len(self.rows(args[-1]))}\n")
        elif args[:3] == ["cooler", "cload", "pairix"]:
            Path(args[-1]).write_text(str(len(self.rows(args[-2]))))
        elif args[:2] == ["cooler", "zoomify"]:
            shutil.copyfile(args[-1], args[args.index("-o")+1])
        elif args[:2] == ["cooler", "info"]:
            Path(output).write_text(json.dumps({"sum": int(Path(args[-1]).read_text())}))
        elif args[0] == "java" and "pre" in args:
            with gzip.open(args[-3], "rt") as handle:
                count = sum(1 for _ in handle)
            Path(args[-2]).write_bytes(b"HIC\0" + str(count).encode())
        elif args[:2] == ["cooler", "ls"] or (args[0] == "java" and "addNorm" in args):
            pass
        else:
            raise AssertionError(args)

    def execution(self):
        import contextlib
        stack = contextlib.ExitStack()
        stack.enter_context(patch("pipeline.tool_versions", return_value={"pairtools": "1.1.3", "cooler": "mock", "java": "mock"}))
        stack.enter_context(patch.object(Runtime, "pipeline", lambda rt, *args, **kwargs: self.fake(rt, *args, **kwargs)))
        stack.enter_context(patch("pipeline.log"))
        return stack

    def test_one_joint_pass_per_library_outputs_and_final_barrier(self):
        with self.execution():
            plan = build_plan(self.c, self.fixture.sheet)
            self.assertEqual(execute(plan, self.root / "config"), 0)
            entities = plan["entities"]
            dedup = [args for args in self.calls if args[:2] == ["pairtools", "dedup"]]
            self.assertEqual(len(dedup), 2)
            self.assertTrue(all("dedup_groups" in args[-1] for args in dedup))
            self.assertEqual(sum(args[0] == "pairix" for args in self.calls), 3)
            self.assertEqual(hic_targets(plan), ["pool:final_selected"])
            self.assertTrue(all(not Path(e["paths"]["norm_mcool"]).exists() for e in entities.values() if e["kind"] == "sample"))
            for n in (1, 2):
                biological, library = entities[f"pool:biological{n}"], entities[f"dedup:library{n}"]
                self.assertEqual(biological["paths"]["accepted"], library["paths"]["accepted"])
                self.assertEqual(int(Path(biological["paths"]["norm_mcool"]).read_text()), 3)
                self.assertFalse(Path(biological["paths"]["norm_hic"]).exists())
            final = entities["pool:final_selected"]
            self.assertEqual(len(self.rows(final["paths"]["accepted"])), 6)
            self.assertEqual(sum(row[2] == "10" and row[4] == "30" for row in self.rows(final["paths"]["accepted"])), 2)
            self.assertEqual(Path(final["paths"]["norm_hic"]).read_bytes(), b"HIC\x006")
            self.assertTrue(all("--balance-args" in args for args in self.calls if args[:2] == ["cooler", "zoomify"]))
            events = [json.loads(line) for line in (self.root / "out/run_metadata/stage_events.jsonl").read_text().splitlines()]
            final_start = next(e["time"] for e in events if e["task"] == "pool:final_selected/accepted" and e["state"] == "running")
            for n in (1, 2):
                finished = next(e["time"] for e in events if e["task"] == f"pool:biological{n}/finished" and e["state"] == "complete")
                self.assertGreaterEqual(final_start, finished)
            self.calls.clear()
            self.assertEqual(execute(build_plan(self.c, self.fixture.sheet), self.root / "config"), 0)
            self.assertEqual(self.calls, [])

    def test_independent_libraries_within_one_biological_are_not_deduplicated_together(self):
        for row in self.fixture.rows:
            row.update(biological_replicate="1", merge_group="one_biological", technical_replicate=str(self.fixture.rows.index(row)+1))
        self.fixture.write_sheet()
        with self.execution():
            plan = build_plan(self.c, self.fixture.sheet)
            self.assertEqual(execute(plan, self.root / "config"), 0)
        biological = plan["entities"]["pool:one_biological"]
        self.assertNotIn("alias_of", biological)
        self.assertEqual(len(self.rows(biological["paths"]["accepted"])), 6)
        self.assertEqual(sum(args[:2] == ["pairtools", "dedup"] for args in self.calls), 2)

    def test_partial_final_membership_rejected(self):
        self.fixture.rows[0]["final_pool"] = ""
        self.fixture.write_sheet()
        with self.assertRaisesRegex(ValueError, "complete biological"):
            build_plan(self.c, self.fixture.sheet)

    def test_single_run_biological_output_and_cross_biological_guard(self):
        self.fixture.rows = self.fixture.rows[:1]
        self.fixture.rows[0]["final_pool"] = ""
        self.fixture.write_sheet()
        plan = build_plan({**self.c, "BUILD_BIOLOGICAL_MATRICES": "true"}, self.fixture.sheet)
        self.assertEqual(plan["pools"], {"biological1": "pool:biological1"})
        self.assertEqual(plan["entities"]["pool:biological1"]["alias_of"], "dedup:library1")
        second = {**self.fixture.rows[0], "sample": "independent_bio", "biological_replicate": "2", "merge_group": "biological2"}
        self.fixture.rows.append(second)
        self.fixture.write_sheet()
        with self.assertRaisesRegex(ValueError, "crosses biological"):
            build_plan(self.c, self.fixture.sheet)

    def test_server_profile_resource_budget(self):
        config = load_config(ROOT / "config/server_140cpu_500gb.conf")
        with patch("stage_scheduler.available_cpus", return_value=140), patch("stage_scheduler.available_memory_gib", return_value=500):
            profile = resource_profile(config)
        self.assertEqual(profile["cpu_budget"], 120)
        self.assertEqual(profile["memory_gib"], 400)
        self.assertEqual(profile["limits"]["alignment"], 2)
        self.assertEqual(config["THREADS_ALIGN"], "32")
        self.assertEqual(config["HIC_OUTPUT_SCOPE"], "final")
        self.assertEqual(size_gib("64g"), 64)
        with patch("stage_scheduler.available_cpus", return_value=140), patch("stage_scheduler.available_memory_gib", return_value=100):
            self.assertEqual(resource_profile(config)["memory_gib"], 80)

    def test_verified_rc4_sorted_checkpoint_reused_without_rewriting_originals(self):
        c = {**self.c, "RUN_FASTP": "false", "RUN_HIC": "false", "REUSE_RC4_OUTPUTS": "false"}
        plan = build_plan(c, self.fixture.sheet)
        old = Runtime(plan, self.root / "config", versions={"pairtools": "1.1.3", "bwa-mem2": "mock"})
        old.code_hashes.update(RC4_HASHES)
        processor = PairProcessor(old)
        with patch.object(Runtime, "pipeline", lambda rt, *args, **kwargs: self.fake(rt, *args, **kwargs)):
            processor.ensure_sorted("sample:row1")
        paths = plan["entities"]["sample:row1"]["paths"]
        originals = [Path(paths["sorted"]), proof_path(Path(paths["sorted"])), Path(paths["done"]) / "sorted.json"]
        before = {str(path): (sha256(path), path.stat().st_mtime_ns) for path in originals}
        json_write(self.root / "out/run_metadata/source_hashes.json", RC4_HASHES)
        current_plan = build_plan({**c, "REUSE_RC4_OUTPUTS": "true"}, self.fixture.sheet)
        current = Runtime(current_plan, self.root / "config", versions=old.versions)
        resumed = PairProcessor(current)
        self.assertTrue(current.rc4_hashes)
        self.assertTrue(resumed.reusable("sample:row1", "sorted"))
        with patch.object(Runtime, "pipeline", side_effect=AssertionError("Must not remap")):
            resumed.ensure_sorted("sample:row1")
        self.assertEqual(before, {str(path): (sha256(path), path.stat().st_mtime_ns) for path in originals})
        changed = Runtime(build_plan({**c, "REUSE_RC4_OUTPUTS": "true", "BWA_EXTRA_ARGS": "-changed"}, self.fixture.sheet), self.root / "config", versions=old.versions)
        self.assertFalse(PairProcessor(changed).reusable("sample:row1", "sorted"))
        marker = originals[-1]
        record = json.loads(marker.read_text())
        json_write(marker, {**record, "status": "failed"})
        self.assertFalse(resumed.reusable("sample:row1", "sorted"))
        json_write(marker, record)
        changed_tools = Runtime(current_plan, self.root / "config", versions={**old.versions, "pairtools": "changed version"})
        self.assertFalse(PairProcessor(changed_tools).reusable("sample:row1", "sorted"))
        self.fixture.fastq.write_text("@changed\nACGT\n+\nIIII\n")
        changed_inputs = Runtime(build_plan({**c, "REUSE_RC4_OUTPUTS": "true"}, self.fixture.sheet), self.root / "config", versions=old.versions)
        self.assertFalse(PairProcessor(changed_inputs).reusable("sample:row1", "sorted"))

    def test_cancellation_does_not_kill_unrelated_process(self):
        plan = build_plan({**self.c, "RUN_HIC": "false"}, self.fixture.sheet)
        rt = Runtime(plan, self.root / "config", versions={"pairtools": "1.1.3"})
        marker, errors = self.root / "started", []
        unrelated = subprocess.Popen([sys.executable, "-B", "-c", "import time; time.sleep(20)"])
        def child():
            try:
                rt.command([sys.executable, "-B", "-c", f"from pathlib import Path; import time; Path({str(marker)!r}).write_text('ready'); time.sleep(20)"], self.root / "child.log")
            except BaseException as exc:
                errors.append(exc)
        worker = threading.Thread(target=child)
        try:
            worker.start()
            deadline = time.monotonic() + 5
            while not marker.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertTrue(marker.exists())
            rt.stop.set()
            worker.join(8)
            self.assertFalse(worker.is_alive())
            self.assertTrue(any(isinstance(exc, InterruptedError) for exc in errors))
            self.assertIsNone(unrelated.poll())
        finally:
            rt.stop.set()
            worker.join(8)
            unrelated.terminate()
            unrelated.wait(timeout=5)


if __name__ == "__main__":
    unittest.main()
