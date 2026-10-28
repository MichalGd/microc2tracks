"""Small refilter/import regressions. Mocked commands test orchestration, not real tools."""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import shutil
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))
import test_integration_features as fixtures
from filter_contact_pairs import filter_contacts
from pairs_processing import PairProcessor, proof_path, validate_pool_inputs
from pipeline import Runtime, execute, tool_versions
from prepare_restriction_fragments import load_digest
from reference_registry import BWA_MEM2_INDEX_SUFFIXES, read_sizes
from refilter_inputs import normalize, source_identity
from run_plan import build_plan
from workflow_common import json_write, sha256, text_open


class RefilterModeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.IntegrationFeatures("test_profiles_generic_dovetail_legacy_and_explicit_enzyme_requirement")
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        self.root = self.fixture.root
        self.c = {**self.fixture.c, "HIC_MIN_CIS_DIST": "10", "RUN_MULTIQC": "false", "RUN_GLOBAL_MULTIQC": "false",
                  "RUN_FASTP": "true", "HIC_OUTPUT_SCOPE": "all", "KEEP_DEDUP_PAIRS": "false", "KEEP_SORTED_PAIRS": "false", "KEEP_TRIMMED_FASTQ": "false"}
        self.jar = self.root / "fake.jar"
        self.jar.write_bytes(b"mock fixture; not a real JAR")
        self.c["JUICER_TOOLS_JAR"] = str(self.jar)
        self.source = self.root / "old_library.dedup.pairs.gz"
        self.records = [fixtures.contact("dangling", 10, 11), fixtures.contact("circle", 12, 13, s1="-", s2="+"),
                        fixtures.contact("short", 1, 3, s1="-", s2="+"), fixtures.contact("edge", 1, 11),
                        fixtures.contact("above", 1, 12), fixtures.contact("trans", 20, 1, c2="plasmidZ"),
                        fixtures.contact("same_coordinates", 1, 11), fixtures.contact("swapped", 11, 1, s1="-", s2="+", kind="UR", q1=40, q2=60),
                        fixtures.contact("low_mapq", 1, 20, q1=29)]
        self.emit_source(self.source, self.records)
        self.index = Path(str(self.source) + ".px2")
        self.index.write_bytes(b"old index: must remain untouched")
        self.sheet = self.root / "refilter.csv"
        self.rows = [{"sample": "existing_library", "assay": "hic", "reference_genome": "toy", "deduplicated_pairs": str(self.source),
                      "input_members": "technical_A;technical_B"}]
        self.write_sheet()

    def emit_source(self, path, rows):
        with gzip.open(path, "wt") as handle:
            handle.write(f"## pairs format v1.0\n#genome_assembly: test\n#chromsize: scaffoldA {len(self.fixture.sequence)}\n#columns: " + " ".join(fixtures.COLS) + "\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")

    def write_sheet(self):
        with self.sheet.open("w", newline="") as handle:
            fields = list(dict.fromkeys(key for row in self.rows for key in row))
            writer = csv.DictWriter(handle, fields, lineterminator="\r\n")
            writer.writeheader(); writer.writerows(self.rows)

    def plan(self, **overrides):
        return build_plan({**self.c, **overrides}, self.sheet, mode="refilter-deduplicated")

    def test_minimal_sheet_needs_no_fastqs_or_bwa_index(self):
        for suffix in BWA_MEM2_INDEX_SUFFIXES:
            Path(str(self.fixture.fasta) + suffix).unlink()
        plan = self.plan()
        self.assertEqual(plan["workflow_mode"], "refilter-deduplicated")
        self.assertEqual(plan["settings"]["WORKFLOW_MODE"], "refilter-deduplicated")
        self.assertEqual(plan["pools"], {})
        self.assertEqual(plan["samples"]["existing_library"]["fastq_r1"], "")
        self.assertNotIn("bwa_index", plan["references"]["toy"]["identity"])
        self.assertEqual(self.plan(DEFAULT_REFERENCE_ID="unused")["samples"]["existing_library"]["reference_genome"], "toy")

    def test_wrong_assembly_missing_mapq_and_pre_dedup_provenance_fail(self):
        with gzip.open(self.source, "rt") as handle:
            original = handle.read()
        for content, message in ((original.replace("assembly: test", "assembly: wrong"), "does not match"),
                                 (original.replace(" mapq1", " absent_mapq"), "required mapping fields"),
                                 (original.replace(f"scaffoldA {len(self.fixture.sequence)}", "scaffoldA 1"), "disagree")):
            with self.subTest(message=message):
                with gzip.open(self.source, "wt") as handle: handle.write(content)
                with self.assertRaisesRegex(ValueError, message): self.plan()
        with gzip.open(self.source, "wt") as handle: handle.write(original)
        json_write(proof_path(self.source), {"stage": "sorted", "pre_dedup": True})
        with self.assertRaisesRegex(ValueError, "pre-dedup"):
            self.plan()
        json_write(proof_path(self.source), {"stage": "dedup", "members": ["technical_A", "technical_B"], "policy": {"min_cis_dist": 1000}})
        self.assertIsNone(self.plan()["samples"]["existing_library"]["input_header"]["previous_filter_policy"])
        self.rows[0]["input_members"] = "other_member"
        self.write_sheet()
        with self.assertRaisesRegex(ValueError, "contradicts known upstream"):
            self.plan()

    def test_dedup_groups_and_duplicate_physical_inputs_are_rejected(self):
        self.rows[0]["dedup_group"] = "do_not_repeat"
        self.write_sheet()
        with self.assertRaisesRegex(ValueError, "never runs deduplication"):
            self.plan()
        self.rows[0]["dedup_group"] = ""
        self.rows.append({**self.rows[0], "sample": "alias"})
        self.write_sheet()
        with self.assertRaisesRegex(ValueError, "Duplicate physical"):
            self.plan()

    def test_inputs_must_be_outside_new_output_roots(self):
        before = sha256(self.source)
        with self.assertRaisesRegex(ValueError, "outside the new output"):
            self.plan(OUTDIR=str(self.root))
        self.assertEqual(sha256(self.source), before)

    def test_pooling_rejects_overlapping_declared_members(self):
        other = self.root / "other.pairs.gz"
        self.emit_source(other, self.records)
        self.rows.append({**self.rows[0], "sample": "other", "deduplicated_pairs": str(other), "input_members": "technical_B;technical_C"})
        self.write_sheet()
        self.plan()  # separate, unpooled comparison products are allowed
        with self.assertRaisesRegex(ValueError, "overlapping input_members"):
            self.plan(RUN_MERGE_ALL_SAMPLES="true")

    def test_normalizer_preserves_multiplicity_swaps_ends_and_removes_stale_fragments(self):
        path = self.root / "annotated.pairs"
        columns = fixtures.COLS + ["rfrag1", "rfrag2", "custom1", "custom2"]
        row = fixtures.contact("reverse", 20, 1, s1="-", s2="+", kind="UR", q1=40, q2=60) + [3, 2, "right", "left"]
        path.write_text("#genome_assembly: test\n#columns: " + " ".join(columns) + "\n" + ("\t".join(map(str,row)) + "\n")*2)
        out = io.StringIO()
        qc = normalize(path, out, read_sizes(self.fixture.sizes, "reference"))
        data = [line.split("\t") for line in out.getvalue().splitlines() if not line.startswith("#")]
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0][2:10], ["1", "scaffoldA", "20", "+", "-", "RU", "60", "40"])
        self.assertEqual(data[0][-2:], ["left", "right"])
        self.assertNotIn("rfrag1", out.getvalue())
        self.assertFalse(qc["deduplication_performed"])
        self.assertEqual(qc["end_swapped_records"], 2)

    def test_input_changes_and_invalid_coordinates_fail_validation(self):
        expected = source_identity(self.source)
        self.emit_source(self.source, [fixtures.contact("invalid", 1, 10000)])
        with self.assertRaisesRegex(ValueError, "changed since planning"):
            normalize(self.source, io.StringIO(), read_sizes(self.fixture.sizes, "ref"), expected)
        with self.assertRaisesRegex(ValueError, "outside configured reference"):
            normalize(self.source, io.StringIO(), read_sizes(self.fixture.sizes, "ref"))

    def test_preflight_does_not_query_aligners(self):
        plan = self.plan()
        queried = []
        def which(name):
            queried.append(name)
            return None if name in {"bwa-mem2", "fastp", "cooltools", "multiqc"} else name
        def version(args, **kwargs):
            return SimpleNamespace(stdout="pairtools, version 1.1.3" if args[0]=="pairtools" else "mock version", stderr="")
        with patch("pipeline.shutil.which", side_effect=which), patch("pipeline.subprocess.run", side_effect=version):
            tool_versions(plan)
        self.assertNotIn("bwa-mem2", queried)
        self.assertNotIn("fastp", queried)

    def test_signature_ignores_unused_alignment_and_dedup_settings(self):
        def keys(plan):
            rt = Runtime(plan, self.root / "config", versions={"pairtools": "1.1.3"})
            return PairProcessor(rt).keys
        original = keys(self.plan())
        self.assertEqual(original, keys(self.plan(DEDUP_MAX_MISMATCH="99", BWA_EXTRA_ARGS="-different", RUN_FASTP="false", PAIRTOOLS_PARSE_MIN_MAPQ="10")))
        self.assertNotEqual(original, keys(self.plan(HIC_MIN_CIS_DIST="11")))

    def fake_pipeline(self, calls):
        def fake(rt, commands, output, log_path, timeout=None):
            calls.extend([list(command) for command in commands])
            for command in commands:
                if command[0] in {"bwa-mem2", "fastp", "samtools"} or command[:2] in (["pairtools", "dedup"], ["pairtools", "parse"]):
                    raise AssertionError(f"Forbidden stage in refilter mode: {command}")
            first = commands[0]
            if len(first)>1 and str(first[1]).endswith("refilter_inputs.py"):
                def get(args, flag): return Path(args[args.index(flag)+1])
                normalized = io.StringIO()
                import_qc = normalize(get(first,"--input"), normalized, read_sizes(get(first,"--chrom-sizes"),"ref"), json.loads(get(first,"--expected-identity").read_text()))
                json_write(get(first,"--qc"), import_qc)
                filtering = commands[1]
                text = io.StringIO()
                with patch("sys.stdin", io.StringIO(normalized.getvalue())):
                    qc = filter_contacts(Path("-"), text, json.loads(get(filtering,"--policy").read_text()),
                        read_sizes(get(filtering,"--chrom-sizes"),"ref"), set(read_sizes(get(filtering,"--matrix-chrom-sizes"),"retained")),
                        load_digest(get(filtering,"--digest")) if "--digest" in filtering else None)
                json_write(get(filtering,"--qc"), qc)
                lines = text.getvalue().splitlines()
                data = [line for line in lines if not line.startswith("#")]
                data.sort(key=lambda line:(line.split("\t")[1],line.split("\t")[3],int(line.split("\t")[2]),int(line.split("\t")[4])))
                with gzip.open(output,"wt") as handle:
                    handle.write("\n".join([line for line in lines if line.startswith("#")] + data) + "\n")
            elif first[0] == "pairix":
                Path(first[-1] + ".px2").write_bytes(b"new mock index")
            elif first[:2] == ["pairtools", "merge"]:
                sources = [Path(value) for value in first[2:] if Path(value).is_file()]
                header, data = [], []
                for index, source in enumerate(sources):
                    with text_open(source) as handle:
                        for line in handle:
                            if line.startswith("#"):
                                if index == 0: header.append(line.rstrip())
                            elif line.strip(): data.append(line.rstrip())
                data.sort(key=lambda line:(line.split("\t")[1],line.split("\t")[3],int(line.split("\t")[2]),int(line.split("\t")[4])))
                with gzip.open(output,"wt") as handle:
                    handle.write("\n".join(header+data)+"\n")
            elif first[:2] == ["pairtools", "stats"]:
                with text_open(Path(first[-1])) as handle: count=sum(not line.startswith("#") for line in handle)
                Path(output).write_text(f"total\t{count}\n")
            elif first[:3] == ["cooler", "cload", "pairix"]:
                with text_open(Path(first[-2])) as handle: count=sum(not line.startswith("#") for line in handle)
                Path(first[-1]).write_text(str(count))
            elif first[:2] == ["cooler", "zoomify"]:
                shutil.copyfile(first[-1], first[first.index("-o")+1])
            elif first[:2] == ["cooler", "info"]:
                Path(output).write_text(json.dumps({"sum": int(Path(first[-1]).read_text())}))
            elif first[0] == "java" and "pre" in first:
                with gzip.open(first[-3],"rt") as handle: data=[line.split() for line in handle]
                self.assertTrue(all(row[3]!=row[7] for row in data))
                Path(first[-2]).write_bytes(b"HIC\0"+str(len(data)).encode())
            elif first[:2] == ["cooler", "ls"] or (first[0] == "java" and "addNorm" in first):
                pass
            else:
                raise AssertionError(first)
        return fake

    def test_runner_skips_alignment_and_dedup_preserves_inputs_and_rebuilds_products(self):
        original, index, original_id = sha256(self.source), sha256(self.index), source_identity(self.source)
        calls = []
        with patch("pipeline.tool_versions", return_value={"pairtools":"1.1.3","cooler":"mock","java":"mock"}), patch.object(Runtime,"pipeline",self.fake_pipeline(calls)), patch("pipeline.log"):
            plan = self.plan()
            self.assertEqual(execute(plan,self.root/"config"),0)
            paths=plan["entities"]["sample:existing_library"]["paths"]
            qc=json.loads(Path(paths["filter_qc"]).read_text())
            self.assertEqual(qc["exclusive_counts"]["accepted"],5)
            self.assertEqual(int(Path(paths["norm_mcool"]).read_text()),5)
            self.assertEqual(Path(paths["norm_hic"]).read_bytes(),b"HIC\x005")
            proof=json.loads(proof_path(Path(paths["accepted"])).read_text())
            self.assertFalse(proof["deduplication_performed"])
            self.assertIsNone(proof["dedup_settings"])
            self.assertEqual(len(proof["import_sources"]),1)
            self.assertFalse((Path(paths["done"])/"dedup.json").exists())
            self.assertFalse((Path(paths["done"])/"sorted.json").exists())
            validate_pool_inputs([Path(paths["accepted"])])
            calls.clear()
            self.assertEqual(execute(self.plan(),self.root/"config"),0)
            self.assertEqual(calls,[])
            self.assertEqual(execute(self.plan(HIC_MIN_CIS_DIST="0"),self.root/"config"),0)
            self.assertEqual(json.loads(Path(paths["filter_qc"]).read_text())["exclusive_counts"]["accepted"],6)
            other = self.root / "independent_library.pairs.gz"
            self.emit_source(other, [["other_"+row[0],*row[1:]] for row in self.records])
            self.rows.append({**self.rows[0], "sample": "independent_library", "deduplicated_pairs": str(other), "input_members": "library_C"})
            self.write_sheet()
            pooled = self.plan(HIC_MIN_CIS_DIST="0",RUN_MERGE_ALL_SAMPLES="true")
            self.assertEqual(execute(pooled,self.root/"config"),0)
            pool_paths = pooled["entities"]["pool:all_samples_merged"]["paths"]
            self.assertEqual(json.loads(Path(pool_paths["filter_qc"]).read_text())["exclusive_counts"]["accepted"],12)
            self.assertEqual(len(json.loads(proof_path(Path(pool_paths["accepted"])).read_text())["import_sources"]),2)
        self.assertEqual(sha256(self.source),original)
        self.assertEqual(sha256(self.index),index)
        self.assertEqual(source_identity(self.source),original_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
