"""Compact scientific, metadata, restart and storage regressions; no production data."""
from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from artifact_state import FileLock, StageStore, StorageGuard, publish
from filter_contact_pairs import classify, filter_contacts
from matrix_products import export_juicer
from pairs_processing import PairProcessor, proof_path, validate_pool_inputs
from pipeline import Runtime, execute, manual_plan, validate_existing_layout
from prepare_restriction_fragments import load_digest, prepare, scan_fasta
from reference_registry import BWA_MEM2_INDEX_SUFFIXES, REGISTRY_COLUMNS
from run_plan import build_plan, resolve_roots
from summarize_run import write_plan_summary
from workflow_common import DEFAULTS, identity, json_write, load_config, registry_data, resolve_policy, signature

COLS = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2".split()


def write_pairs(path, records):
    with path.open("w") as out:
        out.write("## pairs format v1.0\n#genome_assembly: test\n#columns: " + " ".join(COLS) + "\n")
        for row in records:
            out.write("\t".join(map(str, row)) + "\n")


def contact(name, p1=10, p2=20, c1="scaffoldA", c2="scaffoldA", s1="+", s2="-", kind="UU", q1=30, q2=30):
    return [name, c1, p1, c2, p2, s1, s2, kind, q1, q2]


class IntegrationFeatures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.c = {**DEFAULTS, "OUTDIR": str(self.root / "out"), "DEFAULT_REFERENCE_ID": "toy",
                  "MIN_FREE_OUTPUT_GB": "0", "MIN_FREE_SCRATCH_GB": "0", "RUN_PRELIM_DOWNSTREAM": "false",
                  "HIC_PROTOCOL": "restriction_hic", "HIC_RESTRICTION_ENZYMES": "DpnII", "MICROC_MIN_CIS_DIST": "10"}
        self.fasta = self.root / "genome.fa"
        self.sequence = "ttGATCaaCATGaaAAGCTTnn" + "A"*100
        self.fasta.write_text(f">scaffoldA\n{self.sequence}\n>plasmidZ\n{'A'*len(self.sequence)}\n")
        self.sizes = self.root / "sizes.tsv"
        self.sizes.write_text(f"scaffoldA\t{len(self.sequence)}\nplasmidZ\t{len(self.sequence)}\n")
        Path(str(self.fasta) + ".fai").write_text(self.sizes.read_text())
        for suffix in BWA_MEM2_INDEX_SUFFIXES:
            Path(str(self.fasta) + suffix).write_text("fixture")
        self.registry = self.root / "references.tsv"
        with self.registry.open("w", newline="") as out:
            writer = csv.DictWriter(out, REGISTRY_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"reference_id": "toy", "aliases": "", "species": "Synthetic", "assembly": "test", "fasta": str(self.fasta),
                            "bwa_index_prefix": str(self.fasta), "chrom_sizes": str(self.sizes), "canonical_regex": ".*", "browser_preset": "custom"})
        self.c["REFERENCE_REGISTRY"] = str(self.registry)
        self.fastq = self.root / "R.fastq"
        self.fastq.write_text("@r\nACGT\n+\nIIII\n")
        self.sheet = self.root / "samples.csv"
        self.rows = [{"sample": f"row{i}", "assay": "microc", "reference_genome": "toy", "condition": "arbitrary state",
                      "biological_replicate": str((i+1)//2), "technical_replicate": str((i-1)%2+1),
                      "fastq_r1": str(self.fastq), "fastq_r2": str(self.fastq), "merge_group": "",
                      "hic_protocol": "", "restriction_enzymes": "", "dedup_group": "shared" if i<=2 else ""} for i in range(1, 5)]
        self.write_sheet()
        self.registry_data = registry_data(ROOT / "config/hic_protocols.json")

    def tearDown(self):
        self.temp.cleanup()

    def write_sheet(self):
        with self.sheet.open("w", newline="") as out:
            writer = csv.DictWriter(out, list(self.rows[0]), lineterminator="\r\n")
            writer.writeheader()
            writer.writerows(self.rows)

    def policy(self, assay="microc"):
        return resolve_policy({"assay": assay}, self.c)

    def filter(self, rows, policy=None, digest=None):
        source = self.root / "pairs.txt"
        write_pairs(source, rows)
        out = io.StringIO()
        qc = filter_contacts(source, out, policy or self.policy(), {"scaffoldA": len(self.sequence), "plasmidZ": len(self.sequence)}, {"scaffoldA", "plasmidZ"}, digest)
        return [line.split("\t") for line in out.getvalue().splitlines() if not line.startswith("#")], qc

    def test_all_enzyme_cuts_chunk_edges_lowercase_and_combinations(self):
        cuts = {}
        for name, offset in (("DpnII", 0), ("MboI", 0), ("NlaIII", 4), ("HindIII", 1)):
            definition = self.registry_data["enzymes"][name]
            actual = list(scan_fasta(self.fasta, {name: definition}, chunk_size=3))[0][2]
            expected = self.sequence.upper().index(definition["motif"]) + offset
            self.assertEqual(actual, [0, expected, len(self.sequence)])
            cuts[name] = actual
        self.assertEqual(cuts["DpnII"], cuts["MboI"])
        all_cuts = list(scan_fasta(self.fasta, self.registry_data["enzymes"], chunk_size=4))[0][2]
        self.assertEqual(all_cuts, sorted({c for values in cuts.values() for c in values}))

    def test_digest_is_complete_and_cache_reuses_verified_content(self):
        defs = {"DpnII": self.registry_data["enzymes"]["DpnII"]}
        path, meta = prepare(self.fasta, self.sizes, defs, self.root / "cache")
        old = path.stat().st_mtime_ns
        self.assertEqual(prepare(self.fasta, self.sizes, defs, self.root / "cache")[1], meta)
        self.assertEqual(path.stat().st_mtime_ns, old)
        self.assertEqual(list(load_digest(path)["scaffoldA"]), [0, 2, len(self.sequence)])
        with self.assertRaisesRegex(ValueError, "inconsistency"):
            bad = self.root / "bad.sizes"
            bad.write_text("scaffoldA\t1\nplasmidZ\t1\n")
            prepare(self.fasta, bad, defs, self.root / "cache")

    def test_fragment_orientation_classes_and_trans_ordinal_collision(self):
        args = ("scaffoldA", 1, "+", 0, "scaffoldA", 2, "-", 0)
        self.assertEqual(classify(*args), "dangling_end_candidate")
        self.assertEqual(classify("scaffoldA", 1, "-", 0, "scaffoldA", 2, "+", 0), "self_circle_candidate")
        self.assertEqual(classify("scaffoldA", 1, "+", 0, "scaffoldA", 1, "-", 0), "same_fragment_other")
        self.assertEqual(classify("scaffoldA", 1, "+", 0, "scaffoldA", 3, "-", 1), "adjacent_inward_candidate")
        self.assertEqual(classify("scaffoldA", 1, "+", 0, "plasmidZ", 2, "-", 0), "trans")
        self.assertEqual(classify("scaffoldA", 2, "-", 0, "scaffoldA", 1, "+", 0), "dangling_end_candidate")

    def test_cut_coordinates_use_pairs_pos_minus_one(self):
        policy = self.policy("hic")
        policy["remove_same_fragment"] = False
        bed, _ = prepare(self.fasta, self.sizes, policy["enzyme_definitions"], self.root / "cache")
        rows, qc = self.filter([contact("left_right", 2, 3)], policy, load_digest(bed))
        self.assertEqual(rows[0][-2:], ["0", "1"])
        self.assertEqual(qc["exclusive_counts"]["accepted"], 1)

    def test_restriction_defaults_and_opt_in_adjacent_rejection(self):
        p = self.policy("hic")
        bed, _ = prepare(self.fasta, self.sizes, p["enzyme_definitions"], self.root / "cache")
        digest = load_digest(bed)
        records = [contact("dangling", 3, 4), contact("circle", 3, 4, s1="-", s2="+"),
                   contact("other", 3, 4, s1="+", s2="+"), contact("adjacent", 1, 3), contact("trans", 1, 3, c2="plasmidZ")]
        accepted, qc = self.filter(records, p, digest)
        self.assertEqual([r[0] for r in accepted], ["adjacent", "trans"])
        self.assertEqual(sum(v for k,v in qc["exclusive_counts"].items() if k!="input"), 5)
        p["remove_adjacent_inward"] = True
        self.assertEqual([r[0] for r in self.filter(records, p, digest)[0]], ["trans"])

    def test_distance_zero_below_equal_above_and_trans(self):
        records = [contact("below", 10, 19), contact("equal", 10, 20), contact("above", 10, 21), contact("trans", 10, 11, c2="plasmidZ")]
        accepted, _ = self.filter(records)
        self.assertEqual([r[0] for r in accepted], ["equal", "above", "trans"])
        policy = self.policy()
        policy["min_cis_dist"] = 0
        self.assertEqual(len(self.filter(records, policy)[0]), 4)

    def test_microc_contacts_match_legacy_predicate(self):
        records = [contact("good"), contact("close", 10, 11), contact("low", q1=29), contact("wrong", kind="DD"),
                   contact("UR", kind="UR"), contact("trans", 10, 11, c2="plasmidZ")]
        expected = [r[0] for r in records if r[7] in {"UU", "UR", "RU"} and min(r[8:])>=30 and (r[1]!=r[3] or abs(r[2]-r[4])>=10)]
        self.assertEqual([r[0] for r in self.filter(records)[0]], expected)

    def test_mapped_outside_reference_fails_not_silent_acceptance(self):
        with self.assertRaisesRegex(ValueError, "outside reference"):
            self.filter([contact("bad", 0, 20)])

    def test_profiles_generic_dovetail_legacy_and_explicit_enzyme_requirement(self):
        self.assertEqual(self.policy("hic")["enzymes"], ["DpnII"])
        c = {**self.c, "HIC_PROTOCOL": "dovetail_hic", "HIC_RESTRICTION_ENZYMES": ""}
        self.assertEqual(resolve_policy({"assay": "hic"}, c)["enzymes"], ["DpnII"])
        with self.assertRaisesRegex(ValueError, "requires explicit"):
            resolve_policy({"assay": "hic"}, {**self.c, "HIC_RESTRICTION_ENZYMES": ""})
        with self.assertRaisesRegex(ValueError, "original DpnII"):
            resolve_policy({"assay": "hic", "hic_protocol": "dovetail_hic", "restriction_enzymes": "HindIII"}, self.c)
        self.assertEqual(resolve_policy({"assay": "hic"}, {**self.c, "HIC_PROTOCOL": "legacy", "HIC_RESTRICTION_ENZYMES": ""})["enzymes"], [])

    def test_empty_roots_preserve_structure_and_override_only_the_category(self):
        plan = build_plan(self.c, self.sheet)
        e = plan["entities"]["sample:row1"]
        self.assertEqual(e["paths"]["accepted"], str(self.root / "out/row1/03_pairs/row1.valid.mapq30.pairs.gz"))
        moved = build_plan({**self.c, "PAIRS_OUTDIR": str(self.root / "hdd")}, self.sheet)
        self.assertEqual(moved["entities"]["sample:row1"]["paths"]["accepted"], str(self.root / "hdd/row1/03_pairs/row1.valid.mapq30.pairs.gz"))
        self.assertEqual(moved["entities"]["sample:row1"]["paths"]["matrix_dir"], e["paths"]["matrix_dir"])
        self.assertEqual(resolve_roots(self.c)["scratch"], str(self.root / "out/tmp"))

    def test_technical_and_all_pools_use_disjoint_units_once(self):
        plan = build_plan({**self.c, "RUN_MERGE_ALL_SAMPLES": "true"}, self.sheet)
        all_pool = plan["entities"]["pool:all_samples_merged"]
        self.assertEqual(all_pool["units"], ["dedup:shared", "sample:row3", "sample:row4"])
        first = plan["entities"]["pool:arbitrary_state_microc_B1_tech_merged"]
        self.assertEqual(first["units"], ["dedup:shared"])
        self.assertEqual(len(plan["samples"]), 4)

    def test_partial_group_and_unsafe_identifiers_rejected(self):
        self.rows[1]["biological_replicate"] = "2"
        self.rows[1]["technical_replicate"] = "3"
        self.write_sheet()
        with self.assertRaisesRegex(ValueError, "only part"):
            build_plan(self.c, self.sheet)
        self.rows[0]["dedup_group"] = ".."
        self.write_sheet()
        with self.assertRaisesRegex(ValueError, "unsafe dedup"):
            build_plan(self.c, self.sheet)

    def test_restart_invalidates_selection_not_alignment_and_storage_is_not_scientific(self):
        def keys(c):
            rt = Runtime(build_plan(c, self.sheet), self.root / "config", versions={"pairtools": "1.1.3", "bwa-mem2": "2.2.1"})
            return PairProcessor(rt).keys
        original = keys(self.c)
        changed = keys({**self.c, "MICROC_MIN_CIS_DIST": "11"})
        self.assertEqual(original["sample:row1"]["sorted"], changed["sample:row1"]["sorted"])
        self.assertEqual(original["sample:row1"]["dedup"], changed["sample:row1"]["dedup"])
        self.assertNotEqual(original["sample:row1"]["accepted"], changed["sample:row1"]["accepted"])
        moved = keys({**self.c, "PAIRS_OUTDIR": str(self.root / "other")})
        self.assertEqual(original, moved)
        changed_dedup = keys({**self.c, "DEDUP_MAX_MISMATCH": "4"})
        self.assertEqual(original["sample:row1"]["sorted"], changed_dedup["sample:row1"]["sorted"])
        self.assertNotEqual(original["sample:row1"]["dedup"], changed_dedup["sample:row1"]["dedup"])
        self.assertEqual(original, keys({**self.c, "MCOOL_RESOLUTIONS": "2000,10000,100000"}))
        for row in self.rows:
            row["assay"] = "hic"
        self.write_sheet()
        dpnii = keys(self.c)
        hindiii = keys({**self.c, "HIC_RESTRICTION_ENZYMES": "HindIII"})
        self.assertEqual(dpnii["sample:row1"]["sorted"], hindiii["sample:row1"]["sorted"])
        self.assertEqual(dpnii["sample:row1"]["dedup"], hindiii["sample:row1"]["dedup"])
        self.assertNotEqual(dpnii["sample:row1"]["accepted"], hindiii["sample:row1"]["accepted"])

    def store(self):
        roots = resolve_roots(self.c)
        guard = StorageGuard(roots, [], 0, 0)
        return StageStore(self.root / "state", self.root / "scratch/run_unique", [self.root / "outputs"], "destination", guard)

    def test_stage_invalidation_unowned_outputs_and_narrow_cleanup(self):
        store = self.store()
        target = self.root / "outputs/result.txt"
        action = lambda t,w: t["data"].write_text("result")
        store.run("select", "one", {"data": target}, action)
        self.assertTrue(store.complete("select", "one", {"data": target}))
        self.assertFalse(store.complete("select", "two", {"data": target}))
        original = self.root / "outputs/original.fastq"
        original.write_text("do not remove")
        self.assertFalse(store.remove_owned(original))
        self.assertTrue(original.exists())
        with self.assertRaisesRegex(RuntimeError, "unowned"):
            store.run("foreign", "key", {"data": original}, action)
        self.assertTrue(store.remove_owned(target))
        with self.assertRaisesRegex(RuntimeError, "outside"):
            store.remove_owned(self.fastq)

    def test_failed_bundle_not_complete_and_retries_own_published_part(self):
        store = self.store()
        outputs = {"a": self.root / "outputs/a", "b": self.root / "outputs/b"}
        def action(t,w):
            for p in t.values(): p.write_text("data")
        calls = 0
        def fail_second(source, destination):
            nonlocal calls
            calls += 1
            if calls == 2: raise OSError("interrupted transfer")
            return publish(source, destination)
        with patch("artifact_state.publish", side_effect=fail_second):
            with self.assertRaisesRegex(OSError, "interrupted"):
                store.run("bundle", "key", outputs, action)
        self.assertFalse(store.complete("bundle", "key", outputs))
        store.run("bundle", "key", outputs, action)
        self.assertTrue(store.complete("bundle", "key", outputs))

    def test_mount_reserve_and_lock_protections(self):
        with self.assertRaisesRegex(RuntimeError, "mount is absent"):
            StorageGuard(resolve_roots(self.c), [str(self.root / "absent_mount")], 0, 0).check()
        with self.assertRaisesRegex(RuntimeError, "Insufficient"):
            StorageGuard(resolve_roots(self.c), [], 10**10, 0).check()
        path = self.root / "lock"
        with FileLock(path):
            with self.assertRaisesRegex(RuntimeError, "already locked"):
                with FileLock(path): pass
        with FileLock(path): pass

    def test_cross_disk_transfer_checks_checksum_and_keeps_source_on_failure(self):
        source, dest = self.root / "source.bin", self.root / "published/result.bin"
        source.write_bytes(b"original bytes")
        with patch("artifact_state.same_filesystem", return_value=False), patch("artifact_state.sha256", side_effect=["correct", "wrong"]):
            with self.assertRaisesRegex(RuntimeError, "checksum mismatch"):
                publish(source, dest)
        self.assertTrue(source.exists())
        self.assertFalse(dest.exists())
        with patch("artifact_state.same_filesystem", return_value=False):
            publish(source, dest)
        self.assertEqual(dest.read_bytes(), b"original bytes")

    def test_managed_storage_relocation_verifies_copy_without_recomputation(self):
        store = self.store()
        old, new = self.root / "outputs/old/data.bin", self.root / "outputs/new/data.bin"
        store.run("data", "unchanged-science", {"data": old}, lambda t,w: t["data"].write_bytes(b"preserve"))
        def must_not_run(t,w):
            raise AssertionError("storage changes must not rerun computation")
        record = store.run("data", "unchanged-science", {"data": new}, must_not_run)
        self.assertTrue(record["storage_only_republication"])
        self.assertEqual(new.read_bytes(), old.read_bytes())
        self.assertTrue(old.exists())

    def test_juicer_projection_preserves_contacts_and_matches_bin_boundaries(self):
        source = self.root / "accepted.pairs"
        records = [contact("edge", 10, 20), contact("one", 1, 11, c2="plasmidZ")]
        write_pairs(source, records)
        output = self.root / "juicer.gz"
        self.assertEqual(export_juicer(source, output), 2)
        with gzip.open(output, "rt") as handle:
            rows = [line.split() for line in handle]
        self.assertEqual(rows[0], ["0", "scaffoldA", "9", "0", "16", "scaffoldA", "19", "1"])
        self.assertEqual([int(row[2])//10 for row in rows], [(r[2]-1)//10 for r in records])

    def test_report_uses_resolved_paths_and_labels_scopes(self):
        plan = build_plan({**self.c, "MATRIX_OUTDIR": str(self.root / "matrix_hdd")}, self.sheet)
        report = self.root / "report"
        write_plan_summary(plan, report)
        text = (report / "microc2tracks_summary.html").read_text()
        self.assertIn("matrix_hdd", text)
        self.assertIn("need not sum", text)
        self.assertIn("dedup:shared", text)

    def test_legacy_metadata_is_rejected_before_status_or_manifest_overwrite(self):
        plan = build_plan(self.c, self.sheet)
        manifest = self.root / "out/run_metadata/sample_manifest.tsv"
        manifest.parent.mkdir(parents=True)
        manifest.write_text("preserve original metadata")
        status = Path(plan["entities"]["sample:row1"]["paths"]["status"])
        status.parent.mkdir(parents=True)
        status.write_text("preserve original status")
        with self.assertRaisesRegex(RuntimeError, "legacy run metadata"):
            validate_existing_layout(plan)
        self.assertEqual(manifest.read_text(), "preserve original metadata")
        self.assertEqual(status.read_text(), "preserve original status")

    def test_crlf_config_with_spaces_is_read_without_modifying_original(self):
        bash = Path("C:/Program Files/Git/bin/bash.exe") if os.name == "nt" else Path("/bin/bash")
        if not bash.is_file():
            self.skipTest("Bash not installed")
        path = self.root / "config.conf"
        value = str(self.root / "folder with spaces")
        original = b"\xef\xbb\xbf" + (f'OUTDIR="{value}"\r\nTMPDIR=""\r\n').encode()
        path.write_bytes(original)
        with patch.dict(os.environ, {"MICROC2TRACKS_BASH": str(bash)}):
            config = load_config(path)
        self.assertEqual(config["OUTDIR"], value)
        self.assertEqual(config["HIC_PROTOCOL"], "legacy")
        self.assertEqual(path.read_bytes(), original)

    def test_custom_browser_contigs_are_not_forced_to_chr_names(self):
        from package_browser_tracks import canonical_chroms, read_track_rows, SampleSpec, chrom_sizes_candidates
        table = self.root / "track.tsv"
        table.write_text("chrom\tstart\tend\tE1\nscaffoldA\t0\t10\t0.5\n")
        result = read_track_rows(table, "E1", {"scaffoldA", "plasmidZ"})
        self.assertEqual(result[0].chrom, "scaffoldA")
        self.assertEqual(canonical_chroms("custom", "scaffoldA,plasmidZ"), ["scaffoldA", "plasmidZ"])
        self.assertEqual(chrom_sizes_candidates(SampleSpec("x", self.root, self.sizes)), [self.sizes])

    def test_runner_contract_mocked_tools_pooling_restart_and_storage_relocation(self):
        """Orchestration contract only: mocks are NOT scientific/tool validation."""
        calls = []
        def rows(path):
            with text_open(Path(path)) as handle:
                return [line.rstrip().split("\t") for line in handle if not line.startswith("#") and line.strip()]
        def emit(path, data):
            header = "## pairs format v1.0\n#genome_assembly: test\n#columns: " + " ".join(COLS) + "\n"
            with gzip.open(path, "wt") as out:
                out.write(header)
                for row in data: out.write("\t".join(map(str,row)) + "\n")
        def fake(rt, commands, output, log_path, timeout=None):
            args = commands[0]
            calls.append(args[:2])
            if args[0] == "bwa-mem2":
                index = int(log_path.name.split(".")[0][3:])
                emit(output, [contact("common", 10, 30), contact(f"unique{index}", 30+index*10, 50+index*10)])
            elif args[:2] == ["pairtools", "merge"]:
                # Tool flags are fixed by production code; the last values are inputs.
                inputs = [p for p in args[2:] if Path(str(p)).is_file()]
                merged = [row for path in inputs for row in rows(path)]
                merged.sort(key=lambda r: (r[1],r[3],int(r[2]),int(r[4])))
                emit(output, merged)
            elif args[:2] == ["pairtools", "dedup"]:
                data, seen, result = rows(args[-1]), set(), []
                for row in data:
                    key = tuple(row[1:7])
                    if key not in seen: result.append(row); seen.add(key)
                emit(output, result)
                Path(args[args.index("--output-stats")+1]).write_text(f"total\t{len(data)}\ntotal_nodups\t{len(result)}\ntotal_dups\t{len(data)-len(result)}\n")
            elif len(args)>1 and str(args[1]).endswith("filter_contact_pairs.py"):
                get = lambda flag: Path(args[args.index(flag)+1])
                from reference_registry import read_sizes
                with gzip.open(output, "wt") as out:
                    qc = filter_contacts(get("--input"), out, json.loads(get("--policy").read_text()),
                        read_sizes(get("--chrom-sizes"), "ref"), set(read_sizes(get("--matrix-chrom-sizes"), "retained")))
                json_write(get("--qc"), qc)
            elif args[0] == "pairix":
                Path(args[-1] + ".px2").write_bytes(b"fake-index")
            elif args[:2] == ["pairtools", "stats"]:
                Path(output).write_text(f"total\t{len(rows(args[-1]))}\n")
            elif args[:3] == ["cooler", "cload", "pairix"]:
                Path(args[-1]).write_text(str(len(rows(args[-2]))))
            elif args[:2] == ["cooler", "zoomify"]:
                shutil.copyfile(args[-1], args[args.index("-o")+1])
            elif args[:2] == ["cooler", "info"]:
                Path(output).write_text(json.dumps({"sum": int(Path(args[-1]).read_text())}))
            elif args[:2] != ["cooler", "ls"]:
                raise AssertionError(args)
        import shutil
        from workflow_common import text_open
        c = {**self.c, "RUN_FASTP": "false", "RUN_HIC": "false", "RUN_MULTIQC": "false", "RUN_GLOBAL_MULTIQC": "false", "RUN_MERGE_ALL_SAMPLES": "true"}
        with patch("pipeline.tool_versions", return_value={"pairtools": "1.1.3", "cooler": "fake"}), patch.object(Runtime, "pipeline", fake), patch("pipeline.log"):
            plan = build_plan(c, self.sheet)
            self.assertEqual(execute(plan, self.root / "config"), 0)
            pool = plan["entities"]["pool:all_samples_merged"]["paths"]
            self.assertEqual(len(rows(pool["accepted"])), 7)
            self.assertEqual(int(Path(pool["norm_mcool"]).read_text()), 7)
            calls.clear()
            self.assertEqual(execute(build_plan(c, self.sheet), self.root / "config"), 0)
            self.assertEqual(calls, [])
            moved = build_plan({**c, "PAIRS_OUTDIR": str(self.root / "new_pairs"), "MATRIX_OUTDIR": str(self.root / "new_matrices")}, self.sheet)
            self.assertEqual(execute(moved, self.root / "config"), 0)
            self.assertEqual(calls, [])
            self.assertTrue(Path(pool["accepted"]).is_file())  # old copies are preserved
            self.assertEqual(len(rows(moved["entities"]["pool:all_samples_merged"]["paths"]["accepted"])), 7)
            import argparse
            original_plan = self.root / "out/run_metadata/run_plan.json"
            original_bytes = original_plan.read_bytes()
            inputs = [moved["entities"][key]["paths"]["accepted"] for key in ("dedup:shared", "sample:row3")]
            args = argparse.Namespace(name="manual_pool", pairs=",".join(inputs), mode="pool", reference="toy")
            manual = manual_plan(c, args)
            self.assertEqual(execute(manual, self.root / "config"), 0)
            self.assertEqual(len(rows(manual["entities"]["pool:manual_pool"]["paths"]["accepted"])), 5)
            self.assertEqual(original_plan.read_bytes(), original_bytes)
            checked_source = Path(moved["entities"]["dedup:shared"]["paths"]["accepted"])
            proof = json.loads(proof_path(checked_source).read_text())
            marker = Path(proof["completion_record"])
            saved = json.loads(marker.read_text())
            json_write(marker, {**saved, "status": "failed"})
            with self.assertRaisesRegex(ValueError, "incomplete, failed, or stale"):
                validate_pool_inputs([checked_source])
            json_write(marker, saved)
            inputs = [moved["entities"][f"sample:row{i}"]["paths"]["sorted"] for i in (1,2)]
            # Sorted originals can be regenerated/reused independently of accepted outputs.
            # The storage move reused accepted bundles only, leaving original sorted paths.
            inputs = [plan["entities"][f"sample:row{i}"]["paths"]["sorted"] for i in (1,2)]
            args = argparse.Namespace(name="manual_joint", pairs=",".join(inputs), mode="joint-dedup", reference="toy")
            manual = manual_plan(c, args)
            self.assertEqual(execute(manual, self.root / "config"), 0)
            self.assertEqual(len(rows(manual["entities"]["dedup:manual_joint"]["paths"]["accepted"])), 3)
            self.assertEqual(original_plan.read_bytes(), original_bytes)


if __name__ == "__main__":
    unittest.main(verbosity=2)
