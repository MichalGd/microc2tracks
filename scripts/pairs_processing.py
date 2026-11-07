"""Production pair processing shared by rows, joint-dedup groups and pools."""
from __future__ import annotations

import json
import shlex
import shutil
import sys
from pathlib import Path

from prepare_restriction_fragments import prepare
from reference_registry import read_sizes
from refilter_inputs import validate_snapshot
from workflow_common import ROOT, identity, json_write, signature, truth


def proof_path(path: Path) -> Path:
    return Path(str(path) + ".provenance.json")


def read_proof(path: Path) -> dict:
    try:
        return json.loads(proof_path(path).read_text())
    except FileNotFoundError as exc:
        raise ValueError(f"Missing stage provenance for {path}. Joint-dedup requires original sorted streams. Refilter existing deduplicated pairs with microc2tracks.sh --mode refilter-deduplicated; only already accepted legacy pairs may use merge --mode legacy-pool.") from exc


class PairProcessor:
    def __init__(self, runtime):
        self.rt, self.plan, self.c = runtime, runtime.plan, runtime.c
        self.keys = {}
        for name, row in self.plan["samples"].items():
            ref = self.plan["references"][row["reference_genome"]]
            if self.plan.get("workflow_mode") == "refilter-deduplicated":
                self.keys[f"sample:{name}"] = {"input": runtime.key("deduplicated_input", {"pairs": row["input_identity"],
                    "upstream_provenance": row["input_header"]["sidecar_identity"]}, [], ["refilter_inputs.py"])}
                continue
            if any(not Path(row[k]).is_file() for k in ("fastq_r1", "fastq_r2")):
                from checkpoint_resume import sorted_checkpoint
                self.keys[f"sample:{name}"] = {"trim": "missing-raw-inputs", "sorted": sorted_checkpoint(runtime, name)}
                continue
            parameters = {"inputs": [identity(Path(row[k])) for k in ("fastq_r1", "fastq_r2")],
                          "run_fastp": self.c["RUN_FASTP"], "fastp_args": self.c["FASTP_EXTRA_ARGS"]}
            if row.get("preprocessing"):
                from blhic_preprocessing import scientific_tool_identity
                parameters["read_preprocessing"] = row["preprocessing"]
                if row["preprocessing"]["input_state"] == "raw":
                    parameters["trimLinker"] = scientific_tool_identity(runtime.versions.get("trimLinker"))
            trim = runtime.key("trim", parameters, ["fastp"])
            sorted_key = runtime.key("sorted", {"trim": trim, "reference": ref["identity"],
                "bwa_args": self.c["BWA_EXTRA_ARGS"], "parse": {k: self.c[k] for k in self.c if k.startswith("PAIRTOOLS_") and k not in {"PAIRTOOLS_SORT_MEMORY", "PAIRTOOLS_MERGE_MEMORY"}}}, ["bwa-mem2", "pairtools"], ["pairs_processing.py"])
            self.keys[f"sample:{name}"] = {"trim": trim, "sorted": sorted_key}
        for key, entity in self.plan["entities"].items():
            if entity["kind"] in {"external", "external_sorted"}:
                self.keys[key] = {"sorted" if entity["kind"] == "external_sorted" else "accepted": entity["source_signature"]}
                continue
            if entity.get("input_stage") == "deduplicated":
                accepted = runtime.key("refilter_deduplicated", {"input": self.keys[key]["input"], "policy": entity["policy"],
                    "reference": self.plan["references"][entity["reference_id"]]["identity"], "members": entity["members"]},
                    ["pairtools", "pairix"], ["pairs_processing.py", "refilter_inputs.py", "filter_contact_pairs.py", "prepare_restriction_fragments.py"])
                self.keys[key]["accepted"] = accepted
                continue
            if entity["kind"] == "dedup":
                self.keys[key] = {"sorted": runtime.key("joint_sorted", {"members": entity["members"],
                    "sources": [self.keys[f"sample:{name}"]["sorted"] for name in entity["members"]]}, ["pairtools"])}
            if entity["kind"] != "pool":
                dedup = runtime.key("dedup", {"sorted": self.keys[key]["sorted"],
                                             "settings": {k: v for k, v in self.c.items() if k.startswith("DEDUP_")}}, ["pairtools"], ["pairs_processing.py"])
                accepted = runtime.key("accepted", {"dedup": dedup, "policy": entity["policy"],
                    "reference": self.plan["references"][entity["reference_id"]]["identity"],
                    "scope": {"members": entity["members"], "dedup_group": entity["dedup_group"]}}, ["pairtools", "pairix"], ["filter_contact_pairs.py", "prepare_restriction_fragments.py"])
                self.keys[key].update(dedup=dedup, accepted=accepted)
        for key, entity in self.plan["entities"].items():
            if entity["kind"] == "pool":
                self.keys[key] = {"accepted": runtime.key("pool", {"members": entity["members"], "units": entity["units"],
                    "sources": [self.keys[unit]["accepted"] for unit in entity["units"]]}, ["pairtools"], ["pairs_processing.py"])}

    def reusable(self, key, stage):
        entity, keys = self.plan["entities"][key], self.keys[key]
        paths = entity["paths"]
        if stage == "trim":
            row = self.plan["samples"][entity["name"]]
            if self.reusable(key, "sorted") or (not truth(self.c["RUN_FASTP"]) and not row.get("preprocessing")):
                return True
            outputs = self.trim_outputs(entity)
        elif stage in {"sorted", "dedup", "accepted"}:
            outputs = {"pairs": Path(paths[stage]), "provenance": proof_path(Path(paths[stage]))}
            if stage == "dedup":
                outputs["stats"] = Path(paths["dedup_stats"])
            if stage == "accepted":
                outputs.update(index=Path(paths["accepted"] + ".px2"), qc=Path(paths["filter_qc"]))
                if entity.get("input_stage") == "deduplicated":
                    outputs["import_qc"] = Path(paths["import_qc"])
        else:
            return False
        return self.rt.store(entity).complete(stage, keys[stage], outputs)

    def trim_outputs(self, entity):
        paths = entity["paths"]
        result = {"r1": Path(paths["trim_r1"]), "r2": Path(paths["trim_r2"])}
        if truth(self.c["RUN_FASTP"]):
            result.update(html=Path(paths["qc_dir"]) / f"{entity['name']}.fastp.html",
                          json=Path(paths["qc_dir"]) / f"{entity['name']}.fastp.json")
        if self.plan["samples"][entity["name"]].get("preprocessing"):
            result.update(blhic_qc=Path(paths["blhic_qc"]), blhic_stats=Path(paths["blhic_stats"]))
        return result

    def available(self, key, stage):
        """Read-only check for a valid checkpoint, including a storage-only move.

        A move still executes under resource admission; this predicate only
        avoids remapping prerequisites that its verified bundle can replace.
        """
        if self.reusable(key, stage):
            return True
        store = self.rt.store(self.plan["entities"][key])
        record = store.read(stage)
        allowed = {self.keys[key][stage], *self.rt.compatible_keys.get(self.keys[key][stage], [])}
        return (record.get("status") == "success" and record.get("signature") in allowed
                and bool(record.get("outputs")) and all(Path(item["path"]).is_file() and identity(Path(item["path"])) == item
                    for item in record["outputs"].values()))

    def ensure_trim(self, key):
        entity = self.plan["entities"][key]
        row = self.plan["samples"][entity["name"]]
        if any(not Path(row[k]).is_file() for k in ("fastq_r1", "fastq_r2")):
            raise ValueError(f"No reusable sorted/deduplicated checkpoint for {key}; missing original FASTQs are required")
        if (not truth(self.c["RUN_FASTP"]) and not row.get("preprocessing")) or self.reusable(key, "sorted"):
            return Path(row["fastq_r1"]), Path(row["fastq_r2"])
        outputs = self.trim_outputs(entity)
        def trim(t, work):
            def fastp(sources):
                self.rt.command(["fastp", "-i", str(sources[0]), "-I", str(sources[1]), "-o", str(t["r1"]), "-O", str(t["r2"]),
                    "--thread", self.c["THREADS_FASTP"], "--html", str(t["html"]), "--json", str(t["json"]), *shlex.split(self.c["FASTP_EXTRA_ARGS"])],
                    self.rt.log(entity, "fastp"), timeout=int(self.c["FASTP_TIMEOUT_SECONDS"]) or None)
            with self.rt.slots["fastp"]:
                if row.get("preprocessing"):
                    from blhic_preprocessing import prepare
                    prepare(self.rt, entity, row, t, work, fastp)
                else:
                    fastp((row["fastq_r1"], row["fastq_r2"]))
        reserved = f"{entity['id']}/trim" in self.rt.guard.reservations
        try:
            self.rt.store(entity).run("trim", self.keys[key]["trim"], outputs, trim)
        finally:
            if row.get("preprocessing") and not reserved:
                self.rt.guard.release(f"{entity['id']}/trim")
        return outputs["r1"], outputs["r2"]

    def provenance(self, entity, stage, key, data_path):
        result = {"schema": 1, "stage": stage, "signature": key, "reference_id": entity["reference_id"],
                "assembly": self.plan["references"][entity["reference_id"]]["assembly"],
                "reference_identity": self.plan["references"][entity["reference_id"]]["identity"],
                "policy": entity["policy"], "members": entity["members"], "dedup_units": entity.get("provenance_units", entity["units"]),
                "poolable": not entity.get("legacy_pool") and (entity["kind"] != "sample" or not entity["dedup_group"]),
                "pre_dedup": stage == "sorted", "data": str(data_path), "tool_versions": self.rt.versions,
                "qc": entity["paths"]["filter_qc"],
                "parse_settings": {k: v for k, v in self.c.items() if k.startswith("PAIRTOOLS_") and "MEMORY" not in k},
                "dedup_settings": {k: v for k, v in self.c.items() if k.startswith("DEDUP_")}}
        preprocessing = [{"sample": name, "contract": self.plan["samples"][name]["preprocessing"],
                          "receipt": self.plan["entities"][f"sample:{name}"]["paths"]["blhic_qc"]}
                         for name in entity.get("source_rows", entity["members"])
                         if self.plan["samples"].get(name, {}).get("preprocessing")]
        if preprocessing:
            result["read_preprocessing"] = preprocessing
        if self.plan.get("workflow_mode") == "refilter-deduplicated":
            result.update(input_stage="deduplicated", alignment_performed=False, parsing_performed=False, deduplication_performed=False,
                          dedup_settings=None, parse_settings=None,
                          deduplication_status="Assumed from user-declared input; not independently verified or repaired")
            result["import_sources"] = [{"sample": row, "identity": self.plan["samples"][row]["input_identity"],
                "input_members": self.plan["samples"][row]["input_members"], "previous_filter_policy": self.plan["samples"][row]["input_header"]["previous_filter_policy"],
                "read_preprocessing": self.plan["samples"][row]["input_header"].get("read_preprocessing")}
                for row in entity["source_rows"]]
        elif entity.get("import_sources"):
            result.update(import_sources=entity["import_sources"], dedup_settings=None, deduplication_performed=False,
                          deduplication_status="Inherited user-declared deduplication; pooling does not deduplicate")
        return result

    def write_proof(self, temporary, entity, stage, key, destination, extra=None):
        proof = self.provenance(entity, stage, key, destination)
        proof["payload"] = identity(temporary["pairs"])
        proof["payload"]["path"] = str(destination.resolve())
        proof["completion_record"] = str(self.rt.store(entity).record_path(stage).resolve())
        proof.update(extra or {})
        json_write(temporary["provenance"], proof)

    def ensure_sorted(self, key):
        entity, keys = self.plan["entities"][key], self.keys[key]
        if entity.get("input_stage") == "deduplicated":
            raise RuntimeError("Alignment/parse stages are disabled in refilter-deduplicated mode")
        if entity["kind"] == "external_sorted":
            return Path(entity["paths"]["sorted"])
        paths, store = entity["paths"], self.rt.store(entity)
        destination = Path(paths["sorted"])
        outputs = {"pairs": destination, "provenance": proof_path(destination)}
        if store.ensure_reusable("sorted", keys["sorted"], outputs):
            return destination
        if entity["kind"] == "dedup":
            sources = [self.ensure_sorted(f"sample:{name}") for name in entity["members"]]
            def build(t, work):
                self.merge(sources, t["pairs"], entity, work)
                self.write_proof(t, entity, "sorted", keys["sorted"], destination)
            store.run("sorted", keys["sorted"], outputs, build)
            return destination
        row = self.plan["samples"][entity["name"]]
        ref = self.plan["references"][entity["reference_id"]]
        r1, r2 = self.ensure_trim(key)
        if truth(self.c["CHECK_TRIMMED_FASTQ_SYNC"]):
            from check_fastq_pair_names import check_pairs
            check_pairs(r1, r2)
        def align(t, work):
            self.rt.pipeline([
                ["bwa-mem2", "mem", *shlex.split(self.c["BWA_EXTRA_ARGS"]), "-t", self.c["THREADS_ALIGN"], ref["bwa_index_prefix"], str(r1), str(r2)],
                ["pairtools", "parse", "-c", ref["chrom_sizes"], "--assembly", ref["assembly"],
                 "--min-mapq", self.c["PAIRTOOLS_PARSE_MIN_MAPQ"], "--walks-policy", self.c["PAIRTOOLS_WALKS_POLICY"],
                 "--max-inter-align-gap", self.c["PAIRTOOLS_MAX_INTER_ALIGN_GAP"], "--max-molecule-size", self.c["PAIRTOOLS_MAX_MOLECULE_SIZE"],
                 "--report-alignment-end", "5", "--drop-sam", "--drop-seq", "--add-columns", "mapq"],
                ["pairtools", "sort", "--nproc", self.c["THREADS_SORT"], "--memory", self.c["PAIRTOOLS_SORT_MEMORY"], "--tmpdir", str(work)],
                [sys.executable, str(ROOT / "scripts/pairs_receipt.py"), "--report", str(work / "sorted_receipt.json")],
                ["bgzip", "-@", self.c["THREADS_COMPRESS"]]], t["pairs"], self.rt.log(entity, "align_parse_sort"))
            receipt = json.loads((work / "sorted_receipt.json").read_text())
            self.write_proof(t, entity, "sorted", keys["sorted"], destination, {"logical_data_sha256": receipt["data_sha256"], "records": receipt["records"]})
        store.run("sorted", keys["sorted"], outputs, align)
        return destination

    def merge(self, sources, destination, entity, work):
        if len(sources) == 1:
            shutil.copyfile(sources[0], destination)
        else:
            self.rt.pipeline([self.merge_command(sources, work), ["bgzip", "-@", self.c["THREADS_COMPRESS"]]],
                             destination, self.rt.log(entity, "merge"))

    def merge_command(self, sources, work):
        return [sys.executable, str(ROOT / "scripts/checked_merge.py"), "--threads", self.c["THREADS_MERGE"],
                "--memory", self.c["PAIRTOOLS_MERGE_MEMORY"], "--fan-in", self.c["MERGE_FAN_IN"], "--readers", self.c["MERGE_READERS"], "--max-inputs", self.c["MAX_MERGE_INPUTS"],
                "--tmpdir", str(work / "checked_merge"), "--report", str(work / "merge_receipt.json"), *map(str, sources)]

    def streaming(self, key):
        return (self.plan["entities"][key]["kind"] == "dedup" and self.c["MERGE_DEDUP_MODE"] == "stream"
                and not truth(self.c["KEEP_RAW_MERGED_PAIRS"]) and not self.available(key, "sorted"))

    def ensure_dedup(self, key):
        entity, keys = self.plan["entities"][key], self.keys[key]
        if entity.get("input_stage") == "deduplicated":
            raise RuntimeError("Deduplication is disabled in refilter-deduplicated mode")
        paths, store = entity["paths"], self.rt.store(entity)
        destination = Path(paths["dedup"])
        outputs = {"pairs": destination, "stats": Path(paths["dedup_stats"]), "provenance": proof_path(destination)}
        if store.ensure_reusable("dedup", keys["dedup"], outputs):
            return destination
        stream = self.streaming(key)
        sources = [self.ensure_sorted(f"sample:{name}") for name in entity["members"]] if stream else [self.ensure_sorted(key)]
        if stream:
            for source in sources:
                verify_proof_payload(source, read_proof(source))
        def dedup(t, work):
            self.rt.pipeline([self.merge_command(sources, work), ["pairtools", "dedup", "--backend", self.c["DEDUP_BACKEND"], "--method", self.c["DEDUP_METHOD"],
                "--max-mismatch", self.c["DEDUP_MAX_MISMATCH"], "--chunksize", self.c["DEDUP_CHUNKSIZE"], "--carryover", self.c["DEDUP_CARRYOVER"],
                "--output-stats", str(t["stats"])], ["bgzip", "-@", self.c["THREADS_COMPRESS"]]], t["pairs"], self.rt.log(entity, "dedup"))
            receipt = json.loads((work / "merge_receipt.json").read_text())
            totals = dict(line.rstrip().split("\t", 1) for line in t["stats"].read_text().splitlines() if "\t" in line)
            if int(totals.get("total", -1)) != receipt["records"]:
                raise ValueError("Merged input and deduplication total do not reconcile")
            self.write_proof(t, entity, "dedup", keys["dedup"], destination, {"merge_receipt": receipt, "streamed_joint_merge": stream})
        store.run("dedup", keys["dedup"], outputs, dedup)
        return destination

    def ensure_chrom_sizes(self, entity):
        ref, store = self.plan["references"][entity["reference_id"]], self.rt.store(entity)
        import re
        sizes = read_sizes(Path(ref["chrom_sizes"]), "reference")
        retained = {c: size for c, size in sizes.items() if re.fullmatch(entity["policy"]["canonical_regex"], c)}
        if not retained:
            raise ValueError("Configured chromosome policy retains no chromosomes")
        key = self.rt.key("chrom_sizes", {"reference": ref["identity"], "regex": entity["policy"]["canonical_regex"]})
        def write(t, work):
            t["sizes"].write_text("".join(f"{chrom}\t{length}\n" for chrom, length in retained.items()), encoding="utf-8")
        store.run("chrom_sizes", key, {"sizes": Path(entity["paths"]["chrom_sizes"])}, write)

    def ensure_accepted(self, key):
        entity, keys = self.plan["entities"][key], self.keys[key]
        if entity.get("alias_of"):
            source = self.plan["entities"][entity["alias_of"]]
            validate_pool_inputs([Path(source["paths"]["accepted"])])
            return Path(source["paths"]["accepted"])
        paths, store = entity["paths"], self.rt.store(entity)
        destination = Path(paths["accepted"])
        outputs = {"pairs": destination, "index": Path(str(destination) + ".px2"),
                   "qc": Path(paths["filter_qc"]), "provenance": proof_path(destination)}
        imported = entity.get("input_stage") == "deduplicated"
        if imported:
            outputs["import_qc"] = Path(paths["import_qc"])
        self.ensure_chrom_sizes(entity)
        if store.ensure_reusable("accepted", keys["accepted"], outputs):
            return destination
        if entity["kind"] == "pool":
            sources = [Path(self.plan["entities"][unit]["paths"]["accepted"]) for unit in entity["units"]]
            if not entity.get("legacy_pool"):
                validate_pool_inputs(sources)
            def accept(t, work):
                if entity.get("legacy_pool"):
                    raw = work / "legacy_source.pairs.gz"
                    self.merge(sources, raw, entity, work)
                    policy_file = Path(paths["policy"])
                    json_write(policy_file, entity["policy"])
                    ref = self.plan["references"][entity["reference_id"]]
                    self.rt.pipeline([[sys.executable, str(ROOT / "scripts/filter_contact_pairs.py"), "--input", str(raw),
                        "--policy", str(policy_file), "--chrom-sizes", ref["chrom_sizes"], "--matrix-chrom-sizes", paths["chrom_sizes"],
                        "--qc", str(t["qc"]), "--legacy-pool"], ["bgzip", "-@", self.c["THREADS_COMPRESS"]]], t["pairs"], self.rt.log(entity, "legacy_chromosome_filter"))
                else:
                    self.merge(sources, t["pairs"], entity, work)
                    source_qc = [json.loads(Path(self.plan["entities"][unit]["paths"]["filter_qc"]).read_text()) for unit in entity["units"]]
                    count = sum(q["exclusive_counts"]["accepted"] for q in source_qc)
                    json_write(t["qc"], {"policy": entity["policy"], "exclusive_counts": {"input": count, "accepted": count},
                        "accepted_cis": sum(q.get("accepted_cis", 0) for q in source_qc), "accepted_trans": sum(q.get("accepted_trans", 0) for q in source_qc),
                        "disjoint_input_units": entity["units"], "source_qc": [self.plan["entities"][unit]["paths"]["filter_qc"] for unit in entity["units"]]})
                self.rt.command(["pairix", "-f", str(t["pairs"])], self.rt.log(entity, "pairix"))
                self.write_proof(t, entity, "accepted", keys["accepted"], destination)
                if entity.get("legacy_pool"):
                    raw.unlink()
        else:
            source = Path(self.plan["samples"][entity["name"]]["deduplicated_pairs"]) if imported else self.ensure_dedup(key)
            policy_file = Path(paths["policy"])
            json_write(policy_file, entity["policy"])
            ref = self.plan["references"][entity["reference_id"]]
            digest = None
            if entity["policy"]["enzymes"]:
                with self.rt.digest_lock:
                    digest_key = signature({"reference": ref["identity"], "enzymes": entity["policy"]["enzyme_definitions"]})
                    if digest_key not in self.rt.digests:
                        self.rt.digests[digest_key], _ = prepare(Path(ref["fasta"]), Path(ref["chrom_sizes"]), entity["policy"]["enzyme_definitions"],
                                            Path(self.plan["roots"]["cache"]), ref["identity"]["fasta"].get("sha256"), check=self.rt.check)
                    digest = self.rt.digests[digest_key]
            def accept(t, work):
                command = [sys.executable, str(ROOT / "scripts/filter_contact_pairs.py"), "--input", str(source), "--policy", str(policy_file),
                           "--chrom-sizes", ref["chrom_sizes"], "--matrix-chrom-sizes", paths["chrom_sizes"], "--qc", str(t["qc"]),
                           "--batch-records", self.c["FILTER_BATCH_RECORDS"], "--batch-bytes", self.c["FILTER_BATCH_BYTES"], "--engine", self.c["CONTACT_FILTER_ENGINE"]]
                if digest:
                    command += ["--digest", str(digest)]
                commands = [command]
                if imported:
                    row = self.plan["samples"][entity["name"]]
                    validate_snapshot(row)
                    expected = work / "expected_input.json"
                    json_write(expected, row["input_identity"])
                    command[command.index("--input")+1] = "-"
                    commands = [[sys.executable, str(ROOT / "scripts/refilter_inputs.py"), "--input", str(source),
                        "--chrom-sizes", ref["chrom_sizes"], "--expected-identity", str(expected), "--qc", str(t["import_qc"])], command,
                        ["pairtools", "sort", "--nproc", self.c["THREADS_SORT"], "--memory", self.c["PAIRTOOLS_SORT_MEMORY"], "--tmpdir", str(work)]]
                commands.append(["bgzip", "-@", self.c["THREADS_COMPRESS"]])
                self.rt.pipeline(commands, t["pairs"], self.rt.log(entity, "filter"))
                if imported:
                    validate_snapshot(row)
                    qc = json.loads(t["qc"].read_text())
                    qc["import"] = json.loads(t["import_qc"].read_text())
                    qc["previous_filter_policy"] = row["input_header"]["previous_filter_policy"]
                    json_write(t["qc"], qc)
                    expected.unlink()
                if json.loads(t["qc"].read_text())["exclusive_counts"]["accepted"] == 0:
                    raise ValueError(f"No accepted contacts for {entity['id']}; inspect retained QC at {t['qc']} before matrix generation")
                self.rt.command(["pairix", "-f", str(t["pairs"])], self.rt.log(entity, "pairix"))
                self.write_proof(t, entity, "accepted", keys["accepted"], destination, {"digest": str(digest) if digest else None})
        store.run("accepted", keys["accepted"], outputs, accept)
        return destination

    def stats(self, key):
        entity = self.plan["entities"][key]
        output = Path(entity["paths"]["stats"])
        sources = []
        if entity["kind"] == "pool" and not entity.get("legacy_pool") and self.c["POOL_STATS_MODE"] == "merge":
            sources = [self.plan["entities"][unit] for unit in entity["units"]]
            validate_pool_inputs([Path(e["paths"]["accepted"]) for e in sources])
            if any("stats" not in e["paths"] for e in sources):
                sources = []  # standalone external pools retain the direct-scan fallback
            for source in sources:
                self.stats(source["id"])
        command = (["pairtools", "stats", "--merge", *[e["paths"]["stats"] for e in sources]] if sources else
                   ["pairtools", "stats", "--nproc-in", "1", entity["paths"]["accepted"]])
        self.rt.store(entity).run("stats", self.rt.key("stats", self.keys[key]["accepted"], ["pairtools"]), {"stats": output},
            lambda t, work: self.rt.command(command, self.rt.log(entity, "stats"), output=t["stats"]))


def validate_pool_inputs(sources: list[Path]) -> list[dict]:
    units, policies, members, proofs, imported_inputs = set(), set(), set(), [], set()
    for source in sources:
        proof = read_proof(source)
        verify_proof_payload(source, proof)
        if proof.get("stage") != "accepted" or not proof.get("poolable"):
            raise ValueError(f"Pool needs accepted disjoint units, not a pre-dedup stream or a row assigned to joint dedup: {source}")
        current = set(proof["dedup_units"])
        for imported in proof.get("import_sources", []):
            item = imported["identity"]
            source_key = (item["device"], item["inode"]) if item.get("inode") else item["path"]
            if source_key in imported_inputs:
                raise ValueError("Pool inputs overlap: an imported physical source would be counted twice")
            imported_inputs.add(source_key)
        if current & units or set(proof["members"]) & members:
            raise ValueError("Pool inputs overlap: a dedup unit/member would be counted twice")
        units.update(current)
        members.update(proof["members"])
        policies.add(signature({"reference": proof["reference_identity"], "policy": proof["policy"]}))
        proofs.append(proof)
    if len(policies) != 1:
        raise ValueError("Pool inputs have incompatible reference or processing policies")
    return proofs


def verify_proof_payload(source: Path, proof: dict) -> None:
    expected = proof.get("payload", {})
    actual = identity(source)
    # Paths may be deliberately relocated; identity still requires byte length and mtime.
    if any(actual.get(key) != expected.get(key) for key in ("size", "mtime_ns")):
        raise ValueError(f"Pairs/provenance payload identity mismatch: {source}")
    marker = Path(proof.get("completion_record", ""))
    if not marker.is_file():
        raise ValueError(f"Pairs lack a completed stage manifest: {source}")
    record = json.loads(marker.read_text())
    if record.get("status") != "success" or record.get("signature") != proof.get("signature"):
        raise ValueError(f"Pairs stage is incomplete, failed, or stale: {source}")
    if record.get("outputs", {}).get("pairs") != actual:
        raise ValueError(f"Stage manifest does not identify this physical pairs artifact: {source}")
    for item in record["outputs"].values():
        path = Path(item["path"])
        if not path.is_file() or identity(path) != item:
            raise ValueError(f"Pairs bundle is incomplete or has changed: {path}")
