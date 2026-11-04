"""Explicit same-run checkpoint recovery; never invent missing FASTQ identities."""
import json
from pathlib import Path

from workflow_common import identity, truth


def sorted_checkpoint(runtime, name):
    """Recover a known upstream signature only after checking its old contract.

    A retained downstream library checkpoint may replace deleted sorted payloads;
    its own payload and completion bundle must still pass normal reuse checks.
    """
    if not truth(runtime.c["RESUME_FROM_CHECKPOINTS"]):
        raise ValueError(f"FASTQ missing for {name}; use --resume-checkpoints with a managed run, or refilter-deduplicated for legacy pairs")
    old = runtime.previous_plan
    current = runtime.plan
    version = old.get("workflow_version")
    if (not old or version not in {"3.0.0-rc4", "3.0.0-rc5", "3.0.0-rc6"}
            or (version == "3.0.0-rc4" and not runtime.rc4_hashes)
            or (version == "3.0.0-rc5" and not runtime.rc5_hashes)):
        raise ValueError("Checkpoint resume requires a recognized prior run/source fingerprint")
    previous = old.get("samples", {}).get(name)
    row = current["samples"][name]
    if not previous or any(previous.get(k) != row.get(k) for k in ("fastq_r1", "fastq_r2", "reference_genome", "dedup_group")):
        raise ValueError(f"Checkpoint source membership/paths changed for {name}")
    for tool in ("pairtools", "bwa-mem2", "fastp", "numpy", "scipy", "pandas", "pysam"):
        if tool in old.get("tool_versions", {}) and old["tool_versions"][tool] != runtime.versions.get(tool):
            raise ValueError(f"Checkpoint tool/library version changed: {tool}")
    row["fastq_identities"] = {**previous.get("fastq_identities", {}), **row.get("fastq_identities", {})}
    keys = ["RUN_FASTP", "FASTP_EXTRA_ARGS", "BWA_EXTRA_ARGS"] + [k for k in runtime.c if k.startswith("PAIRTOOLS_") and "MEMORY" not in k]
    if any(old["settings"].get(k) != runtime.c[k] for k in keys):
        raise ValueError(f"Alignment/parse/preprocessing contract changed for {name}; FASTQs are required")
    ref = row["reference_genome"]
    if old["references"][ref]["identity"] != current["references"][ref]["identity"]:
        raise ValueError("Reference identity changed; missing FASTQs cannot be remapped")
    saved = old.get("checkpoint_keys", {}).get(f"sample:{name}", {})
    record_path = Path(current["entities"][f"sample:{name}"]["paths"]["done"]) / "sorted.json"
    record = json.loads(record_path.read_text()) if record_path.is_file() else {}
    key = saved.get("sorted") or (record.get("signature") if record.get("status") == "success" else None)
    if not key:
        raise ValueError(f"No recorded sorted-stage lineage for {name}; recover earlier managed manifests or import declared pairs explicitly")
    # If some raw sources remain, do not silently mix a changed file with old lineage.
    for label in ("fastq_r1", "fastq_r2"):
        path = Path(row[label])
        expected = previous.get("fastq_identities", {}).get(label)
        if path.exists() and expected and identity(path) != expected:
            raise ValueError(f"Remaining FASTQ changed: {path}")
    return key
