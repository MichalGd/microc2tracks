"""Small shared configuration, profile and file helpers. No sample metadata is evaluated."""
from __future__ import annotations

import csv
import contextlib
import gzip
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VERSION = "3.0.0-rc2"
DEFAULTS = {
    "WORKFLOW_MODE": "fastq",
    "OUTDIR": "results", "TMPDIR": "", "PAIRS_OUTDIR": "", "MATRIX_OUTDIR": "",
    "TRIMMED_OUTDIR": "", "REFERENCE_CACHE_DIR": "", "LARGE_OUTPUT_STAGING": "destination",
    "REQUIRED_STORAGE_MOUNTS": "", "MIN_FREE_OUTPUT_GB": "10", "MIN_FREE_SCRATCH_GB": "20",
    "STORAGE_CHECK_INTERVAL": "5", "REFERENCE_REGISTRY": str(ROOT / "config/references.tsv"),
    "DEFAULT_REFERENCE_ID": "mm39", "MAPQ_THRESHOLD": "30", "MICROC_MIN_CIS_DIST": "1000",
    "HIC_MIN_CIS_DIST": "0", "HIC_PROTOCOL": "legacy", "HIC_RESTRICTION_ENZYMES": "",
    "HIC_REMOVE_SAME_FRAGMENT": "true", "HIC_REMOVE_ADJACENT_INWARD": "false",
    "HIC_PROTOCOL_REGISTRY": str(ROOT / "config/hic_protocols.json"),
    "FILTER_CANONICAL_CHROMS": "true", "DEDUP_BACKEND": "scipy", "DEDUP_METHOD": "max",
    "DEDUP_MAX_MISMATCH": "3", "DEDUP_CHUNKSIZE": "10000", "DEDUP_CARRYOVER": "100",
    "PAIRTOOLS_PARSE_MIN_MAPQ": "1", "PAIRTOOLS_WALKS_POLICY": "5unique",
    "PAIRTOOLS_MAX_INTER_ALIGN_GAP": "20", "PAIRTOOLS_MAX_MOLECULE_SIZE": "750",
    "BWA_EXTRA_ARGS": "-SP5M", "FASTP_EXTRA_ARGS": "--detect_adapter_for_pe",
    "FASTP_TIMEOUT_SECONDS": "0", "CHECK_TRIMMED_FASTQ_SYNC": "false", "CHECK_FASTQ_GZIP": "false",
    "THREADS_FASTP": "8", "THREADS_ALIGN": "8", "THREADS_SORT": "4", "THREADS_MATRIX": "4",
    "THREADS_MERGE": "4", "THREADS_HIC_NORM": "4", "PAIRTOOLS_SORT_MEMORY": "8G",
    "PAIRTOOLS_MERGE_MEMORY": "8G", "JAVA_HEAP": "8g", "JUICER_TOOLS_JAR": "",
    "MAX_PARALLEL_SAMPLES": "1", "MAX_PARALLEL_FASTP": "1", "MAX_PARALLEL_MATRIX": "1",
    "MAX_PARALLEL_HIC": "1", "MAX_PARALLEL_MERGES": "1", "LOCK_POLL_SECONDS": "5",
    "BASE_RESOLUTION": "1000", "RESOLUTIONS": "1000,5000,10000,100000",
    "MCOOL_RESOLUTIONS": "", "HIC_RESOLUTIONS": "", "HIC_NORMALIZATIONS": "VC,VC_SQRT,KR,SCALE",
    "RUN_FASTP": "true", "RUN_MCOOL": "true", "RUN_HIC": "true", "BUILD_SAMPLE_MATRICES": "true",
    "RUN_MERGE_TECHNICAL_REPLICATES": "true", "RUN_MERGE_ALL_SAMPLES": "false",
    "ALL_SAMPLES_MERGE_NAME": "all_samples_merged", "RUN_PRELIM_DOWNSTREAM": "true",
    "RUN_MULTIQC": "true", "RUN_GLOBAL_MULTIQC": "true", "RUN_FINAL_REPORT": "true",
    "KEEP_SORTED_PAIRS": "true", "KEEP_DEDUP_PAIRS": "true", "KEEP_TRIMMED_FASTQ": "false",
    "KEEP_JUICER_PAIRS": "false", "KEEP_SINGLE_RES_COOL": "false", "KEEP_RAW_MCOOL": "false",
    "KEEP_RAW_HIC": "false", "KEEP_RAW_MERGED_PAIRS": "false", "CLEAN_TMP_ON_SUCCESS": "true",
    "PUBLIC_HIC_DIR": "", "PUBLIC_HIC_BASE_URL": "", "COOLTOOLS_THREADS": "4",
    "TAD_RESOLUTION": "10000", "COMPARTMENT_RESOLUTION": "100000", "LOOP_RESOLUTION": "5000",
    "INSULATION_WINDOWS_BP": "50000 100000 250000", "TAD_BOUNDARY_THRESHOLD": "Li",
    "TAD_BOUNDARY_WINDOW_BP": "100000", "TAD_MIN_SIZE_BP": "40000", "TAD_MIN_BOUNDARY_STRENGTH": "",
    "PHASING_TRACK": "", "REFERENCE_FASTA": "", "BWA_INDEX_PREFIX": "", "CHROM_SIZES": "",
    "CANONICAL_CHROMS_REGEX": "", "GENOME_ASSEMBLY": "",
}


def truth(value: object) -> bool:
    if str(value).lower() in {"true", "1", "yes"}:
        return True
    if str(value).lower() in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"Expected true/false, got {value!r}")


def safe_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Unsafe identifier {value!r}; start with a letter/digit; use letters, digits, ._- only")
    return value


def absolute(value: str | Path, base: Path | None = None) -> str:
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else (base or Path.cwd()) / path).resolve())


def load_config(path: Path) -> dict[str, str]:
    """Configs remain trusted Bash files, as before; never source sample-sheet fields."""
    # The runner also accepts old, clean Bash configs without new keys.
    bash = os.environ.get("MICROC2TRACKS_BASH", "bash")
    from sanitize_text_inputs import normalize_config
    with tempfile.TemporaryDirectory(prefix="microc2tracks_config_") as directory:
        clean = Path(directory) / "config.conf"
        original_bytes = path.read_bytes()
        original = original_bytes.decode("utf-8-sig")
        normalized = normalize_config(original)
        source = path.resolve()
        if original_bytes != normalized.encode("utf-8"):
            if "BASH_SOURCE" in original:
                raise ValueError("Config with BASH_SOURCE needs clean LF text first; run sanitize_text_inputs.py explicitly")
            clean.write_text(normalized, encoding="utf-8", newline="\n")
            source = clean
        result = subprocess.run([bash, "-c", 'set -ae; source "$1"; env -0', "microc2tracks", source.as_posix()],
                                check=True, stdout=subprocess.PIPE)
    exported = dict(item.split("=", 1) for item in result.stdout.decode().split("\0") if "=" in item)
    config = {key: exported.get(key, value) for key, value in DEFAULTS.items()}
    if not exported.get("DEFAULT_REFERENCE_ID") and exported.get("GENOME_ASSEMBLY"):
        config["DEFAULT_REFERENCE_ID"] = exported["GENOME_ASSEMBLY"]
    for key in ("MCOOL_RESOLUTIONS", "HIC_RESOLUTIONS"):
        config[key] = config[key] or config["RESOLUTIONS"]
    for key in ("REFERENCE_REGISTRY", "HIC_PROTOCOL_REGISTRY"):
        value = Path(config[key])
        candidates = [value, path.resolve().parent / value, ROOT / value, ROOT / "config" / value.name]
        config[key] = str(next((p.resolve() for p in candidates if p.is_file()), value))
    validate_config(config)
    return config


def validate_config(c: dict[str, str]) -> None:
    if c["WORKFLOW_MODE"] not in {"fastq", "refilter-deduplicated"}:
        raise ValueError("WORKFLOW_MODE must be fastq/refilter-deduplicated")
    if not c["OUTDIR"].strip():
        raise ValueError("OUTDIR must not be empty")
    for key in c:
        if key.startswith(("RUN_", "KEEP_", "BUILD_", "CHECK_", "CLEAN_", "FILTER_", "HIC_REMOVE_")):
            truth(c[key])
    positive = [k for k in c if k.startswith(("THREADS_", "MAX_PARALLEL_"))]
    positive += ["BASE_RESOLUTION", "DEDUP_CHUNKSIZE", "DEDUP_CARRYOVER", "STORAGE_CHECK_INTERVAL"]
    for key in positive:
        if not c[key].isdigit() or int(c[key]) <= 0:
            raise ValueError(f"{key} must be a positive integer")
    for key in ("MAPQ_THRESHOLD", "MICROC_MIN_CIS_DIST", "HIC_MIN_CIS_DIST", "DEDUP_MAX_MISMATCH",
                "PAIRTOOLS_PARSE_MIN_MAPQ", "PAIRTOOLS_MAX_INTER_ALIGN_GAP", "PAIRTOOLS_MAX_MOLECULE_SIZE", "FASTP_TIMEOUT_SECONDS"):
        if not c[key].isdigit():
            raise ValueError(f"{key} must be a non-negative integer")
    if c["DEDUP_BACKEND"] not in {"scipy", "cython"} or c["DEDUP_METHOD"] not in {"max", "sum"}:
        raise ValueError("DEDUP_BACKEND must be scipy/cython; DEDUP_METHOD must be max/sum")
    if c["PAIRTOOLS_WALKS_POLICY"] not in {"mask", "5any", "5unique", "3any", "3unique", "all"}:
        raise ValueError("Invalid PAIRTOOLS_WALKS_POLICY")
    if c["LARGE_OUTPUT_STAGING"] not in {"destination", "scratch"}:
        raise ValueError("LARGE_OUTPUT_STAGING must be destination/scratch")
    for key in ("MIN_FREE_OUTPUT_GB", "MIN_FREE_SCRATCH_GB"):
        if not 0 <= float(c[key]) < float("inf"):
            raise ValueError(f"{key} must be a finite non-negative number")
    for key in ("MCOOL_RESOLUTIONS", "HIC_RESOLUTIONS"):
        values = c[key].split(",")
        if any(not x.isdigit() or int(x) < int(c["BASE_RESOLUTION"]) or int(x) % int(c["BASE_RESOLUTION"]) for x in values):
            raise ValueError(f"{key} must contain positive multiples of BASE_RESOLUTION")
    if truth(c["RUN_PRELIM_DOWNSTREAM"]) and truth(c["RUN_MCOOL"]):
        for key in ("TAD_RESOLUTION", "COMPARTMENT_RESOLUTION"):
            if c[key] not in c["MCOOL_RESOLUTIONS"].split(","):
                raise ValueError(f"{key} must be present in MCOOL_RESOLUTIONS when preliminary downstream is enabled")


def json_write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def tsv_write(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def signature(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def identity(path: Path, content: bool = False) -> dict:
    stat = path.stat()
    result = {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if content:
        result["sha256"] = sha256(path)
    return result


def text_open(path: Path):
    if str(path) == "-":
        return contextlib.nullcontext(sys.stdin)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix.lower() == ".gz" else path.open(encoding="utf-8")


def registry_data(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    for name, enzyme in data["enzymes"].items():
        safe_id(name)
        motif, offset = enzyme["motif"], enzyme["cut_offset"]
        if not re.fullmatch("[ACGT]+", motif) or motif != motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]:
            raise ValueError(f"Enzyme {name}: this release requires a palindromic A/C/G/T motif")
        if not isinstance(offset, int) or not 0 <= offset <= len(motif):
            raise ValueError(f"Enzyme {name}: invalid cut offset")
    return data


def resolve_policy(row: dict, c: dict) -> dict:
    data = registry_data(Path(c["HIC_PROTOCOL_REGISTRY"]))
    assay = row["assay"]
    profile = (row.get("hic_protocol") or c["HIC_PROTOCOL"]) if assay == "hic" else "microc"
    enzymes: list[str] = []
    if assay == "microc":
        if row.get("hic_protocol") or row.get("restriction_enzymes"):
            raise ValueError("Micro-C rows must leave hic_protocol/restriction_enzymes empty")
    else:
        if profile not in data["profiles"]:
            raise ValueError(f"Unknown Hi-C profile {profile!r}")
        preset = data["profiles"][profile]
        value = row.get("restriction_enzymes") or c["HIC_RESTRICTION_ENZYMES"] or preset.get("enzymes", "")
        if preset["restriction_filtering"]:
            aliases = {key.lower(): key for key in data["enzymes"]}
            try:
                enzymes = sorted({aliases[x.strip().lower()] for x in value.split("+") if x.strip()})
            except KeyError as exc:
                raise ValueError(f"Unsupported enzyme {exc.args[0]!r}; supported: {', '.join(aliases.values())}") from exc
            if not enzymes:
                raise ValueError("Restriction Hi-C requires explicit restriction_enzymes/HIC_RESTRICTION_ENZYMES")
            if profile == "dovetail_hic" and enzymes != ["DpnII"]:
                raise ValueError("dovetail_hic is the original DpnII kit; use restriction_hic for other digests")
        elif value:
            raise ValueError("legacy Hi-C cannot silently ignore an enzyme; select restriction_hic")
    return {"assay": assay, "profile": profile, "enzymes": enzymes,
            "enzyme_definitions": {name: data["enzymes"][name] for name in enzymes},
            "mapq": int(c["MAPQ_THRESHOLD"]),
            "min_cis_dist": int(c["MICROC_MIN_CIS_DIST"] if assay == "microc" else c["HIC_MIN_CIS_DIST"]),
            "remove_same_fragment": bool(enzymes) and truth(c["HIC_REMOVE_SAME_FRAGMENT"]),
            "remove_adjacent_inward": bool(enzymes) and truth(c["HIC_REMOVE_ADJACENT_INWARD"]),
            "filter_canonical": truth(c["FILTER_CANONICAL_CHROMS"])}
