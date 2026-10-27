#!/usr/bin/env python3
"""Read-only access checks as the calling user; not a functional scientific test."""
from __future__ import annotations

import argparse
from dataclasses import replace
import getpass
import os
from pathlib import Path
import shutil
import sys

from reference_registry import Registry, validate_reference_files
from workflow_common import ROOT, VERSION, absolute


def check_access(env_prefix: Path, *, release: Path = ROOT, registry: Path | None = None,
                 reference_ids: list[str] | None = None, jar: Path | None = None,
                 refilter: bool = False, hic: bool = True) -> tuple[list[str], list[str], list[str]]:
    checks, warnings, errors = [], [], []

    def readable(path: Path, label: str):
        try:
            if not path.is_file() or not os.access(path, os.R_OK):
                raise PermissionError(f"not a readable file: {path}")
            with path.open("rb") as handle:
                handle.read(1)
            checks.append(f"{label}: {path}")
        except OSError as exc:
            errors.append(f"{label}: {exc}")

    try:
        for directory in (release / "scripts", release / "config"):
            if not directory.is_dir() or not os.access(directory, os.R_OK | os.X_OK):
                raise PermissionError(f"not a readable/searchable directory: {directory}")
            for path in directory.iterdir():
                if path.is_file() and path.suffix in {".py", ".sh", ".json", ".tsv", ".conf"}:
                    readable(path, "release file")
    except OSError as exc:
        errors.append(f"release: {exc}")

    search_path = str(env_prefix / "bin") + os.pathsep + os.environ.get("PATH", "")
    # Python must come from the chosen shared prefix, never from a user's fallback.
    python = env_prefix / "bin/python"
    if not python.is_file() or not os.access(python, os.R_OK | os.X_OK):
        errors.append(f"shared Python is not readable/executable: {python}")
    else:
        checks.append(f"shared Python: {python}")
    tools = ["pairtools", "bgzip", "pairix", "sort", "cooler"]
    if not refilter:
        tools += ["bwa-mem2", "fastp", "samtools"]
    if hic:
        tools += ["java"]
    for name in tools:
        found = shutil.which(name, path=search_path)
        if not found:
            errors.append(f"required executable is inaccessible: {name}")
        else:
            checks.append(f"executable {name}: {found}")
            resolved = Path(found).resolve().as_posix()
            if resolved.startswith(("/home/", "/root/")):
                warnings.append(f"{name} resolves inside a home directory: {resolved}; check every account's access")
    for name in ("cooltools", "multiqc"):
        found = shutil.which(name, path=search_path)
        if found:
            checks.append(f"optional executable {name}: {found}")
        else:
            warnings.append(f"{name} is inaccessible; required if its analysis/reporting options are enabled")

    if hic:
        if jar is None:
            warnings.append("Juicer JAR was not checked; supply --juicer-jar when .hic output is needed")
        else:
            readable(jar, "Juicer JAR")
    if registry is None:
        warnings.append("Reference paths were not checked; supply --registry and optionally --reference")
    else:
        try:
            refs = Registry(registry)
            selected = [refs.resolve(name) for name in reference_ids] if reference_ids else list(refs.references.values())
            for ref in selected:
                ref = replace(ref, **{key: absolute(getattr(ref, key), registry.resolve().parent)
                    for key in ("fasta", "bwa_index_prefix", "chrom_sizes")},
                    phasing_track=absolute(ref.phasing_track, registry.resolve().parent) if ref.phasing_track else "")
                failures = validate_reference_files(ref, require_bwa=not refilter)
                errors.extend(failures)
                if not failures:
                    checks.append(f"reference {ref.reference_id}: readable existing files; no index created")
        except (OSError, ValueError) as exc:
            errors.append(f"reference registry: {exc}")
    return checks, warnings, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--reference", action="append", help="Check only these registered references (repeatable)")
    parser.add_argument("--juicer-jar", type=Path)
    parser.add_argument("--refilter", action="store_true", help="Skip aligner/FASTQ-tool and BWA-index access checks")
    parser.add_argument("--no-hic", action="store_true", help="Skip Java/Juicer checks if .hic output is disabled")
    args = parser.parse_args()
    if args.reference and not args.registry:
        parser.error("--reference requires --registry")
    user = getpass.getuser()
    if hasattr(os, "geteuid"):
        import pwd
        user = f"{pwd.getpwuid(os.geteuid()).pw_name} (euid={os.geteuid()})"
    print(f"microc2tracks {VERSION}; checking access as {user}")
    checks, warnings, errors = check_access(args.env_prefix, registry=args.registry,
        reference_ids=args.reference, jar=args.juicer_jar, refilter=args.refilter, hic=not args.no_hic)
    for level, rows in (("OK", checks), ("WARNING", warnings), ("ERROR", errors)):
        for row in rows:
            print(f"{level}: {row}")
    if errors:
        return 1
    print("Access checks passed for this account and the checked paths only. "
          "This does not certify every account or replace imports, preflight and real-tool tests.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
