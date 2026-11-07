#!/usr/bin/env python3
"""Build only the pinned optional trimmer; no network or system installation."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

from workflow_common import ROOT, json_write, sha256


def build(destination, compiler="g++"):
    source = ROOT / "vendor/trimLinker"
    manifest = json.loads((source / "source.json").read_text())
    actual = {name: sha256(source / name) for name in manifest["sha256"]}
    if actual != manifest["sha256"]:
        raise ValueError("Bundled trimLinker source checksum differs from the pinned manifest")
    tool = shutil.which(compiler)
    if not tool:
        raise RuntimeError("Building the optional trimmer needs a C++ compiler, OpenMP and zlib headers")
    destination = destination.resolve()
    if destination.exists() or Path(str(destination) + ".build.json").exists():
        raise ValueError("Build destination already exists; use a new versioned path")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".trimLinker-build-", dir=destination.parent) as attempt:
        temporary = Path(attempt) / "trimLinker"
        command = [tool, "-Wall", "-O2", "-Wno-char-subscripts", "-fopenmp", str(source / "trimLinker.cpp"), "-lz", "-o", str(temporary)]
        print(shlex.join(command), flush=True)
        subprocess.run(command, check=True)
        receipt = {"source_commit": manifest["source_commit"], "source_hashes": actual,
                   "binary_sha256": sha256(temporary), "compiler": subprocess.check_output([tool, "--version"], text=True).splitlines()[0],
                   "command": command}
        temporary.chmod(0o755)
        os.replace(temporary, destination)
        json_write(Path(str(destination) + ".build.json"), receipt)
    print(f"Built pinned optional trimmer: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "tools/blhic/trimLinker")
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    build(args.output, args.compiler)
