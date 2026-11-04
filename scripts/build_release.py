#!/usr/bin/env python3
"""Build a manifest-only source archive; never install software or run analyses."""
import argparse
import gzip
import hashlib
from pathlib import Path, PurePosixPath
import tarfile

from workflow_common import ROOT, VERSION


def build(directory, date):
    files = (ROOT / "UPLOAD_MANIFEST.txt").read_text().splitlines()
    if len(files) != len(set(files)):
        raise ValueError("Duplicate release manifest entry")
    for name in files:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or not (ROOT / name).is_file() or (ROOT / name).is_symlink():
            raise ValueError(f"Unsafe/missing source entry: {name}")
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / f"microc2tracks_{VERSION}_{date}.tar.gz"
    with destination.open("xb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
            for name in sorted(files):
                path = ROOT / name
                info = archive.gettarinfo(str(path), arcname="microc2tracks/" + name)
                info.mode, info.uid, info.gid, info.mtime = (0o755 if path.suffix == ".sh" else 0o644), 0, 0, 0
                info.uname = info.gname = "root"
                with path.open("rb") as handle:
                    archive.addfile(info, handle)
    with tarfile.open(destination, "r:gz") as archive:
        members = {m.name: m for m in archive if m.isfile()}
        if set(members) != {"microc2tracks/" + name for name in files}:
            raise ValueError("Archive membership mismatch")
        for name in files:
            if archive.extractfile(members["microc2tracks/" + name]).read() != (ROOT / name).read_bytes():
                raise ValueError(f"Archive payload mismatch: {name}")
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    checksum = destination.with_suffix("").with_suffix(".sha256")
    with checksum.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{digest}  {destination.name}\n")
    print(f"Verified {len(files)} source files: {destination}")
    print(f"SHA256: {digest}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--date", required=True, help="YYYYMMDD build label")
    args = p.parse_args()
    if len(args.date) != 8 or not args.date.isdigit():
        p.error("Date must be YYYYMMDD")
    build(args.output_dir, args.date)
