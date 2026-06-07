#!/usr/bin/env python3
"""Normalize text inputs before Bash and CSV parsing."""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path


SMART_TRANSLATION = str.maketrans(
    {
        "\ufeff": "",
        "\u00a0": " ",
        "\u2018": "'",
        "\u2019": "'",
        "\u201c": '"',
        "\u201d": '"',
    }
)


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def read_utf8(path: Path) -> str:
    try:
        return path.read_bytes().decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        fail(f"{path}: expected UTF-8-compatible text input: {exc}")


def normalize_line_endings(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalize_config(text: str) -> str:
    text = normalize_line_endings(text).translate(SMART_TRANSLATION)
    lines = text.split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(line.rstrip() for line in lines) + "\n"


def normalize_samplesheet(text: str) -> str:
    text = normalize_line_endings(text).translate(SMART_TRANSLATION)
    reader = csv.reader(io.StringIO(text))
    rows = []

    for row in reader:
        if not row:
            continue
        rows.append([field.strip() for field in row])

    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerows(rows)
    return output.getvalue()


def changed_features(original: bytes, original_text: str, new_text: str) -> list[str]:
    features: list[str] = []

    if original.startswith(b"\xef\xbb\xbf") or "\ufeff" in original_text:
        features.append("UTF-8 BOM")
    if b"\r" in original:
        features.append("CRLF/CR line endings")
    if "\u00a0" in original_text:
        features.append("non-breaking spaces")
    if any(char in original_text for char in "\u2018\u2019\u201c\u201d"):
        features.append("smart quotes")
    if original_text != new_text and not features:
        features.append("surrounding whitespace/CSV formatting")

    return features


def sanitize_file(path: Path, kind: str, check_only: bool) -> bool:
    if not path.exists():
        fail(f"{path}: file does not exist")
    if not path.is_file():
        fail(f"{path}: not a regular file")

    original = path.read_bytes()
    original_text = read_utf8(path)

    if kind == "config":
        new_text = normalize_config(original_text)
    elif kind == "samplesheet":
        new_text = normalize_samplesheet(original_text)
    else:
        fail(f"Unsupported kind: {kind}")

    changed = original != new_text.encode("utf-8")
    if not changed:
        print(f"Input hygiene OK: {path}")
        return False

    features = ", ".join(changed_features(original, original_text, new_text))
    if check_only:
        print(f"NEEDS FIX: {path} ({features})", file=sys.stderr)
        return True

    try:
        path.write_text(new_text, encoding="utf-8", newline="")
    except OSError as exc:
        fail(f"{path}: could not write normalized file: {exc}")

    print(f"Sanitized input: {path} ({features})")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize microc2tracks config and sample-sheet text inputs in place."
    )
    parser.add_argument("--kind", choices=["config", "samplesheet"], required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("path", type=Path)
    args = parser.parse_args()

    changed = sanitize_file(args.path, args.kind, args.check_only)
    if changed and args.check_only:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
