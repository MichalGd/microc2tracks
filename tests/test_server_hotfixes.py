#!/usr/bin/env python3
"""Small regressions for the source-identical installed-server hotfixes."""

from __future__ import annotations

import csv
import io
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from summarize_run import html_table  # noqa: E402
from validate_samplesheet import NORMALIZED_COLUMNS, write_normalized  # noqa: E402


class ServerHotfixTests(unittest.TestCase):
    def test_normalized_csv_and_tsv_use_lf_with_empty_last_field(self) -> None:
        row = dict.fromkeys(NORMALIZED_COLUMNS, "")
        row.update(sample="example_sample", assay="microc", reference_genome="mm39")
        with tempfile.TemporaryDirectory() as directory:
            for delimiter, suffix in ((",", "csv"), ("\t", "tsv")):
                with self.subTest(format=suffix):
                    path = Path(directory) / f"normalized.{suffix}"
                    write_normalized(path, [row], delimiter)
                    content = path.read_bytes()
                    self.assertNotIn(b"\r", content)
                    self.assertTrue(content.endswith(delimiter.encode() + b"\n"))
                    parsed = list(csv.DictReader(io.StringIO(content.decode()), delimiter=delimiter))
                    self.assertEqual(parsed, [row])
                    self.assertEqual(parsed[0]["merge_group"], "")

    def test_html_accepts_none_and_missing_cells(self) -> None:
        rendered = html_table("Samples", [{"sample": "one", "optional": None}, {"sample": "two"}])
        self.assertEqual(rendered.count("<td></td>"), 2)
        self.assertNotIn("None", rendered)

    def test_html_still_escapes_special_characters_and_accepts_numbers(self) -> None:
        rendered = html_table("<Samples>", [{"sample": "<x>&", "count": 0}])
        self.assertIn("&lt;Samples&gt;", rendered)
        self.assertIn("<td>&lt;x&gt;&amp;</td>", rendered)
        self.assertIn("<td>0</td>", rendered)


if __name__ == "__main__":
    unittest.main(verbosity=2)
