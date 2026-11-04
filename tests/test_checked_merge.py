"""Real POSIX reader/GNU-sort tests; no production files or scientific tools."""
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from checked_merge import merge, reader


class ReaderTests(unittest.TestCase):
    def test_crc_failure_never_commits_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, fifo, receipt = root / "pairs.gz", root / "out", root / "receipt.json"
            source.write_bytes(gzip.compress(b"#columns: x\na\tb\n")[:-5])
            with self.assertRaises((EOFError, OSError)):
                reader(source, fifo, receipt)
            self.assertFalse(receipt.exists())

    def test_nonterminated_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pairs.gz"
            source.write_bytes(gzip.compress(b"#columns: x\na\tb"))
            with self.assertRaisesRegex(ValueError, "Truncated"):
                reader(source, root / "out", root / "receipt")


@unittest.skipIf(os.name == "nt", "POSIX FIFO/GNU tests run on Linux")
class MergeTests(unittest.TestCase):
    def fixture(self, root, count=11):
        inputs, expected = [], []
        for i in range(count):
            source = root / f"part{i}.pairs.gz"
            data = [f"id{i:03}\tchr1\t10\tchr1\t20\t+\t-\tUU\t{30+i}\t60\n",
                    f"end{i:03}\tchr2\t11\tchr2\t22\t+\t-\tUU\t60\t60\n"]
            source.write_bytes(gzip.compress(("#columns: ignored\n" + "".join(data)).encode()))
            inputs.append(source)
            expected += data
        return inputs, sorted(expected, key=lambda line: (line.split("\t")[1], line.split("\t")[3],
            int(line.split("\t")[2]), int(line.split("\t")[4]), line.split("\t")[7], line))

    def test_more_than_eight_inputs_matches_gnu_tie_order_and_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources, expected = self.fixture(root)
            with (root / "merged").open("wb") as out, patch("checked_merge.headers", return_value=["#test"]), patch("checked_merge.sys.stdout", SimpleNamespace(buffer=out)):
                merge(sources, root / "work", root / "receipt.json", "1M", 2, 8)
            self.assertEqual((root / "merged").read_text().splitlines(True)[1:], expected)
            receipt = json.loads((root / "receipt.json").read_text())
            self.assertEqual(receipt["records"], 22)
            self.assertEqual(len(receipt["reader_receipts"]), 11)
            self.assertFalse(list((root / "work").glob("*.fifo")))

    def test_late_reader_failure_prevents_merge_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources, _ = self.fixture(root, 2)
            sources[1].write_bytes(sources[1].read_bytes()[:-5])
            with (root / "merged").open("wb") as out, patch("checked_merge.headers", return_value=["#test"]), patch("checked_merge.sys.stdout", SimpleNamespace(buffer=out)):
                with self.assertRaisesRegex(RuntimeError, "failed"):
                    merge(sources, root / "work", root / "receipt.json", "1M", 1, 8)
            self.assertFalse((root / "receipt.json").exists())


if __name__ == "__main__":
    unittest.main()
