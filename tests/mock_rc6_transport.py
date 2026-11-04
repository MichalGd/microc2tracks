"""Adapt old orchestration mocks to RC6 transport; not scientific validation."""
import gzip
import json
from pathlib import Path


def adapt(fake):
    def pipeline(rt, commands, output, log_path, timeout=None):
        first = commands[0]
        receipt_command = next((c for c in commands if len(c) > 1 and str(c[1]).endswith("pairs_receipt.py")), None)
        if receipt_command:
            result = fake(rt, [c for c in commands if c is not receipt_command], output, log_path, timeout)
            import hashlib
            with gzip.open(output, "rb") as inp:
                data = b"".join(line for line in inp if not line.startswith(b"#"))
            Path(receipt_command[receipt_command.index("--report") + 1]).write_text(json.dumps({"records": data.count(b"\n"), "data_sha256": hashlib.sha256(data).hexdigest()}))
            return result
        if len(first) > 1 and str(first[1]).endswith("checked_merge.py"):
            sources = first[first.index("--report") + 2:]
            records = 0
            for source in sources:
                with gzip.open(source, "rt") as handle:
                    records += sum(bool(line.strip()) and not line.startswith("#") for line in handle)
            receipt = Path(first[first.index("--report") + 1])
            receipt.parent.mkdir(parents=True, exist_ok=True)
            receipt.write_text(json.dumps({"records": records, "reader_receipts": []}))
            if len(commands) > 1 and commands[1][:2] == ["pairtools", "dedup"]:
                materialized = receipt.parent / "mock_merged.pairs.gz"
                fake(rt, [["pairtools", "merge", *sources]], materialized, log_path, timeout)
                return fake(rt, [[*commands[1], str(materialized)], *commands[2:]], output, log_path, timeout)
            return fake(rt, [["pairtools", "merge", *sources]], output, log_path, timeout)
        if len(first) > 1 and str(first[1]).endswith("juicer_export.py"):
            from matrix_products import export_juicer
            count = export_juicer(Path(first[first.index("--input") + 1]), Path(output))
            Path(first[first.index("--report") + 1]).write_text(json.dumps({"records": count}))
            return
        if len(first) > 1 and str(first[1]).endswith("pool_cooler.py"):
            sources = first[first.index("--expected-total") + 2:]
            total = sum(int(Path(p.split("::")[0]).read_text()) for p in sources)
            Path(first[first.index("--output") + 1]).write_text(str(total))
            return
        if first[:3] == ["pairtools", "stats", "--merge"]:
            total = sum(int(Path(p).read_text().splitlines()[0].split("\t")[1]) for p in first[3:])
            Path(output).write_text(f"total\t{total}\n")
            return
        return fake(rt, commands, output, log_path, timeout)
    return pipeline
