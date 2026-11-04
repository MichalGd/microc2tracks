"""Dependency-free Linux process-tree measurements; RSS is conservative (shared pages sum)."""
from pathlib import Path


def process_tree_metrics(pid):
    records = {}
    for path in Path("/proc").glob("[0-9]*/status"):
        try:
            values = dict(line.split(":", 1) for line in path.read_text().splitlines() if ":" in line)
            records[int(path.parent.name)] = (int(values["PPid"]), int(values.get("VmRSS", "0 kB").split()[0]) * 1024,
                                             int(values.get("Threads", "1")))
        except (OSError, ValueError, KeyError):
            continue
    selected = {pid}
    while True:
        expanded = selected | {p for p, (parent, _, _) in records.items() if parent in selected}
        if expanded == selected:
            break
        selected = expanded
    return {"rss_bytes": sum(records[p][1] for p in selected if p in records),
            "processes": len(selected), "threads": sum(records[p][2] for p in selected if p in records),
            "measurement": "Linux summed process RSS; not PSS or an OS-enforced quota"}
