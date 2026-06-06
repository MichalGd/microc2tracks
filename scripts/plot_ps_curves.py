#!/usr/bin/env python3
"""Plot P(s) curves from one or more cooltools expected-cis tables."""

from __future__ import annotations

import argparse
from pathlib import Path


def choose_column(frame, candidates: list[str]) -> str:
    for column in candidates:
        if column in frame.columns:
            return column
    raise ValueError(f"None of these columns found: {', '.join(candidates)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", nargs="+", required=True, help="Expected TSV files")
    parser.add_argument("--labels", nargs="+", required=True, help="Labels matching --expected")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    import matplotlib.pyplot as plt
    import pandas as pd

    if len(args.expected) != len(args.labels):
        raise SystemExit("--expected and --labels must have the same length")

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 4.5), constrained_layout=True)
    for path, label in zip(args.expected, args.labels):
        frame = pd.read_csv(path, sep="\t")
        dist_col = choose_column(frame, ["dist_bp", "dist"])
        value_col = choose_column(frame, ["balanced.avg.smoothed", "balanced.avg", "count.avg"])
        data = frame[[dist_col, value_col]].dropna()
        data = data[data[dist_col] > 0]
        data = data[data[value_col] > 0]
        ax.loglog(data[dist_col], data[value_col], label=label)

    ax.set_xlabel("genomic separation (bp)")
    ax.set_ylabel("contact frequency")
    ax.legend(frameon=False)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
