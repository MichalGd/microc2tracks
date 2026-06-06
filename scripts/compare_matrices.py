#!/usr/bin/env python3
"""Plot two cooler/mcool matrices side by side for the same genomic region."""

from __future__ import annotations

import argparse
from pathlib import Path


def cooler_uri(path: str, resolution: int) -> str:
    if "::" in path:
        return path
    if path.endswith(".mcool"):
        return f"{path}::resolutions/{resolution}"
    return path


def fetch_matrix(uri: str, region: str, balanced: bool):
    import cooler
    import numpy as np

    clr = cooler.Cooler(uri)
    try:
      matrix = clr.matrix(balance=balanced).fetch(region)
    except Exception:
      if balanced:
          print(f"WARNING: balanced fetch failed for {uri}; using raw counts")
          matrix = clr.matrix(balance=False).fetch(region)
      else:
          raise
    return np.asarray(matrix, dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-a", required=True)
    parser.add_argument("--matrix-b", required=True)
    parser.add_argument("--label-a", default="A")
    parser.add_argument("--label-b", default="B")
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--raw", action="store_true", help="Use raw counts instead of balanced values")
    parser.add_argument("--no-log1p", dest="log1p", action="store_false", help="Plot raw/balanced values without log1p transform")
    parser.set_defaults(log1p=True)
    args = parser.parse_args()

    import matplotlib.pyplot as plt
    import numpy as np

    uri_a = cooler_uri(args.matrix_a, args.resolution)
    uri_b = cooler_uri(args.matrix_b, args.resolution)

    mat_a = fetch_matrix(uri_a, args.region, balanced=not args.raw)
    mat_b = fetch_matrix(uri_b, args.region, balanced=not args.raw)

    if args.log1p:
        mat_a = np.log1p(mat_a)
        mat_b = np.log1p(mat_b)

    finite = np.concatenate([mat_a[np.isfinite(mat_a)], mat_b[np.isfinite(mat_b)]])
    vmax = np.nanpercentile(finite, 99) if finite.size else 1.0

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    for ax, matrix, label in [(axes[0], mat_a, args.label_a), (axes[1], mat_b, args.label_b)]:
        image = ax.imshow(matrix, cmap="magma", origin="lower", vmin=0, vmax=vmax)
        ax.set_title(label)
        ax.set_xlabel(args.region)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(image, ax=axes, shrink=0.75, label="log1p balanced contact" if not args.raw else "log1p raw count")
    fig.suptitle(f"{args.region} at {args.resolution} bp")
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
