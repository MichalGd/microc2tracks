#!/usr/bin/env python3
"""Plot and compare cooltools saddle outputs with a shared color scale."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def sanitize_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value).strip("_") or "sample"


def parse_input(value: str) -> tuple[str, Path]:
    if "=" in value:
        label, path = value.split("=", 1)
        label = label.strip()
    else:
        path = value
        label = Path(path).name.replace(".saddledump.npz", "")
    if not label:
        fail(f"Empty label in input: {value}")
    return label, Path(path)


def load_saddle_matrix(path: Path) -> np.ndarray:
    if not path.exists():
        fail(f"Saddle NPZ not found: {path}")

    with np.load(path) as data:
        if "saddledata" in data:
            matrix = np.asarray(data["saddledata"], dtype=float)
        elif "S" in data and "C" in data:
            with np.errstate(divide="ignore", invalid="ignore"):
                matrix = np.asarray(data["S"], dtype=float) / np.asarray(data["C"], dtype=float)
        else:
            candidates = []
            for key in data.files:
                array = np.asarray(data[key])
                if array.ndim == 2 and array.shape[0] == array.shape[1] and np.issubdtype(array.dtype, np.number):
                    candidates.append((key, array.astype(float)))
            if not candidates:
                fail(f"Could not find a square saddle matrix in {path}; keys: {data.files}")
            matrix = candidates[0][1]

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        fail(f"Saddle matrix is not square in {path}: shape={matrix.shape}")
    return matrix


def prepare_plot_matrix(matrix: np.ndarray, drop_edge_bins: bool) -> np.ndarray:
    if drop_edge_bins and matrix.shape[0] > 2:
        matrix = matrix[1:-1, 1:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(matrix)


def finite_values(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.isfinite(matrix)]


def corner_metrics(matrix: np.ndarray, corner_bins: int) -> dict[str, float]:
    n_bins = matrix.shape[0]
    k = max(1, min(corner_bins, n_bins // 2))

    same = np.concatenate(
        [
            matrix[:k, :k].ravel(),
            matrix[-k:, -k:].ravel(),
        ]
    )
    opposite = np.concatenate(
        [
            matrix[:k, -k:].ravel(),
            matrix[-k:, :k].ravel(),
        ]
    )
    same = same[np.isfinite(same)]
    opposite = opposite[np.isfinite(opposite)]

    same_mean = float(np.nanmean(same)) if same.size else math.nan
    opposite_mean = float(np.nanmean(opposite)) if opposite.size else math.nan
    contrast = same_mean - opposite_mean
    strength_ratio = 2.0**contrast if math.isfinite(contrast) else math.nan

    return {
        "corner_bins": float(k),
        "same_corner_mean_log2oe": same_mean,
        "opposite_corner_mean_log2oe": opposite_mean,
        "corner_contrast_log2": contrast,
        "saddle_strength_ratio": strength_ratio,
    }


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    n = min(a.shape[0], b.shape[0])
    a = a[:n, :n]
    b = b[:n, :n]
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return math.nan
    x = a[mask]
    y = b[mask]
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    n = min(a.shape[0], b.shape[0])
    a = a[:n, :n]
    b = b[:n, :n]
    mask = np.isfinite(a) & np.isfinite(b)
    if not mask.any():
        return math.nan
    return float(np.sqrt(np.nanmean((a[mask] - b[mask]) ** 2)))


def write_matrix_plot(path: Path, label: str, matrix: np.ndarray, vlim: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    image = ax.imshow(
        matrix,
        origin="lower",
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        aspect="equal",
        interpolation="nearest",
    )
    ax.set_box_aspect(1)
    ax.set_xlabel("PC1 rank bins")
    ax.set_ylabel("PC1 rank bins")
    ax.set_title(label)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log2 observed / expected")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def write_panel_plot(path: Path, labels: list[str], matrices: list[np.ndarray], vlim: float, dpi: int) -> None:
    n = len(matrices)
    fig_width = max(5.0 * n, 6.0)
    fig, axes = plt.subplots(1, n, figsize=(fig_width, 5.5), constrained_layout=True)
    if n == 1:
        axes = [axes]

    image = None
    for ax, label, matrix in zip(axes, labels, matrices, strict=True):
        image = ax.imshow(
            matrix,
            origin="lower",
            cmap="coolwarm",
            vmin=-vlim,
            vmax=vlim,
            aspect="equal",
            interpolation="nearest",
        )
        ax.set_box_aspect(1)
        ax.set_xlabel("PC1 rank bins")
        ax.set_ylabel("PC1 rank bins")
        ax.set_title(label)

    if image is not None:
        cbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label("log2 observed / expected")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def write_difference_plot(path: Path, label_a: str, label_b: str, a: np.ndarray, b: np.ndarray, dpi: int) -> None:
    n = min(a.shape[0], b.shape[0])
    diff = a[:n, :n] - b[:n, :n]
    finite = finite_values(diff)
    vlim = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
    if not math.isfinite(vlim) or vlim == 0:
        vlim = 1.0

    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    image = ax.imshow(
        diff,
        origin="lower",
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        aspect="equal",
        interpolation="nearest",
    )
    ax.set_box_aspect(1)
    ax.set_xlabel("PC1 rank bins")
    ax.set_ylabel("PC1 rank bins")
    ax.set_title(f"{label_a} minus {label_b}")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("delta log2 observed / expected")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare cooltools saddle .saddledump.npz files using square plots "
            "with a shared color scale and objective corner/correlation metrics."
        )
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="label=path/to/file.saddledump.npz; repeat for each sample",
    )
    parser.add_argument("--out-dir", required=True, help="output directory for plots and metrics")
    parser.add_argument("--vlim", type=float, default=None, help="fixed symmetric color limit")
    parser.add_argument(
        "--vlim-percentile",
        type=float,
        default=98.0,
        help="percentile of absolute values used for shared color limit when --vlim is absent",
    )
    parser.add_argument(
        "--corner-fraction",
        type=float,
        default=0.2,
        help="fraction of rank bins used for each corner metric",
    )
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument(
        "--keep-edge-bins",
        action="store_true",
        help="keep qrange underflow/overflow bins; default drops the outermost bins",
    )
    args = parser.parse_args()

    items = [parse_input(value) for value in args.input]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = [label for label, _path in items]
    matrices = [
        prepare_plot_matrix(load_saddle_matrix(path), drop_edge_bins=not args.keep_edge_bins)
        for _label, path in items
    ]

    all_values = np.concatenate([finite_values(matrix) for matrix in matrices if finite_values(matrix).size])
    if args.vlim is not None:
        vlim = args.vlim
    elif all_values.size:
        vlim = float(np.nanpercentile(np.abs(all_values), args.vlim_percentile))
    else:
        vlim = 1.0
    if not math.isfinite(vlim) or vlim <= 0:
        vlim = 1.0

    metrics_path = out_dir / "saddle_metrics.tsv"
    with metrics_path.open("w", newline="") as handle:
        fieldnames = [
            "label",
            "path",
            "n_bins",
            "vlim",
            "matrix_mean_log2oe",
            "matrix_min_log2oe",
            "matrix_max_log2oe",
            "corner_bins",
            "same_corner_mean_log2oe",
            "opposite_corner_mean_log2oe",
            "corner_contrast_log2",
            "saddle_strength_ratio",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for (label, path), matrix in zip(items, matrices, strict=True):
            finite = finite_values(matrix)
            corner_bins = max(1, int(round(matrix.shape[0] * args.corner_fraction)))
            row = {
                "label": label,
                "path": str(path),
                "n_bins": matrix.shape[0],
                "vlim": f"{vlim:.8g}",
                "matrix_mean_log2oe": f"{float(np.nanmean(finite)):.8g}" if finite.size else "nan",
                "matrix_min_log2oe": f"{float(np.nanmin(finite)):.8g}" if finite.size else "nan",
                "matrix_max_log2oe": f"{float(np.nanmax(finite)):.8g}" if finite.size else "nan",
            }
            row.update({key: f"{value:.8g}" for key, value in corner_metrics(matrix, corner_bins).items()})
            writer.writerow(row)

    pair_path = out_dir / "saddle_pairwise.tsv"
    with pair_path.open("w", newline="") as handle:
        fieldnames = [
            "label_a",
            "label_b",
            "pearson_direct",
            "pearson_pc1_signflip",
            "pearson_best",
            "rmse_direct_log2",
            "rmse_pc1_signflip_log2",
            "corner_contrast_delta_a_minus_b",
            "saddle_strength_ratio_a_over_b",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        metrics = []
        for matrix in matrices:
            corner_bins = max(1, int(round(matrix.shape[0] * args.corner_fraction)))
            metrics.append(corner_metrics(matrix, corner_bins))

        for i in range(len(matrices)):
            for j in range(i + 1, len(matrices)):
                direct = pearson(matrices[i], matrices[j])
                flipped_matrix = matrices[j][::-1, ::-1]
                flipped = pearson(matrices[i], flipped_matrix)
                direct_rmse = rmse(matrices[i], matrices[j])
                flipped_rmse = rmse(matrices[i], flipped_matrix)
                contrast_delta = metrics[i]["corner_contrast_log2"] - metrics[j]["corner_contrast_log2"]
                strength_ratio = (
                    metrics[i]["saddle_strength_ratio"] / metrics[j]["saddle_strength_ratio"]
                    if metrics[j]["saddle_strength_ratio"]
                    else math.nan
                )
                writer.writerow(
                    {
                        "label_a": labels[i],
                        "label_b": labels[j],
                        "pearson_direct": f"{direct:.8g}",
                        "pearson_pc1_signflip": f"{flipped:.8g}",
                        "pearson_best": f"{max(direct, flipped):.8g}",
                        "rmse_direct_log2": f"{direct_rmse:.8g}",
                        "rmse_pc1_signflip_log2": f"{flipped_rmse:.8g}",
                        "corner_contrast_delta_a_minus_b": f"{contrast_delta:.8g}",
                        "saddle_strength_ratio_a_over_b": f"{strength_ratio:.8g}",
                    }
                )

    for label, matrix in zip(labels, matrices, strict=True):
        write_matrix_plot(out_dir / f"{sanitize_label(label)}.shared_scale.png", label, matrix, vlim, args.dpi)
    write_panel_plot(out_dir / "saddle_shared_scale_panel.png", labels, matrices, vlim, args.dpi)

    if len(matrices) == 2:
        write_difference_plot(
            out_dir / f"{sanitize_label(labels[0])}_minus_{sanitize_label(labels[1])}.png",
            labels[0],
            labels[1],
            matrices[0],
            matrices[1],
            args.dpi,
        )

    print(f"Shared color limit: +/-{vlim:.6g}")
    print(f"Wrote {metrics_path}")
    print(f"Wrote {pair_path}")
    print(f"Wrote plots under {out_dir}")


if __name__ == "__main__":
    main()
