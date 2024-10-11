#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger
from matplotlib import colormaps
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from complexbuilder.common.log import log_setup

mpl.use("Agg")

log_setup(level="DEBUG")
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def load_metrics(json_path: Path) -> dict:
    with json_path.open() as f:
        data = json.load(f)
    return data


def get_pair_value(
    data: dict,
    id1: str,
    id2: str,
    metric: str,
) -> float | None:
    """
    JSON is a key of "id1_id2" format. Since the order may be different, both are searched.
    """
    # Should be lowercase
    id1 = id1.lower()
    id2 = id2.lower()
    key1 = f"{id1}_{id2}"
    key2 = f"{id2}_{id1}"

    def lookup(container: dict) -> float | None:
        if key1 in container and isinstance(container[key1], dict):
            if metric in container[key1]:
                return container[key1][metric]
        if key2 in container and isinstance(container[key2], dict):
            if metric in container[key2]:
                return container[key2][metric]
        return None

    val = lookup(data)
    if val is not None:
        return val

    for sub in data.values():
        if isinstance(sub, dict):
            val = lookup(sub)
            if val is not None:
                return val
    return None


def build_matrix(
    data: dict,
    xs: list[str],
    ys: list[str],
    metric: str,
) -> np.ndarray:
    mat = np.full((len(ys), len(xs)), np.nan, dtype=float)
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            val = get_pair_value(data, x, y, metric)
            if val is not None:
                mat[j, i] = val
    return mat


def draw_highlight_boxes(ax, cells, color="#ffffff", lw=1.0):
    """
    cells: iterable of (row, col) in 0-based index
    """
    for r, c in cells:
        rect = Rectangle(
            (c - 0.5, r - 0.5),
            1,
            1,
            fill=False,
            edgecolor=color,
            linewidth=lw,
            linestyle="dashed",
        )
        ax.add_patch(rect)


def plot_heatmap(
    matrices: list[np.ndarray],
    matrices_titles: list[str],
    xs: list[str],
    label_xs: list[str] | None,
    ys: list[str],
    label_ys: list[str] | None,
    outpath: Path,
    cmap: str = "YlGn",
    dpi: int = 200,
):
    """Plot heatmap from matrices."""
    outpath = Path(outpath)
    fig, axes = plt.subplots(
        1, len(matrices), figsize=(4.0 * len(matrices), 4.0), dpi=dpi
    )

    # axes_list is a unified list of axes
    if len(matrices) == 1:
        axes_list = [axes]
    else:
        axes_list = list(axes.ravel())

    cmap_obj = colormaps.get_cmap(cmap).copy()
    cmap_obj.set_bad(color="#d3d3d3")

    if label_xs is not None and len(label_xs) == len(xs):
        new_label_xs = [
            f"{label}\n({x})" for label, x in zip(label_xs, xs, strict=False)
        ]
    else:
        new_label_xs = xs
    if label_ys is not None and len(label_ys) == len(ys):
        new_label_ys = [
            f"{label}\n({y})" for label, y in zip(label_ys, ys, strict=False)
        ]
    else:
        new_label_ys = ys

    for ax, matrix, title in zip(axes_list, matrices, matrices_titles, strict=False):
        im = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap_obj,
            vmin=0,
            vmax=1,
        )

        ax.set_xticks(range(len(xs)))
        ax.set_yticks(range(len(ys)))
        if label_xs is not None:
            ax.set_xticklabels(new_label_xs, rotation=45, ha="right", fontsize=7)
        else:
            ax.set_xticklabels(xs, rotation=45, ha="right", fontsize=7)
        if label_ys is not None:
            ax.set_yticklabels(new_label_ys, fontsize=7)
        else:
            ax.set_yticklabels(ys, fontsize=7)
        ax.set_aspect("equal")

        ax.set_title(title)
        draw_highlight_boxes(ax, [(0, 0)])
        draw_highlight_boxes(ax, [(1, 1)])
        draw_highlight_boxes(ax, [(2, 2)])
        draw_highlight_boxes(ax, [(3, 3)])
        if matrix.size <= 400:
            for j in range(matrix.shape[0]):
                for i in range(matrix.shape[1]):
                    val = matrix[j, i]
                    if not np.isnan(val):
                        ax.text(
                            i,
                            j,
                            f"{val:.3f}",
                            ha="center",
                            va="center",
                            color="white" if val > 0.5 else "black",
                            fontsize=10,
                        )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=dpi)
    plt.clf()
    plt.close()


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate ipSAE / ipTM heatmap from JSON data."
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input json path including ipTM/ipSAE metrics",
    )
    p.add_argument(
        "-o",
        "--output",
        required=True,
        help="output image file path (png or svg format)",
    )
    p.add_argument(
        "-x",
        "--xids",
        nargs="+",
        required=True,
        help="ID list for X axis (columns)",
    )
    p.add_argument(
        "-lx",
        "--label_xids",
        nargs="+",
        required=False,
        help="Label list for X axis (columns). Should be equal to the number of X IDs.",
    )
    p.add_argument(
        "-y",
        "--yids",
        nargs="+",
        required=True,
        help="ID list for Y axis (rows)",
    )
    p.add_argument(
        "-ly",
        "--label_yids",
        nargs="+",
        required=False,
        help="Label list for Y axis (rows). Should be equal to the number of Y IDs.",
    )
    p.add_argument(
        "-m",
        "--metrics",
        choices=["ipTM+ipSAE+ipSAE_min", "ipTM+ipSAE", "ipSAE", "ipTM"],
        default="ipTM+ipSAE+ipSAE_min",
        help="One or more metrics to plot. Choose from ipTM+ipSAE+ipSAE_min, ipTM+ipSAE, ipTM, or ipSAE.",
    )
    p.add_argument(
        "-c",
        "--colormap",
        default="YlGn",
        help="Colormap to use for the heatmap (default: YlGn)",
    )
    p.add_argument(
        "-d",
        "--dpi",
        type=int,
        default=200,
        help="DPI (dots per inch) for the output image (default: 200)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    # Validate the number of -lx and -ly arguments
    if args.label_xids is not None and (len(args.label_xids) != len(args.xids)):
        raise SystemExit("Number of X labels must match number of X IDs.")
    if args.label_yids is not None and (len(args.label_yids) != len(args.yids)):
        raise SystemExit("Number of Y labels must match number of Y IDs.")

    json_path = Path(args.input)

    data = load_metrics(json_path)
    matrices = []
    if args.metrics == "ipTM+ipSAE+ipSAE_min":
        matrices.append(build_matrix(data, args.xids, args.yids, "ipTM"))
        matrices.append(build_matrix(data, args.xids, args.yids, "ipSAE"))
        matrices.append(build_matrix(data, args.xids, args.yids, "ipSAE_min"))
        matrices_titles = ["ipTM", "ipSAE", "ipSAE_min"]
    elif args.metrics == "ipTM+ipSAE":
        matrices.append(build_matrix(data, args.xids, args.yids, "ipTM"))
        matrices.append(build_matrix(data, args.xids, args.yids, "ipSAE"))
        matrices_titles = ["ipTM", "ipSAE"]
    elif args.metrics == "ipSAE":
        matrices.append(build_matrix(data, args.xids, args.yids, "ipSAE"))
        matrices_titles = ["ipSAE"]
    elif args.metrics == "ipTM":
        matrices.append(build_matrix(data, args.xids, args.yids, "ipTM"))
        matrices_titles = ["ipTM"]

    plot_heatmap(
        matrices,
        matrices_titles,
        args.xids,
        args.label_xids,
        args.yids,
        args.label_yids,
        args.output,
        cmap=args.colormap,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
