import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams["font.family"] = "Arial"
rcParams["font.size"] = 16
rcParams["axes.labelsize"] = 18
rcParams["axes.titlesize"] = 14
rcParams["xtick.labelsize"] = 16
rcParams["ytick.labelsize"] = 16
rcParams["mathtext.fontset"] = "cm"
rcParams["svg.fonttype"] = "none"

METRICS = ["ipTM", "ipSAE", "ipSAE_min", "pDockQ", "pDockQ2", "LIS"]

ROOT = Path(__file__).parent
# MMSEQS2_DIR is the directory containing the results from alphafold3_tools (MMSeqs2)
MMSEQS2_DIR = ROOT / "AF3mmseqs2"
# ORIGINAL_DIR is the directory containing the results from AlphaFold3 (HMMER).
ORIGINAL_DIR = ROOT / "AF3original"

X_LABEL = "alphafold3_tools (MMSeqs2)"
Y_LABEL = "AlphaFold3 (HMMER3)"


def load_all_metrics(root: Path) -> dict[str, dict[str, float]]:
    aggregated: dict[str, dict[str, float]] = {metric: {} for metric in METRICS}
    for json_path in sorted(root.glob("*/complexmetrics.json")):
        subdir = json_path.parent.name
        with open(json_path) as f:
            data = json.load(f)
        for metric, entries in data.items():
            if metric in aggregated:
                for entry in entries:
                    for pair_key, value in entry.items():
                        aggregated[metric][f"{subdir}/{pair_key}"] = value
    return aggregated


def add_second_axis(
    ax,
    xdata,
    ydata,
    color="#0072BC",
    bins=30,
    hist_alpha=0.35,
    show_hist=True,
):
    divider = make_axes_locatable(ax)
    ax_kdex = divider.append_axes("top", size="12%", pad=0.01, sharex=ax)
    ax_kdey = divider.append_axes("right", size="12%", pad=0.01, sharey=ax)

    x = pd.Series(xdata).dropna()
    y = pd.Series(ydata).dropna()

    if show_hist:
        ax_kdex.hist(x, bins=bins, color=color, alpha=hist_alpha, edgecolor="none")
        ax_kdey.hist(
            y,
            bins=bins,
            orientation="horizontal",
            color=color,
            alpha=hist_alpha,
            edgecolor="none",
        )

    for a in (ax_kdex, ax_kdey):
        a.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        a.set_facecolor("none")
        a.grid(False)
        a.set_xlabel("")
        a.set_ylabel("")
        for spine in a.spines.values():
            spine.set_visible(False)


def plot_metric(
    ax: plt.Axes, base: dict[str, float], nomsa: dict[str, float], title: str
) -> None:
    common = sorted(set(base) & set(nomsa))
    x = [base[m] for m in common]
    y = [nomsa[m] for m in common]

    ax.scatter(x, y, s=8, zorder=3, color="steelblue")

    if len(x) >= 2:
        corr = np.corrcoef(x, y)[0, 1]
        ax.text(
            0.05,
            0.95,
            f"r = {corr:.2f}",
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment="top",
        )

    ax.plot(
        [0, 1],
        [0, 1],
        "k--",
        linewidth=1,
        zorder=2,
        label="x = y",
    )
    ax.set(
        xlim=(0, 1),
        xticks=np.linspace(0, 1, 6),
        ylim=(0, 1),
        yticks=np.linspace(0, 1, 6),
    )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel(X_LABEL, fontsize=12)
    ax.set_ylabel(Y_LABEL, fontsize=12)
    ax.tick_params(labelsize=12)
    add_second_axis(ax, x, y)


def plot_diff_boxplot(
    ax: plt.Axes, base: dict[str, float], orig: dict[str, float], title: str
) -> None:
    common = sorted(set(base) & set(orig))
    diffs = [base[k] - orig[k] for k in common]

    if diffs:
        ax.boxplot(
            [diffs],
            vert=False,
            patch_artist=True,
            boxprops=dict(facecolor="steelblue", alpha=0.5),
            medianprops=dict(color="darkred", linewidth=1.5),
            whiskerprops=dict(linewidth=1),
            capprops=dict(linewidth=1),
            flierprops=dict(marker=".", markersize=3, alpha=0.4),
        )

    ax.axvline(0, color="k", linestyle="--", linewidth=1, zorder=2)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks(np.arange(-1.0, 1.01, 0.2))
    ax.set_yticks([])
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("MMSeqs2 - HMMER3", fontsize=10)
    ax.tick_params(labelsize=10)


def export_excel(
    mmseqs2_data: dict[str, dict[str, float]],
    original_data: dict[str, dict[str, float]],
    out: Path,
) -> None:
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for metric in METRICS:
            base = mmseqs2_data[metric]
            hmmer3data = original_data[metric]
            common = sorted(set(base) & set(hmmer3data))
            df = pd.DataFrame(
                {
                    "name": common,
                    X_LABEL: [base[k] for k in common],
                    Y_LABEL: [hmmer3data[k] for k in common],
                }
            )
            df.to_excel(writer, sheet_name=metric, index=False)
    print(f"Saved: {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--excel", action="store_true", help="Export data to xlsx")
    args = parser.parse_args()

    mmseqs2_data = load_all_metrics(MMSEQS2_DIR)
    original_data = load_all_metrics(ORIGINAL_DIR)

    if args.excel:
        export_excel(mmseqs2_data, original_data, ROOT / "complexmetrics_data.xlsx")

    for metric in METRICS:
        base = mmseqs2_data[metric]
        orig = original_data[metric]
        common = sorted(set(base) & set(orig))
        count = sum(1 for k in common if abs(base[k] - orig[k]) > 0.2)
        print(f"{metric}: {count} pairs with |x - y| > 0.2")

    fig = plt.figure(figsize=(15, 16))
    fig.suptitle("Metric comparison", fontsize=13, y=0.995)
    gs = GridSpec(
        4,
        3,
        figure=fig,
        height_ratios=[3, 3, 1, 1],
        hspace=0.7,
        wspace=0.45,
        left=0.07,
        right=0.97,
        top=0.97,
        bottom=0.05,
    )

    scatter_axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    box_axes = [fig.add_subplot(gs[i + 2, j]) for i in range(2) for j in range(3)]

    for ax, metric in zip(scatter_axes, METRICS):
        plot_metric(ax, mmseqs2_data[metric], original_data[metric], metric)

    for ax, metric in zip(box_axes, METRICS):
        plot_diff_boxplot(ax, mmseqs2_data[metric], original_data[metric], metric)
    out = ROOT / "complexmetrics_comparison.svg"
    plt.savefig(out, dpi=200, bbox_inches="tight", format="svg")
    print(f"Saved: {out}")
    plt.clf()
    plt.close()


if __name__ == "__main__":
    main()
