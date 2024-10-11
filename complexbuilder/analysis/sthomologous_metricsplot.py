#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def collect_values(
    data: dict, rmsd_threshold: float
) -> Tuple[List[float], List[float]]:
    iptm_values: List[float] = []
    ipsae_values: List[float] = []
    for _, pairs in data.items():
        if not isinstance(pairs, dict):
            continue
        for _, metrics in pairs.items():
            try:
                rmsd = metrics.get("RMSD")
                if rmsd is None or rmsd > rmsd_threshold:
                    continue
                iptm = metrics.get("ipTM")
                ipsae = metrics.get("ipSAE")
                if iptm is not None:
                    iptm_values.append(float(iptm))
                if ipsae is not None:
                    ipsae_values.append(float(ipsae))
            except Exception:
                continue
    return iptm_values, ipsae_values


def plot_violin(
    iptm_values: List[float],
    ipsae_values: List[float],
    rmsd_threshold: float,
    out: Path,
):
    """
    Plot a seaborn violin plot (ipTM vs ipSAE) without dots.

    Shows distributions for ipTM and ipSAE (RMSD-filtered) side by side.
    Medians are drawn (inner='box').
    Sample sizes (n=) annotated above violins.
    """
    if not iptm_values and not ipsae_values:
        raise ValueError("No data to plot.")

    df_plot = pd.DataFrame(
        {
            "score": iptm_values + ipsae_values,
            "metric": (["ipTM"] * len(iptm_values)) + (["ipSAE"] * len(ipsae_values)),
        }
    )

    fig, ax = plt.subplots(figsize=(4.0, 4.0), dpi=300)
    seaborn.violinplot(
        data=df_plot,
        x="metric",
        y="score",
        hue="metric",
        inner="box",
        legend=False,
        cut=0,
        palette={"ipTM": "#0072BC", "ipSAE": "#FFC740"},
        linewidth=1.0,
        ax=ax,
    )

    # Annotate n
    for cat, grp in df_plot.groupby("metric"):
        x_pos = 0 if cat == "ipTM" else 1
        ax.text(
            x_pos,
            0.98,
            f"$n={len(grp)}$",
            ha="center",
            va="bottom",
            fontsize=8,
            clip_on=False,
        )

    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.set_title(f"ipTM & ipSAE (RMSD ≤ {rmsd_threshold} Å)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)


def main():
    parser = argparse.ArgumentParser(
        description="Plot ipTM / ipSAE violin plots filtered by RMSD. "
        "Usage: .venv/bin/python3.12 -m complexbuilder.analysis.sthomologous_metricsplot -i /path/to/rmsd.json "
        "-t 2.0 -o /path/to/outputplot.svg"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="input JSON file path (must contain RMSD, ipTM, ipSAE).",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--rmsd-threshold",
        type=float,
        default=2.0,
        help="RMSD threshold (Å). Default: %(default)s",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output SVG file path (if not specified, will be saved as out_{threshold}.svg in the same directory as the JSON file)",
    )

    args = parser.parse_args()

    json_path = Path(args.input)
    if not json_path.is_file():
        raise SystemExit(f"Input JSON not found: {json_path}")

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = json_path.parent / f"out_{args.rmsd_threshold}.svg"

    with json_path.open() as f:
        data = json.load(f)

    iptm_values, ipsae_values = collect_values(data, args.rmsd_threshold)

    iptm_median = pd.Series(iptm_values).median() if iptm_values else float("nan")
    print(f"ipTM median: {iptm_median:.2f}")
    ipsae_median = pd.Series(ipsae_values).median() if ipsae_values else float("nan")
    print(f"ipSAE median: {ipsae_median:.2f}")

    if not iptm_values and not ipsae_values:
        raise SystemExit(f"No entries with RMSD <= {args.rmsd_threshold} found.")

    plot_violin(iptm_values, ipsae_values, args.rmsd_threshold, output_path)


if __name__ == "__main__":
    main()
