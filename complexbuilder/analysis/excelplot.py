import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = [
    "Arial",
    "Meiryo",
    "Takao",
]
rcParams["font.size"] = 20
rcParams["axes.labelsize"] = 16
rcParams["axes.titlesize"] = 16
rcParams["xtick.labelsize"] = 16
rcParams["ytick.labelsize"] = 16
rcParams["axes.grid"] = True
rcParams["grid.linestyle"] = "--"
rcParams["grid.linewidth"] = 0.5
rcParams["grid.alpha"] = 0.7

rcParams["axes.edgecolor"] = "black"
rcParams["axes.linewidth"] = 0.8
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"
rcParams["svg.fonttype"] = "none"


def add_second_axis(
    ax,
    xdata,
    ydata,
    color="#0072BC",
    bins=30,
    hist_alpha=0.35,
    show_hist=True,
):
    """Add marginal histograms for x and y around a joint plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Main scatter Axes.
    xdata, ydata : pandas.Series or array-like
        Data for the two axes.
    color : str
        Base color for both KDE and histogram.
    bins : int
        Number of histogram bins (if shown).
    hist_alpha : float
        Alpha value for histogram fill.
    show_hist : bool
        Whether to draw histograms.
    """
    divider = make_axes_locatable(ax)
    ax_kdex = divider.append_axes("top", size="12%", pad=0.02, sharex=ax)
    ax_kdey = divider.append_axes("right", size="12%", pad=0.02, sharey=ax)

    # Clean data
    x = pd.Series(xdata).dropna()
    y = pd.Series(ydata).dropna()

    # Histograms
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

    # Aesthetics: remove ticks/labels but keep transparent background
    for a in (ax_kdex, ax_kdey):
        a.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        a.set_facecolor("none")
        a.grid(False)
        a.set_xlabel("")
        a.set_ylabel("")
        for spine in a.spines.values():
            spine.set_visible(False)


THRESHOLD = 0.55
excelfile = Path(
    "/Users/YoshitakaM/Library/CloudStorage/OneDrive-TheUniversityofTokyo/bgccomplex/sup_data/homocomplexes.xlsx"
)
df = pd.read_excel(excelfile)
ipsae_values = df.iloc[:, 13]  # N column
iptm_values = df.iloc[:, 16]  # Q column
ipsae_homo = df.iloc[:, 18]  # S column
iptm_af_homo = df.iloc[:, 19]  # T column
copy_number = df.iloc[:, 11]  # L column

valid_data = (
    ~(ipsae_values.isna() | iptm_values.isna())  # ipsae and iptm are both present
    & (copy_number == 2)  # and, copy number is 2.
)

ipsae_valid = ipsae_values[valid_data]
iptm_valid = iptm_values[valid_data]

print(f"iptm_valid >= 0.0: {iptm_valid[iptm_valid >= 0].count()}")
print(f"iptm_valid < 0.6: {iptm_valid[iptm_valid < 0.6].count()}")
print(f"iptm_valid >= 0.6: {iptm_valid[iptm_valid >= 0.6].count()}")
print(f"iptm_valid >= 0.7: {iptm_valid[iptm_valid >= 0.7].count()}")
print(f"iptm_valid >= 0.8: {iptm_valid[iptm_valid >= 0.8].count()}")

fig, ax = plt.subplots(1, 3, figsize=(18, 6), dpi=300)

ax[0].scatter(
    iptm_valid, ipsae_valid, alpha=0.7, s=30, c="#0072BC", edgecolors="w", linewidth=0.5
)
ax[0].set_xlabel("Chain pair ipTM")
ax[0].set_ylabel("ipSAE")
# ax[0].axvline(THRESHOLD, c="gray", linestyle="--", linewidth=1.0)
ax[0].grid(True, linestyle="--", alpha=0.7)

ax[0].set_xlim(0, 1)
ax[0].set_ylim(-0.02, 1)
add_second_axis(ax[0], iptm_valid, ipsae_valid, color="#0072BC")

ipsae_homo = df.iloc[:, 18]  # S列
iptm_af_homo = df.iloc[:, 19]  # T列
valid_data = ~(ipsae_homo.isna() | iptm_af_homo.isna())
ipsae_valid = ipsae_values[valid_data]
iptm_af_valid = iptm_values[valid_data]
ipsae_homo_valid = ipsae_homo[valid_data]
iptm_af_homo_valid = iptm_af_homo[valid_data]

for i in range(len(ipsae_valid)):
    x1, y1 = iptm_af_valid.iloc[i], ipsae_valid.iloc[i]  # 始点 (ipTM_af, ipSAE)
    x2, y2 = (
        iptm_af_homo_valid.iloc[i],
        ipsae_homo_valid.iloc[i],
    )

    ax[1].arrow(
        x1,
        y1,
        x2 - x1,
        y2 - y1,
        head_width=0.01,
        head_length=0.02,
        fc="gray",
        ec="gray",
        alpha=0.4,
        length_includes_head=True,
    )

ax[1].scatter(
    iptm_af_valid,
    ipsae_valid,
    alpha=0.7,
    s=30,
    c="#0072BC",
    edgecolors="w",
    linewidth=0.5,
    label="dimer",
)
ax[1].scatter(
    iptm_af_homo_valid,
    ipsae_homo_valid,
    alpha=0.7,
    s=30,
    c="#D25A45",
    label="biological assembly",
    edgecolors="w",
    marker="s",
    linewidth=0.5,
)

ax[1].set_xlim(0, 1)
ax[1].set_ylim(-0.02, 1)
# ax[1].axvline(THRESHOLD, c="gray", linestyle="--", linewidth=1.0)
hans, labs = ax[1].get_legend_handles_labels()
ax[1].legend(handles=hans, labels=labs, fontsize=12)
ax[1].set_xlabel("Chain pair ipTM")
ax[1].set_ylabel("ipSAE")
add_second_axis(ax[1], iptm_af_valid, ipsae_valid, color="#0072BC")
add_second_axis(ax[1], iptm_af_homo_valid, ipsae_homo_valid, color="#D25A45")
excelfile = Path(
    "/Users/YoshitakaM/Library/CloudStorage/OneDrive-TheUniversityofTokyo/bgccomplex/sup_data/positive_hetdimers/heterocomplexes_colored.xlsx"
)
df2 = pd.read_excel(excelfile)
hetero_ipsae_values = df2.iloc[:, 3]  # column D
hetero_iptm_values = df2.iloc[:, 5]  # column F
valid_data = ~(hetero_ipsae_values.isna() | hetero_iptm_values.isna())
hetero_ipsae_valid = hetero_ipsae_values[valid_data]
hetero_iptm_valid = hetero_iptm_values[valid_data]
ax[2].scatter(
    hetero_iptm_valid,
    hetero_ipsae_valid,
    alpha=0.7,
    s=30,
    c="#32C189",
    label="Heterocomplexes in PDB",
    edgecolors="w",
    linewidth=0.5,
)
ax[2].set_xlim(0, 1)
ax[2].set_ylim(-0.02, 1)
ax[2].set_xlabel("Chain pair ipTM", fontsize=16)
ax[2].set_ylabel("ipSAE", fontsize=16)
add_second_axis(ax[2], hetero_iptm_valid, hetero_ipsae_valid, color="#32C189")

plt.tight_layout()
plt.savefig("fig2.svg", format="svg")
plt.show()
plt.clf()
plt.close()
