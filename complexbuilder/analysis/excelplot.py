# %%
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger
from matplotlib import rcParams

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
# %%
excelfile = Path(
    "/Users/YoshitakaM/Library/CloudStorage/OneDrive-TheUniversityofTokyo/bgccomplex/homocomplexes3.xlsx"
)
df = pd.read_excel(excelfile)
ipsae_values = df.iloc[:, 13]  # N列
iptm_values = df.iloc[:, 16]  # Q列

valid_data = ~(ipsae_values.isna() | iptm_values.isna())
ipsae_valid = ipsae_values[valid_data]
iptm_valid = iptm_values[valid_data]
fig, ax = plt.subplots(1, 3, figsize=(18, 6), dpi=300)

ax[0].scatter(
    iptm_valid, ipsae_valid, alpha=0.7, s=30, c="#0072BC", edgecolors="w", linewidth=0.5
)
ax[0].set_title("ipSAE vs ipTM for Homooligomers in PDB")
ax[0].set_xlabel("Chain pair ipTM")
ax[0].set_ylabel("ipSAE")
ax[0].grid(True, linestyle="--", alpha=0.7)

ax[0].set_xlim(0, 1)
ax[0].set_ylim(-0.02, 1)

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
        alpha=0.1,
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
    label="2-mer",
)
ax[1].scatter(
    iptm_af_homo_valid,
    ipsae_homo_valid,
    alpha=0.7,
    s=30,
    c="#D25A45",
    label="oligomeric number of \nGlobal Stoichiometry \ndisplayed in PDB",
    edgecolors="w",
    marker="s",
    linewidth=0.5,
)

ax[1].set_xlim(0, 1)
ax[1].set_ylim(-0.02, 1)
hans, labs = ax[1].get_legend_handles_labels()
ax[1].legend(handles=hans, labels=labs, fontsize=12)
ax[1].set_title("Change in ipSAE and ipTM metrics")
ax[1].set_xlabel("Chain pair ipTM")
ax[1].set_ylabel("ipSAE")

excelfile = Path(
    "/Users/YoshitakaM/Library/CloudStorage/OneDrive-TheUniversityofTokyo/bgccomplex/positive_hetdimers/heterocomplexes_colored.xlsx"
)
df2 = pd.read_excel(excelfile)
hetero_ipsae_values = df2.iloc[:, 3]  # column D
hetero_iptm_values = df2.iloc[:, 6]  # column G
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
ax[2].set_title("Heterocomplexes in PDB")
ax[2].set_xlabel("Chain pair ipTM", fontsize=16)
ax[2].set_ylabel("ipSAE", fontsize=16)
plt.tight_layout()
plt.savefig("fig2.svg", format="svg")
plt.show()
plt.clf()
plt.close()
# %%
