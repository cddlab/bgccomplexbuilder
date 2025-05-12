#!/usr/bin/env python3
# %%
import json
import os

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = [
    "Hiragino Maru Gothic Pro",
    "Yu Gothic",
    "Meiryo",
    "Takao",
    "IPAexGothic",
    "IPAPGothic",
    "VL PGothic",
    "Noto Sans CJK JP",
]

metricsfile = (
    "/Users/YoshitakaM/Desktop/BGC_heteromer/BGC0000028/af3/complexmetrics.json"
)

with open(metricsfile, "r") as f:
    data = json.load(f)

ipSAE_dict = {list(d.keys())[0]: list(d.values())[0] for d in data["ipSAE"]}
ipTM_dict = {list(d.keys())[0]: list(d.values())[0] for d in data["ipTM"]}

# 共通するキーのみを対象に散布図を作成
common_keys = set(ipSAE_dict.keys()) & set(ipTM_dict.keys())

# 散布図用のデータ
x_vals = [ipSAE_dict[k] for k in common_keys]
y_vals = [ipTM_dict[k] for k in common_keys]
# labels = list(common_keys)

# プロット
fig, ax = plt.subplots(dpi=300, figsize=(5, 5))
ax.set_title("Scatter Plot of ipSAE vs ipTM")
ax.set_xlabel("ipSAE")
ax.set_xlim(0, 1)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_ylim(0, 1)
ax.set_ylabel("ipTM")
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.scatter(x_vals, y_vals, s=5, c="blue", alpha=0.5)

plt.tight_layout()
plt.show()

# %%
