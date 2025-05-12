#!/usr/bin/env python3
# %%
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager, rcParams

# Font directory
font_dir = "/home/moriwaki/fonts"
for font_path in font_manager.findSystemFonts(fontpaths=[font_dir]):
    font_manager.fontManager.addfont(font_path)
rcParams["font.family"] = "Arial"
# %%

# BGCディレクトリがある親ディレクトリ
parent_dir = Path("/data2/moriwaki/BGCcomplex/merged")
# 各BGCディレクトリ内のaf3/complexmetrics.jsonのパスを収集
metrics_paths = []
for entry in parent_dir.iterdir():
    if entry.is_dir() and entry.name.startswith("BGC000"):
        json_path = entry / "complexmetrics.json"
        if json_path.exists():
            metrics_paths.append(json_path)

data = {"ipSAE": [], "ipTM": []}

for json_path in metrics_paths:
    with open(json_path, "r") as f:
        current = json.load(f)
    # check if the keys in each section are unique
    for section in ["ipSAE", "ipTM"]:
        seen = set()
        for entry in current.get(section, []):
            # Each entry is assumed to be a dict with one key
            key = list(entry.keys())[0]
            if key in seen:
                sys.stderr.write(
                    f"Warning: Duplicate key '{key}' found in {json_path} for section '{section}'.\n"
                )
            seen.add(key)
    # 全ファイル分のデータを結合
    data["ipSAE"].extend(current.get("ipSAE", []))
    data["ipTM"].extend(current.get("ipTM", []))

ipSAE_dict = {list(d.keys())[0]: list(d.values())[0] for d in data["ipSAE"]}
ipTM_dict = {list(d.keys())[0]: list(d.values())[0] for d in data["ipTM"]}
# %%
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
ax.set_ylim(0.6, 1)
ax.set_ylabel("ipTM")
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.scatter(x_vals, y_vals, s=0.5, c="blue", alpha=0.2)

plt.tight_layout()
plt.savefig("test.eps", format="eps", bbox_inches="tight")
# 487828 points

# %%

# %%
