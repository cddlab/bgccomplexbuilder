#!/usr/bin/env python3
# %%
from pathlib import Path

import pandas as pd

# %%
# blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12_mini1.tsv"
blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12.tsv"

# %%
# TAB-separated data file
df = pd.read_csv(blast_pdbfile, delimiter="\t")
# dfのidentity列の値が95.0以上の行を取得
df_ident95 = df[df["identity"] >= 95.0]

# %%
# df_ident95の1行目の"pdb_id"列の値を取得
# '{5DYV,7PXO}'をリストに変換
df_ident95.iloc[:10]["pdb_id"].str.strip("{}").str.split(", ")

# %%
