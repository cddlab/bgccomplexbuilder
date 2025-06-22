#!/usr/bin/env python3
# %%
import os

import pandas as pd


def make_dataframe(target_dir: str) -> pd.DataFrame:
    """Create a DataFrame from BGC directories.
    Args:
        target_dir (str): The directory containing BGC folders.
    Returns:
        pd.DataFrame: A DataFrame containing the BGC information.
    """
    cols = ["BGC", "proteins", "pdb id", "rmsd"]
    df = pd.DataFrame(columns=cols)
    # find all BGC directories starts with "BGC" in target_dir
    bgc_dirs = [
        d
        for d in os.listdir(target_dir)
        if d.startswith("BGC000") and os.path.isdir(os.path.join(target_dir, d))
    ]

    for dir in bgc_dirs:
        bgcaccession_id = dir.split("_")[0]
        pdb_id = dir.split("_")[1]
        dir_path = os.path.join(target_dir, dir)
        for filename in os.listdir(dir_path):
            if filename.endswith("_af3pae_best.png"):
                filesuffix = filename.split("_af3pae_best.png")[0]
                df.loc[len(df)] = [
                    bgcaccession_id,
                    filesuffix,
                    pdb_id,
                    None,
                ]  # Placeholder for RMSD value

    df.sort_values(by="BGC", inplace=True)
    return df


# %%

target_dir = "/Users/YoshitakaM/Library/CloudStorage/OneDrive-tmd.ac.jp/bgccomplex/positive_hetdimers"
output_sheet = os.path.join(target_dir, "heterocomplexes.xlsx")
df = make_dataframe(target_dir)
df.to_excel(output_sheet, sheet_name="heterocomplexes", index=False)
# %%
