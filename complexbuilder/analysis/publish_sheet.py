#!/usr/bin/env python3

import json
import os

import pandas as pd


def make_heterodataframe(target_dir: str) -> pd.DataFrame:
    """Create a DataFrame from BGC directories.
    Args:
        target_dir (str): The directory containing BGC folders.
    Returns:
        pd.DataFrame: A DataFrame containing the BGC information.
    """
    cols = [
        "BGC",
        "proteins",
        "pdb id",
        "ipSAE",
        "ipSAE_d0chn",
        "ipSAE_d0dom",
        "ipTM_af",
        "ipTM_d0chn",
    ]
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
                metrics_file = os.path.join(dir_path, f"{filesuffix}_ipsae.json")
                with open(metrics_file, "r") as f:
                    metrics = json.load(f)[0]
                df.loc[len(df)] = [
                    bgcaccession_id,
                    filesuffix,
                    pdb_id,
                    metrics["ipSAE"],
                    metrics["ipSAE_d0chn"],
                    metrics["ipSAE_d0dom"],
                    metrics["ipTM_af"],
                    metrics["ipTM_d0chn"],
                ]

    df.sort_values(by="BGC", inplace=True)
    return df


pdbids = [
    "2YJN",
    "6O6E",
    "8K4R",
    "5CZD",
    "1TQY",
    "7THN",
    "4TX3",
    "4TX3",
    "5IG9",
    "8K60",
    "8UC3",
    "6J2U",
    "5DWZ",
    "8GS1",
    "5KP6",
    "8A82",
    "8VSI",
    "8T19",
    "5KP7",
    "6M01",
    "7YN3",
    "6CXT",
    "8HCI",
    "8HK0",
    "6KXD",
    "6M7L",
    "7TCR",
    "6QSP",
    "6QSR",
    "8HCI",
]
for pdbid in pdbids:
    jsonfile = f"/Users/YoshitakaM/Desktop/work/rcsb_pdb_api/rcsb_{pdbid}.json"
    with open(jsonfile, "r") as f:
        rcsbstructdata = json.load(f)

    oligomeric_states = [
        rcsbstructdata["rcsb_struct_symmetry"][i]["oligomeric_state"]
        for i in range(len(rcsbstructdata["rcsb_struct_symmetry"]))
    ]
    print(f"PDB ID {pdbid}: {oligomeric_states}")

target_dir = "/path/to/positive_homomers"
output_sheet = os.path.join(target_dir, "homocomplexes.xlsx")

target_dir = "/path/to/bgccomplex/positive_hetdimers"
output_sheet = os.path.join(target_dir, "heterocomplexes.xlsx")
df = make_heterodataframe(target_dir)
df.to_excel(output_sheet, sheet_name="heterocomplexes", index=False)
