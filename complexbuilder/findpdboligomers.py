#!/usr/bin/env python3
# %%

import json
import os
import time

import requests


def fetch_pdb_assembly(pdb_id, outputdir: str, assembly_id="1") -> None:
    """Fetch assembly information for a given PDB ID from the RCSB PDB API.
    Args:
        pdb_id (str): The PDB ID to fetch assembly information for.
        outputdir (str): The directory to save the output JSON file.
        assembly_id (str): The assembly ID to fetch. Default is "1".
    Returns:
        dict: A dictionary containing the formatted assembly information.
    """
    url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
    else:
        print(f"Failed to fetch data for {pdb_id}")
    # format the data and write to file
    filename = f"{outputdir}/rcsb_{pdb_id.upper()}.json"
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def get_currentallpdbids() -> list[str]:
    """Retrieve the current PDB IDs from the RCSB PDB holdings API.
    url: https://data.rcsb.org/rest/v1/holdings/current/entry_ids
    Returns:
        list[str]: A list of current PDB IDs.
        The IDs are in the format of strings, e.g. ["100D","1ATR","1BU3","1CVB",...]
    """
    url = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
    response = requests.get(url)
    if response.status_code == 200:
        # display the number of current PDB IDs
        print(f"Number of current PDB IDs: {len(response.json())}")
        return response.json()

    else:
        print("Failed to fetch current PDB IDs")
        return []


# %%
pdb_id_list = get_currentallpdbids()
# pdb_id_list = ["100D", "1ATR", "1BU3", "1CVB", "1DWD", "1EW9"]
outputdir = "/Users/YoshitakaM/Desktop/work/rcsb_pdb_api_all"
sleep_duration = 0.3  # seconds
os.makedirs(outputdir, exist_ok=True)
for pdb_id in pdb_id_list:
    if os.path.exists(f"{outputdir}/rcsb_{pdb_id.upper()}.json"):
        print(f"File for {pdb_id} already exists, skipping...")
        continue
    else:
        fetch_pdb_assembly(pdb_id, outputdir, assembly_id="1")
        print(f"Fetched assembly for {pdb_id}")
        time.sleep(sleep_duration)  # Sleep to avoid hitting the API too hard

# %%
