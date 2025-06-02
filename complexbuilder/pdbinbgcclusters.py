#!/usr/bin/env python3
# %%
import json
from pathlib import Path

import pandas as pd

# %%
# blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12_mini1.tsv"
blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12.tsv"

# %%
# TAB-separated data file
df = pd.read_csv(blast_pdbfile, delimiter="\t")
# dfのidentity列の値が95.0以上の行を取得
df2 = df[df["identity"] >= 95.0].copy()


def parse_multi_value(field: str) -> list[str]:
    """
    Parse a string field that may contain multiple values enclosed in {}.
    If the field contains quoted values (e.g. "val1,val2"), it splits accordingly.
    """
    field = field.strip()
    if field.startswith("{") and field.endswith("}"):
        field = field[1:-1]  # remove surrounding braces
        # If there are double quotes, split by '","'
        if '"' in field:
            # Remove any dangling quotes and split
            values = field.split('","')
            return [v.replace('"', "").strip() for v in values]
        else:
            return [x.strip() for x in field.split(",") if x.strip()]
    return [field]


def count_chain_group(chain_group: str) -> int:
    """
    Count number of chain IDs in a comma-separated chain group.
    """
    return len([c.strip() for c in chain_group.split(",") if c.strip()])


def get_max_chain_count(chain_field: str) -> int:
    """
    For a chain_id field that may include multiple groups (separated by commas and braces),
    compute the chain count for each group and return the maximum.
    各行のchain_idについて複数グループがあれば最大値を返す。
    """
    groups = parse_multi_value(chain_field)
    return max(count_chain_group(g) for g in groups)


def is_valid_mibig_accession(accession: str, bgcnumber: int = 2826) -> bool:
    """
    Check if the given mibig_accession value is in the range BGC0000001 to BGC0002826.

    Args:
        accession (str): The accession id string in the format "BGC000XXXX".

    Returns:
        bool: True if the number in accession is between 1 and 2826 (inclusive), False otherwise.
    """
    if not accession.startswith("BGC"):
        return False
    num_part = accession[3:]
    try:
        num = int(num_part)
    except ValueError:
        return False
    return 1 <= num <= bgcnumber


def get_oligomeric_state(pdb_id_field: str):
    """
    Get the oligomeric state from "rcsb_{pdb_Id}.json" file.
    In the json file, the oligomeric state is embedded in the ["rcsb_struct_symmetry"][0]["oligomeric_state"]
    field.
    Args:
        pdb_id_field (str): The PDB ID string. e.g. "{5DYV,7PXO}" or "{8QFU}"
    Returns:
        str: The oligomeric state as a string. If the oligomeric state is not found or invalid, return "Unknown".
        e.g. "2
    """
    # Extract the PDB ID from the field
    pdb_ids = parse_multi_value(pdb_id_field)
    # Define the path to the JSON file
    for pdb_id in pdb_ids:
        json_path = Path(
            f"/Users/YoshitakaM/Desktop/work/rcsb_pdb_api/rcsb_{pdb_id}.json"
        )
        with open(json_path, "r") as f:
            data = json.load(f)

        # Extract the oligomeric state
        try:
            oligomeric_state = data["rcsb_struct_symmetry"][0]["oligomeric_state"]
            return oligomeric_state
        except (KeyError, IndexError):
            return "Unknown"


df2 = df2[df2["mibig_accession"].apply(is_valid_mibig_accession)]

df2["max_chain_count"] = df2["chain_id"].apply(get_max_chain_count)
df2["oligomeric_state"] = df2["pdb_id"].apply(get_oligomeric_state)

# "mibig_accession", "protein_id" ごとに最大の chain_count を求める
max_chain_count = (
    df2.groupby(["mibig_accession", "protein_id"])["max_chain_count"]
    .max()
    .reset_index()
)

# max_chain_countが2以上のものを取得
df2 = df2.merge(
    max_chain_count, on=["mibig_accession", "protein_id"], suffixes=("", "_max")
)

# df3 = df2[df2["max_chain_count_max"] >= 2]
# df3.to_csv("/Users/YoshitakaM/Desktop/pdbinbgcclusters.csv", index=False)
# %%
