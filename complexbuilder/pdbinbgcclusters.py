#!/usr/bin/env python3
# %%
import csv
import io
import json
from itertools import permutations
from pathlib import Path

import pandas as pd

from complexbuilder.common.parser import sanitised_name


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


def parse_chain_ids(chain_field: str) -> list[str]:
    """
    Parse a chain_id field string that may include multiple values
    enclosed in {} and possibly quoted. Commas inside quotes are kept intact.

    For example:
        '{A,"A,B",A}' -> ["A", "A,B", "A"]

    Args:
        chain_field (str): The input chain_id field string.

    Returns:
        list[str]: List of individual chain identifiers.
    """
    chain_field = chain_field.strip()
    if chain_field.startswith("{") and chain_field.endswith("}"):
        # Remove the surrounding braces.
        chain_field = chain_field[1:-1]
    # Use csv reader to correctly handle commas inside quotes.
    reader = csv.reader(io.StringIO(chain_field), skipinitialspace=True)
    return next(reader)


def generate_linked_pairs(strings: list[str]) -> list[str]:
    """
    Generate all concatenated permutation pairs from a list of strings,
    connected by an underscore.

    For example, given:
        ['AAAA', 'BBBB', 'CCCC']
    It returns:
        ["AAAA_BBBB", "AAAA_CCCC", "BBBB_AAAA", "BBBB_CCCC", "CCCC_AAAA", "CCCC_BBBB"]

    Args:
        strings (list[str]): List of strings.

    Returns:
        list[str]: List of concatenated string pairs.
    """
    return [f"{a}_{b}" for a, b in permutations(strings, 2)]


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


def find_shared_pdb_protein_pairs(df: pd.DataFrame) -> list[tuple]:
    """
    Find protein pairs that share the same PDB ID but have different chain IDs within each BGC.

    This function identifies cases where multiple proteins from the same BGC accession
    are found in the same PDB structure but with different chain identifiers. This indicates
    potential protein-protein interactions that have been experimentally validated and
    captured in protein structure databases.

    Args:
        df (pd.DataFrame): A DataFrame containing at minimum the columns:
            - "mibig_accession": BGC identifier
            - "protein_id": Protein identifier
            - "pdb_id": PDB identifier (possibly with multiple values in braces)
            - "chain_id": Chain identifier (possibly with multiple values in braces)

    Returns:
        list[tuple]: A list of tuples, each containing:
            - accession (str): The BGC accession
            - pdb (str): The shared PDB ID
            - protein_set (set): Set of protein IDs associated with this PDB ID
            - chain_set (set): Set of chain IDs associated with these proteins

    Note:
        The function excludes cases where different proteins share a PDB ID but all have
        identical chain IDs, as these likely represent the same physical entity in the structure.
    """
    results = []
    for accession, group in df.groupby("mibig_accession"):
        # Build a dictionary: { pdb_id : list of (protein_id, chain_id) }
        pdb_to_entries = {}
        for _, row in group.iterrows():
            protein = row["protein_id"]
            pdb_ids = parse_multi_value(row["pdb_id"])
            chain_ids = parse_chain_ids(row["chain_id"])  # same length as pdb_ids
            for i in range(len(pdb_ids)):
                pdb = pdb_ids[i]
                chain = chain_ids[i]
                # Append the tuple (protein, chain) for this pdb id
                pdb_to_entries.setdefault(pdb, []).append((protein, chain))

        # For each pdb id, check if there are different protein_ids
        # and that they don't all have the same chain.
        for pdb, entries in pdb_to_entries.items():
            # Get the set of unique protein_ids and unique chains for this pdb id.
            protein_set = {prot for prot, ch in entries}
            chain_set = {ch for prot, ch in entries}

            if len(protein_set) > 1:
                # Exclude if the chain IDs are identical for all occurrences.
                if len(chain_set) == 1:
                    # Same pdbid_chain for all proteins, so ignore.
                    continue
                else:
                    results.append((accession, pdb, protein_set, chain_set))

    return results


# %%
# blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12_mini1.tsv"
blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12.tsv"

df = pd.read_csv(blast_pdbfile, delimiter="\t")
# dfのidentity列の値が95.0以上の行を取得
df2 = df[df["identity"] >= 95.0].copy()
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
df3 = df2.merge(
    max_chain_count, on=["mibig_accession", "protein_id"], suffixes=("", "_max")
)
# %%


# %%
results = find_shared_pdb_protein_pairs(df2)
seen_proteins = {}

for accession, pdb, proteins, _ in results:
    # Convert the proteins set into a frozenset to use as a hashable key.
    protein_key = frozenset(proteins)
    if accession not in seen_proteins:
        seen_proteins[accession] = set()
    # If we've already seen this set of protein IDs for the given accession, skip.
    if protein_key in seen_proteins[accession]:
        continue
    seen_proteins[accession].add(protein_key)
    if len(list(proteins)) == 2:
        linked_pairs = generate_linked_pairs(list(proteins))
        for pair in linked_pairs:
            # print(
            #     f"BGC: {accession}, Linked Pair: {sanitised_name(pair)}, PDB: {pdb}, Chains: {', '.join(chains)}"
            # )
            print(
                f"mkdir -p {accession}_{str(pdb).upper()}\n"
                f"scp -rp yayoi_out:/data2/moriwaki/BGCcomplex/merged/{accession}/{sanitised_name(pair)}/'*' ./{accession}_{pdb}"
            )

    # %%
