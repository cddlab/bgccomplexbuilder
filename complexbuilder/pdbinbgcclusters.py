#!/usr/bin/env python3
# %%
import csv
import io
import json
import os
import re
from itertools import permutations
from pathlib import Path

import pandas as pd

from complexbuilder.common.parser import sanitised_name


def parse_multi_value(field: str) -> list[str]:
    """
    Parse a string field that may contain multiple values enclosed in {}.
    If the field contains quoted values (e.g. "val1,val2"), it splits accordingly.
    e.g. "{A,B,C}" -> ["A", "B", "C"]
    e.g. '{A,"A,B",A}' -> ["A", "A,B", "A"]
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
    If multiple chain groups are present, it returns the maximum count of chains.
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


def _parse_max_homooligomeric_state_to_int(oligomeric_state: str) -> int | None:
    """
    Parse the homooligomeric state string to an integer.
    E.g.
      'Homo 2-mer' -> 2
      'Monomer' -> 1
      'Homo 4-mer' -> 4
      'Hetero 2-mer' -> None
      'Unknown' -> None
      'Homo 6-mer' -> 6
    """
    max_oligomer = None
    if oligomeric_state.startswith("Homo"):
        max_oligomer = int(oligomeric_state.split(" ")[1].split("-")[0])
    elif oligomeric_state == "Monomer":
        max_oligomer = 1
    return max_oligomer


def get_oligomeric_state(pdb_id_field: str, rcsb_data_dir: str):
    """
    Get the oligomeric state from "rcsb_{pdb_Id}.json" file.
    In the json file, the oligomeric state is embedded in the ["rcsb_struct_symmetry"][0]["oligomeric_state"]
    field.
    Args:
        pdb_id_field (str): The PDB ID string. e.g. "{5DYV,7PXO}" or "{8QFU}"
        rcsb_data_dir (str): The directory containing RCSB JSON files.
                             e.g. "/Users/moriwaki/work/rcsb_pdb_api"
    Returns:
        str: The oligomeric state as a string. If the oligomeric state is not found or invalid, return "Unknown".
    """
    # Extract the PDB ID from the field
    pdb_ids = parse_multi_value(pdb_id_field)
    # Define the path to the JSON file
    for pdb_id in pdb_ids:
        json_path = Path(f"{rcsb_data_dir}/rcsb_{pdb_id}.json")
        with open(json_path, "r") as f:
            data = json.load(f)

        # Extract the oligomeric state
        try:
            oligomeric_state = data["rcsb_struct_symmetry"][0]["oligomeric_state"]
            return oligomeric_state
        except (KeyError, IndexError):
            return "Unknown"


def find_pdbid_that_have_different_chain_ids(df: pd.DataFrame) -> list[tuple]:
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


def publish_sheet(
    df: pd.DataFrame,
    target_dir: str = "/Users/YoshitakaM/Desktop/positive_homomers",
    output_file: str = "homocomplexes.xlsx",
) -> None:
    """
    Publish the DataFrame to an Excel file in the specified directory.
    Args:
        df (pd.DataFrame): The DataFrame to be published.
        target_dir (str): The directory where the output file will be saved.
        output_file (str): The name of the output Excel file.
    Returns:
        None
    """
    os.makedirs(target_dir, exist_ok=True)
    output_sheet = os.path.join(target_dir, output_file)
    df.to_excel(output_sheet, sheet_name="homocomplexes", index=False)


def create_proteins_column(row):
    """
    Create a proteins column by concatenating the sanitised protein_id twice with an underscore.
    """
    return f"{sanitised_name(row['protein_id'])}_{sanitised_name(row['protein_id'])}"


def make_homodataframe(target_dir: str) -> pd.DataFrame:
    """Create a DataFrame from BGC directories.
    Args:
        target_dir (str): The directory containing BGC folders.
    Returns:
        pd.DataFrame: A DataFrame containing the BGC information.
    """
    cols = [
        "mibig_accession",
        "proteins",
        "ipSAE",
        "ipSAE_d0chn",
        "ipSAE_d0dom",
        "ipTM_af",
        "ipTM_d0chn",
    ]
    df = pd.DataFrame(columns=cols)

    bgc_dirs = [
        d
        for d in os.listdir(target_dir)
        if d.startswith("BGC000") and os.path.isdir(os.path.join(target_dir, d))
    ]

    for bgcaccession_id in bgc_dirs:
        homocomplex_dirs = [
            d for d in os.listdir(os.path.join(target_dir, bgcaccession_id))
        ]
        for homocomplex_dir in homocomplex_dirs:
            for filename in os.listdir(
                os.path.join(target_dir, bgcaccession_id, homocomplex_dir)
            ):
                if not filename.endswith("_af3pae_best.png"):
                    continue
                if filename.endswith("_af3pae_best.png"):
                    filesuffix = filename.split("_af3pae_best.png")[0]
                    metrics_file = os.path.join(
                        target_dir,
                        bgcaccession_id,
                        homocomplex_dir,
                        f"{filesuffix}_ipsae.json",
                    )
                with open(metrics_file, "r") as f:
                    metrics = json.load(f)[0]
                df.loc[len(df)] = [
                    bgcaccession_id,
                    filesuffix,
                    metrics["ipSAE"],
                    metrics["ipSAE_d0chn"],
                    metrics["ipSAE_d0dom"],
                    metrics["ipTM_af"],
                    metrics["ipTM_d0chn"],
                ]

    df.sort_values(by="mibig_accession", inplace=True)
    return df


def int_to_chain_ids(chain_count: int) -> list[str]:
    """
    Convert an integer chain count to a string of chain IDs.
    For example, 2 -> ["A", "B"], 4 -> ["A", "B", "C", "D"], etc.
    """
    return [chr(65 + i) for i in range(chain_count)]


def transfer_homomer_data(df6: pd.DataFrame, homomer_dir: str) -> None:
    """
    Process homomer data to generate JSON files and shell scripts for transfer.

    Args:
        df6 (pd.DataFrame): DataFrame containing homomer data.
        homomer_dir (str): Directory where the processed files will be stored.
    """
    gpuH100_transfer_dir = os.path.join(homomer_dir, "gpuH100_transfer")
    os.makedirs(gpuH100_transfer_dir, exist_ok=True)
    gpu4090_transfer_dir = os.path.join(homomer_dir, "gpu4090_transfer")
    os.makedirs(gpu4090_transfer_dir, exist_ok=True)

    returnfile_f = os.path.join(homomer_dir, "return_gpu4090.sh")
    returnfile_t = os.path.join(homomer_dir, "return_gpuH100.sh")

    return_handle_f = open(returnfile_f, "w")
    return_handle_f.write("#!/bin/bash\n")
    return_handle_t = open(returnfile_t, "w")
    return_handle_t.write("#!/bin/bash\n")

    for row in df6.itertuples(index=False, name="Pandas"):
        if (
            isinstance(row.parsed_oligomeric_state, float)
            and row.parsed_oligomeric_state > 2.0
        ):
            filename = os.path.join(
                homomer_dir, f"{row.mibig_accession}/{row.proteins}_data.json"
            )
            if not os.path.exists(filename):
                continue

            with open(filename) as f:
                tmp = json.load(f)
                if tmp["sequences"][0]["protein"]["id"] == ["A", "B"]:
                    chain_count = int(row.parsed_oligomeric_state)
                    tmp["sequences"][0]["protein"]["id"] = int_to_chain_ids(chain_count)
                    newname = f"{sanitised_name(row.protein_id)}_{chain_count}mer"
                    tmp["name"] = newname

            seq_length = len(tmp["sequences"][0]["protein"]["sequence"]) * chain_count

            if seq_length < 2001:
                outfile = os.path.join(gpu4090_transfer_dir, f"{newname}.json")
                with open(outfile, "w") as f:
                    json.dump(tmp, f, indent=2, ensure_ascii=False)
                return_handle_f.write(f"mkdir -p ../{row.mibig_accession}\n")
                return_handle_f.write(
                    f"mv {sanitised_name(row.protein_id)}_{chain_count}mer {sanitised_name(row.protein_id)}_{chain_count}mer.json ../{row.mibig_accession}\n"
                )
            elif seq_length < 4001:
                outfile = os.path.join(gpuH100_transfer_dir, f"{newname}.json")
                with open(outfile, "w") as f:
                    json.dump(tmp, f, indent=2, ensure_ascii=False)
                return_handle_t.write(f"mkdir -p ../{row.mibig_accession}\n")
                return_handle_t.write(
                    f"mv {sanitised_name(row.protein_id)}_{chain_count}mer {sanitised_name(row.protein_id)}_{chain_count}mer.json ../{row.mibig_accession}\n"
                )
            else:
                print(
                    f"Skipping {newname} of {row.mibig_accession} as it exceeds 4000 sequence length for transfer."
                )
                continue

    return_handle_f.close()
    return_handle_t.close()


# %%
def _convert_to_repeated_protein_id(protein_id: str) -> str:
    """
    Convert a protein ID to a repeated name.
    e.g. "aam94792.1_4mer" -> "aam94792.1_aam94792.1"
         "aam94793.1_4mer" -> "aam94793.1_aam94793.1"
         "aam70337.1_10mer" -> "aam70337.1_aam70337.1"
         "wp_003722052.1_3mer" -> "wp_003722052.1_wp_003722052.1"
    """
    # remove the suffix "_4mer", "_10mer", etc.
    protein = re.sub(r"(_\d+mer)$", "", protein_id)
    return f"{protein}_{protein}"


def homomer_add_hit_to_df(jsondata: dict) -> pd.DataFrame:
    """
    Convert nested JSON data with BGC entries to a pandas DataFrame.
    Args:
        json_data (dict): Nested JSON data with BGC IDs as top-level keys

    Returns:
        pd.DataFrame: DataFrame with columns: mibig_accession, proteins, ipSAE, ipTM
    """
    rows = []
    for bgc_id, proteins_data in jsondata.items():
        if not proteins_data:
            continue
        for protein_id, metrics in proteins_data.items():
            row = {
                "mibig_accession": bgc_id,
                "proteins": _convert_to_repeated_protein_id(protein_id),
                "ipSAE_homomer": metrics.get("ipSAE"),
                "ipTM_homomer": metrics.get("ipTM"),
            }
            rows.append(row)
    df = pd.DataFrame(rows)
    return df


# %%
blast_pdbfile = "/Users/YoshitakaM/Desktop/blast_pdb24.12.tsv"
df = pd.read_csv(blast_pdbfile, delimiter="\t")
# %%
df2 = df[df["identity"] >= 95.0].copy()
df2 = df2[df2["mibig_accession"].apply(is_valid_mibig_accession)]
df2["oligomeric_state"] = df2["pdb_id"].apply(
    get_oligomeric_state, rcsb_data_dir="/Users/YoshitakaM/Desktop/work/rcsb_pdb_api"
)
df2["parsed_oligomeric_state"] = df2["oligomeric_state"].apply(
    _parse_max_homooligomeric_state_to_int
)
df3 = df2.copy()
max_mask = df3.groupby(["mibig_accession", "protein_id"])[
    "parsed_oligomeric_state"
].transform(lambda x: x == x.max() if x.max() is not None else x.isna())
df3 = df3[max_mask]
df3 = df3.drop_duplicates(["mibig_accession", "protein_id"])
df4 = df3[df3["parsed_oligomeric_state"] >= 2.0].copy()
df4.loc[:, "proteins"] = df4.apply(create_proteins_column, axis=1)
# %%
target_dir = "/Users/YoshitakaM/Desktop/positive_homomers"
output_sheet = os.path.join(target_dir, "homocomplexes23.xlsx")
df5 = make_homodataframe(target_dir)
# %%
df6 = pd.merge(
    df4,
    df5,
    how="left",
    left_on=["mibig_accession", "proteins"],
    right_on=["mibig_accession", "proteins"],
    suffixes=("", "_y"),
)
publish_sheet(df6, target_dir=target_dir, output_file="homocomplexes23.xlsx")


# %%
additional_homomer_file = "/Users/YoshitakaM/Desktop/positive_homomers/homomer_additional/homomer_additional_hitcomplexes.json"
with open(additional_homomer_file, "r") as f:
    additional_homomer_data = json.load(f)
df7 = homomer_add_hit_to_df(additional_homomer_data)
df8 = pd.merge(
    df6,
    df7,
    how="left",
    left_on=["mibig_accession", "proteins"],
    right_on=["mibig_accession", "proteins"],
    suffixes=("", "_y"),
)
publish_sheet(df8, target_dir=target_dir, output_file="homocomplexes3.xlsx")
