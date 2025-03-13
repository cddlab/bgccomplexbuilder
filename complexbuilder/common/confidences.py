import json
import math
import os
import sys

import numpy as np

pae_file_path = sys.argv[1]
pdb_path = sys.argv[2]
pae_cutoff = float(sys.argv[3])
dist_cutoff = float(sys.argv[4])
pae_string = str(int(pae_cutoff))
if pae_cutoff < 10:
    pae_string = "0" + pae_string
dist_string = str(int(dist_cutoff))
if dist_cutoff < 10:
    dist_string = "0" + dist_string


if os.path.splitext(pdb_path)[1] == ".pdb":
    pdb_stem = os.path.splitext(pdb_path)[0]
    path_stem = f"{pdb_stem}_{pae_string}_{dist_string}"
    af2 = True
    af3 = False
    boltz1 = False
    cif = False
elif (
    os.path.splitext(pdb_path)[1] == ".cif"
    and os.path.splitext(pae_file_path)[1] == ".json"
):
    pdb_stem = os.path.splitext(pdb_path)[0]
    path_stem = f"{pdb_stem}_{pae_string}_{dist_string}"
    af2 = False
    af3 = True
    boltz1 = False
    cif = True
elif (
    os.path.splitext(pdb_path)[1] == ".cif"
    and os.path.splitext(pae_file_path)[1] == ".npz"
):
    pdb_stem = os.path.splitext(pdb_path)[0]
    path_stem = f"{pdb_stem}_{pae_string}_{dist_string}"
    af2 = False
    af3 = False
    boltz1 = True
    cif = True
else:
    print("Wrong PDB or PAE file type ", pdb_path)
    sys.exit()

file_path = path_stem + ".txt"
file2_path = path_stem + "_byres.txt"
pml_path = path_stem + ".pml"
OUT = open(file_path, "w")
PML = open(pml_path, "w")
OUT2 = open(file2_path, "w")


def ptm_func(x, d0):
    return 1.0 / (1 + (x / d0) ** 2.0)


ptm_func_vec = np.vectorize(ptm_func)  # vector version


def calc_d0(L: int | float) -> float:
    r"""
    Calculate the d0 value for the number of unique residues in chain B
    that have $PAE_{ij} \lt cutoff$ given the identity of the aligned residue $i$.
    Eq. (15)
    \begin{array}{ll}
        d_0=1.24 \sqrt[3]{L_{PAE<\text { cutoff }}-15}-1.8 & L \geq 27 \\
        d_0=1 & L<27
    \end{array}
    Args:
        L: The number of unique residues in chain B that have PAEij < cutoff
        given the identity of the aligned residue i.
    Returns:
        The d0 value.
    See also:
      "Scoring function for automated assessment of protein structure template quality"
       Zhang and Skolnick, DOI:10.1002/prot.20264
    """
    L = float(L)
    if L < 27:
        return 1.0
    return 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8


def calc_d0_array(L: int | float | np.ndarray) -> np.ndarray:
    """
    Convert L to a NumPy array if it isn't already one
    (enables flexibility in input types)
    Ensure all values of L are at least 19.0
    Calculate d0 using the vectorized operation
    """
    L = np.array(L, dtype=float)
    L = np.maximum(L, 26.523)
    return 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8


def parse_pdb_atom_line(line: str) -> dict:
    """
    Define the parse_atom_line function for PDB lines (by column)
    parsed_line = parse_atom_line(line)
    Example:
        line = "ATOM    123  CA  ALA A  15"
               "  11.111  22.222  33.333  1.00 20.00           C"
        parsed_line = parse_atom_line(line)
        parsed_line == {
            "atom_num": 123,
            "atom_name": "CA",
            "residue_name": "ALA",
            "chain_id": "A",
            "residue_seq_num": 15,
            "x": 11.111,
            "y": 22.222,
            "z": 33.333,
        }
    """
    atom_num = line[6:11].strip()
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    chain_id = line[21].strip()
    residue_seq_num = line[22:26].strip()
    x = line[30:38].strip()
    y = line[38:46].strip()
    z = line[46:54].strip()

    # Convert string numbers to integers or floats as appropriate
    atom_num = int(atom_num)
    residue_seq_num = int(residue_seq_num)
    x = float(x)
    y = float(y)
    z = float(z)

    return {
        "atom_num": atom_num,
        "atom_name": atom_name,
        "residue_name": residue_name,
        "chain_id": chain_id,
        "residue_seq_num": residue_seq_num,
        "x": x,
        "y": y,
        "z": z,
    }


def parse_cif_atom_line(line: str, fielddict) -> dict | None:
    """
    for parsing AF3 and Boltz1 mmCIF files
    ligands do not have residue numbers but modified residues do.
    Return "None" for ligand.
    # AF3 mmcif lines
    # 0      1   2   3     4  5  6 7  8  9  10      11     12      13   14    15 16 17
    # ATOM   1294 N   N     . ARG A 1 159 ? 5.141   -14.096 10.526  1.00 95.62 159 A 1
    # ATOM   1295 C   CA    . ARG A 1 159 ? 4.186   -13.376 11.366  1.00 96.27 159 A 1
    # ATOM   1296 C   C     . ARG A 1 159 ? 2.976   -14.235 11.697  1.00 96.42 159 A 1
    # ATOM   1297 O   O     . ARG A 1 159 ? 2.654   -15.174 10.969  1.00 95.46 159 A 1
    # ...
    # HETATM 1305 N   N     . TPO A 1 160 ? 2.328   -13.853 12.742  1.00 96.42 160 A 1
    # HETATM 1306 C   CA    . TPO A 1 160 ? 1.081   -14.560 13.218  1.00 96.78 160 A 1
    # HETATM 1307 C   C     . TPO A 1 160 ? -2.115  -11.668 12.263  1.00 96.19 160 A 1
    # HETATM 1308 O   O     . TPO A 1 160 ? -1.790  -11.556 11.113  1.00 95.75 160 A 1
    # ...
    # HETATM 2608 P   PG    . ATP C 3 .   ? -6.858  4.182   10.275  1.00 84.94 1   C 1
    # HETATM 2609 O   O1G   . ATP C 3 .   ? -6.178  5.238   11.074  1.00 75.56 1   C 1
    # HETATM 2610 O   O2G   . ATP C 3 .   ? -5.889  3.166   9.748   1.00 75.15 1   C 1
    # ...
    # HETATM 2639 MG  MG    . MG  D 4 .   ? -7.262  2.709   4.825   1.00 91.47 1   D 1
    # HETATM 2640 MG  MG    . MG  E 5 .   ? -4.994  2.251   8.755   1.00 85.96 1   E 1
    """

    linelist = line.split()
    atom_num = linelist[fielddict["id"]]
    atom_name = linelist[fielddict["label_atom_id"]]
    residue_name = linelist[fielddict["label_comp_id"]]
    chain_id = linelist[fielddict["label_asym_id"]]
    residue_seq_num = linelist[fielddict["label_seq_id"]]
    x = linelist[fielddict["Cartn_x"]]
    y = linelist[fielddict["Cartn_y"]]
    z = linelist[fielddict["Cartn_z"]]

    if residue_seq_num == ".":
        return None  # ligand

    # Convert string numbers to integers or floats as appropriate
    atom_num = int(atom_num)
    residue_seq_num = int(residue_seq_num)
    x = float(x)
    y = float(y)
    z = float(z)

    return {
        "atom_num": atom_num,
        "atom_name": atom_name,
        "residue_name": residue_name,
        "chain_id": chain_id,
        "residue_seq_num": residue_seq_num,
        "x": x,
        "y": y,
        "z": z,
    }


def contiguous_ranges(numbers: set[int]) -> str:
    """
    Function for printing out residue numbers in PyMOL scripts
    Args:
        numbers: A list of residue numbers
    Returns:
        A string with residue numbers in contiguous ranges
    Example:

    """
    if len(numbers) == 0:
        return ""

    sorted_numbers = sorted(numbers)
    start = sorted_numbers[0]
    end = start
    ranges = []  # List to store ranges

    def format_range(start, end):
        if start == end:
            return f"{start}"
        else:
            return f"{start}-{end}"

    for number in sorted_numbers[1:]:
        if number == end + 1:
            end = number
        else:
            ranges.append(format_range(start, end))
            start = end = number

    # Append the last range after the loop
    ranges.append(format_range(start, end))

    # Join all ranges with a plus sign and print the result
    string = "+".join(ranges)
    return string


# Load residues from AlphaFold PDB or mmCIF file into lists; each residue is a dictionary
# Read PDB file to get CA coordinates, chainids, and residue numbers
# Convert to np arrays, and calculate distances
residues = []
cb_residues = []
chains = []
atomsitefield_num = 0
atomsitefield_dict = (
    {}
)  # contains order of atom_site fields in mmCIF files; handles any mmCIF field order

# For af3 and boltz1: need mask to identify CA atom tokens in plddt vector and pae matrix;
# Skip ligand atom tokens and non-CA-atom tokens in PTMs (those not in residue_set)
token_mask = list()
residue_set = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

with open(pdb_path, "r") as PDB:
    for line in PDB:
        if line.startswith("_atom_site."):
            line = line.strip()
            (atomsite, fieldname) = line.split(".")
            atomsitefield_dict[fieldname] = atomsitefield_num
            atomsitefield_num += 1

        if line.startswith("ATOM") or line.startswith("HETATM"):
            if cif:
                atom = parse_cif_atom_line(line, atomsitefield_dict)
            else:
                atom = parse_pdb_atom_line(line)

            if atom is None:  # ligand atom
                token_mask.append(0)
                continue

            if atom["atom_name"] == "CA":
                token_mask.append(1)
                residues.append(
                    {
                        "atom_num": atom["atom_num"],
                        "coor": np.array([atom["x"], atom["y"], atom["z"]]),
                        "res": atom["residue_name"],
                        "chainid": atom["chain_id"],
                        "resnum": atom["residue_seq_num"],
                        "residue": f"{atom['residue_name']:3}   "
                        f"{atom['chain_id']:3} "
                        f"{atom['residue_seq_num']:4}",
                    }
                )
                chains.append(atom["chain_id"])

            if atom["atom_name"] == "CB" or (
                atom["residue_name"] == "GLY" and atom["atom_name"] == "CA"
            ):
                cb_residues.append(
                    {
                        "atom_num": atom["atom_num"],
                        "coor": np.array([atom["x"], atom["y"], atom["z"]]),
                        "res": atom["residue_name"],
                        "chainid": atom["chain_id"],
                        "resnum": atom["residue_seq_num"],
                        "residue": f"{atom['residue_name']:3}   "
                        f"{atom['chain_id']:3} "
                        f"{atom['residue_seq_num']:4}",
                    }
                )

            # add nucleic acids and non-CA atoms in PTM residues to tokens (as 0),
            # whether labeled as "HETATM" (af3) or as "ATOM" (boltz1)
            if atom["atom_name"] != "CA" and atom["residue_name"] not in residue_set:
                token_mask.append(0)

# Convert structure information to numpy arrays
numres = len(residues)
CA_atom_num = np.array(
    [res["atom_num"] - 1 for res in residues]
)  # for AF3 atom indexing from 0
CB_atom_num = np.array(
    [res["atom_num"] - 1 for res in cb_residues]
)  # for AF3 atom indexing from 0
coordinates = np.array([res["coor"] for res in cb_residues])
chains = np.array(chains)
unique_chains = np.unique(chains)
token_array = np.array(token_mask)
ntokens = np.sum(token_array)

# Calculate distance matrix using NumPy broadcasting
distances = np.sqrt(
    ((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]) ** 2).sum(axis=2)
)
