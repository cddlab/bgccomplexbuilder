import json
import math
import os
import sys

import numpy as np

pae_file_path = "tests/testfiles/testAB_confidences.json"
pdb_path = "tests/testfiles/testAB_model.cif"
pae_cutoff = 10.0
dist_cutoff = 10.0
pae_string = "10"
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


# Load residues from AlphaFold PDB or mmCIF file into lists;
# each residue is a dictionary
# Read PDB file to get CA coordinates, chainids, and residue numbers
# Convert to np arrays, and calculate distances
residues = []
cb_residues = []
chains = []
atomsitefield_num = 0
atomsitefield_dict = {}  # contains order of atom_site fields in mmCIF files; handles any mmCIF field order

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
# [   1    9   16  23   30   37   44   51   58   65   72   79   86   93]
CB_atom_num = np.array(
    [res["atom_num"] - 1 for res in cb_residues]
)  # for AF3 atom indexing from 0
# [4  12  19  26  33  40  47  54  61  68  75  82  89  96]
coordinates = np.array([res["coor"] for res in cb_residues])

chains = np.array(chains)
unique_chains = np.unique(chains)  # ["A", "B"]
token_array = np.array(token_mask)
ntokens = np.sum(token_array)

# Calculate distance matrix using NumPy broadcasting
distances = np.sqrt(
    ((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]) ** 2).sum(axis=2)
)


if af3:
    if os.path.exists(pae_file_path):
        with open(pae_file_path, "r") as file:
            data = json.load(file)
    else:
        print("AF3 PAE file does not exist: ", pae_file_path)
        sys.exit()

    atom_plddts = np.array(data["atom_plddts"])
    plddt = atom_plddts[CA_atom_num]  # pull out residue plddts from Calpha atoms
    cb_plddt = atom_plddts[
        CB_atom_num
    ]  # pull out residue plddts from Cbeta atoms for pDockQ

    # Get pairwise residue PAE matrix by identifying one token per protein residue.
    # Modified residues have separate tokens for each atom,
    # so need to pull out Calpha atom as token Skip ligands
    if "pae" in data:
        pae_matrix_af3 = np.array(data["pae"])
    else:
        print("no PAE data in AF3 json file; quitting")
        sys.exit()

    # Set pae_matrix for AF3 from subset of full PAE matrix from json file
    token_array = np.array(token_mask)
    pae_matrix = pae_matrix_af3[
        np.ix_(token_array.astype(bool), token_array.astype(bool))
    ]

    # Get iptm matrix from AF3 summary_confidences file
    iptm_af3 = {
        chain1: {chain2: 0 for chain2 in unique_chains if chain1 != chain2}
        for chain1 in unique_chains
    }

    summary_file_path = None
    if "confidences" in pae_file_path:
        summary_file_path = pae_file_path.replace("confidences", "summary_confidences")
    elif "full_data" in pae_file_path:
        summary_file_path = pae_file_path.replace("full_data", "summary_confidences")

    if summary_file_path is not None and os.path.exists(summary_file_path):
        with open(summary_file_path, "r") as file:
            data_summary = json.load(file)
        af3_chain_pair_iptm_data = data_summary["chain_pair_iptm"]
        for chain1 in unique_chains:
            nchain1 = ord(chain1) - ord("A")  # map A,B,C... to 0,1,2...
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue
                nchain2 = ord(chain2) - ord("A")
                iptm_af3[chain1][chain2] = af3_chain_pair_iptm_data[nchain1][nchain2]
    else:
        print("AF3 summary file does not exist: ", summary_file_path)


def init_chainpairdict_zeros(chainlist: list[str]) -> dict:
    """
    Initializes a nested dictionary with all values set to 0
    Args:
        chainlist: A list of chain IDs
    Returns:
        A nested dictionary with all values set to 0
    Example:
    chainlist = ["A", "B", "C"]
    chainpairdict = init_chainpairdict_zeros(chainlist)
    {
        "A": {"B": 0, "C": 0},
        "B": {"A": 0, "C": 0},
        "C": {"A": 0, "B": 0},
    }
    """
    return {
        chain1: {chain2: 0 for chain2 in chainlist if chain1 != chain2}
        for chain1 in chainlist
    }


def init_chainpairdict_npzeros(chainlist: list, arraysize: tuple | int) -> dict:
    """
    Initializes a nested dictionary with NumPy arrays of zeros for each chain pair.

    Each key in the returned dictionary is a chain identifier from the chainlist.
    For each key (chain1), a sub-dictionary is created with keys for every other chain
    (chain2, where chain1 != chain2), each containing a NumPy array of zeros with
    the given shape.

    Args:
        chainlist (list): A list of chain identifiers.
        arraysize (tuple or int): The shape or size to be passed to np.zeros.

    Returns:
        dict: A nested dictionary where for each chain1, each chain2 (chain1 != chain2)
              has an associated NumPy array of zeros.
    """
    return {
        chain1: {
            chain2: np.zeros(arraysize) for chain2 in chainlist if chain1 != chain2
        }
        for chain1 in chainlist
    }


def init_chainpairdict_set(chainlist: list) -> dict:
    """
    Initializes a nested dictionary with empty sets for each chain pair.

    Each key in the returned dictionary is a chain identifier from the chainlist.
    For each key (chain1), a sub-dictionary is created with keys for
    every other chain (chain2, where chain1 != chain2), each containing an empty set.

    Args:
        chainlist (list): A list of chain identifiers.

    Returns:
        dict: A nested dictionary where for each chain1, each chain2 (chain1 != chain2)
              has an associated empty set.
    Examples:
    >>> chainlist = ["A", "B", "C"]
    >>> init_chainpairdict_set(chainlist)
    {
        "A": {"B": set(), "C": set()},
        "B": {"A": set(), "C": set()},
        "C": {"A": set(), "B": set()},
    }
    """
    return {
        chain1: {chain2: set() for chain2 in chainlist if chain1 != chain2}
        for chain1 in chainlist
    }


def compute_pdockq_all(
    unique_chains, numres, chains, distances, cb_plddt, pae_matrix, cutoff=8.0
):
    """
    各チェーン対ごとにpDockQおよびpDockQ2の値を計算する
    共通処理としてpdockq_unique_residuesを算出し、
    pdockqはcb_plddtとnpairsの情報を用い、
    pdockq2はcb_plddtとPAEに基づくptm値を用いる。
    Returns:
        pdockq: dict
        pdockq2: dict
    """
    # 共通で使用するpdockq_unique_residuesの算出
    pdockq_unique_residues = {
        c1: {c2: set() for c2 in unique_chains if c1 != c2} for c1 in unique_chains
    }
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue
            for i in range(numres):
                if chains[i] != chain1:
                    continue
                valid_pairs = (chains == chain2) & (distances[i] <= cutoff)
                if valid_pairs.any():
                    pdockq_unique_residues[chain1][chain2].add(i)
                    for residue in np.where(valid_pairs)[0]:
                        pdockq_unique_residues[chain1][chain2].add(residue)

    # pDockQ, pDockQ2の計算
    pdockq = {c1: {c2: 0.0 for c2 in unique_chains if c1 != c2} for c1 in unique_chains}
    pdockq2 = {
        c1: {c2: 0.0 for c2 in unique_chains if c1 != c2} for c1 in unique_chains
    }
    lis = {c1: {c2: 0.0 for c2 in unique_chains if c1 != c2} for c1 in unique_chains}
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue
            npairs_pdockq = 0
            npairs_pdockq2 = 0
            total_ptm = 0.0
            for i in range(numres):
                if chains[i] != chain1:
                    continue
                valid_pairs = (chains == chain2) & (distances[i] <= cutoff)
                npairs_pdockq += np.sum(valid_pairs)
                if valid_pairs.any():
                    npairs_pdockq2 += np.sum(valid_pairs)
                    pae_list = pae_matrix[i][valid_pairs]
                    # ptm_func_vecはグローバルで定義済みのベクトル版関数
                    pae_list_ptm = ptm_func_vec(pae_list, 10.0)
                    total_ptm += pae_list_ptm.sum()
            if npairs_pdockq > 0:
                mean_plddt = cb_plddt[
                    list(pdockq_unique_residues[chain1][chain2])
                ].mean()
                x = mean_plddt * math.log10(npairs_pdockq)
                pdockq[chain1][chain2] = (
                    0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
                )
            else:
                pdockq[chain1][chain2] = 0.0

            if npairs_pdockq2 > 0:
                mean_plddt = cb_plddt[
                    list(pdockq_unique_residues[chain1][chain2])
                ].mean()
                mean_ptm = total_ptm / npairs_pdockq2
                x = mean_plddt * mean_ptm
                pdockq2[chain1][chain2] = (
                    1.31 / (1 + math.exp(-0.075 * (x - 84.733))) + 0.005
                )
            else:
                pdockq2[chain1][chain2] = 0.0
            # LIS
            mask = (chains[:, None] == chain1) & (chains[None, :] == chain2)
            selected_pae = pae_matrix[mask]
            if selected_pae.size > 0:
                valid_pae = selected_pae[selected_pae <= 12]
                if valid_pae.size > 0:
                    scores = (12 - valid_pae) / 12
                    lis[chain1][chain2] = np.mean(scores)
                else:
                    lis[chain1][chain2] = 0.0
            else:
                lis[chain1][chain2] = 0.0

    return pdockq, pdockq2, lis


def compute_byres_iptm_ipsae(
    unique_chains, numres, chains, pae_matrix, pae_cutoff, dist_cutoff, residues
):
    """
    各チェーン対ごとに by-residue の ipTM/ipSAE を計算し、
    同時に、ユニークな残基や valid pair 数も記録する。
    """
    # 辞書の初期化（既存の初期化関数を使用）
    iptm_d0chn_byres = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0chn_byres = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0dom_byres = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0res_byres = init_chainpairdict_npzeros(unique_chains, numres)
    n0chn = init_chainpairdict_zeros(unique_chains)
    n0res_byres = init_chainpairdict_npzeros(unique_chains, numres)
    d0res_byres = init_chainpairdict_npzeros(unique_chains, numres)
    valid_pair_counts = init_chainpairdict_zeros(unique_chains)
    unique_residues_chain1 = init_chainpairdict_set(unique_chains)
    unique_residues_chain2 = init_chainpairdict_set(unique_chains)
    dist_valid_pair_counts = init_chainpairdict_zeros(unique_chains)
    dist_unique_residues_chain1 = init_chainpairdict_set(unique_chains)
    dist_unique_residues_chain2 = init_chainpairdict_set(unique_chains)

    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue
            # 全体残基数 n0chn, d0chn の算出
            n0chn[chain1][chain2] = np.sum(chains == chain1) + np.sum(chains == chain2)
            d0 = calc_d0(n0chn[chain1][chain2])
            ptm_matrix_d0chn = ptm_func_vec(pae_matrix, d0)

            valid_pairs_iptm = chains == chain2
            valid_pairs_matrix = (chains == chain2) & (pae_matrix < pae_cutoff)
            # Assuming valid_pairs_matrix is already defined
            n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
            d0res_byres_all = calc_d0_array(n0res_byres_all)

            n0res_byres[chain1][chain2] = n0res_byres_all
            d0res_byres[chain1][chain2] = d0res_byres_all
            for i in range(numres):
                if chains[i] != chain1:
                    continue
                valid_pairs_ipsae = valid_pairs_matrix[i]
                iptm_d0chn_byres[chain1][chain2][i] = (
                    ptm_matrix_d0chn[i, valid_pairs_iptm].mean()
                    if valid_pairs_iptm.any()
                    else 0.0
                )
                ipsae_d0chn_byres[chain1][chain2][i] = (
                    ptm_matrix_d0chn[i, valid_pairs_ipsae].mean()
                    if valid_pairs_ipsae.any()
                    else 0.0
                )
                # Track unique residues contributing to the IPSAE for chain1,chain2
                valid_pair_counts[chain1][chain2] += np.sum(valid_pairs_ipsae)
                if valid_pairs_ipsae.any():
                    iresnum = residues[i]["resnum"]
                    unique_residues_chain1[chain1][chain2].add(iresnum)
                    for j in np.where(valid_pairs_ipsae)[0]:
                        unique_residues_chain2[chain1][chain2].add(
                            residues[j]["resnum"]
                        )

                # Track unique residues contributing to iptm in interface
                valid_pairs = (
                    (chains == chain2)
                    & (pae_matrix[i] < pae_cutoff)
                    & (distances[i] < dist_cutoff)
                )
                dist_valid_pair_counts[chain1][chain2] += np.sum(valid_pairs)

                # Track unique residues contributing to the IPTM
                if valid_pairs.any():
                    iresnum = residues[i]["resnum"]
                    dist_unique_residues_chain1[chain1][chain2].add(iresnum)
                    for j in np.where(valid_pairs)[0]:
                        dist_unique_residues_chain2[chain1][chain2].add(
                            residues[j]["resnum"]
                        )

    results = {
        "iptm_d0chn_byres": iptm_d0chn_byres,
        "ipsae_d0chn_byres": ipsae_d0chn_byres,
        "ipsae_d0dom_byres": ipsae_d0dom_byres,
        "ipsae_d0res_byres": ipsae_d0res_byres,
        "n0chn": n0chn,
        "n0res_byres": n0res_byres,
        "d0res_byres": d0res_byres,
        "valid_pair_counts": valid_pair_counts,
        "unique_residues_chain1": unique_residues_chain1,
        "unique_residues_chain2": unique_residues_chain2,
        "dist_valid_pair_counts": dist_valid_pair_counts,
        "dist_unique_residues_chain1": dist_unique_residues_chain1,
        "dist_unique_residues_chain2": dist_unique_residues_chain2,
    }
    return results


def compute_interchain_max(
    unique_chains,
    iptm_d0chn_byres,
    ipsae_d0chn_byres,
    ipsae_d0dom_byres,
    ipsae_d0res_byres,
    n0res_byres,
    d0res_byres,
):
    """
    各チェーン対について、by-residue の ipTM/ipSAE 値の中で
    最大値（asymmetric, max）を求める。
    """
    # 各値の asymmetric および最大値用辞書の初期化（init_chainpairdict_zeros 使用）
    iptm_d0chn_asym = init_chainpairdict_zeros(unique_chains)
    ipsae_d0chn_asym = init_chainpairdict_zeros(unique_chains)
    ipsae_d0dom_asym = init_chainpairdict_zeros(unique_chains)
    ipsae_d0res_asym = init_chainpairdict_zeros(unique_chains)
    n0res = init_chainpairdict_zeros(unique_chains)
    d0res = init_chainpairdict_zeros(unique_chains)

    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue
            # by-residue毎の値から最大値とその residue を算出
            interchain = iptm_d0chn_byres[chain1][chain2]
            max_index = np.argmax(interchain)
            iptm_d0chn_asym[chain1][chain2] = interchain[max_index]

            interchain = ipsae_d0chn_byres[chain1][chain2]
            max_index = np.argmax(interchain)
            ipsae_d0chn_asym[chain1][chain2] = interchain[max_index]

            interchain = ipsae_d0dom_byres[chain1][chain2]
            max_index = np.argmax(interchain)
            ipsae_d0dom_asym[chain1][chain2] = interchain[max_index]

            interchain = ipsae_d0res_byres[chain1][chain2]
            max_index = np.argmax(interchain)
            ipsae_d0res_asym[chain1][chain2] = interchain[max_index]

            n0res[chain1][chain2] = n0res_byres[chain1][chain2][max_index]
            d0res[chain1][chain2] = d0res_byres[chain1][chain2][max_index]

    return {
        "iptm_d0chn_asym": iptm_d0chn_asym,
        "ipsae_d0chn_asym": ipsae_d0chn_asym,
        "ipsae_d0dom_asym": ipsae_d0dom_asym,
        "ipsae_d0res_asym": ipsae_d0res_asym,
        "n0res": n0res,
        "d0res": d0res,
    }


pdockq, pdockq2, lis = compute_pdockq_all(
    unique_chains, numres, chains, distances, cb_plddt, pae_matrix, cutoff=8.0
)

results_byres = compute_byres_iptm_ipsae(
    unique_chains, numres, chains, pae_matrix, pae_cutoff, dist_cutoff, residues
)

interchain_results = compute_interchain_max(
    unique_chains,
    results_byres["iptm_d0chn_byres"],
    results_byres["ipsae_d0chn_byres"],
    results_byres["ipsae_d0dom_byres"],
    results_byres["ipsae_d0res_byres"],
    results_byres["n0res_byres"],
    results_byres["d0res_byres"],
)


def write_byres(
    byres_filepath,
    unique_chains,
    numres,
    chains,
    pae_matrix,
    pae_cutoff,
    plddt,
    residues,
    results_byres,
):
    """
    compute_byres_iptm_ipsae などで計算済みの辞書（results_byres）と、
    グローバル変数（pae_matrix、pae_cutoff、plddt、chains、residuesなど）を用いて、
    各残基ごとの ipTM/ipSAE 出力結果をファイルへ書き出す。

    results_byres には以下のキーが含まれることを前提とする:
      - "iptm_d0chn_byres"
      - "ipsae_d0chn_byres"
      - "ipsae_d0dom_byres"
      - "ipsae_d0res_byres"
      - "n0chn"
      - "n0res_byres"
      - "unique_residues_chain1"
      - "unique_residues_chain2"
    """
    # obtain variables from results_byres
    iptm_d0chn_byres = results_byres["iptm_d0chn_byres"]
    ipsae_d0chn_byres = results_byres["ipsae_d0chn_byres"]
    ipsae_d0dom_byres = results_byres["ipsae_d0dom_byres"]
    ipsae_d0res_byres = results_byres["ipsae_d0res_byres"]
    d0res_byres = results_byres["d0res_byres"]
    n0chn = results_byres["n0chn"]
    n0res_byres = results_byres["n0res_byres"]
    # unique_residues is a nested dict of sets
    unique_residues_chain1 = results_byres["unique_residues_chain1"]
    unique_residues_chain2 = results_byres["unique_residues_chain2"]

    # output header
    with open(byres_filepath, "w") as OUT2:
        OUT2.write(
            "i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      "
            "n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  ipSAE_d0chn "
            "ipSAE_d0dom    ipSAE \n"
        )

        # d0chn, n0dom, d0dom の計算
        d0chn = {c1: {} for c1 in unique_chains}
        n0dom = {c1: {} for c1 in unique_chains}
        d0dom = {c1: {} for c1 in unique_chains}
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue
                # n0chn: sum of the residues belonging to chain1 and chain2
                n0chn_val = np.sum(chains == chain1) + np.sum(chains == chain2)
                # d0chn: calc_d0 applied to n0chn
                d0chn[chain1][chain2] = calc_d0(n0chn_val)
                # n0dom: number of unique residues in chain1 and chain2
                residues_1 = len(unique_residues_chain1[chain1][chain2])
                residues_2 = len(unique_residues_chain2[chain1][chain2])
                n0dom[chain1][chain2] = residues_1 + residues_2
                d0dom[chain1][chain2] = calc_d0(n0dom[chain1][chain2])

        # calculate each chain pair's by-residue output
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue

                ptm_matrix_d0dom = ptm_func_vec(pae_matrix, d0dom[chain1][chain2])
                # valid_pairs_matrix: mask where residues belong to chain2 and pae_matrix < pae_cutoff
                valid_pairs_matrix = (chains == chain2) & (pae_matrix < pae_cutoff)
                # n0res_byres_all: sum of valid pairs for each residue
                n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
                d0res_byres_all = calc_d0_array(n0res_byres_all)

                n0res_byres[chain1][chain2] = n0res_byres_all
                d0res_byres[chain1][chain2] = d0res_byres_all

                # calculate outputs for each residue
                for i in range(numres):
                    if chains[i] != chain1:
                        continue

                    valid_pairs = valid_pairs_matrix[i]
                    ipsae_d0dom_byres[chain1][chain2][i] = (
                        ptm_matrix_d0dom[i, valid_pairs].mean()
                        if valid_pairs.any()
                        else 0.0
                    )

                    ptm_row_d0res = ptm_func_vec(
                        pae_matrix[i], d0res_byres[chain1][chain2][i]
                    )
                    ipsae_d0res_byres[chain1][chain2][i] = (
                        ptm_row_d0res[valid_pairs].mean() if valid_pairs.any() else 0.0
                    )

                    outstring = f"{i+1:<4d}    " + (
                        f"{chain1:4}      "
                        f"{chain2:4}      "
                        f'{residues[i]["resnum"]:4d}           '
                        f'{residues[i]["res"]:3}        '
                        f"{plddt[i]:8.2f}         "
                        f"{int(n0chn[chain1][chain2]):5d}  "
                        f"{int(n0dom[chain1][chain2]):5d}  "
                        f"{int(n0res_byres[chain1][chain2][i]):5d}  "
                        f"{d0chn[chain1][chain2]:8.3f}  "
                        f"{d0dom[chain1][chain2]:8.3f}  "
                        f"{d0res_byres[chain1][chain2][i]:8.3f}   "
                        f"{iptm_d0chn_byres[chain1][chain2][i]:8.4f}    "
                        f"{ipsae_d0chn_byres[chain1][chain2][i]:8.4f}    "
                        f"{ipsae_d0dom_byres[chain1][chain2][i]:8.4f}    "
                        f"{ipsae_d0res_byres[chain1][chain2][i]:8.4f}\n"
                    )
                    OUT2.write(outstring)


write_byres(
    "/Users/YoshitakaM/Desktop/byres.txt",
    unique_chains,
    numres,
    chains,
    pae_matrix,
    pae_cutoff,
    plddt,
    residues,
    results_byres,
)
