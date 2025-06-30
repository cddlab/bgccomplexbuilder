# %%
import json
import os
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger

from complexbuilder.common.log import log_setup

log_setup(level="DEBUG")


def classify_proteins(
    genbank_file: str | Path, max_length: int = 1500, decompose_nrpspks: bool = False
) -> tuple[list[SeqRecord], list[SeqRecord]]:
    """Collect protein sequences from GenBank file and return a list of
    SeqRecord objects.
    The structural domains of the protein belonging to NRPS_PKS are
    returned together as a separate list.

    Args:
        - genbank_file (str | Path): Path to the GenBank file.
        - max_length (int): Maximum length of protein sequences
    Returns:
        - nonproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that do not belong to NRPS_PKS.
        - proteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that belong to NRPS_PKS.
        - decompose_nrpspks (bool): If True, decompose the NRPS_PKS proteins
            into individual domains.
    """

    nonproteins: list[SeqRecord] = []
    proteins: list[SeqRecord] = []
    genbank_file = Path(genbank_file)
    for record in SeqIO.parse(genbank_file, "genbank"):
        logger.info(f"Record ID: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                # protein_id and aa_sequence are required fields
                aa_sequence = feature.qualifiers.get("translation")[0]
                if "protein_id" in feature.qualifiers:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                elif "locus_tag" in feature.qualifiers:
                    # fallback to locus_tag if protein_id is not available
                    protein_id = feature.qualifiers.get("locus_tag")[0]
                else:
                    logger.warning(
                        f"No protein ID and locus tag found in {genbank_file} . "
                        "Use gene ID."
                    )
                    protein_id = feature.qualifiers.get("gene")[0]
                if "product" in feature.qualifiers:
                    product = feature.qualifiers.get("product")[0]
                else:
                    product = ""
                # truncate the protein sequence if it exceeds <clip_length> residues
                if "NRPS_PKS" in feature.qualifiers:
                    pattern = r"Domain: \S+ \((\d+)-(\d+)\).*nrpspksdomains_(\S+)"
                    if len(aa_sequence) < max_length:
                        proteins.append(
                            SeqRecord(
                                Seq(aa_sequence),
                                id=protein_id,
                                description=product,
                            )
                        )
                    elif decompose_nrpspks:
                        for idx in range(len(feature.qualifiers["NRPS_PKS"])):
                            match = re.search(
                                pattern, feature.qualifiers["NRPS_PKS"][idx]
                            )
                            if match:
                                start_res = int(match.group(1))
                                end_res = int(match.group(2))
                                domainname = match.group(3)
                                proteins.append(
                                    SeqRecord(
                                        Seq(aa_sequence[start_res:end_res]),
                                        id=f"{protein_id}@{idx}",
                                        description=domainname,
                                    )
                                )
                            else:
                                raise ValueError(
                                    "Domain residue regions not found for Record ID: "
                                    f"{record.id}"
                                )
                else:
                    nonproteins.append(
                        SeqRecord(Seq(aa_sequence), id=protein_id, description=product)
                    )
    if nonproteins == []:
        logger.warning(
            "No non-NRPS/PKS protein sequences found in "
            f"the GenBank file {genbank_file} ."
        )
    if proteins == []:
        logger.warning(
            f"No NRPS/PKS domains found in the GenBank file {genbank_file} ."
        )
    logger.info(f"Number of non-NRPS/PKS proteins: {len(nonproteins)}")
    logger.info(f"Number of NRPS/PKS proteins: {len(proteins)}")
    return nonproteins, proteins


def classify_proteins2(
    genbank_file: str | Path, max_length: int = 1950, decompose_nrpspks: bool = False
) -> list[SeqRecord]:
    """Collect protein sequences from GenBank file and return a list of
    SeqRecord objects.
    The structural domains of the protein belonging to NRPS_PKS are
    returned together as a separate list.

    Args:
        - genbank_file (str | Path): Path to the GenBank file.
        - max_length (int): Maximum length of protein sequences
    Returns:
        - nonproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that do not belong to NRPS_PKS.
        - proteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that belong to NRPS_PKS.
        - decompose_nrpspks (bool): If True, decompose the NRPS_PKS proteins
            into individual domains.
    """
    proteins: list[SeqRecord] = []
    genbank_file = Path(genbank_file)
    for record in SeqIO.parse(genbank_file, "genbank"):
        logger.info(f"Record ID: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                # protein_id and aa_sequence are required fields
                aa_sequence = feature.qualifiers.get("translation")[0]
                if "protein_id" in feature.qualifiers:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                elif "locus_tag" in feature.qualifiers:
                    # fallback to locus_tag if protein_id is not available
                    protein_id = feature.qualifiers.get("locus_tag")[0]
                else:
                    logger.warning(
                        f"No protein ID and locus tag found in {genbank_file} . "
                        "Use gene ID."
                    )
                    protein_id = feature.qualifiers.get("gene")[0]
                if "product" in feature.qualifiers:
                    product = feature.qualifiers.get("product")[0]
                else:
                    product = ""
                if len(aa_sequence) < max_length:
                    proteins.append(
                        SeqRecord(
                            Seq(aa_sequence),
                            id=protein_id,
                            description=product,
                        )
                    )
                else:
                    logger.warning(
                        f"The protein sequence {protein_id} is too long. "
                        f"Length: {len(aa_sequence)} > {max_length}"
                    )
    logger.info(f"Number of proteins: {len(proteins)}")
    return proteins


def split_proteinids(text: str) -> tuple[str, str]:
    """
    Split protein complexes by an underscore that is followed by a lowercase letter.

    The function looks for underscores ("_") in the input string where the character
    immediately following the underscore is a letter in [a-z]. It then splits the string
    at that underscore into exactly two parts. If no such underscore is found, or if more
    than one valid splitting underscore is found (thus producing more than two parts), an
    error is raised.
    Examples:
        split_proteinids("adakdk.aa_adakdk.aa")
            -> ("adakdk.aa", "adakdk.aa")
        split_proteinids("alsk.01_alsk.01")
            -> ("alsk.01", "alsk.01")
        split_proteinids("trx17522.1_trx20192.1")
            -> ("trx17522.1", "trx20192.1")
        split_proteinids("trx17522.1_fnf07_04290")
            -> ("trx17522.1", "fnf07_04290")
        split_proteinids("ctg1_orf2_ctg1_orf2")
            -> ("ctg1_orf2", "ctg1_orf2")

    Args:
        text (str): The protein complex string to split.

    Returns:
        tuple[str, str]: A tuple with the two parts resulting from the split.

    Raises:
        ValueError: If no valid splitting underscore is found or if an even number of them is detected.
    """
    # Find indices of underscores immediately followed by a lowercase letter.
    # But not followed by "rs" (e.g., "ssgg_rs34700").
    valid_indices = [m.start() for m in re.finditer(r"_(?!rs)(?=[a-z])", text)]
    if not valid_indices:
        raise ValueError(f"No valid splitting underscore found in: {text}")

    if len(valid_indices) == 1:
        split_index = valid_indices[0]
    else:
        # If more than one valid underscore is found, require an odd number.
        if len(valid_indices) % 2 == 0:
            raise ValueError(
                f"Even number of valid splitting underscores found in: {text}"
            )
        # Pick the middle valid underscore.
        mid = len(valid_indices) // 2
        split_index = valid_indices[mid]

    left = text[:split_index]
    right = text[split_index + 1 :]
    return (left, right)


def _check_homo_hetero(dirname: str) -> str:
    left, right = split_proteinids(dirname)
    if left == right:
        return "homo"
    else:
        return "hetero"


def make_hitcomplexlist(
    dir: str | Path,
    ipsae_threshold: float,
    iptm_threshold: float,
    metrics_json: str = "complexmetrics.json",
):
    """
    Create a dictionary of hit complexes based on ipSAE and ipTM thresholds.

    This function scans the specified directory for subdirectories whose names
    start with "BGC000". For each such subdirectory, it looks for a JSON file
    (by default "complexmetrics.json"). If the JSON file exists, the function
    loads the metrics from the "ipSAE" and "ipTM" sections. Each section is expected
    to be a list of single-key dictionaries, mapping a complex identifier to a numeric score.

    A complex is considered a "hit" if its ipSAE and ipTM scores meet or exceed the
    provided thresholds, and it appears in both sections. Additionally, the assembly
    type is determined by the _check_homo_hetero function and is stored as
    "complex_homo_hetero" with values "homo", "hetero", or "otherwise".

    Args:
        dir (str | Path): The path to the directory containing BGC subdirectories.
        ipsae_threshold (float): The minimum ipSAE score required for a complex to be considered a hit.
        iptm_threshold (float): The minimum ipTM score required for a complex to be considered a hit.
        metrics_json (str, optional): The filename of the metrics JSON file in each subdirectory.
            Defaults to "complexmetrics.json".

    Returns:
        dict: A dictionary where each key is a BGC directory name (e.g., "BGC0000001") and each value
              is another dictionary mapping complex identifiers to a dictionary with the following keys:
                  - "ipSAE": The ipSAE score.
                  - "ipTM": The ipTM score.
                  - "complex_homo_hetero": A string indicating the assembly type ("homo", "hetero",
                    or "otherwise").

    Notes:
        If a subdirectory does not contain the specified metrics JSON file, it is skipped with a warning.
    """
    dirname = Path(dir)
    bgcdirs = [
        p.name for p in dirname.iterdir() if p.is_dir() and p.name.startswith("BGC000")
    ]
    # sort by the number in the directory name
    bgcdirs.sort(key=lambda x: int(x.split("BGC")[1]))
    logger.debug(f"Subdirectories in {dirname}: {bgcdirs}")

    # sort by the number in the directory name
    bgcdirs.sort(key=lambda x: int(x.split("BGC")[1]))

    results = {}
    for bgcdir in bgcdirs:
        bgcpath = dirname / bgcdir
        logger.debug(f"Processing {bgcpath}")
        if not os.path.exists(f"{bgcpath}/{metrics_json}"):
            logger.warning(
                f"complexmetrics.json not found in {bgcpath}/{metrics_json} . Skipping."
            )
            continue
        with open(f"{bgcpath}/{metrics_json}") as f:
            complexmetrics = json.load(f)

        ipSAE = {
            list(item.keys())[0]: list(item.values())[0]
            for item in complexmetrics.get("ipSAE", [])
        }
        ipTM = {
            list(item.keys())[0]: list(item.values())[0]
            for item in complexmetrics.get("ipTM", [])
        }

        hitcomplexes = {}
        for key in ipSAE:
            if (
                key in ipTM
                and ipSAE[key] >= ipsae_threshold
                and ipTM[key] >= iptm_threshold
            ):
                hitcomplexes[key] = {
                    "ipSAE": ipSAE[key],
                    "ipTM": ipTM[key],
                    "complex_homo_hetero": _check_homo_hetero(key),
                }
        results[bgcdir] = hitcomplexes
    return results


# %%

dirname = Path("/Users/YoshitakaM/Desktop/BGCcomplex")
ipsae_threshold = 0.6
iptm_threshold = 0.8
print(make_hitcomplexlist(dirname, ipsae_threshold, iptm_threshold))

# %%
