# %%
import argparse
import json
import os
import re
from pathlib import Path

from loguru import logger

from complexbuilder.common.log import log_setup

log_setup(level="SUCCESS")


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
    display_homo_hetero: bool = False,
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
        display_homo_hetero (bool, optional): Whether to include the assembly type in the output.
            Defaults to False.

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
                if display_homo_hetero:
                    hitcomplexes[key] = {
                        "ipSAE": ipSAE[key],
                        "ipTM": ipTM[key],
                        "complex_homo_hetero": _check_homo_hetero(key),
                    }
                else:
                    hitcomplexes[key] = {
                        "ipSAE": ipSAE[key],
                        "ipTM": ipTM[key],
                    }
        results[bgcdir] = hitcomplexes
    return results


def main():
    parser = argparse.ArgumentParser(description="Extract hit complexes from BGCs.")
    parser.add_argument(
        "--bgc_directory",
        metavar="BGC directory",
        default="/data2/moriwaki/BGCcomplex/merged",
        type=str,
        help="Path to root direcotory of BGCs.",
    )
    parser.add_argument(
        "--outputjson",
        metavar="Output JSON file",
        default="hitcomplexes.json",
        type=str,
        help="The filename of the output JSON file.",
    )
    parser.add_argument(
        "--iptm_threshold",
        metavar="ipTM threshold",
        default=0.6,
        type=float,
        help="Minimum ipTM score required for a complex to be considered a hit.",
    )
    parser.add_argument(
        "--ipsae_threshold",
        metavar="ipSAE threshold",
        default=0.0,
        type=float,
        help="Minimum ipSAE score required for a complex to be considered a hit.",
    )
    parser.add_argument(
        "--metrics_json",
        metavar="Metrics JSON file",
        default="complexmetrics.json",
        type=str,
        help="The filename of the metrics JSON file in each subdirectory.",
    )
    parser.add_argument(
        "--display_homo_hetero",
        action="store_true",
        help="Whether to include the assembly type in the output.",
    )
    # args = parser.parse_args()
    args = parser.parse_args(
        [
            "--bgc_directory",
            "/data2/moriwaki/BGCcomplex/homomers_additional",
            "--outputjson",
            "hitcomplexes.json",
            "--iptm_threshold",
            "0.0",
            "--ipsae_threshold",
            "0.0",
            "--metrics_json",
            "complexmetrics.json",
        ]
    )
    make_hitcomplexlist(
        dir=args.bgc_directory,
        ipsae_threshold=args.ipsae_threshold,
        iptm_threshold=args.iptm_threshold,
        metrics_json=args.metrics_json,
        display_homo_hetero=args.display_homo_hetero,
    )


if __name__ == "__main__":
    main()

# %%
