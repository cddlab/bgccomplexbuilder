import argparse
import json
from pathlib import Path

import gemmi
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger
from matplotlib import rcParams

from complexbuilder.analysis.extract_hitcomplexes import split_proteinids
from complexbuilder.common.log import log_setup

log_setup(level="INFO")


def calculate_rmsd_between_two_chains(
    ciffile: str | Path, chain1_id: str = "A", chain2_id: str = "B"
) -> float | None:
    """
    Compare two chains in a CIF file and return the RMSD value.

    Parameters:
        ciffile (str | Path): Path to the CIF file.
        chain1_id (str): ID of the first chain.
        chain2_id (str): ID of the second chain.

    Returns:
        float: RMSD value between the two chains.

    Note:
        If the chain length of either chain is less than 20, return "None"
        because the RMSD calculation is not reliable for short chains.
    """
    ciffile_str = str(ciffile) if isinstance(ciffile, Path) else ciffile
    if not Path(ciffile_str).exists():
        raise FileNotFoundError(f"CIF file not found: {ciffile_str}")
    model = gemmi.read_structure(ciffile_str)[0]
    chain1 = model[chain1_id].get_polymer()
    chain2 = model[chain2_id].get_polymer()

    # Ensure both chains are not empty
    assert chain1 is not None and chain2 is not None, "One of the chains is empty."
    # if the length of the chains is less than 20, return None
    if len(chain1) < 20 or len(chain2) < 20:
        return None

    ptype = chain1.check_polymer_type()
    sr = gemmi.calculate_superposition(
        chain1, chain2, ptype, gemmi.SupSelect.MainChain, trim_cycles=5
    )
    return round(sr.rmsd, 3)


def extract_hetero_complexes(json_data: dict) -> list[tuple[str, str, float, float]]:
    return [
        (bgc_id, complex_key, complex_data.get("ipSAE"), complex_data.get("ipTM"))
        for bgc_id, complexes in json_data.items()
        for complex_key, complex_data in complexes.items()
        if complex_data.get("complex_homo_hetero") == "hetero"
    ]


def draw_histgram(rmsdfile: str | Path, outfigpath: str | Path) -> None:
    """
    Draw a histogram of RMSD values from a JSON file.
    Args:
        rmsdfile (str | Path): Path to the RMSD JSON file.
        outfigpath (str | Path): Path to the output figure file.
    """

    with open(rmsdfile, "r") as f:
        rmsd_data = json.load(f)

    rmsd_values = [
        complex_data["RMSD"]
        for bgc_id, complexes in rmsd_data.items()
        for complex_key, complex_data in complexes.items()
        if complex_data.get("RMSD") is not None
    ]
    df = pd.DataFrame(rmsd_values, columns=["RMSD"])
    # df["RMSD"] <= 2.0の割合を計算
    total = df["RMSD"].count()
    num_below_1 = (df["RMSD"] <= 1.0).sum()
    fraction_below_1 = (df["RMSD"] <= 1.0).mean() * 100
    num_below_2 = (df["RMSD"] <= 2.0).sum()
    fraction_below_2 = (df["RMSD"] <= 2.0).mean() * 100
    logger.info(f"Number of complexes with RMSD <= 1.0: {num_below_1} / {total}")
    logger.info(f"Fraction of complexes with RMSD <= 1.0: {fraction_below_1:.3f} %")
    logger.info(f"Number of complexes with RMSD <= 2.0: {num_below_2} / {total}")
    logger.info(f"Fraction of complexes with RMSD <= 2.0: {fraction_below_2:.3f} %")
    rcParams["font.family"] = "Arial"
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"] = 9
    rcParams["axes.titlesize"] = 9
    rcParams["xtick.labelsize"] = 8
    rcParams["ytick.labelsize"] = 8
    rcParams["mathtext.fontset"] = "cm"

    fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4), dpi=300)
    bin_list = [0 + i * 1 for i in range(0, 101)]
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.hist(df["RMSD"], bins=bin_list, color="#6CD8FD")
    ax.set_xlim(0, 100)
    ax.set_xticks(range(0, 101, 10))
    ax.set_xlabel(r"rmsd [Å]")
    ax.set_ylabel("count")
    ax.set_yticks(range(0, 2000, 500))
    ax.set_ylim(0, 1500)

    axins = ax.inset_axes(
        (0.7, 0.7, 0.27, 0.27),
    )
    bin_list2 = [0 + i * 1 for i in range(0, 11)]
    axins.spines["right"].set_visible(False)
    axins.spines["top"].set_visible(False)
    axins.hist(df["RMSD"], bins=bin_list2, color="#6CD8FD")
    axins.tick_params(labelsize=6)
    axins.set_xlim(0, 10)
    axins.set_xticks(range(0, 11, 2))
    axins.set_xlabel(r"rmsd [Å]", fontsize=6)
    axins.set_ylim(0, 1000)
    axins.set_yticks(range(0, 1200, 200))
    axins.set_ylabel("count", fontsize=6)
    fig.savefig(
        outfigpath,
        bbox_inches="tight",
        pad_inches=0.05,
    )


def write_rmsd_results_to_json(
    bgcdirectory_root: str | Path,
    hitjsonfile: str | Path,
    outputjson: str | Path,
    maxcount: int = 9999,
) -> None:
    """
    Write RMSD results to a JSON file.
    Args:
        bgcdirectory_root (str | Path): Root directory containing BGC data.
        hitjsonfile (str | Path): Path to the input JSON file containing hit complexes.
        outputjson (str | Path): Path to the output JSON file.
        maxcount (int): Maximum number of complexes to process.
    """
    with open(hitjsonfile, "r") as f:
        hit_data = json.load(f)

    count = 0
    results = {}
    with open(outputjson, "w") as f:
        for bgc_id, complex_key, ipSAE, ipTM in extract_hetero_complexes(hit_data):
            logger.info(f"Processing: {bgc_id}, {complex_key}")
            results[bgc_id] = results[bgc_id] if bgc_id in results else {}
            if Path(
                f"{bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
            ).exists():
                rmsd = calculate_rmsd_between_two_chains(
                    f"{bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
                )
                if rmsd is not None:
                    results[bgc_id][complex_key] = {"RMSD": rmsd}
                    results[bgc_id][complex_key]["ipSAE"] = ipSAE
                    results[bgc_id][complex_key]["ipTM"] = ipTM
            else:
                logger.warning(
                    f"File not found: {bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
                )
                results[bgc_id][complex_key] = {"RMSD": None}
            count += 1
            if count >= maxcount:
                logger.info(f"Reached max count of {maxcount}. Stopping.")
                break

    with open(outputjson, "w") as f:
        json.dump(results, f, indent=2)
        logger.info(f"RMSD results saved to {outputjson}")


def make_qualified_ipsae_rmsd_dataframe(
    rmsdfile: str | Path,
    hitjsonfile: str | Path,
    output_suffix: str | Path,
    homoipsae_threshold: float = 0.6,
    heteroipsae_threshold: float = 0.6,
) -> None:
    """
    Create a DataFrame of qualified complexes based on RMSD and ipSAE thresholds
    and save it to a JSON file.
    Args:
        rmsdfile (str | Path): Path to the RMSD JSON file.
        hitjsonfile (str | Path): Path to the hit complexes JSON file.
        output (str | Path): Path to the output JSON file.
        homoipsae_threshold (float): Threshold for homology ipSAE.
        heteroipsae_threshold (float): Threshold for heterology ipSAE.
    """

    with open(rmsdfile, "r") as f:
        data = json.load(f)
    rows = []
    for bgc_id, complexes in data.items():
        for complex_key, metrics in complexes.items():
            prot1, prot2 = split_proteinids(complex_key)
            rows.append(
                {
                    "BGC": bgc_id,
                    "protein1": prot1,
                    "protein2": prot2,
                    "RMSD": metrics.get("RMSD"),
                    "ipSAE": metrics.get("ipSAE"),
                    "ipTM": metrics.get("ipTM"),
                }
            )

    heterodf = pd.DataFrame(
        rows, columns=["BGC", "protein1", "protein2", "RMSD", "ipSAE", "ipTM"]
    )
    heterodf_rmsd = heterodf[heterodf["RMSD"] <= 2.0][heterodf["ipTM"] >= 0.6]
    print(heterodf_rmsd.head())
    with open(hitjsonfile, "r") as f:
        hit_data = json.load(f)
    hitrows = []
    for bgc_id, complexes in hit_data.items():
        for complex_key, metrics in complexes.items():
            hitrows.append(
                {
                    "BGC": bgc_id,
                    "protein1": split_proteinids(complex_key)[0],
                    "protein2": split_proteinids(complex_key)[1],
                    "ipSAE": metrics.get("ipSAE"),
                    "ipTM": metrics.get("ipTM"),
                }
            )
    hitdf = pd.DataFrame(
        hitrows, columns=["BGC", "protein1", "protein2", "ipSAE", "ipTM"]
    )
    homodf = hitdf[hitdf["protein1"] == hitdf["protein2"]].copy()
    homo_ipTM_map = homodf.set_index(["BGC", "protein1"])["ipTM"].to_dict()
    homo_ipSAE_map = homodf.set_index(["BGC", "protein1"])["ipSAE"].to_dict()
    heterodf_rmsd = heterodf_rmsd.copy()
    heterodf_rmsd["prot1_homo_ipTM"] = [
        homo_ipTM_map.get((bgc, p1))
        for bgc, p1 in zip(
            heterodf_rmsd["BGC"], heterodf_rmsd["protein1"], strict=False
        )
    ]
    heterodf_rmsd["prot1_homo_ipSAE"] = [
        homo_ipSAE_map.get((bgc, p1))
        for bgc, p1 in zip(
            heterodf_rmsd["BGC"], heterodf_rmsd["protein1"], strict=False
        )
    ]
    heterodf_rmsd["prot2_homo_ipTM"] = [
        homo_ipTM_map.get((bgc, p2))
        for bgc, p2 in zip(
            heterodf_rmsd["BGC"], heterodf_rmsd["protein2"], strict=False
        )
    ]
    heterodf_rmsd["prot2_homo_ipSAE"] = [
        homo_ipSAE_map.get((bgc, p2))
        for bgc, p2 in zip(
            heterodf_rmsd["BGC"], heterodf_rmsd["protein2"], strict=False
        )
    ]

    qualified = heterodf_rmsd[
        (heterodf_rmsd["prot1_homo_ipSAE"] >= homoipsae_threshold)
        & (heterodf_rmsd["prot2_homo_ipSAE"] >= homoipsae_threshold)
        & (heterodf_rmsd["ipSAE"] > heteroipsae_threshold)
    ]

    outputname = f"{output_suffix}_homoipsae{homoipsae_threshold}_heteroipsae{heteroipsae_threshold}.json"
    with open(outputname, "w") as f:
        json.dump(qualified.to_dict(orient="records"), f, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Process MIBiG GenBank files.")
    parser.add_argument(
        "--bgcdirectoryroot",
        metavar="BGC_DIRECTORY_ROOT",
        type=str,
        help="Directory to predicted BGCs. It should contain subdirectories for each BGC. "
        "e.g. BGC0000001, BGC0000002, ... BGC0002826",
        required=True,
    )
    parser.add_argument(
        "--hitjsonfile",
        metavar="HIT_JSON_FILE",
        type=str,
        help="Path to the hit JSON file.",
        required=True,
    )
    parser.add_argument(
        "--outputrmsdjson",
        metavar="OUTPUT_RMSD_JSON_FILE",
        type=str,
        help="Path to the output RMSD JSON file.",
        default="rmsd.json",
    )
    parser.add_argument(
        "--homoipsae_threshold",
        metavar="HOMO_IPSAE_THRESHOLD",
        type=float,
        default=0.6,
    )
    parser.add_argument(
        "--heteroipsae_threshold",
        metavar="HETERO_IPSAE_THRESHOLD",
        type=float,
        default=0.6,
    )
    # switch whether to compute RMSD or not
    # default is False
    parser.add_argument(
        "--compute_rmsd",
        action="store_true",
        help="If set, compute RMSD and write to output JSON file.",
    )
    args = parser.parse_args()

    if args.compute_rmsd:
        write_rmsd_results_to_json(
            bgcdirectory_root=args.bgcdirectoryroot,
            hitjsonfile=args.hitjsonfile,
            outputjson=args.outputrmsdjson,
        )
    rmsdfile = args.outputrmsdjson
    outfigpath = Path(rmsdfile).stem + "_histogram.svg"
    draw_histgram(rmsdfile, outfigpath)

    make_qualified_ipsae_rmsd_dataframe(
        rmsdfile=args.outputjson,
        hitjsonfile=args.hitjsonfile,
        output_suffix=args.outputjson,
        homoipsae_threshold=args.homoipsae_threshold,
        heteroipsae_threshold=args.heteroipsae_threshold,
    )


if __name__ == "__main__":
    main()
