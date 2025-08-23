#!/usr/bin/env python3
import argparse
import json
import os

from alphafold3tools import paeplot
from alphafold3tools.log import log_setup
from loguru import logger

log_setup(level="INFO")


def run_ipsae_on_af3(
    bgc_directory: str,
    python_binary: str,
    ipsae_py_script: str,
    bgc_number: int,
    pae_cutoff: float = 10,
    dist_cutoff: float = 10,
    skip_paefigs: bool = False,
    overwrite: bool = False,
) -> None:
    target_directory = f"{bgc_directory}/BGC{bgc_number:07d}"
    if not os.path.exists(f"{target_directory}"):
        print(f"af3 directory not found in {target_directory}")
    elif os.path.exists(f"{target_directory}/complexmetrics.json") and not overwrite:
        print(f"'complexmetrics.json' already exists on {target_directory}")
    else:
        ipsaes = []
        ipsae_mins = []
        af3_iptms = []
        af3_subdirs = [
            d
            for d in os.listdir(f"{target_directory}")
            if os.path.isdir(os.path.join(f"{target_directory}", d))
        ]
        for af3_subdir in af3_subdirs:
            if not skip_paefigs:
                logger.info(f"Making PAE figures for {af3_subdir}")
                paeplot.plot_all_paes(
                    f"{target_directory}/{af3_subdir}",
                    "af3pae",
                    dpi=100,
                )
                paeplot.plot_best_pae(
                    f"{target_directory}/{af3_subdir}",
                    "af3pae_best",
                    dpi=100,
                )

            logger.info(f"Running ipSAE on {af3_subdir}")
            os.system(
                f"{python_binary} {ipsae_py_script} "
                f"{target_directory}/{af3_subdir}/{af3_subdir}_confidences.json "
                f"{target_directory}/{af3_subdir}/{af3_subdir}_model.cif "
                f"{pae_cutoff} {dist_cutoff}"
            )

            results = []
            ipSAE_min = float("inf")
            with open(
                f"{target_directory}/{af3_subdir}/{af3_subdir}_ipsae.json", "w"
            ) as f:
                with open(
                    f"{target_directory}/{af3_subdir}/{af3_subdir}_model_{pae_cutoff}_{dist_cutoff}.txt"
                ) as g:
                    lines = g.readlines()

                    header = lines[1].strip().split()
                    type_idx = header.index("Type")
                    ipsae_idx = header.index("ipSAE")
                    ipsae_d0chn_idx = header.index("ipSAE_d0chn")
                    ipsae_d0dom_idx = header.index("ipSAE_d0dom")
                    ipTM_af_idx = header.index("ipTM_af")
                    ipTM_d0chn_idx = header.index("ipTM_d0chn")
                    pDockQ_idx = header.index("pDockQ")
                    pDockQ2_idx = header.index("pDockQ2")
                    LIS_idx = header.index("LIS")
                    n0res_idx = header.index("n0res")
                    n0chn_idx = header.index("n0chn")
                    n0dom_idx = header.index("n0dom")
                    d0res_idx = header.index("d0res")
                    d0chn_idx = header.index("d0chn")
                    d0dom_idx = header.index("d0dom")
                    nres1_idx = header.index("nres1")
                    nres2_idx = header.index("nres2")
                    dist1_idx = header.index("dist1")
                    dist2_idx = header.index("dist2")

                    for line in lines[1:]:
                        if not line.strip():
                            continue
                        parts = line.strip().split()
                        if parts[type_idx] == "asym":
                            ipSAE_min = min(ipSAE_min, float(parts[ipsae_idx]))
                        elif parts[type_idx] == "max":
                            results.append(
                                {
                                    "ipSAE": float(parts[ipsae_idx]),
                                    "ipSAE_min": ipSAE_min,
                                    "ipSAE_d0chn": float(parts[ipsae_d0chn_idx]),
                                    "ipSAE_d0dom": float(parts[ipsae_d0dom_idx]),
                                    "ipTM_af": float(parts[ipTM_af_idx]),
                                    "ipTM_d0chn": float(parts[ipTM_d0chn_idx]),
                                    "pDockQ": float(parts[pDockQ_idx]),
                                    "pDockQ2": float(parts[pDockQ2_idx]),
                                    "LIS": float(parts[LIS_idx]),
                                    "n0res": int(parts[n0res_idx]),
                                    "n0chn": int(parts[n0chn_idx]),
                                    "n0dom": int(parts[n0dom_idx]),
                                    "d0res": float(parts[d0res_idx]),
                                    "d0chn": float(parts[d0chn_idx]),
                                    "d0dom": float(parts[d0dom_idx]),
                                    "nres1": int(parts[nres1_idx]),
                                    "nres2": int(parts[nres2_idx]),
                                    "dist1": float(parts[dist1_idx]),
                                    "dist2": float(parts[dist2_idx]),
                                }
                            )

                json.dump(results, f, indent=2)

            ipsaes.append({f"{af3_subdir}": results[0]["ipSAE"]})
            ipsae_mins.append({f"{af3_subdir}": results[0]["ipSAE_min"]})
            af3_iptms.append({f"{af3_subdir}": results[0]["ipTM_af"]})

        # Sort the results by ipSAE and ipTM
        ipsaes.sort(key=lambda x: list(x.values())[0], reverse=True)
        ipsae_mins.sort(key=lambda x: list(x.values())[0], reverse=True)
        af3_iptms.sort(key=lambda x: list(x.values())[0], reverse=True)
        out = {
            "ipSAE": ipsaes,
            "ipSAE_min": ipsae_mins,
            "ipTM": af3_iptms,
        }
        with open(f"{target_directory}/complexmetrics.json", "w") as f:
            json.dump(out, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Generate bash script to transfer files."
    )
    parser.add_argument(
        "--bgc_directory",
        metavar="BGC directory",
        default="/data2/moriwaki/BGCcomplex/merged",
        type=str,
        help="Path to root direcotory of BGCs.",
    )
    parser.add_argument(
        "--python_binary",
        metavar="Python binary",
        default="/data2/moriwaki/complexbuilder/.venv/bin/python3.12",
        type=str,
        help="Path to the Python binary.",
    )
    parser.add_argument(
        "--ipsae_py_script",
        metavar="ipSAE script",
        default="/data2/moriwaki/complexbuilder/complexbuilder/analysis/ipsae.py",
        type=str,
        help="Path to the ipSAE script.",
    )
    parser.add_argument(
        "-p",
        "--pae_cutoff",
        metavar="PAE cutoff",
        default=10,
        type=float,
        help="PAE cutoff value.",
    )
    parser.add_argument(
        "-d",
        "--dist_cutoff",
        metavar="Distance cutoff",
        default=10,
        type=float,
        help="Distance cutoff value.",
    )
    parser.add_argument(
        "-s",
        "--nstart",
        metavar="start number",
        type=int,
        help="Start number of the BGC.",
    )
    parser.add_argument(
        "-e",
        "--nend",
        metavar="end number",
        type=int,
        help="End number of the BGC.",
    )
    parser.add_argument(
        "--skip_paefigs",
        action="store_true",
        help="Skip making PAE figures.",
    )
    parser.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    args = parser.parse_args()
    bgc_directory = args.bgc_directory
    python_binary = args.python_binary
    ipsae_py_script = args.ipsae_py_script
    pae_cutoff = args.pae_cutoff
    dist_cutoff = args.dist_cutoff
    start = args.nstart
    end = args.nend
    skip_paefigs = args.skip_paefigs
    overwrite = args.overwrite
    for bgc_number in range(start, end + 1):
        run_ipsae_on_af3(
            bgc_directory=bgc_directory,
            python_binary=python_binary,
            ipsae_py_script=ipsae_py_script,
            bgc_number=bgc_number,
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff,
            skip_paefigs=skip_paefigs,
            overwrite=overwrite,
        )


if __name__ == "__main__":
    main()
