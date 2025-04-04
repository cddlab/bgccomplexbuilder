#!/usr/bin/env python3
# %%
import json
import os
import sys

from alphafold3tools import paeplot


def run_ipsae_on_af3(
    bgc_directory: str,
    python_binary: str,
    ipsae_py_script: str,
    bgc_number: int,
    pae_cutoff: float = 10,
    dist_cutoff: float = 10,
    ipsae_cutoff: float = 0.5,
    is_overwrite: bool = False,
) -> None:
    target_directory = f"{bgc_directory}/BGC{bgc_number:07d}"
    if not os.path.exists(f"{target_directory}/af3"):
        print(f"af3 directory not found in {target_directory}")
    elif os.path.exists(f"{target_directory}/af3/ipsaeDONE.txt") and not is_overwrite:
        print(f"ipSAE already ran on {target_directory}")
    else:
        hit_complexes = []
        af3_subdirs = [
            d
            for d in os.listdir(f"{target_directory}/af3")
            if os.path.isdir(os.path.join(f"{target_directory}/af3", d))
        ]
        for af3_subdir in af3_subdirs:
            paeplot.plot_all_paes(
                f"{target_directory}/af3/{af3_subdir}",
                "af3pae",
                dpi=200,
            )
            paeplot.plot_best_pae(
                f"{target_directory}/af3/{af3_subdir}",
                "af3pae_best",
                dpi=200,
            )
            print(f"Running ipSAE on {af3_subdir}")
            os.system(
                f"{python_binary} {ipsae_py_script} "
                f"{target_directory}/af3/{af3_subdir}/{af3_subdir}_confidences.json "
                f"{target_directory}/af3/{af3_subdir}/{af3_subdir}_model.cif "
                f"{pae_cutoff} {dist_cutoff}"
            )

            results = []

            with open(
                f"{target_directory}/af3/{af3_subdir}/{af3_subdir}_ipsae.json", "w"
            ) as f:
                with open(
                    f"{target_directory}/af3/{af3_subdir}/{af3_subdir}_model_{pae_cutoff}_{dist_cutoff}.txt"
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
                        if parts[type_idx] == "max":
                            results.append(
                                {
                                    "ipSAE": float(parts[ipsae_idx]),
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

                json.dump(results, f, indent=4)
            if results[0]["ipSAE"] > ipsae_cutoff:
                hit_complexes.append((af3_subdir, results[0]["ipSAE"]))
        # put DONE marker
        with open(f"{target_directory}/af3/ipsaeDONE.txt", "w") as f:
            f.write(
                f"Hit complexes (ipSAE > {ipsae_cutoff})\n"
                + "\n".join([f"{c[0]}: {c[1]}" for c in hit_complexes])
            )


# %%


ipsae_py_script = (
    "/Volumes/MacintoshHD/workdir/complexbuilder/complexbuilder/common/ipsae.py"
)
bgc_directory = "/Users/YoshitakaM/Desktop/BGC_heteromer"
python_binary = "/Volumes/MacintoshHD/workdir/complexbuilder/.venv/bin/python3.12"
pae_cutoff = 10
dist_cutoff = 10
ipsae_cutoff = 0.5
start = int(sys.argv[1])
end = int(sys.argv[2])

for bgc_number in range(start, end + 1):
    run_ipsae_on_af3(
        bgc_directory=bgc_directory,
        python_binary=python_binary,
        ipsae_py_script=ipsae_py_script,
        bgc_number=bgc_number,
        pae_cutoff=pae_cutoff,
        dist_cutoff=dist_cutoff,
        ipsae_cutoff=ipsae_cutoff,
        is_overwrite=True,
    )
