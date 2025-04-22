#!/usr/bin/env python3
# %%
import argparse
import json
import os
import string
from pathlib import Path

from loguru import logger


def sanitised_name(name) -> str:
    """Returns sanitised version of the name that can be used as a filename."""
    lower_spaceless_name = name.lower().replace(" ", "_")
    allowed_chars = set(string.ascii_lowercase + string.digits + "_-.")
    return "".join(char for char in lower_spaceless_name if char in allowed_chars)


def get_total_length(jsonfile: Path) -> int:
    total_len = 0
    with open(jsonfile) as f:
        data = json.load(f)
    for i in range(len(data["sequences"])):
        id = data["sequences"][i]["protein"]["id"]
        seq_len = len(data["sequences"][i]["protein"]["sequence"])
        total_len += len(id) * seq_len
    return total_len


def run_transfer_and_return_script(
    bgccomplexpath,
    outputfile,
    returnfile_tsubame,
    returnfile_lento,
    nstart,
    nend,
    tsubame_cargo_name="tsubame_cargo",
    lento_cargo_name="lento_cargo",
) -> None:
    """generate a script to transfer the af3 files to tsubame and return the af3 files
    to the original location
    Args:
        bgccomplexpath: path to the BGCcomplex directory
        outputfile: path to the output file
        returnfile: path to the return file
        nstart: start number of the BGC
        nend: end number of the BGC
        tsubame_cargo_name: name of the tsubame cargo directory
        lento_cargo_name: name of the lento cargo directory
    """
    os.remove(outputfile) if os.path.exists(outputfile) else None
    os.remove(returnfile_tsubame) if os.path.exists(returnfile_tsubame) else None
    os.remove(returnfile_lento) if os.path.exists(returnfile_lento) else None
    with open(outputfile, "w") as f:
        f.write(f"mkdir -p {bgccomplexpath}/{tsubame_cargo_name}\n")
        f.write(f"mkdir -p {bgccomplexpath}/{lento_cargo_name}\n")

    for i in range(nstart, nend + 1):
        # e.g. BGCcomplex/BGC0000001/af3
        dirpath = Path(bgccomplexpath) / f"BGC{str(i).zfill(7)}_rest"
        if not dirpath.exists():
            logger.warning(f"{dirpath} does not exist.")
            continue
        else:
            af3_dirpath = dirpath / "af3"
            jsonfiles = af3_dirpath.glob("*.json")
            for jsonfile in jsonfiles:
                basename = jsonfile.stem
                if not os.path.exists(
                    af3_dirpath / sanitised_name(basename) / "TERMS_OF_USE.md"
                ):
                    total_len = get_total_length(jsonfile)
                    logger.info(f"Processing {basename}, total length: {total_len}")
                    with open(outputfile, "a") as f:
                        if total_len > 2000:
                            os.makedirs(
                                f"{bgccomplexpath}/{tsubame_cargo_name}", exist_ok=True
                            )
                            cmd = f"ln -sf {jsonfile} {bgccomplexpath}/{tsubame_cargo_name}"
                        else:
                            os.makedirs(
                                f"{bgccomplexpath}/{lento_cargo_name}", exist_ok=True
                            )
                            cmd = (
                                f"ln -sf {jsonfile} {bgccomplexpath}/{lento_cargo_name}"
                            )
                        f.write(cmd + "\n")
                    if total_len > 2000:
                        with open(returnfile_tsubame, "a") as g:
                            # get the BGC name, e.g. BGC0000136
                            bgc_name = jsonfile.parents[1].name
                            cmd = (
                                f"mv {bgccomplexpath}/{tsubame_cargo_name}/"
                                f"{sanitised_name(basename)} "
                            )
                            cmd += f"{bgccomplexpath}/{bgc_name}/af3\n"
                            g.write(cmd)
                    else:
                        with open(returnfile_lento, "a") as h:
                            bgc_name = jsonfile.parents[1].name
                            cmd = (
                                f"rm -rf {bgccomplexpath}/{bgc_name}/af3/"
                                f"{sanitised_name(basename)}\n"
                                f"mv {bgccomplexpath}/{lento_cargo_name}/"
                                f"{sanitised_name(basename)} "
                            )
                            cmd += f"{bgccomplexpath}/{bgc_name}/af3\n"
                            h.write(cmd)
                    # with open(returnfile, "a") as g:
                    #     # get the BGC name, e.g. BGC0000136
                    #     bgc_name = jsonfile.parents[1].name
                    #     if total_len > 2000:
                    #         cmd = (
                    #             f"mv {bgccomplexpath}/{tsubame_cargo_name}/"
                    #             f"{sanitised_name(basename)} "
                    #         )
                    #         cmd += f"{bgccomplexpath}/{bgc_name}/af3\n"
                    #     else:
                    #         cmd = (
                    #             f"rm -rf {bgccomplexpath}/{bgc_name}/af3/"
                    #             f"{sanitised_name(basename)}\n"
                    #             f"mv {bgccomplexpath}/{lento_cargo_name}/"
                    #             f"{sanitised_name(basename)} "
                    #             f"{bgccomplexpath}/{bgc_name}/af3\n"
                    #         )
                    #     g.write(cmd)


def main():
    parser = argparse.ArgumentParser(
        description="Generate bash script to transfer files."
    )
    parser.add_argument(
        "-c",
        "--bgccomplexpath",
        metavar="BGCcomplex path",
        type=str,
        help="Path to the input GenBank file containing mibig BGC data.",
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
        "-t",
        "--tsubame_cargo_name",
        metavar="tsubame cargo name",
        type=str,
        default="tsubame_cargo5",
        help="Name of the tsubame cargo directory.",
    )
    parser.add_argument(
        "-l",
        "--lento_cargo_name",
        metavar="lento cargo name",
        type=str,
        default="lento_cargo5",
        help="Name of the lento cargo directory.",
    )
    args = parser.parse_args()

    # bgccomplexpath = args.bgccomplexpath
    bgccomplexpath = "/data2/moriwaki/BGCcomplex/restfastas"
    returnpath = "/data2/moriwaki/complexbuilder"
    nstart = args.nstart
    nend = args.nend
    outputfile = f"{returnpath}/af3_transfer{nstart}_{nend}.sh"
    returnfile_tsubame = f"{returnpath}/af3_return{nstart}_{nend}_tsubame.sh"
    returnfile_lento = f"{returnpath}/af3_return{nstart}_{nend}_lento.sh"
    run_transfer_and_return_script(
        bgccomplexpath,
        outputfile,
        returnfile_tsubame,
        returnfile_lento,
        nstart,
        nend,
        args.tsubame_cargo_name,
        args.lento_cargo_name,
    )


if __name__ == "__main__":
    main()
