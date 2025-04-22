#!/usr/bin/env python3
# %%
import os
import shutil
import string

from Bio import SeqIO
from loguru import logger


def sanitised_name(name) -> str:
    """Returns sanitised version of the name that can be used as a filename."""
    lower_spaceless_name = name.lower().replace(" ", "_")
    allowed_chars = set(string.ascii_lowercase + string.digits + "_-.")
    return "".join(char for char in lower_spaceless_name if char in allowed_chars)


def make_list_of_directories(fastafile: str) -> list:
    """Returns a list of directories in the given path.
    Args:
        fastafile (str): The path to the FASTA file.
    Returns:
        list: A list of directories to be tranferred in the given path.
    """
    with open(fastafile, "r") as file:
        records = SeqIO.parse(file, "fasta")
        dirs_to_be_transferred = []
        for record in records:
            name = sanitised_name(record.id)
            if name not in dirs_to_be_transferred:
                dirs_to_be_transferred.append(name)
    return dirs_to_be_transferred


def transfer_from_predicted_dir(
    fastafile: str, predicted_dir: str, restfiles_dir: str, distination: str
) -> None:
    """Transfer matched directories from the predicted directory to the destination.
    Args:
        fastafile (str): The path to the FASTA file.
        predicted_dir (str): The directory containing the predicted directories.
        restfiles_dir (str): The directory containing the restfiles.
        distination (str): The destination directory.
    """
    dirs_to_be_transferred = make_list_of_directories(fastafile)
    # logger.info(f"Num of directories to be transferred: {len(dirs_to_be_transferred)}")
    # get BGC number from the fasta file
    # e.g. "/data2/moriwaki/BGCcomplex/newfastas/newBGC0000001.fasta" -> BGC0000001
    bgc_number = os.path.basename(fastafile).split(".")[0].replace("new", "")
    os.makedirs(distination, exist_ok=True)
    os.makedirs(
        os.path.join(distination, bgc_number),
        exist_ok=True,
    )
    for dir in dirs_to_be_transferred:
        if os.path.isdir(os.path.join(predicted_dir, bgc_number, dir)):
            shutil.move(
                os.path.join(predicted_dir, bgc_number, dir),
                os.path.join(distination, bgc_number),
            )
        elif os.path.isdir(
            os.path.join(restfiles_dir, f"{bgc_number}_rest", "af3", dir)
        ):
            shutil.move(
                os.path.join(restfiles_dir, f"{bgc_number}_rest", "af3", dir),
                os.path.join(distination, bgc_number),
            )
        else:
            print(
                f"Directory {dir} not found in {predicted_dir}/{bgc_number} "
                f"or {restfiles_dir}/{bgc_number}_rest. Skipping."
            )
    # count directories in the destination
    # dirs_in_destination = [
    #     d
    #     for d in os.listdir(os.path.join(distination, bgc_number))
    #     if os.path.isdir(os.path.join(distination, bgc_number, d))
    # ]
    # logger.info(f"Num of directories in destination: {len(dirs_in_destination)}")


# %%
start = 1
end = 1
fastadir = "/data2/moriwaki/BGCcomplex/newfastas"
for i in range(start, end + 1):
    fastafile = os.path.join(fastadir, f"newBGC{str(i).zfill(7)}.fasta")
    if os.path.exists(fastafile):
        predicted_dir = "/data2/moriwaki/BGCcomplex/predicted"
        restfiles_dir = "/data2/moriwaki/BGCcomplex/restfastas"
        distination = "/data2/moriwaki/BGCcomplex/merged"
        transfer_from_predicted_dir(
            fastafile, predicted_dir, restfiles_dir, distination
        )
    else:
        print(f"File {fastafile} does not exist. Skipping.")
# %%
# transfer_from_predicted_dir(
#     "/data2/moriwaki/BGCcomplex/newfastas/newBGC0000001.fasta",
#     "/data2/moriwaki/BGCcomplex/predicted",
#     "/data2/moriwaki/BGCcomplex/restfastas",
#     "/data2/moriwaki/BGCcomplex/merged",
# )
# %%
