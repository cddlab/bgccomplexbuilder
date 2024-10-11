#!/usr/bin/env python3

import argparse
import os

from loguru import logger

from complexbuilder.common.log import log_setup
from complexbuilder.common.parser import classify_proteins2, make_seqrecord_from_fasta
from complexbuilder.common.sequences import (
    concatenate_chunks,
    generate_multimer_input_for_colabfold,
)

parser = argparse.ArgumentParser(
    description="Extract protein sequences from GenBank files."
)
parser.add_argument(
    "-i",
    "--input",
    metavar="input genbank_file",
    type=str,
    help="Path to the input GenBank file containing mibig BGC data.",
)
parser.add_argument(
    "-l",
    "--max_length",
    metavar="maximum length",
    type=int,
    default=1950,
    help="Maximum length of protein sequences to extract.",
)
parser.add_argument(
    "-m",
    "--maxbytes",
    metavar="maximum bytes",
    type=int,
    default=600000,
    help="Maximum bytes of sequences in a single file.",
)
parser.add_argument(
    "--decompose_nrpspks",
    action="store_true",
    help="If True, decompose the NRPS_PKS proteins into individual domains.",
)
parser.add_argument(
    "--use_productname",
    action="store_true",
    help="If True, use the product name as the protein description.",
)
parser.add_argument("--start", type=int, default=1, help="Start index of the BGC ID.")
parser.add_argument("--end", type=int, default=1, help="End index of the BGC ID.")
args = parser.parse_args(
    [
        "--max_length",
        "1950",
        "--maxbytes",
        "6000000",
        "--start",
        "1",
        "--end",
        "2826",
    ]
)

log_setup(level="WARNING")
for i in range(args.start, args.end + 1):
    file = f"/path/to/Downloads/mibig_gbk_4.0/BGC{i:07d}.gbk"
    basename = os.path.basename(file)
    if not os.path.exists(file):
        logger.warning(f"File {file} not found.")
        continue
    proteins = classify_proteins2(file, args.max_length, decompose_nrpspks=False)
    if len(proteins) > 0:
        protein_chunks = generate_multimer_input_for_colabfold(
            proteins,
            extention="fasta",
            use_productname=args.use_productname,
        )
        split_chunks = concatenate_chunks(protein_chunks, args.maxbytes)
        if len(split_chunks) == 1:
            with open(f"new{os.path.splitext(basename)[0]}.fasta", "w") as f:
                f.write(split_chunks[0])
        else:
            for i, chunk in enumerate(split_chunks):
                num = i + 1
                with open(f"new{os.path.splitext(basename)[0]}_{i}.fasta", "w") as f:
                    f.write(chunk)
