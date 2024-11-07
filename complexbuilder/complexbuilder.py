#!/usr/bin/env python3
# %%
import argparse
import os

from complexbuilder.common.log import log_setup
from complexbuilder.common.parser import classify_proteins
from complexbuilder.common.sequences import generate_multimer_input_for_colabfold

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
    "-c",
    "--clip_length",
    metavar="clip_length",
    type=int,
    default=1500,
    help="Maximum length of protein sequences to extract.",
)
log_setup(level="WARNING")
for i in range(87, 101):
    file = f"/Users/YoshitakaM/Downloads/mibig_gbk_3.1/BGC{i:07d}.gbk"
    basename = os.path.basename(file)
    if not os.path.exists(file):
        continue
    nonnrpspksproteins, nrpspksproteins = classify_proteins(
        file,
        clip_length=1500,
    )
    output = generate_multimer_input_for_colabfold(
        nonnrpspksproteins, extention="fasta", use_productname=False
    )
    with open(f"{os.path.splitext(basename)[0]}.fasta", "w") as f:
        f.write(output)

# %%
