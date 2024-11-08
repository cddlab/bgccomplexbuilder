#!/usr/bin/env python3
# %%
import argparse
import os

from loguru import logger

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


# %%
