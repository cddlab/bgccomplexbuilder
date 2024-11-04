#!/usr/bin/env python3
import argparse
from itertools import combinations_with_replacement
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from complexbuilder.common.parser import classify_proteins
from complexbuilder.common.sequences import (
    concatenate_two_sequences,
    generate_seqs_combinations,
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
    "-c",
    "--clip_length",
    metavar="clip_length",
    type=int,
    default=1500,
    help="Maximum length of protein sequences to extract.",
)
args = parser.parse_args()

nonnrpspksproteins, nrpspksproteins = classify_proteins(
    args.input, clip_length=args.clip_length
)
for seq1, seq2 in generate_seqs_combinations(nrpspksproteins):
    concatenated_seq = concatenate_two_sequences(seq1, seq2)
    print(concatenated_seq)
# combine the two lists of SeqRecord objects
