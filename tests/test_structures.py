import os

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from complexbuilder.common.parser import classify_proteins
from complexbuilder.common.sequences import (
    concatenate_two_sequences,
    generate_multimer_input_for_colabfold,
    generate_seqs_combinations,
    get_protein_sequence_from_uniprot,
)


def test_generate_seqs_combinations():
    """test of generate_seqs_combinations"""
    seq1 = SeqRecord(
        Seq("MTEITAAMVKELREST"), id="seq1", description="Example sequence 1"
    )

    seq2 = SeqRecord(Seq("AKAIKES"), id="seq2", description="Example sequence 2")
    seq3 = SeqRecord(Seq("TAKCLSICKSITK"), id="seq3", description="Example sequence 3")
    seqs = [seq1, seq2, seq3]
    generate_seqs_combinations(seqs)
    assert len(generate_seqs_combinations(seqs)) == 6
    test1, test2 = generate_seqs_combinations(seqs)[0]
    assert test1.seq + test2.seq == "MTEITAAMVKELRESTMTEITAAMVKELREST"
    concatenated_seqs: list[SeqRecord] = [
        concatenate_two_sequences(seq1, seq2)
        for seq1, seq2 in generate_seqs_combinations(seqs)
    ]
    assert concatenated_seqs[0].seq == "MTEITAAMVKELREST:MTEITAAMVKELREST"
    assert concatenated_seqs[0].id == "seq1_seq1"
