import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from complexbuilder.common.parser import classify_proteins
from complexbuilder.common.sequences import (
    concatenate_two_sequences,
    generate_multimer_input_for_colabfold,
    generate_seqs_combinations,
)

# def test_generate_seqs_combinations():
#     """test of generate_seqs_combinations"""
#     seq1 = SeqRecord(
#         Seq("MTEITAAMVKELREST"), id="seq1", description="Example sequence 1"
#     )

#     seq2 = SeqRecord(Seq("AKAIKES"), id="seq2", description="Example sequence 2")
#     seq3 = SeqRecord(Seq("TAKCLSICKSITK"), id="seq3", description="Example sequence 3")
#     seqs = [seq1, seq2, seq3]
#     generate_seqs_combinations(seqs)
#     assert len(generate_seqs_combinations(seqs)) == 6
#     test1, test2 = generate_seqs_combinations(seqs)[0]
#     assert test1.seq + test2.seq == "MTEITAAMVKELRESTMTEITAAMVKELREST"
#     concatenated_seqs: list[SeqRecord] = [
#         concatenate_two_sequences(seq1, seq2)
#         for seq1, seq2 in generate_seqs_combinations(seqs)
#     ]
#     assert concatenated_seqs[0].seq == "MTEITAAMVKELREST:MTEITAAMVKELREST"
#     assert concatenated_seqs[0].id == "seq1_seq1"


def test_BGC0000028():
    """test for BGC0000028.gbk."""
    i = 28
    file = f"/Users/YoshitakaM/Downloads/mibig_gbk_4.0/BGC{i:07d}.gbk"
    basename = os.path.basename(file)
    nonnrpspksproteins, nrpspksproteins = classify_proteins(file)
    output = generate_multimer_input_for_colabfold(
        nonnrpspksproteins, extention="fasta", use_productname=False
    )
    with open(f"{os.path.splitext(basename)[0]}.fasta", "w") as f:
        f.write("".join(output))


def test_BGC0000037():
    """test for BGC0000037.gbk. There is only 1 gene that has
    a NRPS/PKS domain. The output fasta should be 0 bytes."""
    i = 37
    file = f"/Users/YoshitakaM/Downloads/mibig_gbk_4.0/BGC{i:07d}.gbk"
    basename = os.path.basename(file)
    nonnrpspksproteins, nrpspksproteins = classify_proteins(file)
    output = generate_multimer_input_for_colabfold(
        nonnrpspksproteins, extention="fasta", use_productname=False
    )
    with open(f"{os.path.splitext(basename)[0]}.fasta", "w") as f:
        f.write("".join(output))


def test_BGC0000053():
    """test for BGC0000053.gbk. There is no product and protein_id.
    Falls back to locus_tag."""
    i = 53
    file = f"/Users/YoshitakaM/Downloads/mibig_gbk_4.0/BGC{i:07d}.gbk"
    nonnrpspksproteins, nrpspksproteins = classify_proteins(file)
    _ = generate_multimer_input_for_colabfold(
        nonnrpspksproteins, extention="fasta", use_productname=False
    )


def test_BGC0000087():
    """test for BGC0000087.gbk. Some proteins do not have protein id and locus tag."""
    i = 87
    file = f"/Users/YoshitakaM/Downloads/mibig_gbk_4.0/BGC{i:07d}.gbk"
    nonnrpspksproteins, nrpspksproteins = classify_proteins(file)
    _ = generate_multimer_input_for_colabfold(
        nonnrpspksproteins, extention="fasta", use_productname=False
    )
