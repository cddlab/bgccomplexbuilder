from itertools import combinations_with_replacement, product

import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_seqs_combinations(
    seqs: tuple[list[SeqRecord], ...],
) -> list[tuple]:
    """
    Generate all possible combinations of two elements from the list,
    allowing for duplicates but not considering (A, B) and (B, A) as different.

    Args:
        - seqs: SeqRecord, the list of elements to combine. If seqs has only one element,
                it is considered as a list of SeqRecord objects.

    Returns:
        - combinations: list of tuples, each tuple contains a pair of elements
    """
    if len(seqs) == 1:
        return list(combinations_with_replacement(seqs[0], 2))
    elif len(seqs) == 2:
        nonnrpspksproteins, nrpspksproteins = seqs
        return list(product(nonnrpspksproteins, nrpspksproteins))
    else:
        raise ValueError("Only two lists of SeqRecord objects are allowed.")


def concatenate_two_sequences(seq1: SeqRecord, seq2: SeqRecord) -> SeqRecord:
    """
    Concatenate the sequences of two SeqRecord objects with a colon.

    Args:
        - seq1: SeqRecord, the first sequence record
        - seq2: SeqRecord, the second sequence record
    Returns:
        - concat_sequence: SeqRecord, a concatenated sequence record with a colon
    Example:
        seq1 = SeqRecord(Seq("MTEITAAMVKELREST"), id="seq1", description="AAA")
        seq2 = SeqRecord(Seq("AKAIKES"), id="seq2", description="BBB")
        concatenate_two_sequences(seq1, seq2) ->
        SeqRecord(Seq("MTEITAAMVKELREST:AKAIKES"),
                  id="seq1_seq2", description="AAA_BBB")
    """
    if not isinstance(seq1, SeqRecord) or not isinstance(seq2, SeqRecord):
        raise ValueError("Both inputs must be SeqRecord objects")
    concat_sequence = f"{seq1.seq}:{seq2.seq}"
    return SeqRecord(
        Seq(concat_sequence),
        id=f"{seq1.id}_{seq2.id}",
        description=f"{seq1.description}_{seq2.description}",
    )


def generate_multimer_input_for_colabfold(
    *proteins: list[SeqRecord],
    extention: str = "csv",
    use_productname: bool = False,
) -> list[str]:
    """
    Generate a string for input of ColabFold from a list of SeqRecord objects.


    Args:
        - proteins: list[SeqRecord], the list of SeqRecord objects
        - extention: str, "csv" or "fasta" are only allowed. Default is "csv"
        - use_productname: bool, if True, use product name defined as the description.
    """
    if extention not in ["csv", "fasta"]:
        raise ValueError("Only 'csv' and 'fasta' are allowed for the extention.")

    concat_seqs = []
    if len(proteins) != 1 and len(proteins) != 2:
        raise ValueError("Only one or two lists of SeqRecord objects are allowed.")
    for seq1, seq2 in generate_seqs_combinations(proteins):
        concat_seqrecord = concatenate_two_sequences(seq1, seq2)
        if use_productname:
            id = concat_seqrecord.description
        else:
            id = concat_seqrecord.id
        if extention == "csv":
            concat_seqs.append(f"{id},{concat_seqrecord.seq}\n")
        elif extention == "fasta":
            concat_seqs.append(f">{id}\n{concat_seqrecord.seq}\n")
    return concat_seqs


def concatenate_chunks(chunks: list[str], max_bytes: int = 500000) -> list[str]:
    """
    Concatenate the chunks of sequences but not exceeding the maximum bytes.
    Args:
        - chunks: list[str], the list of sequence chunks
        - max_bytes: int, the maximum bytes of the concatenated sequence
    Returns:
        - concatenated_chunks: list[str], the list of concatenated sequence chunks
    """
    concatenated_chunks = []
    current_chunk = ""
    for chunk in chunks:
        if len(current_chunk) + len(chunk) <= max_bytes:
            current_chunk += chunk
        else:
            concatenated_chunks.append(current_chunk)
            current_chunk = chunk
    concatenated_chunks.append(current_chunk)
    return concatenated_chunks
