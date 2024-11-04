from itertools import combinations_with_replacement

import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_seqs_combinations(seqs: list[SeqRecord]) -> list[tuple]:
    """
    Generate all possible combinations of two elements from the list,
    allowing for duplicates but not considering (A, B) and (B, A) as different.

    Args:
        - seqs: SeqRecord, the list of elements to combine

    Returns:
        - combinations: list of tuples, each tuple contains a pair of elements
    """
    return list(combinations_with_replacement(seqs, 2))


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
    concat_sequence = seq1.seq + ":" + seq2.seq
    return SeqRecord(
        Seq(concat_sequence),
        id=f"{seq1.id}_{seq2.id}",
        description=f"{seq1.description}_{seq2.description}",
    )


def generate_multimer_input_for_colabfold(
    seqs: list[SeqRecord], extention: str = "csv"
) -> str | ValueError:
    """
    Generate a string for input of ColabFold from a list of SeqRecord objects.

    Args:
        - seqs: list[SeqRecord], the list of SeqRecord objects
        - extention: str, "csv" or "fasta" are only allowed. Default is "csv"
    """
    output: str = ""
    if extention == "csv":
        for seq1, seq2 in generate_seqs_combinations(seqs):
            concat_seqrecord = concatenate_two_sequences(seq1, seq2)
            output += f"{concat_seqrecord.id},{concat_seqrecord.seq}\n"
        return output

    elif extention == "fasta":
        for seq1, seq2 in generate_seqs_combinations(seqs):
            concat_seqrecord = concatenate_two_sequences(seq1, seq2)
            output += f">{concat_seqrecord.id}\n{concat_seqrecord.seq}\n"
        return output
    else:
        return ValueError("The extention must be 'csv' or 'fasta'.")


def get_protein_sequence_from_uniprot(uniprot_id: str) -> str:
    """Retrieve the amino acid sequence from UniProt using a given UniProt ID."""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"

    response = requests.get(url)

    if response.status_code != 200:
        raise ValueError(
            f"Failed to retrieve data for {uniprot_id}."
            "HTTP Status: {response.status_code}"
        )

    fasta_data = response.text
    sequence = "".join(
        line.strip() for line in fasta_data.splitlines() if not line.startswith(">")
    )

    return sequence
