import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_multiple_msa(seqs: dict, copies: int | tuple) -> SeqRecord:
    """
    Make a multiple sequence alignment (MSA) of proteins from a dictionary of sequences
    for ColabFold. The sequences are concatenated with a colon and the number of copies.
    inputs:
    - seqs: dict, dictionary of sequences
    - copies: int | tuple, number of copies of each sequence
    returns:
    - record: SeqRecord, a sequence record
    Example:
    seqs = {"seq1": "MTEITAAMVKELREST", "seq2": "AKAIKES"}
    copies = (1, 2)
    make_multiple_msa(seqs, copies) -> SeqRecord("MTEITAAMVKELREST:AKAIKES:AKAIKES")
    """
    # convert int to tuple
    if isinstance(copies, int):
        copies = (copies,)
    # check if the number of copies is the same as the number of sequences
    if len(copies) != len(seqs):
        raise ValueError(
            "The number of copies must be the same as the number of sequences"
        )

    out_header = "_".join(seqs.keys())
    out_sequence = ""
    for seq, copy in zip(seqs.values(), copies, strict=True):
        for _ in range(copy):
            out_sequence += seq + ":"
    # remove the last colon
    out_sequence = out_sequence[:-1]

    record = SeqRecord(Seq(out_sequence), id=out_header, description="")
    return record


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
