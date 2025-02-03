# %%
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger


def classify_proteins(
    genbank_file: str | Path, max_length: int, decompose_nrpspks: bool = False
) -> tuple[list[SeqRecord], list[SeqRecord]]:
    """Collect protein sequences from GenBank file and return a list of
    SeqRecord objects.
    The structural domains of the protein belonging to NRPS_PKS are
    returned together as a separate list.

    Args:
        - genbank_file (str | Path): Path to the GenBank file.
        - max_length (int): Maximum length of protein sequences
    Returns:
        - nonnrpspksproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that do not belong to NRPS_PKS.
        - nrpspksproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that belong to NRPS_PKS.
        - decompose_nrpspks (bool): If True, decompose the NRPS_PKS proteins
            into individual domains.
    """

    nonnrpspksproteins: list[SeqRecord] = []
    nrpspksproteins: list[SeqRecord] = []
    genbank_file = Path(genbank_file)
    for record in SeqIO.parse(genbank_file, "genbank"):
        logger.info(f"Record ID: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                # protein_id and aa_sequence are required fields
                aa_sequence = feature.qualifiers.get("translation")[0]
                if "protein_id" in feature.qualifiers:
                    protein_id = feature.qualifiers.get("protein_id")[0]
                elif "locus_tag" in feature.qualifiers:
                    # fallback to locus_tag if protein_id is not available
                    protein_id = feature.qualifiers.get("locus_tag")[0]
                else:
                    logger.warning(
                        f"No protein ID and locus tag found in {genbank_file} . "
                        "Use gene ID."
                    )
                    protein_id = feature.qualifiers.get("gene")[0]
                if "product" in feature.qualifiers:
                    product = feature.qualifiers.get("product")[0]
                else:
                    product = ""
                # truncate the protein sequence if it exceeds <clip_length> residues
                if "NRPS_PKS" in feature.qualifiers:
                    pattern = r"Domain: \S+ \((\d+)-(\d+)\).*nrpspksdomains_(\S+)"
                    if len(aa_sequence) < max_length:
                        nrpspksproteins.append(
                            SeqRecord(
                                Seq(aa_sequence),
                                id=protein_id,
                                description=product,
                            )
                        )
                    elif decompose_nrpspks:
                        for idx in range(len(feature.qualifiers["NRPS_PKS"])):
                            match = re.search(
                                pattern, feature.qualifiers["NRPS_PKS"][idx]
                            )
                            if match:
                                start_res = int(match.group(1))
                                end_res = int(match.group(2))
                                domainname = match.group(3)
                                nrpspksproteins.append(
                                    SeqRecord(
                                        Seq(aa_sequence[start_res:end_res]),
                                        id=f"{protein_id}@{idx}",
                                        description=domainname,
                                    )
                                )
                            else:
                                raise ValueError(
                                    "Domain residue regions not found for Record ID: "
                                    f"{record.id}"
                                )
                else:
                    nonnrpspksproteins.append(
                        SeqRecord(Seq(aa_sequence), id=protein_id, description=product)
                    )
    if nonnrpspksproteins == []:
        logger.warning(
            "No non-NRPS/PKS protein sequences found in "
            f"the GenBank file {genbank_file} ."
        )
    if nrpspksproteins == []:
        logger.warning(
            f"No NRPS/PKS domains found in the GenBank file {genbank_file} ."
        )
    logger.info(f"Number of non-NRPS/PKS proteins: {len(nonnrpspksproteins)}")
    logger.info(f"Number of NRPS/PKS proteins: {len(nrpspksproteins)}")
    return nonnrpspksproteins, nrpspksproteins


# %%
