# %%
import json
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def classify_proteins(
    genbank_file: str | Path, clip_length: int = 1500
) -> tuple[list[SeqRecord], list[SeqRecord]]:
    """Collect protein sequences from GenBank file and return a list of
    SeqRecord objects.
    The structural domains of the protein belonging to NRPS_PKS are
    returned together as a separate list.

    Args:
        - genbank_file (str | Path): Path to the GenBank file.
        - clip_length (int): The maximum length of the protein sequence to truncate.
            Default is 1500.
    returns:
        - nonnrpspksproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that do not belong to NRPS_PKS.
        - nrpspksproteins (list[SeqRecord]): List of SeqRecord objects
            for proteins that belong to NRPS_PKS.
    """

    nonnrpspksproteins: list[SeqRecord] = []
    nrpspksproteins: list[SeqRecord] = []
    genbank_file = Path(genbank_file)
    for record in SeqIO.parse(genbank_file, "genbank"):
        print(f"Record ID: {record.id}")
        for feature in record.features:
            if feature.type == "CDS":
                # protein_id and aa_sequence are required fields
                protein_id = feature.qualifiers.get("protein_id")[0]
                aa_sequence = feature.qualifiers.get("translation")[0]
                name = feature.qualifiers.get("product")[0]
                aa_len = len(aa_sequence)
                # truncate the protein sequence if it exceeds <clip_length> residues
                if "NRPS_PKS" in feature.qualifiers and aa_len > clip_length:
                    pattern = r"Domain: \S+ \((\d+)-(\d+)\).*nrpspksdomains_(\S+)"
                    for idx in range(len(feature.qualifiers["NRPS_PKS"])):
                        match = re.search(pattern, feature.qualifiers["NRPS_PKS"][idx])
                        if match:
                            start_res = int(match.group(1))
                            end_res = int(match.group(2))
                            domainname = match.group(3)
                            nrpspksproteins.append(
                                SeqRecord(
                                    Seq(aa_sequence[start_res:end_res]),
                                    id=protein_id,
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
                        SeqRecord(Seq(aa_sequence), id=protein_id, description=name)
                    )
    return nonnrpspksproteins, nrpspksproteins


# %%
genbank_file = "/Users/YoshitakaM/Downloads/mibig_gbk_3.1/BGC0000028.gbk"
nonnrpspksproteins, nrpspksproteins = classify_proteins(genbank_file)
for nrpspksprotein in nrpspksproteins:
    print(nrpspksprotein.id, nrpspksprotein.description, len(nrpspksprotein.seq))
# %%
# json_directory = "/Users/YoshitakaM/Downloads/mibig_json_minimal"
json_directory = "/Users/YoshitakaM/Downloads/mibig_json_3.1"


# すべてのjsonファイルをmibig_dataのdictに変換
def parse_mibig_json(mibig_json: str) -> dict:
    """Parse a MIBiG JSON file and return a dictionary with the parsed data."""
    with open(mibig_json, "r") as f:
        mibig_data = json.load(f)
    return mibig_data


# %%
