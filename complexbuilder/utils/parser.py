# %%
import json
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def parse_mibig_json(mibig_json: str) -> dict:
    """Parse a MIBiG JSON file and return a dictionary with the parsed data."""
    with open(mibig_json, "r") as f:
        mibig_data = json.load(f)
    return mibig_data


# %%
genbank_file = "/Users/YoshitakaM/Downloads/mibig_gbk_3.1/BGC0000028.gbk"
for record in SeqIO.parse(genbank_file, "genbank"):
    print(f"Record ID: {record.id}")

    for feature in record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers:
                translation = feature.qualifiers["translation"][0]
                print(f"Protein translation: {translation}")
            if "NRPS_PKS" in feature.qualifiers:
                pattern = r"Domain: \S+ \((\d+)-(\d+)\)\."
                for idx in range(len(feature.qualifiers["NRPS_PKS"])):
                    match = re.search(pattern, feature.qualifiers["NRPS_PKS"][idx])
                    if match:
                        start_res = int(match.group(1))
                        ending_res = int(match.group(2))
                        print(f"Domain residue region: {start_res}-{ending_res}")
                    else:
                        raise ValueError(
                            f"Domain residue regions not found for Record ID: {record.id}"
                        )

            gene = feature.qualifiers.get("gene", ["N/A"])[0]
            product = feature.qualifiers.get("product", ["N/A"])[0]
            print(f"Gene: {gene}, Product: {product}\n")
# %%
