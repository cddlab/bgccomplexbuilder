# %%
import json
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def combination_of_sequence(genbank_file: str | Path, mibig_json: str | Path):
    """Combine protein records from GenBank to return protein pairs
    for multimer complex building.

    There are seven "biosyn_class" values (Secondary Metabolite record type)
    in MIBiG JSON files:
      1. Polyketide
      2. NRP
      3. RiPP
      4. Terpene
      5. Saccharide
      6. Alkaloid
      7. Other

    inputs:
    - genbank_file (str | Path): Path to the GenBank file.
    - mibig_json (str | Path): Path to the MIBiG JSON file.
    returns:
    - protein_pairs: List of protein pairs.
    """
    pass


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
# json_directory = "/Users/YoshitakaM/Downloads/mibig_json_minimal"
json_directory = "/Users/YoshitakaM/Downloads/mibig_json_3.1"
# すべてのjsonファイルをmibig_dataのdictに変換
mibig_data = {}
for json_file in Path(json_directory).glob("*.json"):
    with open(json_file, "r") as f:
        mibig_data[json_file.stem] = json.load(f)

# %%
# mibig_data[*]["cluster"]["biosyn_class"]の値をすべて表示
all_biosyn_classes = set()
for _, value in mibig_data.items():
    biosyn_classes = value["cluster"]["biosyn_class"]
    all_biosyn_classes.update(biosyn_classes)

for biosyn_class in all_biosyn_classes:
    print(biosyn_class)
# %%
