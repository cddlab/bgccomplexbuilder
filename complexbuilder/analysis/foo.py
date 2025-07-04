# %%
from Bio import SeqIO

input_filename = "/Users/YoshitakaM/Desktop/IPR029058.fasta"
output_filename = "/Users/YoshitakaM/Desktop/IPR029058_modified.fasta"

with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        # ヘッダーIDを'|'で分割し、最初の部分を取得します
        new_id = record.id.split("|")[0]

        # recordのIDとdescriptionを更新します
        record.id = new_id
        record.description = (
            ""  # descriptionを空に設定し、ヘッダーにIDのみが含まれるようにします
        )

        # 更新されたrecordをファイルに書き込みます
        SeqIO.write(record, outfile, "fasta")

# %%
