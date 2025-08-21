#!/usr/bin/env python3
# %%
import json
from pathlib import Path

gbk_directory = "/Users/YoshitakaM/Downloads/mibig_json_4.0"
mibigjsonfiles = list(Path(gbk_directory).rglob("*.json"))
# json_filesについて"status"が"active"のものを抽出
# json_filesについて"status"が"retired"のものを抽出
active_json_files = []
retired_json_files = []

for mibigjson in mibigjsonfiles:
    with open(mibigjson, "r", encoding="utf-8") as f:
        data = json.load(f)
        if data.get("status") == "active":
            # 拡張子を抜いたファイル名を取得
            active_json_files.append(mibigjson.stem)
        elif data.get("status") == "retired":
            # 拡張子を抜いたファイル名を取得
            retired_json_files.append(mibigjson.stem)


for mibigjson in mibigjsonfiles:
    with open(mibigjson, "r", encoding="utf-8") as f:
        data = json.load(f)
        if data.get("status") == "active":
            # 拡張子を抜いたファイル名を取得
            active_json_files.append(mibigjson.stem)
# %%
svg3_directory = "/Users/YoshitakaM/Desktop/svg3"
svg3_files = list(Path(svg3_directory).rglob("*.svg"))
svg3_file_stems = [svg3_file.stem for svg3_file in svg3_files]
# active_json_files側に存在して、svg3_file_stems側に存在しないものを抽出
missing_svg_files = [file for file in active_json_files if file not in svg3_file_stems]
# svg3_file_stems側に存在して、active_json_files側に存在しないものを抽出
extra_svg_files = [file for file in svg3_file_stems if file not in active_json_files]
# retired_json_files側に存在して、svg3_file_stems側に存在しているものを抽出
retired_svg_files = [file for file in retired_json_files if file in svg3_file_stems]
# %%
print("Missing SVG files:")
print(missing_svg_files)
print("Extra SVG files:")
print(extra_svg_files)
print("Retired SVG files:")
print(retired_svg_files)
# %%
