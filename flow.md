
# Flow (To Do)

1. BGC accessionナンバーを取得し、Biosynthetic class、Organism name, Main product, Genes情報を取得する。
2. gbkファイルから、BGCを構成するタンパク質名とアミノ酸配列を取得する
3. 各アミノ酸の配列のホモ二量体またはヘテロ二量体のColabFold向けインプットFASTAファイルを作成する。
   1. このとき、core biosynthetic genesについてはペアを作成しないか、またはオプションでドメインごとに区切って他のドメインとペア作成することを検討する
4. インプットファイルに対するcsvファイルも作成するとよい？
5. ColabFoldで構造予測を行う

NRPS_PKSについて、

1. NRPS_PKSのアミノ酸配列を取得する
2. NRPS_PKSのprotein_idのうち、

    There are seven "biosyn_class" values (Secondary Metabolite record type)
    in MIBiG JSON files:
      1. Polyketide  e.g. BGC0000001, BGC0000028
      2. NRP         e.g. BGC0000034, BGC0000081
      3. RiPP        e.g. BGC0000468, BGC0000472, BGC0000583
      4. Terpene     e.g. BGC0000104, BGC0001078
      5. Saccharide  e.g. BGC0000036, BGC0000055
      6. Alkaloid    e.g. BGC0000017, BGC0000815
      7. Other       e.g. BGC0000283, BGC0000314
    `--only-complete` will pair protein from the BGCs marked as "complete"
    in the MIBiG JSON file. (default: False)

## やるべきこと

計算し終わった複合体の結果は、Protein Data Bankに登録されているタンパク質の複合体情報とどれほど一致しているか？
-
