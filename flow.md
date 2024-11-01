
# Flow (To Do)

1. BGC accessionナンバーを取得し、Biosynthetic class、Organism name, Main product, Genes情報を取得する。
2. gbkファイルから、BGCを構成するタンパク質名とアミノ酸配列を取得する
3. 各アミノ酸の配列のホモ二量体またはヘテロ二量体のColabFold向けインプットFASTAファイルを作成する。
   1. このとき、core biosynthetic genesについてはペアを作成しないか、またはオプションでドメインごとに区切って他のドメインとペア作成することを検討する
4. インプットファイルに対するcsvファイルも作成するとよい？
5. ColabFoldで構造予測を行う
