#!/usr/bin/env python3
# %%
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
from rdkit.Chem.Draw import IPythonConsole

PandasTools.RenderImagesInAllDataFrames(True)


# %%
mibigjsondirectory = Path("/Users/YoshitakaM/Downloads/mibig_json_4.0")
mibig_json_file = mibigjsondirectory / "BGC0000001.json"

with mibig_json_file.open("r") as f:
    mibig_data = json.load(f)

products = [d["structure"] for d in mibig_data["compounds"]]

smiles = products[0]


# %%
def write_html(df, output):
    scripts = """
    <link href="https://cdnjs.cloudflare.com/ajax/libs/foundation/6.4.3/css/foundation.min.css" rel="stylesheet"/>
    <link href="https://cdn.datatables.net/v/zf/jq-3.6.0/dt-1.13.4/b-2.3.6/b-html5-2.3.6/date-1.4.1/fh-3.3.2/sb-1.4.2/datatables.min.css" rel="stylesheet"/>

    <script src="https://cdn.datatables.net/v/zf/jq-3.6.0/dt-1.13.4/b-2.3.6/b-html5-2.3.6/date-1.4.1/fh-3.3.2/sb-1.4.2/datatables.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/foundation/6.4.3/js/foundation.min.js"></script>

    <script>
        $(document).ready(function() {$('.my-table').DataTable({
            select: true,
            displayLength: 25,
            buttons: ['copy'],
            fixedHeader: true,
            dom: 'iQrtBlp',
        });})
    </script>
    """

    html = df.to_html(classes="my-table")
    html = scripts + html
    with open(output, mode="w") as f:
        f.write(html)


df = pd.DataFrame({"SMILES": ["CCO", "c1ccccc1", "CC(=O)O"]})

# SMILES列からMolオブジェクトを生成して "ROMol" という列に追加
PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES", molCol="ROMol")
df.head()
write_html(df, "hoge.html")
# %%
