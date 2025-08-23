#!/usr/bin/env python3
import json
from pathlib import Path

import pandas as pd
from loguru import logger
from rdkit.Chem import PandasTools
from rdkit.Chem.PandasTools import ChangeMoleculeRendering

from complexbuilder.common.log import log_setup

log_setup(level="SUCCESS")


def write_html(df, output):
    scripts = """
    <link href="https://cdnjs.cloudflare.com/ajax/libs/foundation/6.9.0/css/foundation.min.css" rel="stylesheet"/>
    <link href="https://cdn.datatables.net/v/zf/jq-3.6.0/dt-1.13.4/b-2.3.6/b-html5-2.3.6/date-1.4.1/fh-3.3.2/sb-1.4.2/datatables.min.css" rel="stylesheet"/>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/foundation/6.9.0/js/foundation.min.js"></script>
    <script src="https://cdn.datatables.net/v/zf/jq-3.6.0/dt-1.13.4/b-2.3.6/b-html5-2.3.6/date-1.4.1/fh-3.3.2/sb-1.4.2/datatables.min.js"></script>


    <script>
        $(document).ready(function() {$('.my-table').DataTable({
            select: true,
            displayLength: 25,
            buttons: ['copy'],
            fixedHeader: true,
            dom: 'iQrtBlp',
            language: {
                searchBuilder: {
                    title: 'BGC PPI network visualizer'
                }
            },
        });})
    </script>
    """
    style_block = """
    <style type="text/css">
      .test-container {
        text-align: center;
        padding: 20px;
      }
    </style>
    """
    div_start = '<div class="test-container">'
    div_end = "</div>"

    html = df.to_html(classes="my-table", escape=False, index=False)
    caption = """
    (2025.08.21) Updated. Structurally homologous (RMSD &lt; 2.0 Å) protein pairs are highlighted in blue.<br>
    The raw JSON file can be downloaded from <a href="hitcomplex_iptm0.55_ipsae0.0_all.json">here (4.2 MB)</a>. <br>You can download these contents from <a href="networks.zip">here (54MB)</a> and view them offline.<br>Please be patient while the contents are being loaded...
    """
    html = scripts + style_block + div_start + caption + html + div_end
    with open(output, mode="w") as f:
        f.write(html)


def add_data(
    mibig_json_file: Path, svgdirectory: Path, pre_df: pd.DataFrame | None = None
) -> pd.DataFrame:
    """
    Adds data from a MIBiG JSON file to a DataFrame.
    """
    with mibig_json_file.open("r") as f:
        data = json.load(f)

    compounds = [d["structure"] for d in data["compounds"] if "structure" in d]
    if not compounds:
        compounds = [None]
    compoundnames = [d["name"] for d in data["compounds"]]
    accession_id = data["accession"]
    version = data["version"]
    classes = [d["class"] for d in data["biosynthesis"]["classes"] if "class" in d]
    acc_ver = f"{accession_id}.{version}"
    taxname = data["taxonomy"]["name"]
    networksvgpath = Path(svgdirectory, f"{accession_id}.svg")
    network = f'<a href="{networksvgpath}" target="_blank" rel="noopener noreferrer"><img src="{networksvgpath}" alt="Network" /></a>'

    compoundname = "<br>".join(compoundnames if compoundnames else "No compound name")
    classes_str = "<br>".join(classes if classes else "No class")
    logger.info(f"Processing {mibig_json_file}: {acc_ver} {taxname} {compoundname}")
    df = pd.DataFrame(
        {
            "Accession": [acc_ver],
            "Taxonomy": [taxname],
            "Class": [classes_str],
            "SMILES": [compounds[0]] if compounds else None,
            "Compound Names": [compoundname],
            "Network": [network],
        }
    )
    # Italic font
    df["Taxonomy"] = df["Taxonomy"].apply(lambda x: f"<i>{x}</i>")
    # Hyperlink to MIBiG entry
    df["Accession"] = df["Accession"].apply(
        lambda x: f'<a href="https://mibig.secondarymetabolites.org/repository/{x}" target="_blank" rel="noopener noreferrer">{x}</a>'
    )

    if df["SMILES"].notnull().any():
        PandasTools.AddMoleculeColumnToFrame(
            df, smilesCol="SMILES", molCol="Representative Structure"
        )
    else:
        df["Representative Structure"] = None
    if pre_df is not None:
        df = pd.concat([pre_df, df])
    return df


mibigjsondirectory = Path("/Users/YoshitakaM/Downloads/mibig_json_4.0")
svgdirectory = Path("svg")

for i in range(1, 2827):
    mibig_json_file = mibigjsondirectory / f"BGC000{i:04d}.json"
    if not mibig_json_file.exists():
        continue
    with open(mibig_json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    if data.get("status") != "active":
        continue
    if i == 1:
        print(f"Processing {mibig_json_file}")
        df = add_data(mibig_json_file, svgdirectory)
    else:
        if mibig_json_file.exists():
            df = add_data(mibig_json_file, svgdirectory, df)
ChangeMoleculeRendering(df)
df.drop(columns=["SMILES"], inplace=True)
write_html(df, "publish.html")
