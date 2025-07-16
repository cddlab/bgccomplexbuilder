#!/usr/bin/env python3
# %%
import json
from pathlib import Path

import pandas as pd
from loguru import logger
from rdkit.Chem import PandasPatcher, PandasTools
from rdkit.Chem.PandasTools import ChangeMoleculeRendering

from complexbuilder.common.log import log_setup

log_setup(level="DEBUG")


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

    html = df.to_html(classes="my-table", escape=False)
    html = scripts + html
    with open(output, mode="w") as f:
        f.write(html)


# %%
def add_data(mibig_json_file: Path, pre_df: pd.DataFrame | None = None) -> pd.DataFrame:
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
    acc_ver = f"{accession_id}.{version}"
    taxname = data["taxonomy"]["name"]

    rep_compoundname = compoundnames[0] if compoundnames else "Unknown"
    logger.info(f"Processing {mibig_json_file}: {acc_ver} {taxname} {rep_compoundname}")
    df = pd.DataFrame(
        {
            "Accession": [acc_ver],
            "Taxonomy": [taxname],
            "SMILES": [compounds[0]] if compounds else None,
            "Representative Compound Name": [rep_compoundname],
        }
    )
    # Hyperlink to MIBiG entry
    df["Accession"] = df["Accession"].apply(
        lambda x: f'<a href="https://mibig.secondarymetabolites.org/repository/{x}">{x}</a>'
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

for i in range(1, 201):
    mibig_json_file = mibigjsondirectory / f"BGC000{i:04d}.json"
    if i == 1:
        print(f"Processing {mibig_json_file}")
        df = add_data(mibig_json_file)
    else:
        if mibig_json_file.exists():
            df = add_data(mibig_json_file, df)
ChangeMoleculeRendering(df)
df.drop(columns=["SMILES"], inplace=True)
write_html(df, "hoge.html")
# %%
