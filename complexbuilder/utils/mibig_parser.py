# %%
import json


def parse_mibig_json(mibig_json: str) -> dict:
    """Parse a MIBiG JSON file and return a dictionary with the parsed data."""
    with open(mibig_json, "r") as f:
        mibig_data = json.load(f)
    return mibig_data


mibig_data = parse_mibig_json(
    "/Users/YoshitakaM/Downloads/mibig_json_3.1/BGC0000028.json"
)

# %%
finalproduct_name = mibig_data["cluster"]["compounds"][0]["name"]
