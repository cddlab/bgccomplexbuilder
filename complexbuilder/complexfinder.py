#!/usr/bin/env python3
# %%
import argparse
import json
import os
from pathlib import Path

import numpy as np
from loguru import logger


def parse_jsonfile(file: str | Path) -> dict:
    """
    Parse a chain_iptm tuple from "summary_confidences.json" file in the
    specified directory.
    Args:
        file (str): Path to "summary_confidences.json" file.
    Returns:
    """
    with open(os.path.join(file)) as f:
        data = json.load(f)
    return data


pae_directory = Path("/Users/YoshitakaM/Desktop/BGC0001296/af3")


subdirectory = "bat51063.1_bat51064.1"
summaryfile = pae_directory / subdirectory / f"{subdirectory}_summary_confidences.json"
summary_data = parse_jsonfile(summaryfile)

print(summary_data["iptm"])

# %%
confidencefile = pae_directory / subdirectory / f"{subdirectory}_confidences.json"
with open(confidencefile) as f:
    confidencedata = json.load(f)
# plot pae as numpy array
pae = np.array(confidencedata["pae"])
# show all keys in confidence data
print(confidencedata.keys())
# %%
