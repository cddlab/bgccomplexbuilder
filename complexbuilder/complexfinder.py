#!/usr/bin/env python3
# %%
import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from loguru import logger
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.morphology import remove_small_objects, skeletonize


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


def show_valley_mask(subdirectory: str) -> np.ndarray:
    """
    Visualize the valley mask.
    Args:
        subdirectory (str): Subdirectory name.
    Returns:
        labeled_valleys (np.ndarray): Labeled valley mask.
    """
    print(subdirectory)
    confidencefile = pae_directory / subdirectory / f"{subdirectory}_confidences.json"
    with open(confidencefile) as f:
        confidencedata = json.load(f)
    pae = np.array(confidencedata["pae"])
    pae_rev = 32.75 - pae

    data_smooth = ndimage.gaussian_filter(pae, sigma=1)
    threshold = threshold_otsu(data_smooth)
    valley_mask = data_smooth < threshold
    valley_mask = remove_small_objects(valley_mask, min_size=20)
    labeled_valleys = label(valley_mask, return_num=False)
    plt.figure(figsize=(6, 6))
    plt.imshow(pae, cmap="Greens_r")
    plt.colorbar(label="Value")
    plt.contour(labeled_valleys, colors="blue", linewidths=1)
    plt.title("Detected Valley Regions")
    plt.show()

    return labeled_valleys


# %%
pae_directory = Path("/Users/YoshitakaM/Desktop/BGC0000028/af3")
a = show_valley_mask("adc79625.1_adc79625.1")
# %%
b = show_valley_mask("adc79628.1_adc79628.1")
# %%
