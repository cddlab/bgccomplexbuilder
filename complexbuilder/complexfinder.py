#!/usr/bin/env python3
# %%
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from loguru import logger
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.morphology import remove_small_objects


def get_chain_ids_and_lengths(token_chain_ids: list[str]) -> dict:
    """Get the chain IDs and lengths from the token chain IDs.
    Args:
        token_chain_ids (list[str]): List of token chain IDs
    Returns:
        count_dict (dict): Dictionary containing the chain IDs and lengths
    """
    unique_chain_ids = list(dict.fromkeys(token_chain_ids))
    count_dict = {
        chain_id: token_chain_ids.count(chain_id) for chain_id in unique_chain_ids
    }
    return count_dict


def create_valley_mask(
    confidencefile: str | Path, threshold: float = -1.0
) -> np.ndarray:
    """
    Create a valley mask from the PAE matrix in a confidence file.
    Args:
        confidencefile (str | Path): Path to the confidence file.
    Returns:
        labeled_valleys (np.ndarray): Labeled valley mask.
    """

    with open(confidencefile) as f:
        confidencedata = json.load(f)
    pae = np.array(confidencedata["pae"])

    data_smooth = ndimage.gaussian_filter(pae, sigma=1)
    threshold = threshold_otsu(data_smooth) if threshold < 0 else threshold
    logger.debug(f"Threshold: {threshold}")
    valley_mask = remove_small_objects(data_smooth < threshold, min_size=100)
    return np.asarray(valley_mask, dtype=int)


def list_subdirectories(path: str | Path) -> list[str]:
    return [p.name for p in Path(path).iterdir() if p.is_dir()]


def map_with_colorbar(
    fig,
    ax,
    data,
    labeled_valleys,
    chain_ids_and_lengths,
    model_name="",
    cmap="Greens_r",
    vmin=0,
    vmax=31.75,
    **kwargs,
):
    """Add a colorbar to the plot.

    Args:
        mappable (plt.cm.ScalarMappable): ScalarMappable object
        ax: Axes object
    """
    ax.set_title(model_name)
    ax.set_xlabel("Scored Residue")
    ax.set_ylabel("Aligned Residue")
    mappable: plt.cm.ScalarMappable = ax.imshow(
        data["pae"],
        label=model_name,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    # add black trace between the boundaries of the chains
    pos = 0
    positions = []
    labels = []
    for chain_id, chain_len in chain_ids_and_lengths.items():
        pos += chain_len
        ax.axvline(x=pos, color="black", linewidth=0.5)
        ax.axhline(y=pos, color="black", linewidth=0.5)
        labels.append(chain_id)
        positions.append(pos)
    ax.set_xlim(0, pos)
    ax.set_ylim(pos, 0)
    ax.set_aspect("equal")
    ax.contour(labeled_valleys, colors="blue", linewidths=1)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.5)
    cbar = fig.colorbar(mappable, cax=cax, ax=ax, **kwargs)
    return cbar


def create_interchain_mask(chain_ids_and_lengths: dict[str, int]) -> np.ndarray:
    """
    Create a mask for the interchain regions.
    Args:
        chain_ids_and_lengths (dict[str, int]): Dictionary containing the chain IDs
        and lengths
    Returns:
        interchain_mask (np.ndarray[int, int]): Mask for the interchain regions.
        Regions from the same chain are set to 0, otherwise 1.
    """
    total_length = sum(chain_ids_and_lengths.values())
    interchain_mask = np.ones((total_length, total_length), dtype=int)
    pos = 0
    # set the interchain regions to 0
    for _, chain_len in chain_ids_and_lengths.items():
        interchain_mask[pos : pos + chain_len, pos : pos + chain_len] = 0
        pos += chain_len
    return interchain_mask


# %%

logger.remove()
af3directory = Path("/Users/YoshitakaM/Desktop/BGC0001296/af3")
subdirectories = list_subdirectories(af3directory)
cmap = "Greens_r"
subdirectory = "bat51058.1_bat51062.1"
confidencefile = af3directory / subdirectory / f"{subdirectory}_confidences.json"
with open(confidencefile) as f:
    data = json.load(f)
chain_ids_and_lengths = get_chain_ids_and_lengths(data["token_chain_ids"])
lowpae_mask = create_valley_mask(confidencefile)
interchain_mask = create_interchain_mask(chain_ids_and_lengths)
interchain_valley_mask = lowpae_mask * interchain_mask
aligned_residues = np.sum(interchain_valley_mask, axis=0)
scored_residues = np.sum(interchain_valley_mask, axis=1)

# %%
for subdirectory in subdirectories:
    confidencefile = af3directory / subdirectory / f"{subdirectory}_confidences.json"
    with open(confidencefile) as f:
        data = json.load(f)
    chain_ids_and_lengths = get_chain_ids_and_lengths(data["token_chain_ids"])
    labeled_valleys = create_valley_mask(confidencefile)
    fig, ax = plt.subplots(figsize=(3.6, 4.2), dpi=100)
    map_with_colorbar(
        fig,
        ax,
        data,
        labeled_valleys,
        chain_ids_and_lengths,
        model_name=subdirectory,
        cmap=cmap,
        orientation="horizontal",
        pad=0.2,
        label="Expected Position Error (Ångströms)",
    )
    plt.tight_layout()
    plt.savefig(f"{af3directory}/{subdirectory}/{subdirectory}_valleys.png")
    # print iptm value
    summaryfile = (
        af3directory / subdirectory / f"{subdirectory}_summary_confidences.json"
    )
    # with open(summaryfile) as f:
    #     summarydata = json.load(f)
    # print(f"dir: {subdirectory}, iptm: {summarydata['iptm']}")
# %%
