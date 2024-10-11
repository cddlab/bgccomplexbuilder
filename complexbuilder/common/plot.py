# %%
import json
from pathlib import Path
from string import ascii_lowercase, ascii_uppercase

import matplotlib.pyplot as plt
import numpy as np

alphabet_list = list(ascii_uppercase + ascii_lowercase)


def plot_ticks(Ls, axes=None):
    if axes is None:
        axes = plt.gca()
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color="black")
        plt.plot([L, L], [0, Ln], color="black")
    ticks = np.cumsum([0] + Ls)
    ticks = (ticks[1:] + ticks[:-1]) / 2
    axes.set_yticks(ticks)
    axes.set_yticklabels(alphabet_list[: len(ticks)])


def plot_paes(paes, Ls=None, dpi=100, fig=True):
    num_models = len(paes)
    if fig:
        plt.figure(figsize=(3 * num_models, 2), dpi=dpi)
    for n, pae in enumerate(paes):
        axes = plt.subplot(1, num_models, n + 1)
        plot_pae(pae, axes, caption=f"rank_{n + 1}", Ls=Ls)
    return plt


def plot_pae(pae, axes, caption="PAE", caption_pad=None, Ls=None, colorkey_size=1.0):
    axes.set_title(caption, pad=caption_pad)
    Ln = pae.shape[0]
    image = axes.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
    if Ls is not None and len(Ls) > 1:
        plot_ticks(Ls, axes=axes)
    plt.colorbar(mappable=image, ax=axes, shrink=colorkey_size)


def _json_to_pae(file: str | Path) -> np.ndarray:
    """create a predicted aligned error (PAE) map from a json file of ColabFold"""
    file = "/Users/YoshitakaM/Desktop/Orf1196_Orf1198T/Orf1196_Orf1198_T_predicted_aligned_error_v1.json"  # noqa
    with open(file, "r") as f:
        data = json.load(f)
    pae = data["predicted_aligned_error"]

    pae_np = np.asarray(pae)
    return pae_np


# paes = [
#     _json_to_pae(
#         "/Users/YoshitakaM/Desktop/Orf1196_Orf1198T/
# Orf1196_Orf1198_T_predicted_aligned_error_v1.json"
#     )
# ]
# paes_plot = plot_paes(paes, Ls=[599, 130], dpi=300)
# paes_plot.savefig(str("foo.png"), bbox_inches="tight")
# paes_plot.close()
# %%
