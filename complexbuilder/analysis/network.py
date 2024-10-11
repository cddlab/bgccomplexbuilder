#!/usr/bin/env python3

import argparse
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from Bio import SeqIO
from loguru import logger

from complexbuilder.analysis.extract_hitcomplexes import split_proteinids
from complexbuilder.common.log import log_setup
from complexbuilder.common.parser import sanitised_name

log_setup(level="DEBUG")
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def get_bgcgenes(dataInt: dict[str, dict], bgc_id: str) -> list[str]:
    """
    Extract a unique list of gene names for a given BGC from the input dataInt.

    For the specified BGC key in dataInt, each key (e.g. protein header)
    is split using split_proteinids() into a tuple of gene names.
    The function then returns a list of unique gene names (bgcgenes).

    Args:
        dataInt (dict[str, dict]):
            A dictionary containing BGC data where each key is a BGC identifier.
        bgc_id (str):
            The BGC identifier for which the gene list is to be extracted.

    Returns:
        list[str]: A list of unique gene names for the specified BGC.
    """
    bgcgenes = []
    for header in dataInt[bgc_id]:
        geneNames = split_proteinids(header)
        for name in geneNames:
            if name not in bgcgenes:
                bgcgenes.append(name)
    return bgcgenes


def extract_cds_info(mibiggbkdir: str, bgc_id: str):
    """
    Retrieve gene information for a given BGC ID.
    This function reads a GenBank file and extracts information about
    each CDS feature.

    Args:
        gbk_file (str):
    """
    cds_info_list = []
    gbk_file = os.path.join(mibiggbkdir, f"{bgc_id}.gbk")
    with open(gbk_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers
                    # each CDS feature may have multiple qualifiers
                    # protein_id, gene, locus_tag, product
                    # are common qualifiers for CDS features
                    # but may not be present in all records
                    # so we use .get() to avoid KeyError
                    protein_id = qualifiers.get("protein_id", [None])[0]
                    gene = qualifiers.get("gene", [None])[0]
                    locus_tag = qualifiers.get("locus_tag", [None])[0]
                    product = qualifiers.get("product", [None])[0]
                    cds_info = {
                        "protein_id": protein_id,
                        "gene": gene,
                        "locus_tag": locus_tag,
                        "product": product,
                    }
                    cds_info_list.append(cds_info)
    return cds_info_list


def make_node_desc_dict(
    bgcgenes: list[str],
    cds_info_list: list[dict[str, str]],
) -> dict[str, str]:
    """
    Return a dictionary mapping each node in bgcgenes to a concatenated string of available labels.

    For each node (e.g. a protein identifier, gene name, or locus_tag) in bgcgenes,
    the function searches cds_info_list for a CDS feature where one of the keys
    "protein_id", "gene", or "locus_tag" matches the node (after applying sanitised_name()).
    When a match is found, the description is built by concatenating (with newline separators)
    the values of protein_id, gene, locus_tag and product from that CDS feature.

    Example:
        bgcgenes = ["aek75497.1", "aek75507.1"]
        cds_info_list = [
            {
                "protein_id": "AEK75497.1",
                "gene": "gene1",
                "locus_tag": "locus1",
                "product": "product1",
            },
            {
                "protein_id": "AEK75507.1",
                "gene": None,
                "locus_tag": None,
                "product": "product2",
            },
            {
                "protein_id": None,
                "gene": "abyU",
                "locus_tag": None,
                "product": "product3",
            },
        ]

        -> {
             "aek75497.1": "AEK75497.1\ngene1\nlocus1\nproduct1",
             "aek75507.1": "AEK75507.1\nproduct2",
             "abyU": "abyU\nproduct3",
           }

    Args:
        bgcgenes: List of gene identifiers (nodes) from the BGC.
        cds_info_list: List of dictionaries, each containing keys "protein_id", "gene",
                       "locus_tag" and "product" for a CDS feature.

    Returns:
        A dictionary mapping each node to its description string.
    """
    node_desc_dict = {}
    for node in bgcgenes:
        for cds_info in cds_info_list:
            for label in ["protein_id", "gene", "locus_tag"]:
                value = cds_info.get(label)
                if value is None:
                    continue
                if sanitised_name(value) == sanitised_name(node):
                    parts = []
                    if cds_info.get("protein_id"):
                        parts.append(cds_info["protein_id"])
                    if cds_info.get("gene"):
                        parts.append(cds_info["gene"])
                    if cds_info.get("locus_tag"):
                        parts.append(cds_info["locus_tag"])
                    if cds_info.get("product"):
                        parts.append(cds_info["product"])
                    node_desc_dict[node] = "\n".join(parts)
                    break
            if node in node_desc_dict:
                break
    return node_desc_dict


def fold_description(description: str, width: int = 20) -> str:
    """
    Wrap a text string so that each line does not exceed width characters.
    The function attempts to break lines at spaces rather than cutting off
    exactly at 'width'.

    Args:
        description (str): The description to fold.
        width (int): The maximum line width.

    Returns:
        str: The folded description.
    """
    original_lines = description.splitlines()
    wrapped_segments = []
    for segment in original_lines:
        # Wrap each segment separately.
        words = segment.split()
        lines = []
        current_line = ""
        for word in words:
            if current_line:
                if len(current_line) + 1 + len(word) > width:
                    lines.append(current_line)
                    current_line = word
                else:
                    current_line += " " + word
            else:
                current_line = word
        if current_line:
            lines.append(current_line)
        # If the segment was empty, preserve the empty line.
        wrapped_segments.append("\n".join(lines) if lines else "")
    # Join the wrapped segments with newline characters.
    return "\n".join(wrapped_segments)


def make_network_svg(
    mibiggbkdir: str,
    hitcomplexesPath: str,
    rmsdPath: str,
    outputdir: str,
    rmsd_threshold: float = 2.0,
    max_count: int = 3000,
) -> None:
    """
    Create a network SVG for each BGC in the hitcomplexes data.

    Args:
        mibiggbkdir (str): Path to the directory containing MIBiG GenBank files.
        hitcomplexesPath (str): Path to the JSON file containing hit complexes.
        rmsdPath (str): Path to the JSON file containing RMSD data.
        outputdir (str): Directory where the SVG files will be saved.
        max_count (int): Maximum number of BGCs to process. Default is 3000.
    """
    if not os.path.exists(hitcomplexesPath):
        logger.error(f"Hit complexes file not found: {hitcomplexesPath}")
        return

    count = 0

    with open(hitcomplexesPath, "r") as f:
        dataInt = json.load(f)
    with open(rmsdPath, "r") as f:
        rmsd_data = json.load(f)
    for bgc_id in dataInt:
        bgcgenes = get_bgcgenes(dataInt, bgc_id)
        cds_info_list = extract_cds_info(mibiggbkdir, bgc_id)
        node_desc = make_node_desc_dict(bgcgenes, cds_info_list)
        G = nx.Graph()
        for gene in bgcgenes:
            # add nodes with descriptions
            G.add_node(
                gene,
                description=node_desc.get(gene, "No description available"),
            )
        for protein_dimer in dataInt[bgc_id]:
            # protein dimer is like "aek75497.1_aek75507.1"
            protein_a, protein_b = split_proteinids(protein_dimer)
            iptm_score = dataInt[bgc_id][protein_dimer]["ipTM"]
            rmsd_value = (
                rmsd_data.get(bgc_id, {})
                .get(protein_dimer, {})
                .get("RMSD", float("inf"))
            )
            G.add_edge(protein_a, protein_b, iptm=iptm_score, rmsd=rmsd_value)
        iptms = [G[u][v]["iptm"] for u, v in G.edges()]
        for u, v, data in G.edges(data=True):
            w = data.get("iptm")
            if not isinstance(w, (int, float)):
                print(f"problematic iptm: edge=({u}, {v}), iptm={w}, type={type(w)}")

        ### Draw Network ###
        num_nodes = len(G.nodes)
        # proportional size based on number of nodes.
        base_width, base_height = 8.0, 6.0
        scale_factor_w, scale_factor_h = 0.80, 0.60
        fig_width = base_width + scale_factor_w * num_nodes
        fig_height = base_height + scale_factor_h * num_nodes
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=300)
        pos = nx.circular_layout(G, scale=1)
        # edge_color is set to a colormap based on iptms
        # iptms are normalized to the range [0, 1] for colormap
        cmap = plt.get_cmap("coolwarm")
        rmsd_values = [G[u][v]["rmsd"] for u, v in G.edges()]
        # RMSD < rmsd_thresholdの場合、紫色に設定
        edge_colors = [
            cmap(w) if rmsd > rmsd_threshold else (0, 0, 1.0, 1.0)
            for w, rmsd in zip(iptms, rmsd_values, strict=False)
        ]
        # widths
        widths = [w * 8 for w in iptms]
        # node size is proportional to the length of the description
        node_sizes = [len(G.nodes[node]["description"]) * 120 for node in G.nodes]
        nx.draw(
            G,
            pos,
            ax=ax,
            width=widths,
            edge_color=edge_colors,
            edge_cmap=cmap,
            alpha=0.5,
            node_color="lightblue",
            node_size=node_sizes,
            font_weight="bold",
        )

        edge_labels = nx.get_edge_attributes(G, "iptm")
        node_labels = {
            node: fold_description(G.nodes[node]["description"], width=20)
            for node in G.nodes
        }

        # draw edge labels and node labels
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax)
        nx.draw_networkx_labels(
            G, pos, labels=node_labels, font_size=10, font_family="Arial", ax=ax
        )
        ax.axis("off")
        os.makedirs(outputdir, exist_ok=True)
        fig.savefig(os.path.join(outputdir, f"{bgc_id}.svg"))
        plt.close()
        plt.clf()
        count += 1
        print(f"Processed {count} / {len(dataInt)}: {bgc_id}")
        if count > max_count:
            print(f"More than {max_count} BGCs processed, stopping to avoid overload.")
            break


def main():
    parser = argparse.ArgumentParser(description="Draw network SVGs for hit complexes.")
    parser.add_argument(
        "--mibiggbkdir",
        metavar="MIBiG GenBank directory",
        type=str,
        help="Directory containing MIBiG GenBank files downloaded from MIBiG 4.0."
        " e.g. /path/to/Downloads/mibig_gbk_4.0",
    )
    parser.add_argument(
        "--hitcomplexes_path",
        metavar="The filename of the output JSON file.",
        type=str,
        help="Input JSON file containing hit complexes. This file is generated by extract_hitcomplexes.py."
        " e.g. hitcomplex_iptm0.55_ipsae0.0_all.json",
    )
    parser.add_argument(
        "--rmsdjson_path",
        metavar="The filename of the RMSD JSON file.",
        type=str,
        help="Input JSON file containing RMSD data. This file is generated by structhomologues.py."
        " e.g. rmsd.json",
    )
    parser.add_argument(
        "--outputdir",
        metavar="Output directory",
        type=str,
        help="Directory to save output files.",
    )
    parser.add_argument(
        "--rmsd_threshold",
        metavar="RMSD threshold",
        type=float,
        default=2.0,
        help="RMSD threshold for edge coloring. Edges with RMSD below this value will be colored blue.",
    )
    args = parser.parse_args()
    make_network_svg(
        args.mibiggbkdir,
        args.hitcomplexes_path,
        args.rmsdjson_path,
        args.outputdir,
        args.rmsd_threshold,
    )


if __name__ == "__main__":
    main()
