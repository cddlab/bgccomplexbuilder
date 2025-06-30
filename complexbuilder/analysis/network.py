#!/usr/bin/env python3
# %%
import json
import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
from loguru import logger

from complexbuilder.common.log import log_setup
from complexbuilder.common.parser import split_proteinids

log_setup(level="DEBUG")


# %%
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
        logger.debug(f"geneNames from header '{header}': {geneNames}")
        for name in geneNames:
            if name not in bgcgenes:
                bgcgenes.append(name)
    return bgcgenes


mibigJsonPath = "/Users/YoshitakaM/Downloads/mibig_json_4.0"
hitcomplexesPath = "/Users/YoshitakaM/Desktop/hitcomplexes.json"

# %%
cmap = cm.get_cmap("coolwarm")

noGeneInfoList = []

with open(hitcomplexesPath, "r") as f:
    dataInt = json.load(f)


get_bgcgenes(dataInt, "BGC0000001")  # Example call to the function
# %%

count = 0
for bgc_id in dataInt:
    logger.debug(f"Processing {bgc_id}")
    with open(os.path.join(mibigJsonPath, bgc_id + ".json"), "r") as f:
        dataBGC = json.load(f)
    G = nx.Graph()
    bgcgenes = get_bgcgenes(dataInt, bgc_id)
    for nodeC in bgcgenes:
        nodeC_desc = "no info"
        if "genes" not in dataBGC:
            noGeneInfoList.append(bgc_id)
            print(f"No gene info for {bgc_id}, skipping node {nodeC}")
            continue
        for geneInfo in dataBGC["genes"]:
            for geneInfo2 in dataBGC["genes"][geneInfo]:
                if geneInfo2["id"].lower() == nodeC:
                    if "name" in geneInfo2 and geneInfo2["name"]:
                        if nodeC_desc == "no info":
                            nodeC_desc = ""
                        nodeC_desc += geneInfo2["name"] + "\n"
                    if "product" in geneInfo2 and geneInfo2["product"]:
                        if nodeC_desc == "no info":
                            nodeC_desc = ""
                        nodeC_desc += geneInfo2["product"] + "\n"
        G.add_node(nodeC, description=nodeC_desc)
    for j in dataInt[bgc_id]:
        geneName_list = split_proteinids(j)
        complesValue = [dataInt[bgc_id][j]["ipTM"]]
        G.add_edge(geneName_list[0], geneName_list[1], weight=complesValue[0])
    # ノード位置
    print("getting dataInt is finished")
    # 重み取得
    weights = [G[u][v]["weight"] for u, v in G.edges()]
    for u, v, data in G.edges(data=True):
        w = data.get("weight")
        if not isinstance(w, (int, float)):
            print(f"problematic weight: edge=({u}, {v}), weight={w}, type={type(w)}")
    pos = nx.shell_layout(G)
    # 描画
    edge_colors = [cmap(w) for w in weights]
    widths = [w * 10 for w in weights]
    plt.figure(figsize=(12, 8))
    nx.draw(
        G,
        pos,
        # with_labels=True,
        width=widths,
        # width=5,
        edge_color=edge_colors,
        edge_cmap=cmap,
        node_color="lightblue",
        node_size=1000,
        font_weight="bold",
    )

    # エッジラベル（重み表示）
    edge_labels = nx.get_edge_attributes(G, "weight")

    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    # labelの追加
    custom_labels = {
        node: f"{node}\n{G.nodes[node]['description']}" for node in G.nodes
    }
    nx.draw_networkx_labels(G, pos, labels=custom_labels, font_size=10)
    plt.axis("off")
    plt.tight_layout()
    os.makedirs("/Users/YoshitakaM/Desktop/svg2/", exist_ok=True)
    plt.savefig("/Users/YoshitakaM/Desktop/svg2/" + bgc_id + ".svg")  # ← ここで保存
    plt.close()  # 表示せず終了（表示したい場合は plt.show() を使ってもOK）
    plt.clf()
    count += 1
    print(f"Processed {count} / {len(dataInt)}: {bgc_id}")
    if count > 30:
        print("30個以上のBGCを処理しました。")
        break

# In[ ]:
