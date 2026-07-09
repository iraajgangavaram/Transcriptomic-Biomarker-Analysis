"""
Phase 6
STRING Protein-Protein Interaction Network Analysis
"""

import os
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


INPUT_FILE = "results/differential_expression_gene_symbols.csv"

RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(
    RESULT_DIR,
    exist_ok=True
)

os.makedirs(
    FIGURE_DIR,
    exist_ok=True
)


STRING_URL = (
    "https://string-db.org/api/tsv/network"
)


def get_string_network(genes):

    params = {
        "identifiers": "%0d".join(genes),
        "species": 9606,
        "required_score": 200
    }

    response = requests.get(
        STRING_URL,
        params=params
    )

    response.raise_for_status()

    return response.text


def main():

    data = pd.read_csv(
        INPUT_FILE
    )


    # Select top genes

    genes = (
        data
        .sort_values("padj")
        .head(200)["Gene_Symbol"]
        .tolist()
    )


    print(
        "Genes submitted:",
        len(genes)
    )


    print(
        "Querying STRING..."
    )


    network_text = get_string_network(
        genes
    )


    output = os.path.join(
        RESULT_DIR,
        "string_network.csv"
    )


    with open(output, "w") as file:
        file.write(network_text)


    print(
        "Saved:",
        output
    )


    # Read network

    network = pd.read_csv(
        output,
        sep="\t"
    )


    # Build graph

    G = nx.from_pandas_edgelist(
        network,
        "preferredName_A",
        "preferredName_B",
        edge_attr="score"
    )


    print(
        "Nodes:",
        G.number_of_nodes()
    )

    print(
        "Edges:",
        G.number_of_edges()
    )


    # Plot

    plt.figure(
        figsize=(10,8)
    )


    pos = nx.spring_layout(
        G,
        seed=42
    )


    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=500,
        font_size=8
    )


    plt.title(
        "STRING Protein Interaction Network"
    )


    plt.tight_layout()


    figure = os.path.join(
        FIGURE_DIR,
        "interaction_network.png"
    )


    plt.savefig(
        figure,
        dpi=300
    )


    plt.close()


    print(
        "Saved:",
        figure
    )


if __name__ == "__main__":
    main()