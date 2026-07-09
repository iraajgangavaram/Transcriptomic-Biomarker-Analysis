"""
Phase 6.1
Hub Gene Identification
"""

import os
import pandas as pd
import networkx as nx


INPUT_FILE = "results/string_network.csv"

OUTPUT_FILE = "results/hub_genes.csv"


def main():

    network = pd.read_csv(
        INPUT_FILE,
        sep="\t"
    )


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


    # Centrality calculations

    degree = nx.degree_centrality(G)

    betweenness = nx.betweenness_centrality(G)


    hub_data = pd.DataFrame(
        {
            "Gene": list(degree.keys()),
            "Degree_Centrality": list(degree.values()),
            "Betweenness_Centrality": [
                betweenness[g]
                for g in degree.keys()
            ]
        }
    )


    hub_data = hub_data.sort_values(
        by="Degree_Centrality",
        ascending=False
    )


    hub_data.to_csv(
        OUTPUT_FILE,
        index=False
    )


    print(
        "Saved:",
        OUTPUT_FILE
    )


    print(
        "\nTop Hub Genes:"
    )

    print(
        hub_data.head(10)
    )


if __name__ == "__main__":
    main()