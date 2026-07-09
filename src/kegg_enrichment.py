"""
Phase 5.3
KEGG Pathway Enrichment Analysis
"""

import os
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import math


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


def main():

    data = pd.read_csv(
        INPUT_FILE
    )


    # Use top ranked genes because of small sample size

    top_genes = (
        data
        .sort_values("padj")
        .head(200)
    )


    gene_list = top_genes[
        "Gene_Symbol"
    ].tolist()


    print(
        "Genes used:",
        len(gene_list)
    )


    print(
        "Running KEGG enrichment..."
    )


    enrichment = gp.enrichr(
        gene_list=gene_list,
        gene_sets=[
            "KEGG_2021_Human"
        ],
        organism="human",
        outdir=None
    )


    results = enrichment.results


    output = os.path.join(
        RESULT_DIR,
        "kegg_enrichment.csv"
    )


    results.to_csv(
        output,
        index=False
    )


    print(
        "Saved:",
        output
    )


    # Plot top 10 pathways

    top = results.head(10)


    plt.figure(figsize=(10,6))


    plt.barh(
        top["Term"],
        -top["Adjusted P-value"].apply(
            math.log10
        )
    )


    plt.xlabel(
        "-log10 Adjusted P-value"
    )


    plt.title(
        "Top KEGG Pathways"
    )


    plt.tight_layout()


    figure = os.path.join(
        FIGURE_DIR,
        "kegg_enrichment_barplot.png"
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