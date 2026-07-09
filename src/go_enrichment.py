"""
Phase 5.2
GO Biological Process Enrichment Analysis
"""

import os
import pandas as pd
import gseapy as gp
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


def main():

    data = pd.read_csv(
        INPUT_FILE
    )


    # Select significant genes

    top_genes = (
    data
    .sort_values("padj")
    .head(200)
)

    significant_genes = top_genes[
    "Gene_Symbol"
    ].tolist()

    print(
        "Significant genes:",
        len(significant_genes)
    )


    if len(significant_genes) == 0:
        print("No significant genes found")
        return


    print(
        "Running GO enrichment..."
    )


    enrichment = gp.enrichr(
        gene_list=significant_genes,
        gene_sets=[
            "GO_Biological_Process_2023"
        ],
        organism="human",
        outdir=None
    )


    results = enrichment.results


    output = os.path.join(
        RESULT_DIR,
        "go_enrichment.csv"
    )


    results.to_csv(
        output,
        index=False
    )


    print(
        "Saved:",
        output
    )


    # Plot top pathways

    top = results.head(10)


    plt.figure(figsize=(10,6))

    plt.barh(
        top["Term"],
        -top["Adjusted P-value"].apply(
            lambda x: __import__("math").log10(x)
        )
    )


    plt.xlabel(
        "-log10 Adjusted P-value"
    )

    plt.title(
        "Top GO Biological Processes"
    )


    plt.tight_layout()


    figure = os.path.join(
        FIGURE_DIR,
        "go_enrichment_barplot.png"
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