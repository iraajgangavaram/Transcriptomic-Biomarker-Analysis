"""
Phase 4.2
Heatmap of Top Differentially Expressed Genes
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


EXPRESSION_FILE = "data/processed/expression_matrix.csv"
RESULT_FILE = "results/differential_expression_results.csv"

FIGURE_DIR = "figures"

os.makedirs(
    FIGURE_DIR,
    exist_ok=True
)


def main():

    # Load data

    expression = pd.read_csv(
        EXPRESSION_FILE,
        index_col=0
    )

    results = pd.read_csv(
        RESULT_FILE
    )


    # Select top 50 significant genes

    top_genes = (
        results
        .sort_values(
            "padj"
        )
        .head(50)["Gene_ID"]
        .tolist()
    )


    # Extract expression

    heatmap_data = expression.loc[
        top_genes
    ]


    # Log transform

    heatmap_data = np.log2(
        heatmap_data + 1
    )


    # Z-score by gene

    heatmap_data = (
        heatmap_data
        .sub(heatmap_data.mean(axis=1), axis=0)
        .div(heatmap_data.std(axis=1), axis=0)
    )


    # Plot

    plt.figure(figsize=(10,12))

    plt.imshow(
        heatmap_data,
        aspect="auto"
    )

    plt.colorbar(
        label="Z-score"
    )


    plt.xticks(
        range(len(heatmap_data.columns)),
        heatmap_data.columns,
        rotation=45
    )


    plt.yticks(
        range(len(heatmap_data.index)),
        heatmap_data.index,
        fontsize=6
    )


    plt.title(
        "Top Differentially Expressed Genes"
    )


    plt.tight_layout()


    output = os.path.join(
        FIGURE_DIR,
        "top_genes_heatmap.png"
    )


    plt.savefig(
        output,
        dpi=300
    )


    plt.close()


    print(
        "✓ Saved:",
        output
    )


if __name__ == "__main__":
    main()