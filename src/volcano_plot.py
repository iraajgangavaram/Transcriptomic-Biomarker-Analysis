"""
Phase 4.1
Volcano Plot for Differential Expression
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


RESULT_FILE = "results/differential_expression_results.csv"

FIGURE_DIR = "figures"

os.makedirs(
    FIGURE_DIR,
    exist_ok=True
)


def main():

    results = pd.read_csv(
        RESULT_FILE
    )

    results["-log10_pvalue"] = (
        -np.log10(results["pvalue"])
    )


    results["Significance"] = "Not Significant"


    results.loc[
        (results["padj"] < 0.05)
        &
        (results["log2FoldChange"] > 1),
        "Significance"
    ] = "Upregulated"


    results.loc[
        (results["padj"] < 0.05)
        &
        (results["log2FoldChange"] < -1),
        "Significance"
    ] = "Downregulated"



    plt.figure(figsize=(10,6))


    for group in results["Significance"].unique():

        subset = results[
            results["Significance"] == group
        ]

        plt.scatter(
            subset["log2FoldChange"],
            subset["-log10_pvalue"],
            s=10,
            label=group
        )


    plt.axvline(
        1,
        linestyle="--"
    )

    plt.axvline(
        -1,
        linestyle="--"
    )


    plt.axhline(
        -np.log10(0.05),
        linestyle="--"
    )


    plt.xlabel(
        "log2 Fold Change"
    )

    plt.ylabel(
        "-log10(p-value)"
    )


    plt.title(
        "Volcano Plot: Alzheimer's Disease vs Control"
    )


    plt.legend()


    plt.tight_layout()


    output = os.path.join(
        FIGURE_DIR,
        "volcano_plot.png"
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