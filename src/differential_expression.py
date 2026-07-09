"""
Phase 4
Differential Expression Analysis

Compares:
Alzheimer's disease vs Control samples

Outputs:
- differential_expression_results.csv
- significant_genes.csv
"""

import os
import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


# -----------------------------
# Paths
# -----------------------------

EXPRESSION_FILE = "data/processed/expression_matrix.csv"
METADATA_FILE = "data/metadata/sample_groups.csv"

RESULT_DIR = "results"

os.makedirs(
    RESULT_DIR,
    exist_ok=True
)


# -----------------------------
# Load Data
# -----------------------------

def load_data():

    expression = pd.read_csv(
        EXPRESSION_FILE,
        index_col=0
    )

    metadata = pd.read_csv(
        METADATA_FILE
    )

    return expression, metadata


# -----------------------------
# Differential Expression
# -----------------------------

def differential_expression(expression, metadata):

    print("Running differential expression...")


    # log transform
    log_expression = np.log2(
        expression + 1
    )


    # sample groups

    disease_samples = metadata[
        metadata["Group"] == "Alzheimer's disease"
    ]["Sample"].tolist()


    control_samples = metadata[
        metadata["Group"] == "Control"
    ]["Sample"].tolist()


    results = []


    for gene in log_expression.index:

        disease_values = log_expression.loc[
            gene,
            disease_samples
        ]

        control_values = log_expression.loc[
            gene,
            control_samples
        ]


        # Fold change

        log2FC = (
            disease_values.mean()
            -
            control_values.mean()
        )


        # Statistical test

        statistic, pvalue = ttest_ind(
            disease_values,
            control_values,
            equal_var=False
        )


        results.append(
            [
                gene,
                log2FC,
                pvalue
            ]
        )


    results = pd.DataFrame(
        results,
        columns=[
            "Gene_ID",
            "log2FoldChange",
            "pvalue"
        ]
    )


    # Multiple testing correction

    results["padj"] = multipletests(
        results["pvalue"],
        method="fdr_bh"
    )[1]


    return results


# -----------------------------
# Main
# -----------------------------

def main():

    expression, metadata = load_data()


    print(
        "Expression matrix:",
        expression.shape
    )


    results = differential_expression(
        expression,
        metadata
    )


    output = os.path.join(
        RESULT_DIR,
        "differential_expression_results.csv"
    )


    results.to_csv(
        output,
        index=False
    )


    print(
        "Saved:",
        output
    )


    # Significant genes

    significant = results[
        (results["padj"] < 0.05)
        &
        (abs(results["log2FoldChange"]) > 1)
    ]


    significant_file = os.path.join(
        RESULT_DIR,
        "significant_genes.csv"
    )


    significant.to_csv(
        significant_file,
        index=False
    )


    print(
        "Significant genes:",
        len(significant)
    )


    print(
        "Saved:",
        significant_file
    )


if __name__ == "__main__":
    main()