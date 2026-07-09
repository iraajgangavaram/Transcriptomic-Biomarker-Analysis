"""
Phase 3
Quality Control for RNA-seq Data

Generates:
1. Raw count distribution
2. Log2(count + 1) distribution
3. Boxplot of log2-transformed counts
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# -----------------------------
# File Paths
# -----------------------------

EXPRESSION_FILE = "data/processed/expression_matrix.csv"
METADATA_FILE = "data/metadata/sample_groups.csv"

FIGURE_DIR = "figures"
RESULTS_DIR = "results"

os.makedirs(FIGURE_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_data():
    expression = pd.read_csv(
        EXPRESSION_FILE,
        index_col=0
    )

    metadata = pd.read_csv(
        METADATA_FILE
    )

    return expression, metadata


def plot_raw_distribution(expression):

    plt.figure(figsize=(10,6))

    plt.hist(
        expression.values.flatten(),
        bins=100
    )

    plt.title("Raw RNA-seq Count Distribution")
    plt.xlabel("Read Count")
    plt.ylabel("Frequency")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGURE_DIR,
            "expression_distribution.png"
        ),
        dpi=300
    )

    plt.close()

    print("✓ expression_distribution.png")


def plot_log_distribution(log_expression):

    plt.figure(figsize=(10,6))

    plt.hist(
        log_expression.values.flatten(),
        bins=100
    )

    plt.title("Log2 RNA-seq Count Distribution")
    plt.xlabel("log2(count + 1)")
    plt.ylabel("Frequency")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGURE_DIR,
            "log_expression_distribution.png"
        ),
        dpi=300
    )

    plt.close()

    print("✓ log_expression_distribution.png")


def plot_boxplot(log_expression):

    plt.figure(figsize=(10,6))

    plt.boxplot(
        log_expression.values,
        tick_labels=log_expression.columns,
        showfliers=False
    )

    plt.xticks(rotation=45)

    plt.title("Sample Expression Boxplot")
    plt.ylabel("log2(count + 1)")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGURE_DIR,
            "boxplot_log_counts.png"
        ),
        dpi=300
    )

    plt.close()

    print("✓ boxplot_log_counts.png")
def plot_pca(log_expression, metadata):
    """
    Perform PCA and save PCA plot.
    """

    # Samples should be rows
    X = log_expression.T

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)

    # Create dataframe
    pca_df = pd.DataFrame(
        pcs,
        columns=["PC1", "PC2"]
    )

    pca_df["Sample"] = X.index

    # Merge with metadata
    pca_df = pca_df.merge(
        metadata,
        on="Sample"
    )

    # Plot
    plt.figure(figsize=(8, 6))

    for group in pca_df["Group"].unique():

        subset = pca_df[pca_df["Group"] == group]

        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            s=100,
            label=group
        )

    # Label samples
    for _, row in pca_df.iterrows():
        plt.text(
            row["PC1"],
            row["PC2"],
            row["Sample"],
            fontsize=8
        )

    plt.xlabel(
        f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
    )

    plt.ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
    )

    plt.title("Principal Component Analysis")

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGURE_DIR,
            "pca_plot.png"
        ),
        dpi=300
    )

    plt.close()

    print("✓ pca_plot.png")

def plot_correlation_heatmap(log_expression):

    correlation = log_expression.corr(method="pearson")

    plt.figure(figsize=(8, 6))

    plt.imshow(correlation, aspect="auto")

    plt.colorbar(label="Pearson Correlation")

    plt.xticks(
        range(len(correlation.columns)),
        correlation.columns,
        rotation=45
    )

    plt.yticks(
        range(len(correlation.columns)),
        correlation.columns
    )

    plt.title("Sample Correlation Heatmap")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGURE_DIR,
            "sample_correlation_heatmap.png"
        ),
        dpi=300
    )

    plt.close()

    print("✓ sample_correlation_heatmap.png")

def main():

    expression, metadata = load_data()

    print("Expression matrix:", expression.shape)
    print("Metadata:", metadata.shape)

    log_expression = np.log2(expression + 1)

    plot_raw_distribution(expression)
    plot_log_distribution(log_expression)
    plot_boxplot(log_expression)
    plot_pca(log_expression, metadata)
    plot_correlation_heatmap(log_expression)


if __name__ == "__main__":
    main()