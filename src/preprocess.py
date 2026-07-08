"""
Phase 2.1:
Preprocess RNA-seq count matrix.

Tasks:
- Load raw count data
- Clean gene identifiers
- Remove unwanted columns
- Save processed expression matrix
"""

import os
import pandas as pd


# -----------------------------
# File paths
# -----------------------------

input_file = "data/raw/GSE163877_counts.txt.gz"

output_dir = "data/processed"

output_file = os.path.join(
    output_dir,
    "expression_matrix.csv"
)


os.makedirs(
    output_dir,
    exist_ok=True
)


# -----------------------------
# Load count matrix
# -----------------------------

print("Loading count matrix...")


counts = pd.read_csv(
    input_file,
    sep="\t",
    compression="gzip"
)


print("Original shape:")
print(counts.shape)


# -----------------------------
# Inspect columns
# -----------------------------

print("\nColumns:")
print(counts.columns.tolist())


# -----------------------------
# Clean gene identifiers
# -----------------------------

# First column contains gene IDs
gene_column = counts.columns[0]


counts = counts.rename(
    columns={
        gene_column: "Gene_ID"
    }
)


# Remove duplicate genes if present
counts = counts.drop_duplicates(
    subset="Gene_ID"
)


# Set gene IDs as index

counts = counts.set_index(
    "Gene_ID"
)


# -----------------------------
# Save processed matrix
# -----------------------------

counts.to_csv(
    output_file
)


print("\nProcessed matrix saved:")
print(output_file)


print("\nFinal shape:")
print(counts.shape)