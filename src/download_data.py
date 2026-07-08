"""
Phase 2:
Download RNA-seq dataset and metadata from GEO.

Dataset:
GSE163877 - Alzheimer's disease RNA-seq study
"""

import os
import requests
import pandas as pd


# -----------------------------
# Configuration
# -----------------------------

DATASET = "GSE163877"

RAW_DIR = "data/raw"
METADATA_DIR = "data/metadata"


os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(METADATA_DIR, exist_ok=True)


# -----------------------------
# Download count matrix
# -----------------------------

counts_url = (
    "https://ftp.ncbi.nlm.nih.gov/"
    "geo/series/GSE163nnn/"
    "GSE163877/suppl/"
    "GSE163877_VBB_Counts.txt.gz"
)


counts_file = os.path.join(
    RAW_DIR,
    "GSE163877_counts.txt.gz"
)


print("Downloading count matrix...")


response = requests.get(
    counts_url,
    stream=True
)

response.raise_for_status()


with open(counts_file, "wb") as file:
    for chunk in response.iter_content(
        chunk_size=8192
    ):
        file.write(chunk)


print(
    "Saved:",
    counts_file
)


# -----------------------------
# Download GEO series matrix
# (contains sample metadata)
# -----------------------------

metadata_url = (
    "https://ftp.ncbi.nlm.nih.gov/"
    "geo/series/GSE163nnn/"
    "GSE163877/matrix/"
    "GSE163877_series_matrix.txt.gz"
)


metadata_file = os.path.join(
    METADATA_DIR,
    "GSE163877_metadata.txt.gz"
)


print("Downloading metadata...")


response = requests.get(
    metadata_url,
    stream=True
)

response.raise_for_status()


with open(metadata_file, "wb") as file:
    for chunk in response.iter_content(
        chunk_size=8192
    ):
        file.write(chunk)


print(
    "Saved:",
    metadata_file
)


# -----------------------------
# Quick validation
# -----------------------------

print("\nChecking count matrix...")


counts = pd.read_csv(
    counts_file,
    sep="\t",
    compression="gzip",
    nrows=5
)


print(counts.head())


print("\nDownload complete.")