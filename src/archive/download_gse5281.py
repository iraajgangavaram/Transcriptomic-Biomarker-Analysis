"""
Phase 1:
Download gene expression data from GEO.

Functions:
- Download GEO dataset
- Extract expression matrix
- Extract sample metadata
- Save files for downstream analysis
"""


import os
import GEOparse
import pandas as pd



# -----------------------------
# Dataset configuration
# -----------------------------

GEO_ACCESSION = "GSE5281"

RAW_FOLDER = "data/raw"

METADATA_FOLDER = "data/metadata"



# -----------------------------
# Create folders
# -----------------------------

os.makedirs(
    RAW_FOLDER,
    exist_ok=True
)

os.makedirs(
    METADATA_FOLDER,
    exist_ok=True
)



# -----------------------------
# Download GEO dataset
# -----------------------------

print(
    f"Downloading dataset {GEO_ACCESSION}"
)


gse = GEOparse.get_GEO(
    geo=GEO_ACCESSION,
    destdir=RAW_FOLDER
)


print(
    "Dataset downloaded successfully"
)



# -----------------------------
# Extract expression data
# -----------------------------

expression_tables = []


for sample_id, sample in gse.gsms.items():

    print(
        "Processing:",
        sample_id
    )

    table = sample.table

    table["Sample"] = sample_id

    expression_tables.append(
        table
    )



expression_matrix = pd.concat(
    expression_tables,
    ignore_index=True
)



expression_matrix.to_csv(
    "data/raw/expression_matrix.csv",
    index=False
)


print(
    "Expression matrix saved"
)



# -----------------------------
# Extract metadata
# -----------------------------

metadata = []


for sample_id, sample in gse.gsms.items():

    metadata.append(
        {
            "Sample": sample_id,

            "Title":
                sample.metadata["title"][0],

            "Organism":
                sample.metadata["organism_ch1"][0],

            "Source":
                sample.metadata["source_name_ch1"][0]
        }
    )



metadata_df = pd.DataFrame(
    metadata
)



metadata_df.to_csv(
    "data/metadata/sample_metadata.csv",
    index=False
)


print(
    "Metadata saved"
)



print(
    "Phase 1 complete!"
)