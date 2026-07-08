"""
Phase 1: GEO Dataset Acquisition Pipeline

Downloads a public gene expression dataset from NCBI GEO,
extracts expression data and sample metadata,
and saves them for downstream transcriptomic analysis.
"""


import os
import GEOparse
import pandas as pd



# =====================================================
# Configuration
# =====================================================

GEO_ACCESSION = "GSE5281"

RAW_FOLDER = "data/raw"
METADATA_FOLDER = "data/metadata"



# =====================================================
# Create required folders
# =====================================================

os.makedirs(
    RAW_FOLDER,
    exist_ok=True
)

os.makedirs(
    METADATA_FOLDER,
    exist_ok=True
)



# =====================================================
# Download GEO dataset
# =====================================================

print("\nDownloading GEO dataset:", GEO_ACCESSION)


gse = GEOparse.get_GEO(
    geo=GEO_ACCESSION,
    destdir=RAW_FOLDER,
    silent=False
)


print("\nDataset downloaded successfully")



# =====================================================
# Extract expression data
# =====================================================

print("\nExtracting expression data...")


expression_list = []


for sample_id, sample in gse.gsms.items():

    print("Processing sample:", sample_id)


    table = sample.table.copy()

    table["Sample"] = sample_id


    expression_list.append(table)



expression_matrix = pd.concat(
    expression_list,
    ignore_index=True
)



expression_output = os.path.join(
    RAW_FOLDER,
    "expression_matrix.csv"
)


expression_matrix.to_csv(
    expression_output,
    index=False
)


print(
    "\nExpression matrix saved:",
    expression_output
)



# =====================================================
# Extract sample metadata
# =====================================================

print("\nExtracting metadata...")


metadata_list = []


for sample_id, sample in gse.gsms.items():

    metadata = sample.metadata


    metadata_list.append(
        {

            "Sample":
                sample_id,


            "Title":
                metadata.get(
                    "title",
                    ["Unknown"]
                )[0],


            "Source":
                metadata.get(
                    "source_name_ch1",
                    ["Unknown"]
                )[0],


            "Organism":
                metadata.get(
                    "organism_ch1",
                    ["Unknown"]
                )[0],


            "Characteristics":
                "; ".join(
                    metadata.get(
                        "characteristics_ch1",
                        ["Unknown"]
                    )
                )

        }
    )



metadata_df = pd.DataFrame(
    metadata_list
)



metadata_output = os.path.join(
    METADATA_FOLDER,
    "sample_metadata.csv"
)



metadata_df.to_csv(
    metadata_output,
    index=False
)



print(
    "Metadata saved:",
    metadata_output
)



# =====================================================
# Summary
# =====================================================

print("\n==============================")
print("PHASE 1 COMPLETE")
print("==============================")

print(
    "Number of samples:",
    len(gse.gsms)
)

print(
    "Expression shape:",
    expression_matrix.shape
)

print("\nFiles created:")

print(
    expression_output
)

print(
    metadata_output
)
