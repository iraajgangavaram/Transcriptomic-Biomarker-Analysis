"""
Phase 5.1
Convert Ensembl Gene IDs to Gene Symbols
"""

import os
import pandas as pd
import mygene


INPUT_FILE = "results/differential_expression_results.csv"

OUTPUT_DIR = "results"

OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    "differential_expression_gene_symbols.csv"
)


def main():

    results = pd.read_csv(
        INPUT_FILE
    )


    gene_ids = results["Gene_ID"].tolist()


    print(
        "Annotating genes..."
    )


    mg = mygene.MyGeneInfo()


    annotations = mg.querymany(
        gene_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human"
    )


    annotation_df = pd.DataFrame(
        annotations
    )


    annotation_df = annotation_df[
        ["query", "symbol"]
    ]


    annotation_df.columns = [
        "Gene_ID",
        "Gene_Symbol"
    ]


    results = results.merge(
        annotation_df,
        on="Gene_ID",
        how="left"
    )


    results = results.dropna(
        subset=["Gene_Symbol"]
    )


    results.to_csv(
        OUTPUT_FILE,
        index=False
    )


    print(
        "Saved:",
        OUTPUT_FILE
    )


    print(
        "Genes annotated:",
        len(results)
    )


if __name__ == "__main__":
    main()