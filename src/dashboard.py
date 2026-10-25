"""
Streamlit dashboard for the GSE163877 analysis.

Run from the repository root:
    streamlit run src/dashboard.py
"""

from pathlib import Path

import pandas as pd
import streamlit as st

RESULTS = Path("results")
FIGURES = Path("figures")


def show_image(path, caption=None):
    path = Path(path)
    if path.exists():
        st.image(str(path), caption=caption)
    else:
        st.info(f"{path} has not been generated.")


def optional_table(path):
    path = Path(path)
    return pd.read_csv(path) if path.exists() else None


st.set_page_config(page_title="Transcriptomic Biomarker Analysis", layout="wide")
st.title("Integrative Transcriptomic Analysis of Alzheimer's Disease")
st.write(
    "Middle temporal gyrus, 3 Alzheimer's disease and 4 control samples. "
    "With this few samples no gene passes FDR < 0.05, so everything below "
    "is exploratory."
)

de = pd.read_csv(RESULTS / "differential_expression_results.csv")
tested = de["padj"].notna()

st.header("Differential expression")
col1, col2, col3 = st.columns(3)
col1.metric("Genes tested", int(tested.sum()))
col2.metric("Unadjusted p < 0.05", int((de["pvalue"] < 0.05).sum()))
col3.metric("FDR < 0.05", int((de["padj"] < 0.05).sum()))
st.caption(f"Smallest adjusted p-value: {de['padj'].min():.2f}")
show_image(FIGURES / "volcano_plot.png")

st.header("Quality control")
show_image(FIGURES / "pca_plot.png", "PCA of samples")
show_image(FIGURES / "sample_correlation_heatmap.png", "Sample correlation")

go = optional_table(RESULTS / "go_enrichment.csv")
kegg = optional_table(RESULTS / "kegg_enrichment.csv")
hubs = optional_table(RESULTS / "hub_genes.csv")

if go is None and kegg is None and hubs is None:
    st.header("Enrichment and network analysis")
    st.info(
        "These outputs were removed because they were built from a gene "
        "ranking produced before the FDR bug was fixed. Rerun "
        "gene_annotation.py, go_enrichment.py, kegg_enrichment.py and "
        "string_network.py to regenerate them as exploratory results."
    )
else:
    st.header("Functional enrichment (exploratory)")
    if go is not None:
        show_image(FIGURES / "go_enrichment_barplot.png")
        st.dataframe(go.head(20))
    if kegg is not None:
        show_image(FIGURES / "kegg_enrichment_barplot.png")
        st.dataframe(kegg.head(20))
    if hubs is not None:
        st.header("Protein interaction network (exploratory)")
        show_image(FIGURES / "interaction_network.png")
        st.dataframe(hubs.head(20))
