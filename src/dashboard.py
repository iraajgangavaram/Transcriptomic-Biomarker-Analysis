"""
Phase 7
Interactive Transcriptomic Analysis Dashboard
"""

import streamlit as st
import pandas as pd
from PIL import Image


# -----------------------------
# Page Setup
# -----------------------------

st.set_page_config(
    page_title="Transcriptomic Biomarker Analysis",
    layout="wide"
)


# -----------------------------
# Title
# -----------------------------

st.title(
    "Integrative Transcriptomic Analysis of Alzheimer's Disease"
)

st.write(
    """
    This dashboard presents results from RNA-seq analysis including:
    
    - Quality control
    - Differential expression analysis
    - Functional enrichment
    - Protein interaction network analysis
    """
)


# -----------------------------
# Load Results
# -----------------------------

@st.cache_data
def load_data():

    de_results = pd.read_csv(
        "results/differential_expression_results.csv"
    )

    go_results = pd.read_csv(
        "results/go_enrichment.csv"
    )

    kegg_results = pd.read_csv(
        "results/kegg_enrichment.csv"
    )

    hub_results = pd.read_csv(
        "results/hub_genes.csv"
    )

    return (
        de_results,
        go_results,
        kegg_results,
        hub_results
    )


de, go, kegg, hubs = load_data()


# -----------------------------
# Differential Expression
# -----------------------------

st.header(
    "Differential Expression"
)


col1, col2 = st.columns(2)


with col1:

    st.metric(
        "Genes analysed",
        len(de)
    )


with col2:

    st.metric(
        "Significant genes",
        len(
            de[
                de["padj"] < 0.05
            ]
        )
    )


st.subheader(
    "Volcano Plot"
)


st.image(
    "figures/volcano_plot.png"
)


# -----------------------------
# PCA
# -----------------------------

st.header(
    "Sample Separation"
)


st.image(
    "figures/pca_plot.png"
)


# -----------------------------
# Heatmap
# -----------------------------

st.header(
    "Top Differentially Expressed Genes"
)


st.image(
    "figures/top_genes_heatmap.png"
)


# -----------------------------
# Enrichment
# -----------------------------

st.header(
    "Functional Enrichment"
)


tab1, tab2 = st.tabs(
    [
        "GO Biological Processes",
        "KEGG Pathways"
    ]
)


with tab1:

    st.image(
        "figures/go_enrichment_barplot.png"
    )

    st.dataframe(
        go.head(20)
    )


with tab2:

    st.image(
        "figures/kegg_enrichment_barplot.png"
    )

    st.dataframe(
        kegg.head(20)
    )


# -----------------------------
# Network Analysis
# -----------------------------

st.header(
    "Protein Interaction Network"
)


st.image(
    "figures/interaction_network.png"
)


st.subheader(
    "Top Hub Genes"
)


st.dataframe(
    hubs.head(20)
)


st.success(
    "Dashboard loaded successfully"
)