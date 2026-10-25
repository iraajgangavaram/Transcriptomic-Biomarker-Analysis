# An End-to-End Python Transcriptomics Pipeline for Alzheimer's Disease Biomarker Discovery

## Overview

This project is a Python workflow for analysing publicly available RNA-seq data
to look for genes, biological pathways and interaction networks associated with
Alzheimer's disease. It compares post-mortem brain tissue from people with
Alzheimer's disease against healthy controls, using data from the NCBI Gene
Expression Omnibus (GEO).

The pipeline covers every stage from raw counts to interpretation:

- downloading and preprocessing RNA-seq count data from GEO
- quality control with expression distributions, sample correlation and PCA
- differential expression analysis with multiple-testing correction
- gene annotation from Ensembl IDs to gene symbols
- Gene Ontology and KEGG pathway enrichment
- protein-protein interaction network analysis with the STRING database
- hub gene identification from network centrality
- an interactive Streamlit dashboard for exploring the results

## Biological Question

**Which genes and regulatory pathways are altered in Alzheimer's disease compared with healthy brain tissue?**

Alzheimer's disease is the most common cause of dementia, and transcriptomic
studies of brain tissue are one way to find the molecular processes that change
with disease. This project builds a reproducible workflow for asking that
question with public data, and for generating candidate genes and pathways for
further study.

## Dataset

| | |
|---|---|
| GEO accession | [GSE163877](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163877) |
| Organism | *Homo sapiens* |
| Tissue | Post-mortem middle temporal gyrus |
| Samples | 3 Alzheimer's disease, 4 healthy control |
| Data | Simultaneous coding and non-coding RNA sequencing counts |

The middle temporal gyrus is one of the cortical regions affected early in
Alzheimer's disease, which makes it a useful tissue for studying disease-related
expression changes.

## Analysis Workflow

```text
GEO RNA-seq counts
        |
        v
Preprocessing  ->  Quality control
        |
        v
Differential expression
        |
        v
Gene annotation
        |
        +------------------+------------------+
        v                  v                  v
GO enrichment       KEGG enrichment     STRING network
                                              |
                                              v
                                     Hub gene analysis
        |
        v
Interactive dashboard
```

## Methods

### 1. Data processing

Raw counts are downloaded from GEO with `src/download_data.py`. In
`src/preprocess.py`, gene identifiers are cleaned, duplicate genes are removed
and the data is saved as a gene-by-sample expression matrix.

### 2. Quality control

`src/quality_control.py` and `notebooks/02_quality_control.ipynb` check the data
before any testing:

- raw and log-transformed expression distributions
- per-sample box plots of log counts
- a Pearson correlation heatmap between samples
- principal component analysis (PCA) of all samples

### 3. Differential expression analysis

`src/differential_expression.py` compares Alzheimer's disease with control
samples for each gene:

- log2 transformation of counts (log2(count + 1))
- log2 fold change between group means
- Welch's two-sample t-test
- Benjamini-Hochberg false discovery rate (FDR) correction

Genes with identical values in every sample, and genes with zero variance in
both groups, cannot be tested meaningfully and are excluded before the FDR
correction.

### 4. Gene annotation

`src/gene_annotation.py` maps Ensembl gene IDs to gene symbols with MyGene, so
the ranked genes can be used in enrichment and network tools.

### 5. Functional enrichment analysis

The top-ranked genes are analysed with gseapy:

- **Gene Ontology (GO) biological processes** (`src/go_enrichment.py`), to find
  processes linked to the expression changes
- **KEGG pathways** (`src/kegg_enrichment.py`), to find molecular and signalling
  pathways involved

### 6. Protein interaction network analysis

`src/string_network.py` retrieves protein-protein interactions between the
top-ranked genes from the STRING database. `src/hub_analysis.py` then builds the
network with NetworkX and ranks genes by degree and betweenness centrality to
identify hub genes, which are candidates for a regulatory role.

### 7. Interactive dashboard

`src/dashboard.py` is a Streamlit app that brings the results together:
summary statistics, the volcano plot, QC plots, and the enrichment and network
results when they have been generated.

## Results

### Quality control

All samples are highly correlated with one another (Pearson r of about 0.94 or
higher), so no sample stands out as an outlier. In the PCA, the first two
components explain 20.6% and 18.9% of the variance, and the Alzheimer's disease
and control samples do not form separate clusters. This suggests that disease
status is not the largest source of variation in these samples, which is common
in post-mortem brain tissue.

![PCA Plot](screenshots/pca_plot.png)

### Differential expression

| | Genes |
|---|---|
| Tested | 52,253 |
| Unadjusted p < 0.05 | 1,548 |
| FDR < 0.05 | 0 |

The volcano plot shows the distribution of fold changes and p-values across all
tested genes. No gene reaches FDR < 0.05 (the smallest adjusted p-value is
0.91). With 3 Alzheimer's disease and 4 control samples, the study has power to
detect only very large effects, so the ranked gene list is best treated as a
source of candidates rather than confirmed biomarkers.

![Volcano Plot](screenshots/volcano_plot.png)

My project [count-de-benchmark](https://github.com/iraajgangavaram/count-de-benchmark)
reanalyses this dataset with median-of-ratios normalisation and an
empirical-Bayes moderated t-test, and reaches the same conclusion. It also uses
simulations to show how detection power changes with the number of replicates.

### Enrichment and network analysis

Running the annotation, enrichment and STRING scripts on the ranked gene list
produces GO and KEGG enrichment tables and bar plots, an interaction network
figure and a ranked list of hub genes. These steps call online services
(MyGene, Enrichr and STRING), so they need internet access.

## Biological Interpretation

The quality control results show consistent, high-quality samples, and the
differential expression analysis gives a full ranking of genes by evidence of
change between Alzheimer's disease and control tissue. Because no gene passes
FDR correction in this small dataset, the enrichment and network results are
hypothesis-generating: they point to pathways and hub genes that are worth
testing in larger cohorts, rather than establishing them.

The same workflow can be applied directly to larger Alzheimer's disease
RNA-seq datasets, where the statistical power would be much higher.

## Running the Pipeline

Install the requirements:

```bash
pip install -r requirements.txt
```

Run the scripts from the repository root:

```bash
python src/download_data.py
python src/preprocess.py
python src/quality_control.py
python src/differential_expression.py
python src/volcano_plot.py
python src/gene_annotation.py
python src/go_enrichment.py
python src/kegg_enrichment.py
python src/string_network.py
python src/hub_analysis.py
python src/heatmap.py
```

Launch the dashboard:

```bash
streamlit run src/dashboard.py
```

## Repository Structure

```text
Transcriptomic-Biomarker-Analysis/
├── data/
│   ├── metadata/     # sample metadata and group labels
│   ├── processed/    # processed expression matrix
│   └── raw/          # downloaded counts (not tracked; see src/download_data.py)
├── figures/          # QC, differential expression, enrichment and network plots
├── notebooks/        # quality control notebook
├── results/          # differential expression results and other output tables
├── screenshots/      # images used in this README
├── src/              # pipeline scripts and Streamlit dashboard
│   └── archive/      # earlier exploratory scripts
├── requirements.txt
└── README.md
```

## Technologies Used

- **Programming:** Python
- **Data analysis:** pandas, NumPy, SciPy, statsmodels, scikit-learn
- **Bioinformatics:** NCBI GEO, gseapy (Enrichr), STRING database, MyGene
- **Visualisation:** Matplotlib, seaborn, NetworkX, Streamlit

## Limitations

- The dataset has seven samples, which limits statistical power.
- Counts are log-transformed without library-size normalisation.
- Post-mortem tissue carries variation from factors such as age, post-mortem
  interval and RNA quality, which are not modelled here.

## Future Extensions

- Applying the pipeline to larger cohorts, such as the AMP-AD studies (ROSMAP,
  Mayo, MSBB)
- Using DESeq2, edgeR or limma-voom for count-based differential expression
- Validating candidate genes across independent datasets
- Integrating clinical metadata
- Machine learning classification of disease status
- Single-cell transcriptomics integration

## Author

Iraaj Gangavaram

## Licence

MIT. See [LICENSE](LICENSE).
