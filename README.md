# An End-to-End Python Transcriptomics Pipeline for Alzheimer's Disease Biomarker Discovery

A Python pipeline for analysing public RNA-seq data from Alzheimer's disease and
control brain tissue (GEO GSE163877): quality control, differential expression,
gene annotation, pathway enrichment, protein interaction networks and an
interactive dashboard.

## Dataset

| | |
|---|---|
| GEO accession | [GSE163877](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163877) |
| Organism | *Homo sapiens* |
| Tissue | Post-mortem middle temporal gyrus |
| Samples | 3 Alzheimer's disease, 4 control |
| Data | Coding and non-coding RNA-seq counts |

## Result

| | Genes |
|---|---|
| Tested | 52,253 |
| Unadjusted p < 0.05 | 1,548 |
| Benjamini-Hochberg FDR < 0.05 | 0 |

The smallest adjusted p-value is 0.91. With 52,253 tests you would expect
around 2,600 genes below p = 0.05 by chance alone, so the unadjusted count is
consistent with chance. With 3 versus 4 samples the study has power to detect
only very large effects, so the results here are exploratory.

![Volcano plot](screenshots/volcano_plot.png)

My separate project
[count-de-benchmark](https://github.com/iraajgangavaram/count-de-benchmark)
reanalyses this dataset with median-of-ratios normalisation and an
empirical-Bayes moderated t-test, with the same outcome: no gene passes 10% FDR and the
p-value histogram is close to flat. That project also uses simulations to show
how many replicates you need before a moderated test starts finding true
changes.

## Statistical notes

- **Undefined tests.** 6,534 genes have identical values in all samples, so
  their t-test p-value is undefined. They are excluded before the
  Benjamini-Hochberg correction.
- **Zero-variance genes.** Genes with zero variance in both groups (for example
  counts 2, 2, 2 against 0, 0, 0, 0) would get a p-value of exactly 0. These 38
  genes are treated as untestable.
- **Enrichment and network steps.** GO, KEGG, STRING and hub-gene analyses take
  the top-ranked genes from the differential expression step. As no gene passes
  FDR, their output is best read as hypothesis-generating.

## Pipeline

| Step | Script | Output |
|---|---|---|
| Download | `src/download_data.py` | Raw counts in `data/raw/` (not tracked) |
| Preprocess | `src/preprocess.py` | `data/processed/expression_matrix.csv` |
| Quality control | `src/quality_control.py`, `notebooks/02_quality_control.ipynb` | Distribution, correlation and PCA plots in `figures/` |
| Differential expression | `src/differential_expression.py` | `results/differential_expression_results.csv` |
| Volcano plot | `src/volcano_plot.py` | `figures/volcano_plot.png` |
| Gene symbols | `src/gene_annotation.py` | Gene-symbol table for the steps below |
| Enrichment | `src/go_enrichment.py`, `src/kegg_enrichment.py` | GO and KEGG tables and bar plots |
| Network | `src/string_network.py`, `src/hub_analysis.py` | STRING interactions and hub genes |
| Heatmap | `src/heatmap.py` | Top-ranked genes across samples |
| Dashboard | `src/dashboard.py` | Streamlit summary of whatever outputs exist |

The differential expression step uses Welch's t-test on log2(count + 1) for
each gene, with Benjamini-Hochberg correction. The annotation, enrichment and
STRING steps call online services, so they need internet access.

![PCA plot](screenshots/pca_plot.png)

## Running it

```bash
pip install -r requirements.txt
python src/download_data.py
python src/preprocess.py
python src/differential_expression.py
python src/volcano_plot.py
streamlit run src/dashboard.py
```

Run the scripts from the repository root.

## Limitations

- Seven samples limit statistical power.
- Counts are log-transformed without library-size normalisation; size factors
  or TMM (as in DESeq2, edgeR and count-de-benchmark) are a natural extension.
- Variance estimates from 3 versus 4 samples are noisy; moderated tests such as
  limma or DESeq2 address this.
- Post-mortem tissue brings confounders (age, post-mortem interval, RNA
  quality) that this analysis does not model.

## Future work

Larger AD RNA-seq cohorts such as the AMP-AD
studies (ROSMAP, Mayo, MSBB) have hundreds of samples, and analysing one of
them with DESeq2 or limma-voom would be the natural next step.

## Repository layout

```text
data/metadata/     sample groups and GEO metadata
data/processed/    expression matrix
figures/           QC plots and volcano plot
notebooks/         quality control notebook
results/           differential expression tables
screenshots/       images used in this README
src/               pipeline scripts and dashboard
src/archive/       earlier exploratory scripts
```

## Licence

MIT. See [LICENSE](LICENSE).
