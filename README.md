# Immune Deconvolution & Gene-Immune Correlation Explorer

A Shiny app for exploring the immune cell composition of tumor samples from gene expression data. It runs multiple immune-cell deconvolution methods, compares results across cancer types and methods, and correlates gene expression with immune cell fractions — all through an interactive UI.

## Features

- **Multiple data sources**
  - Synthetic demo data (instant, no setup required)
  - Real TCGA RNA-seq data via [`TCGAbiolinks`](https://bioconductor.org/packages/TCGAbiolinks/)
  - Your own uploaded expression matrix (CSV/TSV, including cBioPortal-style RSEM files)

- **Nine deconvolution methods** via [`immunedeconv`](https://github.com/icbi-lab/immunedeconv):
  - quanTIseq, TIMER, CIBERSORT, CIBERSORT (absolute), MCP-counter, xCell, EPIC, ABIS, ConsensusTME, ESTIMATE
  - Falls back to a small built-in NNLS demo method if `immunedeconv` isn't installed

- **Cross-cancer and cross-method comparison**
  - Stacked bar charts of average immune composition by cancer type and method
  - Boxplots comparing a given cell type across cancer types
  - Method-agreement scatter plots

- **Gene ↔ immune correlation analysis**
  - Enter any gene(s) and compute Spearman correlations against selected immune cell fractions
  - Auto-generated scatter plots with correlation coefficients and p-values

- **Export results** to Excel (`.xlsx`) — both deconvolution results and correlation tables

## Installation

```r
install.packages(c(
  "shiny", "DT", "dplyr", "tidyr", "tibble", "stringr",
  "purrr", "ggplot2", "openxlsx", "shinyjs"
))

remotes::install_github("icbi-lab/immunedeconv")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
```

### Optional: CIBERSORT

CIBERSORT requires registration at [cibersortx.stanford.edu](https://cibersortx.stanford.edu). Once you have access:

1. Download `CIBERSORT.R` and `LM22.txt`
2. Place both files in a `cibersort_files/` folder next to `app.R`

Without these files, CIBERSORT and CIBERSORT (absolute) will be skipped automatically.

## Usage

```r
shiny::runApp("app.R")
```

1. Choose a data source (demo data works out of the box).
2. Select one or more cancer types and up to 3 deconvolution methods, then click **Run deconvolution**.
3. Explore results in **Deconvolution results** and **Compare cancer types & methods**.
4. Go to **Gene ↔ Immune correlations**, enter genes (e.g. `CCL2`), pick cell types, and click **Compute correlations**.

## Notes

- The synthetic demo dataset is for illustration only — it is not real biological expression data.
- TIMER and ConsensusTME require a per-sample cancer-type label using supported lowercase TCGA codes (e.g. `brca`, `luad`, `skcm`).
- Uploads without a metadata file are grouped under a `"Not specified"` placeholder; TIMER and ConsensusTME are skipped in that case.

## Citations

- Newman AM, et al. Robust enumeration of cell subsets from tissue expression profiles (CIBERSORT). *Nat Methods*. 2015;12(5):453-457.
- Sturm G, et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology (immunedeconv). *Bioinformatics*. 2019;35(14):i436-i445.
- Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature*. 2012;490(7418):61-70. Data from the [TCGA Research Network](https://www.cancer.gov/tcga).

See the **Method notes** tab in the app for the full list of method citations.
