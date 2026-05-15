# MASLD-META analysis code

This repository contains the analysis code used for the MASLD-META transcriptomic meta-analysis. The analysis integrates publicly available bulk RNA-seq liver biopsy datasets with histological metadata to identify genes, pathways, and regulatory programs associated with MASLD progression across NAFLD Activity Score (NAS) and fibrosis stage.

The code is organized as a set of R Markdown and R scripts covering sample grouping, differential expression analysis, regression analysis, gene disease pattern analysis, gene scoring, over-representation analysis, temporal pathway analysis, and sex-stratified analysis.

## Project overview

The study analyzes nine GEO bulk RNA-seq datasets with MASLD liver biopsy samples and available histological annotations. The main analyses compare gene expression across NAS groups and fibrosis stages, then integrate evidence from multiple approaches to prioritize progression-associated genes.

The main analytical layers are:

1. NAS grouping using PCA and PERMANOVA
2. Individual-study differential gene expression analysis with DESeq2
3. Meta-analysis of differential expression p-values with metaRNASeq
4. Global differential expression analysis after ComBat-seq batch correction
5. Linear regression against NAS group or fibrosis score
6. Gene disease pattern analysis using GAMMA rank correlation and self-organizing maps
7. Frequency-based integration of all evidence into a total gene score
8. MSigDB and DoRothEA over-representation analysis
9. Temporal and sex-stratified enrichment analyses

## Repository structure

```text
MASLD META/
├── NAS Groupings/
│   └── NAS_Clustering.Rmd
├── DGEA/
│   ├── Individual/
│   │   ├── NAS_individual_analysis.Rmd
│   │   └── Fibrosis_individual_analysis.Rmd
│   └── Global/
│       ├── NAS_global_analysis.Rmd
│       └── Fibrosis_global_analysis.Rmd
├── Regression/
│   ├── NAS_Regression.Rmd
│   └── Fibrosis_Regression.Rmd
├── GeneDiseasePaterns/
│   ├── genediseasePatterns_workflow_NAS.Rmd
│   └── genediseasePatterns_workflow_Fibrosis.Rmd
├── Scoring/
│   ├── Merge_Analysis.Rmd
│   ├── Score.Rmd
│   ├── pubmed_gtex.R
│   └── tcga.R
├── ORA/
│   ├── NAS_ORA.Rmd
│   ├── Fibrosis_ORA.Rmd
│   ├── Sex_Stratified_ORA.Rmd
│   └── Temporal_ORA_network.R
└── Sex Stratified Analysis/
    ├── NAS_gender.Rmd
    └── Fibrosis_gender.Rmd
```

## Expected input files

The repository assumes that several input files are available locally in the working directory. These files are not included in this code-only repository.

### Raw count matrices

Each GEO dataset is expected as a tab-delimited count matrix named as:

```text
GSE135251.tsv
GSE193066.tsv
GSE225740.tsv
GSE130970.tsv
GSE192959.tsv
GSE207310.tsv
GSE162694.tsv
GSE174478.tsv
GSE185051.tsv
```

Most scripts expect genes as rows and samples as columns. The first column should contain gene identifiers.

### GEO metadata

Metadata are retrieved inside the scripts with:

```r
GEOquery::getGEO(id, GSEMatrix = TRUE)
```

Several scripts use dataset-specific metadata column indices for NAS, fibrosis stage, or sex. These indices should be checked if GEO metadata formatting changes.

## R package requirements

The scripts use the following main R packages:

```r
set.seed(42)

# Data handling
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

# File input and output
library(readxl)
library(writexl)
library(rmarkdown)

# GEO and expression data infrastructure
library(GEOquery)
library(Biobase)
library(TCGAbiolinks)

# Gene annotation and ID conversion
library(org.Hs.eg.db)
library(biomaRt)

# RNA seq differential expression and batch correction
library(DESeq2)
library(edgeR)
library(limma)
library(sva)

# Meta analysis
library(metaRNASeq)

# Disease progression and pattern analysis
library(rococo)
library(kohonen)
library(phenopath)

# Statistical analysis
library(vegan)
library(transport)

# Parallelization and progress tracking
library(parallel)
library(pbapply)

# Functional enrichment and pathway analysis
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(dorothea)
library(decoupleR)
library(viper)
library(OmnipathR)

# Literature mining
library(easyPubMed)

# Network analysis and visualization
library(igraph)
library(ggraph)

# Plotting and visualization
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(plotly)
library(pheatmap)
library(Rtsne)
```

## Citation

Manuscript not available yet.

## License

License not available yet.

## Contact

For questions about the analysis workflow, contact the corresponding or first author listed in the manuscript or open an issue in the GitHub repository.

