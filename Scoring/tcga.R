if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(DESeq2)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)

query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
  )

options(timeout = 1000)

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 10
)


data <- GDCprepare(query)

saveRDS(data, "data.RDS")


########
data <- readRDS("data.RDS")


counts <- assay(data)

sample_type <- data$shortLetterCode
table(sample_type)

keep_samples <- sample_type != "TR"

# subset WITHOUT changing order
counts <- counts[, keep_samples]

sample_type <- sample_type[keep_samples]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(sample_type),
  design = ~ sample_type
)

# remove low counts
dds <- dds[rowSums(counts(dds) >= 10) >= 50, ]

dds <- DESeq(dds)

res <- results(dds, contrast=c("sample_type","TP","NT"))

res <- as.data.frame(res)

gene_names <- mapIds(
  org.Hs.eg.db,
  keys = sub("\\..*", "", rownames(res)),
  column = "SYMBOL", 
  keytype = "ENSEMBL",
  multiVals = "first"
)

res$SYMBOL <- gene_names

res_filtered <- res[res$padj < 0.05,]


####### our lists
gene_names <- mapIds(
  org.Hs.eg.db,
  keys = sub("\\..*", "", rownames(counts)),
  column = "SYMBOL", 
  keytype = "ENSEMBL",
  multiVals = "first"
)


nas <- read_xlsx("patent_merged_nas.xlsx")
fibrosis <- read_xlsx("patent_merged_fibrosis.xlsx")
available_genes_tcga_nas <- nas[nas$Gene %in% gene_names,]$Gene
not_available_genes_tcga_nas <- nas[!nas$Gene %in% gene_names,]$Gene
available_genes_tcga_fibrosis <-  fibrosis[fibrosis$Gene %in% gene_names,]$Gene
not_available_genes_tcga_fibrosis <-  fibrosis[!fibrosis$Gene %in% gene_names,]$Gene
  
res_filtered_nas <- res[res$SYMBOL %in% nas$Gene, c("log2FoldChange","padj","SYMBOL")]
res_filtered_fibrosis <- res[res$SYMBOL %in% fibrosis$Gene, c("log2FoldChange","padj","SYMBOL")]

merged_df_nas <- nas %>%
  left_join(
    res_filtered_nas,
    by = c("Gene" = "SYMBOL")
  ) %>%
  mutate(
    across(
      -Gene,
      ~ ifelse(is.na(.), "Not Found in TCGA data", as.character(.))
    )
  )

colnames(merged_df_nas)[c(17,18)] <- c("TCGA_log2FC","TCGA_padj")

write_xlsx(merged_df_nas,"patent_merged_nas_tcga.xlsx")

merged_df_fibrosis <- fibrosis %>%
  left_join(
    res_filtered_fibrosis,
    by = c("Gene" = "SYMBOL")
  ) %>%
  mutate(
    across(
      -Gene,
      ~ ifelse(is.na(.), "Not Found in TCGA data", as.character(.))
    )
  )

colnames(merged_df_fibrosis)[c(17,18)] <- c("TCGA_log2FC","TCGA_padj")

write_xlsx(merged_df_fibrosis,"patent_merged_fibrosis_tcga.xlsx")




