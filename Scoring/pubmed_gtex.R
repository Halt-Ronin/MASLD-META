#### Libraries

library(easyPubMed)
library(data.table)
library(dplyr)
library(biomaRt)
library(readxl)
library(writexl)
library(org.Hs.eg.db)
library(ggplot2)

#### Data

fibrosis_top_score <- read_xlsx("fibrosis_top_scoring_genes.xlsx")
nas_top_score <- read_xlsx("nas_top_scoring_genes.xlsx")


################### PUBMED


# Function to check PubMed hits for multiple MASLD-related terms
check_nafld_related_pubmed <- function(gene) {
  query <- paste0(
    gene, 
    " AND (NAFLD OR MASLD OR MASH OR NASH OR ",
    "\"nonalcoholic fatty liver disease\" OR ",
    "\"metabolic dysfunction-associated steatotic liver disease\" OR ",
    "\"steatohepatitis\" OR \"fatty liver\" OR ", 
    "\"hepatic steatosis\" OR \"liver steatosis\" OR ",
    "\"liver fibrosis\" OR \"hepatic inflammation\")"
  )
  
  res <- tryCatch({
    search_result <- epm_query(query)
    search_result@meta$exp_count > 0
  }, error = function(e) FALSE)
  
  return(res)
}

# Apply to your gene list
associated_with_nafld_related_nas <- sapply(nas_top_score$Gene, check_nafld_related_pubmed)
associated_with_nafld_related_fibrosis <- sapply(fibrosis_top_score$Gene, check_nafld_related_pubmed)

# Get genes not associated with any of the NAFLD-related terms
non_nafld_related_genes_nas <- nas_top_score[!associated_with_nafld_related_nas,]
non_nafld_related_genes_fibrosis <- fibrosis_top_score[!associated_with_nafld_related_fibrosis,]

nas_top_score$pubmed_mentioned <- associated_with_nafld_related_nas
fibrosis_top_score$pubmed_mentioned <- associated_with_nafld_related_fibrosis


write_xlsx(nas_top_score, "pubmed_column_nas.xlsx")
write_xlsx(fibrosis_top_score, "pubmed_colum_fibrosis.xlsx")

####################################

file_path <- "GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
gtex_df <- fread(file_path, skip = 2)

gtex_df$gene_id <- sub("\\..*", "", gtex_df$Name)
expr_data <- as.data.frame(gtex_df[, -(1:2)])  # Drop Name + Description

expr_mat <- aggregate(expr_data, by = list(ENSG = gtex_df$gene_id), FUN = mean)
rownames(expr_mat) <- expr_mat$ENSG
expr_mat <- expr_mat[, -1]

annots <- select(org.Hs.eg.db, keys=rownames(expr_mat), 
                 columns="SYMBOL", keytype="ENSEMBL")
na.omit(annots)
annots <- annots[!duplicated(annots$ENSEMBL),] 
annots <- annots[!duplicated(annots$SYMBOL),] 
annots <- annots[!is.na(annots$SYMBOL),] 

expr_mat <- expr_mat[rownames(expr_mat) %in% annots$ENSEMBL,]
annots <- annots[annots$ENSEMBL %in% rownames(expr_mat),]
sum(rownames(expr_mat) != annots$ENSEMBL)

rownames(expr_mat) <- annots$SYMBOL


liver_col <- grep("^Liver", colnames(expr_mat), value = TRUE)
if (length(liver_col) == 0) stop("Liver tissue not found in GTEx columns.")



#genes <- fibrosis_top_score$Gene
genes <- nas_top_score$Gene

result_df <- data.frame(
  gene = genes,
  category = NA_integer_
)

for (i in seq_along(genes)) {
  gene <- genes[i]
  
  if (!(gene %in% rownames(expr_mat))) {
    result_df$category[i] <- NA  # Not found in GTEx
    next
  }
  
  expr_values <- as.numeric(expr_mat[gene, ])
  names(expr_values) <- colnames(expr_mat)
  liver_expr <- expr_values[liver_col]
  
  if (all(is.na(liver_expr)) || sum(liver_expr, na.rm = TRUE) == 0) {
    result_df$category[i] <- 0  # Not expressed in liver
  } else {
    expressed_tissues <- names(expr_values)[expr_values > 0]
    max_tissue <- names(which.max(expr_values))
    
    if (length(expressed_tissues) == 1 && expressed_tissues %in% liver_col) {
      result_df$category[i] <- 3  # Only expressed in liver
    } else if (max_tissue %in% liver_col) {
      result_df$category[i] <- 2  # Highest in liver
    } else {
      result_df$category[i] <- 1  # Expressed in liver, but not highest
    }
  }
}

result_df <- result_df %>%
  mutate(category = case_when(
    category == 3 ~ "Only expressed in liver",
    category == 2 ~ "Highest in Liver",
    category == 1 ~ "Expressed in Liver",
    category == 0 ~ "Not Expressed in Liver",
    is.na(category) ~ "Not found"
  ))

#print(sum(result_df$Gene != fibrosis_top_score$Gene))
#fibrosis_top_score$GTEX <- result_df$category

print(sum(result_df$Gene != nas_top_score$Gene))
nas_top_score$GTEX <- result_df$category

#########

write_xlsx(nas_top_score,"nas_pubmed_gtex.xlsx")
write_xlsx(fibrosis_top_score,"fibross_pubmed_gtex.xlsx")


########
