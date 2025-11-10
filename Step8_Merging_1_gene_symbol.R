# Load necessary library
library(dplyr)
library(tidyverse)
BiocManager::install("biomaRt")
library(biomaRt)
library(data.table)

#-------------------Read and combine all QTL files-------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results"
qtl_files <- list.files(folder_path, pattern = "_QTL\\.csv$", full.names = TRUE)
coloc_QTL <- lapply(qtl_files, read.csv, stringsAsFactors = FALSE) %>%
  bind_rows()
head(coloc_QTL)
#-------------------Add gene information-------------------
# Create gene ID column without version number
coloc_QTL$ensembl_id_no_version <- sub("\\..*", "", coloc_QTL$Gene_Ensembl)

# Connect to BioMart database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query gene information
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = unique(coloc_QTL$ensembl_id_no_version),
  mart = mart
)

# Rename columns for matching
names(genes)[names(genes) == "ensembl_gene_id"] <- "ensembl_id_no_version"

# Merge data (keep all original rows)
coloc_QTL <- merge(
  coloc_QTL, 
  genes[, c('ensembl_id_no_version', 'hgnc_symbol', 'entrezgene_id')],
  by = "ensembl_id_no_version",
  all.x = TRUE
)

# Rename new column to Gene_symbol
names(coloc_QTL)[names(coloc_QTL) == "hgnc_symbol"] <- "Gene_symbol"

# Remove temporary column
coloc_QTL$ensembl_id_no_version <- NULL

# View results
head(coloc_QTL[, c('Gene_Ensembl', 'Gene_symbol')])

write.csv(coloc_QTL, 'D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/coloc_QTL_symbol.csv')
