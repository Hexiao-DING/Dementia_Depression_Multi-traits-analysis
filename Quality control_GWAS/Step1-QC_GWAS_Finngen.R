# ============================
# Step 0. Install and load required packages
# ============================
install.packages("httr")
library(data.table)
library(dplyr)
library(httr)
library(rtracklayer)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

setwd("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS")

# ===============================================================
# Step 1. Process Undefined Dementia (UD) dataset (hg19)
# ===============================================================

ud <- fread("Downloaded_original_data/summary_stats_release_finngen_R12_F5_DEMNAS.gz")
colSums(is.na(ud))
colnames(ud)

# Standardize column names
name_mapping <- c(
  rsids = "SNP",
  '#chrom' = "CHR",
  pos = "BP",
  beta = "beta",
  sebeta = "se",
  pval = "pval",
  af_alt = "FRQ"
)
names(ud) <- name_mapping[names(ud)]

# Add sample size and calculate derived columns
ud$N <- 500348
ud$MAF <- pmin(ud$FRQ, 1 - ud$FRQ)
ud$Z <- ud$beta / ud$se

# Keep only relevant columns
ud <- ud %>%
  dplyr::select(SNP, CHR, BP, alt, ref, beta, se, Z, pval, N, FRQ, MAF)

# Format summary statistics and save as GRCh38
format_sumstats(
  path = ud,
  ref_genome = "GRCh38",
  dbSNP = 144,
  save_path = "Format_sumstats_data/UD38.tsv.gz"
)

# Reload formatted data
ud <- read.table("Format_sumstats_data/UD38.tsv.gz", header = TRUE, sep = "\t")
colSums(is.na(ud))

# ============================
# Step 1.1. LiftOver: hg38 to hg19
# ============================

chain <- import.chain("D:/hg38ToHg19.over.chain/hg38ToHg19.over.chain")

gr_hg38 <- GRanges(
  seqnames = paste0("chr", ud$CHR),
  ranges = IRanges(start = ud$BP, end = ud$BP),
  strand = "*",
  mcols = ud
)

gr_hg19 <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(gr_hg19)

converted_data <- as.data.frame(gr_hg19) %>%
  dplyr::rename(
    CHR = seqnames,
    BP_hg19 = start,
    BP_hg38 = mcols.BP
  ) %>%
  mutate(CHR = as.integer(gsub("chr", "", CHR)))

final_ud <- ud %>%
  inner_join(converted_data, by = c("CHR", "BP" = "BP_hg38")) %>%
  select(-BP) %>%
  dplyr::rename(BP = "BP_hg19")

failed_ud <- ud %>%
  anti_join(converted_data, by = c("CHR", "BP" = "BP_hg38"))

colnames(final_ud)
colSums(is.na(final_ud))

# Rename columns for downstream analysis
final_ud <- dplyr::rename(final_ud, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

# Prepare MTAG and LDSC formatted datasets
ud_mtag <- final_ud %>%
  dplyr::select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
ud_ldsc <- final_ud %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

# Export results
write.table(ud_mtag, "MTAG_input_data/MTAGUDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ud_ldsc, "LDSC_input_data/LDSCUDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ud_ldsc, "Final_processed_data/UDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)



# ===============================================================
# Step 2. Process Vascular Dementia (VD) dataset (hg38)
# ===============================================================

vd <- fread("Downloaded_original_data/summary_stats_release_finngen_R12_F5_VASCDEM.gz")
colSums(is.na(vd))
colnames(vd)

# Standardize names
name_mapping <- c(
  rsids = "SNP",
  '#chrom' = "CHR",
  pos = "BP",
  beta = "beta",
  sebeta = "se",
  pval = "pval",
  af_alt = "FRQ"
)
names(vd) <- name_mapping[names(vd)]

vd$N <- 500348
vd$MAF <- pmin(vd$FRQ, 1 - vd$FRQ)
vd$Z <- vd$beta / vd$se

vd <- vd %>%
  dplyr::select(SNP, CHR, BP, alt, ref, beta, se, Z, pval, N, FRQ, MAF)

# Format and save as GRCh38
format_sumstats(
  path = vd,
  ref_genome = "GRCh38",
  dbSNP = 144,
  save_path = "Format_sumstats_data/VD38.tsv.gz"
)

vd <- read.table("Format_sumstats_data/VD38.tsv.gz", header = TRUE, sep = "\t")

# LiftOver: hg38 to hg19
gr_hg38 <- GRanges(
  seqnames = paste0("chr", vd$CHR),
  ranges = IRanges(start = vd$BP, end = vd$BP),
  strand = "*",
  mcols = vd
)
gr_hg19 <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(gr_hg19)

converted_data <- as.data.frame(gr_hg19) %>%
  dplyr::rename(CHR = seqnames, BP_hg19 = start, BP_hg38 = mcols.BP) %>%
  mutate(CHR = as.integer(gsub("chr", "", CHR)))

final_vd <- vd %>%
  inner_join(converted_data, by = c("CHR", "BP" = "BP_hg38")) %>%
  select(-BP) %>%
  dplyr::rename(BP = "BP_hg19")

final_vd <- dplyr::rename(final_vd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

# Prepare and export
vd_mtag <- final_vd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
vd_ldsc <- final_vd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(vd_mtag, "MTAG_input_data/MTAGVDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(vd_ldsc, "LDSC_input_data/LDSCVDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(vd_ldsc, "Final_processed_data/VDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)



# ===============================================================
# Step 3. Process Overall Dementia dataset (Dementia, hg38)
# ===============================================================

dementia <- fread("Downloaded_original_data/summary_stats_release_finngen_R12_F5_DEMENTIA.gz")
colSums(is.na(dementia))
colnames(dementia)

# Rename columns
name_mapping <- c(
  rsids = "SNP",
  '#chrom' = "CHR",
  pos = "BP",
  beta = "beta",
  sebeta = "se",
  pval = "pval",
  af_alt = "FRQ"
)
names(dementia) <- name_mapping[names(dementia)]

dementia$N <- 500348
dementia$MAF <- pmin(dementia$FRQ, 1 - dementia$FRQ)
dementia$Z <- dementia$beta / dementia$se

dementia <- dementia %>%
  dplyr::select(SNP, CHR, BP, alt, ref, beta, se, Z, pval, N, FRQ, MAF)

# Format and save
format_sumstats(
  path = dementia,
  ref_genome = "GRCh38",
  dbSNP = 144,
  save_path = "Format_sumstats_data/DEM38.tsv.gz"
)

dementia <- read.table("Format_sumstats_data/DEM38.tsv.gz", header = TRUE, sep = "\t")

# LiftOver to hg19
gr_hg38 <- GRanges(
  seqnames = paste0("chr", dementia$CHR),
  ranges = IRanges(start = dementia$BP, end = dementia$BP),
  strand = "*",
  mcols = dementia
)
gr_hg19 <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(gr_hg19)

converted_data <- as.data.frame(gr_hg19) %>%
  dplyr::rename(CHR = seqnames, BP_hg19 = start, BP_hg38 = mcols.BP) %>%
  mutate(CHR = as.integer(gsub("chr", "", CHR)))

final_dementia <- dementia %>%
  inner_join(converted_data, by = c("CHR", "BP" = "BP_hg38")) %>%
  select(-BP) %>%
  dplyr::rename(BP = "BP_hg19")

final_dementia <- dplyr::rename(final_dementia, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

# Export formatted data
dementia_mtag <- final_dementia %>%
  dplyr::select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
dementia_ldsc <- final_dementia %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(dementia_mtag, "MTAG_input_data/MTAGDEMsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dementia_ldsc, "LDSC_input_data/LDSCDEMsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dementia_ldsc, "Final_processed_data/DEMsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
