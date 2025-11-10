# ============================
# Step 0. Install and load required packages
# ============================
install.packages("httr")
library(rtracklayer)
library(data.table)
library(dplyr)
library(httr)

setwd("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS")

# ============================
# Step 1. Process DD dataset (GCST90476922, build hg38)
# ============================

dd <- fread("Downloaded_original_data/GCST90476922.tsv.gz")
colnames(dd)

# Check missing base pair locations
View(dd[is.na(dd$base_pair_location), ])

# Compute beta and standard error from odds ratio and confidence interval
dd$beta <- log(dd$odds_ratio)
dd$standard_error <- (log(dd$ci_upper) - log(dd$ci_lower)) / (2 * 1.96)  # Assuming CI for OR

# Rename columns to standardized format
name_mapping <- c(
  rsid = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ",
  n = "N"
)
names(dd) <- name_mapping[names(dd)]

# Calculate minor allele frequency (MAF) and Z-score
dd$MAF <- pmin(dd$FRQ, 1 - dd$FRQ)
dd$Z <- dd$beta / dd$se
colnames(dd)

View(dd[is.na(dd$BP), ])

# Keep only required columns
dd <- dd %>% 
  dplyr::select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

# Format summary statistics and save as GRCh38
format_sumstats(
  path = dd,
  ref_genome = "GRCh38",
  dbSNP = 144,
  save_path = "Format_sumstats_data/DD38.tsv.gz"
)

# Reload formatted file
dd <- read.table("Format_sumstats_data/DD38.tsv.gz", header = TRUE, sep = "\t")

# ============================
# Step 1.1. LiftOver: hg38 to hg19
# ============================

chain <- import.chain("D:/hg38ToHg19.over.chain/hg38ToHg19.over.chain")

# Create GRanges object for conversion
gr_hg38 <- GRanges(
  seqnames = paste0("chr", dd$CHR),
  ranges = IRanges(start = dd$BP, end = dd$BP),
  strand = "*",
  mcols = dd
)

# Coordinate transformation
gr_hg19 <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(gr_hg19)

# Extract converted coordinates
converted_data <- as.data.frame(gr_hg19) %>%
  dplyr::rename(
    CHR = seqnames,
    BP_hg19 = start,
    BP_hg38 = mcols.BP
  ) %>%
  mutate(CHR = as.integer(gsub("chr", "", CHR)))

# Merge converted positions back to original data
final_dd <- dd %>%
  inner_join(converted_data, by = c("CHR", "BP" = "BP_hg38")) %>%
  select(-BP) %>%
  dplyr::rename(BP = "BP_hg19")

# Get unmapped SNPs
failed_dd <- dd %>%
  anti_join(converted_data, by = c("CHR", "BP" = "BP_hg38"))

colnames(final_dd)

# Rename for downstream tools (A1 = effect allele)
final_dd <- dplyr::rename(final_dd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

# Export MTAG and LDSC formatted files
dd_mtag <- final_dd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
dd_ldsc <- final_dd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(dd_mtag, "MTAG_input_data/MTAGDDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dd_ldsc, "LDSC_input_data/LDSCDDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dd_ldsc, "Final_processed_data/DDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ============================
# Step 2. Process MADD dataset (GCST90225526, build hg19)
# ============================

madd <- fread("Downloaded_original_data/GCST90225526_buildGRCh37.tsv")
colSums(is.na(madd))

name_mapping <- c(
  variant_id = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ"
)
names(madd) <- name_mapping[names(madd)]

madd$N <- 502480
madd$MAF <- pmin(madd$FRQ, 1 - madd$FRQ)
madd$Z <- madd$beta / madd$se

colnames(madd)

madd <- madd %>%
  dplyr::select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

# Format and save as GRCh37
format_sumstats(
  path = madd,
  ref_genome = "GRCh37",
  dbSNP = 144,
  save_path = "Format_sumstats_data/MADD37.tsv.gz"
)

madd <- read.table("Format_sumstats_data/MADD37.tsv.gz", header = TRUE, sep = "\t")
colSums(is.na(madd))

# Rename and prepare for export
madd <- dplyr::rename(madd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

madd_mtag <- madd %>%
  select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
madd_ldsc <- madd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(madd_mtag, "MTAG_input_data/MTAGMADDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(madd_ldsc, "LDSC_input_data/LDSCMADDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(madd_ldsc, "Final_processed_data/MADDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# ============================
# Step 3. Process MDD dataset (GCST90468123, build hg19)
# ============================

mdd <- fread("Downloaded_original_data/GCST90468123.tsv.gz")
colSums(is.na(mdd))
colnames(mdd)

name_mapping <- c(
  variant_id = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ",
  n = "N"
)
names(mdd) <- name_mapping[names(mdd)]

mdd$MAF <- pmin(mdd$FRQ, 1 - mdd$FRQ)
mdd$Z <- mdd$beta / mdd$se

mdd <- mdd %>%
  dplyr::select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

# Format and save as GRCh37
format_sumstats(
  path = mdd,
  ref_genome = "GRCh37",
  dbSNP = 144,
  save_path = "Format_sumstats_data/MDD37.tsv.gz"
)

mdd <- read.table("Format_sumstats_data/MDD37.tsv.gz", header = TRUE, sep = "\t")
colSums(is.na(mdd))

# Check SNP duplicates
has_duplicates <- anyDuplicated(mdd$SNP) > 0

# Rename and prepare for export
mdd <- dplyr::rename(mdd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

mdd_mtag <- mdd %>%
  select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
mdd_ldsc <- mdd %>%
  dplyr::select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(mdd_mtag, "MTAG_input_data/MTAGMDDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mdd_ldsc, "LDSC_input_dataLDSCMDDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mdd_ldsc, "Final_processed_data/MDDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
