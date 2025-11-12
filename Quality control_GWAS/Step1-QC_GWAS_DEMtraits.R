# ===============================================================
# Step 0. Install and load required packages
# ===============================================================

if(!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("SNPlocs.Hsapiens.dbSNP144.GRCh37",
                       "SNPlocs.Hsapiens.dbSNP144.GRCh38",
                       "BSgenome.Hsapiens.NCBI.GRCh38",
                       "BSgenome.Hsapiens.1000genomes.hs37d5",
                       "rtracklayer"))

install.packages(c("data.table", "dplyr", "MungeSumstats"))

library(data.table)
library(dplyr)
library(MungeSumstats)
library(rtracklayer)
setwd("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS")

# ===============================================================
# Step 1. Alzheimer's disease (GCST90012877, build hg19)
# ===============================================================

ad <- fread("Downloaded_original_data/GCST90012877_buildGRCh37.tsv.gz")

# Rename columns
name_mapping <- c(
  variant_id = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  GWAS_BETA = "beta",
  GWAS_SE = "se",
  GWAS_P = "pval",
  effect_allele_frequency = "FRQ"
)
names(ad) <- name_mapping[names(ad)]

# Add derived variables
ad$N <- 472868
ad$MAF <- pmin(ad$FRQ, 1 - ad$FRQ)
ad$Z <- ad$beta / ad$se

# Keep only required columns
ad <- ad %>%
  dplyr::select(SNP, CHR, BP, effect_allele, other_allele,
                beta, se, Z, pval, N, FRQ, MAF)

# Format and save GRCh37
format_sumstats(
  path = ad,
  ref_genome = "GRCh37",
  dbSNP = 144,
  save_path = "Format_sumstats_data/AD37.tsv.gz"
)

ad <- read.table("Format_sumstats_data/AD37.tsv.gz", header = TRUE, sep = "\t")

# Prepare for MTAG / LDSC
ad <- dplyr::rename(ad, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

ad_mtag <- ad %>% select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
ad_ldsc <- ad %>% select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(ad_mtag, "MTAG_input_data/MTAGADsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ad_ldsc, "LDSC_input_data/LDSCADsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ad_ldsc, "Format_sumstats_data/ADsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ===============================================================
# Step 2. Cognitive Performance (GCST006572, build hg19)
# ===============================================================

cp <- fread("Downloaded_original_data/30038396-GCST006572-EFO_0008354-build37.f.tsv.gz")
colSums(is.na(cp)) 
colnames(cp)
name_mapping <- c(
  variant_id = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ"
)
names(cp) <- name_mapping[names(cp)]

cp$N <- 257841
cp$MAF <- pmin(cp$FRQ, 1 - cp$FRQ)
cp$Z <- cp$beta / cp$se

cp <- cp %>%
  select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

format_sumstats(
  path = cp,
  ref_genome = "GRCh37",
  dbSNP = 144,
  save_path = "Format_sumstats_data/CP37.tsv.gz"
)

cp <- read.table("Format_sumstats_data/CP37.tsv.gz", header = TRUE, sep = "\t")

cp <- dplyr::rename(cp, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

cp_mtag <- cp %>% select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
cp_ldsc <- cp %>% select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(cp_mtag, "MTAG_input_data/MTAGCPsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cp_ldsc, "LDSC_input_data/LDSCCPsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cp_ldsc, "Final_processed_data/CPsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ===============================================================
# Step 3. Frontotemporal Dementia (GCST90558311, build hg19)
# ===============================================================

fd <- fread("Downloaded_original_data/GCST90558311.tsv")

name_mapping <- c(
  rsid = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ",
  Z_STAT = "Z",
  n = "N"
)
names(fd) <- name_mapping[names(fd)]

fd$MAF <- pmin(fd$FRQ, 1 - fd$FRQ)

fd <- fd %>%
  select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

format_sumstats(
  path = fd,
  ref_genome = "GRCh37",
  dbSNP = 144,
  save_path = "Format_sumstats_data/FD37.tsv.gz"
)

fd <- read.table("Format_sumstats_data/FD37.tsv.gz", header = TRUE, sep = "\t")

fd <- dplyr::rename(fd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

fd_mtag <- fd %>% select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
fd_ldsc <- fd %>% select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(fd_mtag, "MTAG_input_data/MTAGFDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fd_ldsc, "LDSC_input_data/LDSCFDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fd_ldsc, "Final_processed_data/FDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ===============================================================
# Step 4. Lewy Body Dementia (GCST90001390, build hg38 to hg19)
# ===============================================================

lbd <- fread("Downloaded_original_data/GCST90001390_buildGRCh38.tsv.gz")

name_mapping <- c(
  variant_id = "SNP",
  chromosome = "CHR",
  base_pair_location = "BP",
  beta = "beta",
  standard_error = "se",
  p_value = "pval",
  effect_allele_frequency = "FRQ"
)
names(lbd) <- name_mapping[names(lbd)]

lbd$N <- 6618
lbd$MAF <- pmin(lbd$FRQ, 1 - lbd$FRQ)
lbd$Z <- lbd$beta / lbd$se

lbd <- lbd %>%
  select(SNP, CHR, BP, effect_allele, other_allele, beta, se, Z, pval, N, FRQ, MAF)

# Format and convert (hg38 to hg19)
format_sumstats(
  path = lbd,
  ref_genome = "GRCh38",
  dbSNP = 144,
  save_path = "Format_sumstats_data/LBD38.tsv.gz"
)
lbd <- read.table("Format_sumstats_data/LBD38.tsv.gz", header = TRUE, sep = "\t")

# LiftOver
chain <- import.chain("D:/hg38ToHg19.over.chain/hg38ToHg19.over.chain")
gr_hg38 <- GRanges(seqnames = paste0("chr", lbd$CHR),
                   ranges = IRanges(start = lbd$BP, end = lbd$BP),
                   strand = "*", mcols = lbd)
gr_hg19 <- unlist(liftOver(gr_hg38, chain))

converted_data <- as.data.frame(gr_hg19) %>%
  rename(CHR = seqnames, BP_hg19 = start, BP_hg38 = mcols.BP) %>%
  mutate(CHR = as.integer(gsub("chr", "", CHR)))

final_lbd <- lbd %>%
  inner_join(converted_data, by = c("CHR", "BP" = "BP_hg38")) %>%
  select(-BP) %>%
  rename(BP = "BP_hg19")

# Prepare for MTAG/LDSC
final_lbd <- rename(final_lbd, A2 = A1, A1 = A2, beta = BETA, se = SE, pval = P)

lbd_mtag <- final_lbd %>% select(SNP, CHR, BP, A2, A1, pval, N, FRQ, Z)
lbd_ldsc <- final_lbd %>% select(SNP, CHR, BP, A2, A1, beta, se, pval, N, FRQ, MAF)

write.table(lbd_mtag, "MTAG_input_data/MTAGLBDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lbd_ldsc, "LDSC_input_data/LDSCLBDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lbd_ldsc, "Final_processed_data/LBDsummary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

