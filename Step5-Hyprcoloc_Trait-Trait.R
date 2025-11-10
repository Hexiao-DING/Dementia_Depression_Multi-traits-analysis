library(readr)
library(dplyr)
library(data.table)
library(devtools)
library(RcppEigen)
library(hyprcoloc)
?hyprcoloc

# Function to get top SNPs from a single file
get_target_snps_single <- function(file, snp_col = "SNP") {
  df <- read.csv(file, stringsAsFactors = FALSE)
  top_snps <- df[df$Signal.type == "Top", ]
  unique_snps <- unique(top_snps[[snp_col]])
}

# Function to get top SNPs from two files
get_target_snps_double <- function(file1, file2, snp_col = "SNP") {
  df1 <- read.csv(file1, stringsAsFactors = FALSE)
  df2 <- read.csv(file2, stringsAsFactors = FALSE)
  
  if (!identical(names(df1), names(df2))) {
    stop("Column names of the two datasets are not identical, please check files:", file1, " and ", file2)
  }
  
  combined_df <- rbind(df1, df2)
  top_snps <- combined_df[combined_df$Signal.type == "Top", ]
  unique_snps <- unique(top_snps[[snp_col]])
}

# Function to run hyprcoloc analysis
run_hyprcoloc_analysis <- function(folder, origin_gwas1, origin_gwas2, target_snps, traits) {
  
  # Read GWAS result files
  gwas1 <- fread(file.path(folder, "mtag_results_trait_1.txt"))
  gwas2 <- fread(file.path(folder, "mtag_results_trait_2.txt"))
  
  # Initialize a list to store significant results
  all_significant_results <- list()
  
  # Loop through each target SNP
  for (target_snp in target_snps) {
    # Get SNP information
    snp_info1 <- gwas1[SNP == target_snp, .(CHR, BP)]
    snp_info2 <- gwas2[SNP == target_snp, .(CHR, BP)]
    
    # Extract regional SNPs (+/- 200kb)
    gwas1_f <- gwas1[CHR == snp_info1$CHR & 
                       BP >= (snp_info1$BP - 200000) & 
                       BP <= (snp_info1$BP + 200000)]
    
    gwas2_f <- gwas2[CHR == snp_info2$CHR & 
                       BP >= (snp_info2$BP - 200000) & 
                       BP <= (snp_info2$BP + 200000)]
    
    # Union of SNPs
    union_snps <- unique(union(gwas1_f$SNP, gwas2_f$SNP))
    
    # Create full data frame
    full_data <- data.frame(SNP = union_snps, stringsAsFactors = FALSE)
    
    # Merge GWAS summary statistics
    origin_gwas_check1 <- merge(full_data, 
                                origin_gwas1[, .(SNP, pval)], 
                                by = "SNP", all.x = TRUE)
    origin_gwas_check2 <- merge(full_data, 
                                origin_gwas2[, .(SNP, pval)], 
                                by = "SNP", all.x = TRUE)
    
    # Skip if no SNP reaches significance (P < 5e-4)
    if (min(origin_gwas_check1$pval, na.rm = TRUE) > 5e-4 && 
        min(origin_gwas_check2$pval, na.rm = TRUE) > 5e-4) {
      message("Skipping ", target_snp, ": No significant SNP (P < 5e-4) in both datasets")
      next
    }
    
    # Merge MTAG GWAS data
    gwas1_merge <- merge(full_data, 
                         gwas1[, .(SNP, beta1 = mtag_beta, se1 = mtag_se)], 
                         by = "SNP", all.x = TRUE)
    
    full_data <- merge(gwas1_merge, 
                       gwas2[, .(SNP, beta2 = mtag_beta, se2 = mtag_se)], 
                       by = "SNP", all.x = TRUE)
    
    # Create beta and se matrices
    beta_matrix <- as.matrix(full_data[, c("beta1", "beta2")])
    se_matrix <- as.matrix(full_data[, c("se1", "se2")])
    
    rownames(beta_matrix) <- full_data$SNP
    rownames(se_matrix) <- full_data$SNP
    colnames(beta_matrix) <- traits
    colnames(se_matrix) <- traits
    
    # Run hyprcoloc
    hyprcoloc_result <- hyprcoloc(
      beta_matrix,
      se_matrix,
      trait.names = traits,
      snp.id = rownames(beta_matrix),
      binary.outcomes = c(1,1),
      prior.1 = 1e-4,
      prior.c = 0.02
    )
    
    # Extract results
    result_df <- as.data.frame(hyprcoloc_result$results)
    
    # Filter posterior probability > 0.5
    if (any(result_df$posterior_prob >= 0.5 & !is.na(result_df$posterior_prob))) {
      sig_idx <- which(result_df$posterior_prob >= 0.5)
      
      sig_df <- data.frame(
        target_snp = target_snp,
        traits = result_df$traits[sig_idx],
        candidate_snps = result_df$candidate_snp[sig_idx],
        posterior_probs = result_df$posterior_prob[sig_idx],
        PP_explained = result_df$posterior_explained_by_snp[sig_idx],
        regional_prob = result_df$regional_prob[sig_idx],
        stringsAsFactors = FALSE
      )
      
      all_significant_results[[target_snp]] <- sig_df
    }
  }
  
  # Combine all significant results
  if (length(all_significant_results) > 0) {
    final_df <- do.call(rbind, all_significant_results)
    final_df <- dplyr::rename(final_df, 'SNP' = 'target_snp')
  } else {
    final_df <- data.frame()
  }
  
  return(final_df)
}

# --------------------- DD_DEM ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_DEM"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DEM_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DEM_DD.csv")
traits <- c('DD', 'DEM')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Dementia, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_DEM_DD.csv", row.names = FALSE)

# --------------------- MDD_DEM ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_DEM"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DEM_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DEM_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_DEM.csv")
traits <- c('MDD', 'DEM')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Dementia, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_DEM_MDD.csv", row.names = FALSE)

# --------------------- MADD_DEM ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_DEM"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DEM_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DEM_MADD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MADD_DEM.csv")
traits <- c('MADD', 'DEM')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Dementia, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_DEM_MADD.csv", row.names = FALSE)

# --------------------- DD_AD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_AD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/AD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_AD_DD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DD_AD.csv")
traits <- c('DD', 'AD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'AD, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_AD_DD.csv", row.names = FALSE)

# --------------------- MDD_AD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_AD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/AD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_AD_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_AD.csv")
traits <- c('MDD', 'AD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'AD, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_AD_MDD.csv", row.names = FALSE)

# --------------------- MADD_AD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_AD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/AD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_AD_MADD.csv")
traits <- c('MADD', 'AD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'AD, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_AD_MADD.csv", row.names = FALSE)

# --------------------- DD_CP ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_CP"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/CP_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_CP_DD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DD_CP.csv")
traits <- c('DD', 'CP')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Cognitive Performance, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_DD.csv", row.names = FALSE)

# --------------------- MDD_CP ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_CP"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/CP_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_CP_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_CP.csv")
traits <- c('MDD', 'CP')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Cognitive Performance, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_MDD.csv", row.names = FALSE)

# --------------------- MADD_CP ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_CP"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/CP_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_CP_MADD.csv")
traits <- c('MADD', 'CP')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Cognitive Performance, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_MADD.csv", row.names = FALSE)

# --------------------- DD_FD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_FD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/FD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_FD_DD.csv")
target_snps
traits <- c('DD', 'FD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'FD, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_FD_DD.csv", row.names = FALSE)

# --------------------- MDD_FD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_FD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/FD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_FD_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_FD.csv")
traits <- c('MDD', 'FD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'FD, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_FD_MDD.csv", row.names = FALSE)

# --------------------- DD_UD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_UD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/UD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_UD_DD.csv")
target_snps
traits <- c('DD', 'UD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'UD, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_UD_DD.csv", row.names = FALSE)

# --------------------- MDD_UD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_UD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/UD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_UD_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_UD.csv")
target_snps
traits <- c('MDD', 'UD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'UD, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_UD_MDD.csv", row.names = FALSE)

# --------------------- MADD_UD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_UD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/UD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_UD_MADD.csv")
target_snps
traits <- c('MADD', 'UD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
# Print significant colocalized SNPs
if (nrow(results) > 0) {
  print("Significant colocalized SNPs found:")
  for (i in 1:nrow(results)) {
    cat("\nTarget SNP:", results$SNP[i])
    cat("\nCandidate SNP:", results$candidate_snps[i])
    cat("\nPosterior Probability:", results$posterior_probs[i], "\n")
  }
} else {
  print("No colocalization results with posterior probability > 0.5")
}
results$Traits <- 'UD, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_UD_MADD.csv", row.names = FALSE)

# --------------------- DD_VD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_VD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/VD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_VD_DD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DD_VD.csv")
target_snps
traits <- c('DD', 'VD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
# Print significant colocalized SNPs
if (nrow(results) > 0) {
  print("Significant colocalized SNPs found:")
  for (i in 1:nrow(results)) {
    cat("\nTarget SNP:", results$SNP[i])
    cat("\nCandidate SNP:", results$candidate_snps[i])
    cat("\nPosterior Probability:", results$posterior_probs[i], "\n")
  }
} else {
  print("No colocalization results with posterior probability > 0.5")
}
results$Traits <- 'VD, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_VD_DD.csv", row.names = FALSE)

# --------------------- MDD_VD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_VD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/VD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_VD_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_VD.csv")
target_snps
traits <- c('MDD', 'VD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
# Print significant colocalized SNPs
if (nrow(results) > 0) {
  print("Significant colocalized SNPs found:")
  for (i in 1:nrow(results)) {
    cat("\nTarget SNP:", results$SNP[i])
    cat("\nCandidate SNP:", results$candidate_snps[i])
    cat("\nPosterior Probability:", results$posterior_probs[i], "\n")
  }
} else {
  print("No colocalization results with posterior probability > 0.5")
}
results$Traits <- 'VD, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_VD_MDD.csv", row.names = FALSE)

# --------------------- MADD_VD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_VD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/VD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_VD_MADD.csv")
target_snps
traits <- c('MADD', 'VD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'VD, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_VD_MADD.csv", row.names = FALSE)

# --------------------- DD_LBD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/DD_LBD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/LBD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_LBD_DD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_DD_LBD.csv")
traits <- c('DD', 'LBD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Lewy Body Dementia, Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_LBD_DD.csv", row.names = FALSE)

# --------------------- MDD_LBD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MDD_LBD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/LBD_summary.txt")
target_snps <- get_target_snps_double(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_LBD_MDD.csv",
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_MDD_LBD.csv")
traits <- c('MDD', 'LBD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Lewy Body Dementia, Major Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_LBD_MDD.csv", row.names = FALSE)

# --------------------- MADD_LBD ---------------------
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results/MADD_LBD"
origin_gwas1 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADD_summary.txt")
origin_gwas2 <- fread("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/LBD_summary.txt")
target_snps <- get_target_snps_single(
  "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Results/FUMA/FUMA_results_LBD_MADD.csv")
traits <- c('MADD', 'LBD')
results <- run_hyprcoloc_analysis(folder_path, origin_gwas1, origin_gwas2, target_snps, traits)
results$Traits <- 'Lewy Body Dementia, Mixed Anxiety and Depression'
write.csv(results, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_LBD_MADD.csv", row.names = FALSE)

