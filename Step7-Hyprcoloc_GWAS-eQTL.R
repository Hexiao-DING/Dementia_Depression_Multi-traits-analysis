# Hyprcoloc pipeline that pairs MTAG GWAS results with GTEx v7 eQTL panels.

library(dplyr)
library(tidyr)
library(MungeSumstats)
library(data.table)
library(stringr)
library(duckdb)
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ------Fixed Step: Load EQTL data and process in batch-------------------------
load_qtl_data <- function(folder_path) {
  # List all files that match the pattern
  files <- list.files(folder_path, pattern = "\\.v7\\.signif_variant_gene_pairs\\.txt\\.gz$", full.names = TRUE)
  
  # Initialize list
  qtl_data_list <- list()
  
  # Loop and read files
  for (file in files) {
    # Extract tissue name
    filename <- basename(file)
    tissue <- str_extract(filename, "^[^.]+")
    
    # Read data
    data <- fread(file)
    
    # Store in list
    qtl_data_list[[tissue]] <- data
  }
  
  # Return list
  return(qtl_data_list)
}

folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/eQTL/GTEx_Analysis_v7_eQTL"
qtl_data_list <- load_qtl_data(folder_path)
length(qtl_data_list)

# ------Process each EQTL dataset-------------------------
process_eqtl_data_list <- function(data_list) {
  
  # Process each dataset
  processed_list <- lapply(names(data_list), function(name) {
    df <- data_list[[name]]
    
    # Check for required columns
    required_cols <- c("variant_id", "maf", "slope", "slope_se", "pval_nominal", "gene_id")
    missing_cols <- setdiff(required_cols, colnames(df))
    
    if (length(missing_cols) > 0) {
      stop("Dataset '", name, "' is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # 1. Split variant_id column
    df_processed <- df %>%
      tidyr::separate(variant_id, 
                      into = c("CHR", "BP", "A1", "A2", "Build"),
                      sep = "_", remove = FALSE) %>%
      dplyr::rename(EF = maf, 
                    BETA = slope, 
                    std = slope_se, 
                    P = pval_nominal, 
                    GN = gene_id)
    
    # 2. Format summary statistics
    df_formatted <- MungeSumstats::format_sumstats(
      df_processed[, c("CHR", "BP", "A1", "A2", "BETA", "P", "std", "EF", "GN")],
      ref_genome = "GRCh37",
      nThread = 22,
      dbSNP = 144,
      check_dups = FALSE,
      return_data = TRUE
    )
    
    # 3. Calculate MAF and rename columns
    df_final <- df_formatted %>%
      dplyr::mutate(MAF = pmin(EF, 1 - EF)) %>%
      dplyr::rename(beta = BETA, 
                    EAF = EF, 
                    se = STD,
                    Gene = GN) %>%
      as.data.frame()
    
    return(df_final)
  })
  
  # Preserve original list names
  names(processed_list) <- names(data_list)
  
  return(processed_list)
}

qtl_list <- process_eqtl_data_list(qtl_data_list)
qtl_list
length(unique(qtl_list[['Whole_Blood']]$SNP))

#------------------Fixed Functions-------------------------
create_gwas_list <- function(trait1, trait2) {
  folder_path <- file.path("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results", 
                           paste0(trait1, "_", trait2))
  
  # Check if folder exists
  if (!dir.exists(folder_path)) {
    stop("Folder does not exist: ", folder_path)
  }
  
  # Construct file paths
  file1 <- file.path(folder_path, "mtag_results_trait_1.txt")
  file2 <- file.path(folder_path, "mtag_results_trait_2.txt")
  
  # Read data
  trait1_data <- data.table::fread(file1)
  trait2_data <- data.table::fread(file2)
  
  # Create named list
  gwas_list <- list(
    trait1 = trait1_data,
    trait2 = trait2_data
  )
  
  # Rename list elements according to input traits
  names(gwas_list) <- c(trait1, trait2)
  
  return(gwas_list)
}

get_target_snps <- function(file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_Trait-Trait.csv", 
                            filter_col = "Traits", filter_trait, pp_col = "PP", snp_col = "Candidate.Causal.SNP") {
  df <- read.csv(file, stringsAsFactors = FALSE)
  filtered_df <- df[
    df[[filter_col]] == filter_trait & 
      df[[pp_col]] >= 0.8,
  ]
  unique_snps <- unique(filtered_df[[snp_col]])
  return(unique_snps)
}

run_coloc_analysis <- function(gwas_files, qtl_files, target_snps, 
                               window_size = 200000, gwas_threshold = 5e-4, qtl_threshold = 5e-4) {
  
  # Initialize result list
  results_list <- list()
  result_count <- 1
  
  # Loop through each target SNP
  for (target_snp in target_snps) {
    message("\n===== Processing target SNP: ", target_snp, " =====")
    
    # Loop through each GWAS dataset
    for (gwas_name in names(gwas_files)) {
      GWASdata <- gwas_files[[gwas_name]]
      
      # Loop through each QTL dataset
      for (qtl_name in names(qtl_files)) {
        QTLdata <- qtl_files[[qtl_name]]
        QTLdata$SNP_Gene <- paste(QTLdata$SNP, QTLdata$Gene, sep = "_")
        
        # Get target SNP info in GWAS
        snp_info1 <- GWASdata[GWASdata$SNP == target_snp, c("CHR", "BP")]
        chr <- snp_info1$CHR
        bp_center <- snp_info1$BP
        
        GWASdata_f <- GWASdata[GWASdata$CHR == chr & 
                                 GWASdata$BP >= (bp_center - window_size) & 
                                 GWASdata$BP <= (bp_center + window_size), ]
        
        # Get target SNP info in QTL
        snp_info2 <- QTLdata[QTLdata$SNP == target_snp, c("CHR", "BP")]
        
        # Check if target SNP exists in QTL
        if(nrow(snp_info2) > 0) {
          # Use QTL position
          QTLdata_f <- QTLdata[QTLdata$CHR == snp_info2$CHR & 
                                 QTLdata$BP >= (snp_info2$BP - window_size) & 
                                 QTLdata$BP <= (snp_info2$BP + window_size), ]
        } else {
          # Use GWAS position
          message("Note: ", target_snp, " not found in QTL, using GWAS position")
          QTLdata_f <- QTLdata[QTLdata$CHR == chr & 
                                 QTLdata$BP >= (bp_center - window_size) & 
                                 QTLdata$BP <= (bp_center + window_size), ]
        }
        
        # Get union and intersect of SNPs
        union_snps <- unique(union(GWASdata_f$SNP, QTLdata_f$SNP))
        common_snps <- intersect(GWASdata$SNP[GWASdata$SNP %in% union_snps],
                                 QTLdata$SNP[QTLdata$SNP %in% union_snps])
        
        if (length(common_snps) == 0) {
          message("Skipping ", gwas_name, " - ", qtl_name, ": No common SNPs")
          next
        }
        
        # Create full dataframe
        full_data <- data.frame(SNP = common_snps, stringsAsFactors = FALSE)
        
        # Merge GWAS data
        gwas_merge <- merge(full_data, 
                            GWASdata[, c("SNP", "mtag_beta", "mtag_se", "mtag_pval")], 
                            by = "SNP", all.x = TRUE)
        colnames(gwas_merge)[(ncol(gwas_merge)-2):ncol(gwas_merge)] <- c("beta1", "se1", "pval1")
        
        # Merge QTL data
        full_data <- merge(gwas_merge, 
                           QTLdata[, c("SNP", "beta", "se", "P", "SNP_Gene")], 
                           by = "SNP", all.x = TRUE)
        colnames(full_data)[(ncol(full_data)-3):ncol(full_data)] <- c("beta2", "se2", "pval2", "SNP_Gene")
        
        # Skip if no significant SNPs
        if (min(gwas_merge$pval1, na.rm = TRUE) > gwas_threshold && 
            min(full_data$pval2, na.rm = TRUE) > qtl_threshold) {
          message("Skipping ", gwas_name, " - ", qtl_name, ": No significant SNPs in both datasets")
          next
        }
        
        # Create beta and SE matrices
        beta_matrix <- as.matrix(full_data[, c("beta1", "beta2")])
        se_matrix <- as.matrix(full_data[, c("se1", "se2")])
        
        # Set row and column names
        traits <- c('GWAS', 'QTL')
        rownames(beta_matrix) <- full_data$SNP_Gene
        colnames(beta_matrix) <- traits
        rownames(se_matrix) <- full_data$SNP_Gene
        colnames(se_matrix) <- traits
        
        # Perform coloc analysis
        hyprcoloc_result <- tryCatch({
          hyprcoloc(
            beta_matrix,
            se_matrix,
            trait.names = traits,
            binary.outcomes = c(1,0),
            snp.id = rownames(beta_matrix),
            prior.1 = 1e-4,
            prior.c = 0.02
          )
        }, error = function(e) {
          message("Coloc analysis failed: ", e$message)
          return(NULL)
        })
        
        if (!is.null(hyprcoloc_result)) {
          # Create result dataframe
          result_df <- data.frame(
            GWAS = gwas_name,
            QTL = qtl_name,
            Target_SNP = target_snp,
            PP = hyprcoloc_result$results$posterior_prob,
            PP_explained = hyprcoloc_result$results$posterior_explained_by_snp,
            regional_prob = hyprcoloc_result$results$regional_prob,
            Candidate_SNPs = hyprcoloc_result$results$candidate_snp,
            stringsAsFactors = FALSE
          )
          
          # Add to result list
          results_list[[result_count]] <- result_df
          result_count <- result_count + 1
          
          message("Completed: ", target_snp, " | ", gwas_name, " - ", qtl_name)
        }
      }
    }
  }
  
  # Combine all results
  if (length(results_list) > 0) {
    final_results <- do.call(rbind, results_list)
    return(final_results)
  } else {
    message("No successful analysis completed")
    return(NULL)
  }
}

#------------------Clean and Expand Coloc Results-------------------------
clean_and_expand_coloc <- function(coloc_results, trait_map) {
  
  coloc_clean <- coloc_results[!is.na(coloc_results$PP) & coloc_results$PP >= 0.5, ]
  
  if ("Candidate_SNPs" %in% names(coloc_clean)) {
    library(tidyr)
    coloc_clean <- coloc_clean %>%
      separate(
        Candidate_SNPs,
        into = c("Causal_SNP", "Gene_Ensembl"),
        sep = "_",
        remove = FALSE,
        extra = "merge"
      )
  }
  
  coloc_clean$Traits <- trait_map
  
  return(coloc_clean)
}
#------------------DD_DEM---------------------
gwas_files <- create_gwas_list("DD", "DEM")
target_snps <- get_target_snps(
  filter_trait = "Dementia, Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'DD_DEM')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_DEM_DD_QTL.csv", row.names = FALSE)

#------------------MADD_DEM---------------------
gwas_files <- create_gwas_list("MADD", "DEM")
target_snps <- get_target_snps(
  filter_trait = "Dementia, Mixed Anxiety and Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MADD_DEM')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_DEM_MADD_QTL.csv", row.names = FALSE)

#------------------DD_CP---------------------
gwas_files <- create_gwas_list("DD", "CP")
target_snps <- get_target_snps(
  filter_trait = "CP, Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'DD_CP')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_DD_QTL.csv", row.names = FALSE)

#------------------MDD_CP---------------------
gwas_files <- create_gwas_list("MDD", "CP")
target_snps <- get_target_snps(
  filter_trait = "CP, Major Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MDD_CP')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_MDD_QTL.csv", row.names = FALSE)

#------------------MADD_CP---------------------
gwas_files <- create_gwas_list("MADD", "CP")
target_snps <- get_target_snps(
  filter_trait = "CP, Mixed Anxiety and Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MADD_CP')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_CP_MADD_QTL.csv", row.names = FALSE)

#------------------MADD_UD---------------------
gwas_files <- create_gwas_list("MADD", "UD")
target_snps <- get_target_snps(
  filter_trait = "UD, Mixed Anxiety and Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MADD_UD')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_UD_MADD_QTL.csv", row.names = FALSE)

#------------------DD_VD---------------------
gwas_files <- create_gwas_list("DD", "VD")
target_snps <- get_target_snps(
  filter_trait = "VD, Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)# No colocalization results for this trait pair
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'DD_VD')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_VD_DD_QTL.csv", row.names = FALSE)

#------------------MADD_VD---------------------
gwas_files <- create_gwas_list("MADD", "VD")
target_snps <- get_target_snps(
  filter_trait = "VD, Mixed Anxiety and Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)# No colocalization results for this trait pair
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MADD_VD')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_VD_MADD_QTL.csv", row.names = FALSE)

#------------------MADD_LBD---------------------
gwas_files <- create_gwas_list("MADD", "LBD")
target_snps <- get_target_snps(
  filter_trait = "LBD, Mixed Anxiety and Depression")
target_snps
coloc_results <- run_coloc_analysis(
  gwas_files = gwas_files,
  qtl_files = qtl_list,
  target_snps = target_snps
)# No colocalization results for this trait pair
coloc_results_clean <- clean_and_expand_coloc(coloc_results, trait_map = 'MADD_LBD')
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/Raw_coloc_results/coloc_LBD_MADD_QTL.csv", row.names = FALSE)
