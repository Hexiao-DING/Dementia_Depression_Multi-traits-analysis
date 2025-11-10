# Hyprcoloc workflow that integrates MTAG GWAS outputs with MetaBrain eQTL datasets.

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(duckdb)
library(DBI)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(hyprcoloc)
# ------ Fixed Steps: Import eQTL Data and Batch Processing -------------------------
load_qtl_data <- function(folder_path) {
  # List all matching files
  files <- list.files(folder_path, pattern = ".*chr[0-9]+\\.txt\\.gz$", full.names = TRUE)
  
  # Initialize list
  qtl_data_list <- list()
  
  # Iterate and read files
  for (file in files) {
    
    # Extract tissue name and related information
    filename <- basename(file)
    # Use regex to extract basalganglia-EUR-30PCs-chr19 part
    tissue <- str_extract(filename, "(?<=-)[a-zA-Z]+-[A-Z]+-[0-9]+PCs-chr[0-9]+(?=\\.txt\\.gz)")
    
    # Read data
    data <- fread(file, header = TRUE, sep="\t", quote="")
    
    # Store in list
    qtl_data_list[[tissue]] <- data
  }
  
  # Return list
  return(qtl_data_list)
}

# Cortex-EUR/Hippocampus-EUR/Spinalcord-EUR
folder_path <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/eQTL/Cortex-EUR"
qtl_data_list <- load_qtl_data(folder_path)

length(qtl_data_list)

process_eqtl_data_list <- function(data_list) {
  purrr::imap(data_list, function(df, name) {
    # First add non-effect allele column
    df_processed <- df %>%
      dplyr::mutate(
        A2 = purrr::map2_chr(SNPAlleles, SNPEffectAllele, ~ {
          alleles <- unlist(strsplit(.x, "/"))
          setdiff(alleles, .y)
        }),
        MAF = pmin(SNPEffectAlleleFreq, 1 - SNPEffectAlleleFreq)
      ) %>%
      # Then separate SNP column (outside mutate)
      tidyr::separate(SNP, into = c("CHR", "BP", "SNP", "A1_A2"), sep = ":", remove = FALSE) %>%
      # Finally rename other columns
      dplyr::rename(
        EAF = SNPEffectAlleleFreq,
        beta = MetaBeta,
        se = MetaSE,
        pval = MetaP,
        N = MetaPN,
        A1 = SNPEffectAllele
      )
    
    # Format and return
    MungeSumstats::format_sumstats(
      dplyr::select(df_processed, CHR, BP, A1, A2, beta, se, pval, MAF, Gene),
      ref_genome = "GRCh37",
      nThread = 22,
      dbSNP = 144,
      check_dups = FALSE,
      return_data = TRUE
    )
  })
}

qtl_list <- process_eqtl_data_list(qtl_data_list)

qtl_list

#------------------Fixed Functions------------------
create_gwas_list <- function(trait1, trait2) {
  folder_path <- file.path("D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results", 
                           paste0(trait1, "_", trait2))
  
  # Check if folder exists
  if (!dir.exists(folder_path)) {
    stop("Folder does not exist: ", folder_path)
  }
  
  # Build file paths
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
  
  # Use input trait names as list element names
  names(gwas_list) <- c(trait1, trait2)
  
  return(gwas_list)
}

get_target_snps <- function(file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/coloc_Trait-Trait.csv", 
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
  
  # Initialize results list
  results_list <- list()
  result_count <- 1
  
  # Loop through each target SNP
  for (target_snp in target_snps) {
    message("\n===== Processing target SNP: ", target_snp, " =====")
    
    # Loop through each GWAS file
    for (gwas_name in names(gwas_files)) {
      # Read GWAS data
      GWASdata <- gwas_files[[gwas_name]]
      
      # Loop through each QTL dataset
      for (qtl_name in names(qtl_files)) {
        QTLdata <- qtl_files[[qtl_name]]
        QTLdata$SNP_Gene <- paste(QTLdata$SNP, QTLdata$GENE, sep = "_")
        
        # Get current target SNP information in GWAS
        snp_info1 <- GWASdata[GWASdata$SNP == target_snp, c("CHR", "BP")]
        chr <- snp_info1$CHR
        bp_center <- snp_info1$BP
        
        GWASdata_f <- GWASdata[GWASdata$CHR == chr & 
                                 GWASdata$BP >= (bp_center - 200000) & 
                                 GWASdata$BP <= (bp_center + 200000), ]
        
        # Try to get target SNP information in QTL
        snp_info2 <- QTLdata[QTLdata$SNP == target_snp, c("CHR", "BP")]
        
        # Check if target SNP exists in QTL
        if(nrow(snp_info2) > 0) {
          # If exists, use QTL position information to extract region
          QTLdata_f <- QTLdata[QTLdata$CHR == snp_info2$CHR & 
                                 QTLdata$BP >= (snp_info2$BP - 200000) & 
                                 QTLdata$BP <= (snp_info2$BP + 200000), ]
        } else {
          # If not exists, use GWAS position information to extract QTL region
          message("Note: ", target_snp, " does not exist in QTL data, using GWAS position information")
          QTLdata_f <- QTLdata[QTLdata$CHR == chr & 
                                 QTLdata$BP >= (bp_center - 200000) & 
                                 QTLdata$BP <= (bp_center + 200000), ]
        }
        
        # Get union SNPs
        union_snps <- unique(union(GWASdata_f$SNP, QTLdata_f$SNP))
        common_snps <- intersect(GWASdata$SNP[GWASdata$SNP %in% union_snps],
                                 QTLdata$SNP[QTLdata$SNP %in% union_snps])
        
        if (length(common_snps) == 0) {
          message("Skip ", gwas_name, " - ", qtl_name, 
                  ": No common SNPs")
          next
        }
        
        # Create complete data frame
        full_data <- data.frame(SNP = common_snps, stringsAsFactors = FALSE)
        
        # Merge GWAS data
        gwas_merge <- merge(full_data, 
                            GWASdata[, c("SNP", "mtag_beta", "mtag_se", "mtag_pval")], 
                            by = "SNP", all.x = TRUE)
        colnames(gwas_merge)[(ncol(gwas_merge)-2):ncol(gwas_merge)] <- c("beta1", "se1", "pval1")
        
        full_data <- merge(gwas_merge, 
                           QTLdata[, c("SNP", "BETA", "SE", "P", "SNP_Gene")], 
                           by = "SNP", all.x = TRUE)
        colnames(full_data)[(ncol(full_data)-3):ncol(full_data)] <- c("beta2", "se2", "pval2", "SNP_Gene")
        
        # Check if there are significant SNPs
        if (min(gwas_merge$pval1, na.rm = TRUE) > gwas_threshold && 
            min(full_data$pval2, na.rm = TRUE) > qtl_threshold) {
          message("Skip ", gwas_name, " - ", qtl_name, 
                  ": No significant SNPs in both datasets")
          next
        }
        
        # Create matrices
        beta_matrix <- as.matrix(full_data[, c("beta1", "beta2")])
        se_matrix <- as.matrix(full_data[, c("se1", "se2")])
        
        # Set row and column names
        traits <- c('GWAS', 'QTL')
        rownames(beta_matrix) <- full_data$SNP_Gene  # Use SNP as row names
        colnames(beta_matrix) <- traits
        rownames(se_matrix) <- full_data$SNP_Gene
        colnames(se_matrix) <- traits
        
        # Perform colocalization analysis
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
          message("Colocalization analysis failed: ", e$message)
          return(NULL)
        })
        
        if (!is.null(hyprcoloc_result)) {
          # Create result data frame
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
          
          # Add to results list
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
    message("No successfully completed analyses")
    return(NULL)
  }
}

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_results/coloc_DEM_DD_basalganglia.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_results/coloc_DEM_MADD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/coloc_CP_DD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/coloc_CP_MDD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/cortex_CP_MADD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/coloc_UD_MADD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/coloc_VD_DD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_result/coloc_VD_MADD_QTL.csv", row.names = FALSE)

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
write.csv(coloc_results_clean, file = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/Coloc_results/MetaBrain_coloc_resultc/coloc_LBD_MADD_QTL.csv", row.names = FALSE)
