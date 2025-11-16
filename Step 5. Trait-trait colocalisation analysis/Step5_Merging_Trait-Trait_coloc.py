# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:26:27 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:26:27 2025

@author: user
"""

import pandas as pd
import os

def merge_datasets(coloc_files, snp_files, output_path):
    """
    Merge multiple datasets and perform data processing.
    
    Parameters:
    coloc_files: List of full paths to coloc result files
    snp_files: List of full paths to SNP information files
    output_path: Final output file path
    """
    # 1. Read all coloc datasets
    coloc_dfs = [pd.read_csv(file) for file in coloc_files]
    
    # Combine all coloc data
    combined_coloc = pd.concat(coloc_dfs, ignore_index=True)
    
    # 2. Rename statistical columns
    combined_coloc = combined_coloc.rename(columns={
        "candidate_snps": "Candidate Causal SNP",
        "posterior_probs": "PP",
        "PP_explained": "%PP explained",
        "regional_prob": "RP"
    })
    
    # 3. Read all SNP information datasets
    snp_dfs = [pd.read_csv(file) for file in snp_files]
    
    # Combine all SNP information
    combined_snp = pd.concat(snp_dfs, ignore_index=True)
    
    # 4. Merge coloc results with SNP information
    merged = pd.merge(
        combined_coloc[["SNP", 'Candidate Causal SNP', 'PP', '%PP explained', 'RP', 'Traits']],
        combined_snp[["SNP", "Trait", "Chr", "Pos", "Gene"]],
        on=["SNP"],
        how="left"
    )
    
    # 5. Adjust column order
    final_columns = [
        "SNP", "Chr", "Pos", "Gene", "Trait", "Traits",
        "PP", "Candidate Causal SNP", "%PP explained", "RP"
    ]
    
    # 6. Remove duplicates and save results
    result = merged[final_columns].drop_duplicates()
    result.to_csv(output_path, index=False)
    
    print(f"Processing completed! Merged {len(coloc_files)} coloc files and {len(snp_files)} SNP files.")
    print(f"The final dataset contains {len(result)} unique SNP records.")
    return result

# ---------------- Example usage ----------------
if __name__ == "__main__":
    # List of full paths for coloc files
    coloc_files = [
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_DEM_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_DEM_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_CP_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_CP_MDD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_CP_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_UD_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_VD_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_VD_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\Raw_coloc_results\coloc_LBD_MADD.csv"
    ]
    
    # List of full paths for SNP files
    snp_files = [
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_DEM_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_MADD_DEM.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_DEM_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_CP_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_MDD_CP.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_CP_MDD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_MADD_CP.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_CP_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_UD_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_DD_VD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_VD_DD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_VD_MADD.csv",
        r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\1_FUMA_results_table\FUMA_results_LBD_MADD.csv"
    ]
    
    # Output path
    output_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results\coloc_Trait-Trait.csv"
    
    # Execute merge
    result_df = merge_datasets(coloc_files, snp_files, output_path)
    
    # View results
    print(result_df.head())
