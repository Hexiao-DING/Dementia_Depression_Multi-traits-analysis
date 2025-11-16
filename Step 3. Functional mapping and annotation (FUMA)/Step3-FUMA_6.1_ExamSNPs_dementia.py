# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 01:41:28 2025

@author: user
"""
import pandas as pd

def merge_mtag_results(main_trait, depression_traits, gws_threshold=5e-8):
    """
    Merge MTAG results for main trait and multiple depression traits
    
    Parameters:
    main_trait (str): Abbreviation of main trait (e.g., "DEM")
    depression_traits (list): List of depression trait abbreviations (e.g., ["DD", "MDD", "MADD"])
    gws_threshold (float): Genome-wide significance threshold, default is 5e-8
    """
    main_trait_file = f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_{main_trait}.csv"
    base_dir = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/MTAG_results"
    output_path = f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/6_ExamSNPs/Examination_SNPs_{main_trait}_MTAG.csv"

    # 1. Load main trait dataset
    try:
        main_df = pd.read_csv(main_trait_file)
        print(f"Successfully loaded main trait data: {len(main_df)} rows")
    except Exception as e:
        print(f"Error loading main trait data: {str(e)}")
        return None
    
    # Rename main trait statistic columns
    main_df = main_df.rename(columns={
        "Beta": f"{main_trait}_Beta",
        "P-value": f"{main_trait}_P"
    })
    
    # 2. Load and merge depression trait datasets
    merged_df = main_df.copy()
    
    for trait in depression_traits:
        try:
            trait_file = f"{base_dir}/{trait}_{main_trait}/mtag_results_trait_1.txt"
            trait_df = pd.read_csv(trait_file, sep='\t')
            print(f"Successfully loaded {trait} data: {len(trait_df)} rows")
            
            # Rename depression trait statistic columns
            trait_df = trait_df.rename(columns={
                "mtag_beta": f"{trait}_Beta",
                "mtag_pval": f"{trait}_P"
            })
            
            # Merge into main dataset
            merged_df = pd.merge(
                merged_df,
                trait_df[["SNP", f"{trait}_Beta", f"{trait}_P"]],
                on="SNP",
                how="left"
            )
            
        except Exception as e:
            print(f"Error loading {trait} data: {str(e)}")
            continue
    
    # 3. Calculate CV_GWS_assoc column (check for genome-wide significant associations in depression traits)
    depression_p_cols = [f"{trait}_P" for trait in depression_traits]
    
    # Initialize column as "no"
    merged_df["CV_GWS_assoc"] = "no"
    
    # Check P-values for each depression trait
    for p_col in depression_p_cols:
        if p_col in merged_df.columns:
            # Check if p-value is below threshold and not missing
            trait_significant = merged_df[p_col].notna() & (merged_df[p_col] < gws_threshold)
            merged_df.loc[trait_significant, "CV_GWS_assoc"] = "yes"
    
    # 4. Adjust column order to target format
    # Determine available columns
    available_columns = merged_df.columns.tolist()
    
    # Base columns
    base_columns = ["Trait", "Unique", "SNP", "Chr", "Pos", "EA", "OA", "MAF", "Genomic Locus", "Gene", "CV_GWS_assoc"]
    
    # Statistic columns
    beta_columns = [f"{main_trait}_Beta"] + [f"{trait}_Beta" for trait in depression_traits]
    p_columns = [f"{main_trait}_P"] + [f"{trait}_P" for trait in depression_traits]
    
    # Keep only existing columns
    final_columns = []
    for col in base_columns + beta_columns + p_columns:
        if col in available_columns:
            final_columns.append(col)
    
    # Add any other columns that might exist
    other_columns = [col for col in available_columns if col not in final_columns]
    final_columns.extend(other_columns)
    
    result_df = merged_df[final_columns]
    
    # 5. Save results
    try:
        result_df.to_csv(output_path, index=False)
        print(f"Results saved to: {output_path}")
        
        # Validate results
        print("Result dataset preview:")
        print(result_df.head())
        print(f"\nCV_GWS_assoc distribution:")
        print(result_df["CV_GWS_assoc"].value_counts())
        
        return result_df
        
    except Exception as e:
        print(f"Error saving results: {str(e)}")
        return None

# Usage examples
if __name__ == "__main__":
    main_trait = "DEM"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "AD"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "CP"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "FD"
    depression_traits = ["DD", "MDD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "UD"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "VD"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )

if __name__ == "__main__":
    main_trait = "LBD"
    depression_traits = ["DD", "MDD", "MADD"]
    result = merge_mtag_results(
        main_trait=main_trait,
        depression_traits=depression_traits
    )
