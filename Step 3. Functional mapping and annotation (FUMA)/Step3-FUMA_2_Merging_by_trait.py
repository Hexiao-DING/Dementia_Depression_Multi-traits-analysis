# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 19:23:31 2025

@author: user
"""

import pandas as pd

# ---------------------------------------------
# Part 1: Combine results for single-trait analysis
# ---------------------------------------------

traits = ["DEM", "AD", "UD", "VD", "LBD", "CP", "FD"]

for trait in traits:
    all_snp_info_dfs = []  # store dataframes for current trait
    
    # Subfolder structure per trait
    subfolders = [f"{trait}_DD", f"{trait}_MDD", f"{trait}_MADD"]
    
    for subfolder in subfolders:
        # Updated input path
        snp_info_file = (
            f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/"
            f"Files/FUMA_results/1_FUMA_results_table/FUMA_results_{subfolder}.csv"
        )
        
        try:
            snp_info_df = pd.read_csv(snp_info_file)
            print(f"Loaded {subfolder}: {len(snp_info_df)} rows")
            all_snp_info_dfs.append(snp_info_df)
        except FileNotFoundError:
            print(f"Warning: file not found: {snp_info_file}")
            continue
        except Exception as e:
            print(f"Warning: error reading {snp_info_file}: {str(e)}")
            continue
    
    if not all_snp_info_dfs:
        print(f"Warning: no valid files found for {trait}")
        continue
    
    # Combine all dataframes
    combined_df = pd.concat(all_snp_info_dfs, axis=0, ignore_index=True)
    print(f"Total merged rows for {trait}: {len(combined_df)}")
    
    # Output path remains the same
    output_file = (
        f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/"
        f"Files/FUMA_results/2_Merging_results/FUMA_results_{trait}.csv"
    )
    combined_df.to_csv(output_file, index=False)
    print(f"Saved combined file: {output_file}\n")


# ---------------------------------------------
# Part 2: Combine results for dual-trait analysis
# ---------------------------------------------

traits = ["DD", "MADD", "MDD"]

for trait in traits:
    all_snp_info_dfs = []  # store dataframes for current trait
    
    subfolders = [
        f"{trait}_DEM", f"{trait}_AD", f"{trait}_UD", f"{trait}_VD",
        f"{trait}_LBD", f"{trait}_CP", f"{trait}_FD"
    ]
    
    for subfolder in subfolders:
        # Updated input path
        snp_info_file = (
            f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/"
            f"Files/FUMA_results/1_FUMA_results_table/FUMA_results_{subfolder}.csv"
        )
        
        try:
            snp_info_df = pd.read_csv(snp_info_file)
            print(f"Loaded {subfolder}: {len(snp_info_df)} rows")
            all_snp_info_dfs.append(snp_info_df)
        except FileNotFoundError:
            print(f"Warning: file not found: {snp_info_file}")
            continue
        except Exception as e:
            print(f"Warning: error reading {snp_info_file}: {str(e)}")
            continue
    
    if not all_snp_info_dfs:
        print(f"Warning: no valid files found for {trait}")
        continue
    
    # Combine all dataframes
    combined_df = pd.concat(all_snp_info_dfs, axis=0, ignore_index=True)
    print(f"Total merged rows for {trait}: {len(combined_df)}")
    
    # Output path remains the same
    output_file = (
        f"D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/"
        f"Files/FUMA_results/2_Merging_results/FUMA_results_{trait}.csv"
    )
    combined_df.to_csv(output_file, index=False)
    print(f"Saved combined file: {output_file}\n")
