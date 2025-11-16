# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 11:19:10 2025

@author: user
"""
import pandas as pd
import os


def remove_matching_rows(data2, reference_file, output_path=None):
    """
    Remove rows from data2 where both 'SNP' and 'Trait' match entries in the reference file.

    Parameters:
        data2 (pd.DataFrame): Target DataFrame to be filtered.
        reference_file (str): Path to the reference CSV file containing failed SNPs.
        output_path (str, optional): Path to save the filtered DataFrame (only if not empty).

    Returns:
        pd.DataFrame: Filtered version of data2 (rows not matching reference entries).
    """
    # Load reference file (data1)
    data1 = pd.read_csv(reference_file)

    # Define matching keys
    keys = ["SNP", "Trait"]

    # Keep only rows in data2 not found in data1
    mask = ~data2[keys].apply(tuple, axis=1).isin(data1[keys].apply(tuple, axis=1))
    result = data2[mask].copy()

    # Summary of filtering
    removed_count = len(data2) - len(result)
    print(f"Removed {removed_count} rows; {len(result)} rows remain.")

    # Save file only if data remains
    if len(result) > 0 and output_path:
        result.to_csv(output_path, index=False)
        print(f"Saved filtered file: {output_path}")
    else:
        print("Warning: no remaining data after filtering. File not saved.")

    return result


def process_selected_files():
    """
    Process selected FUMA results files by removing rows matching the reference failed SNPs file.
    """
    base_dir = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results"
    input_dir1 = os.path.join(base_dir, "1_FUMA_results_table")
    input_dir2 = os.path.join(base_dir, "2_Merging_results")
    output_dir = os.path.join(base_dir, "3_Robust_filter_results")
    reference_file = os.path.join(output_dir, "FUMA_results_MADD_failed_snps.csv")

    os.makedirs(output_dir, exist_ok=True)

    # Files to process
    file_groups = {
        "Table": ["MADD_VD", "MADD_LBD", "MADD_DEM", "MADD_UD", "MADD_CP"],
        "Merged": ["MADD"]
    }

    # Process files in both directories
    for dir_label, subfolders in file_groups.items():
        input_dir = input_dir1 if dir_label == "Table" else input_dir2

        for name in subfolders:
            input_path = os.path.join(input_dir, f"FUMA_results_{name}.csv")
            output_path = os.path.join(output_dir, f"cleaned_FUMA_results_{name}.csv")

            if not os.path.exists(input_path):
                print(f"Warning: file not found: {input_path}")
                continue

            print(f"\nProcessing {dir_label} file: {name}")

            data2 = pd.read_csv(input_path)
            remove_matching_rows(data2, reference_file, output_path)

    print("\nAll selected files processed successfully!")


if __name__ == "__main__":
    process_selected_files()

