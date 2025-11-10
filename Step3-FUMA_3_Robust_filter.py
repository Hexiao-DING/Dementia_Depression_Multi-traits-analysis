# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 22:11:54 2025

@author: user
"""

import pandas as pd
import numpy as np
import os


def filter_robust_snps(mtag_df, gwas_df):
    """
    Filter robust MTAG SNPs that satisfy two conditions:
    1. P-value < 0.01 in the original GWAS.
    2. Effect directions between MTAG and GWAS are consistent.

    Parameters:
        mtag_df : pd.DataFrame
            MTAG results containing ['SNP', 'A1', 'A2', 'beta', 'pval']
        gwas_df : pd.DataFrame
            GWAS results containing ['SNP', 'A1', 'A2', 'beta', 'pval']

    Returns:
        pd.DataFrame
            Filtered SNPs that passed both conditions.
    """
    # Merge MTAG and GWAS results by SNP
    merged = pd.merge(
        mtag_df,
        gwas_df,
        on='SNP',
        suffixes=('_MTAG', '_GWAS'),
        how='inner'
    )
    if merged.empty:
        return pd.DataFrame()

    # Condition 1: GWAS p-value < 0.01
    condition1 = merged['pval_GWAS'] <= 0.01

    # Identify allele alignment
    same_allele = (
        (merged['A1_MTAG'] == merged['A1_GWAS']) &
        (merged['A2_MTAG'] == merged['A2_GWAS'])
    )
    flipped_allele = (
        (merged['A1_MTAG'] == merged['A2_GWAS']) &
        (merged['A2_MTAG'] == merged['A1_GWAS'])
    )

    # Condition 2: Consistent effect directions
    direction_ok = (
        (same_allele & (merged['beta_MTAG'] * merged['beta_GWAS'] > 0)) |
        (flipped_allele & (merged['beta_MTAG'] * merged['beta_GWAS'] < 0))
    )

    # Apply both filters
    filtered_df = merged[condition1 & direction_ok].copy()

    # Mark allele orientation (same or flipped)
    filtered_df['allele_orientation'] = np.where(
        same_allele[condition1 & direction_ok],
        'same',
        'flipped'
    )

    # Select relevant columns for output
    result_cols = [
        'Trait', 'SNP', 'A1_MTAG', 'A2_MTAG', 'beta_MTAG', 'pval_MTAG',
        'A1_GWAS', 'A2_GWAS', 'beta_GWAS', 'pval_GWAS', 'allele_orientation'
    ]
    return filtered_df[result_cols]


def get_failure_reason(row, gwas_df):
    """
    Determine the reason why a SNP failed the robust filtering.

    Parameters:
        row : pd.Series
            One SNP row from MTAG results.
        gwas_df : pd.DataFrame
            Corresponding GWAS data.

    Returns:
        str : Reason for failure.
    """
    rsid = row['SNP']
    gwas_match = gwas_df[gwas_df['SNP'] == rsid]

    if gwas_match.empty:
        return "Not found in GWAS"

    gwas_row = gwas_match.iloc[0]

    if gwas_row['pval'] > 0.01:
        return "GWAS P-value > 0.01"

    same_allele = (row['A1'] == gwas_row['A1']) and (row['A2'] == gwas_row['A2'])
    flipped_allele = (row['A1'] == gwas_row['A2']) and (row['A2'] == gwas_row['A1'])

    if same_allele and (row['beta'] * gwas_row['beta'] < 0):
        return "Effect direction mismatch (same allele)"
    if flipped_allele and (row['beta'] * gwas_row['beta'] > 0):
        return "Effect direction mismatch (flipped allele)"

    return "Unknown reason"


if __name__ == "__main__":
    # Define dataset file paths
    DATASETS = {
        'DD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_DD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DDsummary.txt",
            'output_file': "FUMA_results_DD_failed_snps.csv"
        },
        'MADD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_MADD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MADDsummary.txt",
            'output_file': "FUMA_results_MADD_failed_snps.csv"
        },
        'MDD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_MDD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/MDDsummary.txt",
            'output_file': "FUMA_results_MDD_failed_snps.csv"
        },
        'AD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_AD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/ADsummary.txt",
            'output_file': "FUMA_results_AD_failed_snps.csv"
        },
        'DEM': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_DEM.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/DEMsummary.txt",
            'output_file': "FUMA_results_DEM_failed_snps.csv"
        },
        'UD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_UD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/UDsummary.txt",
            'output_file': "FUMA_results_UD_failed_snps.csv"
        },
        'VD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_VD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/VDsummary.txt",
            'output_file': "FUMA_results_VD_failed_snps.csv"
        },
        'CP': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_CP.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/CPsummary.txt",
            'output_file': "FUMA_results_CP_failed_snps.csv"
        },
        'LBD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_LBD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/LBDsummary.txt",
            'output_file': "FUMA_results_LBD_failed_snps.csv"
        },
        'FD': {
            'mtag_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/2_Merging_results/FUMA_results_FD.csv",
            'gwas_file': "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Data/GWAS/Final_processed_data/FDsummary.txt",
            'output_file': "FUMA_results_FD_failed_snps.csv"
        }
    }

    # Column mapping from MTAG output to unified naming
    COLUMN_MAP = {
        'Chr': 'CHR',
        'Pos': 'POS',
        'EA': 'A2',
        'OA': 'A1',
        'Beta': 'beta',
        'SE': 'se',
        'P-value': 'pval'
    }

    # Columns to include in the final failed SNP output
    DISPLAY_COLS = ['Trait', 'SNP', 'CHR', 'POS', 'A1', 'A2', 'beta', 'pval', 'failure_reason']

    # Output directory
    output_dir = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/3_Robust_filter_results"
    os.makedirs(output_dir, exist_ok=True)

    # Process each dataset
    for name, config in DATASETS.items():
        print(f"\n{'='*50}")
        print(f"Processing dataset: {name}")
        print(f"{'='*50}")

        # Load MTAG and GWAS data
        mtag_data = pd.read_csv(config['mtag_file'])
        gwas_df = pd.read_csv(config['gwas_file'], sep='\t')

        # Rename columns for consistency
        mtag_df = mtag_data.rename(columns=COLUMN_MAP)

        # Ensure Trait column exists; use the one from file if available
        if 'Trait' not in mtag_df.columns:
            mtag_df['Trait'] = name

        # Run robust filtering
        robust_snps = filter_robust_snps(mtag_df, gwas_df)

        # Determine failed SNPs
        all_snps = set(mtag_df['SNP'])
        passed_snps = set(robust_snps['SNP']) if not robust_snps.empty else set()
        failed_snps = all_snps - passed_snps

        if not failed_snps:
            print("All SNPs passed filtering.")
            continue

        failed_df = mtag_df[mtag_df['SNP'].isin(failed_snps)].copy()

        # Annotate failure reasons
        failed_df['failure_reason'] = failed_df.apply(
            lambda row: get_failure_reason(row, gwas_df), axis=1
        )

        # Summary statistics of failure reasons
        reason_counts = failed_df['failure_reason'].value_counts()
        print(f"Failed SNPs summary ({len(failed_df)} total):")
        for reason, count in reason_counts.items():
            print(f"  - {reason}: {count} SNPs")

        # Save failed SNPs to output file
        output_path = os.path.join(output_dir, config['output_file'])
        failed_df[DISPLAY_COLS].to_csv(output_path, index=False)
        print(f"Results saved to: {output_path}")

