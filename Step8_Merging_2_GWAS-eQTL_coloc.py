# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 17:14:32 2025

@author: user
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict

# ========== General Optimization Functions ==========
def create_traits_cache(base_path, traits_list):
    """Create file cache for all Traits to avoid repeated reading"""
    cache = defaultdict(dict)
    
    for traits in set(traits_list):
        trait_dir = os.path.join(base_path, traits)
        file1 = os.path.join(trait_dir, "mtag_results_trait_1.txt")
        file2 = os.path.join(trait_dir, "mtag_results_trait_2.txt")
        
        try:
            # Read both files at once and create SNP index
            df1 = pd.read_csv(file1, sep='\t').set_index('SNP')
            df2 = pd.read_csv(file2, sep='\t').set_index('SNP')
            cache[traits] = {
                'df1': df1[['CHR', 'BP', 'mtag_pval']],
                'df2': df2[['CHR', 'BP', 'mtag_pval']]
            }
        except Exception as e:
            print(f"Error loading {traits}: {str(e)}")
            cache[traits] = None
    
    return cache

def process_row(row, cache, snp_col='SNP'):
    """Process single row data using cache for efficiency"""
    traits = row['Traits']
    gwas = row['GWAS']
    
    # Skip invalid Traits cache
    if cache.get(traits) is None:
        return {f'pval1': np.nan, f'pval2': np.nan, 'CHR': np.nan, 'BP': np.nan}
    
    # Get target SNP
    snp = row[snp_col]
    trait_parts = traits.split('_')
    
    # Determine file mapping relationship
    if gwas == trait_parts[0]:
        target_df = cache[traits]['df1']
        other_df = cache[traits]['df2']
    elif gwas == trait_parts[1]:
        target_df = cache[traits]['df2']
        other_df = cache[traits]['df1']
    else:
        return {f'pval1': np.nan, f'pval2': np.nan, 'CHR': np.nan, 'BP': np.nan}
    
    # Use index for fast lookup (.loc[] is 100x faster than linear search)
    result = {'CHR': np.nan, 'BP': np.nan, 'pval1': np.nan, 'pval2': np.nan}
    
    if snp in target_df.index:
        target_data = target_df.loc[snp]
        result['CHR'] = target_data['CHR']
        result['BP'] = target_data['BP']
        result['pval1'] = target_data['mtag_pval']
    
    if snp in other_df.index:
        result['pval2'] = other_df.loc[snp]['mtag_pval']
    
    return result

# ========== Optimized Main Functions ==========
def match_snp_data(df):
    rename_dict = {
        'Target_SNP': 'SNP',
        'Gene_Ensembl': 'Gene (Ensembl)',
        'Gene_symbol': 'Gene symbol',
        'QTL': 'Tissue',
        'Causal_SNP': 'Causal SNP',
        'PP_explained': '%PP explained',
        'regional_prob': 'RP'
    }
    df = df.rename(columns=rename_dict)
    
    base_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\MTAG_results"
    
    # Step 1: Create file cache for all Traits
    traits_list = df['Traits'].tolist()
    cache = create_traits_cache(base_path, traits_list)
    
    # Step 2: Batch process all rows
    results = []
    for _, row in df.iterrows():
        results.append(process_row(row, cache, snp_col='SNP'))
    
    # Step 3: Merge results
    result_df = pd.DataFrame(results)
    df = pd.concat([df, result_df], axis=1)
    
    # Step 4: Calculate flag column
    df['MTAG_targetSNP_sig'] = df.apply(
        lambda row: 'yes' if not pd.isna(row['pval1']) and not pd.isna(row['pval2']) 
                   and row['pval1'] < 0.01 and row['pval2'] < 0.01 
            else 'no' if not pd.isna(row['pval1']) and not pd.isna(row['pval2']) 
                    and row['pval1'] >= 0.01 and row['pval2'] >= 0.01 
            else 'uncertain', axis=1
    )
    
    final_columns = [
        "SNP", "CHR", "BP", "Traits", "Gene (Ensembl)", "Gene symbol", 
        "Tissue", "PP", "Causal SNP", "%PP explained", "RP", 'MTAG_targetSNP_sig', 'GWAS'
    ]
    
    return df[final_columns]

def match_sig_col(df):
    base_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\MTAG_results"
    
    # Step 1: Create file cache for all Traits
    traits_list = df['Traits'].tolist()
    cache = create_traits_cache(base_path, traits_list)
    
    # Step 2: Batch process all rows
    results = []
    for _, row in df.iterrows():
        # Use same process_row function, but specify snp_col as 'Causal SNP'
        res = process_row(row, cache, snp_col='Causal SNP')
        results.append({'pval1': res['pval1'], 'pval2': res['pval2']})
    
    # Step 3: Merge results
    result_df = pd.DataFrame(results)
    df = pd.concat([df, result_df], axis=1)
    
    # Step 4: Calculate flag column
    df['MTAG_Causal_SNP_sig'] = df.apply(
        lambda row: 'yes' if not pd.isna(row['pval1']) and not pd.isna(row['pval2']) 
                   and row['pval1'] < 0.01 and row['pval2'] < 0.01 
            else 'no' if not pd.isna(row['pval1']) and not pd.isna(row['pval2']) 
                    and row['pval1'] >= 0.01 and row['pval2'] >= 0.01 
            else 'uncertain', axis=1
    )
        
    return df

# Main program
if __name__ == "__main__":
    base_dir = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\Coloc_results"
    file = 'coloc_QTL_symbol.csv'
    coloc_df = pd.read_csv(os.path.join(base_dir, file))
    #coloc_df = coloc_df[coloc_df['QTL'].str.contains('Brain', case=False, na=False)]
    enhanced_df = match_snp_data(coloc_df)
    df = match_sig_col(enhanced_df)
    #filtered_df = df[df['Tissue'].str.contains('Brain', case=False, na=False)]
    output_path = os.path.join(base_dir, "coloc_GWAS-QTL.csv")
    df.to_csv(output_path, index=False)