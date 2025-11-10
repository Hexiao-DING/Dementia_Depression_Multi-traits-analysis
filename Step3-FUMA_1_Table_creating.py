# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 02:54:00 2025

@author: user
"""


import os
import pandas as pd
import re


def process_fuma_data(job_dir, job_name="AF"):
    """
    Process FUMA output data and generate standardized summary tables.
    
    Parameters:
        job_dir: str, directory containing FUMA output files
        job_name: str, name of the FUMA analysis job (default = "AF")
    
    Returns:
        pandas.DataFrame: formatted summary table
    """
    try:
        file_paths = {
            'lead_snps': os.path.join(job_dir, 'leadSNPs.txt'),
            'ind_sig_snps': os.path.join(job_dir, 'IndSigSNPs.txt'),
            'snps': os.path.join(job_dir, 'snps.txt'),
            'annovar': os.path.join(job_dir, 'annov.txt'),
            'genomic_loci': os.path.join(job_dir, 'GenomicRiskLoci.txt'),
            'mapped_genes': os.path.join(job_dir, 'genes.txt')
        }

        print("Reading input files...")
        lead_snps = pd.read_csv(file_paths['lead_snps'], sep='\t', dtype={'rsID': str})
        ind_sig_snps = pd.read_csv(file_paths['ind_sig_snps'], sep='\t', dtype={'rsID': str})
        snps = pd.read_csv(file_paths['snps'], sep='\t', dtype={'rsID': str})
        annovar = pd.read_csv(file_paths['annovar'], sep='\t')
        genomic_loci = pd.read_csv(file_paths['genomic_loci'], sep='\t')
        mapped_genes = pd.read_csv(file_paths['mapped_genes'], sep='\t')

        # Ensure chromosome columns are standardized
        for df in [ind_sig_snps, snps, annovar]:
            if 'chr' in df.columns:
                df['chr'] = df['chr'].astype(str).str.replace('chr', '').str.strip()
            if 'pos' in df.columns:
                df['pos'] = df['pos'].astype(int)

        # Map SNPs to genomic risk loci
        print("Building genomic risk locus mapping...")
        snp_to_locus = {}
        for _, row in genomic_loci.iterrows():
            if pd.notna(row['IndSigSNPs']):
                for snp in str(row['IndSigSNPs']).split(';'):
                    snp = snp.strip()
                    if snp:
                        snp_to_locus[snp] = row['GenomicLocus']

        # Process ANNOVAR annotation data - keep nearest gene per SNP
        print("Processing ANNOVAR annotation data...")
        annovar['abs_dist'] = annovar['dist'].abs()
        annovar_nearest = annovar.sort_values(['chr', 'pos', 'abs_dist']).drop_duplicates(['chr', 'pos'])

        # Merge SNP base information
        print("Merging SNP data...")
        result_df = ind_sig_snps[['rsID', 'chr', 'pos', 'p']].copy()
        result_df = pd.merge(
            result_df,
            snps[['rsID', 'chr', 'pos', 'effect_allele', 'non_effect_allele', 'MAF', 'beta', 'se']],
            on=['rsID', 'chr', 'pos'],
            how='left'
        )

        # Add gene annotation info
        print("Adding gene annotation info...")
        result_df = pd.merge(
            result_df,
            annovar_nearest[['chr', 'pos', 'symbol', 'dist', 'annot']].rename(columns={
                'symbol': 'gene_name',
                'dist': 'distance',
                'annot': 'Function'
            }),
            on=['chr', 'pos'],
            how='left'
        )

        # Add genomic locus and signal type info
        print("Adding genomic risk locus info...")
        result_df['GenomicLocus'] = result_df['rsID'].map(snp_to_locus)
        result_df['Signal type'] = result_df['rsID'].isin(lead_snps['rsID']).map({True: 'Top', False: 'Secondary'})

        # Handle distance direction
        result_df['Direction'] = result_df['distance'].apply(
            lambda x: 'Upstream' if pd.notna(x) and x < 0 else ('Downstream' if pd.notna(x) and x > 0 else '')
        )
        result_df['abs_distance'] = result_df['distance'].abs()

        # Build final output DataFrame
        print("Building final result table...")
        final_df = pd.DataFrame({
            'Trait': job_name,
            'SNP': result_df['rsID'],
            'Chr': result_df['chr'],
            'Pos': result_df['pos'],
            'EA': result_df['effect_allele'],
            'OA': result_df['non_effect_allele'],
            'MAF': result_df['MAF'].round(4),
            'Beta': result_df['beta'].round(6),
            'SE': result_df['se'].round(6),
            'P-value': result_df['p'].apply(lambda x: f"{float(x):.2e}"),
            'Signal type': result_df['Signal type'],
            'Genomic Locus': result_df['GenomicLocus'],
            'Gene': result_df['gene_name'].fillna('Intergenic'),
            'Distance': result_df.apply(
                lambda r: f"{r['abs_distance']}bp {r['Direction']}".strip()
                if pd.notna(r['abs_distance'])
                else 'Intergenic',
                axis=1
            ),
            'Function': result_df['Function'].fillna('Intergenic')
        })

        # Handle missing values
        for col in ['MAF', 'Beta', 'SE']:
            final_df[col] = final_df[col].apply(lambda x: 'NA' if pd.isna(x) else x)

        final_df.fillna({
            'Genomic Locus': 'NA',
            'Gene': 'Intergenic',
            'Distance': 'Intergenic',
            'Function': 'Intergenic'
        }, inplace=True)

        # Sort and finalize
        final_df = final_df.sort_values(['Chr', 'Pos'])
        final_df.reset_index(drop=True, inplace=True)
        
        print("Data processing completed successfully!")
        return final_df

    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")
        return None


# Batch processing of multiple FUMA job folders
if __name__ == "__main__":
    job_dirs = {
    # Trait pairings anchored on depression disorder cohorts
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/DD_VD": "DD_VD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/DD_AD": "DD_AD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/AD_DD": "AD_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/DEM_DD": "DEM_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/UD_DD": "UD_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/VD_DD": "VD_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/CP_DD": "CP_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/LBD_DD": "LBD_DD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/FD_DD": "FD_DD",

    # Trait pairings anchored on mixed anxiety and depression cohorts
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MADD_VD": "MADD_VD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MADD_LBD": "MADD_LBD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MADD_DEM": "MADD_DEM",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MADD_UD": "MADD_UD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MADD_CP": "MADD_CP",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/AD_MADD": "AD_MADD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/DEM_MADD": "DEM_MADD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/UD_MADD": "UD_MADD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/VD_MADD": "VD_MADD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/CP_MADD": "CP_MADD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/LBD_MADD": "LBD_MADD",

    # Trait pairings anchored on major depressive disorder cohorts
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_VD": "MDD_VD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_AD": "MDD_AD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_FD": "MDD_FD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_LBD": "MDD_LBD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_CP": "MDD_CP",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_UD": "MDD_UD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/MDD_DEM": "MDD_DEM",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/AD_MDD": "AD_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/DEM_MDD": "DEM_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/UD_MDD": "UD_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/VD_MDD": "VD_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/CP_MDD": "CP_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/LBD_MDD": "LBD_MDD",
    "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/Raw_FUMA_results/FD_MDD": "FD_MDD"
}

    
    all_results = []
    
    output_main_dir = "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/Files/FUMA_results/1_FUMA_results_table"
    os.makedirs(output_main_dir, exist_ok=True)
    
    for job_dir, job_name in job_dirs.items():
        print(f"\n{'='*50}")
        print(f"Processing FUMA job: {job_name}\nDirectory: {job_dir}")
        
        job_result = process_fuma_data(job_dir=job_dir, job_name=job_name)
        
        if job_result is not None:
            all_results.append(job_result)
            output_file = os.path.join(output_main_dir, f"FUMA_results_{job_name}.csv")
            job_result.to_csv(output_file, index=False)
            
            print(f"\nSaved results to {output_file}")
            print("\nSummary:")
            print(f"- Total SNPs: {len(job_result)}")
            print(f"- Top signals: {sum(job_result['Signal type'] == 'Top')}")
            print(f"- Genes involved: {job_result['Gene'].nunique()}")
            print(f"- Genomic loci: {job_result['Genomic Locus'].nunique()}")
            print("\nPreview:")
            print(job_result.head(5).to_string(index=False))
        else:
            print(f"Error processing job {job_name}")

        
        

