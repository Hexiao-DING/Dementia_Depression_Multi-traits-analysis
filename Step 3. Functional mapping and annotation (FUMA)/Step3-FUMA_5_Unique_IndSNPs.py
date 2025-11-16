# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 16:35:04 2025

@author: user
"""

import os
import pandas as pd


def process_fuma_data(job_dir, job_name="AF"):
    """
    Process FUMA output files and generate standardized SNP summary tables.

    Parameters:
        job_dir (str): Directory containing FUMA output files.
        job_name (str): Name of the FUMA analysis (default = "AF").

    Returns:
        pandas.DataFrame or None: Processed result DataFrame, or None if failed.
    """
    try:
        # Define input file paths
        file_paths = {
            "lead_snps": os.path.join(job_dir, "leadSNPs.txt"),
            "ind_sig_snps": os.path.join(job_dir, "IndSigSNPs.txt"),
            "snps": os.path.join(job_dir, "snps.txt"),
            "annovar": os.path.join(job_dir, "annov.txt"),
            "genomic_loci": os.path.join(job_dir, "GenomicRiskLoci.txt"),
            "mapped_genes": os.path.join(job_dir, "genes.txt"),
        }

        print("Reading FUMA output files...")
        lead_snps = pd.read_csv(file_paths["lead_snps"], sep="\t", dtype={"rsID": str})
        ind_sig_snps = pd.read_csv(file_paths["ind_sig_snps"], sep="\t", dtype={"rsID": str})
        snps = pd.read_csv(file_paths["snps"], sep="\t", dtype={"rsID": str})
        annovar = pd.read_csv(file_paths["annovar"], sep="\t")
        genomic_loci = pd.read_csv(file_paths["genomic_loci"], sep="\t")

        # Clean chromosome and position columns
        for df in [ind_sig_snps, snps, annovar]:
            if "chr" in df.columns:
                df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False).str.strip()
            if "pos" in df.columns:
                df["pos"] = df["pos"].astype(int)

        # Map SNPs to genomic loci
        print("Mapping SNPs to genomic loci...")
        snp_to_locus = {}
        for _, row in genomic_loci.iterrows():
            if pd.notna(row["IndSigSNPs"]):
                for snp in str(row["IndSigSNPs"]).split(";"):
                    snp = snp.strip()
                    if snp:
                        snp_to_locus[snp] = row["GenomicLocus"]

        # Process ANNOVAR annotations (nearest gene)
        print("Processing ANNOVAR annotations...")
        annovar["abs_dist"] = annovar["dist"].abs()
        annovar_nearest = (
            annovar.sort_values(["chr", "pos", "abs_dist"])
            .drop_duplicates(["chr", "pos"])
        )

        # Merge SNP-level info
        print("Merging SNP-level data...")
        result_df = ind_sig_snps[["rsID", "chr", "pos", "p"]].copy()
        result_df = pd.merge(
            result_df,
            snps[["rsID", "chr", "pos", "effect_allele", "non_effect_allele", "MAF", "beta", "se"]],
            on=["rsID", "chr", "pos"],
            how="left",
        )

        # Add gene annotation
        result_df = pd.merge(
            result_df,
            annovar_nearest[["chr", "pos", "symbol", "dist", "annot"]].rename(
                columns={
                    "symbol": "gene_name",
                    "dist": "distance",
                    "annot": "Function",
                }
            ),
            on=["chr", "pos"],
            how="left",
        )

        # Add locus and signal type
        result_df["Genomic Locus"] = result_df["rsID"].map(snp_to_locus)
        result_df["Signal type"] = result_df["rsID"].isin(lead_snps["rsID"]).map(
            {True: "Top", False: "Secondary"}
        )

        # Compute gene distance direction
        result_df["Direction"] = result_df["distance"].apply(
            lambda x: "Upstream" if pd.notna(x) and x < 0 else ("Downstream" if pd.notna(x) and x > 0 else "")
        )
        result_df["abs_distance"] = result_df["distance"].abs()

        # Build formatted table
        print("Formatting final result table...")
        final_df = pd.DataFrame({
            "Trait": job_name,
            "SNP": result_df["rsID"],
            "Chr": result_df["chr"],
            "Pos": result_df["pos"],
            "EA": result_df["effect_allele"],
            "OA": result_df["non_effect_allele"],
            "MAF": result_df["MAF"].round(4),
            "Beta": result_df["beta"].round(6),
            "SE": result_df["se"].round(6),
            "P-value": result_df["p"].apply(lambda x: f"{float(x):.2e}"),
            "Signal type": result_df["Signal type"],
            "Genomic Locus": result_df["Genomic Locus"].fillna("NA"),
            "Gene": result_df["gene_name"].fillna("Intergenic"),
            "Distance": result_df.apply(
                lambda r: f"{r['abs_distance']}bp {r['Direction']}".strip()
                if pd.notna(r["abs_distance"]) else "Intergenic",
                axis=1,
            ),
            "Function": result_df["Function"].fillna("Intergenic"),
        })

        # Replace missing numeric values with "NA"
        for col in ["MAF", "Beta", "SE"]:
            final_df[col] = final_df[col].apply(lambda x: "NA" if pd.isna(x) else x)

        # Sort for readability
        final_df = final_df.sort_values(["Chr", "Pos"]).reset_index(drop=True)

        print(f"Successfully processed {job_name}: {len(final_df)} SNPs")
        return final_df

    except Exception as e:
        print(f"Error processing {job_name}: {e}")
        return None


if __name__ == "__main__":
    # Define job directories and names
    base_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld"
    job_dirs = {
        "DEM": os.path.join(base_path, "DEM"),
        "AD": os.path.join(base_path, "AD"),
        "CP": os.path.join(base_path, "CP"),
        "FD": os.path.join(base_path, "FD"),
        "UD": os.path.join(base_path, "UD"),
        "VD": os.path.join(base_path, "VD"),
        "LBD": os.path.join(base_path, "LBD"),
        "DD": os.path.join(base_path, "DD"),
        "MDD": os.path.join(base_path, "MDD"),
        "MADD": os.path.join(base_path, "MADD"),
    }

    output_dir = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\5_Final_IndSNPs"
    os.makedirs(output_dir, exist_ok=True)

    all_results = []

    # Process each FUMA job
    for job_name, job_dir in job_dirs.items():
        print(f"\n{'='*60}\nProcessing: {job_name}\nPath: {job_dir}")
        df = process_fuma_data(job_dir, job_name)
        if df is not None:
            all_results.append(df)
            output_file = os.path.join(output_dir, f"final_independent_snps_{job_name}.csv")
            df.to_csv(output_file, index=False)
            print(f"Saved: {output_file}")
        else:
            print(f"Skipped {job_name} due to errors.")

    # Optionally merge all results
    # if all_results:
    #     combined_file = os.path.join(output_dir, "FUMA_results_combined.csv")
    #     pd.concat(all_results).to_csv(combined_file, index=False)
    #     print(f"\nCombined results saved to: {combined_file}")
