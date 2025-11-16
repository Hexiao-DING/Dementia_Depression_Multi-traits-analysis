# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 16:57:29 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 19:43:55 2025

@author: 86198
"""

import pandas as pd

def list_locus_categories_for_alzheimers(file_path, dis):
    """
    List all categories in the 'GenomicLocus' column when the 'Trait' column 
    contains a specified disease keyword (e.g., "Alzheimer's disease").

    Parameters:
    ----------
    file_path : str
        Path to the CSV or TXT file (tab-delimited).
    dis : str
        Disease keyword to search for (e.g., "Alzheimer's disease").

    Returns:
    -------
    list
        A list of unique categories from the 'GenomicLocus' column.
    """
    try:
        # Read file (assuming it's tab-separated)
        df = pd.read_csv(file_path, sep='\t')

        # Filter rows containing the disease keyword (case-insensitive)
        disease_df = df[df['Trait'].str.contains(dis, case=False, na=False)]

        # Extract unique values from 'GenomicLocus' column
        locus_categories = disease_df['GenomicLocus'].unique().tolist()

        # Print results
        print(f"File: {file_path}")
        print(f"Rows containing '{dis}': {len(disease_df)}")
        print(f"Number of unique Locus categories: {len(locus_categories)}")
        print("Locus category list:")
        for i, category in enumerate(locus_categories, 1):
            print(f"{i}. {category}")

        return locus_categories

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
        return []
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return []


# ==================== Example Usage ====================

# Alzheimer's disease (D)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_D\gwascatalog.txt"
dis = "Alzheimer's disease"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# D Novel Loci: 1 5 11 12 13 15

# Vascular dementia (VD)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_VD\gwascatalog.txt"
dis = "Dementia"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# AD/VD/LBD/MDD no Novel

# Cognitive performance (CP)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_CP\gwascatalog.txt"
dis = "Cognitive"  # "Cognitive performance" or "Cognitive ability"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# CP Novel Loci: 50 76 158

# Unipolar depression (UD)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_UD\gwascatalog.txt"
dis = "Alzheimer's disease"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# UD Novel Loci: 1 3

# Lewy body dementia (LBD)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_LBD\gwascatalog.txt"
dis = "Dementia"  # "Dementia with Lewy bodies"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)

# Depression disorder (DD)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_DD\gwascatalog.txt"
dis = "Depression"  # "Depressive disorder"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# DD Novel Loci: 1

# Mixed depression/anxiety disorder (MD)
file_path = r"D:\Projects_data&code\Stage1_bioinformatics_ADandDepression\Files\FUMA_results\4_FUMA_for_ld\MADD_MD\gwascatalog.txt"
dis = "anxiety"  # "Depressive disorder" or "Anxiety"
locus_categories = list_locus_categories_for_alzheimers(file_path, dis)
# MD Novel Loci: 1 2
