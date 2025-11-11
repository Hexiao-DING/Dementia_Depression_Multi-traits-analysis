# Dementia-Depression Multi-trait Analysis Pipeline

This repository contains a comprehensive, production-ready workflow for investigating the shared genetic architecture between dementia endophenotypes and depression-spectrum traits. The code base covers the complete arc from genome-wide association study (GWAS) preprocessing, through multi-trait meta-analysis (MTAG) and FUMA-based functional annotation, to HyPrColoc colocalisation against tissue-resolved eQTL resources. 

## Table of Contents

- [Background](#background)
- [Trait Catalogue and GWAS Sources](#trait-catalogue-and-gwas-sources)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Initial Setup](#initial-setup)
- [Data Preparation](#data-preparation)
- [Analytical Modules](#analytical-modules)
- [Key Outputs](#key-outputs)
- [Reproducibility Guidelines](#reproducibility-guidelines)
- [Acknowledgements](#acknowledgements)

## Background

We integrate multiple GWAS cohorts spanning dementias (Alzheimer's disease, vascular dementia, Lewy body dementia, frontotemporal dementia, undefined dementia, and composite FinnGen endpoints) and depression-related traits (depressive disorders, major depressive disorder, mixed anxiety and depression). The primary goal is to characterise pleiotropic loci and regulatory mechanisms that underpin the comorbidity between neurodegeneration and psychiatric phenotypes. MTAG is used to boost discovery power across correlated traits, FUMA supplies locus-level annotation and gene prioritisation, and HyPrColoc quantifies the evidence for shared causal variants between MTAG signals and tissue-specific eQTL panels.

## Trait Catalogue and GWAS Sources

All analyses rely on publicly available or FinnGen-released GWAS summary statistics. File names correspond to the archives stored under `Data/GWAS/Downloaded_original_data/`. Please review and respect each data provider's terms of use before redistribution.

### Dementia Traits (7 outcomes)

| Trait | Summary Statistics | Source | Primary Citation |
| ----- | ------------------ | ------ | ---------------- |
| Dementia (FinnGen endpoint F5_DEMENTIA) | `Finn-b-F5_DEMENTIA.tsv.gz` | [FinnGen](https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA) | Kurki, M. I. *et al.* (2023) *Nature* 613, 508-518. |
| Alzheimer's disease | `GCST90012877.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589840) | Schwartzentruber, J. *et al.* (2021) *Nat. Genet.* 53, 392-402. |
| Cognitive performance | `GCST006572.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/30038396) | Lee, J. J. *et al.* (2018) *Nat. Genet.* 50, 1112-1121. |
| Vascular dementia (FinnGen endpoint F5_VASCDEM) | `Finn-b-F5_VASCDEM.tsv.gz` | [FinnGen](https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM) | Kurki, M. I. *et al.* (2023) *Nature* 613, 508-518. |
| Lewy body dementia | `GCST90001390.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589841) | Chia, R. *et al.* (2021) *Nat. Genet.* 53, 294-303. |
| Frontotemporal dementia | `GCST90558311.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40280976) | Pottier, C. *et al.* (2025) *Nat. Commun.* 16, 3914. |
| Undefined dementia (FinnGen endpoint F5_DEMNAS) | `Finn_b_F5_DEMNAS.tsv.gz` | [FinnGen](https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS) | Kurki, M. I. *et al.* (2023) *Nature* 613, 508-518. |

### Depression Traits (3 outcomes)

| Trait | Summary Statistics | Source | Primary Citation |
| ----- | ------------------ | ------ | ---------------- |
| Depressive disorders | `GCST90476922.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39024449) | Verma, A. *et al.* (2024) *Science* 385, eadj1182. |
| Major depressive disorder | `GCST90468123.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39789286) | Loya, H. *et al.* (2025) *Nat. Genet.* 57, 461-468. |
| Mixed anxiety and depressive disorder | `GCST90225526.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37259642) | Brasher, M. S. *et al.* (2023) *Genes Brain Behav.* 22, e12851. |

## Repository Structure

- `R/` – Quality control, reference-genome liftover, and HyPrColoc analysis scripts written in R.
- `python/` – Python utilities for FUMA result aggregation, filtering, and MTAG-eQTL integration steps.
- `docs/` – Environment setup notes and user guides (MTAG installation on Linux, SMR replication instructions).
- `figs/` – Placeholder directory for derived figures or summary plots.
- `Others(delete)/` – Legacy working directory retained from the source project (safe to remove if not required).
- `README.md` – This documentation file.

## Prerequisites

### R Environment

- R (>= 4.2 recommended).
- Core packages: `data.table`, `dplyr`, `tidyr`, `stringr`, `duckdb`, `DBI`, `MungeSumstats`, `rtracklayer`, `BiocManager`.
- Bioconductor data dependencies: `SNPlocs.Hsapiens.dbSNP144.GRCh37`, chain files for genomic liftover (for example `hg38ToHg19.over.chain`), HyPrColoc (`hyprcoloc`).

### Python Environment

- Python (>= 3.8).
- Packages: `pandas`, `numpy`, `pyarrow` (optional but improves I/O), plus any modules referenced at the top of each script.

### External Tools

- MTAG (Jon Jala et al.) for multi-trait GWAS meta-analysis.
- FUMA (Functional Mapping and Annotation) web outputs (scripts assume results downloaded locally).
- SMR (Summary-data-based Mendelian Randomization) for replication analyses (steps detailed in `docs/Step4-NovelLoci_2_replication_SMR.txt`).
- PLINK reference data (1KG EUR) and brain eQTL summary statistics (for example MetaBrain or GTEx) as specified by the scripts.

## Initial Setup

1. **Clone the repository**
   ```bash
   git clone <your-fork-or-remote>
   ```
2. **Create conda environments** (optional but recommended)
   ```bash
   conda create -n gwas_r r-base=4.2 -y
   conda create -n gwas_py python=3.10 -y
   ```
3. **Install R dependencies** within the R environment. Example:
   ```r
   install.packages(c("data.table", "dplyr", "tidyr", "stringr", "duckdb", "DBI", "httr"))
   BiocManager::install(c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP144.GRCh37", "AnnotationHub"))
   install.packages("hyprcoloc", repos = "https://cloud.r-project.org")
   ```
4. **Install Python dependencies**
   ```bash
   conda activate gwas_py
   pip install pandas numpy pyarrow
   ```
5. **Review documentation** under `docs/` for tool-specific setup (e.g., MTAG on Linux, SMR replication workflow).

## Data Preparation

The scripts reference absolute Windows-style paths rooted at `D:/Projects_data&code/Stage1_bioinformatics_ADandDepression`. Adjust these paths to match your environment by either:

- Mirroring the original directory hierarchy, or
- Editing the path constants inside the scripts to point to your local data directories.

Required data assets include:

- Raw GWAS summary statistics (`Data/GWAS/Downloaded_original_data/`).
- Intermediate directories for formatted summary statistics (`Format_sumstats_data/`, `Final_processed_data/`, `MTAG_input_data/`, `LDSC_input_data/`).
- MTAG output folders under `Files/MTAG_results/`.
- FUMA outputs organised as `Files/FUMA_results/Raw_FUMA_results/`.
- Colocalisation results under `Files/Coloc_results/`.
- Brain eQTL panels (e.g., MetaBrain or GTEx v7) stored in `Data/eQTL/`.

Ensure these directories exist prior to running any script; most scripts will not attempt to create missing parent directories.

## Analytical Modules

The pipeline is divided into modular steps, each designed to solve a specific analytical question. You can run the full sequence or select individual components depending on your study design.

1. **GWAS Quality Control and Harmonisation (`R/Step1-QC_GWAS_*`)**
   - Objective: Standardise public GWAS summary statistics into a harmonised format suitable for joint analysis.
   - Logic: Apply allele orientation checks, compute beta and standard error (from odds ratios where needed), impute minor allele frequency, filter for data completeness, and perform liftover between GRCh38 and GRCh37 builds. Outputs feed MTAG and LDSC.
2. **MTAG Configuration and Execution (`docs/Step2-MTAG_linux.txt`)**
   - Objective: Provide reproducible instructions for running MTAG on Linux, ensuring consistent sample-size handling and trait alignment.
   - Logic: Create Python 2.7 environment, install MTAG dependencies, stage harmonised input files, and launch MTAG with parameters optimised for correlating dementia and depression traits.
3. **FUMA Post-processing Suite (`python/Step3-*`)**
   - `Step3-FUMA_1_Table_creating.py`: Batch-extract per-locus tables from raw FUMA downloads to generate harmonised CSV outputs.
   - `Step3-FUMA_2_Merging_by_trait.py`: Merge locus annotations across trait pairs to compare shared signals and facilitate downstream prioritisation.
   - `Step3-FUMA_3_Robust_filter.py`: Evaluate MTAG signals against original GWAS for robustness, flagging SNPs that fail concordance checks.
   - `Step3-FUMA_4_Delete_failed.py`: Remove FUMA entries that correspond to non-robust SNPs, keeping only high-confidence loci.
   - `Step3-FUMA_5_Unique_IndSNPs.py`: Derive a deduplicated list of independent lead SNPs per trait combination.
   - `Step3-FUMA_6.1_ExamSNPs_dementia.py` and `Step3-FUMA_6.2_ExamSNPs_depression.py`: Generate trait-specific inspection tables for manual review.
4. **Novel Locus Nomination (`python/Step4-NovelLoci_1_definition.py`)**
   - Objective: Identify novel loci by integrating filtered FUMA results with additional heuristics (distance to known loci, multi-trait overlap).
   - Logic: Flag loci absent from established catalogues, annotate gene context, and collate candidate lists for replication.
5. **Trait-Trait Colocalisation (`R/Step5-Hyprcoloc_Trait-Trait.R`)**
   - Objective: Quantify whether pleiotropic signals across pairs of traits can be explained by shared causal variants.
   - Logic: Combine MTAG outputs with sentinel SNP lists, construct effect matrices, run HyPrColoc, and store posterior probabilities plus candidate SNP sets.
6. **Colocalisation Result Integration (`python/Step6_Merging_Trait-Trait_coloc.py`)**
   - Objective: Aggregate HyPrColoc outputs across all trait pairs to facilitate comparative interpretation and visualisation.
   - Logic: Merge per-pair CSV files, harmonise trait labels, and export consolidated summaries for plotting or reporting.
7. **MTAG vs eQTL Colocalisation (`R/Step7-Hyprcoloc_GWAS-eQTL.R` and `R/Step9-Hyprcoloc_3_GWAS-eQTLBrain.R`)**
   - Objective: Evaluate whether MTAG-identified variants colocalise with expression quantitative trait loci in disease-relevant brain tissues (e.g., MetaBrain, GTEx).
   - Logic: Load brain-region eQTL panels, harmonise alleles with MTAG hits, define genomic windows around target SNPs, and run HyPrColoc to quantify shared regulation.
8. **Replication with SMR (`docs/Step4-NovelLoci_2_replication_SMR.txt`)**
   - Objective: Provide a replicable workflow for validating candidate genes using SMR with 1000 Genomes LD references.
   - Logic: Convert eQTL panels to BESD format, run SMR against MTAG or GWAS summary statistics, and export results for comparative analysis.

Each script contains detailed parameter annotations; please consult the header comments before adapting to new datasets.

## Key Outputs

- Processed GWAS summary statistics in `Data/GWAS/Final_processed_data/`.
- MTAG meta-analysis results in `Files/MTAG_results/`.
- FUMA-derived locus annotations in `Files/FUMA_results/`.
- Trait-trait and GWAS-eQTL colocalisation summaries in `Files/Coloc_results/`.
- Optional SMR replication outputs under the SMR working directory described in the documentation.

## Acknowledgements

- **Version**: 1.0 
- **Last Updated**: 2025-11-10
- **Contributor**: Hexiao Ding (PolyU); Na Li (SCU)
- **Supervisor**: Dr. Jung Sun Yoo (PolyU)
- **Institute**: The Hong Kong Polytechnic University, Sichuan University
