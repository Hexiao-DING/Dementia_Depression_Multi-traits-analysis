# Dementia_Depression_Multi-traits-analysis
This is a comperhensive R code for MTAG analysis in Dementia_Depression comorbidity. This repository consolidates the analytical scripts used to interrogate the genetic architecture shared by Alzheimer-related dementia phenotypes and depression-spectrum traits. The code base covers the complete workflow from genome-wide association study (GWAS) preprocessing to multi-trait meta-analysis, FUMA-driven functional annotation, and HyPrColoc-based colocalisation against eQTL data sets. 

## Table of Contents

- [Background](#background)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Initial Setup](#initial-setup)
- [Data Preparation](#data-preparation)
- [Workflow Summary](#workflow-summary)
- [Key Outputs](#key-outputs)
- [Reproducibility Guidelines](#reproducibility-guidelines)
- [Contributing](#contributing)
- [License](#license)

## Background

The project integrates multiple GWAS cohorts (dementia, vascular dementia, Lewy body dementia, depression, mixed anxiety and depression, and related traits) to characterise pleiotropic genetic loci. Multi-Trait Analysis of GWAS (MTAG) is used to boost statistical power, while FUMA provides downstream gene prioritisation and locus annotation. Colocalisation analyses are performed with HyPrColoc to evaluate shared genetic signals between MTAG-derived GWAS results and brain-region eQTL data sets.

## Repository Structure

- `R/` – Quality control, reference-genome liftover, and HyPrColoc analysis scripts written in R.
- `python/` – Python utilities for FUMA result aggregation, filtering, and MTAG-eQTL integration steps.
- `docs/` – Environment setup notes and user guides (MTAG installation on Linux, SMR replication instructions).
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

## Workflow Summary

1. **GWAS Quality Control (`R/Step1-*` scripts)**
   - Standardise column names, compute effect sizes, perform allele harmonisation.
   - Run genomic liftover between GRCh38 and GRCh37 as required.
   - Export MTAG-ready and LDSC-ready summary statistics.
2. **MTAG Analysis (`docs/Step2-MTAG_linux.txt`)**
   - Set up MTAG on a Linux environment.
   - Execute multi-trait meta-analyses per trait pairing.
3. **FUMA Aggregation (`python/Step3-*` scripts)**
   - Extract locus-level tables from FUMA result archives.
   - Merge trait-specific outputs and generate robust SNP sets.
4. **Novel Loci Definition (`python/Step4-NovelLoci_1_definition.py`)**
   - Integrate FUMA results with additional criteria to nominate novel loci.
5. **Trait-Trait Colocalisation (`R/Step5-Hyprcoloc_Trait-Trait.R`)**
   - Apply HyPrColoc to MTAG trait pairs to identify shared causal signals.
6. **MTAG-eQTL Colocalisation (`R/Step7-Hyprcoloc_GWAS-eQTL.R`, `R/Step9-Hyprcoloc_3_GWAS-eQTLBrain.R`)**
   - Load brain-region eQTL panels, match to MTAG results, and run HyPrColoc.
7. **Supporting Utilities**
   - `python/Step6_Merging_Trait-Trait_coloc.py` consolidates colocalisation summaries.
   - `docs/Step4-NovelLoci_2_replication_SMR.txt` outlines SMR replication of candidate loci.

Each script is self-contained; read the header comments for detailed parameter descriptions and file expectations.

## Key Outputs

- Processed GWAS summary statistics in `Data/GWAS/Final_processed_data/`.
- MTAG meta-analysis results in `Files/MTAG_results/`.
- FUMA-derived locus annotations in `Files/FUMA_results/`.
- Trait-trait and GWAS-eQTL colocalisation summaries in `Files/Coloc_results/`.
- Optional SMR replication outputs under the SMR working directory described in the documentation.

## Reproducibility Guidelines

- **Version control**: commit only scripts and configuration files. Large matrix or summary data sets should be excluded (consider `.gitignore` rules).
- **Metadata tracking**: capture software versions, LD reference panels, and eQTL releases used in each run within a lab notebook or an additional markdown log.
- **Random seeds**: set deterministic seeds where scripts perform random sampling or bootstrapping.
- **Documentation**: update the `docs/` folder with any deviations from the default pipeline to keep collaborators aligned.

## Acknowledgements

- **Version**: 1.0 
- **Last Updated**: 2025-11-10
- **Contributor**: Hexiao Ding (PolyU); Na Li (SCU)
- **Superviosr**: Dr. Jung Sun Yoo (PolyU)
- **Institute**: The Hong Kong Polytechnic University, Sichuan University
