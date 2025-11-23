<div align="center">

# 🧬 Dementia–Depression Multi‑Trait Analysis Pipeline

**A comprehensive, end-to-end workflow for pleiotropy analysis across dementia and depression traits**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%3E%3D%203.8-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20Windows-lightgrey)]()

[Documentation](#-overview) • [Quick Start](#-quick-start) • [Workflow](#-workflow-steps) • [Citation](#-citation-and-references)

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Trait Catalogue](#-trait-catalogue-and-gwas-sources)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Data Preparation](#-data-preparation)
- [Workflow Steps](#-workflow-steps)
- [Outputs](#-outputs)
- [Troubleshooting](#-troubleshooting)
- [FAQ](#-faq)
- [Citation and References](#-citation-and-references)
- [Acknowledgements](#-acknowledgements)

---

## 🎯 Overview

This pipeline provides a **robust, modular workflow** to study pleiotropy between dementia phenotypes and depression-spectrum traits. It integrates harmonised GWAS summary statistics with multi-trait meta-analysis (MTAG), functional annotation (FUMA), and colocalisation analysis (HyPrColoc) against brain eQTL resources.

### Key Objectives

- 🔍 **Identify** shared genetic architecture between dementia and depression
- 📊 **Harmonise** public GWAS summary statistics to a uniform schema
- 🔬 **Annotate** loci with functional mapping and gene context
- 🎯 **Evaluate** shared genetic signals using colocalisation methods
- ✅ **Ensure** reproducibility through clear environment specifications and deterministic processing

### Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Pipeline Workflow                                   │
└─────────────────────────────────────────────────────────────────────────────┘

  Raw GWAS Data
       │
       ▼
  ┌────────────────────────────────────────────────────────────┐
  │  Step 1: Quality Control & Harmonisation (R)               │
  │  • Allele checks & orientation                             │
  │  • β/SE computation from OR                                │
  │  • MAF imputation                                          │
  │  • Build liftover (GRCh37 ↔ GRCh38)                        │
  └────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────┐
  │  Step 2: Multi-Trait Meta-Analysis (MTAG, Linux)           │
  │  • Trait alignment                                         │
  │  • Joint effect estimation                                 │
  └────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────┐
  │  Step 3: Functional Annotation (FUMA, Web + Python)        │
  │  • Locus identification                                    │
  │  • Gene mapping & annotation                               │
  │  • Robustness filtering                                    │
  └────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
  ┌────────────────────────────────────────────────────────────┐
  │  Steps 4-7: Colocalisation Analysis (R)                    │
  │  • Trait-trait colocalisation (HyPrColoc)                  │
  │  • GWAS-eQTL colocalisation (MetaBrain, GTEx)              │
  │  • Result integration & visualization                      │
  └────────────────────────────────────────────────────────────┘
```

---

## ✨ Features

- 🔄 **Modular Design**: Execute individual steps or run the complete pipeline
- 🛡️ **Robust QC**: Comprehensive quality control and harmonisation procedures
- 📈 **Multi-Trait Analysis**: Leverage MTAG for improved power and precision
- 🧬 **Functional Annotation**: Integrate FUMA for gene mapping and context
- 🎯 **Colocalisation**: Identify shared causal variants across traits and eQTLs
- 🔬 **Reproducible**: Version-controlled environments and deterministic processing
- 📚 **Well-Documented**: Comprehensive documentation and examples

---

## 📊 Trait Catalogue and GWAS Sources

This pipeline harmonises publicly available GWAS summary statistics. All data should be placed under `Data/GWAS/Downloaded_original_data/`. Please **respect data usage licenses** for each provider.

> **⚠️ Important**: If your file names differ, update them consistently across all steps (QC → MTAG → FUMA → colocalisation).

### Dementia-Related Traits (7 outcomes)

| Trait | Summary Statistics | Source | Primary Citation |
|:------|:-------------------|:-------|:-----------------|
| **Dementia**<br>(FinnGen F5_DEMENTIA) | `Finn-b-F5_DEMENTIA.tsv.gz` | [FinnGen Risteys](https://r10.risteys.finngen.fi/endpoints/F5_DEMENTIA) | Kurki *et al.* (2023) *Nature* **613**, 508-518 |
| **Alzheimer's disease** | `GCST90012877.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589840) | Schwartzentruber *et al.* (2021) *Nat. Genet.* **53**, 392-402 |
| **Cognitive performance** | `GCST006572.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/30038396) | Lee *et al.* (2018) *Nat. Genet.* **50**, 1112-1121 |
| **Vascular dementia**<br>(FinnGen F5_VASCDEM) | `Finn-b-F5_VASCDEM.tsv.gz` | [FinnGen Risteys](https://r9.risteys.finngen.fi/endpoints/F5_VASCDEM) | Kurki *et al.* (2023) *Nature* **613**, 508-518 |
| **Lewy body dementia** | `GCST90001390.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/33589841) | Chia *et al.* (2021) *Nat. Genet.* **53**, 294-303 |
| **Frontotemporal dementia** | `GCST90558311.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/40280976) | Pottier *et al.* (2025) *Nat. Commun.* **16**, 3914 |
| **Undefined dementia**<br>(FinnGen F5_DEMNAS) | `Finn_b_F5_DEMNAS.tsv.gz` | [FinnGen Risteys](https://r12.risteys.finngen.fi/endpoints/F5_DEMNAS) | Kurki *et al.* (2023) *Nature* **613**, 508-518 |

### Depression-Spectrum Traits (3 outcomes)

| Trait | Summary Statistics | Source | Primary Citation |
|:------|:-------------------|:-------|:-----------------|
| **Depressive disorders** | `GCST90476922.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39024449) | Verma *et al.* (2024) *Science* **385**, eadj1182 |
| **Major depressive disorder** | `GCST90468123.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/39789286) | Loya *et al.* (2025) *Nat. Genet.* **57**, 461-468 |
| **Mixed anxiety and depressive disorder** | `GCST90225526.tsv.gz` | [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/37259642) | Brasher *et al.* (2023) *Genes Brain Behav.* **22**, e12851 |

### Genome Build Notes

The pipeline supports both **GRCh37** and **GRCh38** with automatic liftover capabilities. Prefer the latest builds as required by downstream tools and harmonise alleles to a consistent reference. See the QC step for detailed build liftover and allele orientation handling.

---

## 💻 Installation

### Prerequisites

- **R** ≥ 4.2
- **Python** ≥ 3.8
- **Conda** (recommended for environment management)
- **Linux** environment (required for MTAG execution)

### Step 1: Clone Repository

```bash
git clone https://github.com/Hexiao-DING/Dementia_Depression_Multi-traits-analysis.git
cd Dementia_Depression_Multi-traits-analysis
```

### Step 2: Create Conda Environments

```bash
# Create R environment
conda create -n gwas_r r-base=4.2 -y

# Create Python environment
conda create -n gwas_py python=3.10 -y
```

### Step 3: Install R Dependencies

```bash
conda activate gwas_r
```

```r
# Install CRAN packages
install.packages(c(
  "data.table", "dplyr", "tidyr", "stringr",
  "duckdb", "DBI", "httr"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "MungeSumstats",
  "rtracklayer",
  "SNPlocs.Hsapiens.dbSNP144.GRCh37",
  "AnnotationHub"
))

# Install HyPrColoc
install.packages("hyprcoloc", repos = "https://cloud.r-project.org")
```

### Step 4: Install Python Dependencies

```bash
conda activate gwas_py
pip install pandas numpy pyarrow
```

### Step 5: External Tools

- **MTAG**: Follow instructions in `docs/Step2-MTAG_linux.txt` for Linux setup
- **FUMA**: Web-based tool ([fuma.ctglab.nl](https://fuma.ctglab.nl/))
- **SMR**: Optional replication tool (see `docs/Step4-NovelLoci_2_replication_SMR.txt`)
- **PLINK**: Reference data (1000 Genomes EUR) for LD-based tools
- **eQTL Panels**: MetaBrain and GTEx brain eQTL data

---

## 🚀 Quick Start

### Minimal Example

```bash
# 1. Activate R environment
conda activate gwas_r

# 2. Run QC for a single trait
Rscript R/Step1-QC_GWAS_dementia.R

# 3. Activate Python environment
conda activate gwas_py

# 4. Process FUMA results
python python/Step3-FUMA_1_Table_creating.py \
  --in Files/FUMA_results/Raw_FUMA_results \
  --out Files/FUMA_results/tables
```

### Full Pipeline

See [Workflow Steps](#-workflow-steps) for detailed instructions on running the complete analysis pipeline.

---

## 📁 Data Preparation

### Directory Structure

Scripts reference paths under `D:/Projects_data&code/Stage1_bioinformatics_ADandDepression`. You can either:

1. **Mirror the hierarchy locally**, or
2. **Edit path constants** in scripts to match your local directories

### Required Data Directories

```
Data/
├── GWAS/
│   ├── Downloaded_original_data/      # Raw GWAS summary statistics
│   ├── Format_sumstats_data/           # Intermediate harmonised data
│   ├── Final_processed_data/           # Final QC'd GWAS
│   ├── MTAG_input_data/                # MTAG-formatted inputs
│   └── LDSC_input_data/                # LDSC-formatted inputs
├── eQTL/                               # Brain eQTL panels (MetaBrain, GTEx)
└── Files/
    ├── MTAG_results/                   # MTAG outputs
    ├── FUMA_results/
    │   ├── Raw_FUMA_results/           # Downloaded FUMA outputs
    │   ├── tables/                     # Extracted tables
    │   ├── robust/                     # Robust filtered results
    │   └── novel/                      # Novel loci
    └── Coloc_results/                  # Colocalisation results
```

### Project Structure

```
.
├── R/                                  # R scripts
│   ├── Step1-QC_GWAS_*.R              # GWAS QC and harmonisation
│   ├── Step5-Hyprcoloc_Trait-Trait.R  # Trait-trait colocalisation
│   ├── Step7-Hyprcoloc_GWAS-eQTL.R    # GWAS-eQTL colocalisation
│   └── Step9-Hyprcoloc_3_GWAS-eQTLBrain.R
├── python/                             # Python scripts
│   ├── Step3-FUMA_*.py                 # FUMA post-processing
│   ├── Step4-NovelLoci_*.py           # Novel locus identification
│   └── Step6_Merging_*.py             # Result integration
├── docs/                               # Documentation
│   ├── Step2-MTAG_linux.txt           # MTAG setup guide
│   └── Step4-NovelLoci_2_replication_SMR.txt
└── figs/                               # Figures and plots
```

> **Note**: Ensure parent folders exist before running scripts; most scripts do not create directories automatically.

### Dataset Readiness Checklist

Before starting the pipeline, verify:

- [ ] All raw GWAS archives present under `Data/GWAS/Downloaded_original_data/`
- [ ] Reference build known (GRCh37/GRCh38) for each study
- [ ] Chain files available for liftover (if needed)
- [ ] eQTL panels downloaded and documented (MetaBrain/GTEx)
- [ ] PLINK reference specified (1000 Genomes EUR) for LD-based tools
- [ ] Disk space verified for intermediate outputs (recommend ≥100 GB free)

---

## 🔬 Workflow Steps

### Step 1: GWAS QC and Harmonisation

**Scripts**: `R/Step1-QC_GWAS_*.R`

**Purpose**: Standardises public GWAS to a uniform schema for MTAG/LDSC.

**Key Actions**:
- ✅ Allele checks and orientation
- ✅ Compute β/SE from OR when needed
- ✅ MAF imputation
- ✅ Missingness filters
- ✅ GRCh38↔GRCh37 liftover

**Example**:

```r
source("R/Step1-QC_GWAS_dementia.R")
# Configure input/output paths at top of the script
run_qc_for_trait("dementia")
```

---

### Step 2: MTAG Execution

**Documentation**: `docs/Step2-MTAG_linux.txt`

**Purpose**: Multi-trait meta-analysis to identify shared genetic signals.

**Requirements**: 
- Linux environment
- Python 2.7-compatible MTAG installation

**Example**:

```bash
# After preparing harmonised inputs in MTAG format
python mtag.py \
  --sumstats trait1.sumstats.gz,trait2.sumstats.gz,... \
  --out Files/MTAG_results/mtag_joint \
  --n_min 50000 \
  --force
```

---

### Step 3: FUMA Post-Processing

**Scripts**: `python/Step3-FUMA_*.py`

**Purpose**: Extract, merge, and filter FUMA annotation results.

**Script Overview**:

| Script | Purpose |
|:-------|:--------|
| `Step3-FUMA_1_Table_creating.py` | Extract per-locus tables from raw FUMA zips to harmonised CSV |
| `Step3-FUMA_2_Merging_by_trait.py` | Merge locus-level annotations across trait pairs |
| `Step3-FUMA_3_Robust_filter.py` | Concordance/robustness checks vs. original GWAS |
| `Step3-FUMA_4_Delete_failed.py` | Drop non-robust entries, preserve high-confidence loci |
| `Step3-FUMA_5_Unique_IndSNPs.py` | Deduplicate lead SNPs by trait combination |
| `Step3-FUMA_6.1_ExamSNPs_dementia.py` | Generate dementia-specific inspection tables |
| `Step3-FUMA_6.2_ExamSNPs_depression.py` | Generate depression-specific inspection tables |

**Examples**:

```bash
conda activate gwas_py

# Extract FUMA tables
python python/Step3-FUMA_1_Table_creating.py \
  --in Files/FUMA_results/Raw_FUMA_results \
  --out Files/FUMA_results/tables

# Robust filtering
python python/Step3-FUMA_3_Robust_filter.py \
  --mtag Files/MTAG_results \
  --gwas Data/GWAS/Final_processed_data \
  --out Files/FUMA_results/robust
```

---

### Step 4: Novel Locus Definition

**Script**: `python/Step4-NovelLoci_1_definition.py`

**Purpose**: Flags loci absent from major catalogues, annotates gene context, aggregates candidates.

**Example**:

```bash
python python/Step4-NovelLoci_1_definition.py \
  --in Files/FUMA_results/robust \
  --out Files/FUMA_results/novel
```

---

### Step 5: Trait–Trait Colocalisation

**Script**: `R/Step5-Hyprcoloc_Trait-Trait.R`

**Purpose**: Constructs effect matrices from MTAG and sentinel SNPs, runs HyPrColoc, saves posterior probabilities and credible sets.

**Example**:

```r
source("R/Step5-Hyprcoloc_Trait-Trait.R")
run_hyprcoloc_trait_trait(
  mtag_dir = "Files/MTAG_results",
  sentinels = "Files/FUMA_results/novel/sentinels.csv",
  out_dir = "Files/Coloc_results/trait_trait"
)
```

---

### Step 6: Coloc Result Integration

**Script**: `python/Step6_Merging_Trait-Trait_coloc.py`

**Purpose**: Merges per-pair HyPrColoc results into consolidated summary tables.

**Example**:

```bash
python python/Step6_Merging_Trait-Trait_coloc.py \
  --in Files/Coloc_results/trait_trait \
  --out Files/Coloc_results/trait_trait_merged.csv
```

---

### Step 7: GWAS–eQTL Colocalisation

**Scripts**: `R/Step7-Hyprcoloc_GWAS-eQTL.R`, `R/Step9-Hyprcoloc_3_GWAS-eQTLBrain.R`

**Purpose**: Tests MTAG loci for colocalisation with brain eQTL panels (MetaBrain/GTEx), using harmonised allele windows around targets.

**Example**:

```r
source("R/Step7-Hyprcoloc_GWAS-eQTL.R")
run_hyprcoloc_eqtl(
  mtag_hits = "Files/MTAG_results/mtag_joint.signals.csv",
  eqtl_root = "Data/eQTL",
  out_dir = "Files/Coloc_results/gwas_eqtl"
)
```

---

### Step 8: Optional Replication with SMR

**Documentation**: `docs/Step4-NovelLoci_2_replication_SMR.txt`

**Purpose**: Convert eQTL to BESD format, run SMR against MTAG/GWAS, export for comparison.

---

## 📤 Outputs

The pipeline generates outputs in the following directories:

| Output Type | Location | Description |
|:------------|:---------|:------------|
| **Harmonised GWAS** | `Data/GWAS/Final_processed_data/` | QC'd and harmonised summary statistics |
| **MTAG Results** | `Files/MTAG_results/` | Multi-trait meta-analysis outputs |
| **FUMA Loci** | `Files/FUMA_results/` | Functional annotation results |
| **Colocalisation** | `Files/Coloc_results/` | Trait-trait and GWAS-eQTL colocalisation summaries |
| **SMR Outputs** | As configured in `docs/` | Optional replication results |

---

## 🔧 Troubleshooting

### Common Issues and Solutions

| Issue | Solution |
|:------|:---------|
| **Missing columns during harmonisation** | Check delimiter, header case, required fields (SNP/CHR/BP/EA/NEA/BETA/SE/OR/P) |
| **Build mismatch (hg19 vs hg38)** | Confirm chain files and `rtracklayer` liftover paths are correct |
| **MTAG alignment errors** | Verify allele orientation and sample-size metadata; ensure identical SNP IDs across traits |
| **FUMA parsing errors** | Ensure raw FUMA directory structure is intact; unzip files before running step 3.1 |
| **HyPrColoc failures** | Check MAF/effect alignment, region window sizes, and missingness filters |

---

## ❓ FAQ

### Do I need Linux for MTAG?

Yes, MTAG requires a Linux environment. Follow `docs/Step2-MTAG_linux.txt` to create a Python 2.7-compatible MTAG environment on Linux.

### Which genome build should I use?

The pipeline supports both GRCh37/GRCh38 with liftover capabilities. Keep a consistent build across steps or perform liftover during the QC step.

### Can I run only FUMA and colocalisation?

Yes. Provide harmonised GWAS/MTAG inputs in the expected formats and start from step 3 (FUMA post-processing) or step 5 (colocalisation).

### Where do I put FUMA downloads?

Place downloaded FUMA outputs under `Files/FUMA_results/Raw_FUMA_results/` before running step 3 scripts. Ensure the directory structure matches FUMA's output format.

### How long does the full pipeline take?

Processing time varies by dataset size and computational resources. Typical runtime:
- **Step 1 (QC)**: 1-2 hours per trait
- **Step 2 (MTAG)**: 2-4 hours for all traits
- **Step 3 (FUMA post-processing)**: 30 minutes - 1 hour
- **Steps 4-7 (Colocalisation)**: 1-3 hours depending on number of loci

---

## 📚 Citation and References

### Source Repository

- **Repository**: [https://github.com/Hexiao-DING/Dementia_Depression_Multi-traits-analysis](https://github.com/Hexiao-DING/Dementia_Depression_Multi-traits-analysis)

### Key Methodological References

**MTAG**  
Turley P, Walters RK, Maghzian O, *et al*. Multi-trait analysis of genome-wide association summary statistics using MTAG. *Nat Genet*. 2018;50:229-237.  
[DOI: 10.1038/s41588-017-0009-4](https://doi.org/10.1038/s41588-017-0009-4)

**FUMA**  
Watanabe K, Taskesen E, van Bochoven A, *et al*. Functional mapping and annotation of genetic associations with FUMA. *Nat Commun*. 2017;8:1826.  
[DOI: 10.1038/s41467-017-01261-5](https://doi.org/10.1038/s41467-017-01261-5)

**HyPrColoc**  
Foley CN, Staley JR, Breen PG, *et al*. A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. *PLoS Genet*. 2021;17:e1009440.  
[DOI: 10.1371/journal.pgen.1009440](https://doi.org/10.1371/journal.pgen.1009440)

---


## 🙏 Acknowledgements

| Role | Name | Affiliation |
|:-----|:-----|:------------|
| **Contributors** | Hexiao Ding | The Hong Kong Polytechnic University |
| **Contributors** | Na Li | Sichuan University |
| **Supervisor** | Dr. Jung Sun Yoo | The Hong Kong Polytechnic University |

**Institutions**:
- The Hong Kong Polytechnic University
- Sichuan University

**Version**: 1.2  
**Last Updated**: 2025-11-23

