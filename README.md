# Dementia–Depression Multi‑Trait Analysis
# GWAS → MTAG → FUMA → Colocalisation

> Robust, end‑to‑end pipeline for pleiotropy analysis across dementia and depression traits.

A complete, modular workflow to study pleiotropy between dementia phenotypes and depression-spectrum traits. The pipeline harmonises public GWAS summary statistics, runs multi-trait meta-analysis (MTAG), annotates loci with FUMA, and evaluates shared genetic signals using HyPrColoc against brain eQTL resources. It is designed for robustness, reproducibility, and selective execution of individual steps.

## Table of Contents

- [Overview](#overview)
- [Trait Catalogue and GWAS Sources](#trait-catalogue-and-gwas-sources)
- [Directory Layout](#directory-layout)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Data Preparation](#data-preparation)
- [Workflow Steps](#workflow-steps)
- [Outputs](#outputs)
- [Configuration Tips](#configuration-tips)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Citation and References](#citation-and-references)
- [Acknowledgements](#acknowledgements)

## Overview

- **Goal**: Identify and interpret shared genetic architecture for dementia and depression by integrating harmonised GWAS with MTAG, functional annotation (FUMA), and colocalisation (HyPrColoc with brain eQTL).
- **Design**: Each step is a standalone module. You can execute the full pipeline or run only the parts you need.
- **Reproducibility**: Clear environment specs, strict file conventions, and deterministic processing where applicable.

Flow at a glance:

```
GWAS sumstats  ──▶  QC + Harmonise (R)  ──▶  MTAG (Linux)  ──▶  FUMA annotate (web)
      │                    │                         │               │
      ▼                    ▼                         ▼               ▼
 LDSC inputs         Final_processed            MTAG_results     FUMA tables
      └───────────────────────────────────────────────────────────────────┐
                                                                          ▼
                                 HyPrColoc trait–trait (R) + GWAS–eQTL (R) ──▶ Coloc_results
```

## Trait Catalogue and GWAS Sources

This pipeline harmonises publicly available GWAS summary statistics stored under `Data/GWAS/Downloaded_original_data/`. Below are the core traits used in the dementia–depression analyses, with expected filenames, primary sources, and key references. Please respect data usage licenses for each provider.

Note: If your file names differ, update them consistently across steps (QC → MTAG → FUMA → colocalisation).

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

Where available, prefer the latest builds (GRCh37/GRCh38 as required by downstream tools) and harmonise alleles to a consistent reference. See the QC step for build liftover and allele orientation handling.

## Directory Layout

- `R/` – GWAS QC/harmonisation, trait–trait colocalisation, GWAS–eQTL colocalisation scripts.
- `python/` – FUMA post-processing, robust filtering, merging utilities.
- `docs/` – How-tos for MTAG on Linux and SMR replication instructions.
- `figs/` – Optional figures and plots.
- `Others(delete)/` – Legacy workspace from the source project (remove if not needed).

## Requirements

### R

- R (≥ 4.2)
- CRAN: `data.table`, `dplyr`, `tidyr`, `stringr`, `duckdb`, `DBI`, `httr`
- Bioconductor: `MungeSumstats`, `rtracklayer`, `SNPlocs.Hsapiens.dbSNP144.GRCh37`, `AnnotationHub`
- Colocalisation: `hyprcoloc`

### Python

- Python (≥ 3.8)
- Packages: `pandas`, `numpy`, `pyarrow` (optional but faster I/O)

### External Tools

- MTAG (multi-trait GWAS meta-analysis)
- FUMA (Functional Mapping and Annotation) – run on web, download outputs locally
- SMR (Summary-data-based Mendelian Randomization) – optional replication
- PLINK reference data (e.g., 1KG EUR) and brain eQTL panels (MetaBrain, GTEx)

## Quick Start

1) Clone and create environments

```bash
git clone <your-remote>
conda create -n gwas_r r-base=4.2 -y
conda create -n gwas_py python=3.10 -y
```

2) Install dependencies

```r
# R
install.packages(c("data.table","dplyr","tidyr","stringr","duckdb","DBI","httr"))
BiocManager::install(c("MungeSumstats","SNPlocs.Hsapiens.dbSNP144.GRCh37","AnnotationHub"))
install.packages("hyprcoloc", repos = "https://cloud.r-project.org")
```

```bash
# Python
conda activate gwas_py
pip install pandas numpy pyarrow
```

3) Read the notes in `docs/` (MTAG on Linux, SMR replication).

## Data Preparation

Scripts reference Windows-style roots under `D:/Projects_data&code/Stage1_bioinformatics_ADandDepression`. Either:

- Mirror the same hierarchy locally, or
- Edit path constants in scripts to your local data directories.

Required assets:

- Raw GWAS: `Data/GWAS/Downloaded_original_data/`
- Harmonised/derived: `Format_sumstats_data/`, `Final_processed_data/`, `MTAG_input_data/`, `LDSC_input_data/`
- MTAG outputs: `Files/MTAG_results/`
- FUMA downloads: `Files/FUMA_results/Raw_FUMA_results/`
- Colocalisation: `Files/Coloc_results/`
- Brain eQTL panels: `Data/eQTL/`

Ensure parent folders exist; most scripts do not create them automatically.

### Dataset readiness checklist

- [ ] All raw GWAS archives present under `Data/GWAS/Downloaded_original_data/`
- [ ] Reference build known (GRCh37/GRCh38) for each study
- [ ] Chain files available for liftover (if needed)
- [ ] eQTL panels downloaded and documented (MetaBrain/GTEx)
- [ ] PLINK reference specified (1KG EUR) for LD-based tools
- [ ] Disk space verified for intermediate outputs

## Workflow Steps

1) GWAS QC and Harmonisation (`R/Step1-QC_GWAS_*`)

- Standardises public GWAS to a uniform schema for MTAG/LDSC.
- Actions: allele checks, compute β/SE (from OR when needed), MAF imputation, missingness filters, GRCh38↔GRCh37 liftover.

Example (R):

```r
source("R/Step1-QC_GWAS_<trait>.R")
# configure input/output paths at top of the script
run_qc_for_trait("<trait_code>")
```

2) MTAG Execution (see `docs/Step2-MTAG_linux.txt`)

- Creates a Python 2.7-compatible MTAG environment on Linux.
- Aligns traits and launches MTAG with recommended parameters for dementia–depression analyses.

Example (Linux shell):

```bash
# after preparing harmonised inputs in MTAG format
python mtag.py \
  --sumstats trait1.sumstats.gz,trait2.sumstats.gz,... \
  --out Files/MTAG_results/mtag_joint \
  --n_min 50000 --force
```

3) FUMA Post-processing (`python/Step3-*`)

- `Step3-FUMA_1_Table_creating.py`: Extract per-locus tables from raw FUMA zips to harmonised CSV.
- `Step3-FUMA_2_Merging_by_trait.py`: Merge locus-level annotations across trait pairs.
- `Step3-FUMA_3_Robust_filter.py`: Concordance/robustness checks vs. original GWAS.
- `Step3-FUMA_4_Delete_failed.py`: Drop non-robust entries, preserve high-confidence loci.
- `Step3-FUMA_5_Unique_IndSNPs.py`: Deduplicate lead SNPs by trait combination.
- `Step3-FUMA_6.1_ExamSNPs_dementia.py`, `Step3-FUMA_6.2_ExamSNPs_depression.py`: Generate trait-specific inspection tables.

Examples (Python):

```bash
conda activate gwas_py
python python/Step3-FUMA_1_Table_creating.py --in Files/FUMA_results/Raw_FUMA_results --out Files/FUMA_results/tables
python python/Step3-FUMA_3_Robust_filter.py --mtag Files/MTAG_results --gwas Data/GWAS/Final_processed_data --out Files/FUMA_results/robust
```

4) Novel Locus Definition (`python/Step4-NovelLoci_1_definition.py`)

- Flags loci absent from major catalogues, annotates gene context, aggregates candidates.

Example:

```bash
python python/Step4-NovelLoci_1_definition.py --in Files/FUMA_results/robust --out Files/FUMA_results/novel
```

5) Trait–Trait Colocalisation (`R/Step5-Hyprcoloc_Trait-Trait.R`)

- Constructs effect matrices from MTAG and sentinel SNPs, runs HyPrColoc, saves posterior probabilities and credible sets.

Example (R):

```r
source("R/Step5-Hyprcoloc_Trait-Trait.R")
run_hyprcoloc_trait_trait(
  mtag_dir = "Files/MTAG_results",
  sentinels = "Files/FUMA_results/novel/sentinels.csv",
  out_dir = "Files/Coloc_results/trait_trait"
)
```

6) Coloc Result Integration (`python/Step6_Merging_Trait-Trait_coloc.py`)

- Merges per-pair HyPrColoc results into consolidated summary tables.

Example:

```bash
python python/Step6_Merging_Trait-Trait_coloc.py \
  --in Files/Coloc_results/trait_trait \
  --out Files/Coloc_results/trait_trait_merged.csv
```

7) GWAS–eQTL Colocalisation (`R/Step7-Hyprcoloc_GWAS-eQTL.R`, `R/Step9-Hyprcoloc_3_GWAS-eQTLBrain.R`)

- Tests MTAG loci for colocalisation with brain eQTL panels (MetaBrain/GTEx), using harmonised allele windows around targets.

Example (R):

```r
source("R/Step7-Hyprcoloc_GWAS-eQTL.R")
run_hyprcoloc_eqtl(
  mtag_hits = "Files/MTAG_results/mtag_joint.signals.csv",
  eqtl_root = "Data/eQTL",
  out_dir = "Files/Coloc_results/gwas_eqtl"
)
```

8) Optional Replication with SMR (`docs/Step4-NovelLoci_2_replication_SMR.txt`)

- Convert eQTL to BESD, run SMR against MTAG/GWAS, export for comparison.

## Outputs

- Harmonised GWAS: `Data/GWAS/Final_processed_data/`
- MTAG results: `Files/MTAG_results/`
- FUMA loci: `Files/FUMA_results/`
- Colocalisation summaries: `Files/Coloc_results/`
- Optional SMR outputs as configured in `docs/`

## Troubleshooting

- Missing columns during harmonisation: check delimiter, header case, required fields (SNP/CHR/BP/EA/NEA/BETA/SE/OR/P).
- Build mismatch (hg19 vs hg38): confirm chain files and `rtracklayer` liftover paths.
- MTAG alignment errors: verify allele orientation and sample-size metadata; ensure identical SNP IDs across traits.
- FUMA parsing: ensure raw FUMA directory structure is intact; unzip before running step 3.1.
- HyPrColoc failures: check MAF/effect alignment, region window sizes, and missingness filters.


## Citation and References

- MTAG: Turley P, et al. Multi-trait analysis of genome-wide association summary statistics using MTAG. Nat Genet. 2018.
- FUMA: Watanabe K, et al. Functional mapping and annotation of genetic associations. Nat Commun. 2017.
- HyPrColoc: Foley CN, et al. A fast and robust method for colocalisation. PLoS Genet. 2021.

## Acknowledgements

- Version: 1.1
- Last Updated: 2025-11-16
- Contributors: Hexiao Ding (PolyU); Na Li (SCU)
- Supervisor: Dr. Jung Sun Yoo (PolyU)
- Institutes: The Hong Kong Polytechnic University; Sichuan University

---

License: MIT (see `LICENSE` in the source repository).  
Contributions and issue reports are welcome. For questions, please open an issue or contact the contributors listed above.
