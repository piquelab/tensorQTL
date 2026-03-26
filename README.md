# eQTL Mapping Pipeline — tensorQTL

A full eQTL mapping pipeline from genotype VCF to significant eGene calls, using [tensorQTL](https://github.com/broadinstitute/tensorqtl) for GPU-accelerated cis-QTL mapping.

Originally developed for the GxP project (GxP_20250730), but structured to be portable to any study with a similar design: multiple treatment/timepoint conditions, RNA-seq phenotypes, and a shared genotype panel.

---

## Overview

```
VCF (filtered SNPs)
      │
      ▼
  PLINK2 files (.pgen/.pvar/.psam)          [00 / 01]
      │
      ▼
  Sorted phenotype BEDs + covariate files   [01]
      │
      ▼
  tensorQTL cis_nominal (SLURM array)       [02]
      │
      ├─▶  .parquet → .txt conversion       [03]  ← use for mashr / downstream
      │
      └─▶  tensorQTL cis (permutation)      [02,  ← use for calling eGenes]
                │
                ▼
           FDR correction + eGene summary   [04]
```

---

## Scripts

| # | Script | Description |
|---|--------|-------------|
| 00 | `00_setup_tensorQTL.sh` | One-time conda environment setup (rpy2, R qvalue) |
| 01 | `01_tensorQTL.sh` | VCF → PLINK2; build condition manifest; submit SLURM array |
| 02 | `02_tensorQTL.sbatch` | SLURM array job — one condition per task, run both cis and cis_nominal mode |
| 03 | `03_convert_parquet_to_txt.sh` | Convert cis_nominal parquet output to tab-delimited text |
| 04 | `04_tQTL_fdr_analysis.R` | Load cis (permuted) results; apply FDR; plots + summary |

---

## Requirements

### Software

| Tool | Version | Notes |
|------|---------|-------|
| conda / mamba | any | environment management |
| PLINK2 | 2.0 | `module load plink/2.0` on cluster |
| tensorQTL | ≥ 0.6 | installed in conda env (see setup) |
| Python | 3.11 | via conda env |
| R | ≥ 4.1 | for FDR analysis script |
| pandas + pyarrow | any | parquet conversion |
| rpy2 | ≤ 3.5.12 | Python–R bridge for qvalue |

### R packages

```r
# CRAN
install.packages(c("data.table", "ggplot2", "dplyr", "tidyr", "viridis"))

# Bioconductor
BiocManager::install("qvalue")
```

### Input files

| File | Description |
|------|-------------|
| `*.vcf.gz` | Filtered, biallelic SNP VCF (no `chr` prefix on contigs) |
| `prefix_<COND>_suffix.bed.gz` | Quantile-normalised, batch-corrected phenotype BED per condition |
| `covariates/prefix_<COND>_suffix.txt` | Covariate matrix per condition (genotype PCs + PEER/SVs) |

---

## Setup (one-time)
Recommend running interactively.

```bash
bash scripts/00_setup_tensorQTL.sh
```

This will:
1. Create the `tensorqtl_p3.11_env` conda environment from a base YAML
2. Install `rpy2 ≤ 3.5.12` via pip
3. Install `r-essentials` via conda
4. Install the Bioconductor `qvalue` package inside R
5. Verify the full Python → rpy2 → R → qvalue chain

After setup, activate the environment before running any other scripts:

```bash
conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
```

> **Note:** The `LD_LIBRARY_PATH` line is is for the pique group R/4.3.2 which is known to be compatible with all analyses. 

---

## Usage

### Step 1 — Prepare inputs and submit jobs

```bash
bash scripts/01_tensorQTL.sh
```

What it does:
- Converts the input VCF to PLINK2 format (autosomes, biallelic only)
- Discovers all condition-level phenotype BEDs matching the naming pattern
- Creates sorted BEDs if they don't already exist
- Validates that covariate files exist for each condition
- Writes `tensorqtl_conditions.tsv` (the per-job manifest)
- Submits `02_tensorQTL.sbatch` as a SLURM array sized to the number of conditions

**Edit before running:** Set the paths and **naming patterns** at the top of the script:

| Variable | Description | Example |
|----------|-------------|---------|
| `BASE` | Root directory containing BEDs and `covariates/` | `/path/to/eQTL_mapping` |
| `VCF` | Filtered biallelic SNP VCF (no `chr` prefix) | `study_genotypes.vcf.gz` |
| `BED_PREFIX` | Filename prefix before the condition label | `"study_"` |
| `BED_SUFFIX` | Filename suffix after the condition label | `"_qnorm_corrected.bed.gz"` |
| `COV_PREFIX` | Filename prefix before the condition label | `"study_"` |
| `COV_SUFFIX` | Filename suffix after the condition label | `"_10-SVs.txt"` |
| `SORTED_SUFFIX` | Suffix for sorted BED copies (usually matches `BED_SUFFIX` with `.sorted` inserted) | `"_qnorm_corrected.sorted.bed.gz"` |

The condition label is extracted automatically by stripping `BED_PREFIX` and `BED_SUFFIX` from each BED filename.

---

### Step 2 — Submit tensorQTL in cis_nominal (all SNP-gene pairs) and cis (permutations) modes via SLURM

`02_tensorQTL.sbatch` runs automatically as a SLURM array from Step 1.  
Each task maps to one row in `tensorqtl_conditions.tsv` and runs:

```
python3 -m tensorqtl <plink_prefix> <pheno_bed> <out_prefix> \
  --covariates <cov_file> \
  --mode cis_nominal \
  --fdr 0.1 \
  --window 100000
```
AND
```
python3 -m tensorqtl <plink_prefix> <pheno_bed> <out_prefix> \
  --covariates <cov_file> \
  --mode cis \
  --fdr 0.1 \
  --window 100000
```

**Key parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| `--mode` | `cis_nominal` | All SNP–gene pairs within window; no permutation |
| `--window` | 100,000 bp | ±100 kb from TSS |
| `--fdr` | 0.1 | Passed to tensorQTL's internal FDR step |

**Outputs** (one set per condition in `tensorqtl_output_cis-nominal_SV15_100kb/`):
```
<CONDITION>_cis.cis_nominal_pairs.<chunk>.parquet
<CONDITION>_cis.cis_qtl_thresholds.txt.gz
```

To monitor jobs:
```bash
squeue -u $USER
```

To re-run a single condition manually (e.g., for debugging):
```bash
SLURM_ARRAY_TASK_ID=3 bash scripts/02_tensorQTL.sbatch
```

---

### Step 3 — Convert parquet to text (for mashr / downstream)

```bash
bash scripts/03_convert_parquet_to_txt.sh
```

Converts all `.parquet` files in the cis_nominal output directory to tab-delimited `.txt` files. Output columns include `phenotype_id`, `variant_id`, `tss_distance`, `pval_nominal`, `slope`, `slope_se`.

---

### Step 4 — FDR correction and eGene summary (cis permutation results)

```bash
/path/to/Rscript scripts/04_tQTL_fdr_analysis.R
```

This script operates on **cis permutation** output (not cis_nominal). It:
1. Loads `<CONDITION>_cis.cis_qtl.txt.gz` for each condition
2. Selects the best available p-value (`pval_beta` > `pval_perm` > `pval_nominal`)
3. Applies cross-gene FDR via `p.adjust()`
4. Counts eGenes at multiple thresholds (FDR 0.05/0.10/0.20, Bonferroni, nominal)
5. Writes QQ plots, method comparison bar charts, and a plain-text summary

**Outputs** (in `tensorqtl_output_cis_SV15_100kb/analysis_results_fdr/`):
```
eQTL_summary_fdr.txt                    — per-condition eGene counts
FDR_ANALYSIS_SUMMARY.txt                — overall summary
figures/QQ_fdr_<CONDITION>.png
figures/eGenes_method_comparison.png
figures/pvalue_distributions_<COND>.png
```
A list of unique significant eGenes (FDR < 0.10, any condition) is written to:
```
tensorqtl_output_cis_SV15_100kb/unique_significant_egenes.txt
```

---

## Condition manifest format

`tensorqtl_conditions.tsv` — generated by `01_tensorQTL.sh`:

```
condition       bed                                          covariates
BPA_100nM_24    /path/to/GxP-eQTL_BPA_100nM_24_...sorted.bed.gz    /path/to/covariates/...txt
BPA_100nM_6     /path/to/GxP-eQTL_BPA_100nM_6_...sorted.bed.gz     /path/to/covariates/...txt
...
```

SLURM array task IDs map to rows: task 1 = first data row, task 2 = second, etc.

---

## Notes on FDR strategy

tensorQTL's cis permutation mode corrects for multiple testing **within each gene** using a beta-approximation to the permutation p-value (`pval_beta`). The FDR script applies a **second-pass FDR** across all genes within a condition — this is an additional step and deliberately conservative.

| Step | What it corrects for |
|------|----------------------|
| tensorQTL permutation | Multiple SNPs per gene |
| `p.adjust(pval_beta, "fdr")` | Multiple genes per condition |

For exploratory comparisons across conditions, `pval_nominal < 0.001` or FDR < 0.10 are reasonable thresholds; for a final significant eGene list, use FDR < 0.05.

---

## Directory structure

```
eQTL_mapping/
├── tensor/
│   ├── GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.{pgen,pvar,psam}
│   ├── tensorqtl_conditions.tsv
│   ├── tensorqtl_output_cis-nominal_SV15_100kb/
│   │   ├── <CONDITION>_cis.cis_nominal_pairs.*.parquet
│   │   └── txt_files/
│   │       └── <CONDITION>_cis.cis_nominal_pairs.*.txt
│   └── tensorqtl_output_cis_SV15_100kb/
│       ├── <CONDITION>_cis.cis_qtl.txt.gz
│       └── analysis_results_fdr/
│           ├── eQTL_summary_fdr.txt
│           ├── FDR_ANALYSIS_SUMMARY.txt
│           └── figures/
├── GxP-eQTL_<COND>_qnorm_genotypePC_COMBATNone_corrected.bed.gz
├── GxP-eQTL_<COND>_qnorm_genotypePC_COMBATNone_corrected.sorted.bed.gz
└── covariates/
    └── GxP-eQTL_<COND>-SV1-15.covariates-COMBATNone-FastQTL.txt
```

---

## Reference

Taylor-Weiner A, et al. (2019). Scaling computational genomics to millions of individuals with GPUs. *Genome Biology*, 20, 228. https://doi.org/10.1186/s13059-019-1836-7
