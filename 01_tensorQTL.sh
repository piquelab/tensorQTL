#!/usr/bin/env bash
# =============================================================================
# 01_tensorQTL.sh
# Input preparation for tensorQTL cis-eQTL mapping
#
# Steps:
#   1. Convert filtered VCF to PLINK2 format (.pgen/.pvar/.psam)
#   2. Build the condition manifest (tensorqtl_conditions.tsv) from phenotype
#      BEDs and covariate files
#   3. Sort BED files (preserving header) if not already sorted
#   4. Submit the per-condition tensorQTL array job (02_tensorQTL.sbatch)
#
# Usage:
#   bash 01_tensorQTL.sh
#
# Expected inputs (set BASE path below):
#   - Filtered, biallelic SNP VCF (gzipped, no "chr" prefix on contigs)
#   - Per-condition phenotype BEDs:
#       GxP-eQTL_<CONDITION>_qnorm_genotypePC_COMBATNone_corrected.bed.gz
#   - Per-condition covariate files:
#       covariates/GxP-eQTL_<CONDITION>-SV1-15.covariates-COMBATNone-FastQTL.txt
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION — edit variables in this block before running
# ---------------------------------------------------------------------------

# Root directory containing phenotype BEDs and covariates/ subdirectory
BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping"

# Where PLINK2 files and the condition manifest will be written
TENSOR_DIR="${BASE}/tensor"

# Input VCF (filtered, biallelic SNPs, no 'chr' contig prefix)
VCF="${BASE}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz"

# Prefix for output PLINK2 files (.pgen / .pvar / .psam)
OUT_PREFIX="${TENSOR_DIR}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"

# Output manifest TSV (consumed by 02_tensorQTL.sbatch)
OUT_TSV="${TENSOR_DIR}/tensorqtl_conditions.tsv"

# ---------------------------------------------------------------------------
# NAMING PATTERNS — edit these to match YOUR file naming conventions
#
# BED_PATTERN   : glob passed to bash to find all per-condition phenotype BEDs.
#                 Must be in ${BASE}. The condition label will be extracted by
#                 stripping BED_PREFIX and BED_SUFFIX from each filename.
#
#   Example:  files named  study_<COND>_qnorm_corrected.bed.gz
#             BED_PREFIX="study_"
#             BED_SUFFIX="_qnorm_corrected.bed.gz"
#
# COV_PATTERN   : same as BED_PATTERN
#
#   Example:  covariates/study_<COND>_PCs.covariates.txt
#
# SORTED_SUFFIX : suffix appended to the base BED name for the sorted copy.
#                 Usually safe to leave as-is.
# ---------------------------------------------------------------------------
BED_PREFIX="GxP-eQTL_"
BED_SUFFIX="_qnorm_genotypePC_COMBATNone_corrected.bed.gz"
BED_PATTERN="${BASE}/${BED_PREFIX}*${BED_SUFFIX}"

SORTED_SUFFIX="_qnorm_genotypePC_COMBATNone_corrected.sorted.bed.gz"

COV_PREFIX="GxP-eQTL_"
COV_SUFFIX="-SV1-15.covariates-COMBATNone-FastQTL.txt"
COV_PATTERN="${BASE}/covariates/${COV_PREFIX}<COND>${COV_SUFFIX}"

mkdir -p "${TENSOR_DIR}"

# ===========================================================================
# STEP 1: VCF → PLINK2 format
# ===========================================================================
echo "[$(date)] Converting VCF to PLINK2 format..."

module load plink/2.0

plink2 \
  --vcf "${VCF}" \
  --max-alleles 2 \
  --set-missing-var-ids @:#:\$r:\$a \
  --chr 1-22 \
  --make-pgen \
  --out "${OUT_PREFIX}"

echo "[$(date)] PLINK2 files written to: ${OUT_PREFIX}.{pgen,pvar,psam}"

# ===========================================================================
# STEP 2: Build condition manifest
# ===========================================================================
echo "[$(date)] Building condition manifest: ${OUT_TSV}"

shopt -s nullglob
beds=( ${BED_PATTERN} )

if (( ${#beds[@]} == 0 )); then
  echo "ERROR: No phenotype BED files found matching: ${BED_PATTERN}" >&2
  echo "       Check BASE, BED_PREFIX, and BED_SUFFIX in the CONFIGURATION block." >&2
  exit 1
fi

echo -e "condition\tbed\tcovariates" > "${OUT_TSV}"

for bed in "${beds[@]}"; do
  fname=$(basename "${bed}")

  # Extract condition label by stripping prefix and suffix
  cond="${fname#${BED_PREFIX}}"
  cond="${cond%${BED_SUFFIX}}"

  # Path for sorted BED (required by tensorQTL)
  sorted_bed="${BASE}/${BED_PREFIX}${cond}${SORTED_SUFFIX}"

  # Sort if not already done
  if [[ ! -s "${sorted_bed}" ]]; then
    echo "  Sorting BED for ${cond}..."
    (zcat "${bed}" | head -1; zcat "${bed}" | tail -n +2 | sort -k1,1 -k2,2n) \
      | gzip > "${sorted_bed}"
  else
    echo "  Sorted BED already exists for ${cond}, skipping."
  fi

  # Resolve covariate file path (replace <COND> placeholder)
  cov="${COV_PATTERN//<COND>/${cond}}"
  if [[ ! -s "${cov}" ]]; then
    echo "WARNING: Missing covariate file for ${cond}: ${cov}" >&2
    continue
  fi

  echo -e "${cond}\t${sorted_bed}\t${cov}" >> "${OUT_TSV}"
done

# ===========================================================================
# STEP 3: Submit SLURM array job
# ===========================================================================
N=$(( $(wc -l < "${OUT_TSV}") - 1 ))  # subtract header line
echo "[$(date)] Manifest complete. Conditions found: ${N}"
echo "  Written to: ${OUT_TSV}"
wc -l "${OUT_TSV}"

if (( N == 0 )); then
  echo "ERROR: No valid conditions written to manifest. Aborting." >&2
  exit 1
fi

echo "[$(date)] Submitting SLURM array job (1-${N})..."
sbatch --array=1-"${N}" 02_tensorQTL.sbatch

echo "[$(date)] Done. Monitor jobs with: squeue -u \$USER"
