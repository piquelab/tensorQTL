#!/bin/bash
# =============================================================================
# tensorQTL Pre-processing Script
# Configure the variables in the USER CONFIGURATION section below, then run.
# =============================================================================

# =============================================================================
# USER CONFIGURATION
# =============================================================================

### Base directory
BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping"

### Genotype VCF file (full path)
VCF="${BASE}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr.vcf.gz"

### Output prefix for PLINK2 files (no extension)
OUT_PREFIX="${BASE}/tensor/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"

### BED file pattern
# Use a glob pattern that matches all your phenotype BED files.
# The script will extract the condition name from the filename using
# BED_PREFIX and BED_SUFFIX (the parts before and after the condition name).
BED_DIR="${BASE}"
BED_PREFIX="GxP-eQTL_"
BED_SUFFIX="_qnorm_genotypePC_COMBATNone_corrected.bed.gz"

### Covariate file pattern
# The condition name will be inserted between COV_PREFIX and COV_SUFFIX.
COV_DIR="${BASE}/covariates"
COV_PREFIX="GxP-eQTL_"
COV_SUFFIX="-SV1-15.covariates-COMBATNone-FastQTL.txt"

### Output TSV (condition table)
OUT_TSV="${BASE}/tensor/tensorqtl_conditions.tsv"

### Downstream sbatch script to submit
SBATCH_SCRIPT="02_tensorqtl.sbatch"

### Chromosomes to include in PLINK2 (e.g. "1-22" for autosomes)
CHROMOSOMES="1-22"

# =============================================================================
# END OF USER CONFIGURATION — do not edit below unless you know what you're doing
# =============================================================================

set -euo pipefail

### 1 — Set up output directory
mkdir -p "${BASE}/tensor"

### 2 — VCF to PLINK2
echo "=== Step 1: Converting VCF to PLINK2 format ==="
module load plink/2.0

plink2 \
  --vcf "${VCF}" \
  --max-alleles 2 \
  --set-missing-var-ids @:#:\$r:\$a \
  --chr "${CHROMOSOMES}" \
  --make-pgen \
  --out "${OUT_PREFIX}"

echo "Done. PLINK2 prefix: ${OUT_PREFIX}"

### 3 — Build condition table
echo ""
echo "=== Step 2: Building condition table ==="

shopt -s nullglob
beds=( "${BED_DIR}/${BED_PREFIX}"*"${BED_SUFFIX}" )

if (( ${#beds[@]} == 0 )); then
  echo "ERROR: No phenotype BEDs found matching pattern:" >&2
  echo "       ${BED_DIR}/${BED_PREFIX}*${BED_SUFFIX}" >&2
  exit 1
fi

echo -e "condition\tbed\tcovariates" > "${OUT_TSV}"

for bed in "${beds[@]}"; do
  fname=$(basename "${bed}")

  # Extract condition name by stripping prefix and suffix
  cond="${fname#${BED_PREFIX}}"
  cond="${cond%${BED_SUFFIX}}"

  # Define sorted BED path
  sorted_bed="${BED_DIR}/${BED_PREFIX}${cond}${BED_SUFFIX%.bed.gz}.sorted.bed.gz"

  # Create sorted BED if it doesn't already exist
  if [[ ! -s "${sorted_bed}" ]]; then
    echo "  Creating sorted BED for: ${cond}"
    (zcat "${bed}" | head -1; zcat "${bed}" | tail -n +2 | sort -k1,1 -k2,2n) \
      | gzip > "${sorted_bed}"
  else
    echo "  Sorted BED already exists for: ${cond}"
  fi

  # Define covariate file path
  cov="${COV_DIR}/${COV_PREFIX}${cond}${COV_SUFFIX}"

  if [[ ! -s "${cov}" ]]; then
    echo "  WARNING: Missing covariates for ${cond} — skipping" >&2
    echo "           Expected: ${cov}" >&2
    continue
  fi

  echo -e "${cond}\t${sorted_bed}\t${cov}" >> "${OUT_TSV}"
done

echo ""
echo "Wrote condition table: ${OUT_TSV}"
wc -l "${OUT_TSV}"

# Count data rows (excluding header)
N=$(( $(wc -l < "${OUT_TSV}") - 1 ))
echo "Conditions found: N = ${N}"

if (( N == 0 )); then
  echo "ERROR: No valid conditions written to TSV. Check BED and covariate paths." >&2
  exit 1
fi

### 4 — Submit array job
echo ""
echo "=== Step 3: Submitting sbatch array (1-${N}) ==="
sbatch --array=1-"${N}" "${SBATCH_SCRIPT}"
