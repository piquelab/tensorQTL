#!/usr/bin/env bash
# =============================================================================
# tensorQTL cis-nominal SBATCH Script
# Configure the variables in the USER CONFIGURATION section below.
# =============================================================================
#SBATCH -J tQTL
#SBATCH -q primary
#SBATCH -n 4
#SBATCH --mem=24G
#SBATCH -t 08:00:00
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err

# =============================================================================
# USER CONFIGURATION
# =============================================================================

### Base tensor directory (where PLINK files and condition TSV live)
BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"

### PLINK2 prefix (no extension — expects .pgen/.pvar/.psam)
PLINK_PREFIX="${BASE}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"

### Condition table produced by the pre-processing script
COND_TSV="${BASE}/tensorqtl_conditions.tsv"

### Output directory for tensorQTL results
OUT_DIR="${BASE}/tensorqtl_output"

### tensorQTL run parameters
MODE="cis_nominal"        # tensorQTL mode (e.g. cis, cis_nominal, cis_independent, interaction, etc.)
FDR=0.1                   # FDR threshold
WINDOW=100000             # cis-window in bp

OUT_DIR="${OUT_DIR}/${MODE}"

### Conda environment name
CONDA_ENV="tensorqtl_p3.11"

### Extra LD_LIBRARY_PATH entries needed (space-separated, or leave empty)
LD_PATHS="/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib"

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================

set -euo pipefail

### Output directory
mkdir -p "${OUT_DIR}"

### Environment
module purge
module load plink/2.0
source activate "${CONDA_ENV}"

if [[ -n "${EXTRA_LD_PATHS}" ]]; then
  export LD_LIBRARY_PATH="${EXTRA_LD_PATHS}:${LD_LIBRARY_PATH:-}"
fi

python3 -m tensorqtl --help >/dev/null 2>&1 \
  || { echo "ERROR: tensorQTL not found in environment '${CONDA_ENV}'"; exit 1; }

### Array task — skip header; get the SLURM_ARRAY_TASK_ID-th data line
line=$(awk -v n="${SLURM_ARRAY_TASK_ID}" \
  'NR==1{next} NR==n+1{print; exit}' "${COND_TSV}") || true

if [[ -z "${line:-}" ]]; then
  echo "No line for task ${SLURM_ARRAY_TASK_ID} in ${COND_TSV}" >&2
  exit 0
fi

IFS=$'\t' read -r CONDITION PHENO_BED COV_FILE <<< "${line}"

### Output prefix for this condition
OUT_PREFIX="${CONDITION}_${MODE}"

### Log run details
echo "[$(date)] Starting tensorQTL ${MODE}"
echo "Condition   : ${CONDITION}"
echo "Genotypes   : ${PLINK_PREFIX}.[pgen|pvar|psam]"
echo "Phenotypes  : ${PHENO_BED}"
echo "Covariates  : ${COV_FILE}"
echo "Mode        : ${MODE}"
echo "FDR         : ${FDR}"
echo "Window (bp) : ${WINDOW}"
echo "Out dir     : ${OUT_DIR}"
echo "Out prefix  : ${OUT_PREFIX}"

### Run tensorQTL
python3 -m tensorqtl \
  "${PLINK_PREFIX}" \
  "${PHENO_BED}" \
  "${OUT_PREFIX}" \
  --covariates "${COV_FILE}" \
  --mode "${MODE}" \
  --fdr "${FDR}" \
  --window "${WINDOW}" \
  -o "${OUT_DIR}"

echo "[$(date)] Done."
