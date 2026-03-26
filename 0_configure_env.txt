#!/usr/bin/env bash
# =============================================================================
# 00_setup_tensorQTL.sh
# Environment setup for tensorQTL eQTL mapping
#
# Run once to create the conda environment before launching any analysis.
# Adapts an existing environment and installs the additional packages
# needed for q-value computation via R/rpy2.
# =============================================================================
 
set -euo pipefail
 
ENV_NAME="tensorqtl_p3.11_env"
BASE_YAML="/wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml"
 
# ---------------------------------------------------------------------------
# 1. Create environment from collaborator YAML
# ---------------------------------------------------------------------------
echo "Creating conda environment '${ENV_NAME}' from ${BASE_YAML}..."
conda env create --name "${ENV_NAME}" --file="${BASE_YAML}"
 
conda activate "${ENV_NAME}"
 
# Point to the correct R shared library (using R 4.3.2)
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
 
# Verify tensorQTL loads
echo "Verifying tensorQTL installation..."
python3 -m tensorqtl --help | head
# Expected warning at this stage:
#   Warning: 'rfunc' cannot be imported.
#   R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values.
 
# ---------------------------------------------------------------------------
# 2. Install rpy2 and R dependencies for q-value computation
# ---------------------------------------------------------------------------
echo "Installing rpy2..."
pip3 install 'rpy2<=3.5.12'
 
echo "Configuring conda channels for R packages..."
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y r-essentials
 
# ---------------------------------------------------------------------------
# 3. Install the Bioconductor qvalue package inside R
# ---------------------------------------------------------------------------
echo "Installing Bioconductor qvalue package in R..."
Rscript - <<'EOF'
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)
cat("qvalue version:", as.character(packageVersion("qvalue")), "\n")
EOF
 
# ---------------------------------------------------------------------------
# 4. Final verification
# ---------------------------------------------------------------------------
echo "Final verification..."
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
 
python3 -c "
from rpy2.robjects.packages import importr
qvalue = importr('qvalue')
print('rpy2 <-> R bridge OK')
print('qvalue package accessible from Python')
"
 
python3 -m tensorqtl --help
echo ""
echo "Setup complete. Activate the environment with:"
echo "  conda activate ${ENV_NAME}"
echo "  export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:\$LD_LIBRARY_PATH"
 
