#!/bin/bash
# Convert TensorQTL parquet files to tab-delimited text files
# =============================================================================
# USER CONFIGURATION
# =============================================================================

### Base tensor directory
BASE_DIR="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"

### Directory containing .parquet files to convert
INPUT_DIR="${BASE_DIR}/tensorqtl_output_cis-nominal"

### Directory to write .txt output files
OUTPUT_DIR="${INPUT_DIR}/txt_files"

### Output file delimiter (use $'\t' for tab, ',' for CSV, etc.)
DELIMITER=$'\t'

### Output file extension (e.g. txt, tsv, csv)
OUT_EXT="txt"

### Number of header rows to preview per file in the summary (0 to disable)
PREVIEW_ROWS=5

### Max files to list in the summary (0 to disable)
SUMMARY_LIST_MAX=10

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================

set -euo pipefail

mkdir -p "${OUTPUT_DIR}"

echo "[$(date)] Starting parquet to ${OUT_EXT} conversion..."
echo "Input directory  : ${INPUT_DIR}"
echo "Output directory : ${OUTPUT_DIR}"
echo "Delimiter        : $(cat -A <<< "${DELIMITER}")"

### Collect parquet files
shopt -s nullglob
parquet_files=( "${INPUT_DIR}"/*.parquet )

if (( ${#parquet_files[@]} == 0 )); then
  echo "ERROR: No .parquet files found in ${INPUT_DIR}" >&2
  exit 1
fi

echo "Found ${#parquet_files[@]} parquet file(s) to convert"

# =============================================================================
# CONVERT EACH FILE
# =============================================================================

converted=0
failed=0

for parquet_file in "${parquet_files[@]}"; do
  base_name=$(basename "${parquet_file}" .parquet)
  out_file="${OUTPUT_DIR}/${base_name}.${OUT_EXT}"

  echo ""
  echo "Converting: $(basename "${parquet_file}")"

  python3 - "${parquet_file}" "${out_file}" "${DELIMITER}" <<'PYEOF'
import sys
import pandas as pd

parquet_path, out_path, delimiter = sys.argv[1], sys.argv[2], sys.argv[3]

try:
    df = pd.read_parquet(parquet_path)
    print(f"  Shape   : {df.shape[0]:,} rows x {df.shape[1]} columns")
    print(f"  Columns : {', '.join(df.columns)}")
    df.to_csv(out_path, sep=delimiter, index=False)
    print(f"  Saved   : {out_path}")
except Exception as e:
    print(f"  ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF

  if [[ $? -eq 0 ]]; then
    (( converted++ ))
    echo "  ✓ SUCCESS  ($(du -h "${out_file}" | cut -f1))"
  else
    (( failed++ ))
    echo "  ✗ FAILED"
  fi

done

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=============================================================="
echo "[$(date)] Conversion complete!"
echo "=============================================================="
echo "Successfully converted : ${converted} file(s)"
echo "Failed                 : ${failed} file(s)"
echo "Output directory       : ${OUTPUT_DIR}"

if (( converted > 0 )); then

  if (( SUMMARY_LIST_MAX > 0 )); then
    echo ""
    echo "Sample of converted files:"
    ls -lh "${OUTPUT_DIR}/"*".${OUT_EXT}" | head -"${SUMMARY_LIST_MAX}"
  fi

  if (( PREVIEW_ROWS > 0 )); then
    first_out=$(ls "${OUTPUT_DIR}/"*".${OUT_EXT}" | head -1)
    echo ""
    echo "Preview of: $(basename "${first_out}") (first ${PREVIEW_ROWS} rows)"
    head -"$(( PREVIEW_ROWS + 1 ))" "${first_out}"
  fi

fi

if (( failed > 0 )); then
  echo ""
  echo "WARNING: ${failed} file(s) failed to convert. Check output above for details." >&2
  exit 1
fi

echo ""
echo "All done!"
