#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
OUT="$PROJECT/06_feasibility/stage3"
LOG="$PROJECT/logs/stage3_output_validation.log"

/usr/bin/time -v -o "$OUT/stage3_output_validation.time.txt" \
  Rscript "$PROJECT/04_scripts/stage3_validate_outputs.R" >"$LOG" 2>&1
Rscript "$PROJECT/04_scripts/stage3_compile_runtime_memory.R" >>"$LOG" 2>&1
gzip -t "$PROJECT/02_processed_data/malignant_epithelial/10X_counts/matrix.mtx.gz"
gzip -t "$PROJECT/02_processed_data/malignant_epithelial/target_counts/matrix.mtx.gz"
echo "STAGE3_OUTPUT_VALIDATION_OK"
