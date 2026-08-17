#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
SCRIPT_DIR="${PROJECT}/04_scripts"
MAL_DIR="${PROJECT}/02_processed_data/malignant_epithelial"
STAGE3_OUT="${PROJECT}/06_feasibility/stage3"
LOG_DIR="${PROJECT}/logs"
SOURCE="${SCRIPT_DIR}/stage3_stream_extract_malignant.cpp"
BINARY="${PROJECT}/00_tools/bin/stage3_stream_extract_malignant"
INPUT="${PROJECT}/02_processed_data/All_single_cells/counts/10X_counts/matrix.mtx.gz"

mkdir -p "${MAL_DIR}/10X_counts" "${MAL_DIR}/target_counts" "${STAGE3_OUT}" "${LOG_DIR}" "${PROJECT}/00_tools/bin"

g++ -O3 -std=c++17 -Wall -Wextra -Wpedantic "${SOURCE}" -lz -o "${BINARY}"

/usr/bin/time -v -o "${STAGE3_OUT}/stage3_prepare_manifest.time.txt" \
  Rscript "${SCRIPT_DIR}/stage3_prepare_manifest.R" \
  > "${LOG_DIR}/stage3_prepare_manifest.log" 2>&1

/usr/bin/time -v -o "${STAGE3_OUT}/stage3_stream_extraction.time.txt" \
  "${BINARY}" \
    "${INPUT}" \
    "${MAL_DIR}/full_to_malignant_col.int32.bin" \
    "${MAL_DIR}/full_to_target_row.int32.bin" \
    32137 1117245 96942 9 \
    "${MAL_DIR}/10X_counts/matrix.mtx.gz" \
    "${MAL_DIR}/target_counts/matrix.mtx.gz" \
    "${MAL_DIR}/malignant_cell_library_sizes.tsv" \
    "${STAGE3_OUT}/stage3_stream_extraction_audit.txt" \
    "${MAL_DIR}/stage3_extract_tmp" \
  > "${LOG_DIR}/stage3_stream_extraction.log" 2>&1

printf 'STAGE3_EXTRACTION_PIPELINE_OK\n'
