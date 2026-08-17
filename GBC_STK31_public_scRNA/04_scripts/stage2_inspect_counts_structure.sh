#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
COUNTS="${PROJECT}/02_processed_data/All_single_cells/counts/10X_counts"
OUT="${PROJECT}/06_feasibility/counts_structure.txt"

{
  printf 'inspection_scope: headers, feature names, barcodes, dimensions; matrix body not loaded\n'
  printf 'matrix_format: gzip-compressed Matrix Market coordinate sparse matrix\n'
  printf 'matrix_header:\n'
  gzip -cd "${COUNTS}/matrix.mtx.gz" | head -5 || true
  printf 'feature_examples:\n'
  gzip -cd "${COUNTS}/features.tsv.gz" | head -5 || true
  printf 'feature_rows: '
  gzip -cd "${COUNTS}/features.tsv.gz" | wc -l
  printf 'barcode_examples:\n'
  gzip -cd "${COUNTS}/barcodes.tsv.gz" | head -5 || true
  printf 'barcode_rows: '
  gzip -cd "${COUNTS}/barcodes.tsv.gz" | wc -l
  printf 'file_details:\n'
  find "${COUNTS}" -maxdepth 1 -type f -printf '%f\t%s bytes\n' | sort
  printf 'file_types:\n'
  file "${COUNTS}"/*
} > "${OUT}"

cat "${OUT}"
