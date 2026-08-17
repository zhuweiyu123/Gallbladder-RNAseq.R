#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
OUT="${PROJECT}/06_feasibility/stage2_environment.txt"

{
  date --iso-8601=seconds
  uname -a
  printf '\nDISK\n'
  df -h "${PROJECT}"
  df -B1 "${PROJECT}"
  printf '\nMEMORY\n'
  free -h
  free -b
  printf '\nR_PACKAGES\n'
  Rscript -e 'pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table", "readxl"); cat("R=", R.version.string, "\n", sep=""); for (p in pkgs) cat(p, "=", if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT_INSTALLED", "\n", sep="")'
  printf '\nDOWNLOAD_MD5\n'
  md5sum "${PROJECT}/01_raw_data/scRNA_processed_data_all.zip" "${PROJECT}/01_raw_data/NK.RDS"
  printf '\nPROJECT_FILE_SIZES\n'
  stat --printf='%n\t%s bytes\n' \
    "${PROJECT}/01_raw_data/scRNA_processed_data_all.zip" \
    "${PROJECT}/01_raw_data/NK.RDS" \
    "${PROJECT}/02_processed_data/All_single_cells/GBC_counts.rar" \
    "${PROJECT}/02_processed_data/All_single_cells/GBC_Metadata.txt" \
    "${PROJECT}/02_processed_data/All_single_cells/sample_info.xlsx"
} > "${OUT}"

cat "${OUT}"
