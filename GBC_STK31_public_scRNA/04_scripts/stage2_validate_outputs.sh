#!/usr/bin/env bash
set -euo pipefail

PROJECT="/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
cd "${PROJECT}"

required=(
  metadata_summary.csv
  celltype_counts.csv
  sample_celltype_counts.csv
  patient_sample_mapping.csv
  sample_site_mapping.csv
  target_gene_presence.csv
  nk_rds_structure.txt
  counts_structure.txt
  stage2_summary.md
)

for name in "${required[@]}"; do
  test -s "06_feasibility/${name}"
  printf 'OK\t06_feasibility/%s\n' "${name}"
done

for script in 04_scripts/*.sh; do
  bash -n "${script}"
done
printf 'BASH_PARSE_OK\n'

Rscript -e 'fs <- list.files("04_scripts", pattern="[.]R$", full.names=TRUE); for (f in fs) parse(file=f); cat("R_PARSE_OK", length(fs), "\n")'
Rscript -e 'library(data.table); fs <- list.files("06_feasibility", pattern="[.]csv$", full.names=TRUE); for (f in fs) fread(f, nrows=2); cat("CSV_READ_OK", length(fs), "\n")'

md5sum 01_raw_data/scRNA_processed_data_all.zip 01_raw_data/NK.RDS

if pgrep -af 'Rscript .*GBC_STK31|stage2_inspect|stage2_audit' | grep -v stage2_validate_outputs; then
  printf 'WARNING_ACTIVE_STAGE2_PROCESS\n' >&2
  exit 4
fi
printf 'NO_ACTIVE_STAGE2_ANALYSIS_PROCESS\n'
printf 'STAGE2_VALIDATION_OK\n'
