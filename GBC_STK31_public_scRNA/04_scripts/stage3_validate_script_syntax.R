#!/usr/bin/env Rscript

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
scripts <- file.path(project, "04_scripts", c(
  "stage3_prepare_manifest.R",
  "stage3_malignant_detection.R",
  "stage3_nk_detection.R",
  "stage3_primary_matching.R"
))
for (script in scripts) parse(file = script)
cat("STAGE3_R_PARSE_OK", length(scripts), "\n")
