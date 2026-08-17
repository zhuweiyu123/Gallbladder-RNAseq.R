#!/usr/bin/env Rscript

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
scripts <- file.path(project, "04_scripts", c(
  "stage4A_current_paper_minimal.R"
))
for (script in scripts) parse(file = script)
cat("STAGE4A_R_PARSE_OK", length(scripts), "\n")
