#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
input <- file.path(project, "02_processed_data", "All_single_cells", "sample_info.xlsx")
out <- file.path(project, "06_feasibility", "sample_info_audit.txt")

x <- as.data.table(read_excel(input, skip = 1, .name_repair = "unique"))
x <- x[!is.na(`Sample ID`) & trimws(as.character(`Sample ID`)) != ""]

lines <- c(
  paste0("rows: ", nrow(x)),
  paste0("columns: ", ncol(x)),
  paste0("column_names: ", paste(names(x), collapse = " | "))
)

for (nm in names(x)) {
  vals <- unique(as.character(x[[nm]]))
  vals <- vals[!is.na(vals)]
  lines <- c(
    lines,
    paste0("field.", nm, ".n_unique_nonmissing: ", length(vals)),
    paste0("field.", nm, ".values: ", paste(vals, collapse = " | "))
  )
}

writeLines(lines, out, useBytes = TRUE)
cat(paste(lines, collapse = "\n"), "\n")
