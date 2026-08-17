#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
out <- file.path(project, "06_feasibility", "stage3")

components <- data.table(
  component = c("manifest_and_index", "stream_extract_malignant", "malignant_detection", "NK_detection", "primary_matching", "output_validation"),
  file = file.path(out, c(
    "stage3_prepare_manifest.time.txt",
    "stage3_stream_extraction.time.txt",
    "stage3_malignant_detection.time.txt",
    "stage3_nk_detection.time.txt",
    "stage3_primary_matching.time.txt",
    "stage3_output_validation.time.txt"
  ))
)

extract_field <- function(lines, field) {
  prefix <- paste0(trimws(field), ":")
  normalized <- trimws(lines)
  hit <- normalized[startsWith(normalized, prefix)]
  if (!length(hit)) return(NA_character_)
  trimws(substring(hit[1], nchar(prefix) + 1L))
}

runtime <- components[, {
  if (!file.exists(file)) {
    .(command = NA_character_, elapsed = NA_character_, max_rss_kb = NA_real_, exit_status = NA_integer_)
  } else {
    lines <- readLines(file, warn = FALSE)
    .(
      command = extract_field(lines, "Command being timed"),
      elapsed = extract_field(lines, "Elapsed (wall clock) time (h:mm:ss or m:ss)"),
      max_rss_kb = as.numeric(extract_field(lines, "Maximum resident set size (kbytes)")),
      exit_status = as.integer(extract_field(lines, "Exit status"))
    )
  }
}, by = .(component, file)]
runtime[, `:=`(
  max_rss_mib = max_rss_kb / 1024,
  max_rss_gib = max_rss_kb / 1024^2
)]
fwrite(runtime, file.path(out, "stage3_runtime_memory.csv"))

peak_index <- which.max(runtime$max_rss_kb)
lines <- c(
  "Stage 3 runtime and memory audit",
  paste0("generated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R_version=", R.version.string),
  paste0("platform=", R.version$platform),
  paste0("peak_component=", runtime$component[peak_index]),
  paste0("peak_max_rss_kb=", format(runtime$max_rss_kb[peak_index], scientific = FALSE)),
  paste0("peak_max_rss_gib=", sprintf("%.3f", runtime$max_rss_gib[peak_index])),
  "",
  apply(runtime, 1L, function(row) paste(row, collapse = "\t")),
  "",
  "Notes:",
  "- The full 1,117,245-cell matrix was never materialized in R.",
  "- The full malignant matrix remains sparse Matrix Market; no dense conversion was performed.",
  "- The first NK attempt failed before output generation due to a script-length check and is preserved separately as *.failed_attempt.* in logs/time files."
)
writeLines(lines, file.path(out, "stage3_runtime_memory.txt"), useBytes = TRUE)
cat(paste(lines[1:7], collapse = "\n"), "\n")
