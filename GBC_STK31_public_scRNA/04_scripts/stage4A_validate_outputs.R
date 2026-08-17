#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
fig_dir <- file.path(project, "11_figures", "current_paper_minimal")
table_dir <- file.path(project, "12_tables", "current_paper_minimal")
technical_dir <- file.path(table_dir, "technical")
report <- file.path(project, "06_feasibility", "stage4A_current_paper_minimal_summary.md")

checks <- list()
add_check <- function(check, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = check, passed = isTRUE(passed), detail = detail)
  if (!isTRUE(passed)) stop(sprintf("Stage4A validation failed: %s (%s)", check, detail), call. = FALSE)
}

figure_stems <- c(
  "Figure1_STK31_celltype_localization",
  "Figure2_malignant_epithelial_HLAI",
  "Figure3_NK_inhibitory_receptor_subtypes",
  "Figure4_primary_STK31_HLAI_pseudobulk_trends",
  "Supplementary_FigureS1_all_malignant_STK31_HLAI_pseudobulk_trends",
  "Supplementary_FigureS2_KIR3DL1_subtype_detection"
)
for (stem in figure_stems) {
  for (extension in c("png", "pdf")) {
    path <- file.path(fig_dir, paste0(stem, ".", extension))
    add_check(paste0("figure_", basename(path)), file.exists(path) && file.info(path)$size > 1000, path)
  }
}

major <- fread(file.path(table_dir, "Figure1_STK31_major_celltype_stats.csv"))
hla <- fread(file.path(table_dir, "Figure2_malignant_HLAI_expression_stats.csv"))
nk <- fread(file.path(table_dir, "Figure3_NK_inhibitory_receptor_subtype_stats.csv"))
pb <- fread(file.path(table_dir, "Figure4_malignant_TMM_pseudobulk_target_genes.csv"))
correlations <- fread(file.path(table_dir, "Figure4_STK31_HLAI_Spearman_correlations.csv"))
audit <- fread(file.path(technical_dir, "stage4A_malignant_pseudobulk_stream_audit.tsv"))

add_check("figure1_major_celltype_rows", nrow(major) == 7L, nrow(major))
add_check("figure1_detection_bounds", all(major$STK31_detected_cell_n <= major$cell_n) && all(major$STK31_RNA_detected_fraction >= 0 & major$STK31_RNA_detected_fraction <= 1))
add_check("figure2_HLA_main_rows", nrow(hla) == 4L && setequal(as.character(hla$gene), c("HLA-A", "HLA-B", "HLA-C", "B2M")), nrow(hla))
add_check("figure2_malignant_denominator", all(hla$malignant_cell_n == 96942L) && all(hla$RNA_detected_n <= hla$malignant_cell_n))
add_check("figure3_NK_subtype_receptor_grid", nrow(nk) == 70L && uniqueN(nk$nk_subtype) == 14L && uniqueN(nk$gene) == 5L)
add_check("figure3_KIR2DL2_absent", !"KIR2DL2" %in% nk$gene)
add_check("figure3_detection_bounds", all(nk$RNA_detected_n <= nk$NK_cell_n) && all(nk$RNA_detected_fraction >= 0 & nk$RNA_detected_fraction <= 1))
add_check("pseudobulk_sample_n", nrow(pb) == 58L && uniqueN(pb$sample_name) == 58L)
add_check("pseudobulk_primary_n", pb[standardized_site == "Primary tumor", .N] == 11L)
add_check("pseudobulk_required_targets", all(paste0(c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M"), "_log2CPM") %in% names(pb)))
add_check("pseudobulk_stream_integrity", audit$declared_nnz == audit$actual_nnz && audit$output_rows == 32137L && audit$output_cols == 58L)
add_check("correlation_table_shape", nrow(correlations) == 8L && setequal(correlations$scope, c("Primary tumor", "All malignant-containing samples")))
add_check("primary_correlations_n", all(correlations[scope == "Primary tumor", n] == 11L))
add_check("summary_exists", file.exists(report) && file.info(report)$size > 1000, report)

result <- rbindlist(checks)
fwrite(result, file.path(table_dir, "stage4A_validation_checks.csv"))
writeLines(c(
  "STAGE4A_VALIDATION_OK",
  paste0("checks_total=", nrow(result)),
  paste0("checks_passed=", sum(result$passed)),
  paste0("checks_failed=", sum(!result$passed))
), file.path(table_dir, "stage4A_validation_summary.txt"), useBytes = TRUE)
cat("STAGE4A_VALIDATION_OK", nrow(result), "\n")
