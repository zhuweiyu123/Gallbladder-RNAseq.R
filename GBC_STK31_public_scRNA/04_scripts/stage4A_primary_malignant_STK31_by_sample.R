#!/usr/bin/env Rscript

# Per-sample descriptive STK31 summary for author-annotated malignant epithelial
# cells from the 11 Primary tumor samples. Reuses Stage 3 and Stage 4A tables;
# no expression matrix is read and no statistical test is performed.

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
stage3_path <- file.path(project, "06_feasibility", "stage3", "malignant_sample_gene_detection.csv")
tmm_path <- file.path(project, "12_tables", "current_paper_minimal", "Figure4_malignant_TMM_pseudobulk_target_genes.csv")
output_path <- file.path(project, "12_tables", "current_paper_minimal", "Primary_tumor_malignant_STK31_by_sample.csv")

stage3 <- fread(stage3_path, showProgress = FALSE)
sample_stk <- stage3[
  standardized_site == "Primary tumor" & gene == "STK31",
  .(sample_name, patient_id, patient_id_source, site, standardized_site,
    TNM_stage, histological_type, metastatic_group,
    malignant_epithelial_cell_n = malignant_cell_n,
    STK31_RNA_detected_cell_n = RNA_detected_n,
    STK31_RNA_detected_fraction = RNA_detected_fraction,
    STK31_raw_UMI = pseudobulk_raw_count,
    STK31_raw_library_CPM = pseudobulk_CPM)
]
stopifnot(nrow(sample_stk) == 11L, !anyDuplicated(sample_stk$sample_name))

tmm <- fread(tmm_path, showProgress = FALSE)[
  standardized_site == "Primary tumor",
  .(sample_name, STK31_TMM_CPM, STK31_log2CPM)
]
stopifnot(nrow(tmm) == 11L, !anyDuplicated(tmm$sample_name))

result <- merge(sample_stk, tmm, by = "sample_name", all.x = TRUE, sort = FALSE)
stopifnot(nrow(result) == 11L, !anyNA(result$STK31_TMM_CPM))
setorder(result, sample_name)
fwrite(result, output_path)
message("STAGE4A_PRIMARY_MALIGNANT_STK31_BY_SAMPLE_OK n=", nrow(result), " total_malignant_cells=", sum(result$malignant_epithelial_cell_n))
