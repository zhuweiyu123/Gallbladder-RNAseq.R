#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

root <- "/home/zhuweiyu/codex-r/results/merged_cnv_binary_annotation"
rds_path <- file.path(root, "merged_cnv_binary_annotated.rds")
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

obj <- readRDS(rds_path)
meta <- as.data.table(obj[[]], keep.rownames = "barcode")
malignant_statuses <- c("malignant_high_confidence", "malignant_likely")
normal_statuses <- c("cnv_normal_like", "uncertain", "unresolved_low_quality")
is_epi <- as.character(meta$celltype) == "Epithelial"
is_malignant <- meta$epithelial_binary_label %in% "Malignant epithelial cells"
is_normal <- meta$epithelial_binary_label %in% "Normal epithelial cells"

required_files <- c(
  "tables/cell_level_annotation_with_cnv_binary.csv.gz",
  "tables/epithelial_binary_counts_by_patient_status.tsv",
  "tables/final_celltype_counts_by_patient.tsv",
  "tables/annotation_transfer_summary.tsv",
  "differential_expression/markers_malignant_vs_normal_epithelial_celllevel_wilcox.csv",
  "differential_expression/markers_stk31_high_vs_low_malignant_epithelial.csv",
  "differential_expression/markers_stk31_high_malignant_epithelial_vs_refined_nk.csv",
  "differential_expression/markers_stk31_high_malignant_epithelial_vs_all_other.csv",
  "differential_expression/pseudobulk_edgeR_paired_malignant_vs_normal.csv",
  "differential_expression/pseudobulk_sample_metadata.tsv",
  "figures/UMAP_final_celltype_cnv_binary.png",
  "figures/UMAP_epithelial_binary_by_patient.png",
  "figures/UMAP_STK31_high_low_malignant_epithelial.png",
  "figures/Volcano_pseudobulk_malignant_vs_normal_epithelial.png",
  "README.md", "sessionInfo.txt"
)

checks <- data.table(
  check = c(
    "rds_cells", "barcode_unique", "umap_present", "old_annotation_matched",
    "cnv_matched", "epithelial_total", "malignant_total", "normal_total",
    "malignant_status_mapping_exact", "normal_status_mapping_exact",
    "missing_cnv_epithelial_normal", "missing_cnv_source_explicit",
    "non_epithelial_not_binary_epithelial", "no_unsplit_epithelial_final_label",
    "old_epithelial_conflicts_resolved", "stk31_high_malignant",
    "stk31_low_malignant", "all_required_files_present_nonempty"
  ),
  observed = as.character(c(
    ncol(obj), !anyDuplicated(meta$barcode), "umap" %in% names(obj@reductions),
    sum(meta$previous_annotation_matched), sum(meta$cnv_metadata_matched),
    sum(is_epi), sum(is_malignant, na.rm = TRUE), sum(is_normal, na.rm = TRUE),
    all(meta$cnv_consensus_status[is_malignant] %in% malignant_statuses),
    all(meta$cnv_consensus_status[is_normal & !is.na(meta$cnv_consensus_status)] %in% normal_statuses),
    sum(is_epi & is.na(meta$cnv_consensus_status) & is_normal, na.rm = TRUE),
    all(meta$analysis_celltype_cnv_binary_source[
      is_epi & is.na(meta$cnv_consensus_status)
    ] == "new_epithelial_without_cnv_assigned_normal_by_user_rule"),
    all(!meta$analysis_celltype_cnv_binary[!is_epi] %in%
          c("Malignant epithelial cells", "Normal epithelial cells")),
    !any(meta$analysis_celltype_cnv_binary == "Epithelial"),
    sum(meta$analysis_celltype_cnv_binary_source ==
          "new_annotation_resolved_old_epithelial_conflict"),
    sum(meta$stk31_malignant_group == "STK31-high malignant epithelial", na.rm = TRUE),
    sum(meta$stk31_malignant_group == "STK31-low malignant epithelial", na.rm = TRUE),
    all(file.exists(file.path(root, required_files)) & file.info(file.path(root, required_files))$size > 0)
  )),
  expected = as.character(c(
    27108, TRUE, TRUE, 25821, 4304, 4763, 1047, 3716,
    TRUE, TRUE, 459, TRUE, TRUE, TRUE, 89, 117, 930, TRUE
  ))
)
checks[, pass := observed == expected]
fwrite(checks, file.path(root, "tables", "independent_readback_validation.tsv"), sep = "\t")
print(checks)
stop_if_not(all(checks$pass), "Independent read-back validation failed")

pb_meta <- fread(file.path(root, "differential_expression", "pseudobulk_sample_metadata.tsv"))
stop_if_not(
  identical(sort(unique(as.character(pb_meta$patient_id))), c("P1", "P2", "P4")),
  "Unexpected pseudobulk paired patients"
)
stop_if_not(nrow(pb_meta) == 6L, "Unexpected pseudobulk sample count")

message("INDEPENDENT_READBACK_PASS")
