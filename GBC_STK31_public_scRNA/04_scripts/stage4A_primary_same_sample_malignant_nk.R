#!/usr/bin/env Rscript

# Same-primary-sample malignant epithelial–NK pairing.
# Descriptive only; no cross-site matching and no expression-based filtering.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
stage3_dir <- file.path(project, "06_feasibility", "stage3")
table_dir <- file.path(project, "12_tables", "current_paper_minimal")
fig_dir <- file.path(project, "11_figures", "current_paper_minimal")
report_dir <- file.path(project, "06_feasibility")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
}

expected_exact <- c(
  "GBC_005_P", "GBC_023_P", "GBC_033_P", "GBC_034_P", "GBC_038_P",
  "GBC_042_P", "GBC_047_P", "GBC_056_P", "GBC_073_P", "GBC_075_P", "GBC_086_P"
)
expected_qc <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")

pairs <- fread(file.path(stage3_dir, "primary_matched_stk31_hla_nk.csv"), showProgress = FALSE)
stop_if_not(nrow(pairs) == 11L && !anyDuplicated(pairs$sample_name), "Expected 11 unique exact primary pairs")
setorder(pairs, sample_name)
stop_if_not(identical(pairs$sample_name, expected_exact), "Exact-pair sample list differs from the locked primary list")
stop_if_not(all(pairs$standardized_site == "Primary tumor"), "A non-primary sample entered the pair table")
stop_if_not(all(pairs$match_scope == "exact_same_primary_sample"), "Pairing is not exact same-sample primary matching")
stop_if_not(all(pairs$has_malignant_compartment & pairs$has_NK_compartment), "A paired sample lacks a required compartment")

pairs[, `:=`(
  pass_malignant_100 = malignant_cell_n >= 100L,
  pass_NK_100 = NK_cell_n >= 100L
)]
pairs[, pass_paired_QC := pass_malignant_100 & pass_NK_100]
stop_if_not(identical(pairs[pass_paired_QC == TRUE, sample_name], expected_qc), "Paired QC set differs from the locked four-sample set")

# Use the already completed Primary-8-only TMM normalization for malignant genes.
mal_tmm <- fread(file.path(table_dir, "Primary8_STK31_HLAI_TMM_CPM.csv"), showProgress = FALSE)
mal_keep <- mal_tmm[, .(
  sample_name,
  malignant_STK31_TMM_CPM = STK31_TMM_CPM,
  malignant_STK31_log2CPM = STK31_log2CPM,
  malignant_HLA.A_TMM_CPM = `HLA-A_TMM_CPM`,
  malignant_HLA.B_TMM_CPM = `HLA-B_TMM_CPM`,
  malignant_HLA.C_TMM_CPM = `HLA-C_TMM_CPM`,
  malignant_B2M_TMM_CPM = B2M_TMM_CPM
)]
pairs <- merge(pairs, mal_keep, by = "sample_name", all.x = TRUE, sort = FALSE)
setorder(pairs, sample_name)
stop_if_not(!anyNA(pairs[pass_paired_QC == TRUE, malignant_STK31_TMM_CPM]), "A paired-QC sample lacks Primary-8 malignant TMM values")

pair_output_columns <- c(
  "sample_name", "patient_id", "standardized_site", "TNM_stage", "histological_type",
  "malignant_cell_n", "NK_cell_n", "pass_malignant_100", "pass_NK_100", "pass_paired_QC",
  "malignant_STK31_TMM_CPM", "malignant_HLA.A_TMM_CPM", "malignant_HLA.B_TMM_CPM",
  "malignant_HLA.C_TMM_CPM", "malignant_B2M_TMM_CPM",
  "STK31_RNA_detected_n", "STK31_RNA_detected_fraction",
  "KIR3DL1_RNA_detected_n", "KIR3DL1_RNA_detected_fraction", "KIR3DL1_pseudobulk_CPM",
  "KIR2DL1_RNA_detected_fraction", "KIR2DL3_RNA_detected_fraction",
  "KLRC1_RNA_detected_fraction", "KLRD1_RNA_detected_fraction",
  "NKG7_RNA_detected_fraction", "GNLY_RNA_detected_fraction", "PRF1_RNA_detected_fraction",
  "GZMB_RNA_detected_fraction", "FGFBP2_RNA_detected_fraction", "FCGR3A_RNA_detected_fraction",
  "XCL1_RNA_detected_fraction", "XCL2_RNA_detected_fraction", "match_scope"
)
stop_if_not(all(pair_output_columns %in% names(pairs)), "A required pair-table column is absent")
fwrite(pairs[, ..pair_output_columns], file.path(table_dir, "Primary_same_sample_malignant_NK_pairs.csv"))

qc <- pairs[pass_paired_QC == TRUE]
stop_if_not(nrow(qc) == 4L, "Paired QC analysis must contain four samples")
cor_specs <- data.table(
  malignant_variable = c(
    "malignant_STK31_TMM_CPM", "malignant_HLA.B_TMM_CPM",
    "malignant_STK31_TMM_CPM", "malignant_HLA.B_TMM_CPM",
    "malignant_STK31_TMM_CPM"
  ),
  NK_variable = c(
    "KIR3DL1_RNA_detected_fraction", "KIR3DL1_RNA_detected_fraction",
    "KLRC1_RNA_detected_fraction", "KLRC1_RNA_detected_fraction",
    "NKG7_RNA_detected_fraction"
  )
)
correlations <- cor_specs[, {
  x <- qc[[malignant_variable]]
  y <- qc[[NK_variable]]
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  .(
    n = length(x),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    interpretation = "exploratory_n4_not_for_inference"
  )
}, by = .(malignant_variable, NK_variable)]
fwrite(correlations, file.path(table_dir, "Primary_same_sample_malignant_NK_correlations.csv"))

plot_pairs <- copy(pairs)
plot_pairs[, sample_name := factor(sample_name, levels = expected_exact)]
plot_pairs[, QC_label := fifelse(pass_paired_QC, "Both compartments >=100 cells", "Low-count descriptive only")]
pair_colors <- c("Both compartments >=100 cells" = "#1f78b4", "Low-count descriptive only" = "#9e9e9e")

count_long <- melt(
  plot_pairs[, .(sample_name, QC_label, malignant_cell_n, NK_cell_n)],
  id.vars = c("sample_name", "QC_label"),
  variable.name = "compartment", value.name = "cell_n"
)
count_long[, compartment := factor(compartment, levels = c("malignant_cell_n", "NK_cell_n"), labels = c("Malignant epithelial", "NK"))]
panel_a <- ggplot(count_long, aes(sample_name, cell_n, color = QC_label, shape = compartment, group = compartment)) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_color_manual(values = pair_colors, name = NULL) +
  scale_shape_manual(values = c("Malignant epithelial" = 16, "NK" = 17), name = "Compartment") +
  labs(x = NULL, y = "Cell number (log10 scale)", title = "A. Same-primary-sample compartment sizes") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

panel_b <- ggplot(qc, aes(malignant_HLA.B_TMM_CPM, KIR3DL1_RNA_detected_fraction, label = sample_name)) +
  geom_point(size = 3.2, color = "#1f78b4") +
  geom_text(vjust = -0.8, size = 3.1, check_overlap = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Malignant HLA-B TMM-CPM (Primary-8 normalization)", y = "NK KIR3DL1 RNA detection fraction",
       title = "B. HLA-B vs KIR3DL1 (paired QC, n=4)") +
  theme_classic(base_size = 11)

panel_c <- ggplot(qc, aes(malignant_STK31_TMM_CPM, KIR3DL1_RNA_detected_fraction, label = sample_name)) +
  geom_point(size = 3.2, color = "#b2182b") +
  geom_text(vjust = -0.8, size = 3.1, check_overlap = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Malignant STK31 TMM-CPM (Primary-8 normalization)", y = "NK KIR3DL1 RNA detection fraction",
       title = "C. STK31 vs KIR3DL1 (paired QC, n=4)") +
  theme_classic(base_size = 11)

pair_figure <- panel_a / (panel_b + panel_c) + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(title = "Exact same-sample malignant epithelial-NK pairing: Primary tumor only")
save_figure(pair_figure, "Primary_same_sample_malignant_NK_pairing", 12, 8.5)

rho_hlab <- correlations[malignant_variable == "malignant_HLA.B_TMM_CPM" & NK_variable == "KIR3DL1_RNA_detected_fraction"]
rho_stk <- correlations[malignant_variable == "malignant_STK31_TMM_CPM" & NK_variable == "KIR3DL1_RNA_detected_fraction"]
summary_lines <- c(
  "# Stage 4A：原发灶同一样本恶性上皮–NK 配对",
  "",
  "## 配对规则",
  "",
  "只纳入 `standardized_site == \"Primary tumor\"`，并要求恶性上皮与 NK 来自完全相同的 `sample_name`；未进行跨部位或跨样本配对。",
  "",
  "## 样本数量",
  "",
  "- 11 个原发灶样本同时存在恶性上皮和 NK，可进入完整描述性配对表。",
  "- 预设结构阈值为恶性上皮和 NK 均 >=100 个细胞，仅 4 个样本合格：GBC_033_P、GBC_047_P、GBC_056_P、GBC_073_P。",
  "- 其余 7 个样本仅作描述，不因 STK31、HLA-B 或 KIR3DL1 表达状态而筛除。",
  "",
  "## 描述性结果",
  "",
  sprintf("- 配对 QC 样本中，恶性 HLA-B TMM-CPM 与 NK KIR3DL1 检出率：Spearman rho=%.3f，p=%.4f，n=4。", rho_hlab$spearman_rho, rho_hlab$p_value),
  sprintf("- 配对 QC 样本中，恶性 STK31 TMM-CPM 与 NK KIR3DL1 检出率：Spearman rho=%.3f，p=%.4f，n=4。", rho_stk$spearman_rho, rho_stk$p_value),
  "- n=4 极小，且 KIR3DL1 scRNA 检出稀疏；上述数值只能用于透明描述，不能验证 STK31–HLA-B–KIR3DL1 机制。",
  "",
  "## 使用建议",
  "",
  "该配对结果最多作为补充材料中的可行性/样本结构图和表，不建议放入正文作为关联或因果证据。"
)
writeLines(summary_lines, file.path(report_dir, "stage4A_primary_same_sample_malignant_nk_summary.md"), useBytes = TRUE)
message("STAGE4A_PRIMARY_SAME_SAMPLE_MALIGNANT_NK_OK exact_pairs=11 paired_qc=4")
