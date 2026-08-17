#!/usr/bin/env Rscript

# Stage 4A: primary-tumor-only malignant epithelial pseudobulk TMM.
# Descriptive validation only. TMM is independently estimated within the
# predefined Primary QC (>=100 malignant cells) and Primary all (11 samples)
# sets, using the full transcriptome pseudobulk count matrix.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
fig_dir <- file.path(project, "11_figures", "current_paper_minimal")
table_dir <- file.path(project, "12_tables", "current_paper_minimal")
technical_dir <- file.path(table_dir, "technical")
report_dir <- file.path(project, "06_feasibility")
for (path in c(fig_dir, table_dir, report_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
fmt_pct <- function(x, digits = 2L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
}

target_genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M")
primary_expected <- c(
  "GBC_005_P", "GBC_023_P", "GBC_033_P", "GBC_034_P", "GBC_038_P",
  "GBC_042_P", "GBC_047_P", "GBC_056_P", "GBC_073_P", "GBC_075_P", "GBC_086_P"
)
primary_qc_expected <- c(
  "GBC_005_P", "GBC_033_P", "GBC_034_P", "GBC_042_P",
  "GBC_047_P", "GBC_056_P", "GBC_073_P", "GBC_086_P"
)

matrix_path <- file.path(technical_dir, "stage4A_malignant_pseudobulk_counts.mtx.gz")
audit_path <- file.path(technical_dir, "stage4A_malignant_pseudobulk_stream_audit.tsv")
features_path <- file.path(project, "02_processed_data", "malignant_epithelial", "10X_counts", "features.tsv.gz")
sample_meta_path <- file.path(table_dir, "Figure4_malignant_pseudobulk_sample_metadata.csv")
stage3_detection_path <- file.path(project, "06_feasibility", "stage3", "malignant_sample_gene_detection.csv")
old58_path <- file.path(table_dir, "Figure4_malignant_TMM_pseudobulk_target_genes.csv")

stop_if_not(file.exists(matrix_path) && file.exists(features_path) && file.exists(sample_meta_path), "Required Stage 4A pseudobulk inputs are absent")
audit <- fread(audit_path, showProgress = FALSE)
stop_if_not(nrow(audit) == 1L && audit$output_rows == 32137L && audit$output_cols == 58L && audit$declared_nnz == audit$actual_nnz,
            "Pseudobulk stream audit is inconsistent")

sample_meta <- fread(sample_meta_path, showProgress = FALSE)
stop_if_not(nrow(sample_meta) == 58L && !anyDuplicated(sample_meta$sample_name), "Expected 58 unique pseudobulk sample metadata rows")

features <- fread(cmd = paste("gzip -cd", shQuote(features_path)), header = FALSE, showProgress = FALSE)
stop_if_not(nrow(features) == 32137L && ncol(features) >= 2L, "Unexpected feature table shape")
gene_symbol <- as.character(features[[2]])
target_rows <- vapply(target_genes, function(gene) {
  hit <- which(gene_symbol == gene)
  stop_if_not(length(hit) == 1L, paste("Target gene must occur exactly once:", gene))
  hit
}, integer(1))

counts_sparse <- readMM(gzfile(matrix_path))
counts_sparse <- as(counts_sparse, "CsparseMatrix")
stop_if_not(identical(dim(counts_sparse), c(32137L, 58L)), "Pseudobulk matrix must be 32137 genes x 58 samples")
stop_if_not(all(is.finite(counts_sparse@x)) && all(counts_sparse@x >= 0) && all(abs(counts_sparse@x - round(counts_sparse@x)) < 1e-8),
            "Pseudobulk counts must be finite nonnegative integers")
counts <- as.matrix(counts_sparse)
rm(counts_sparse)
storage.mode(counts) <- "double"
rownames(counts) <- make.unique(gene_symbol)
colnames(counts) <- sample_meta$sample_name
stop_if_not(identical(colnames(counts), sample_meta$sample_name), "Pseudobulk column order differs from sample metadata")

stage3 <- fread(stage3_detection_path, showProgress = FALSE)
library_manifest <- unique(stage3[, .(sample_name, malignant_pseudobulk_library_size)])
stop_if_not(nrow(library_manifest) == 58L && !anyDuplicated(library_manifest$sample_name), "Invalid Stage 3 malignant pseudobulk library manifest")
expected_library_sizes <- library_manifest[match(sample_meta$sample_name, library_manifest$sample_name), malignant_pseudobulk_library_size]
stop_if_not(!anyNA(expected_library_sizes) && all(colSums(counts) == expected_library_sizes),
            "Pseudobulk matrix columns do not match the saved 58-row sample metadata order")
primary_detail <- stage3[
  standardized_site == "Primary tumor" & gene == "STK31",
  .(sample_name, patient_id, patient_id_source, TNM_stage, histological_type,
    malignant_epithelial_cell_n = malignant_cell_n,
    STK31_RNA_detected_cell_n = RNA_detected_n,
    STK31_RNA_detected_fraction = RNA_detected_fraction,
    STK31_raw_UMI = pseudobulk_raw_count)
]
setorder(primary_detail, sample_name)
stop_if_not(identical(primary_detail$sample_name, primary_expected), "Primary sample set/order differs from prespecified list")
stop_if_not(sum(primary_detail$malignant_epithelial_cell_n) == 36910L, "Primary malignant-cell total must be 36,910")
primary_detail[, pass_malignant_cell_threshold := malignant_epithelial_cell_n >= 100L]
stop_if_not(identical(primary_detail[pass_malignant_cell_threshold == TRUE, sample_name], primary_qc_expected), "Primary QC set differs from prespecified samples")

run_tmm <- function(sample_names, label) {
  stop_if_not(all(sample_names %in% colnames(counts)), paste("Missing samples for", label))
  counts_use <- counts[, sample_names, drop = FALSE]
  keep_gene <- rowSums(counts_use) > 0
  stop_if_not(sum(keep_gene) > 0L, paste("No nonzero genes in", label))
  dge <- DGEList(counts = counts_use[keep_gene, , drop = FALSE])
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_cpm <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE, prior.count = 0)
  log2_cpm <- cpm(dge, log = TRUE, normalized.lib.sizes = TRUE, prior.count = 1)
  factors <- as.data.table(dge$samples, keep.rownames = "sample_name")
  setnames(factors, "lib.size", "raw_library_size")
  factors[, effective_library_size := raw_library_size * norm.factors]
  factors[, analysis_set := label]
  stop_if_not(all(is.finite(factors$norm.factors) & factors$norm.factors > 0), paste("Invalid TMM factors in", label))
  stop_if_not(all(is.finite(factors$effective_library_size) & factors$effective_library_size > 0), paste("Invalid effective library size in", label))
  target_keys <- make.unique(gene_symbol)[target_rows]
  stop_if_not(all(target_keys %in% rownames(tmm_cpm)), paste("Target gene lost in", label))
  for (gene in target_genes) {
    key <- make.unique(gene_symbol)[target_rows[[gene]]]
    zero_raw <- counts[key, sample_names] == 0
    stop_if_not(all(tmm_cpm[key, sample_names][zero_raw] == 0), paste("Zero raw counts have nonzero TMM CPM for", gene, label))
  }
  list(tmm_cpm = tmm_cpm, log2_cpm = log2_cpm, factors = factors, keep_gene_n = sum(keep_gene))
}

primary11_names <- primary_detail$sample_name
primary8_names <- primary_detail[pass_malignant_cell_threshold == TRUE, sample_name]
tmm11 <- run_tmm(primary11_names, "Primary_all_11")
tmm8 <- run_tmm(primary8_names, "Primary_QC_8")

target_value_table <- function(result, sample_names, prefix = "") {
  out <- data.table(sample_name = sample_names)
  for (gene in target_genes) {
    key <- make.unique(gene_symbol)[target_rows[[gene]]]
    out[, (paste0(gene, "_raw_count")) := as.numeric(counts[key, sample_names])]
    out[, (paste0(gene, "_TMM_CPM")) := as.numeric(result$tmm_cpm[key, sample_names])]
    out[, (paste0(gene, "_log2CPM")) := as.numeric(result$log2_cpm[key, sample_names])]
  }
  out
}

primary11_target <- target_value_table(tmm11, primary11_names)
primary8_target <- target_value_table(tmm8, primary8_names)
primary11_descriptive <- merge(primary_detail, primary11_target[, .(sample_name, STK31_primary11_raw_count = STK31_raw_count,
  STK31_primary11_TMM_CPM = STK31_TMM_CPM, STK31_primary11_log2CPM = STK31_log2CPM)], by = "sample_name", all.x = TRUE, sort = FALSE)
setorder(primary11_descriptive, sample_name)
stop_if_not(all(primary11_descriptive$STK31_raw_UMI == primary11_descriptive$STK31_primary11_raw_count), "Stage 3 and pseudobulk STK31 raw counts disagree")
fwrite(primary11_descriptive, file.path(table_dir, "Primary11_STK31_descriptive.csv"))

primary8_factors <- merge(primary_detail[pass_malignant_cell_threshold == TRUE, .(sample_name, patient_id, malignant_epithelial_cell_n)],
                          tmm8$factors[, .(sample_name, raw_library_size, norm.factors, effective_library_size, analysis_set)],
                          by = "sample_name", all.x = TRUE, sort = FALSE)
setorder(primary8_factors, sample_name)
fwrite(primary8_factors, file.path(table_dir, "Primary8_TMM_factors.csv"))

primary8_results <- merge(primary_detail[pass_malignant_cell_threshold == TRUE], primary8_target, by = "sample_name", all.x = TRUE, sort = FALSE)
primary8_results <- merge(primary8_results, primary8_factors[, .(sample_name, norm.factors, effective_library_size)], by = "sample_name", all.x = TRUE, sort = FALSE)
setorder(primary8_results, sample_name)
fwrite(primary8_results, file.path(table_dir, "Primary8_STK31_HLAI_TMM_CPM.csv"))

analysis_sets <- list(Primary_QC_8 = primary8_names, Primary_all_11_sensitivity = primary11_names)
correlations <- rbindlist(lapply(names(analysis_sets), function(set_name) {
  samples <- analysis_sets[[set_name]]
  values <- if (set_name == "Primary_QC_8") primary8_target else primary11_target
  rbindlist(lapply(target_genes[-1], function(gene) {
    test <- suppressWarnings(cor.test(values$STK31_log2CPM, values[[paste0(gene, "_log2CPM")]], method = "spearman", exact = FALSE))
    data.table(
      analysis_set = set_name,
      gene = gene,
      n = nrow(values),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      STK31_raw_count_nonzero_sample_n = sum(values$STK31_raw_count > 0),
      analysis_label = "exploratory_small_n_noncausal"
    )
  }))
}))
fwrite(correlations, file.path(table_dir, "Primary_only_STK31_HLAI_correlations.csv"))

# Compare re-estimated Primary-8 values with the previously stored all-58 TMM results.
old58 <- fread(old58_path, showProgress = FALSE)
old58_primary8 <- old58[sample_name %in% primary8_names, .(sample_name, STK31_all58_TMM_CPM = STK31_TMM_CPM,
  HLA.A_all58_TMM_CPM = `HLA-A_TMM_CPM`, HLA.B_all58_TMM_CPM = `HLA-B_TMM_CPM`,
  HLA.C_all58_TMM_CPM = `HLA-C_TMM_CPM`, B2M_all58_TMM_CPM = B2M_TMM_CPM)]
comparison <- merge(primary8_results[, .(sample_name, STK31_primary8_TMM_CPM = STK31_TMM_CPM,
  HLA.A_primary8_TMM_CPM = `HLA-A_TMM_CPM`, HLA.B_primary8_TMM_CPM = `HLA-B_TMM_CPM`,
  HLA.C_primary8_TMM_CPM = `HLA-C_TMM_CPM`, B2M_primary8_TMM_CPM = B2M_TMM_CPM)], old58_primary8, by = "sample_name", all.x = TRUE)
stop_if_not(nrow(comparison) == 8L && !anyNA(comparison), "Could not compare Primary-8 and all-58 TMM values")
fwrite(comparison, file.path(table_dir, "Primary8_vs_all58_TMM_comparison.csv"))

# Descriptive plot: retains all 11 samples and preserves sample-name order.
plot_data <- copy(primary11_descriptive)
plot_data[, threshold_label := fifelse(pass_malignant_cell_threshold, ">=100 malignant cells", "<100 cells: descriptive only")]
plot_data[, sample_name := factor(sample_name, levels = primary_expected)]
palette <- c(">=100 malignant cells" = "#1f78b4", "<100 cells: descriptive only" = "#9e9e9e")
panel_a <- ggplot(plot_data, aes(x = sample_name, y = STK31_RNA_detected_fraction, fill = threshold_label)) +
  geom_col(width = 0.72, color = "white") +
  geom_text(aes(label = paste0(STK31_RNA_detected_cell_n, " / ", malignant_epithelial_cell_n)), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "STK31 RNA detection fraction", title = "A. STK31 detection in malignant epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
panel_b <- ggplot(plot_data, aes(x = sample_name, y = STK31_primary11_TMM_CPM, color = threshold_label, shape = threshold_label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
  geom_point(size = 3.3) +
  scale_color_manual(values = palette, name = NULL) +
  scale_shape_manual(values = c(">=100 malignant cells" = 16, "<100 cells: descriptive only" = 1), name = NULL) +
  labs(x = NULL, y = "STK31 TMM-CPM (Primary-all-11 normalization)", title = "B. Sample-level STK31 pseudobulk expression") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
sample_plot <- panel_a / panel_b + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Primary tumor malignant epithelial cells: descriptive STK31 by sample")
save_figure(sample_plot, "Primary_STK31_by_sample", width = 10.5, height = 8)

qc_hlab <- correlations[analysis_set == "Primary_QC_8" & gene == "HLA-B"]
all_hlab <- correlations[analysis_set == "Primary_all_11_sensitivity" & gene == "HLA-B"]
comparison_ratio <- comparison$STK31_primary8_TMM_CPM / comparison$STK31_all58_TMM_CPM
comparison_ratio <- comparison_ratio[is.finite(comparison_ratio) & comparison_ratio > 0]
summary_lines <- c(
  "# Stage 4A：原发灶恶性上皮细胞 pseudobulk–TMM",
  "",
  "## 分析范围",
  "",
  "本任务仅为描述性公共单细胞验证。未进行差异表达、STK31 高/低分组、CellChat、discovery 分支分析，也未重新扫描单细胞表达矩阵。",
  "",
  "## 输入核验",
  "",
  "- 复用恶性上皮细胞的全转录组 pseudobulk 原始计数矩阵：32,137 个基因 × 58 个样本。",
  "- 仅在逐列全转录组 library size 与 Stage 3 中对应样本的 library-size manifest 完全一致后，才将矩阵列映射到保存的 58 行样本元数据顺序；5 个目标基因均在 `features.tsv.gz` 第 2 列中恰好出现一次。",
  "- 矩阵条目均已核验为有限、非负整数；流式提取审计中声明的非零条目数与实际读取数一致。",
  "",
  "## 预设的原发灶样本集",
  "",
  "- Primary all set：11 个样本、共 36,910 个恶性上皮细胞；TMM 仅在这 11 个样本内部重新计算。",
  "- Primary QC set：8 个样本，每个样本 ≥100 个恶性上皮细胞；TMM 仅在这 8 个样本内部独立重新计算。",
  "- 低细胞数、仅作描述性展示的样本：GBC_023_P（7 个细胞）、GBC_038_P（2 个细胞）、GBC_075_P（59 个细胞）。其纳入或排除只由预设细胞数阈值决定，与 STK31 表达无关。",
  "",
  "## 11 个原发灶样本中的 STK31 检出",
  "",
  paste0("STK31 原始计数在 ", sum(primary11_target$STK31_raw_count > 0), "/11 个样本中非零。逐样本计数和 Primary-all-11 TMM 值见 `Primary11_STK31_descriptive.csv` 及配套图。"),
  "原始计数为 0 时，非 log 的 TMM-CPM 必为 0；log2CPM 列仅以 prior.count=1 保证数值稳定，不能解释为检出表达。",
  "",
  "## 探索性 STK31–HLA-B 关联",
  "",
  sprintf("- Primary QC set（n=8）：Spearman rho=%.3f，p=%.4f；STK31 原始计数非零样本=%d/8。", qc_hlab$spearman_rho, qc_hlab$p_value, qc_hlab$STK31_raw_count_nonzero_sample_n),
  sprintf("- Primary all-set 敏感性分析（n=11）：Spearman rho=%.3f，p=%.4f；STK31 原始计数非零样本=%d/11。", all_hlab$spearman_rho, all_hlab$p_value, all_hlab$STK31_raw_count_nonzero_sample_n),
  "- 上述结果仅为小样本探索性描述，不能作为 STK31 调控 HLA-I 的证据。",
  "",
  "## 与既往 58 样本 TMM 的比较",
  "",
  paste0("Primary-8 的 TMM 为重新估计，未沿用混合不同部位的 58 样本计算结果。对 STK31 原始计数非零的 Primary-8 样本，Primary-8/all-58 TMM-CPM 比值范围为 ",
         sprintf("%.3f", min(comparison_ratio)), "–", sprintf("%.3f", max(comparison_ratio)), "；详情见 `Primary8_vs_all58_TMM_comparison.csv`。两种分析的 TMM 参考组成不同，不能互换使用。"),
  "",
  "## 图的定位",
  "",
  "`Primary_STK31_by_sample` 最适合放在补充材料中作为描述性图：它透明呈现零表达和低细胞数样本，但不足以支撑正文的机制性主张。",
  "",
  "Stage 4A 的 Primary-only TMM 分析到此停止。"
)
writeLines(summary_lines, file.path(report_dir, "stage4A_primary_only_tmm_summary.md"), useBytes = TRUE)
message("STAGE4A_PRIMARY_ONLY_TMM_OK primary11=11 primary8=8")
