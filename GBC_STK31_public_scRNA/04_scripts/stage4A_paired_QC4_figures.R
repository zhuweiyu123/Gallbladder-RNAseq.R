#!/usr/bin/env Rscript

# Rebuild the four minimal-paper public-scRNA figures using only the four
# exact same-primary-sample malignant epithelial-NK pairs with >=100 cells
# in both compartments. Author annotations are retained; no reclustering.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
scripts_dir <- file.path(project, "04_scripts")
fig_dir <- file.path(project, "11_figures", "current_paper_minimal", "paired_QC4")
table_dir <- file.path(project, "12_tables", "current_paper_minimal", "paired_QC4")
technical_dir <- file.path(table_dir, "technical")
report_dir <- file.path(project, "06_feasibility")
log_dir <- file.path(project, "logs")
for (path in c(fig_dir, table_dir, technical_dir, report_dir, log_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_gz_tsv <- function(path, ...) fread(cmd = paste("gzip -cd", shQuote(path)), showProgress = FALSE, ...)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
}
run_checked <- function(command, args, log_file) {
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file, useBytes = TRUE)
  status <- attr(output, "status")
  stop_if_not(is.null(status) || status == 0L, paste("Command failed:", command))
}

qc4 <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")
pair_table <- fread(file.path(project, "12_tables", "current_paper_minimal", "Primary_same_sample_malignant_NK_pairs.csv"))
stop_if_not(identical(pair_table[pass_paired_QC == TRUE, sample_name], qc4), "QC4 sample set differs from the locked exact-pair set")

major_levels <- c(
  "Normal epithelial cells", "Malignant epithelial cells", "T cells", "NK cells",
  "B cells", "Myeloid cells", "Fibroblasts", "Endothelial cells"
)
major_labels <- c(
  "Normal\nepithelial", "Malignant\nepithelial", "T cells", "NK cells",
  "B cells", "Myeloid cells", "Fibroblasts", "Endothelial\ncells"
)
hla_main <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRD1")

# Figure 1: full-cell STK31 localisation in QC4 only, via low-memory stream.
metadata_path <- file.path(project, "02_processed_data", "All_single_cells", "GBC_Metadata.txt")
barcodes_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "barcodes.tsv.gz")
features_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "features.tsv.gz")
matrix_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "matrix.mtx.gz")
metadata <- fread(metadata_path, showProgress = FALSE)[name != "type"]
barcodes <- read_gz_tsv(barcodes_path, header = FALSE, col.names = "name")
features <- read_gz_tsv(features_path, header = FALSE)
stop_if_not(nrow(metadata) == 1117245L && identical(metadata$name, barcodes$name), "Full metadata-barcode mapping failed")
stk31_row <- which(features[[2]] == "STK31")
stop_if_not(length(stk31_row) == 1L, "STK31 feature must be unique")

metadata[, major_celltype := NA_character_]
metadata[subtype == "Normal epithelial cells", major_celltype := "Normal epithelial cells"]
metadata[subtype == "Malignant epithelial cells", major_celltype := "Malignant epithelial cells"]
metadata[celltype %in% c("CD8+ T cell", "CD4+ T cells"), major_celltype := "T cells"]
metadata[celltype == "NK cells", major_celltype := "NK cells"]
metadata[celltype == "B cells", major_celltype := "B cells"]
metadata[celltype %in% c("Monocytes & Macrophages", "Neutrophils", "Dendritic cells"), major_celltype := "Myeloid cells"]
metadata[celltype == "Mesenchymal cells" & grepl("^F_", subtype), major_celltype := "Fibroblasts"]
metadata[celltype == "Endothelial cells", major_celltype := "Endothelial cells"]
metadata[, group_index := match(major_celltype, major_levels)]
metadata[is.na(group_index) | !sample_name %in% qc4, group_index := 0L]
cell_counts <- metadata[group_index > 0L, .(cell_n = .N), by = .(group_index, major_celltype)]
setorder(cell_counts, group_index)
stop_if_not(nrow(cell_counts) == length(major_levels), "A requested QC4 major cell type has zero cells")
fwrite(cell_counts, file.path(table_dir, "QC4_Figure1_major_celltype_cell_numbers.csv"))

map_path <- file.path(technical_dir, "QC4_full_to_major_type.int32.bin")
writeBin(as.integer(metadata$group_index), map_path, size = 4L, endian = .Platform$endian)
rm(metadata, barcodes)
gc()
cpp <- file.path(scripts_dir, "stage4A_extract_stk31_major_types.cpp")
bin <- file.path(technical_dir, "QC4_extract_stk31_major_types")
stats_path <- file.path(table_dir, "QC4_Figure1_STK31_major_celltype_stats.csv")
needs_stream <- !file.exists(stats_path) || nrow(fread(stats_path, showProgress = FALSE)) != length(major_levels)
if (needs_stream) {
  run_checked("g++", c("-O3", "-std=c++17", cpp, "-lz", "-o", bin), file.path(log_dir, "stage4A_QC4_compile_stk31.log"))
  run_checked(
    "/usr/bin/time",
    c("-v", "-o", file.path(log_dir, "stage4A_QC4_stk31_stream.time.txt"), bin,
      matrix_path, map_path, "32137", "1117245", as.character(stk31_row), as.character(length(major_levels)), stats_path),
    file.path(log_dir, "stage4A_QC4_stk31_stream.log")
  )
}
stats <- fread(stats_path)
stop_if_not(nrow(stats) == length(major_levels) && all(stats$declared_nnz == stats$actual_nnz), "QC4 STK31 stream failed integrity check")
if (!"major_celltype" %in% names(stats)) {
  stats <- merge(cell_counts, stats, by = "group_index", sort = FALSE)
  setorder(stats, group_index)
  stats[, STK31_RNA_detected_fraction := STK31_detected_cell_n / cell_n]
  stats[, STK31_pseudobulk_CPM := STK31_raw_UMI / group_library_UMI * 1e6]
  stats[, log1p_STK31_pseudobulk_CPM := log1p(STK31_pseudobulk_CPM)]
} else {
  stop_if_not(all(stats$cell_n == cell_counts$cell_n), "Reused QC4 Figure 1 cell counts differ")
}
stats[, major_celltype := factor(major_celltype, levels = major_levels)]
fwrite(stats, stats_path)

fig1a <- ggplot(stats, aes(major_celltype, "STK31")) +
  geom_point(aes(size = STK31_RNA_detected_fraction, color = log1p_STK31_pseudobulk_CPM)) +
  scale_x_discrete(labels = major_labels) +
  scale_size_continuous(name = "% RNA detected", range = c(2.5, 13), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "A. STK31 across major cell types (paired QC4)") +
  theme_classic(base_size = 10.5) +
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())
mal_stk <- stats[major_celltype == "Malignant epithelial cells"]
fig1b <- ggplot(mal_stk, aes("Malignant epithelial", STK31_RNA_detected_fraction)) +
  geom_col(width = 0.55, fill = "#8c2d04") +
  geom_text(aes(label = paste0(STK31_detected_cell_n, " / ", cell_n, " cells\n", sprintf("%.2f%%", 100 * STK31_RNA_detected_fraction))), vjust = -0.35, size = 3.4) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "STK31 RNA detection fraction", title = "B. STK31 in malignant epithelial cells") +
  theme_classic(base_size = 10.5)
fig1 <- fig1a + fig1b + plot_layout(widths = c(1.8, 1))
save_figure(fig1, "QC4_Figure1_STK31_celltype_localization", 12, 4.8)

# Figure 2: HLA-I in QC4 malignant epithelial cells.
target_counts <- readRDS(file.path(project, "02_processed_data", "malignant_epithelial", "target_counts", "malignant_target_raw_counts.rds"))
mal_meta <- fread(file.path(project, "02_processed_data", "malignant_epithelial", "malignant_cell_metadata.csv.gz"), showProgress = FALSE)
mal_library <- fread(file.path(project, "02_processed_data", "malignant_epithelial", "malignant_cell_library_sizes.tsv"), showProgress = FALSE)
stop_if_not(identical(colnames(target_counts), mal_meta$barcode), "Malignant target matrix and metadata order differ")
stop_if_not(identical(as.integer(mal_library$malignant_col), as.integer(mal_meta$malignant_col)), "Malignant library rows differ from metadata")
mal_idx <- which(mal_meta$sample_name %in% qc4)
stop_if_not(length(mal_idx) == 13368L, "QC4 malignant epithelial cell total must be 13,368")
qc_target <- target_counts[, mal_idx, drop = FALSE]
qc_library <- mal_library$all_gene_raw_umi[mal_idx]
hla_stats <- rbindlist(lapply(hla_main, function(gene) {
  values <- as.numeric(qc_target[gene, ])
  data.table(
    gene = gene, malignant_cell_n = length(values), RNA_detected_n = sum(values > 0),
    RNA_detected_fraction = mean(values > 0), total_raw_UMI = sum(values),
    pseudobulk_CPM = sum(values) / sum(qc_library) * 1e6
  )
}))
hla_stats[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
hla_stats[, gene := factor(gene, levels = hla_main)]
fwrite(hla_stats, file.path(table_dir, "QC4_Figure2_malignant_HLAI_stats.csv"))
fig2 <- ggplot(hla_stats, aes(gene, "Malignant epithelial")) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(name = "% RNA detected", range = c(5, 16), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL, title = "HLA-I RNA in paired-QC4 malignant epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_figure(fig2, "QC4_Figure2_malignant_epithelial_HLAI", 6.4, 4.8)
rm(target_counts, qc_target)
gc()

# Figure 3: NK inhibitory receptors by author subtype in QC4 only.
nk <- readRDS(file.path(project, "01_raw_data", "NK.RDS"))
stop_if_not(inherits(nk, "Seurat") && DefaultAssay(nk) == "RNA", "NK.RDS is not the expected RNA Seurat object")
nk_counts_all <- LayerData(nk[["RNA"]], layer = "counts")
nk_meta_all <- as.data.table(nk[[]], keep.rownames = "barcode")
stop_if_not(identical(colnames(nk_counts_all), nk_meta_all$barcode), "NK counts and metadata order differ")
nk_idx <- which(nk_meta_all$orig.ident %in% qc4)
stop_if_not(length(nk_idx) == 1135L, "QC4 NK cell total must be 1,135")
nk_counts <- nk_counts_all[, nk_idx, drop = FALSE]
nk_meta <- nk_meta_all[nk_idx]
nk_meta[, `:=`(nk_subtype = as.character(celltype), all_gene_raw_umi = as.numeric(Matrix::colSums(nk_counts)))]
stop_if_not(all(nk_receptors %in% rownames(nk_counts)) && all(nk_meta$all_gene_raw_umi > 0), "QC4 NK target features/library sizes invalid")
subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^.*_C([0-9]+)_.*$", "\\1", nk_subtype))]
setorder(subtype_sizes, cluster_number, nk_subtype)
nk_dot <- rbindlist(lapply(nk_receptors, function(gene) {
  working <- nk_meta[, .(nk_subtype, all_gene_raw_umi)]
  working[, raw_UMI := as.numeric(nk_counts[gene, ])]
  working[, .(
    NK_cell_n = .N, RNA_detected_n = sum(raw_UMI > 0), RNA_detected_fraction = mean(raw_UMI > 0),
    total_raw_UMI = sum(raw_UMI), mean_logCP10K = mean(log1p(raw_UMI / all_gene_raw_umi * 1e4))
  ), by = nk_subtype][, gene := gene]
}))
nk_dot <- merge(nk_dot, subtype_sizes[, .(nk_subtype, NK_cell_n, cluster_number)], by = c("nk_subtype", "NK_cell_n"), all.x = TRUE, sort = FALSE)
nk_dot[, `:=`(
  gene = factor(gene, levels = rev(nk_receptors)),
  subtype_display = paste0(nk_subtype, "\n(n=", NK_cell_n, ")")
)]
subtype_levels <- nk_dot[gene == "KLRC1"][order(cluster_number, nk_subtype), subtype_display]
nk_dot[, subtype_display := factor(subtype_display, levels = subtype_levels)]
fwrite(nk_dot, file.path(table_dir, "QC4_Figure3_NK_receptor_subtype_stats.csv"))
fig3 <- ggplot(nk_dot, aes(subtype_display, gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = mean_logCP10K)) +
  scale_size_continuous(name = "% RNA detected", range = c(1.4, 10), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "mean log1p(CP10K)", low = "#e0f3db", high = "#08589e") +
  labs(x = "Author NK/NKT subtype (QC4 cell number shown)", y = NULL, title = "NK inhibitory receptors in paired-QC4 primary samples") +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1), axis.ticks = element_blank())
fig3_width <- max(10, 0.85 * uniqueN(nk_dot$nk_subtype))
save_figure(fig3, "QC4_Figure3_NK_inhibitory_receptor_subtypes", fig3_width, 5.8)
rm(nk, nk_counts_all, nk_counts, nk_meta_all, nk_meta)
gc()

# Figure 4: sample-level STK31-HLA-I trends, QC4 only.
pb <- fread(file.path(project, "12_tables", "current_paper_minimal", "Primary8_STK31_HLAI_TMM_CPM.csv"))
pb <- pb[sample_name %in% qc4]
setorder(pb, sample_name)
stop_if_not(identical(pb$sample_name, qc4), "QC4 Primary-8 TMM subset differs from the locked sample order")
trend_cor <- rbindlist(lapply(hla_main, function(gene) {
  test <- suppressWarnings(cor.test(pb$STK31_log2CPM, pb[[paste0(gene, "_log2CPM")]], method = "spearman", exact = FALSE))
  data.table(gene = gene, n = 4L, spearman_rho = unname(test$estimate), p_value = test$p.value,
             interpretation = "descriptive_n4_not_for_inference")
}))
fwrite(trend_cor, file.path(table_dir, "QC4_Figure4_STK31_HLAI_correlations.csv"))
trend_long <- rbindlist(lapply(hla_main, function(gene) data.table(
  sample_name = pb$sample_name, gene = gene, STK31_log2CPM = pb$STK31_log2CPM,
  HLA_log2CPM = pb[[paste0(gene, "_log2CPM")]]
)))
annotation <- copy(trend_cor)
annotation[, facet_label := sprintf("%s\nn=4, rho=%.2f, p=%.3g", gene, spearman_rho, p_value)]
trend_long <- merge(trend_long, annotation[, .(gene, facet_label)], by = "gene", all.x = TRUE, sort = FALSE)
annotation[, facet_label := factor(facet_label, levels = annotation[match(hla_main, gene), facet_label])]
trend_long[, facet_label := factor(facet_label, levels = levels(annotation$facet_label))]
fig4 <- ggplot(trend_long, aes(STK31_log2CPM, HLA_log2CPM, label = sample_name)) +
  geom_point(size = 3, color = "#1f78b4") +
  geom_text(vjust = -0.75, size = 3, check_overlap = TRUE) +
  facet_wrap(~facet_label, scales = "free_y", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.18))) +
  labs(x = "Malignant STK31 log2 TMM-CPM", y = "HLA-I log2 TMM-CPM",
       title = "Paired-QC4 primary samples: descriptive STK31-HLA-I trends") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", color = "grey70"))
save_figure(fig4, "QC4_Figure4_STK31_HLAI_pseudobulk_trends", 12, 3.8)

summary_lines <- c(
  "# Stage 4A：4 个同一样本配对 QC 原发灶样本重绘",
  "",
  "仅使用 GBC_033_P、GBC_047_P、GBC_056_P、GBC_073_P。四个样本均满足恶性上皮和 NK 各 >=100 个细胞，且两类细胞来自完全相同的 Primary tumor `sample_name`。",
  "",
  paste0("- QC4 恶性上皮细胞总数：", format(length(mal_idx), big.mark = ","), "。"),
  paste0("- QC4 NK 细胞总数：1,135；作者 NK/NKT 亚群数：", uniqueN(nk_dot$nk_subtype), "。"),
  paste0("- QC4 恶性上皮 STK31 检出：", mal_stk$STK31_detected_cell_n, "/", mal_stk$cell_n,
         "（", sprintf("%.2f%%", 100 * mal_stk$STK31_RNA_detected_fraction), "）。"),
  "- Figure 1–3 可用于观察 QC4 子集中的定位和检出结构，但不能代表整个原发胆囊癌队列。",
  "- Figure 4 仅有 n=4，任何 Spearman 结果均为描述性，不得用于关联或因果结论。",
  "",
  "输出位于 `11_figures/current_paper_minimal/paired_QC4/` 和 `12_tables/current_paper_minimal/paired_QC4/`。"
)
writeLines(summary_lines, file.path(report_dir, "stage4A_paired_QC4_figures_summary.md"), useBytes = TRUE)
message("STAGE4A_PAIRED_QC4_FIGURES_OK samples=4 malignant=13368 NK=1135")
