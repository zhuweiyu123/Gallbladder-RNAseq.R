#!/usr/bin/env Rscript

# Primary-tumor-only companion version of Stage 4A Figure 1.
# Uses author annotations and a streaming Matrix Market scan; no Seurat object
# or dense all-cell matrix is created.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
scripts_dir <- file.path(project, "04_scripts")
fig_dir <- file.path(project, "11_figures", "current_paper_minimal")
table_dir <- file.path(project, "12_tables", "current_paper_minimal")
technical_dir <- file.path(table_dir, "technical")
log_dir <- file.path(project, "logs")
for (path in c(fig_dir, table_dir, technical_dir, log_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
run_checked <- function(command, args, log_file) {
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file, useBytes = TRUE)
  status <- attr(output, "status")
  stop_if_not(is.null(status) || status == 0L, paste("Command failed:", command))
}
read_gz_tsv <- function(path, ...) fread(cmd = paste("gzip -cd", shQuote(path)), showProgress = FALSE, ...)
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
}

major_levels <- c(
  "Normal epithelial cells", "Malignant epithelial cells", "T cells", "NK cells", "B cells",
  "Myeloid cells", "Fibroblasts", "Endothelial cells"
)
major_labels <- c(
  "Normal\nepithelial", "Malignant\nepithelial", "T cells", "NK cells", "B cells",
  "Myeloid cells", "Fibroblasts", "Endothelial\ncells"
)
metadata_path <- file.path(project, "02_processed_data", "All_single_cells", "GBC_Metadata.txt")
barcodes_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "barcodes.tsv.gz")
features_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "features.tsv.gz")
matrix_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "matrix.mtx.gz")
malignant_sample_path <- file.path(project, "06_feasibility", "stage3", "malignant_sample_gene_detection.csv")

metadata <- fread(metadata_path, showProgress = FALSE)[name != "type"]
stop_if_not(nrow(metadata) == 1117245L, "Unexpected all-cell metadata row count")
stop_if_not(identical(names(metadata), c("name", "celltype", "subtype", "sample_name")), "Unexpected all-cell metadata schema")
barcodes <- read_gz_tsv(barcodes_path, header = FALSE, col.names = "name")
stop_if_not(identical(metadata$name, barcodes$name), "Metadata and barcode order differ")
features <- read_gz_tsv(features_path, header = FALSE)
stk31_row <- which(features[[2]] == "STK31")
stop_if_not(length(stk31_row) == 1L, "STK31 feature row is not unique")

sample_sites <- unique(fread(malignant_sample_path, showProgress = FALSE)[, .(sample_name, standardized_site)])
stop_if_not(!anyDuplicated(sample_sites$sample_name), "Duplicate sample-site mapping")
primary_samples <- sample_sites[standardized_site == "Primary tumor", sample_name]
stop_if_not(length(primary_samples) == 11L, "Expected 11 primary tumor malignant-cell samples")
stop_if_not(all(primary_samples %in% metadata$sample_name), "A primary sample is absent from full metadata")

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
metadata[is.na(group_index) | !sample_name %in% primary_samples, group_index := 0L]

cell_counts <- metadata[group_index > 0L, .N, by = .(group_index, major_celltype)]
setorder(cell_counts, group_index)
stop_if_not(nrow(cell_counts) == length(major_levels), "A requested major primary-tumor cell type has no cells")
fwrite(cell_counts, file.path(table_dir, "Figure1_primary_tumor_STK31_major_celltype_cell_numbers.csv"))

map_path <- file.path(technical_dir, "stage4A_primary_full_to_major_type.int32.bin")
writeBin(as.integer(metadata$group_index), map_path, size = 4L, endian = .Platform$endian)
rm(metadata, barcodes)
gc()

cpp <- file.path(scripts_dir, "stage4A_extract_stk31_major_types.cpp")
bin <- file.path(technical_dir, "stage4A_extract_stk31_primary_major_types")
stats_path <- file.path(table_dir, "Figure1_primary_tumor_STK31_major_celltype_stats.csv")
needs_stream <- !file.exists(stats_path) || nrow(fread(stats_path, showProgress = FALSE)) != length(major_levels)
if (needs_stream) {
  run_checked("g++", c("-O3", "-std=c++17", cpp, "-lz", "-o", bin), file.path(log_dir, "stage4A_compile_primary_stk31_stream.log"))
  run_checked(
    "/usr/bin/time",
    c("-v", "-o", file.path(log_dir, "stage4A_primary_stk31_major_type_stream.time.txt"), bin,
      matrix_path, map_path, "32137", "1117245", as.character(stk31_row),
      as.character(length(major_levels)), stats_path),
    file.path(log_dir, "stage4A_primary_stk31_major_type_stream.log")
  )
}

stats <- fread(stats_path)
stop_if_not(nrow(stats) == length(major_levels), "Unexpected primary STK31 statistic rows")
stop_if_not(all(stats$declared_nnz == stats$actual_nnz), "Primary stream did not consume declared nnz")
if (!"major_celltype" %in% names(stats)) {
  stats <- merge(cell_counts, stats, by = "group_index", sort = FALSE)
  setorder(stats, group_index)
  stop_if_not(all(stats$N == cell_counts$N), "Primary cell count merge failed")
  setnames(stats, "N", "cell_n")
  stats[, major_celltype := factor(major_celltype, levels = major_levels)]
  stats[, STK31_RNA_detected_fraction := STK31_detected_cell_n / cell_n]
  stats[, STK31_pseudobulk_CPM := STK31_raw_UMI / group_library_UMI * 1e6]
  stats[, log1p_STK31_pseudobulk_CPM := log1p(STK31_pseudobulk_CPM)]
} else {
  stop_if_not(identical(as.integer(stats$cell_n), as.integer(cell_counts$N)), "Reused primary statistics differ from metadata cell counts")
  stats[, major_celltype := factor(major_celltype, levels = major_levels)]
}
fwrite(stats, stats_path)

fig_a <- ggplot(stats, aes(x = major_celltype, y = "STK31")) +
  geom_point(aes(size = STK31_RNA_detected_fraction, color = log1p_STK31_pseudobulk_CPM)) +
  scale_size_continuous(name = "% RNA detected", range = c(2.5, 13), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  scale_x_discrete(labels = major_labels) +
  labs(x = NULL, y = NULL, title = "STK31 across major cell types: Primary tumor only") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.ticks.y = element_blank(), axis.text.y = element_blank())

malignant <- stats[major_celltype == "Malignant epithelial cells"]
fig_b <- ggplot(malignant, aes(x = "Malignant epithelial cells", y = STK31_RNA_detected_fraction)) +
  geom_col(width = 0.55, fill = "#8c2d04") +
  geom_text(aes(label = paste0(format(STK31_detected_cell_n, big.mark = ","), " / ", format(cell_n, big.mark = ","), " cells\n", fmt_pct(STK31_RNA_detected_fraction, 2))), vjust = -0.35, size = 3.4) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "STK31 RNA detection fraction", title = "STK31 in primary malignant cells") +
  theme_classic(base_size = 11)

figure <- fig_a + fig_b + plot_layout(widths = c(1.7, 1)) + plot_annotation(tag_levels = "A")
save_figure(figure, "Figure1_primary_tumor_STK31_celltype_localization", width = 11, height = 4.8)
message("STAGE4A_PRIMARY_FIGURE1_OK samples=", length(primary_samples), " malignant_cells=", malignant$cell_n, " STK31_detected=", malignant$STK31_detected_cell_n)
