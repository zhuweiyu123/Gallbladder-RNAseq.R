#!/usr/bin/env Rscript

# Current-paper minimal public scRNA validation.
# Scope: four core figures plus two direct supplementary displays only.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
scripts_dir <- file.path(project, "04_scripts")
fig_dir <- file.path(project, "11_figures", "current_paper_minimal")
table_dir <- file.path(project, "12_tables", "current_paper_minimal")
technical_dir <- file.path(table_dir, "technical")
report_dir <- file.path(project, "06_feasibility")
log_dir <- file.path(project, "logs")
for (path in c(fig_dir, table_dir, technical_dir, report_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
safe_fraction <- function(numerator, denominator) ifelse(denominator > 0, numerator / denominator, NA_real_)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
}
run_checked <- function(command, args, log_file) {
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(output, log_file, useBytes = TRUE)
  status <- attr(output, "status")
  stop_if_not(is.null(status) || status == 0L, paste("Command failed:", command, paste(args, collapse = " ")))
}
read_gz_tsv <- function(path, ...) {
  fread(cmd = paste("gzip -cd", shQuote(path)), showProgress = FALSE, ...)
}
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)

major_levels <- c(
  "Malignant epithelial cells", "T cells", "NK cells", "B cells",
  "Myeloid cells", "Fibroblasts", "Endothelial cells"
)
hla_main <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
hla_extended <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2")
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRD1")

# ---- Figure 1: all-cell STK31 localisation without a full Seurat object ----
metadata_path <- file.path(project, "02_processed_data", "All_single_cells", "GBC_Metadata.txt")
full_barcodes_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "barcodes.tsv.gz")
full_features_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "features.tsv.gz")
full_matrix_path <- file.path(project, "02_processed_data", "All_single_cells", "counts", "10X_counts", "matrix.mtx.gz")

message("Reading author metadata for the all-cell STK31 location panel")
metadata <- fread(metadata_path, showProgress = FALSE)
metadata <- metadata[name != "type"]
stop_if_not(nrow(metadata) == 1117245L, "Unexpected all-cell metadata row count")
stop_if_not(identical(names(metadata), c("name", "celltype", "subtype", "sample_name")), "Unexpected metadata schema")
full_barcodes <- read_gz_tsv(full_barcodes_path, header = FALSE, col.names = "name")
stop_if_not(identical(metadata$name, full_barcodes$name), "Full metadata and barcode order differ")
full_features <- read_gz_tsv(full_features_path, header = FALSE)
stop_if_not(nrow(full_features) == 32137L && ncol(full_features) >= 2L, "Unexpected full feature table")
stk31_row <- which(full_features[[2]] == "STK31")
stop_if_not(length(stk31_row) == 1L, "STK31 must occur exactly once in full feature symbols")

metadata[, major_celltype := NA_character_]
metadata[subtype == "Malignant epithelial cells", major_celltype := "Malignant epithelial cells"]
metadata[celltype %in% c("CD8+ T cell", "CD4+ T cells"), major_celltype := "T cells"]
metadata[celltype == "NK cells", major_celltype := "NK cells"]
metadata[celltype == "B cells", major_celltype := "B cells"]
metadata[celltype %in% c("Monocytes & Macrophages", "Neutrophils", "Dendritic cells"), major_celltype := "Myeloid cells"]
metadata[celltype == "Mesenchymal cells" & grepl("^F_", subtype), major_celltype := "Fibroblasts"]
metadata[celltype == "Endothelial cells", major_celltype := "Endothelial cells"]
metadata[, group_index := match(major_celltype, major_levels)]
metadata[is.na(group_index), group_index := 0L]
major_cell_counts <- metadata[group_index > 0L, .N, by = .(group_index, major_celltype)]
setorder(major_cell_counts, group_index)
stop_if_not(nrow(major_cell_counts) == length(major_levels), "A requested major cell type has no cells")
fwrite(major_cell_counts, file.path(table_dir, "Figure1_STK31_major_celltype_cell_numbers.csv"))

major_map_path <- file.path(technical_dir, "stage4A_full_to_major_type.int32.bin")
writeBin(as.integer(metadata$group_index), major_map_path, size = 4L, endian = .Platform$endian)
rm(full_barcodes)
gc()

stk_cpp <- file.path(scripts_dir, "stage4A_extract_stk31_major_types.cpp")
stk_bin <- file.path(technical_dir, "stage4A_extract_stk31_major_types")
stk_stats_path <- file.path(table_dir, "Figure1_STK31_major_celltype_stats.csv")
force_stream <- identical(Sys.getenv("STAGE4A_FORCE_STREAM"), "true")
if (force_stream || !file.exists(stk_stats_path)) {
  run_checked(
    "g++",
    c("-O3", "-std=c++17", stk_cpp, "-lz", "-o", stk_bin),
    file.path(log_dir, "stage4A_compile_stk31_stream.log")
  )
  run_checked(
    "/usr/bin/time",
    c("-v", "-o", file.path(log_dir, "stage4A_stk31_major_type_stream.time.txt"), stk_bin,
      full_matrix_path, major_map_path, "32137", "1117245", as.character(stk31_row),
      as.character(length(major_levels)), stk_stats_path),
    file.path(log_dir, "stage4A_stk31_major_type_stream.log")
  )
} else {
  message("Reusing existing, previously completed all-cell STK31 stream output")
}
major_stk <- fread(stk_stats_path)
stop_if_not(nrow(major_stk) == length(major_levels), "Unexpected STK31 major-type statistic rows")
stop_if_not(all(major_stk$declared_nnz == major_stk$actual_nnz), "STK31 stream did not consume declared nnz")
if (!"major_celltype" %in% names(major_stk)) {
  major_stk <- merge(major_cell_counts, major_stk, by = "group_index", sort = FALSE)
  setorder(major_stk, group_index)
  major_stk[, `:=`(
    major_celltype = factor(major_celltype, levels = major_levels),
    STK31_RNA_detected_fraction = safe_fraction(STK31_detected_cell_n, N)
  )]
  major_stk[, STK31_pseudobulk_CPM := STK31_raw_UMI / group_library_UMI * 1e6]
  major_stk[, mean_raw_UMI_per_cell := STK31_raw_UMI / N]
  major_stk[, log1p_STK31_pseudobulk_CPM := log1p(STK31_pseudobulk_CPM)]
  setnames(major_stk, "N", "cell_n")
} else {
  stop_if_not(all(major_stk$cell_n == major_cell_counts$N), "Reused Figure 1 table differs from metadata cell counts")
  major_stk[, major_celltype := factor(major_celltype, levels = major_levels)]
}
fwrite(major_stk, stk_stats_path)

fig1a <- ggplot(major_stk, aes(x = major_celltype, y = "STK31")) +
  geom_point(aes(size = STK31_RNA_detected_fraction, color = log1p_STK31_pseudobulk_CPM)) +
  scale_size_continuous(name = "% RNA detected", range = c(2.5, 13), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "STK31 across author-annotated major cell types") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.ticks.y = element_blank(), axis.text.y = element_blank())

malignant_stk <- major_stk[major_celltype == "Malignant epithelial cells"]
fig1b <- ggplot(malignant_stk, aes(x = "Malignant epithelial cells", y = STK31_RNA_detected_fraction)) +
  geom_col(width = 0.55, fill = "#8c2d04") +
  geom_text(aes(label = paste0(format(STK31_detected_cell_n, big.mark = ","), " / ", format(cell_n, big.mark = ","), " cells\n", fmt_pct(STK31_RNA_detected_fraction, 2))), vjust = -0.35, size = 3.4) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "STK31 RNA detection fraction", title = "STK31 in malignant epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

fig1 <- fig1a + fig1b + plot_layout(widths = c(1.7, 1)) + plot_annotation(tag_levels = "A")
save_figure(fig1, "Figure1_STK31_celltype_localization", width = 11, height = 4.8)

# ---- Figure 2: HLA-I in author-labelled malignant epithelial cells ----
mal_target_path <- file.path(project, "02_processed_data", "malignant_epithelial", "target_counts", "malignant_target_raw_counts.rds")
mal_library_path <- file.path(project, "02_processed_data", "malignant_epithelial", "malignant_cell_library_sizes.tsv")
target_counts <- readRDS(mal_target_path)
mal_library <- fread(mal_library_path)
stop_if_not(inherits(target_counts, "sparseMatrix") && identical(dim(target_counts), c(9L, 96942L)), "Unexpected malignant target sparse matrix")
stop_if_not(nrow(mal_library) == ncol(target_counts) && all(mal_library$all_gene_raw_umi > 0), "Invalid malignant cell library sizes")
stop_if_not(all(hla_main %in% rownames(target_counts)), "A primary HLA-I gene is missing")

mal_hla <- rbindlist(lapply(hla_main, function(gene) {
  values <- as.numeric(target_counts[gene, ])
  data.table(
    gene = gene,
    malignant_cell_n = length(values),
    RNA_detected_n = sum(values > 0),
    RNA_detected_fraction = mean(values > 0),
    total_raw_UMI = sum(values),
    pseudobulk_CPM = sum(values) / sum(mal_library$all_gene_raw_umi) * 1e6,
    mean_raw_UMI_per_cell = mean(values)
  )
}))
mal_hla[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
mal_hla[, gene := factor(gene, levels = rev(hla_main))]
fwrite(mal_hla, file.path(table_dir, "Figure2_malignant_HLAI_expression_stats.csv"))

fig2 <- ggplot(mal_hla, aes(x = "Malignant epithelial cells", y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(name = "% RNA detected", range = c(5, 16), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL, title = "HLA-I-related genes in malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(fig2, "Figure2_malignant_epithelial_HLAI", width = 6.4, height = 4.8)

# ---- Figure 3: author NK/NKT subtype inhibitory receptor spectrum ----
nk_path <- file.path(project, "01_raw_data", "NK.RDS")
message("Reading NK.RDS for the author-defined subtype receptor dot plot")
nk <- readRDS(nk_path)
stop_if_not(inherits(nk, "Seurat") && DefaultAssay(nk) == "RNA", "NK.RDS must be an RNA Seurat object")
nk_counts <- LayerData(nk[["RNA"]], layer = "counts")
nk_meta <- as.data.table(nk[[]], keep.rownames = "barcode")
stop_if_not(identical(colnames(nk_counts), nk_meta$barcode), "NK counts and metadata differ in order")
stop_if_not(all(nk_receptors %in% rownames(nk_counts)), "An NK receptor feature is missing")
stop_if_not(!"KIR2DL2" %in% rownames(nk_counts), "KIR2DL2 unexpectedly present")
nk_meta[, nk_subtype := as.character(celltype)]
nk_meta[, all_gene_raw_umi := as.numeric(Matrix::colSums(nk_counts))]
stop_if_not(all(nk_meta$all_gene_raw_umi > 0), "NK contains zero-library cells")
subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^.*_C([0-9]+)_.*$", "\\1", nk_subtype))]
setorder(subtype_sizes, cluster_number, nk_subtype)
stop_if_not(nrow(subtype_sizes) == 14L, "Expected 14 author NK/NKT subtypes")

nk_dot <- rbindlist(lapply(nk_receptors, function(gene) {
  values <- as.numeric(nk_counts[gene, ])
  working <- nk_meta[, .(nk_subtype, all_gene_raw_umi)]
  working[, raw_UMI := values]
  working[, .(
    NK_cell_n = .N,
    RNA_detected_n = sum(raw_UMI > 0),
    RNA_detected_fraction = mean(raw_UMI > 0),
    total_raw_UMI = sum(raw_UMI),
    mean_logCP10K = mean(log1p(raw_UMI / all_gene_raw_umi * 1e4))
  ), by = nk_subtype][, gene := gene]
}))
nk_dot <- merge(nk_dot, subtype_sizes[, .(nk_subtype, NK_cell_n, cluster_number)], by = c("nk_subtype", "NK_cell_n"), all.x = TRUE, sort = FALSE)
nk_dot[, `:=`(
  gene = factor(gene, levels = rev(nk_receptors)),
  subtype_display = paste0(nk_subtype, "\n(n=", format(NK_cell_n, big.mark = ","), ")")
)]
subtype_display_levels <- nk_dot[gene == "KLRC1"][order(cluster_number, nk_subtype), subtype_display]
nk_dot[, subtype_display := factor(subtype_display, levels = subtype_display_levels)]
setorder(nk_dot, cluster_number, nk_subtype, gene)
fwrite(nk_dot, file.path(table_dir, "Figure3_NK_inhibitory_receptor_subtype_stats.csv"))

fig3 <- ggplot(nk_dot, aes(x = subtype_display, y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = mean_logCP10K)) +
  scale_size_continuous(name = "% RNA detected", range = c(1.4, 10), labels = function(x) paste0(round(100 * x, 1), "%")) +
  scale_color_gradient(name = "mean log1p(CP10K)", low = "#e0f3db", high = "#08589e") +
  labs(x = "Author NK/NKT subtype (cell number shown)", y = NULL, title = "NK inhibitory receptor RNA detection across author-defined subtypes") +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1), axis.ticks = element_blank())
save_figure(fig3, "Figure3_NK_inhibitory_receptor_subtypes", width = 13.5, height = 5.8)

kir3dl1_subtype <- nk_dot[gene == "KIR3DL1"]
figS2 <- ggplot(kir3dl1_subtype[order(RNA_detected_fraction)], aes(x = RNA_detected_fraction, y = reorder(nk_subtype, RNA_detected_fraction))) +
  geom_col(fill = "#6a51a3") +
  geom_text(aes(label = paste0(RNA_detected_n, " / ", NK_cell_n)), hjust = -0.15, size = 3) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x, 1), "%"), expand = expansion(mult = c(0, 0.2))) +
  labs(x = "KIR3DL1 RNA detection fraction", y = NULL, title = "Supplementary: KIR3DL1 detection by author NK/NKT subtype") +
  theme_classic(base_size = 11)
save_figure(figS2, "Supplementary_FigureS2_KIR3DL1_subtype_detection", width = 7.8, height = 5.8)
rm(nk, nk_counts, nk_meta)
gc()

# ---- Figure 4: malignant sample pseudobulk with full-transcriptome TMM factors ----
mal_counts_path <- file.path(project, "02_processed_data", "malignant_epithelial", "10X_counts", "matrix.mtx.gz")
mal_features_path <- file.path(project, "02_processed_data", "malignant_epithelial", "10X_counts", "features.tsv.gz")
mal_cell_metadata_path <- file.path(project, "02_processed_data", "malignant_epithelial", "malignant_cell_metadata.csv.gz")
mal_sample_long_path <- file.path(project, "06_feasibility", "stage3", "malignant_sample_gene_detection.csv")
mal_cells <- fread(mal_cell_metadata_path, showProgress = FALSE)
sample_long <- fread(mal_sample_long_path, showProgress = FALSE)
sample_meta <- unique(sample_long[, .(
  sample_name, patient_id, patient_id_source, site, standardized_site,
  TNM_stage, histological_type, metastatic_group
)])
setorder(sample_meta, sample_name)
stop_if_not(nrow(sample_meta) == 58L && nrow(mal_cells) == 96942L, "Unexpected malignant sample/cell counts")
stop_if_not(!anyDuplicated(sample_meta$sample_name), "Duplicate malignant sample metadata")
cell_to_sample <- match(mal_cells$sample_name, sample_meta$sample_name)
stop_if_not(!anyNA(cell_to_sample), "A malignant cell has no sample mapping")
sample_map_path <- file.path(technical_dir, "stage4A_malignant_cell_to_sample.int32.bin")
writeBin(as.integer(cell_to_sample), sample_map_path, size = 4L, endian = .Platform$endian)
fwrite(sample_meta, file.path(table_dir, "Figure4_malignant_pseudobulk_sample_metadata.csv"))

pb_cpp <- file.path(scripts_dir, "stage4A_malignant_pseudobulk.cpp")
pb_bin <- file.path(technical_dir, "stage4A_malignant_pseudobulk")
pb_mtx_path <- file.path(technical_dir, "stage4A_malignant_pseudobulk_counts.mtx.gz")
pb_audit_path <- file.path(technical_dir, "stage4A_malignant_pseudobulk_stream_audit.tsv")
if (force_stream || !file.exists(pb_mtx_path) || !file.exists(pb_audit_path)) {
  run_checked(
    "g++",
    c("-O3", "-std=c++17", pb_cpp, "-lz", "-o", pb_bin),
    file.path(log_dir, "stage4A_compile_pseudobulk_stream.log")
  )
  run_checked(
    "/usr/bin/time",
    c("-v", "-o", file.path(log_dir, "stage4A_malignant_pseudobulk_stream.time.txt"), pb_bin,
      mal_counts_path, sample_map_path, "32137", "96942", as.character(nrow(sample_meta)), pb_mtx_path, pb_audit_path),
    file.path(log_dir, "stage4A_malignant_pseudobulk_stream.log")
  )
} else {
  stop_if_not(system2("gzip", c("-t", pb_mtx_path)) == 0L, "Existing pseudobulk matrix fails gzip integrity test")
  message("Reusing existing, previously completed malignant pseudobulk stream output")
}
pb_audit <- fread(pb_audit_path)
stop_if_not(pb_audit$declared_nnz == pb_audit$actual_nnz && pb_audit$output_rows == 32137L && pb_audit$output_cols == 58L,
            "Malformed malignant pseudobulk stream output")
pb <- readMM(gzfile(pb_mtx_path))
pb <- as(pb, "CsparseMatrix")
stop_if_not(identical(dim(pb), c(32137L, 58L)), "Unexpected pseudobulk matrix dimensions")
mal_features <- read_gz_tsv(mal_features_path, header = FALSE)
stop_if_not(nrow(mal_features) == nrow(pb), "Malignant feature count differs from pseudobulk matrix")
feature_symbol <- as.character(mal_features[[2]])
target_rows <- vapply(c("STK31", hla_main), function(gene) {
  index <- which(feature_symbol == gene)
  stop_if_not(length(index) == 1L, paste("Expected exactly one pseudobulk feature for", gene))
  index
}, integer(1))
pb_dense <- as.matrix(pb)
rownames(pb_dense) <- make.unique(feature_symbol)
colnames(pb_dense) <- sample_meta$sample_name
keep_gene <- rowSums(pb_dense) > 0
dge <- DGEList(counts = pb_dense[keep_gene, , drop = FALSE])
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
tmm_cpm <- cpm(dge, log = FALSE, prior.count = 0)
target_unique_rows <- match(make.unique(feature_symbol[target_rows]), rownames(log_cpm))
names(target_unique_rows) <- names(target_rows)
stop_if_not(!anyNA(target_unique_rows), "Target pseudobulk row mapping failed")

pb_target <- copy(sample_meta)
for (gene in names(target_rows)) {
  row_index <- target_unique_rows[[gene]]
  pb_target[, (paste0(gene, "_raw_count")) := as.numeric(pb_dense[rownames(log_cpm)[row_index], ])]
  pb_target[, (paste0(gene, "_TMM_CPM")) := as.numeric(tmm_cpm[row_index, ])]
  pb_target[, (paste0(gene, "_log2CPM")) := as.numeric(log_cpm[row_index, ])]
}
fwrite(pb_target, file.path(table_dir, "Figure4_malignant_TMM_pseudobulk_target_genes.csv"))
fwrite(as.data.table(dge$samples, keep.rownames = "sample_name"), file.path(table_dir, "Figure4_TMM_normalization_factors.csv"))

correlation_table <- rbindlist(lapply(c("Primary tumor", "All malignant-containing samples"), function(scope) {
  scope_data <- if (scope == "Primary tumor") pb_target[standardized_site == "Primary tumor"] else copy(pb_target)
  rbindlist(lapply(hla_main, function(gene) {
    x <- scope_data$STK31_log2CPM
    y <- scope_data[[paste0(gene, "_log2CPM")]]
    test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    data.table(
      scope = scope,
      gene = gene,
      n = length(x),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      STK31_nonzero_sample_n = sum(scope_data$STK31_raw_count > 0),
      gene_nonzero_sample_n = sum(scope_data[[paste0(gene, "_raw_count")]] > 0),
      analysis_label = if (scope == "Primary tumor") "exploratory_primary_patient_samples" else "supplementary_nonindependent_sample_level"
    )
  }))
}))
fwrite(correlation_table, file.path(table_dir, "Figure4_STK31_HLAI_Spearman_correlations.csv"))

make_trend_plot <- function(scope_name, title_text) {
  scope_data <- if (scope_name == "Primary tumor") pb_target[standardized_site == "Primary tumor"] else copy(pb_target)
  long <- rbindlist(lapply(hla_main, function(gene) data.table(
    sample_name = scope_data$sample_name,
    gene = gene,
    STK31_log2CPM = scope_data$STK31_log2CPM,
    HLA_log2CPM = scope_data[[paste0(gene, "_log2CPM")]]
  )))
  annotation <- correlation_table[scope == scope_name, .(gene, n, spearman_rho, p_value)]
  ranges <- long[, .(x = min(STK31_log2CPM) + 0.03 * diff(range(STK31_log2CPM)),
                     y = max(HLA_log2CPM) - 0.04 * diff(range(HLA_log2CPM))), by = gene]
  annotation <- merge(annotation, ranges, by = "gene", all.x = TRUE)
  annotation[, label := sprintf("n=%d\nSpearman rho=%.2f\np=%.3g", n, spearman_rho, p_value)]
  ggplot(long, aes(x = STK31_log2CPM, y = HLA_log2CPM)) +
    geom_point(size = 2.5, alpha = 0.78, color = "#1f78b4") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, color = "#b2182b", linetype = "dashed") +
    geom_text(data = annotation, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.1) +
    facet_wrap(~gene, scales = "free_y", nrow = 1) +
    labs(x = "STK31 log2(TMM-CPM + prior count)", y = "HLA-I gene log2(TMM-CPM + prior count)", title = title_text) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", color = "grey70"))
}

fig4 <- make_trend_plot("Primary tumor", "Primary-tumor malignant pseudobulk: exploratory STK31-HLA-I trends")
save_figure(fig4, "Figure4_primary_STK31_HLAI_pseudobulk_trends", width = 12, height = 3.7)
figS1 <- make_trend_plot("All malignant-containing samples", "Supplementary: all malignant-containing samples (sample-level, non-independent)")
save_figure(figS1, "Supplementary_FigureS1_all_malignant_STK31_HLAI_pseudobulk_trends", width = 12, height = 3.7)

# ---- Required report: figure hierarchy and interpretation limits ----
primary_hlab <- correlation_table[scope == "Primary tumor" & gene == "HLA-B"]
primary_n <- primary_hlab$n
trend_judgement <- if (abs(primary_hlab$spearman_rho) >= 0.5) {
  "可作为明确标注为探索性的样本级趋势图展示；n较小，不能表述为因果或稳健患者级证据。"
} else {
  "趋势幅度有限；建议保留为补充或用作简短描述，而不作为正文的关键关联证据。"
}
summary_lines <- c(
  "# Stage 4A：当前湿实验文章的精简公共单细胞验证",
  "",
  "## 范围",
  "",
  "本任务生成3张主图和3张直接补充图；未进行全基因差异分析、通路富集、CellChat/NicheNet、重新聚类、轨迹、inferCNV、生存或发现型网络分析。共享Stage 1–3输入保持只读，未写入`14_discovery/`。",
  "",
  "## 推荐图组",
  "",
  "| 图 | 定位 | 准确figure message |",
  "|---|---|---|",
  "| Figure 1 | 主图 | 作者注释的多种主要细胞类型均可检测STK31 RNA；在所比较的细胞类型中，恶性上皮细胞的STK31检出比例最高，但单细胞信号仍稀疏。 |",
  "| Figure 2 | 主图 | 恶性上皮细胞中HLA-A、HLA-B、HLA-C和B2M具有广泛RNA检出，为HLA-I相关湿实验提供临床组织背景。 |",
  "| Figure 3 | 主图 | 作者定义的NK/NKT亚群中可检测到KLRC1、KIR3DL1、KIR2DL1、KIR2DL3和KLRD1；KIR3DL1为低检出候选受体，不能解释为生物学缺失。 |",
  "| Figure 4 | 补充（探索性） | 在11个原发肿瘤恶性细胞pseudobulk样本中，未观察到清晰的STK31–HLA-B样本级趋势；该图用于透明报告，不构成因果或稳健关联证据。 |",
  "| Supplementary Figure S1 | 补充 | 所有含恶性细胞样本的STK31–HLA-I样本级趋势供透明性检查；该集合含多部位和重复患者样本，不能代替原发患者分析。 |",
  "| Supplementary Figure S2 | 补充 | KIR3DL1在作者NK/NKT亚群中呈稀疏、异质的RNA检出；图中同时给出阳性数和细胞数。 |",
  "",
  "## STK31–HLA-B主要趋势",
  "",
  sprintf("原发肿瘤分析使用%d个样本（每个为一个派生患者的Primary sample），全转录组恶性pseudobulk经edgeR TMM标准化。STK31对HLA-B的Spearman rho=%.3f，p=%.4g。%s", primary_n, primary_hlab$spearman_rho, primary_hlab$p_value, trend_judgement),
  "",
  "所有相关性在不筛除STK31低/零样本的条件下计算；p值仅作为探索性描述。",
  "",
  "## KIR3DL1是否保留在主图",
  "",
  "保留在Figure 3的受体DotPlot中，但不单列为机制主图或以其低检出作生物学结论。它服务于湿实验提出的HLA-Bw4–KIR3DL1候选轴；总HLA-B RNA不能替代HLA-Bw4分型，KIR3DL1的scRNA稀疏还可能受掉零、同源性和多态性影响。",
  "",
  "## 不应写入正文结论的内容",
  "",
  "- 不应把STK31–HLA-B的样本级趋势写成STK31调控HLA-B的因果证据。",
  "- 不应把KIR3DL1低检出写成NK缺失KIR3DL1。",
  "- 不应把HLA-B总RNA写成HLA-Bw4或直接受体–配体互作证据。",
  "- 原发精确配对中结构可评估样本只有4个，因此不应在公共数据中声称患者级STK31–HLA-I–KIR3DL1机制已被验证。",
  "",
  "## 输出",
  "",
  "- 图：`11_figures/current_paper_minimal/`",
  "- 图表数据：`12_tables/current_paper_minimal/`",
  "- 流式技术文件：`12_tables/current_paper_minimal/technical/`",
  "- 本报告由`04_scripts/stage4A_current_paper_minimal.R`生成。",
  "",
  "Stage 4A到此停止；未进入深度挖掘。"
)
writeLines(summary_lines, file.path(report_dir, "stage4A_current_paper_minimal_summary.md"), useBytes = TRUE)
fwrite(data.table(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Supplementary Figure S1", "Supplementary Figure S2"),
  role = c("main", "main", "main", "supplementary_exploratory", "supplementary", "supplementary"),
  file_stem = c(
    "Figure1_STK31_celltype_localization", "Figure2_malignant_epithelial_HLAI",
    "Figure3_NK_inhibitory_receptor_subtypes", "Figure4_primary_STK31_HLAI_pseudobulk_trends",
    "Supplementary_FigureS1_all_malignant_STK31_HLAI_pseudobulk_trends",
    "Supplementary_FigureS2_KIR3DL1_subtype_detection"
  )
), file.path(table_dir, "figure_manifest_current_paper_minimal.csv"))

message("STAGE4A_CURRENT_PAPER_MINIMAL_OK")
