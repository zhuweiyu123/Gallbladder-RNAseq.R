#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
figure_dir <- file.path(project, "11_figures", "current_paper_minimal", "paired_QC4")
table_dir <- file.path(project, "12_tables", "current_paper_minimal", "paired_QC4")
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

qc4 <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")

# QC4 malignant epithelial HLA-I figure with HLA-E.
hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
target_counts <- readRDS(file.path(
  project, "02_processed_data", "malignant_epithelial", "target_counts",
  "malignant_target_raw_counts.rds"
))
mal_meta <- fread(file.path(
  project, "02_processed_data", "malignant_epithelial",
  "malignant_cell_metadata.csv.gz"
), showProgress = FALSE)
mal_library <- fread(file.path(
  project, "02_processed_data", "malignant_epithelial",
  "malignant_cell_library_sizes.tsv"
), showProgress = FALSE)
stop_if_not(identical(colnames(target_counts), mal_meta$barcode),
            "Malignant target matrix and metadata order differ")
stop_if_not(
  identical(as.integer(mal_library$malignant_col), as.integer(mal_meta$malignant_col)),
  "Malignant library rows differ from metadata"
)
stop_if_not(all(hla_genes %in% rownames(target_counts)),
            "A requested QC4 HLA-I target is absent")
mal_idx <- which(mal_meta$sample_name %in% qc4)
stop_if_not(length(mal_idx) == 13368L,
            "QC4 malignant epithelial cell total must be 13,368")
qc_target <- target_counts[, mal_idx, drop = FALSE]
qc_library <- mal_library$all_gene_raw_umi[mal_idx]
hla_stats <- rbindlist(lapply(hla_genes, function(gene) {
  values <- as.numeric(qc_target[gene, ])
  data.table(
    gene = gene,
    malignant_cell_n = length(values),
    RNA_detected_n = sum(values > 0),
    RNA_detected_fraction = mean(values > 0),
    total_raw_UMI = sum(values),
    pseudobulk_CPM = sum(values) / sum(qc_library) * 1e6
  )
}))
hla_stats[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  gene = factor(gene, levels = rev(hla_genes))
)]
fwrite(
  hla_stats,
  file.path(table_dir, "QC4_Figure2_malignant_HLAI_stats_plus_HLAE.csv")
)

p_hla <- ggplot(hla_stats, aes("Malignant epithelial cells", gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL,
       title = "HLA-I-related genes in paired-QC4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(
  p_hla, "QC4_Figure2_malignant_epithelial_HLAI_plus_HLAE",
  6.8, 5.4
)

rm(target_counts, qc_target)
gc()

# QC4 public NK/NKT subtype figure with KIR2DL2.
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRD1")
nk <- readRDS(file.path(project, "01_raw_data", "NK.RDS"))
stop_if_not(inherits(nk, "Seurat") && DefaultAssay(nk) == "RNA",
            "NK.RDS is not the expected RNA Seurat object")
nk_counts_all <- LayerData(nk[["RNA"]], layer = "counts")
nk_meta_all <- as.data.table(nk[[]], keep.rownames = "barcode")
stop_if_not(identical(colnames(nk_counts_all), nk_meta_all$barcode),
            "NK counts and metadata order differ")
nk_idx <- which(nk_meta_all$orig.ident %in% qc4)
stop_if_not(length(nk_idx) == 1135L, "QC4 NK cell total must be 1,135")
nk_counts <- nk_counts_all[, nk_idx, drop = FALSE]
nk_meta <- nk_meta_all[nk_idx]
nk_meta[, `:=`(
  nk_subtype = as.character(celltype),
  all_gene_raw_umi = as.numeric(Matrix::colSums(nk_counts))
)]
missing_nk_receptors <- setdiff(nk_receptors, rownames(nk_counts))
stop_if_not(identical(missing_nk_receptors, "KIR2DL2"),
            paste("Unexpected missing QC4 NK receptors:",
                  paste(missing_nk_receptors, collapse = ", ")))
available_nk_receptors <- setdiff(nk_receptors, missing_nk_receptors)
stop_if_not(all(nk_meta$all_gene_raw_umi > 0), "QC4 NK library sizes are invalid")

subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^.*_C([0-9]+)_.*$", "\\1", nk_subtype))]
setorder(subtype_sizes, cluster_number, nk_subtype)
nk_dot <- rbindlist(lapply(available_nk_receptors, function(gene) {
  working <- nk_meta[, .(nk_subtype, all_gene_raw_umi)]
  working[, raw_UMI := as.numeric(nk_counts[gene, ])]
  working[, .(
    NK_cell_n = .N,
    RNA_detected_n = sum(raw_UMI > 0),
    RNA_detected_fraction = mean(raw_UMI > 0),
    total_raw_UMI = sum(raw_UMI),
    mean_logCP10K = mean(log1p(raw_UMI / all_gene_raw_umi * 1e4))
  ), by = nk_subtype][, gene := gene]
}))
nk_dot <- merge(
  nk_dot, subtype_sizes[, .(nk_subtype, NK_cell_n, cluster_number)],
  by = c("nk_subtype", "NK_cell_n"), all.x = TRUE, sort = FALSE
)
nk_unavailable <- copy(subtype_sizes)
nk_unavailable[, `:=`(
  RNA_detected_n = NA_integer_,
  RNA_detected_fraction = NA_real_,
  total_raw_UMI = NA_real_,
  mean_logCP10K = NA_real_,
  gene = "KIR2DL2"
)]
nk_dot <- rbindlist(list(nk_dot, nk_unavailable), use.names = TRUE, fill = TRUE)
nk_dot[, `:=`(
  gene = factor(gene, levels = rev(nk_receptors)),
  subtype_display = paste0(nk_subtype, "\n(n=", NK_cell_n, ")")
)]
subtype_levels <- nk_dot[gene == "KLRC1"][order(cluster_number, nk_subtype), subtype_display]
nk_dot[, subtype_display := factor(subtype_display, levels = subtype_levels)]
fwrite(
  nk_dot,
  file.path(table_dir, "QC4_Figure3_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)

p_nk <- ggplot(nk_dot, aes(subtype_display, gene)) +
  geom_point(
    data = nk_dot[!is.na(RNA_detected_fraction)],
    aes(size = RNA_detected_fraction, color = mean_logCP10K)
  ) +
  geom_point(
    data = nk_dot[is.na(RNA_detected_fraction)],
    shape = 4, color = "#777777", size = 3, stroke = 1
  ) +
  scale_size_continuous(
    name = "% RNA detected", range = c(1.4, 10),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "mean log1p(CP10K)", low = "#e0f3db", high = "#08589e") +
  labs(
    x = "Author NK/NKT subtype (QC4 cell number shown)", y = NULL,
    title = "NK inhibitory receptors in paired-QC4 primary samples",
    caption = "KIR2DL2 is absent from the QC4 NK count feature matrix; × denotes unavailable, not zero."
  ) +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        axis.ticks = element_blank())
fig3_width <- max(10, 0.85 * uniqueN(nk_dot$nk_subtype))
save_figure(
  p_nk, "QC4_Figure3_NK_inhibitory_receptor_subtypes_plus_KIR2DL2",
  fig3_width, 6.5
)

required <- c(
  file.path(figure_dir, paste0(
    "QC4_Figure2_malignant_epithelial_HLAI_plus_HLAE.", c("png", "pdf")
  )),
  file.path(figure_dir, paste0(
    "QC4_Figure3_NK_inhibitory_receptor_subtypes_plus_KIR2DL2.", c("png", "pdf")
  )),
  file.path(table_dir, "QC4_Figure2_malignant_HLAI_stats_plus_HLAE.csv"),
  file.path(table_dir, "QC4_Figure3_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "An updated public QC4 output is absent or empty")
message("PUBLIC_QC4_HLAE_KIR2DL2_FIGURES_PASS")
print(hla_stats)
print(nk_dot[gene == "KIR2DL2"])
