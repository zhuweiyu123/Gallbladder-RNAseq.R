#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

project <- "/home/zhuweiyu/codex-r"
local_root <- file.path(project, "results", "merged_cnv_binary_annotation")
object_path <- file.path(local_root, "merged_cnv_binary_annotated.rds")
local_figure_dir <- file.path(local_root, "figures")
local_table_dir <- file.path(local_root, "tables")

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
save_figure <- function(plot, directory, stem, width, height) {
  ggsave(file.path(directory, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# P1/P2/P4 epithelial figures from the final integrated object.
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj, assay = "RNA")
}
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "RNA counts must be sparse")
target_genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
stop_if_not(all(target_genes %in% rownames(counts)), "A requested epithelial target gene is absent")

patients <- c("P1", "P2", "P4")
meta_all <- as.data.table(obj[[]], keep.rownames = "cell_id")
meta <- meta_all[
  patient_id %in% patients &
    epithelial_binary_label %in% c("Normal epithelial cells", "Malignant epithelial cells")
]
meta[, display_group := paste(patient_id, sub(" epithelial cells$", "", epithelial_binary_label))]
group_levels <- as.vector(rbind(
  paste(patients, "Normal"), paste(patients, "Malignant")
))
meta[, display_group := factor(display_group, levels = group_levels)]
expected_n <- c(
  "P1 Normal" = 509L, "P1 Malignant" = 251L,
  "P2 Normal" = 1485L, "P2 Malignant" = 697L,
  "P4 Normal" = 47L, "P4 Malignant" = 99L
)
stop_if_not(
  identical(as.integer(table(meta$display_group)), as.integer(expected_n[group_levels])),
  "Unexpected P1/P2/P4 epithelial group sizes"
)

group_stats <- function(gene) {
  rbindlist(lapply(group_levels, function(group) {
    cells <- meta[display_group == group, cell_id]
    values <- as.numeric(counts[gene, cells, drop = TRUE])
    library_umi <- sum(counts[, cells, drop = FALSE])
    data.table(
      gene = gene,
      display_group = group,
      patient_id = sub(" .*", "", group),
      epithelial_group = sub("^[^ ]+ ", "", group),
      cell_n = length(cells),
      RNA_detected_n = sum(values > 0),
      RNA_detected_fraction = mean(values > 0),
      total_raw_UMI = sum(values),
      group_library_UMI = library_umi,
      pseudobulk_CPM = sum(values) / library_umi * 1e6,
      mean_raw_UMI_per_cell = mean(values)
    )
  }))
}

stk <- group_stats("STK31")
stk[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  display_group = factor(display_group, levels = group_levels)
)]
fwrite(stk, file.path(local_table_dir, "Figure_P1_P2_P4_STK31_epithelial_stats.csv"))

p_stk_a <- ggplot(stk, aes(x = display_group, y = "STK31")) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(3, 13),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "STK31 in P1/P2/P4 epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

stk_malignant <- stk[epithelial_group == "Malignant"]
stk_malignant[, patient_id := factor(patient_id, levels = patients)]
p_stk_b <- ggplot(stk_malignant, aes(x = patient_id, y = RNA_detected_fraction)) +
  geom_col(width = 0.58, fill = "#8c2d04") +
  geom_text(
    aes(label = paste0(RNA_detected_n, " / ", cell_n, " cells\n",
                       fmt_pct(RNA_detected_fraction, 2))),
    vjust = -0.35, size = 3.2
  ) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "STK31 RNA detection fraction",
       title = "STK31 in malignant epithelial cells") +
  theme_classic(base_size = 11)

p_stk <- p_stk_a + p_stk_b + plot_layout(widths = c(1.75, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(
  p_stk, local_figure_dir, "Figure_P1_P2_P4_STK31_epithelial_localization",
  12.5, 5.0
)

hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
hla <- rbindlist(lapply(hla_genes, group_stats))[epithelial_group == "Malignant"]
hla[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  gene = factor(gene, levels = rev(hla_genes)),
  display_group = factor(display_group, levels = paste(patients, "Malignant"))
)]
fwrite(hla, file.path(local_table_dir, "Figure_P1_P2_P4_malignant_HLAI_plus_HLAE_stats.csv"))

p_hla <- ggplot(hla, aes(x = display_group, y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL,
       title = "HLA-I-related genes in P1/P2/P4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(
  p_hla, local_figure_dir, "Figure_P1_P2_P4_malignant_epithelial_HLAI_plus_HLAE",
  8.4, 5.4
)

# Local P1/P2/P4 NK inhibitory-receptor figure using the frozen local NK annotation.
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRD1")
missing_nk_receptors <- setdiff(nk_receptors, rownames(counts))
stop_if_not(identical(missing_nk_receptors, "KIR2DL2"),
            paste("Unexpected missing local NK receptors:",
                  paste(missing_nk_receptors, collapse = ", ")))
available_nk_receptors <- setdiff(nk_receptors, missing_nk_receptors)
nk_meta <- meta_all[
  patient_id %in% patients &
    analysis_celltype_previous_or_new == "NK_cell" &
    as.character(previous_tnk_subcluster) %in% c("4", "7", "12")
]
stop_if_not(nrow(nk_meta) == 1129L, "Unexpected strict P1/P2/P4 local NK total")
nk_counts <- counts[, nk_meta$cell_id, drop = FALSE]
nk_meta[, `:=`(
  nk_subtype = paste0("NK_C", as.character(previous_tnk_subcluster)),
  all_gene_raw_umi = as.numeric(Matrix::colSums(nk_counts))
)]
stop_if_not(all(nk_meta$all_gene_raw_umi > 0), "A local NK library size is invalid")

subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^NK_C", "", nk_subtype))]
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
  file.path(local_table_dir, "Figure_P1_P2_P4_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
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
    x = "Strict NK subcluster (P1/P2/P4 cell number shown)", y = NULL,
    title = "NK inhibitory receptors in local P1/P2/P4 strict NK cells",
    caption = "KIR2DL2 is absent from the local count feature matrix; × denotes unavailable, not zero."
  ) +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        axis.ticks = element_blank())
fig3_width <- max(10, 0.85 * uniqueN(nk_dot$nk_subtype))
save_figure(
  p_nk, local_figure_dir, "Figure_P1_P2_P4_NK_inhibitory_receptor_subtypes_plus_KIR2DL2",
  fig3_width, 6.5
)

required <- c(
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_STK31_epithelial_localization.", c("png", "pdf")
  )),
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_malignant_epithelial_HLAI_plus_HLAE.", c("png", "pdf")
  )),
  file.path(local_table_dir, "Figure_P1_P2_P4_STK31_epithelial_stats.csv"),
  file.path(local_table_dir, "Figure_P1_P2_P4_malignant_HLAI_plus_HLAE_stats.csv"),
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_NK_inhibitory_receptor_subtypes_plus_KIR2DL2.", c("png", "pdf")
  )),
  file.path(local_table_dir, "Figure_P1_P2_P4_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "An updated figure or table is absent or empty")
message("LOCAL_P1_P2_P4_STRICT_NK_HLAE_AND_KIR2DL2_FIGURES_PASS")
print(stk)
print(hla)
print(nk_dot[gene == "KIR2DL2"])
