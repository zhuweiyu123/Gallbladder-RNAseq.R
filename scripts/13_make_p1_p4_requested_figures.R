#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

root <- "/home/zhuweiyu/codex-r/results/merged_cnv_binary_annotation"
object_path <- file.path(root, "merged_cnv_binary_annotated.rds")
figure_dir <- file.path(root, "figures")
table_dir <- file.path(root, "tables")

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
safe_fraction <- function(numerator, denominator) ifelse(denominator > 0, numerator / denominator, NA_real_)
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj, assay = "RNA")
}
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "RNA counts must be sparse")
genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M")
stop_if_not(all(genes %in% rownames(counts)), "A requested gene is absent")

meta <- as.data.table(obj[[]], keep.rownames = "cell_id")
meta <- meta[patient_id %in% c("P1", "P4") &
  epithelial_binary_label %in% c("Normal epithelial cells", "Malignant epithelial cells")]
meta[, display_group := paste(patient_id, sub(" epithelial cells$", "", epithelial_binary_label))]
group_levels <- c("P1 Normal", "P1 Malignant", "P4 Normal", "P4 Malignant")
meta[, display_group := factor(display_group, levels = group_levels)]
expected_n <- c("P1 Normal" = 509L, "P1 Malignant" = 251L,
                "P4 Normal" = 47L, "P4 Malignant" = 99L)
observed_n <- table(meta$display_group)
stop_if_not(identical(as.integer(observed_n), as.integer(expected_n)),
            "Unexpected P1/P4 epithelial group sizes")

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
stk[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
stk[, display_group := factor(display_group, levels = group_levels)]
fwrite(stk, file.path(table_dir, "Figure_P1_P4_STK31_epithelial_stats.csv"))

p_stk_a <- ggplot(stk, aes(x = display_group, y = "STK31")) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(3, 13),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "STK31 in P1/P4 epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

stk_malignant <- stk[epithelial_group == "Malignant"]
p_stk_b <- ggplot(stk_malignant, aes(x = patient_id, y = RNA_detected_fraction)) +
  geom_col(width = 0.58, fill = "#8c2d04") +
  geom_text(
    aes(label = paste0(RNA_detected_n, " / ", cell_n, " cells\n",
                       fmt_pct(RNA_detected_fraction, 2))),
    vjust = -0.35, size = 3.4
  ) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "STK31 RNA detection fraction",
       title = "STK31 in malignant epithelial cells") +
  theme_classic(base_size = 11)

p_stk <- p_stk_a + p_stk_b + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(p_stk, "Figure_P1_P4_STK31_epithelial_localization", 11, 4.8)

hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
hla <- rbindlist(lapply(hla_genes, group_stats))
hla <- hla[epithelial_group == "Malignant"]
hla[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
hla[, gene := factor(gene, levels = rev(hla_genes))]
hla[, display_group := factor(display_group, levels = c("P1 Malignant", "P4 Malignant"))]
fwrite(hla, file.path(table_dir, "Figure_P1_P4_malignant_HLAI_stats.csv"))

p_hla <- ggplot(hla, aes(x = display_group, y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL, title = "HLA-I-related genes in P1/P4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(p_hla, "Figure_P1_P4_malignant_epithelial_HLAI", 7.2, 4.8)

required <- c(
  file.path(figure_dir, paste0("Figure_P1_P4_STK31_epithelial_localization.", c("png", "pdf"))),
  file.path(figure_dir, paste0("Figure_P1_P4_malignant_epithelial_HLAI.", c("png", "pdf"))),
  file.path(table_dir, "Figure_P1_P4_STK31_epithelial_stats.csv"),
  file.path(table_dir, "Figure_P1_P4_malignant_HLAI_stats.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "A requested P1/P4 output is absent or empty")
message("P1_P4_REQUESTED_FIGURES_PASS")
print(stk)
print(hla)
