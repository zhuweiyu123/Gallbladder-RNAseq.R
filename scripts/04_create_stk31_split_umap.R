suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

base_dir <- "/home/zhuweiyu/codex-r"
target_gene <- "STK31"
tumor_epithelial_celltype_label <- "Epithelial"
high_label <- "STK31_high_tumor_epithelial"
low_label <- "STK31_low_tumor_epithelial"

read_summary_value <- function(path, item) {
  summary <- read.csv(path, stringsAsFactors = FALSE)
  value <- summary$value[summary$item == item]
  if (length(value) == 0) stop("Missing summary item: ", item, " in ", path)
  value[[1]]
}

get_gene_expr <- function(object, gene) {
  assay <- DefaultAssay(object)
  mat <- tryCatch(
    GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) GetAssayData(object, assay = assay, slot = "data")
  )
  if (!gene %in% rownames(mat)) stop(gene, " was not found in assay ", assay)
  as.numeric(mat[gene, ])
}

add_split_celltype <- function(object, cutoff) {
  if (!"manual_celltype" %in% colnames(object@meta.data)) {
    stop("manual_celltype is missing from object metadata")
  }
  stk31_expr <- get_gene_expr(object, target_gene)
  tumor_epithelial_idx <- as.character(object$manual_celltype) == tumor_epithelial_celltype_label
  if (sum(tumor_epithelial_idx) == 0) {
    stop("No tumor epithelial cells found with manual_celltype == ", tumor_epithelial_celltype_label)
  }

  group <- rep("Other", length(stk31_expr))
  group[tumor_epithelial_idx & stk31_expr > cutoff] <- high_label
  group[tumor_epithelial_idx & stk31_expr <= cutoff] <- low_label
  names(group) <- colnames(object)

  object$tumor_epithelial_stk31_group <- group
  object$manual_celltype_stk31_tumor_epithelial_split <- as.character(object$manual_celltype)
  object$manual_celltype_stk31_tumor_epithelial_split[group == high_label] <- high_label
  object$manual_celltype_stk31_tumor_epithelial_split[group == low_label] <- low_label
  object$manual_celltype_stk31_tumor_epithelial_split <- factor(object$manual_celltype_stk31_tumor_epithelial_split)
  object
}

write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

plot_split_umap <- function(object, path, title) {
  write_plot(
    path,
    DimPlot(
      object,
      reduction = "umap",
      group.by = "manual_celltype_stk31_tumor_epithelial_split",
      label = TRUE,
      repel = TRUE
    ) + ggtitle(title)
  )
  message("Wrote ", path)
}

message("Creating merged manual celltype UMAP with STK31 tumor epithelial split")
merged_out_dir <- file.path(base_dir, "results", "merged_stk31_nk_analysis")
merged_obj <- readRDS(file.path(base_dir, "results", "merged_basic_seurat", "gallbladder_cancer_merged_basic_seurat.rds"))
merged_annotation <- read.csv(file.path(merged_out_dir, "cluster_celltype_annotation.csv"), stringsAsFactors = FALSE)
cluster_to_manual <- setNames(merged_annotation$manual_celltype, as.character(merged_annotation$cluster))
merged_obj$manual_celltype <- unname(cluster_to_manual[as.character(merged_obj$seurat_clusters)])
merged_obj$manual_celltype[is.na(merged_obj$manual_celltype) | !nzchar(merged_obj$manual_celltype)] <- "Unassigned"
merged_obj$manual_celltype <- factor(merged_obj$manual_celltype)
merged_cutoff <- as.numeric(read_summary_value(file.path(merged_out_dir, "analysis_summary.csv"), "tumor_epithelial_stk31_high_cutoff"))
merged_obj <- add_split_celltype(merged_obj, merged_cutoff)
plot_split_umap(
  merged_obj,
  file.path(merged_out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  "Manual annotation with STK31 high/low tumor epithelial split (merged)"
)

message("Creating tissue2 manual celltype UMAP with STK31 tumor epithelial split")
tissue2_out_dir <- file.path(base_dir, "results", "tissue2_stk31_nk_analysis")
tissue2_obj <- readRDS(file.path(tissue2_out_dir, "tissue2_stk31_nk_annotated_seurat.rds"))
tissue2_cutoff <- as.numeric(read_summary_value(file.path(tissue2_out_dir, "analysis_summary.csv"), "tumor_epithelial_stk31_high_cutoff"))
tissue2_obj <- add_split_celltype(tissue2_obj, tissue2_cutoff)
plot_split_umap(
  tissue2_obj,
  file.path(tissue2_out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  "Manual annotation with STK31 high/low tumor epithelial split"
)

validation_path <- file.path(tissue2_out_dir, "validation", "tissue2_stk31_nk_validated_seurat.rds")
if (file.exists(validation_path)) {
  message("Creating validation UMAP with STK31 tumor epithelial split")
  validation_obj <- readRDS(validation_path)
  validation_obj <- add_split_celltype(validation_obj, tissue2_cutoff)
  validation_obj$validation_celltype_stk31_tumor_epithelial_split <- validation_obj$manual_celltype_stk31_tumor_epithelial_split
  write_plot(
    file.path(tissue2_out_dir, "validation", "umap_validation_celltypes_with_stk31_tumor_epithelial_split.pdf"),
    DimPlot(
      validation_obj,
      reduction = "umap",
      group.by = "validation_celltype_stk31_tumor_epithelial_split",
      label = TRUE,
      repel = TRUE
    ) + ggtitle("Validation cell types with STK31 high/low tumor epithelial split")
  )
}

message("Done")
