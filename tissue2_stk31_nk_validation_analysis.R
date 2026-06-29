suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260625)

base_result_dir <- "/home/zhuweiyu/codex-r/results/tissue2_stk31_nk_analysis"
input_rds <- file.path(base_result_dir, "tissue2_stk31_nk_annotated_seurat.rds")
out_dir <- "/home/zhuweiyu/codex-r/results/tissue2_stk31_nk_validation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"

analysis_summary_file <- file.path(base_result_dir, "analysis_summary.csv")
manual_annotation_file <- file.path(base_result_dir, "manual_celltype_annotation.csv")
candidate_lr_file <- file.path(base_result_dir, "candidate_ligand_receptor_pathways.csv")
candidate_pathway_summary_file <- file.path(base_result_dir, "candidate_pathway_summary.csv")
de_stk31_vs_nk_file <- file.path(base_result_dir, "markers_stk31_high_vs_nk.csv")
de_stk31_vs_all_file <- file.path(base_result_dir, "markers_stk31_high_vs_all_other.csv")
de_nk_vs_all_file <- file.path(base_result_dir, "markers_nk_vs_all_other.csv")

identity_marker_sets <- list(
  Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "ACTA2", "TAGLN"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "ESAM"),
  T_cell = c("CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "CCR7"),
  NK_cell = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCN1", "TYROBP", "CST3"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "CD68"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "BANK1", "CD74"),
  Plasma_cell = c("MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT")
)

nk_function_sets <- list(
  NK_cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH"),
  NK_activation = c("IFNG", "TNF", "CD69", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5"),
  NK_checkpoint = c("TIGIT", "HAVCR2", "LAG3", "PDCD1", "CTLA4", "TOX"),
  NK_migration = c("CXCR3", "CCR5", "XCR1", "CX3CR1", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2"),
  TGFbeta_response = c("TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "TAGLN", "ACTA2")
)

pathway_panel_genes <- unique(c(
  "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
  "NECTIN2", "PVR", "TIGIT", "HAVCR2", "LAG3", "PDCD1", "CD274",
  "KLRK1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "KIR2DL1", "KIR3DL1",
  "IL18", "IL18R1", "IL15", "IL2RB", "IFNG", "IFNGR1",
  "TNFSF10", "TNFRSF10B", "FASLG", "FAS", "ICAM1", "ITGAL",
  "CXCL9", "CXCL10", "CXCL11", "CCL5", "XCL1", "XCL2"
))

validation_feature_panel <- unique(c(
  target_gene,
  unlist(identity_marker_sets, use.names = FALSE),
  unlist(nk_function_sets, use.names = FALSE),
  pathway_panel_genes
))

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

parse_semicolon_item <- function(summary_df, item_name, fallback = character(0)) {
  if (is.null(summary_df) || !all(c("item", "value") %in% colnames(summary_df))) {
    return(fallback)
  }
  idx <- which(summary_df$item == item_name)
  if (length(idx) == 0) {
    return(fallback)
  }
  value <- summary_df$value[idx[1]]
  if (is.na(value) || !nzchar(value)) {
    return(fallback)
  }
  unlist(strsplit(value, ";", fixed = TRUE))
}

available_genes <- function(object, genes) {
  intersect(genes, rownames(object))
}

fetch_gene_matrix <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) {
    return(NULL)
  }
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

mean_nonzero <- function(x) {
  mean(x > 0)
}

write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

write_multi_page_pdf <- function(path, plots, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  for (plt in plots) {
    print(plt)
  }
  dev.off()
}

summarize_gene_by_group <- function(object, genes, group_vec, group_name = "group") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) {
    return(data.frame())
  }
  expr <- fetch_gene_matrix(object, genes)
  group_vec <- as.character(group_vec)
  groups <- sort(unique(group_vec))

  out <- do.call(rbind, lapply(groups, function(grp) {
    idx <- group_vec == grp
    data.frame(
      group = grp,
      gene = genes,
      avg_expr = Matrix::rowMeans(expr[, idx, drop = FALSE]),
      pct_expr = Matrix::rowMeans(expr[, idx, drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )
  }))
  names(out)[1] <- group_name
  out
}

top_de_genes <- function(markers, n = 25) {
  if (is.null(markers) || nrow(markers) == 0) {
    return(data.frame())
  }
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC), markers$gene), ]
  head(markers, n)
}

run_optional_go <- function(markers, prefix) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Skipping GO enrichment for ", prefix, " because clusterProfiler/org.Hs.eg.db is not installed.")
    return(invisible(NULL))
  }

  if (!all(c("gene", "avg_log2FC") %in% colnames(markers))) {
    return(invisible(NULL))
  }

  if ("p_val_adj" %in% colnames(markers)) {
    sig <- markers$gene[markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25]
  } else {
    sig <- markers$gene[abs(markers$avg_log2FC) > 0.25]
  }
  sig <- unique(sig[!is.na(sig)])
  if (length(sig) < 10) {
    message("Skipping GO enrichment for ", prefix, " because fewer than 10 genes are significant.")
    return(invisible(NULL))
  }

  converted <- clusterProfiler::bitr(
    sig,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  if (nrow(converted) < 10) {
    message("Skipping GO enrichment for ", prefix, " because fewer than 10 genes were converted.")
    return(invisible(NULL))
  }

  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )

  write.csv(as.data.frame(ego), file.path(out_dir, paste0(prefix, "_go_bp.csv")), row.names = FALSE)
}

stop_if_missing(input_rds)
message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

if (!target_gene %in% rownames(obj)) {
  stop(target_gene, " was not found in the Seurat object.")
}

analysis_summary <- safe_read_csv(analysis_summary_file)
manual_annotation <- safe_read_csv(manual_annotation_file)

stk31_high_clusters <- parse_semicolon_item(analysis_summary, "stk31_high_clusters", c("14", "0", "8", "7"))
nk_clusters <- parse_semicolon_item(analysis_summary, "nk_candidate_clusters", c("2", "3"))

if (!is.null(manual_annotation) && "cluster" %in% colnames(manual_annotation) && "manual_celltype" %in% colnames(manual_annotation)) {
  cluster_to_manual <- setNames(manual_annotation$manual_celltype, as.character(manual_annotation$cluster))
  obj$validation_celltype <- unname(cluster_to_manual[as.character(obj$seurat_clusters)])
} else if ("manual_celltype" %in% colnames(obj@meta.data)) {
  obj$validation_celltype <- as.character(obj$manual_celltype)
} else if ("broad_celltype" %in% colnames(obj@meta.data)) {
  obj$validation_celltype <- as.character(obj$broad_celltype)
} else {
  obj$validation_celltype <- as.character(obj$seurat_clusters)
}

obj$validation_celltype[is.na(obj$validation_celltype) | !nzchar(obj$validation_celltype)] <- "Unassigned"
obj$validation_celltype[as.character(obj$seurat_clusters) %in% nk_clusters] <- "NK_cell"
obj$validation_celltype <- factor(obj$validation_celltype)

obj$stk31_high_group <- ifelse(as.character(obj$seurat_clusters) %in% stk31_high_clusters, "STK31_high_cluster", "Other")
obj$nk_group <- ifelse(as.character(obj$seurat_clusters) %in% nk_clusters, "NK_cell", "Other")
obj$validation_group <- ifelse(obj$stk31_high_group == "STK31_high_cluster", "STK31_high_cluster", ifelse(obj$nk_group == "NK_cell", "NK_cell", "Other"))
obj$validation_group <- factor(obj$validation_group, levels = c("STK31_high_cluster", "NK_cell", "Other"))

cluster_annotation <- unique(data.frame(
  cluster = as.character(obj$seurat_clusters),
  validation_celltype = as.character(obj$validation_celltype),
  stk31_high_group = as.character(obj$stk31_high_group),
  nk_group = as.character(obj$nk_group),
  stringsAsFactors = FALSE
))
cluster_annotation <- cluster_annotation[order(as.numeric(cluster_annotation$cluster)), ]
write.csv(cluster_annotation, file.path(out_dir, "cluster_validation_annotation.csv"), row.names = FALSE)

celltype_counts <- as.data.frame(table(obj$validation_celltype), stringsAsFactors = FALSE)
colnames(celltype_counts) <- c("validation_celltype", "cell_count")
celltype_counts$pct_cells <- celltype_counts$cell_count / sum(celltype_counts$cell_count)
write.csv(celltype_counts, file.path(out_dir, "validation_celltype_counts.csv"), row.names = FALSE)

message("Writing validation UMAPs")
write_plot(
  file.path(out_dir, "umap_validation_celltypes.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "validation_celltype", label = TRUE, repel = TRUE) +
    ggtitle("Validation cell types")
)
write_plot(
  file.path(out_dir, "umap_validation_group.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "validation_group", label = TRUE, repel = TRUE) +
    ggtitle("STK31-high and NK_cell groups")
)
write_plot(
  file.path(out_dir, "umap_stk31_expression.pdf"),
  FeaturePlot(obj, features = target_gene, reduction = "umap", order = TRUE) +
    ggtitle("STK31 expression")
)

validation_feature_genes <- available_genes(obj, validation_feature_panel)
if (length(validation_feature_genes) > 0) {
  feature_plots <- lapply(validation_feature_genes, function(gene) {
    FeaturePlot(obj, features = gene, reduction = "umap", order = TRUE) + ggtitle(gene)
  })
  write_multi_page_pdf(
    file.path(out_dir, "validation_feature_umaps.pdf"),
    feature_plots,
    width = 7,
    height = 6
  )
}

identity_dotplot_genes <- unique(unlist(identity_marker_sets, use.names = FALSE))
identity_dotplot_genes <- available_genes(obj, identity_dotplot_genes)
if (length(identity_dotplot_genes) > 0) {
  write_plot(
    file.path(out_dir, "identity_marker_dotplot_by_celltype.pdf"),
    DotPlot(obj, features = identity_dotplot_genes, group.by = "validation_celltype") +
      RotatedAxis() +
      ggtitle("Identity markers by validation cell type"),
    width = 16,
    height = 7
  )
}

identity_marker_table <- summarize_gene_by_group(obj, identity_dotplot_genes, obj$validation_celltype, group_name = "validation_celltype")
if (nrow(identity_marker_table) > 0) {
  write.csv(identity_marker_table, file.path(out_dir, "identity_marker_expression_by_celltype.csv"), row.names = FALSE)
}

message("Scoring NK function signatures")
score_names <- character(0)
for (score_name in names(nk_function_sets)) {
  genes <- available_genes(obj, nk_function_sets[[score_name]])
  if (length(genes) < 2) {
    next
  }
  obj <- AddModuleScore(obj, features = list(genes), name = score_name, assay = DefaultAssay(obj))
  score_names <- c(score_names, paste0(score_name, "1"))
}

if (length(score_names) > 0) {
  score_rename <- setNames(names(nk_function_sets)[seq_along(score_names)], score_names)
  names(score_rename) <- score_names
  for (old_name in names(score_rename)) {
    obj[[score_rename[[old_name]]]] <- obj[[old_name]]
  }
}

nk_score_columns <- intersect(names(nk_function_sets), colnames(obj@meta.data))
if (length(nk_score_columns) > 0) {
  score_feature_plots <- lapply(nk_score_columns, function(score_name) {
    FeaturePlot(obj, features = score_name, reduction = "umap", order = TRUE) + ggtitle(score_name)
  })
  write_multi_page_pdf(
    file.path(out_dir, "nk_signature_umaps.pdf"),
    score_feature_plots,
    width = 7,
    height = 6
  )

  score_violin_plots <- lapply(nk_score_columns, function(score_name) {
    VlnPlot(obj, features = score_name, group.by = "validation_celltype", pt.size = 0.01) +
      RotatedAxis() +
      ggtitle(score_name)
  })
  write_multi_page_pdf(
    file.path(out_dir, "nk_signature_violin_by_celltype.pdf"),
    score_violin_plots,
    width = 10,
    height = 5
  )

  nk_score_summary <- do.call(rbind, lapply(nk_score_columns, function(score_name) {
    data.frame(
      validation_celltype = as.character(obj$validation_celltype),
      score_name = score_name,
      score_value = obj@meta.data[[score_name]],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(nk_score_summary, file.path(out_dir, "nk_signature_scores_by_cell.csv"), row.names = FALSE)

  nk_score_summary_by_type <- aggregate(
    score_value ~ validation_celltype + score_name,
    data = nk_score_summary,
    FUN = mean
  )
  write.csv(nk_score_summary_by_type, file.path(out_dir, "nk_signature_scores_by_celltype.csv"), row.names = FALSE)
}

message("Writing pathway evidence plots")
pathway_genes <- available_genes(obj, pathway_panel_genes)
if (length(pathway_genes) > 0) {
  write_plot(
    file.path(out_dir, "pathway_gene_dotplot_by_celltype.pdf"),
    DotPlot(obj, features = pathway_genes, group.by = "validation_celltype") +
      RotatedAxis() +
      ggtitle("Candidate pathway genes by validation cell type"),
    width = 18,
    height = 8
  )

  pathway_expression <- summarize_gene_by_group(obj, pathway_genes, obj$validation_celltype, group_name = "validation_celltype")
  write.csv(pathway_expression, file.path(out_dir, "pathway_gene_expression_by_celltype.csv"), row.names = FALSE)
}

if (file.exists(candidate_pathway_summary_file)) {
  pathway_summary <- read.csv(candidate_pathway_summary_file, stringsAsFactors = FALSE)
  write.csv(pathway_summary, file.path(out_dir, "candidate_pathway_summary_copy.csv"), row.names = FALSE)
}

if (file.exists(candidate_lr_file)) {
  candidate_lr <- read.csv(candidate_lr_file, stringsAsFactors = FALSE)
  write.csv(candidate_lr, file.path(out_dir, "candidate_ligand_receptor_pathways_copy.csv"), row.names = FALSE)

  top_lr <- candidate_lr[order(-candidate_lr$interaction_score), ]
  top_lr <- head(top_lr, min(20, nrow(top_lr)))
  if (nrow(top_lr) > 0) {
    write.csv(top_lr, file.path(out_dir, "top_candidate_interactions.csv"), row.names = FALSE)
  }
}

message("Reading differential expression tables")
de_stk31_vs_nk <- safe_read_csv(de_stk31_vs_nk_file)
de_stk31_vs_all <- safe_read_csv(de_stk31_vs_all_file)
de_nk_vs_all <- safe_read_csv(de_nk_vs_all_file)

if (!is.null(de_stk31_vs_nk) && nrow(de_stk31_vs_nk) > 0) {
  write.csv(top_de_genes(de_stk31_vs_nk, 40), file.path(out_dir, "top_stk31_vs_nk_genes.csv"), row.names = FALSE)
  run_optional_go(de_stk31_vs_nk, "stk31_high_vs_nk")
}
if (!is.null(de_stk31_vs_all) && nrow(de_stk31_vs_all) > 0) {
  write.csv(top_de_genes(de_stk31_vs_all, 40), file.path(out_dir, "top_stk31_vs_all_genes.csv"), row.names = FALSE)
  run_optional_go(de_stk31_vs_all, "stk31_high_vs_all_other")
}
if (!is.null(de_nk_vs_all) && nrow(de_nk_vs_all) > 0) {
  write.csv(top_de_genes(de_nk_vs_all, 40), file.path(out_dir, "top_nk_vs_all_genes.csv"), row.names = FALSE)
  run_optional_go(de_nk_vs_all, "nk_vs_all_other")
}

summary_table <- data.frame(
  item = c(
    "target_gene",
    "stk31_high_clusters",
    "stk31_high_cells",
    "nk_candidate_clusters",
    "nk_candidate_cells",
    "validation_groups",
    "validation_celltypes"
  ),
  value = c(
    target_gene,
    paste(stk31_high_clusters, collapse = ";"),
    sum(obj$stk31_high_group == "STK31_high_cluster"),
    paste(nk_clusters, collapse = ";"),
    sum(obj$nk_group == "NK_cell"),
    paste(levels(obj$validation_group), collapse = ";"),
    paste(levels(obj$validation_celltype), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)

saveRDS(obj, file.path(out_dir, "tissue2_stk31_nk_validation_seurat.rds"))

message("Done. Key outputs:")
message("- ", file.path(out_dir, "validation_summary.csv"))
message("- ", file.path(out_dir, "cluster_validation_annotation.csv"))
message("- ", file.path(out_dir, "identity_marker_dotplot_by_celltype.pdf"))
message("- ", file.path(out_dir, "nk_signature_umaps.pdf"))
message("- ", file.path(out_dir, "pathway_gene_dotplot_by_celltype.pdf"))
message("- ", file.path(out_dir, "top_stk31_vs_nk_genes.csv"))
