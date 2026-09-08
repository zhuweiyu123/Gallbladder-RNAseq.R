#!/usr/bin/env Rscript

# Complete the missing CNV-binary STK31 figures that require the full Seurat
# object and R/CellChat/clusterProfiler environment.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

pick_project_dir <- function() {
  candidates <- unique(c(
    Sys.getenv("CODEX_R_PROJECT_DIR", unset = NA_character_),
    getwd(),
    "/home/zhuweiyu/codex-r"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(file.path(
      candidate, "results", "merged_cnv_binary_annotation",
      "merged_cnv_binary_annotated.rds"
    ))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot find merged_cnv_binary_annotated.rds. Set CODEX_R_PROJECT_DIR.")
}

project_dir <- pick_project_dir()
root_dir <- file.path(project_dir, "results", "merged_cnv_binary_annotation")
object_path <- file.path(root_dir, "merged_cnv_binary_annotated.rds")
de_dir <- file.path(root_dir, "differential_expression")
figure_dir <- file.path(root_dir, "figures_stk31_de_go_cellchat")
table_dir <- file.path(root_dir, "tables_stk31_de_go_cellchat")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

padj_cut <- 0.05
logfc_cut <- 0.25
min_cells_cellchat <- as.integer(Sys.getenv("CELLCHAT_MIN_CELLS", "10"))
max_cells_cellchat <- as.integer(Sys.getenv("CELLCHAT_MAX_CELLS_PER_GROUP", "500"))

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

fetch_layer <- function(object, layer_name = "data", genes = NULL, cells = NULL) {
  mat <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = "RNA", layer = layer_name)
  } else {
    GetAssayData(object, assay = "RNA", slot = layer_name)
  }
  if (!is.null(genes)) genes <- intersect(genes, rownames(mat))
  if (is.null(genes)) genes <- rownames(mat)
  if (is.null(cells)) cells <- colnames(mat)
  mat[genes, cells, drop = FALSE]
}

zscore_rows <- function(mat) {
  t(scale(t(as.matrix(mat))))
}

message("Reading CNV-binary annotated Seurat object: ", object_path)
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj, assay = "RNA")
}
meta <- obj[[]]
need_cols <- c("patient_id", "epithelial_binary_label", "analysis_celltype_cnv_binary")
stop_if_not(all(need_cols %in% colnames(meta)),
            paste("Missing metadata columns:", paste(setdiff(need_cols, colnames(meta)), collapse = ", ")))

counts <- fetch_layer(obj, "counts")
data_mat <- fetch_layer(obj, "data")
stop_if_not("STK31" %in% rownames(counts), "STK31 missing from RNA counts")

cells <- colnames(obj)
stk31_counts <- as.numeric(counts["STK31", ])
malignant <- !is.na(meta$epithelial_binary_label) &
  meta$epithelial_binary_label == "Malignant epithelial cells"
normal <- !is.na(meta$epithelial_binary_label) &
  meta$epithelial_binary_label == "Normal epithelial cells"
stk31_high <- malignant & stk31_counts > 0
stk31_low <- malignant & stk31_counts == 0
nk_cells <- if ("analysis_celltype_previous_or_new" %in% colnames(meta)) {
  !is.na(meta$analysis_celltype_previous_or_new) &
    meta$analysis_celltype_previous_or_new == "NK_cell"
} else {
  !is.na(meta$analysis_celltype_cnv_binary) &
    meta$analysis_celltype_cnv_binary == "NK_cell"
}

validation <- data.table(
  check = c("malignant_epithelial", "normal_epithelial", "stk31_high_malignant",
            "stk31_low_malignant", "nk_cells", "p2_likely_in_malignant_rule"),
  observed = c(sum(malignant, na.rm = TRUE), sum(normal, na.rm = TRUE),
               sum(stk31_high, na.rm = TRUE), sum(stk31_low, na.rm = TRUE),
               sum(nk_cells, na.rm = TRUE), "malignant_likely encoded as malignant"),
  expected = c(1047, 3716, 117, 930, 1479, "malignant_likely encoded as malignant")
)
validation[, pass := observed == expected]
fwrite(validation, file.path(table_dir, "R_object_missing_figure_validation.tsv"), sep = "\t")
stop_if_not(all(validation$pass), "Object validation failed")

group_label <- rep(NA_character_, ncol(obj))
group_label[stk31_high] <- "STK31_high_malignant_epithelial"
group_label[stk31_low] <- "STK31_low_malignant_epithelial"
group_label[normal] <- "Normal_epithelial"
group_label[nk_cells] <- "NK_cell"
obj$stk31_cnv_binary_plot_group <- group_label

# ------------------------------------------------------------
# Object-level immune target dotplot
# ------------------------------------------------------------
immune_genes <- unique(c(
  "STK31", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "TAP1", "TAP2", "NLRC5",
  "TIGIT", "PVR", "NECTIN2", "CD274", "PDCD1LG2",
  "TGFB1", "TGFBR1", "TGFBR2", "IFNG", "IFNGR1", "IFNGR2",
  "NKG7", "GNLY", "PRF1", "GZMB",
  "KLRC1", "KLRD1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3",
  "KLRK1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3"
))
immune_genes <- intersect(immune_genes, rownames(data_mat))
plot_groups <- c(
  "STK31_high_malignant_epithelial", "STK31_low_malignant_epithelial",
  "Normal_epithelial", "NK_cell"
)
dot_rows <- rbindlist(lapply(plot_groups, function(group) {
  idx <- which(obj$stk31_cnv_binary_plot_group == group)
  if (length(idx) == 0) return(NULL)
  rbindlist(lapply(immune_genes, function(gene) {
    values_data <- as.numeric(data_mat[gene, idx])
    values_counts <- as.numeric(counts[gene, idx])
    data.table(
      group = group,
      gene = gene,
      cells = length(idx),
      detected_fraction = mean(values_counts > 0),
      mean_log_normalized = mean(values_data)
    )
  }))
}))
fwrite(dot_rows, file.path(table_dir, "Dotplot_object_immune_targets_table.csv"))
dot_rows[, group := factor(group, levels = plot_groups)]
dot_rows[, gene := factor(gene, levels = rev(immune_genes))]
p_dot <- ggplot(dot_rows, aes(x = group, y = gene)) +
  geom_point(aes(size = detected_fraction, color = mean_log_normalized)) +
  scale_size_continuous(range = c(0.8, 8),
                        labels = function(x) paste0(round(100 * x), "%")) +
  scale_color_gradient(low = "#d9f0f3", high = "#b2182b") +
  theme_classic(base_size = 10) +
  labs(
    title = "Immune-related target genes in CNV-binary STK31 groups",
    subtitle = "P2 malignant_likely is included as malignant; STK31-high = raw count > 0",
    x = NULL, y = NULL, size = "% detected", color = "Mean log-normalized"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank())
save_both(p_dot, "Dotplot_object_immune_targets_cnv_binary", 9, 7.5)

# ------------------------------------------------------------
# Top expression heatmap for STK31-high vs low malignant cells
# ------------------------------------------------------------
de_hl <- fread(file.path(de_dir, "markers_stk31_high_vs_low_malignant_epithelial.csv"))
top_up <- head(de_hl[avg_log2FC > 0][order(p_val_adj, -abs(avg_log2FC)), gene], 15)
top_down <- head(de_hl[avg_log2FC < 0][order(p_val_adj, -abs(avg_log2FC)), gene], 15)
heat_genes <- intersect(unique(c(top_up, top_down)), rownames(data_mat))
avg_rows <- rbindlist(lapply(unique(meta$patient_id[malignant]), function(patient) {
  rbindlist(lapply(c("STK31_high", "STK31_low"), function(status) {
    status_mask <- if (status == "STK31_high") stk31_high else stk31_low
    idx <- which(meta$patient_id == patient & status_mask)
    if (length(idx) == 0) return(NULL)
    values <- Matrix::rowMeans(data_mat[heat_genes, idx, drop = FALSE])
    data.table(patient_id = patient, stk31_group = status, gene = heat_genes, mean_log_normalized = values)
  }))
}))
fwrite(avg_rows, file.path(table_dir, "Heatmap_STK31_high_low_malignant_patient_average_table.csv"))
avg_wide <- dcast(avg_rows, gene ~ patient_id + stk31_group, value.var = "mean_log_normalized")
hm <- as.matrix(avg_wide[, -1, with = FALSE])
rownames(hm) <- avg_wide$gene
hm_z <- zscore_rows(hm)
hm_long <- as.data.table(as.table(hm_z))
colnames(hm_long) <- c("gene", "sample_group", "z")
hm_long[, gene := factor(gene, levels = rev(heat_genes))]
hm_long[, sample_group := factor(sample_group, levels = colnames(hm))]
p_hm <- ggplot(hm_long, aes(x = sample_group, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, na.value = "#f0f0f0") +
  theme_classic(base_size = 10) +
  labs(
    title = "Top STK31-high vs low malignant DE genes",
    subtitle = "Patient-level average log-normalized expression; row z-score",
    x = NULL, y = NULL, fill = "Row z"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.ticks = element_blank())
save_both(p_hm, "Heatmap_object_STK31_high_low_malignant_top_DE", 8.5, 8)

# ------------------------------------------------------------
# Immune module scores
# ------------------------------------------------------------
modules <- list(
  Antigen_presentation_MHCI = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "TAP1", "TAP2", "NLRC5"),
  IFN_response = c("IFNG", "IFNGR1", "IFNGR2", "STAT1", "STAT2", "IRF1", "IRF7", "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"),
  TIGIT_PVR_checkpoint = c("TIGIT", "PVR", "NECTIN2", "CD226", "CD96", "PVRIG"),
  TGFb_axis = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7"),
  NK_cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CTSW"),
  NK_inhibitory_receptors = c("KLRC1", "KLRD1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "LILRB1"),
  NKG2D_ligands = c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6")
)
score_rows <- rbindlist(lapply(names(modules), function(module_name) {
  genes <- intersect(modules[[module_name]], rownames(data_mat))
  if (length(genes) == 0) {
    return(data.table(
      module = character(), group = character(), patient_id = character(),
      cell = character(), score = numeric(), n_genes = integer(),
      genes_used = character()
    ))
  }
  eligible <- which(!is.na(obj$stk31_cnv_binary_plot_group))
  scores <- Matrix::colMeans(data_mat[genes, eligible, drop = FALSE])
  data.table(
    module = module_name,
    group = obj$stk31_cnv_binary_plot_group[eligible],
    patient_id = meta$patient_id[eligible],
    cell = colnames(obj)[eligible],
    score = as.numeric(scores),
    n_genes = length(genes),
    genes_used = paste(genes, collapse = ";")
  )
}), fill = TRUE)
fwrite(score_rows, file.path(table_dir, "ModuleScore_immune_pathways_by_group_table.csv"))
score_rows[, group := factor(group, levels = plot_groups)]
p_score <- ggplot(score_rows, aes(x = group, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.2, fill = "white") +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(
    "STK31_high_malignant_epithelial" = "#b2182b",
    "STK31_low_malignant_epithelial" = "#ef8a62",
    "Normal_epithelial" = "#2166ac",
    "NK_cell" = "#053061"
  ), drop = FALSE) +
  theme_classic(base_size = 9.5) +
  labs(
    title = "Immune pathway/module scores by CNV-binary STK31 group",
    subtitle = "Mean log-normalized expression of available genes in each module",
    x = NULL, y = "Module score"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
save_both(p_score, "ModuleScore_immune_pathways_by_group", 11, 9)

# ------------------------------------------------------------
# GO / GSEA
# ------------------------------------------------------------
write_go_skip <- function(prefix, note) {
  fwrite(data.table(note = note), file.path(table_dir, paste0(prefix, "_GO_GSEA_skipped.csv")))
}

run_go_gsea <- function(de_path, prefix, fc_col, p_col, padj_col) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    write_go_skip(prefix, "clusterProfiler/org.Hs.eg.db not installed")
    return(invisible(NULL))
  }
  de <- fread(de_path)
  de <- de[is.finite(get(fc_col)) & is.finite(get(p_col))]
  de[, rank_score := sign(get(fc_col)) * -log10(pmax(get(p_col), 1e-300))]
  de <- de[order(-abs(rank_score))]
  de <- de[!duplicated(gene)]
  converted <- suppressMessages(clusterProfiler::bitr(
    de$gene, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  ranked <- merge(de[, .(gene, rank_score)], converted,
                  by.x = "gene", by.y = "SYMBOL", all = FALSE)
  ranked <- ranked[!duplicated(ENTREZID)]
  stats <- ranked$rank_score
  names(stats) <- ranked$ENTREZID
  stats <- sort(stats, decreasing = TRUE)
  fwrite(ranked, file.path(table_dir, paste0(prefix, "_ranked_genes_entrez.csv")))

  gsea <- tryCatch(
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
      ont = "BP", pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(gsea, "error")) {
    write_go_skip(prefix, paste("gseGO failed:", conditionMessage(gsea)))
  } else {
    gsea_df <- as.data.frame(gsea)
    fwrite(gsea_df, file.path(table_dir, paste0(prefix, "_gseGO_BP.csv")))
    if (nrow(gsea_df) > 0) {
      top <- head(gsea_df[order(gsea_df$p.adjust), ], 20)
      top$Description <- factor(top$Description, levels = rev(top$Description))
      p <- ggplot(top, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(low = "#b2182b", high = "#2166ac", trans = "log10") +
        theme_classic(base_size = 10) +
        labs(title = paste0("GSEA GO BP: ", prefix), x = "NES", y = NULL,
             color = "FDR", size = "Set size")
      save_both(p, paste0("GSEA_GO_BP_", prefix), 9.5, max(5, 0.28 * nrow(top) + 2))
    }
  }

  sig <- de[get(padj_col) < padj_cut & abs(get(fc_col)) > logfc_cut, gene]
  fwrite(data.table(prefix = prefix, n_sig = length(sig), padj_cut = padj_cut,
                    logfc_cut = logfc_cut),
         file.path(table_dir, paste0(prefix, "_ORA_input_summary.csv")))
  if (length(sig) < 10) {
    write_go_skip(paste0(prefix, "_ORA"), "fewer than 10 FDR/logFC significant genes")
    return(invisible(NULL))
  }
  sig_conv <- suppressMessages(clusterProfiler::bitr(
    unique(sig), fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  ora <- clusterProfiler::enrichGO(
    gene = unique(sig_conv$ENTREZID), OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  ora_df <- as.data.frame(ora)
  fwrite(ora_df, file.path(table_dir, paste0(prefix, "_enrichGO_BP.csv")))
  if (nrow(ora_df) > 0) {
    top <- head(ora_df[order(ora_df$p.adjust), ], 20)
    top$Description <- factor(top$Description, levels = rev(top$Description))
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#b2182b", high = "#2166ac", trans = "log10") +
      theme_classic(base_size = 10) +
      labs(title = paste0("ORA GO BP: ", prefix), x = "-log10(FDR)", y = NULL,
           color = "FDR", size = "Genes")
    save_both(p, paste0("ORA_GO_BP_", prefix), 9.5, max(5, 0.28 * nrow(top) + 2))
  }
  invisible(NULL)
}

run_go_gsea(
  file.path(de_dir, "pseudobulk_edgeR_paired_malignant_vs_normal.csv"),
  "paired_malignant_vs_normal_epithelial", "logFC", "PValue", "FDR"
)
run_go_gsea(
  file.path(de_dir, "markers_stk31_high_vs_low_malignant_epithelial.csv"),
  "stk31_high_vs_low_malignant_epithelial", "avg_log2FC", "p_val", "p_val_adj"
)

# ------------------------------------------------------------
# CellChat
# ------------------------------------------------------------
axis_genes <- list(
  MHC_I_KIR = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "KIR", "KLRC1", "KLRD1"),
  TIGIT_PVR_NECTIN = c("TIGIT", "PVR", "NECTIN", "CD226", "CD96"),
  TGFb = c("TGFB", "TGFBR"),
  NKG2D = c("NKG2D", "KLRK1", "MICA", "MICB", "ULBP"),
  IFNG_IFN = c("IFNG", "IFNGR", "IFN")
)

summarize_cellchat_axes <- function(comm) {
  if (nrow(comm) == 0) return(data.table())
  text_cols <- intersect(c("pathway_name", "pathway", "ligand", "receptor", "interaction_name"), names(comm))
  comm$axis_text <- apply(comm[, text_cols, drop = FALSE], 1, paste, collapse = " ")
  target_pairs <- rbind(
    data.table(source = "STK31_high_malignant_epithelial", target = "NK_cell", direction = "high_to_NK"),
    data.table(source = "STK31_low_malignant_epithelial", target = "NK_cell", direction = "low_to_NK"),
    data.table(source = "NK_cell", target = "STK31_high_malignant_epithelial", direction = "NK_to_high"),
    data.table(source = "NK_cell", target = "STK31_low_malignant_epithelial", direction = "NK_to_low")
  )
  out <- rbindlist(lapply(seq_len(nrow(target_pairs)), function(i) {
    sub <- comm[comm$source == target_pairs$source[i] & comm$target == target_pairs$target[i], ]
    rbindlist(lapply(names(axis_genes), function(axis) {
      pattern <- paste(axis_genes[[axis]], collapse = "|")
      hit <- grepl(pattern, sub$axis_text, ignore.case = TRUE)
      data.table(
        direction = target_pairs$direction[i],
        axis = axis,
        total_probability = sum(sub$prob[hit], na.rm = TRUE),
        ligand_receptor_rows = sum(hit)
      )
    }))
  }))
  out
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  fwrite(data.table(note = "CellChat not installed"),
         file.path(table_dir, "CellChat_cnv_binary_skipped.csv"))
} else {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan("sequential")
    options(future.globals.maxSize = 8 * 1024^3)
  }
  cc_group <- as.character(meta$analysis_celltype_cnv_binary)
  cc_group[stk31_high] <- "STK31_high_malignant_epithelial"
  cc_group[stk31_low] <- "STK31_low_malignant_epithelial"
  cc_group[normal] <- "Normal_epithelial"
  cc_group[nk_cells] <- "NK_cell"
  cc_group[cc_group %in% c("T_NK_ambiguous", "T_NK_unresolved_new", "Ambiguous", "Unassigned")] <-
    "Other_immune_or_unassigned"
  obj$cellchat_cnv_binary_group <- cc_group
  group_counts <- as.data.table(table(group = obj$cellchat_cnv_binary_group))
  group_counts <- group_counts[N >= min_cells_cellchat]
  keep_cells <- cells[obj$cellchat_cnv_binary_group %in% group_counts$group]
  split_cells <- split(keep_cells, obj$cellchat_cnv_binary_group[match(keep_cells, cells)])
  keep_cells <- unlist(lapply(split_cells, function(x) {
    if (length(x) > max_cells_cellchat) sample(x, max_cells_cellchat) else x
  }), use.names = FALSE)
  cc_data <- fetch_layer(obj, "data", cells = keep_cells)
  cc_meta <- data.frame(
    labels = factor(obj$cellchat_cnv_binary_group[match(keep_cells, cells)]),
    row.names = keep_cells,
    stringsAsFactors = FALSE
  )
  fwrite(as.data.table(table(group = cc_meta$labels)),
         file.path(table_dir, "CellChat_cnv_binary_group_counts.csv"))

  message("Running CellChat on CNV-binary STK31 groups")
  cellchat <- CellChat::createCellChat(object = cc_data, meta = cc_meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_cellchat)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(table_dir, "CellChat_cnv_binary_object.rds"))
  comm <- CellChat::subsetCommunication(cellchat)
  fwrite(as.data.table(comm), file.path(table_dir, "CellChat_cnv_binary_communications.csv"))

  group_size <- as.numeric(table(cc_meta$labels))
  names(group_size) <- names(table(cc_meta$labels))
  pdf(file.path(figure_dir, "CellChat_cnv_binary_circle_count.pdf"), width = 8, height = 8)
  CellChat::netVisual_circle(cellchat@net$count, vertex.weight = group_size,
                             weight.scale = TRUE, label.edge = FALSE,
                             title.name = "Number of interactions")
  dev.off()
  pdf(file.path(figure_dir, "CellChat_cnv_binary_circle_weight.pdf"), width = 8, height = 8)
  CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = group_size,
                             weight.scale = TRUE, label.edge = FALSE,
                             title.name = "Interaction weights")
  dev.off()
  pdf(file.path(figure_dir, "CellChat_cnv_binary_heatmap_weight.pdf"), width = 8, height = 7)
  print(CellChat::netVisual_heatmap(cellchat, measure = "weight"))
  dev.off()

  focus <- intersect(
    c("STK31_high_malignant_epithelial", "STK31_low_malignant_epithelial", "NK_cell"),
    levels(cc_meta$labels)
  )
  if (length(focus) >= 2) {
    pdf(file.path(figure_dir, "CellChat_cnv_binary_focused_STK31_NK_bubble.pdf"),
        width = 11, height = 8)
    print(CellChat::netVisual_bubble(cellchat, sources.use = focus, targets.use = focus,
                                     remove.isolate = FALSE))
    dev.off()
  }

  axis_df <- summarize_cellchat_axes(comm)
  fwrite(axis_df, file.path(table_dir, "CellChat_cnv_binary_candidate_axis_summary.csv"))
  if (nrow(axis_df) > 0) {
    p_axis <- ggplot(axis_df, aes(x = direction, y = axis, size = ligand_receptor_rows,
                                  color = total_probability)) +
      geom_point() +
      scale_color_gradient(low = "#d9f0f3", high = "#b2182b") +
      theme_classic(base_size = 10) +
      labs(
        title = "CellChat candidate immune-axis summary",
        subtitle = "Focused on STK31-high/low malignant epithelial cells and NK cells",
        x = NULL, y = NULL, size = "LR rows", color = "Total probability"
      ) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_both(p_axis, "CellChat_cnv_binary_candidate_axis_summary", 8.5, 5.5)
  }
}

writeLines(c(
  "# CNV-binary missing figure rerun",
  "",
  paste0("Date: ", Sys.time()),
  paste0("Project directory: ", project_dir),
  "",
  "P2 malignant_likely is included as malignant through the frozen CNV-binary object.",
  "Cell-level Wilcoxon and CellChat outputs remain exploratory because cells are not independent patient replicates.",
  "For malignant vs normal epithelial, paired pseudobulk has no FDR<0.05 hits; ranked GSEA is the preferred enrichment display."
), file.path(root_dir, "STK31_CNV_BINARY_DE_GO_CELLCAT_RUN.md"))

message("Done. Outputs written under: ", figure_dir, " and ", table_dir)
