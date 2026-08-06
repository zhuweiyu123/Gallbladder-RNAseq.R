# ============================================================
# 06_refine_t_nk_and_rerun_stk31_nk.R
# ============================================================
# PI task: strict T/NK re-annotation, then NK downstream + CellChat
# high-vs-low differential interaction export.
#
# Reuses small helper patterns from scripts/02_merged_basic_analysis.R
# (fetch_gene_matrix, STK31 high/low assignment, LR table, module scores).
# Does NOT overwrite existing result directories.
#
# Stage gate: Stage 2/3 run only if Stage 1 QC passes.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260717)

# -------------------- paths / params --------------------
project_dir <- "/home/zhuweiyu/codex-r"
input_rds <- file.path(
  project_dir, "results/merged_basic_seurat/gallbladder_cancer_merged_basic_seurat.rds"
)
old_annotation_file <- file.path(
  project_dir, "results/merged_stk31_nk_analysis/cluster_celltype_annotation.csv"
)
old_summary_file <- file.path(
  project_dir, "results/merged_stk31_nk_analysis/analysis_summary.csv"
)

stage1_dir <- file.path(project_dir, "results/merged_tnk_refined_annotation")
stage2_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
stage3_dir <- file.path(project_dir, "results/merged_cellchat_refined_high_low_nk")
dir.create(stage1_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stage2_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stage3_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
tumor_epithelial_celltype_label <- "Epithelial"
tumor_epithelial_stk31_high_quantile <- 0.75
min_cells_for_de <- 10
min_cells_per_cellchat_group <- 10
max_cells_per_cellchat_group <- 1500
run_exact_wilcox <- TRUE
# Speed control for large DE (without presto). Descriptive only.
max_cells_per_ident_de <- as.integer(Sys.getenv("STK31_MAX_CELLS_DE", "500"))
# Resume control: "1" full run; "2" skip stage1 if QC passed; "3" stage3 only
stage_start <- as.integer(Sys.getenv("STK31_STAGE_START", "1"))

# Subclustering
tnk_resolution <- 0.4
tnk_n_variable <- 2000
tnk_n_pcs <- 30
tnk_use_dims <- 1:20
nk_score_cluster_quantile <- 0.85
score_delta <- 0.05

t_markers <- c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD2", "CD247")
nk_markers <- c(
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "KLRD1", "KLRF1",
  "FCGR3A", "NCAM1", "TYROBP", "CTSW", "CST7", "XCL1", "XCL2"
)
core_nk_markers <- c("NKG7", "GNLY", "PRF1", "KLRD1", "KLRF1", "FCGR3A", "NCAM1")
core_t_markers <- c("CD3D", "CD3E", "TRAC")
exclude_markers <- c(
  "MS4A1", "CD79A", "LYZ", "S100A8", "S100A9",
  "EPCAM", "KRT8", "KRT18", "COL1A1", "PECAM1"
)

nk_function_sets <- list(
  NK_cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH"),
  NK_activation = c("IFNG", "TNF", "CSF2", "CCL3", "CCL4", "CCL5"),
  NK_checkpoint = c("TIGIT", "HAVCR2", "LAG3", "PDCD1", "CD96"),
  NK_migration = c("CXCR3", "CXCR4", "CCR5", "S1PR5", "ITGAL"),
  TGFbeta_response = c("TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SKIL")
)

lr_reference <- data.frame(
  ligand = c(
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "XCL1", "XCL2", "IL15", "IL18",
    "IL12A", "IL12B", "TGFB1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
    "HLA-A", "HLA-B", "HLA-C", "CD274", "LGALS9", "NECTIN2", "PVR", "ICAM1",
    "TNFSF10", "FASLG", "IFNG", "LTA", "CSF2"
  ),
  receptor = c(
    "CXCR3", "CXCR3", "CXCR3", "CCR5", "XCR1", "XCR1", "IL2RB", "IL18R1",
    "IL12RB1", "IL12RB1", "TGFBR2", "KLRK1", "KLRK1", "KLRK1", "KLRK1", "KLRK1",
    "KIR2DL1", "KIR3DL1", "KIR2DL1", "PDCD1", "HAVCR2", "TIGIT", "TIGIT", "ITGAL",
    "TNFRSF10B", "FAS", "IFNGR1", "TNFRSF1A", "CSF2RA"
  ),
  pathway = c(
    "CXCR3 chemotaxis", "CXCR3 chemotaxis", "CXCR3 chemotaxis", "CCR5 migration",
    "XCR1 recruitment", "XCR1 recruitment", "IL15 NK activation", "IL18 IFNG induction",
    "IL12 cytotoxic activation", "IL12 cytotoxic activation", "TGF beta suppression",
    "NKG2D stress ligand", "NKG2D stress ligand", "NKG2D stress ligand", "NKG2D stress ligand",
    "NKG2D stress ligand", "MHC I inhibitory KIR", "MHC I inhibitory KIR", "MHC I inhibitory KIR",
    "PD1 checkpoint", "TIM3 checkpoint", "TIGIT checkpoint", "TIGIT checkpoint", "LFA1 adhesion",
    "TRAIL apoptosis", "FAS apoptosis", "IFNG response", "TNF family", "GM-CSF signaling"
  ),
  stringsAsFactors = FALSE
)

mechanism_axes <- list(
  "MHC-I / HLA axis" = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "KIR2DL1", "KIR3DL1", "KIR2DL3"),
  "TIGIT / NECTIN axis" = c("NECTIN2", "NECTIN3", "PVR", "TIGIT", "CD96"),
  "TGFb axis" = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2"),
  "IFNG / IFN axis" = c("IFNG", "IFNGR1", "IFNGR2", "STAT1"),
  "NKG2D axis" = c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "KLRK1")
)

mhc_i_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M")

# -------------------- helpers (from 02 patterns) --------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

available_genes <- function(object, genes) intersect(genes, rownames(object))

fetch_gene_matrix <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NULL)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

write_base_pdf <- function(path, expr, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  force(expr)
  dev.off()
}

module_score <- function(object, genes) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(rep(NA_real_, ncol(object)))
  as.numeric(Matrix::colMeans(fetch_gene_matrix(object, genes)))
}

gene_mean <- function(object, genes, cells = NULL) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NA_real_)
  mat <- fetch_gene_matrix(object, genes)
  if (!is.null(cells)) mat <- mat[, intersect(cells, colnames(mat)), drop = FALSE]
  if (ncol(mat) == 0) return(NA_real_)
  mean(as.numeric(Matrix::colMeans(mat)))
}

read_summary_item <- function(path, item_name, fallback = character(0)) {
  if (!file.exists(path)) return(fallback)
  summary_df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("item", "value") %in% colnames(summary_df))) return(fallback)
  idx <- which(summary_df$item == item_name)
  if (length(idx) == 0) return(fallback)
  value <- summary_df$value[idx[1]]
  if (is.na(value) || !nzchar(as.character(value))) return(fallback)
  unlist(strsplit(as.character(value), ";", fixed = TRUE))
}

assign_tumor_epithelial_stk31_group <- function(object, stk31_expr, celltype_col, celltype_label, high_quantile = 0.75) {
  tumor_epithelial_idx <- as.character(object[[celltype_col]][, 1]) == celltype_label
  if (sum(tumor_epithelial_idx) == 0) {
    stop("No tumor epithelial cells found with ", celltype_col, " == ", celltype_label)
  }
  tumor_epithelial_expr <- stk31_expr[tumor_epithelial_idx]
  cutoff <- as.numeric(stats::quantile(tumor_epithelial_expr, high_quantile, na.rm = TRUE))
  if (!is.finite(cutoff)) cutoff <- 0
  if (cutoff > 0) {
    high_idx <- tumor_epithelial_idx & stk31_expr >= cutoff
  } else {
    high_idx <- tumor_epithelial_idx & stk31_expr > 0
  }
  low_idx <- tumor_epithelial_idx & !high_idx
  group <- rep("Other", length(stk31_expr))
  group[high_idx] <- "STK31_high_tumor_epithelial"
  group[low_idx] <- "STK31_low_tumor_epithelial"
  names(group) <- colnames(object)
  list(
    group = group, cutoff = cutoff,
    tumor_epithelial_cells = sum(tumor_epithelial_idx),
    high_cells = sum(high_idx), low_cells = sum(low_idx)
  )
}

run_marker_test <- function(object, cells.1, cells.2, output_file) {
  if (length(cells.1) < min_cells_for_de || length(cells.2) < min_cells_for_de) {
    warning("Too few cells for marker test: ", output_file)
    empty <- data.frame()
    write.csv(empty, output_file, row.names = FALSE)
    return(empty)
  }
  object$.comparison_group <- "unused"
  object$.comparison_group[colnames(object) %in% cells.1] <- "group_1"
  object$.comparison_group[colnames(object) %in% cells.2] <- "group_2"
  markers <- FindMarkers(
    object, ident.1 = "group_1", ident.2 = "group_2",
    group.by = ".comparison_group", logfc.threshold = 0.1, min.pct = 0.1,
    test.use = "wilcox",
    max.cells.per.ident = max_cells_per_ident_de
  )
  markers$gene <- rownames(markers)
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
  write.csv(markers, output_file, row.names = FALSE)
  markers
}

average_expression_for_genes <- function(object, cells, genes) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0 || length(cells) == 0) {
    return(data.frame(gene = character(0), avg_expr = numeric(0), pct_expr = numeric(0)))
  }
  expr <- fetch_gene_matrix(object, genes)[genes, cells, drop = FALSE]
  data.frame(
    gene = genes,
    avg_expr = as.numeric(Matrix::rowMeans(expr)),
    pct_expr = as.numeric(Matrix::rowMeans(expr > 0)),
    stringsAsFactors = FALSE
  )
}

infer_lr_links <- function(object, high_cells, nk_cells, lr_table) {
  lr_genes <- unique(c(lr_table$ligand, lr_table$receptor))
  high_expr <- average_expression_for_genes(object, high_cells, lr_genes)
  nk_expr <- average_expression_for_genes(object, nk_cells, lr_genes)
  names(high_expr)[names(high_expr) != "gene"] <- paste0(
    "tumor_epithelial_stk31_high_", names(high_expr)[names(high_expr) != "gene"]
  )
  names(nk_expr)[names(nk_expr) != "gene"] <- paste0(
    "nk_group_", names(nk_expr)[names(nk_expr) != "gene"]
  )
  forward <- merge(lr_table, high_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  forward <- merge(forward, nk_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  forward$direction <- "STK31_high_tumor_epithelial_to_NK_cell"
  reverse <- merge(lr_table, nk_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  reverse <- merge(reverse, high_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  reverse$direction <- "NK_cell_to_STK31_high_tumor_epithelial"
  reverse <- reverse[, names(forward)]
  links <- rbind(forward, reverse)
  links[is.na(links)] <- 0
  links$interaction_score <- sqrt(
    links$tumor_epithelial_stk31_high_avg_expr * links$nk_group_avg_expr
  )
  links[order(-links$interaction_score), ]
}

label_evidence <- function(prob_high, prob_low, similar_ratio = 1.25, min_prob = 1e-6) {
  if ((is.na(prob_high) || prob_high < min_prob) && (is.na(prob_low) || prob_low < min_prob)) {
    return("not_detected")
  }
  ph <- ifelse(is.na(prob_high), 0, prob_high)
  pl <- ifelse(is.na(prob_low), 0, prob_low)
  if (ph >= pl * similar_ratio && ph > pl) return("higher_in_high")
  if (pl >= ph * similar_ratio && pl > ph) return("higher_in_low")
  "similar"
}

write_session_info <- function(path) {
  sink(path)
  cat("timestamp:", as.character(Sys.time()), "\n")
  cat("script: 06_refine_t_nk_and_rerun_stk31_nk.R\n")
  cat("project_dir:", project_dir, "\n")
  cat("input_rds:", input_rds, "\n")
  cat("tnk_resolution:", tnk_resolution, "\n")
  cat("score_delta:", score_delta, "\n")
  cat("tumor_epithelial_stk31_high_quantile:", tumor_epithelial_stk31_high_quantile, "\n")
  print(sessionInfo())
  sink()
}

list_output_files <- function(dirs) {
  files <- unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, recursive = TRUE, full.names = TRUE)
  }))
  data.frame(file = files, stringsAsFactors = FALSE)
}

# -------------------- start --------------------
message("=== 06 STK31 T/NK refine + rerun ===")
message("Start: ", Sys.time())
message("stage_start=", stage_start, " max_cells_per_ident_de=", max_cells_per_ident_de)
write_session_info(file.path(stage1_dir, "sessionInfo_start.txt"))
stop_if_missing(input_rds)
stop_if_missing(old_annotation_file)

cluster_annotation <- read.csv(old_annotation_file, stringsAsFactors = FALSE)
cluster_to_manual <- setNames(cluster_annotation$manual_celltype, as.character(cluster_annotation$cluster))
cluster_to_broad <- setNames(cluster_annotation$broad_celltype, as.character(cluster_annotation$cluster))
old_nk_clusters <- read_summary_item(old_summary_file, "nk_candidate_clusters", fallback = c("0", "1"))

resume_rds <- file.path(stage1_dir, "refined_merged_object.rds")
qc_file <- file.path(stage1_dir, "tnk_refinement_qc_summary.csv")

if (stage_start >= 2) {
  stop_if_missing(resume_rds)
  stop_if_missing(qc_file)
  message("Resuming from Stage ", stage_start, " using ", resume_rds)
  obj <- readRDS(resume_rds)
  DefaultAssay(obj) <- "RNA"
  old_nk_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% old_nk_clusters]
  qc_exist <- read.csv(qc_file, stringsAsFactors = FALSE)
  stage1_pass <- tolower(qc_exist$value[qc_exist$item == "stage1_pass"][1]) == "true"
  high_nk_cells <- colnames(obj)[obj$analysis_celltype == "NK_cell"]
  message("Loaded refined object. stage1_pass=", stage1_pass, " high_conf_nk=", length(high_nk_cells))
  if (!stage1_pass) {
    stop("Cannot resume Stage 2/3 because stage1_pass is FALSE")
  }
} else {
message("Loading merged Seurat object")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  if ("layers" %in% slotNames(obj[["RNA"]])) {
    try(obj <- JoinLayers(obj, assay = "RNA"), silent = TRUE)
  }
}

obj$legacy_manual_celltype <- unname(cluster_to_manual[as.character(obj$seurat_clusters)])
obj$legacy_broad_celltype <- unname(cluster_to_broad[as.character(obj$seurat_clusters)])
obj$legacy_manual_celltype[is.na(obj$legacy_manual_celltype)] <- "Unassigned"
obj$legacy_broad_celltype[is.na(obj$legacy_broad_celltype)] <- "Unassigned"

old_nk_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% old_nk_clusters]
message("Legacy NK clusters: ", paste(old_nk_clusters, collapse = ","),
        " n=", length(old_nk_cells))

# ============================================================
# STAGE 1: strict T/NK re-annotation
# ============================================================
message("=== STAGE 1: T/NK refinement ===")

# Score NK-ness on all clusters to expand candidates beyond 0/1/5
obj$global_nk_score <- module_score(obj, nk_markers)
obj$global_t_score <- module_score(obj, t_markers)

cluster_ids <- sort(unique(as.character(obj$seurat_clusters)))
cluster_nk_tbl <- do.call(rbind, lapply(cluster_ids, function(cl) {
  idx <- as.character(obj$seurat_clusters) == cl
  data.frame(
    cluster = cl,
    cells = sum(idx),
    mean_nk_score = mean(obj$global_nk_score[idx], na.rm = TRUE),
    mean_t_score = mean(obj$global_t_score[idx], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
cluster_nk_tbl <- cluster_nk_tbl[order(-cluster_nk_tbl$mean_nk_score), ]
write.csv(cluster_nk_tbl, file.path(stage1_dir, "global_cluster_tnk_scores.csv"), row.names = FALSE)

nk_cut <- as.numeric(stats::quantile(cluster_nk_tbl$mean_nk_score, nk_score_cluster_quantile, na.rm = TRUE))
candidate_clusters <- unique(c(
  "0", "1", "5",
  cluster_nk_tbl$cluster[cluster_nk_tbl$mean_nk_score >= nk_cut],
  old_nk_clusters
))
candidate_clusters <- candidate_clusters[candidate_clusters %in% cluster_ids]
message("Candidate T/NK clusters: ", paste(candidate_clusters, collapse = ","))

candidate_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% candidate_clusters]
if (length(candidate_cells) < 100) stop("Too few candidate T/NK cells: ", length(candidate_cells))

tnk <- subset(obj, cells = candidate_cells)
# Reprocess subset
message("Reclustering candidate T/NK cells: ", ncol(tnk))
tnk <- NormalizeData(tnk, verbose = FALSE)
tnk <- FindVariableFeatures(tnk, nfeatures = tnk_n_variable, verbose = FALSE)
tnk <- ScaleData(tnk, verbose = FALSE)
tnk <- RunPCA(tnk, npcs = tnk_n_pcs, verbose = FALSE)
tnk <- FindNeighbors(tnk, dims = tnk_use_dims, verbose = FALSE)
tnk <- FindClusters(tnk, resolution = tnk_resolution, verbose = FALSE)
tnk <- RunUMAP(tnk, dims = tnk_use_dims, verbose = FALSE)
tnk$tnk_subcluster <- as.character(Idents(tnk))

tnk$t_score <- module_score(tnk, t_markers)
tnk$nk_score <- module_score(tnk, nk_markers)
tnk$t_minus_nk <- tnk$t_score - tnk$nk_score
tnk$nk_minus_t <- tnk$nk_score - tnk$t_score
tnk$exclude_score <- module_score(tnk, exclude_markers)
tnk$core_t_score <- module_score(tnk, core_t_markers)
tnk$core_nk_score <- module_score(tnk, core_nk_markers)

sub_ids <- sort(unique(tnk$tnk_subcluster))
sub_scores <- do.call(rbind, lapply(sub_ids, function(sc) {
  idx <- tnk$tnk_subcluster == sc
  data.frame(
    tnk_subcluster = sc,
    cells = sum(idx),
    mean_t_score = mean(tnk$t_score[idx], na.rm = TRUE),
    mean_nk_score = mean(tnk$nk_score[idx], na.rm = TRUE),
    mean_t_minus_nk = mean(tnk$t_minus_nk[idx], na.rm = TRUE),
    mean_nk_minus_t = mean(tnk$nk_minus_t[idx], na.rm = TRUE),
    mean_core_t = mean(tnk$core_t_score[idx], na.rm = TRUE),
    mean_core_nk = mean(tnk$core_nk_score[idx], na.rm = TRUE),
    mean_exclude = mean(tnk$exclude_score[idx], na.rm = TRUE),
    mean_CD3D = gene_mean(tnk, "CD3D", colnames(tnk)[idx]),
    mean_CD3E = gene_mean(tnk, "CD3E", colnames(tnk)[idx]),
    mean_TRAC = gene_mean(tnk, "TRAC", colnames(tnk)[idx]),
    mean_NKG7 = gene_mean(tnk, "NKG7", colnames(tnk)[idx]),
    mean_GNLY = gene_mean(tnk, "GNLY", colnames(tnk)[idx]),
    mean_PRF1 = gene_mean(tnk, "PRF1", colnames(tnk)[idx]),
    mean_KLRD1 = gene_mean(tnk, "KLRD1", colnames(tnk)[idx]),
    stringsAsFactors = FALSE
  )
}))

# Global thresholds from subcluster means (conservative)
med_nk <- stats::median(sub_scores$mean_nk_score, na.rm = TRUE)
med_t <- stats::median(sub_scores$mean_t_score, na.rm = TRUE)
med_ex <- stats::median(sub_scores$mean_exclude, na.rm = TRUE)
q75_ex <- as.numeric(stats::quantile(sub_scores$mean_exclude, 0.75, na.rm = TRUE))

assign_subcluster_label <- function(row) {
  # Contaminant if exclusion markers dominate and lymphoid scores weak
  if (is.finite(row$mean_exclude) && row$mean_exclude >= q75_ex &&
      row$mean_core_nk < med_nk && row$mean_core_t < med_t) {
    return("Other_or_contaminant")
  }
  nk_high <- row$mean_nk_score >= med_nk && row$mean_core_nk >= med_nk * 0.9
  t_high <- row$mean_t_score >= med_t && row$mean_core_t >= med_t * 0.9
  nk_pref <- row$mean_nk_minus_t >= score_delta
  t_pref <- row$mean_t_minus_nk >= score_delta
  # High-confidence NK: NK markers high, TCR/CD3 axis not co-high
  if (nk_high && nk_pref && row$mean_core_t < row$mean_core_nk &&
      row$mean_core_t <= med_t) {
    return("NK_cell")
  }
  # Clear T
  if (t_high && t_pref && row$mean_core_t > row$mean_core_nk) {
    return("T_cell")
  }
  # Both high or unclear
  if (nk_high && t_high) return("T_NK_ambiguous")
  if (nk_high && !t_pref) return("T_NK_ambiguous")
  if (t_high) return("T_cell")
  "T_NK_ambiguous"
}

sub_scores$subcluster_refined_label <- vapply(seq_len(nrow(sub_scores)), function(i) {
  assign_subcluster_label(sub_scores[i, ])
}, character(1))
write.csv(sub_scores, file.path(stage1_dir, "tnk_subcluster_marker_scores.csv"), row.names = FALSE)

# Cell-level refined labels: start from subcluster label, then cell-level veto
sub_map <- setNames(sub_scores$subcluster_refined_label, sub_scores$tnk_subcluster)
tnk$refined_celltype <- unname(sub_map[tnk$tnk_subcluster])

# Cell-level corrections
cell_nk <- tnk$nk_score
cell_t <- tnk$t_score
cell_core_nk <- tnk$core_nk_score
cell_core_t <- tnk$core_t_score
cell_ex <- tnk$exclude_score
med_cell_nk <- stats::median(cell_nk, na.rm = TRUE)
med_cell_t <- stats::median(cell_t, na.rm = TRUE)
q90_ex_cell <- as.numeric(stats::quantile(cell_ex, 0.90, na.rm = TRUE))

refined <- tnk$refined_celltype
# Promote only clear NK cells if subcluster was ambiguous but cell is NK-like
promote_nk <- refined == "T_NK_ambiguous" &
  cell_core_nk > cell_core_t + score_delta &
  cell_nk >= med_cell_nk &
  cell_core_t <= med_cell_t
# Demote NK if cell has high CD3/TCR
demote_nk <- refined == "NK_cell" & (
  cell_core_t > cell_core_nk |
    (cell_core_t >= med_cell_t & cell_core_nk < cell_core_t + score_delta)
)
# Contaminants
contam <- cell_ex >= q90_ex_cell & cell_core_nk < med_cell_nk & cell_core_t < med_cell_t

refined[promote_nk] <- "NK_cell"
refined[demote_nk] <- "T_NK_ambiguous"
refined[contam] <- "Other_or_contaminant"
# Final high-confidence NK: require core NK > core T and NK score above median
final_nk_ok <- refined == "NK_cell" & cell_core_nk > cell_core_t & cell_nk >= med_cell_nk
refined[refined == "NK_cell" & !final_nk_ok] <- "T_NK_ambiguous"
tnk$refined_celltype <- refined

# analysis_celltype inside subset: only high-confidence NK as NK_cell
tnk$analysis_celltype <- tnk$refined_celltype
tnk$analysis_celltype[tnk$refined_celltype != "NK_cell" & tnk$refined_celltype != "T_cell"] <-
  ifelse(tnk$refined_celltype[tnk$refined_celltype != "NK_cell" & tnk$refined_celltype != "T_cell"] ==
           "Other_or_contaminant", "Other_or_contaminant", "T_NK_ambiguous")

# Transfer to full object
obj$tnk_subcluster <- NA_character_
obj$t_score <- NA_real_
obj$nk_score <- NA_real_
obj$t_minus_nk <- NA_real_
obj$nk_minus_t <- NA_real_
obj$refined_celltype <- "Not_in_TNK_candidate"
obj$analysis_celltype <- obj$legacy_manual_celltype

common <- intersect(colnames(obj), colnames(tnk))
obj$tnk_subcluster[common] <- tnk$tnk_subcluster[common]
obj$t_score[common] <- tnk$t_score[common]
obj$nk_score[common] <- tnk$nk_score[common]
obj$t_minus_nk[common] <- tnk$t_minus_nk[common]
obj$nk_minus_t[common] <- tnk$nk_minus_t[common]
obj$refined_celltype[common] <- tnk$refined_celltype[common]

# analysis_celltype on full object:
# - candidate cells: refined T/NK/ambiguous/contaminant
# - non-candidate: legacy manual (but force legacy NK clusters not in refined NK -> T_cell or keep non-NK)
obj$analysis_celltype <- as.character(obj$legacy_manual_celltype)
obj$analysis_celltype[common] <- tnk$refined_celltype[common]
# Never keep legacy NK label for cells that failed high-confidence NK
legacy_nk_mask <- as.character(obj$seurat_clusters) %in% old_nk_clusters
not_new_nk <- obj$analysis_celltype != "NK_cell"
# If was legacy "NK" via cluster but not high-conf, set to refined or T_cell/ambiguous
obj$analysis_celltype[legacy_nk_mask & not_new_nk & obj$refined_celltype == "Not_in_TNK_candidate"] <- "T_NK_ambiguous"
# Epithelial and others from legacy stay
obj$analysis_celltype[obj$legacy_manual_celltype == tumor_epithelial_celltype_label &
                        obj$refined_celltype == "Not_in_TNK_candidate"] <- tumor_epithelial_celltype_label

# Cell-level annotation export
cell_annot <- data.frame(
  cell = colnames(obj),
  sample = obj$sample,
  seurat_clusters = as.character(obj$seurat_clusters),
  legacy_broad_celltype = obj$legacy_broad_celltype,
  legacy_manual_celltype = obj$legacy_manual_celltype,
  tnk_subcluster = obj$tnk_subcluster,
  t_score = obj$t_score,
  nk_score = obj$nk_score,
  t_minus_nk = obj$t_minus_nk,
  nk_minus_t = obj$nk_minus_t,
  refined_celltype = obj$refined_celltype,
  analysis_celltype = obj$analysis_celltype,
  stringsAsFactors = FALSE
)
write.csv(cell_annot, file.path(stage1_dir, "refined_celltype_annotation_by_cell.csv"), row.names = FALSE)

# By original cluster summary
by_cluster <- do.call(rbind, lapply(cluster_ids, function(cl) {
  idx <- as.character(obj$seurat_clusters) == cl
  tab <- table(factor(obj$refined_celltype[idx], levels = c(
    "NK_cell", "T_cell", "T_NK_ambiguous", "Other_or_contaminant", "Not_in_TNK_candidate"
  )))
  data.frame(
    seurat_clusters = cl,
    cells = sum(idx),
    n_NK_cell = as.integer(tab["NK_cell"]),
    n_T_cell = as.integer(tab["T_cell"]),
    n_T_NK_ambiguous = as.integer(tab["T_NK_ambiguous"]),
    n_Other_or_contaminant = as.integer(tab["Other_or_contaminant"]),
    n_Not_in_TNK_candidate = as.integer(tab["Not_in_TNK_candidate"]),
    legacy_broad = cluster_to_broad[cl],
    legacy_manual = cluster_to_manual[cl],
    stringsAsFactors = FALSE
  )
}))
write.csv(by_cluster, file.path(stage1_dir, "refined_celltype_annotation_by_cluster.csv"), row.names = FALSE)

# Subcluster counts by sample
sub_sample <- as.data.frame(table(
  sample = tnk$sample,
  tnk_subcluster = tnk$tnk_subcluster,
  refined_celltype = tnk$refined_celltype
), stringsAsFactors = FALSE)
colnames(sub_sample)[4] <- "cells"
write.csv(sub_sample, file.path(stage1_dir, "tnk_subcluster_counts_by_sample.csv"), row.names = FALSE)

# Plots
write_plot(
  file.path(stage1_dir, "umap_tnk_subclusters.pdf"),
  DimPlot(tnk, group.by = "tnk_subcluster", label = TRUE) + ggtitle("T/NK candidate subclusters"),
  width = 8, height = 6
)
write_plot(
  file.path(stage1_dir, "umap_refined_celltype.pdf"),
  DimPlot(tnk, group.by = "refined_celltype", label = TRUE) + ggtitle("Refined T/NK labels"),
  width = 8, height = 6
)

plot_genes <- available_genes(tnk, c(core_t_markers, core_nk_markers, "FCGR3A", "NCAM1", "GZMB"))
if (length(plot_genes) > 0) {
  write_plot(
    file.path(stage1_dir, "dotplot_t_vs_nk_markers.pdf"),
    DotPlot(tnk, features = plot_genes, group.by = "refined_celltype") +
      RotatedAxis() + ggtitle("T vs NK markers by refined label"),
    width = 10, height = 5
  )
  write_plot(
    file.path(stage1_dir, "featureplot_core_t_nk_markers.pdf"),
    FeaturePlot(tnk, features = intersect(plot_genes, c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "PRF1", "KLRD1")),
                ncol = 3),
    width = 12, height = 10
  )
  write_plot(
    file.path(stage1_dir, "violin_core_t_nk_by_refined.pdf"),
    VlnPlot(tnk, features = intersect(plot_genes, c("CD3D", "TRAC", "NKG7", "GNLY", "PRF1", "KLRD1")),
            group.by = "refined_celltype", pt.size = 0, ncol = 3),
    width = 12, height = 8
  )
}

# QC summary / acceptance
high_nk_cells <- colnames(obj)[obj$analysis_celltype == "NK_cell"]
high_t_cells <- colnames(obj)[obj$analysis_celltype == "T_cell"]
# For T comparison use refined T in candidates + legacy T outside
t_compare_cells <- colnames(obj)[obj$refined_celltype == "T_cell" |
                                   (obj$analysis_celltype == "T_cell" & obj$refined_celltype == "Not_in_TNK_candidate")]
if (length(t_compare_cells) < 20) t_compare_cells <- high_t_cells

split_cluster <- function(cl) {
  idx <- as.character(obj$seurat_clusters) == cl
  tab <- table(factor(obj$refined_celltype[idx], levels = c(
    "NK_cell", "T_cell", "T_NK_ambiguous", "Other_or_contaminant", "Not_in_TNK_candidate"
  )))
  as.list(as.integer(tab))
}
c0 <- split_cluster("0")
c1 <- split_cluster("1")
names(c0) <- c("NK_cell", "T_cell", "T_NK_ambiguous", "Other_or_contaminant", "Not_in_TNK_candidate")
names(c1) <- c("NK_cell", "T_cell", "T_NK_ambiguous", "Other_or_contaminant", "Not_in_TNK_candidate")

nk_cd3 <- gene_mean(obj, core_t_markers, high_nk_cells)
t_cd3 <- gene_mean(obj, core_t_markers, t_compare_cells)
nk_cyto <- gene_mean(obj, c("KLRD1", "GNLY", "NKG7", "PRF1"), high_nk_cells)
t_cyto <- gene_mean(obj, c("KLRD1", "GNLY", "NKG7", "PRF1"), t_compare_cells)

nk_by_sample <- do.call(rbind, lapply(sort(unique(obj$sample)), function(s) {
  idx <- obj$sample == s
  n_nk <- sum(idx & obj$analysis_celltype == "NK_cell")
  data.frame(
    sample = s,
    total_cells = sum(idx),
    high_confidence_nk = n_nk,
    pct_high_confidence_nk = n_nk / sum(idx),
    stringsAsFactors = FALSE
  )
}))

stage1_pass_cd3 <- is.finite(nk_cd3) && is.finite(t_cd3) && nk_cd3 < t_cd3
stage1_pass_nkmark <- is.finite(nk_cyto) && is.finite(t_cyto) && nk_cyto > t_cyto
stage1_pass_n <- length(high_nk_cells) >= 50
stage1_pass <- stage1_pass_cd3 && stage1_pass_nkmark && stage1_pass_n

qc_rows <- list(
  c("stage1_pass", as.character(stage1_pass)),
  c("stage1_pass_cd3_lower_in_nk_than_t", as.character(stage1_pass_cd3)),
  c("stage1_pass_nk_markers_higher_in_nk_than_t", as.character(stage1_pass_nkmark)),
  c("stage1_pass_min_nk_n", as.character(stage1_pass_n)),
  c("candidate_clusters", paste(candidate_clusters, collapse = ";")),
  c("n_candidate_cells", as.character(length(candidate_cells))),
  c("n_tnk_subclusters", as.character(length(sub_ids))),
  c("cluster0_n_NK", as.character(c0[["NK_cell"]])),
  c("cluster0_n_T", as.character(c0[["T_cell"]])),
  c("cluster0_n_ambiguous", as.character(c0[["T_NK_ambiguous"]])),
  c("cluster0_n_contaminant", as.character(c0[["Other_or_contaminant"]])),
  c("cluster1_n_NK", as.character(c1[["NK_cell"]])),
  c("cluster1_n_T", as.character(c1[["T_cell"]])),
  c("cluster1_n_ambiguous", as.character(c1[["T_NK_ambiguous"]])),
  c("cluster1_n_contaminant", as.character(c1[["Other_or_contaminant"]])),
  c("high_confidence_nk_total", as.character(length(high_nk_cells))),
  c("high_confidence_t_total", as.character(length(t_compare_cells))),
  c("legacy_nk_total", as.character(length(old_nk_cells))),
  c("mean_CD3_TRAC_axis_in_high_conf_NK", as.character(nk_cd3)),
  c("mean_CD3_TRAC_axis_in_T", as.character(t_cd3)),
  c("mean_KLRD1_GNLY_NKG7_PRF1_in_high_conf_NK", as.character(nk_cyto)),
  c("mean_KLRD1_GNLY_NKG7_PRF1_in_T", as.character(t_cyto)),
  c("score_delta", as.character(score_delta)),
  c("nk_score_cluster_quantile", as.character(nk_score_cluster_quantile)),
  c("note", "Only analysis_celltype==NK_cell used as high-confidence NK downstream")
)
qc_df <- as.data.frame(do.call(rbind, qc_rows), stringsAsFactors = FALSE)
colnames(qc_df) <- c("item", "value")
# append per-sample NK
for (i in seq_len(nrow(nk_by_sample))) {
  qc_df <- rbind(
    qc_df,
    data.frame(
      item = paste0("sample_", nk_by_sample$sample[i], "_high_conf_nk"),
      value = paste0(nk_by_sample$high_confidence_nk[i], " (",
                     signif(nk_by_sample$pct_high_confidence_nk[i], 3), ")"),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(qc_df, file.path(stage1_dir, "tnk_refinement_qc_summary.csv"), row.names = FALSE)
write.csv(nk_by_sample, file.path(stage1_dir, "high_confidence_nk_by_sample.csv"), row.names = FALSE)

saveRDS(obj, file.path(stage1_dir, "refined_merged_object.rds"))
# Also keep tnk subset for debugging
saveRDS(tnk, file.path(stage1_dir, "tnk_candidate_subcluster_object.rds"))

message("Stage 1 pass = ", stage1_pass)
message("High-conf NK n = ", length(high_nk_cells), " ; legacy NK n = ", length(old_nk_cells))
if (!stage1_pass) {
  message("STAGE 1 FAILED QC. Stopping before Stage 2/3.")
  inv <- list_output_files(stage1_dir)
  write.csv(inv, file.path(stage1_dir, "output_file_inventory.csv"), row.names = FALSE)
  write_session_info(file.path(stage1_dir, "sessionInfo_end.txt"))
  quit(save = "no", status = 2)
}
} # end stage_start < 2 (full stage 1)

# ============================================================
# STAGE 2: NK downstream with refined labels
# ============================================================
if (stage_start > 2) {
  message("Skipping Stage 2 because stage_start=", stage_start)
} else {
message("=== STAGE 2: refined STK31/NK analysis ===")

# STK31 high/low on epithelial (legacy epithelial label preserved for non-candidates)
stk31_expr <- as.numeric(fetch_gene_matrix(obj, target_gene)[target_gene, ])
obj$stk31_expr <- stk31_expr
# Ensure epithelial present in analysis_celltype
obj$analysis_celltype[obj$legacy_manual_celltype == tumor_epithelial_celltype_label &
                        !(obj$analysis_celltype %in% c("NK_cell", "T_cell", "T_NK_ambiguous"))] <-
  tumor_epithelial_celltype_label

epi_assign <- assign_tumor_epithelial_stk31_group(
  obj, stk31_expr,
  celltype_col = "analysis_celltype",
  celltype_label = tumor_epithelial_celltype_label,
  high_quantile = tumor_epithelial_stk31_high_quantile
)
obj$tumor_epithelial_stk31_group <- epi_assign$group

# Sensitivity definitions
epi_idx <- obj$analysis_celltype == tumor_epithelial_celltype_label
epi_cells <- colnames(obj)[epi_idx]
# 1) counts > 0 on data layer already >0 approx; also try counts layer if available
stk31_counts <- tryCatch({
  as.numeric(fetch_gene_matrix(obj, target_gene, slot = "counts")[target_gene, ])
}, error = function(e) stk31_expr)
def_counts <- epi_idx & stk31_counts > 0
def_quantile <- obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"
# 2 already quantile global epithelial
# 3 per-sample epithelial quantile
def_sample <- rep(FALSE, ncol(obj))
names(def_sample) <- colnames(obj)
for (s in unique(obj$sample)) {
  s_epi <- epi_idx & obj$sample == s
  if (sum(s_epi) < 5) next
  cut_s <- as.numeric(stats::quantile(stk31_expr[s_epi], tumor_epithelial_stk31_high_quantile, na.rm = TRUE))
  if (!is.finite(cut_s) || cut_s == 0) {
    def_sample[s_epi & stk31_expr > 0] <- TRUE
  } else {
    def_sample[s_epi & stk31_expr >= cut_s] <- TRUE
  }
}

sens <- data.frame(
  definition = c("counts_gt_0", "epithelial_global_quantile", "epithelial_per_sample_quantile"),
  n_high = c(sum(def_counts), sum(def_quantile), sum(def_sample)),
  n_epithelial = rep(sum(epi_idx), 3),
  pct_of_epithelial = c(
    mean(def_counts[epi_idx]),
    mean(def_quantile[epi_idx]),
    mean(def_sample[epi_idx])
  ),
  cutoff_note = c(
    "stk31 counts > 0 within all cells then intersect epithelial for n_high count of epi&pos",
    paste0("global epithelial quantile q=", tumor_epithelial_stk31_high_quantile,
           " cutoff=", epi_assign$cutoff),
    paste0("per-sample epithelial quantile q=", tumor_epithelial_stk31_high_quantile)
  ),
  stringsAsFactors = FALSE
)
# fix counts definition to epithelial only high
sens$n_high[1] <- sum(epi_idx & stk31_counts > 0)
sens$pct_of_epithelial[1] <- mean(stk31_counts[epi_idx] > 0)
write.csv(sens, file.path(stage2_dir, "stk31_definition_sensitivity.csv"), row.names = FALSE)

high_cells <- colnames(obj)[obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"]
low_cells <- colnames(obj)[obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"]
nk_cells <- colnames(obj)[obj$analysis_celltype == "NK_cell"]
other_cells <- setdiff(colnames(obj), nk_cells)

analysis_summary <- data.frame(
  item = c(
    "target_gene", "tumor_epithelial_celltype_label", "tumor_epithelial_stk31_high_cutoff",
    "tumor_epithelial_cells", "stk31_high_tumor_epithelial_cells", "stk31_low_tumor_epithelial_cells",
    "high_confidence_nk_cells", "legacy_nk_cells", "samples_analyzed",
    "stk31_definition", "nk_definition"
  ),
  value = c(
    target_gene, tumor_epithelial_celltype_label, as.character(epi_assign$cutoff),
    as.character(epi_assign$tumor_epithelial_cells),
    as.character(epi_assign$high_cells), as.character(epi_assign$low_cells),
    as.character(length(nk_cells)), as.character(length(old_nk_cells)),
    paste(sort(unique(obj$sample)), collapse = ";"),
    paste0("epithelial_global_quantile_", tumor_epithelial_stk31_high_quantile,
           " (fallback counts>0 if cutoff==0)"),
    "analysis_celltype==NK_cell high-confidence refined"
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_summary, file.path(stage2_dir, "analysis_summary.csv"), row.names = FALSE)

# Old vs refined NK overlap
overlap <- data.frame(
  cell = colnames(obj),
  legacy_nk = colnames(obj) %in% old_nk_cells,
  refined_nk = colnames(obj) %in% nk_cells,
  stringsAsFactors = FALSE
)
write.csv(overlap, file.path(stage2_dir, "old_vs_refined_nk_overlap.csv"), row.names = FALSE)
overlap_summary <- data.frame(
  metric = c(
    "legacy_nk_n", "refined_nk_n", "intersection_n",
    "legacy_only_n", "refined_only_n",
    "jaccard", "fraction_legacy_kept_as_refined"
  ),
  value = c(
    length(old_nk_cells), length(nk_cells), length(intersect(old_nk_cells, nk_cells)),
    length(setdiff(old_nk_cells, nk_cells)), length(setdiff(nk_cells, old_nk_cells)),
    length(intersect(old_nk_cells, nk_cells)) /
      length(union(old_nk_cells, nk_cells)),
    length(intersect(old_nk_cells, nk_cells)) / max(1, length(old_nk_cells))
  ),
  stringsAsFactors = FALSE
)
write.csv(overlap_summary, file.path(stage2_dir, "old_vs_refined_nk_summary.csv"), row.names = FALSE)

message("Running DE: NK vs all other")
run_marker_test(obj, nk_cells, other_cells, file.path(stage2_dir, "markers_nk_vs_all_other.csv"))
message("Running DE: STK31-high epithelial vs NK")
run_marker_test(
  obj, high_cells, nk_cells,
  file.path(stage2_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv")
)

message("Candidate L-R")
lr_links <- infer_lr_links(obj, high_cells, nk_cells, lr_reference)
write.csv(lr_links, file.path(stage2_dir, "candidate_ligand_receptor_pathways.csv"), row.names = FALSE)
pathway_summary <- aggregate(interaction_score ~ direction + pathway, data = lr_links, FUN = max)
pathway_summary <- pathway_summary[order(-pathway_summary$interaction_score), ]
write.csv(pathway_summary, file.path(stage2_dir, "candidate_pathway_summary.csv"), row.names = FALSE)

# Sample-level STK31-NK relationship
sample_rel <- do.call(rbind, lapply(sort(unique(as.character(obj$sample))), function(s) {
  idx <- obj$sample == s
  data.frame(
    sample = s,
    total_cells = sum(idx),
    epithelial_cells = sum(idx & obj$analysis_celltype == tumor_epithelial_celltype_label),
    stk31_high_epithelial = sum(idx & obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"),
    pct_epithelial_stk31_high = mean(obj$tumor_epithelial_stk31_group[idx & obj$analysis_celltype == tumor_epithelial_celltype_label] ==
                                       "STK31_high_tumor_epithelial"),
    high_confidence_nk = sum(idx & obj$analysis_celltype == "NK_cell"),
    pct_nk = mean(obj$analysis_celltype[idx] == "NK_cell"),
    stringsAsFactors = FALSE
  )
}))
write.csv(sample_rel, file.path(stage2_dir, "sample_level_stk31_nk_relationship.csv"), row.names = FALSE)

# NK function scores by sample (high-conf NK only)
for (score_name in names(nk_function_sets)) {
  obj[[paste0("score_", score_name)]] <- module_score(obj, nk_function_sets[[score_name]])
}
nk_score_rows <- list()
for (score_name in names(nk_function_sets)) {
  col <- paste0("score_", score_name)
  for (s in sort(unique(as.character(obj$sample)))) {
    idx <- obj$sample == s & obj$analysis_celltype == "NK_cell"
    vals <- obj[[col]][idx, 1]
    nk_score_rows[[length(nk_score_rows) + 1]] <- data.frame(
      sample = s, score_name = score_name,
      cells = sum(idx),
      mean_score = ifelse(length(vals) > 0, mean(vals, na.rm = TRUE), NA_real_),
      median_score = ifelse(length(vals) > 0, median(vals, na.rm = TRUE), NA_real_),
      pct_score_positive = ifelse(length(vals) > 0, mean(vals > 0, na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  }
}
nk_scores_by_sample <- do.call(rbind, nk_score_rows)
write.csv(nk_scores_by_sample, file.path(stage2_dir, "nk_function_scores_by_sample.csv"), row.names = FALSE)

# Plots
obj$relationship_group <- ifelse(
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial",
  ifelse(
    obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial", "STK31_low_tumor_epithelial",
    ifelse(obj$analysis_celltype == "NK_cell", "NK_cell", "Other")
  )
)
if ("umap" %in% names(obj@reductions) || "UMAP" %in% names(obj@reductions)) {
  write_plot(
    file.path(stage2_dir, "umap_analysis_celltype.pdf"),
    DimPlot(obj, group.by = "analysis_celltype", label = TRUE, repel = TRUE) +
      ggtitle("Analysis celltype (refined NK)"),
    width = 10, height = 7
  )
  write_plot(
    file.path(stage2_dir, "umap_relationship_group.pdf"),
    DimPlot(obj, group.by = "relationship_group") +
      ggtitle("STK31 high/low epithelial and refined NK"),
    width = 9, height = 6
  )
}
if (length(nk_cells) > 0) {
  write_plot(
    file.path(stage2_dir, "dotplot_nk_markers_by_analysis_celltype.pdf"),
    DotPlot(obj, features = available_genes(obj, c(core_t_markers, core_nk_markers)),
            group.by = "analysis_celltype") + RotatedAxis(),
    width = 11, height = 6
  )
}
if (nrow(nk_scores_by_sample) > 0) {
  write_plot(
    file.path(stage2_dir, "nk_function_scores_by_sample.pdf"),
    ggplot(nk_scores_by_sample, aes(x = sample, y = mean_score, fill = score_name)) +
      geom_col(position = "dodge") + theme_bw() +
      ggtitle("NK function module scores (high-confidence NK)") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)),
    width = 10, height = 5
  )
}

message("Stage 2 done")
} # end stage 2

# ============================================================
# STAGE 3: CellChat high vs low differential
# ============================================================
message("=== STAGE 3: CellChat high-vs-low ===")
# Ensure STK31 groups exist if resumed into stage 3 only
if (!"tumor_epithelial_stk31_group" %in% colnames(obj@meta.data)) {
  stk31_expr <- as.numeric(fetch_gene_matrix(obj, target_gene)[target_gene, ])
  obj$analysis_celltype[obj$legacy_manual_celltype == tumor_epithelial_celltype_label &
                          !(obj$analysis_celltype %in% c("NK_cell", "T_cell", "T_NK_ambiguous"))] <-
    tumor_epithelial_celltype_label
  epi_assign <- assign_tumor_epithelial_stk31_group(
    obj, stk31_expr,
    celltype_col = "analysis_celltype",
    celltype_label = tumor_epithelial_celltype_label,
    high_quantile = tumor_epithelial_stk31_high_quantile
  )
  obj$tumor_epithelial_stk31_group <- epi_assign$group
}
if (!"analysis_celltype" %in% colnames(obj@meta.data)) {
  stop("analysis_celltype missing; cannot run Stage 3")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  warning("CellChat not installed; writing empty Stage 3 placeholders and exiting successfully after Stage 1/2.")
  write.csv(data.frame(note = "CellChat not installed"), file.path(stage3_dir, "cellchat_skipped.csv"), row.names = FALSE)
} else {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan("sequential")
    options(future.globals.maxSize = 8 * 1024^3)
  }

  # Build cellchat groups
  obj$cellchat_group <- as.character(obj$analysis_celltype)
  obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"] <-
    "STK31_high_tumor_epithelial"
  obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"] <-
    "STK31_low_tumor_epithelial"
  obj$cellchat_group[obj$analysis_celltype == "NK_cell"] <- "NK_cell"
  # Collapse rare ambiguous labels
  obj$cellchat_group[obj$cellchat_group %in% c("T_NK_ambiguous", "Other_or_contaminant", "Unassigned")] <-
    "Other_immune_or_unassigned"

  group_counts <- sort(table(obj$cellchat_group), decreasing = TRUE)
  keep_groups <- names(group_counts[group_counts >= min_cells_per_cellchat_group])
  # Always try to keep key groups if present with enough cells
  key_groups <- c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell")
  for (kg in key_groups) {
    if (!kg %in% keep_groups && sum(obj$cellchat_group == kg) >= 3) {
      message("Warning: key group ", kg, " has only ", sum(obj$cellchat_group == kg), " cells")
    }
  }
  cc_obj <- subset(obj, cells = colnames(obj)[obj$cellchat_group %in% keep_groups])
  cc_obj$cellchat_group <- factor(as.character(cc_obj$cellchat_group))

  before <- as.data.frame(table(cc_obj$cellchat_group), stringsAsFactors = FALSE)
  colnames(before) <- c("cellchat_group", "cells_before_sampling")

  cells_by_group <- split(colnames(cc_obj), cc_obj$cellchat_group)
  keep_cells <- unlist(lapply(cells_by_group, function(cells) {
    if (length(cells) > max_cells_per_cellchat_group) sample(cells, max_cells_per_cellchat_group) else cells
  }), use.names = FALSE)
  cc_obj <- subset(cc_obj, cells = keep_cells)
  cc_obj$cellchat_group <- droplevels(factor(as.character(cc_obj$cellchat_group)))

  after <- as.data.frame(table(cc_obj$cellchat_group), stringsAsFactors = FALSE)
  colnames(after) <- c("cellchat_group", "cells_after_sampling")
  group_counts_out <- merge(before, after, by = "cellchat_group", all = TRUE)
  group_counts_out$cells_after_sampling[is.na(group_counts_out$cells_after_sampling)] <- 0
  write.csv(group_counts_out, file.path(stage3_dir, "cellchat_group_counts.csv"), row.names = FALSE)

  required_for_diff <- c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell")
  if (!all(required_for_diff %in% levels(cc_obj$cellchat_group))) {
    warning("Missing required CellChat groups after filtering: ",
            paste(setdiff(required_for_diff, levels(cc_obj$cellchat_group)), collapse = ", "))
  }

  message("Running CellChat (this may take a while)")
  data_input <- fetch_gene_matrix(cc_obj, rownames(cc_obj), slot = "data")
  meta <- data.frame(labels = cc_obj$cellchat_group, row.names = colnames(cc_obj), stringsAsFactors = FALSE)
  cellchat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_per_cellchat_group)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(stage3_dir, "merged_cellchat_object.rds"))

  communication <- CellChat::subsetCommunication(cellchat)
  write.csv(communication, file.path(stage3_dir, "merged_cellchat_communications.csv"), row.names = FALSE)

  # Helper to extract pair probs
  get_pair_prob <- function(comm, source, target) {
    sub <- comm[comm$source == source & comm$target == target, , drop = FALSE]
    if (nrow(sub) == 0) return(sub)
    sub
  }

  high_to_nk <- get_pair_prob(communication, "STK31_high_tumor_epithelial", "NK_cell")
  low_to_nk <- get_pair_prob(communication, "STK31_low_tumor_epithelial", "NK_cell")
  nk_to_high <- get_pair_prob(communication, "NK_cell", "STK31_high_tumor_epithelial")
  nk_to_low <- get_pair_prob(communication, "NK_cell", "STK31_low_tumor_epithelial")

  # Aggregate L-R level
  make_diff_table <- function(high_df, low_df, direction_label) {
    # key by pathway + ligand + receptor if columns exist
    key_cols <- intersect(c("pathway_name", "ligand", "receptor", "interaction_name"), colnames(communication))
    if (!"pathway_name" %in% colnames(communication)) {
      # CellChat sometimes uses pathway_name
      if ("pathway" %in% colnames(communication)) {
        communication$pathway_name <<- communication$pathway
      }
    }
    # rebuild keys safely
    add_key <- function(df) {
      if (nrow(df) == 0) {
        df$key <- character(0)
        return(df)
      }
      pn <- if ("pathway_name" %in% names(df)) df$pathway_name else if ("pathway" %in% names(df)) df$pathway else "unknown"
      lig <- if ("ligand" %in% names(df)) df$ligand else "NA"
      rec <- if ("receptor" %in% names(df)) df$receptor else "NA"
      df$pathway_name <- pn
      df$ligand <- lig
      df$receptor <- rec
      df$key <- paste(pn, lig, rec, sep = "||")
      df
    }
    high_df <- add_key(high_df)
    low_df <- add_key(low_df)
    keys <- unique(c(high_df$key, low_df$key))
    if (length(keys) == 0) {
      return(data.frame(
        pathway_name = character(0), ligand = character(0), receptor = character(0),
        prob_high_to_nk = numeric(0), prob_low_to_nk = numeric(0),
        delta_high_minus_low = numeric(0), ratio_high_over_low = numeric(0),
        evidence_label = character(0), direction = character(0),
        stringsAsFactors = FALSE
      ))
    }
    do.call(rbind, lapply(keys, function(k) {
      h <- high_df[high_df$key == k, , drop = FALSE]
      l <- low_df[low_df$key == k, , drop = FALSE]
      ph <- if (nrow(h)) sum(h$prob, na.rm = TRUE) else 0
      pl <- if (nrow(l)) sum(l$prob, na.rm = TRUE) else 0
      parts <- strsplit(k, "||", fixed = TRUE)[[1]]
      data.frame(
        pathway_name = parts[1],
        ligand = parts[2],
        receptor = parts[3],
        prob_high = ph,
        prob_low = pl,
        delta_high_minus_low = ph - pl,
        ratio_high_over_low = ifelse(pl > 0, ph / pl, ifelse(ph > 0, Inf, NA_real_)),
        evidence_label = label_evidence(ph, pl),
        direction = direction_label,
        stringsAsFactors = FALSE
      )
    }))
  }

  # Fix column naming for required files
  high_low_to_nk <- make_diff_table(high_to_nk, low_to_nk, "epithelial_to_NK")
  if (nrow(high_low_to_nk) > 0) {
    names(high_low_to_nk)[names(high_low_to_nk) == "prob_high"] <- "prob_high_to_nk"
    names(high_low_to_nk)[names(high_low_to_nk) == "prob_low"] <- "prob_low_to_nk"
  } else {
    high_low_to_nk <- data.frame(
      pathway_name = character(0), ligand = character(0), receptor = character(0),
      prob_high_to_nk = numeric(0), prob_low_to_nk = numeric(0),
      delta_high_minus_low = numeric(0), ratio_high_over_low = numeric(0),
      evidence_label = character(0), direction = character(0)
    )
  }
  write.csv(high_low_to_nk, file.path(stage3_dir, "high_vs_low_to_nk_cellchat_diff.csv"), row.names = FALSE)

  nk_to_high_low <- make_diff_table(nk_to_high, nk_to_low, "NK_to_epithelial")
  if (nrow(nk_to_high_low) > 0) {
    names(nk_to_high_low)[names(nk_to_high_low) == "prob_high"] <- "prob_nk_to_high"
    names(nk_to_high_low)[names(nk_to_high_low) == "prob_low"] <- "prob_nk_to_low"
  }
  write.csv(nk_to_high_low, file.path(stage3_dir, "nk_to_high_vs_low_cellchat_diff.csv"), row.names = FALSE)

  # Mechanism axis comparison
  axis_rows <- list()
  for (axis_name in names(mechanism_axes)) {
    genes <- mechanism_axes[[axis_name]]
    pick_axis <- function(df, prob_col1, prob_col2) {
      if (nrow(df) == 0) return(c(0, 0))
      hit <- grepl(paste(genes, collapse = "|"), paste(df$ligand, df$receptor, df$pathway_name), ignore.case = TRUE)
      # also MHC pathway name heuristics
      if (axis_name == "MHC-I / HLA axis") {
        hit <- hit | grepl("MHC-I|HLA|KIR", df$pathway_name, ignore.case = TRUE) |
          grepl("^HLA-", df$ligand) | grepl("^HLA-", df$receptor)
      }
      if (axis_name == "TIGIT / NECTIN axis") {
        hit <- hit | grepl("TIGIT|NECTIN|PVR|CD226", df$pathway_name, ignore.case = TRUE)
      }
      if (axis_name == "TGFb axis") {
        hit <- hit | grepl("TGF", df$pathway_name, ignore.case = TRUE)
      }
      if (axis_name == "IFNG / IFN axis") {
        hit <- hit | grepl("IFN|INF", df$pathway_name, ignore.case = TRUE)
      }
      if (axis_name == "NKG2D axis") {
        hit <- hit | grepl("NKG2D|MICA|MICB|ULBP", df$pathway_name, ignore.case = TRUE)
      }
      c(sum(df[[prob_col1]][hit], na.rm = TRUE), sum(df[[prob_col2]][hit], na.rm = TRUE))
    }
    if (nrow(high_low_to_nk) > 0) {
      probs <- pick_axis(high_low_to_nk, "prob_high_to_nk", "prob_low_to_nk")
    } else {
      probs <- c(0, 0)
    }
    axis_rows[[length(axis_rows) + 1]] <- data.frame(
      mechanism_axis = axis_name,
      direction = "epithelial_to_NK",
      total_prob_high = probs[1],
      total_prob_low = probs[2],
      delta_high_minus_low = probs[1] - probs[2],
      evidence_label = label_evidence(probs[1], probs[2]),
      stringsAsFactors = FALSE
    )
  }
  axis_df <- do.call(rbind, axis_rows)
  write.csv(axis_df, file.path(stage3_dir, "mechanism_axis_high_low_comparison.csv"), row.names = FALSE)

  # MHC-I special check
  if (nrow(high_low_to_nk) > 0) {
    mhc_hit <- grepl("MHC-I|HLA|KIR", high_low_to_nk$pathway_name, ignore.case = TRUE) |
      grepl("^HLA-", high_low_to_nk$ligand) |
      high_low_to_nk$ligand %in% mhc_i_genes
    mhc_sub <- high_low_to_nk[mhc_hit, , drop = FALSE]
    prob_high_mhc <- sum(mhc_sub$prob_high_to_nk, na.rm = TRUE)
    prob_low_mhc <- sum(mhc_sub$prob_low_to_nk, na.rm = TRUE)
  } else {
    mhc_sub <- high_low_to_nk
    prob_high_mhc <- 0
    prob_low_mhc <- 0
  }
  if (prob_high_mhc <= prob_low_mhc) {
    mhc_conclusion <- "MHC-I/KIR is not higher in STK31-high epithelial to NK compared with STK31-low epithelial to NK in the current CellChat output."
  } else {
    mhc_conclusion <- "MHC-I/KIR total probability is higher in STK31-high epithelial to NK than STK31-low in the current CellChat output; still observational/candidate only, not causal."
  }
  mhc_check <- data.frame(
    item = c(
      "prob_high_to_nk_mhc_i", "prob_low_to_nk_mhc_i", "delta_high_minus_low",
      "n_interactions_high_or_low", "conclusion"
    ),
    value = c(
      as.character(prob_high_mhc), as.character(prob_low_mhc),
      as.character(prob_high_mhc - prob_low_mhc),
      as.character(nrow(mhc_sub)), mhc_conclusion
    ),
    stringsAsFactors = FALSE
  )
  write.csv(mhc_check, file.path(stage3_dir, "mhc_i_high_low_check.csv"), row.names = FALSE)
  write.csv(mhc_sub, file.path(stage3_dir, "mhc_i_high_low_interactions_detail.csv"), row.names = FALSE)

  # Plots
  if (nrow(high_low_to_nk) > 0) {
    plot_df <- high_low_to_nk
    plot_df$pair <- paste(plot_df$ligand, plot_df$receptor, sep = "->")
    plot_df <- plot_df[order(-pmax(plot_df$prob_high_to_nk, plot_df$prob_low_to_nk)), ]
    plot_df <- head(plot_df, 40)
    long <- rbind(
      data.frame(pair = plot_df$pair, group = "STK31_high->NK", prob = plot_df$prob_high_to_nk,
                 pathway = plot_df$pathway_name, stringsAsFactors = FALSE),
      data.frame(pair = plot_df$pair, group = "STK31_low->NK", prob = plot_df$prob_low_to_nk,
                 pathway = plot_df$pathway_name, stringsAsFactors = FALSE)
    )
    long$pair <- factor(long$pair, levels = rev(unique(plot_df$pair)))
    write_plot(
      file.path(stage3_dir, "cellchat_high_low_nk_bubble.pdf"),
      ggplot(long, aes(x = group, y = pair, size = prob, color = pathway)) +
        geom_point(alpha = 0.85) + theme_bw() +
        labs(title = "CellChat: STK31-high vs low epithelial to refined NK",
             x = NULL, y = "Ligand->Receptor", size = "Prob") +
        theme(axis.text.x = element_text(angle = 20, hjust = 1)),
      width = 10, height = max(6, 0.22 * length(unique(long$pair)) + 2)
    )
  } else {
    write_plot(
      file.path(stage3_dir, "cellchat_high_low_nk_bubble.pdf"),
      ggplot() + theme_void() + ggtitle("No high/low->NK interactions detected"),
      width = 6, height = 4
    )
  }

  write_plot(
    file.path(stage3_dir, "cellchat_mechanism_axis_high_low_barplot.pdf"),
    {
      ad <- axis_df
      ad_long <- rbind(
        data.frame(axis = ad$mechanism_axis, group = "high->NK", prob = ad$total_prob_high, stringsAsFactors = FALSE),
        data.frame(axis = ad$mechanism_axis, group = "low->NK", prob = ad$total_prob_low, stringsAsFactors = FALSE)
      )
      ggplot(ad_long, aes(x = axis, y = prob, fill = group)) +
        geom_col(position = "dodge") + theme_bw() +
        labs(title = "Mechanism axis: high vs low epithelial to NK", x = NULL, y = "Total CellChat probability") +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    },
    width = 10, height = 5
  )

  # Heatmap of top pairs high vs low
  if (nrow(high_low_to_nk) > 0 && requireNamespace("pheatmap", quietly = TRUE)) {
    top <- head(high_low_to_nk[order(-pmax(high_low_to_nk$prob_high_to_nk, high_low_to_nk$prob_low_to_nk)), ], 30)
    mat <- as.matrix(top[, c("prob_high_to_nk", "prob_low_to_nk")])
    rownames(mat) <- paste(top$ligand, top$receptor, sep = "->")
    pdf(file.path(stage3_dir, "cellchat_high_low_nk_heatmap.pdf"), width = 6, height = max(5, 0.25 * nrow(mat) + 2))
    pheatmap::pheatmap(
      mat, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA,
      main = "CellChat prob: high vs low epithelial -> NK",
      color = colorRampPalette(c("white", "#fee08b", "#d73027"))(50)
    )
    dev.off()
  } else if (nrow(high_low_to_nk) > 0) {
    top <- head(high_low_to_nk[order(-pmax(high_low_to_nk$prob_high_to_nk, high_low_to_nk$prob_low_to_nk)), ], 30)
    top$pair <- paste(top$ligand, top$receptor, sep = "->")
    long <- rbind(
      data.frame(pair = top$pair, group = "high", prob = top$prob_high_to_nk, stringsAsFactors = FALSE),
      data.frame(pair = top$pair, group = "low", prob = top$prob_low_to_nk, stringsAsFactors = FALSE)
    )
    write_plot(
      file.path(stage3_dir, "cellchat_high_low_nk_heatmap.pdf"),
      ggplot(long, aes(x = group, y = reorder(pair, prob), fill = prob)) +
        geom_tile() + scale_fill_gradient(low = "white", high = "#d73027") +
        theme_bw() + labs(title = "CellChat high vs low -> NK", y = NULL),
      width = 6, height = max(5, 0.25 * length(unique(long$pair)) + 2)
    )
  } else {
    write_plot(
      file.path(stage3_dir, "cellchat_high_low_nk_heatmap.pdf"),
      ggplot() + theme_void() + ggtitle("No interactions for heatmap"),
      width = 6, height = 4
    )
  }

  # Optional circle plots
  if (!is.null(cellchat@net$count)) {
    group_size <- as.numeric(table(cellchat@idents))
    write_base_pdf(
      file.path(stage3_dir, "cellchat_circle_count.pdf"),
      CellChat::netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE,
                                 label.edge = FALSE, title.name = "Number of interactions"),
      width = 8, height = 8
    )
  }

  message("MHC-I conclusion: ", mhc_conclusion)
  message("Stage 3 done")
}

# -------------------- inventory --------------------
inv <- list_output_files(c(stage1_dir, stage2_dir, stage3_dir))
write.csv(inv, file.path(stage1_dir, "output_file_inventory.csv"), row.names = FALSE)
write_session_info(file.path(stage1_dir, "sessionInfo_end.txt"))
message("=== ALL DONE ===")
message("Stage1: ", stage1_dir)
message("Stage2: ", stage2_dir)
message("Stage3: ", stage3_dir)
message("End: ", Sys.time())
