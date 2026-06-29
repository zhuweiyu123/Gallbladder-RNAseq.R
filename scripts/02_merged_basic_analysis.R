suppressPackageStartupMessages({
  # 加载单细胞分析和绘图需要的 R 包；suppressPackageStartupMessages 用来隐藏启动提示。
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# 固定随机种子，让 PCA/UMAP/聚类等带随机性的步骤尽量可重复。
set.seed(20260624)

# 服务器路径：这个脚本设计为在 aiserver 上运行，不是在本地 Windows 上运行。
# base_dir 是 10x 原始矩阵所在目录；out_dir 是四个样本合并分析的结果输出目录。
base_dir <- "/home/zhuweiyu/codex-r/gallbladder_cancer"
out_dir <- "/home/zhuweiyu/codex-r/results/merged_basic_seurat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 本脚本会直接合并分析四个样本。
samples <- c("tissue1", "tissue2", "tissue4", "tissue5")

# 第一轮质控参数，先用保守阈值跑通流程；看完 QC 图后可以再调整。
# min_cells_per_gene：基因至少在多少细胞中出现，低于该值的基因会被过滤。
# min/max_features_per_cell：每个细胞检测到的基因数下限/上限，用来去掉低质量细胞和疑似 doublet。
# max_percent_mt：线粒体基因比例上限，过高通常提示细胞状态较差。
min_cells_per_gene <- 3
min_features_per_cell <- 200
max_features_per_cell <- 6000
max_percent_mt <- 20

# 降维和聚类参数。
# n_variable_features：挑选多少个高变基因用于 PCA。
# n_pcs：计算主成分的总数。
# use_dims：实际用于聚类和 UMAP 的前多少个主成分，后面 PCs 噪音大，取 1:20 更稳。
# cluster_resolution：聚类分辨率，数值越高通常分出更多 cluster。
n_variable_features <- 2000
n_pcs <- 30
use_dims <- 1:20
cluster_resolution <- 0.5

read_sample <- function(sample_id) {
  # 每个样本的 10x filtered_feature_bc_matrix 目录。
  data_dir <- file.path(base_dir, sample_id, "filtered_feature_bc_matrix")
  if (!dir.exists(data_dir)) {
    stop("Missing 10x directory: ", data_dir)
  }

  # 读取 10x 数据。部分 10x 输出会返回 list，这里优先取 Gene Expression 矩阵。
  counts <- Read10X(data.dir = data_dir)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  # 创建 Seurat 对象，并记录样本名，方便后面按 sample 分组画图。
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells_per_gene,
    min.features = min_features_per_cell
  )
  obj$sample <- sample_id

  # 计算每个细胞线粒体基因比例。人类基因名通常以 MT- 开头。
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

qc_table <- function(object_list) {
  # 汇总每个样本的细胞数、基因数、UMI 数和线粒体比例，输出成 CSV 方便检查。
  do.call(rbind, lapply(object_list, function(obj) {
    data.frame(
      sample = unique(obj$sample),
      cells = ncol(obj),
      genes = nrow(obj),
      median_nFeature_RNA = median(obj$nFeature_RNA),
      median_nCount_RNA = median(obj$nCount_RNA),
      median_percent_mt = median(obj$percent.mt)
    )
  }))
}

message("Reading samples from: ", base_dir)
objects <- lapply(samples, read_sample)
names(objects) <- samples

# 保存过滤前的 QC 概况，作为后面判断过滤是否过严/过松的基准。
write.csv(
  qc_table(objects),
  file.path(out_dir, "qc_before_filter.csv"),
  row.names = FALSE
)

# 按前面设定的阈值过滤细胞。
filtered <- lapply(objects, function(obj) {
  subset(
    obj,
    subset = nFeature_RNA >= min_features_per_cell &
      nFeature_RNA <= max_features_per_cell &
      percent.mt <= max_percent_mt
  )
})

# 保存过滤后的 QC 概况，用来和过滤前对比。
write.csv(
  qc_table(filtered),
  file.path(out_dir, "qc_after_filter.csv"),
  row.names = FALSE
)

message("Merging samples")

# 把四个样本合并成一个 Seurat 对象；add.cell.ids 会给细胞条形码加样本前缀，避免重名。
combined <- merge(
  x = filtered[[1]],
  y = filtered[-1],
  add.cell.ids = samples,
  project = "gallbladder_cancer"
)

# 画 QC 小提琴图：基因数、UMI 数、线粒体比例，并按样本分组。
pdf(file.path(out_dir, "qc_violin.pdf"), width = 12, height = 6)
print(VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  ncol = 3,
  pt.size = 0.01
))
dev.off()

message("Running normalization and dimensional reduction")

# NormalizeData：对每个细胞做标准化，降低测序深度差异的影响。
combined <- NormalizeData(combined)

# FindVariableFeatures：找高变基因，后续 PCA 主要基于这些信息量更高的基因。
combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = n_variable_features
)

# ScaleData：对高变基因做中心化/标准化。第一轮只 scale 高变基因，可以明显降低内存占用。
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)

# RunPCA：把高维表达矩阵压缩成主成分，便于后续建邻居图和可视化。
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = n_pcs, verbose = FALSE)

# FindNeighbors + FindClusters：根据 PCA 空间中的相似性构建细胞图，并进行聚类。
combined <- FindNeighbors(combined, dims = use_dims)
combined <- FindClusters(combined, resolution = cluster_resolution)

# RunUMAP：把细胞降到二维，方便观察 cluster 和样本分布。
combined <- RunUMAP(combined, dims = use_dims)

# 按样本来源上色的 UMAP，用来观察不同样本是否混合、是否有明显批次差异。
pdf(file.path(out_dir, "umap_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(combined, reduction = "umap", group.by = "sample"))
dev.off()

# 按聚类结果上色的 UMAP，并在图上标出 cluster 编号。
pdf(file.path(out_dir, "umap_by_cluster.pdf"), width = 8, height = 6)
print(DimPlot(combined, reduction = "umap", label = TRUE))
dev.off()

# 统计每个 cluster 中来自各样本的细胞数，帮助判断 cluster 是否由某个样本主导。
cluster_counts <- as.data.frame.matrix(table(combined$seurat_clusters, combined$sample))
cluster_counts$cluster <- rownames(cluster_counts)
cluster_counts <- cluster_counts[, c("cluster", samples)]
write.csv(
  cluster_counts,
  file.path(out_dir, "cluster_counts_by_sample.csv"),
  row.names = FALSE
)

# 保存本次分析的关键参数和总览数字，方便以后回看。
analysis_summary <- data.frame(
  metric = c(
    "total_cells_after_filter",
    "total_genes",
    "clusters",
    "variable_features",
    "pcs",
    "use_dims",
    "cluster_resolution"
  ),
  value = c(
    ncol(combined),
    nrow(combined),
    length(levels(combined$seurat_clusters)),
    n_variable_features,
    n_pcs,
    paste(use_dims, collapse = "-"),
    cluster_resolution
  )
)
write.csv(
  analysis_summary,
  file.path(out_dir, "analysis_summary.csv"),
  row.names = FALSE
)

# Seurat v5 合并多个样本后会保留 data.tissue* 等多 layer；后续表达汇总和 marker 分析需要先合并成单层。
if (packageVersion("SeuratObject") >= "5.0.0") {
  combined <- JoinLayers(combined, assay = "RNA")
}

# 保存完整 Seurat 对象，后续可以直接读取它继续做 marker、注释、差异分析等。
saveRDS(combined, file.path(out_dir, "gallbladder_cancer_merged_basic_seurat.rds"))

message("Basic analysis done. Outputs written to: ", out_dir)

# ============================================================
# STK31-NK 完整分析（在合并数据上）
# 用途：胆囊癌合并数据（tissue1/2/4/5），聚焦 STK31
#   [1] 定位 STK31 高表达细胞群
#   [2] 鉴定 NK 细胞候选群（NK marker score）
#   [3] 差异分析（STK31-high vs NK、vs all，NK vs all）
#   [4] 按 sample 分组比较 STK31 表达
#   [5] 候选 ligand-receptor 互作初筛
#   [6] 细胞身份验证（身份 marker DotPlot / UMAP）
#   [7] NK 功能模块评分（细胞毒/活化/耗竭/迁移/TGF-beta）
#   [8] 通路关键基因证据图
# ============================================================

stk31_out_dir <- "/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis"
dir.create(stk31_out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 配置参数 ----------
target_gene <- "STK31"
min_pct_for_stk31_cluster <- 0.01
top_stk31_quantile <- 0.75
run_exact_wilcox <- TRUE
manual_nk_clusters <- character(0)

# NK 标志基因列表
nk_marker_genes <- c(
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "KLRD1", "KLRF1",
  "FCGR3A", "NCAM1", "TYROBP", "CTSW", "CST7", "SPON2", "XCL1", "XCL2"
)

# 粗略细胞类型标志基因
broad_celltype_markers <- list(
  Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7"),
  NK_cell = nk_marker_genes,
  B_cell = c("MS4A1", "CD79A", "CD79B", "BANK1", "CD74"),
  Plasma_cell = c("MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"),
  Myeloid = c("LYZ", "S100A8", "S100A9", "FCN1", "LST1", "CST3"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "CD68", "APOE"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "ESAM"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT")
)

# 候选 ligand-receptor 参考表
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

# ---------- 工具函数 ----------

available_genes <- function(object, genes) {
  intersect(genes, rownames(object))
}

fetch_gene_matrix <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NULL)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

mean_nonzero <- function(x) mean(x > 0)

summarize_gene_by_cluster <- function(object, gene) {
  expr <- as.numeric(fetch_gene_matrix(object, gene)[gene, ])
  clusters <- as.character(object$seurat_clusters)
  summary <- do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
    idx <- clusters == cluster_id
    data.frame(
      cluster = cluster_id, cells = sum(idx),
      avg_expr = mean(expr[idx]), median_expr = median(expr[idx]),
      pct_expr = mean_nonzero(expr[idx]), stringsAsFactors = FALSE
    )
  }))
  summary[order(-summary$avg_expr, -summary$pct_expr), ]
}

score_nk_clusters <- function(object, markers) {
  markers <- available_genes(object, markers)
  if (length(markers) == 0) stop("None of the NK marker genes were found.")
  expr <- fetch_gene_matrix(object, markers)
  cell_score <- Matrix::colMeans(expr)
  clusters <- as.character(object$seurat_clusters)
  score_table <- do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
    idx <- clusters == cluster_id
    data.frame(
      cluster = cluster_id, cells = sum(idx),
      nk_score = mean(cell_score[idx]),
      nk_marker_pct = mean(cell_score[idx] > 0),
      detected_markers = paste(markers, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  score_table[order(-score_table$nk_score, -score_table$nk_marker_pct), ]
}

score_broad_celltypes <- function(object, marker_list) {
  clusters <- as.character(object$seurat_clusters)
  scores <- do.call(rbind, lapply(names(marker_list), function(cell_type) {
    markers <- available_genes(object, marker_list[[cell_type]])
    if (length(markers) == 0) return(NULL)
    expr <- fetch_gene_matrix(object, markers)
    cell_score <- Matrix::colMeans(expr)
    do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
      idx <- clusters == cluster_id
      data.frame(
        cluster = cluster_id, cell_type = cell_type,
        score = mean(cell_score[idx]),
        marker_pct = mean(cell_score[idx] > 0),
        detected_markers = paste(markers, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))
  }))
  best <- do.call(rbind, lapply(split(scores, scores$cluster), function(x) {
    x[order(-x$score, -x$marker_pct), ][1, ]
  }))
  rownames(best) <- NULL
  list(long = scores, best = best[order(as.numeric(as.character(best$cluster))), ])
}

choose_stk31_clusters <- function(summary_table) {
  cutoff <- as.numeric(stats::quantile(summary_table$avg_expr, top_stk31_quantile, na.rm = TRUE))
  selected <- summary_table$cluster[
    summary_table$avg_expr >= cutoff & summary_table$pct_expr >= min_pct_for_stk31_cluster
  ]
  if (length(selected) == 0) selected <- summary_table$cluster[1]
  selected
}

choose_nk_clusters <- function(score_table) {
  if (length(manual_nk_clusters) > 0) return(manual_nk_clusters)
  cutoff <- as.numeric(stats::quantile(score_table$nk_score, 0.90, na.rm = TRUE))
  selected <- score_table$cluster[score_table$nk_score >= cutoff & score_table$nk_marker_pct > 0]
  if (length(selected) == 0) selected <- score_table$cluster[1]
  selected
}

write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

run_marker_test <- function(object, cells.1, cells.2, output_file) {
  if (length(cells.1) < 3 || length(cells.2) < 3) {
    warning("Too few cells for marker test: ", output_file)
    return(data.frame())
  }
  if (run_exact_wilcox) {
    object$.comparison_group <- "unused"
    object$.comparison_group[colnames(object) %in% cells.1] <- "group_1"
    object$.comparison_group[colnames(object) %in% cells.2] <- "group_2"
    markers <- FindMarkers(
      object, ident.1 = "group_1", ident.2 = "group_2",
      group.by = ".comparison_group", logfc.threshold = 0, min.pct = 0.05, test.use = "wilcox"
    )
    markers$gene <- rownames(markers)
    markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
  } else {
    expr <- fetch_gene_matrix(object, rownames(object))
    expr1 <- expr[, cells.1, drop = FALSE]
    expr2 <- expr[, cells.2, drop = FALSE]
    avg1 <- Matrix::rowMeans(expm1(expr1))
    avg2 <- Matrix::rowMeans(expm1(expr2))
    pct1 <- Matrix::rowMeans(expr1 > 0)
    pct2 <- Matrix::rowMeans(expr2 > 0)
    markers <- data.frame(
      gene = rownames(expr),
      avg_log2FC = log2((avg1 + 1e-9) / (avg2 + 1e-9)),
      avg_expr_group1 = avg1, avg_expr_group2 = avg2,
      pct.1 = pct1, pct.2 = pct2, pct_diff = pct1 - pct2,
      score = abs(log2((avg1 + 1e-9) / (avg2 + 1e-9))) * pmax(pct1, pct2),
      stringsAsFactors = FALSE
    )
    markers <- markers[pmax(markers$pct.1, markers$pct.2) >= 0.05, ]
    markers <- markers[order(-markers$score, -abs(markers$avg_log2FC)), ]
  }
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
    avg_expr = Matrix::rowMeans(expr),
    pct_expr = Matrix::rowMeans(expr > 0),
    stringsAsFactors = FALSE
  )
}

infer_lr_links <- function(object, stk31_cells, nk_cells, lr_table) {
  lr_genes <- unique(c(lr_table$ligand, lr_table$receptor))
  stk31_expr <- average_expression_for_genes(object, stk31_cells, lr_genes)
  nk_expr <- average_expression_for_genes(object, nk_cells, lr_genes)
  names(stk31_expr)[names(stk31_expr) != "gene"] <- paste0("stk31_group_", names(stk31_expr)[names(stk31_expr) != "gene"])
  names(nk_expr)[names(nk_expr) != "gene"] <- paste0("nk_group_", names(nk_expr)[names(nk_expr) != "gene"])
  forward <- merge(lr_table, stk31_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  forward <- merge(forward, nk_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  forward$direction <- "STK31_high_cell_to_NK_cell"
  reverse <- merge(lr_table, nk_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  reverse <- merge(reverse, stk31_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  reverse$direction <- "NK_cell_to_STK31_high_cell"
  reverse <- reverse[, names(forward)]
  links <- rbind(forward, reverse)
  links[is.na(links)] <- 0
  links$interaction_score <- sqrt(links$stk31_group_avg_expr * links$nk_group_avg_expr)
  links <- links[order(-links$interaction_score, -links$stk31_group_pct_expr, -links$nk_group_pct_expr), ]
  links
}

summarize_gene_by_group <- function(object, genes, group_vec, group_name = "group") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(data.frame())
  expr <- fetch_gene_matrix(object, genes)
  group_vec <- as.character(group_vec)
  groups <- sort(unique(group_vec))
  out <- do.call(rbind, lapply(groups, function(grp) {
    idx <- group_vec == grp
    data.frame(
      group = grp, gene = genes,
      avg_expr = Matrix::rowMeans(expr[, idx, drop = FALSE]),
      pct_expr = Matrix::rowMeans(expr[, idx, drop = FALSE] > 0),
      stringsAsFactors = FALSE
    )
  }))
  names(out)[1] <- group_name
  out
}

run_go_if_available <- function(markers, output_prefix) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("clusterProfiler/org.Hs.eg.db not installed; skip GO enrichment.")
    return(invisible(NULL))
  }
  if ("p_val_adj" %in% colnames(markers)) {
    sig <- markers$gene[markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25]
  } else {
    sig <- markers$gene[abs(markers$avg_log2FC) > 0.5 & abs(markers$pct_diff) > 0.05]
  }
  sig <- unique(sig[!is.na(sig)])
  if (length(sig) < 10) {
    message("Fewer than 10 significant genes; skip GO enrichment for ", output_prefix)
    return(invisible(NULL))
  }
  converted <- clusterProfiler::bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  if (nrow(converted) < 10) {
    message("Fewer than 10 converted genes; skip GO enrichment for ", output_prefix)
    return(invisible(NULL))
  }
  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  write.csv(as.data.frame(ego), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
}

# ---------- 主分析开始 ----------

if (!target_gene %in% rownames(combined)) {
  stop(target_gene, " was not found in the Seurat object.")
}

message("=== STK31-NK analysis on merged data ===")

# --- 1. STK31 表达概况 ---
message("Summarizing STK31 expression by cluster")
stk31_cluster_summary <- summarize_gene_by_cluster(combined, target_gene)
write.csv(stk31_cluster_summary, file.path(stk31_out_dir, "stk31_expression_by_cluster.csv"), row.names = FALSE)

stk31_high_clusters <- choose_stk31_clusters(stk31_cluster_summary)

# 按样本分组汇总 STK31 表达（合并数据的独有优势）
stk31_by_sample <- summarize_gene_by_group(combined, target_gene, combined$sample, group_name = "sample")
write.csv(stk31_by_sample, file.path(stk31_out_dir, "stk31_expression_by_sample.csv"), row.names = FALSE)

# --- 2. NK cluster 自动筛选 ---
message("Scoring NK-like clusters")
nk_cluster_scores <- score_nk_clusters(combined, nk_marker_genes)
write.csv(nk_cluster_scores, file.path(stk31_out_dir, "nk_marker_score_by_cluster.csv"), row.names = FALSE)
nk_clusters <- choose_nk_clusters(nk_cluster_scores)

# --- 3. 自动细胞类型注释 ---
message("Scoring broad cell types for cluster annotation hints")
broad_scores <- score_broad_celltypes(combined, broad_celltype_markers)
write.csv(broad_scores$long, file.path(stk31_out_dir, "broad_celltype_scores_by_cluster.csv"), row.names = FALSE)
write.csv(broad_scores$best, file.path(stk31_out_dir, "broad_celltype_best_guess_by_cluster.csv"), row.names = FALSE)

cluster_to_celltype <- setNames(as.character(broad_scores$best$cell_type), as.character(broad_scores$best$cluster))
combined$broad_celltype <- unname(cluster_to_celltype[as.character(combined$seurat_clusters)])
combined$broad_celltype[is.na(combined$broad_celltype)] <- "Unassigned"
combined$broad_celltype <- as.factor(combined$broad_celltype)

combined$manual_celltype <- as.character(combined$broad_celltype)
combined$manual_celltype[as.character(combined$seurat_clusters) %in% nk_clusters] <- "NK_cell"
combined$manual_celltype <- as.factor(combined$manual_celltype)
combined$cluster_celltype_label <- as.character(combined$manual_celltype)

cluster_celltype_annotation <- unique(data.frame(
  cluster = as.character(combined$seurat_clusters),
  broad_celltype = as.character(combined$broad_celltype),
  manual_celltype = as.character(combined$manual_celltype),
  nk_candidate_by_marker_score = as.character(combined$seurat_clusters) %in% nk_clusters,
  cluster_celltype_label = combined$cluster_celltype_label,
  stringsAsFactors = FALSE
))
cluster_celltype_annotation <- cluster_celltype_annotation[order(as.numeric(cluster_celltype_annotation$cluster)), ]
write.csv(cluster_celltype_annotation, file.path(stk31_out_dir, "cluster_celltype_annotation.csv"), row.names = FALSE)

# --- 4. 定义分析用分组 ---
combined$stk31_group <- ifelse(as.character(combined$seurat_clusters) %in% stk31_high_clusters, "STK31_high_cluster", "Other")
combined$nk_group <- ifelse(as.character(combined$seurat_clusters) %in% nk_clusters, "NK_cell", "Other")
combined$relationship_group <- ifelse(
  combined$stk31_group == "STK31_high_cluster", "STK31_high_cluster",
  ifelse(combined$nk_group == "NK_cell", "NK_cell", "Other")
)

stk31_cells <- colnames(combined)[combined$stk31_group == "STK31_high_cluster"]
nk_cells <- colnames(combined)[combined$nk_group == "NK_cell"]
if (length(stk31_cells) == 0) stop("No STK31 high cells found.")
if (length(nk_cells) == 0) stop("No NK cells found.")
other_cells <- colnames(combined)[combined$relationship_group == "Other"]

# --- STK31 阳性细胞分布 ---
stk31_expr <- as.numeric(fetch_gene_matrix(combined, target_gene)[target_gene, ])
stk31_positive_distribution <- as.data.frame(table(combined$seurat_clusters[stk31_expr > 0]))
colnames(stk31_positive_distribution) <- c("cluster", "stk31_positive_cells")
stk31_positive_distribution$pct_of_all_stk31_positive <- stk31_positive_distribution$stk31_positive_cells / sum(stk31_positive_distribution$stk31_positive_cells)
stk31_positive_distribution <- merge(
  stk31_positive_distribution,
  stk31_cluster_summary[, c("cluster", "cells", "avg_expr", "pct_expr")],
  by = "cluster", all.x = TRUE
)
stk31_positive_distribution <- stk31_positive_distribution[order(-stk31_positive_distribution$stk31_positive_cells), ]
write.csv(stk31_positive_distribution, file.path(stk31_out_dir, "stk31_positive_cell_distribution.csv"), row.names = FALSE)

# --- 样本级别的 STK31-positive 分布 ---
stk31_pos_by_sample <- as.data.frame(table(
  combined$sample[stk31_expr > 0]
))
colnames(stk31_pos_by_sample) <- c("sample", "stk31_positive_cells")
total_by_sample <- as.data.frame(table(combined$sample))
colnames(total_by_sample) <- c("sample", "total_cells")
stk31_pos_by_sample <- merge(stk31_pos_by_sample, total_by_sample, by = "sample")
stk31_pos_by_sample$pct_stk31_positive <- stk31_pos_by_sample$stk31_positive_cells / stk31_pos_by_sample$total_cells
write.csv(stk31_pos_by_sample, file.path(stk31_out_dir, "stk31_positive_by_sample.csv"), row.names = FALSE)

# 输出分析总览表
summary_table <- data.frame(
  item = c(
    "target_gene", "stk31_high_clusters", "stk31_high_cells", "nk_candidate_clusters",
    "nk_candidate_cells", "other_cells", "detected_nk_markers", "samples_analyzed"
  ),
  value = c(
    target_gene,
    paste(stk31_high_clusters, collapse = ";"),
    length(stk31_cells),
    paste(nk_clusters, collapse = ";"),
    length(nk_cells),
    length(other_cells),
    paste(available_genes(combined, nk_marker_genes), collapse = ";"),
    paste(unique(combined$sample), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(stk31_out_dir, "analysis_summary.csv"), row.names = FALSE)

# ---------- 生成可视化 ----------
message("Writing STK31 and NK visualization")

write_plot(
  file.path(stk31_out_dir, "umap_broad_celltype_annotation.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "broad_celltype", label = TRUE, repel = TRUE) +
    ggtitle("Marker-based broad cell type annotation (merged)")
)
write_plot(
  file.path(stk31_out_dir, "umap_manual_celltype_annotation.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "manual_celltype", label = TRUE, repel = TRUE) +
    ggtitle("Manual annotation with NK marker score candidates (merged)")
)
write_plot(
  file.path(stk31_out_dir, "umap_cluster_celltype_labels.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "cluster_celltype_label", label = TRUE, repel = TRUE) +
    ggtitle("Cluster labels with manual cell type annotation (merged)"),
  width = 10, height = 7
)
write_plot(
  file.path(stk31_out_dir, "umap_stk31_expression.pdf"),
  FeaturePlot(combined, features = target_gene, reduction = "umap", order = TRUE) + ggtitle("STK31 expression (merged)")
)
write_plot(
  file.path(stk31_out_dir, "umap_stk31_high_and_nk_groups.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "relationship_group", label = FALSE) +
    ggtitle("STK31 high clusters and NK candidate clusters (merged)")
)
write_plot(
  file.path(stk31_out_dir, "umap_by_sample_merged.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "sample", label = FALSE) +
    ggtitle("Merged UMAP by sample")
)
write_plot(
  file.path(stk31_out_dir, "stk31_violin_by_cluster.pdf"),
  VlnPlot(combined, features = target_gene, group.by = "seurat_clusters", pt.size = 0.01) +
    ggtitle("STK31 by cluster (merged)"),
  width = 10, height = 5
)
write_plot(
  file.path(stk31_out_dir, "nk_marker_dotplot.pdf"),
  DotPlot(combined, features = available_genes(combined, nk_marker_genes), group.by = "seurat_clusters") +
    RotatedAxis() + ggtitle("NK marker expression by cluster (merged)"),
  width = 12, height = 6
)
write_plot(
  file.path(stk31_out_dir, "broad_celltype_marker_dotplot.pdf"),
  DotPlot(combined, features = unique(unlist(lapply(broad_celltype_markers, function(x) available_genes(combined, x)))), group.by = "seurat_clusters") +
    RotatedAxis() + ggtitle("Broad cell type markers by cluster (merged)"),
  width = 16, height = 7
)

# --- UMAP 中心距离：STK31-high vs NK_cell ---
if ("umap" %in% names(combined@reductions)) {
  umap <- Embeddings(combined, "umap")
  centroids <- aggregate(umap, by = list(group = combined$relationship_group), FUN = mean)
  write.csv(centroids, file.path(stk31_out_dir, "relationship_group_umap_centroids.csv"), row.names = FALSE)
  if (all(c("STK31_high_cluster", "NK_cell") %in% centroids$group)) {
    umap_cols <- setdiff(colnames(centroids), "group")[1:2]
    stk31_center <- as.numeric(centroids[centroids$group == "STK31_high_cluster", umap_cols])
    nk_center <- as.numeric(centroids[centroids$group == "NK_cell", umap_cols])
    distance_table <- data.frame(
      comparison = "STK31_high_cluster_vs_NK_cell",
      umap_centroid_distance = sqrt(sum((stk31_center - nk_center)^2))
    )
    write.csv(distance_table, file.path(stk31_out_dir, "stk31_nk_umap_centroid_distance.csv"), row.names = FALSE)
  }
}

# ---------- 差异分析 ----------
message("Running differential expression tests")
stk31_vs_nk <- run_marker_test(
  combined, cells.1 = stk31_cells, cells.2 = nk_cells,
  output_file = file.path(stk31_out_dir, "markers_stk31_high_vs_nk.csv")
)
stk31_vs_other <- run_marker_test(
  combined, cells.1 = stk31_cells, cells.2 = c(nk_cells, other_cells),
  output_file = file.path(stk31_out_dir, "markers_stk31_high_vs_all_other.csv")
)
nk_vs_other <- run_marker_test(
  combined, cells.1 = nk_cells, cells.2 = c(stk31_cells, other_cells),
  output_file = file.path(stk31_out_dir, "markers_nk_vs_all_other.csv")
)

# ---------- 候选 ligand-receptor 互作 ----------
message("Inferring candidate ligand-receptor pathways")
lr_links <- infer_lr_links(combined, stk31_cells, nk_cells, lr_reference)
write.csv(lr_links, file.path(stk31_out_dir, "candidate_ligand_receptor_pathways.csv"), row.names = FALSE)

pathway_summary <- aggregate(
  interaction_score ~ direction + pathway,
  data = lr_links, FUN = max
)
pathway_summary <- pathway_summary[order(-pathway_summary$interaction_score), ]
write.csv(pathway_summary, file.path(stk31_out_dir, "candidate_pathway_summary.csv"), row.names = FALSE)

# --- GO 富集（可选）---
message("Trying optional GO enrichment")
run_go_if_available(stk31_vs_nk, file.path(stk31_out_dir, "stk31_high_vs_nk"))
run_go_if_available(stk31_vs_other, file.path(stk31_out_dir, "stk31_high_vs_all_other"))
run_go_if_available(nk_vs_other, file.path(stk31_out_dir, "nk_vs_all_other"))

# --- 每个 cluster 的 marker 基因（用于验证 STK31-high cluster 的身份）---
message("Finding all cluster markers (this may take a while on merged data)")
all_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  all_markers,
  file.path(stk31_out_dir, "all_cluster_markers.csv"),
  row.names = FALSE
)

message("STK31-NK analysis done. Key outputs:")
message("- ", file.path(stk31_out_dir, "all_cluster_markers.csv"))
message("- ", file.path(stk31_out_dir, "stk31_expression_by_cluster.csv"))
message("- ", file.path(stk31_out_dir, "stk31_expression_by_sample.csv"))
message("- ", file.path(stk31_out_dir, "nk_marker_score_by_cluster.csv"))
message("- ", file.path(stk31_out_dir, "markers_stk31_high_vs_nk.csv"))
message("- ", file.path(stk31_out_dir, "candidate_ligand_receptor_pathways.csv"))
