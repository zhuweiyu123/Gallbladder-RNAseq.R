suppressPackageStartupMessages({
  # 加载单细胞分析和绘图需要的 R 包；suppressPackageStartupMessages 用来隐藏启动提示。
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# 固定随机种子，让 PCA/UMAP/聚类等带随机性的步骤尽量可重复。
set.seed(20260624)

# 服务器路径：这个脚本设计为在 aiserver 上运行，不是在本地 Windows 上运行。
# base_dir 是 10x 原始矩阵所在目录；out_dir 是本次 tissue2 分析结果输出目录。
base_dir <- "/home/zhuweiyu/codex-r/gallbladder_cancer"
out_dir <- "/home/zhuweiyu/codex-r/results/tissue2_basic_seurat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 本脚本只分析 tissue2 一个样本。
samples <- c("tissue2")

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
# 单样本分析时不需要 merge；保留这个判断是为了脚本以后扩展到多样本时也能用。
if (length(filtered) == 1) {
  combined <- filtered[[1]]
} else {
  combined <- merge(
    x = filtered[[1]],
    y = filtered[-1],
    add.cell.ids = samples,
    project = "gallbladder_cancer"
  )
}

# 画 QC 小提琴图：基因数、UMI 数、线粒体比例。
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

# 按样本来源上色的 UMAP。单样本时主要用于确认图正常生成。
pdf(file.path(out_dir, "umap_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(combined, reduction = "umap", group.by = "sample"))
dev.off()

# 按聚类结果上色的 UMAP，并在图上标出 cluster 编号。
pdf(file.path(out_dir, "umap_by_cluster.pdf"), width = 8, height = 6)
print(DimPlot(combined, reduction = "umap", label = TRUE))
dev.off()

# 统计每个 cluster 中来自各样本的细胞数。单样本时就是 tissue2 的 cluster 细胞数。
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

# 保存完整 Seurat 对象，后续可以直接读取它继续做 marker、注释、差异分析等。
saveRDS(combined, file.path(out_dir, "gallbladder_cancer_basic_seurat.rds"))

message("Done. Outputs written to: ", out_dir)
