# ========================================================================
# 【中文阅读指南】单样本 tissue2：从原始计数到聚类图
# 输入：服务器 gallbladder_cancer/tissue2 下的 10x 矩阵；矩阵的行是基因、列是细胞。
# 流程：读取数据 → 质控过滤 → 标准化 → 高变基因 → PCA → 邻居图/聚类 → UMAP。
# 输出：results/tissue2_basic_seurat 下的 QC 表、图和 gallbladder_cancer_basic_seurat.rds。
# 阅读重点：先看质控阈值，再看过滤前后细胞数；完成后可接着阅读脚本 03。
# 聚类编号只是算法分组，不能直接当作细胞类型；后续要结合标志基因判断。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  # 加载单细胞分析和绘图需要的 R 包；suppressPackageStartupMessages 用来隐藏启动提示。
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# 固定随机种子，让 PCA/UMAP/聚类等带随机性的步骤尽量可重复。
# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260624)

# 服务器路径：这个脚本设计为在 aiserver 上运行，不是在本地 Windows 上运行。
# base_dir 是 10x 原始矩阵所在目录；out_dir 是本次 tissue2 分析结果输出目录。
base_dir <- "/home/zhuweiyu/codex-r/gallbladder_cancer"
out_dir <- "/home/zhuweiyu/codex-r/results/tissue2_basic_seurat"
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
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

# 【函数：read_sample】读取一个 sample_id 的 10x 计数并建立 Seurat 对象，记录样本来源及线粒体比例。
# 返回值是该样本的细胞对象；样本标签用于后续分组和条形码区分。
read_sample <- function(sample_id) {
  # 每个样本的 10x filtered_feature_bc_matrix 目录。
  data_dir <- file.path(base_dir, sample_id, "filtered_feature_bc_matrix")
  if (!dir.exists(data_dir)) {
    stop("Missing 10x directory: ", data_dir)
  }

  # 读取 10x 数据。部分 10x 输出会返回 list，这里优先取 Gene Expression 矩阵。
  # 【10x 数据】读取 matrix、features 和 barcodes；通常返回基因×细胞的稀疏原始计数矩阵。
  counts <- Read10X(data.dir = data_dir)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  # 创建 Seurat 对象，并记录样本名，方便后面按 sample 分组画图。
  # 【创建对象】把 counts 与细胞信息装进 Seurat；对象同时保存表达矩阵、元数据及后续降维结果。
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells_per_gene,
    min.features = min_features_per_cell
  )
  obj$sample <- sample_id

  # 计算每个细胞线粒体基因比例。人类基因名通常以 MT- 开头。
  # 【线粒体比例】按匹配的线粒体基因统计每个细胞的计数占比；percent.mt 通常按 0–100 的百分数保存。
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

# 【函数：qc_table】输入多个样本的对象列表，逐样本统计细胞数和质控指标。
# 返回可写入 CSV 的表，用过滤前/后两张表检查数据损失。
qc_table <- function(object_list) {
  # 汇总每个样本的细胞数、基因数、UMI 数和线粒体比例，输出成 CSV 方便检查。
  # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
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
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
write.csv(
  qc_table(objects),
  file.path(out_dir, "qc_before_filter.csv"),
  row.names = FALSE
)

# 按前面设定的阈值过滤细胞。
filtered <- lapply(objects, function(obj) {
  # 【取细胞子集】按 cells 或条件保留需要的细胞；这一步改变本次分析范围，要同步核对后面的分母。
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
# 【标准化】默认 LogNormalize 将每个细胞按总计数缩放再取 log1p，减少测序深度差异的影响。
combined <- NormalizeData(combined)

# FindVariableFeatures：找高变基因，后续 PCA 主要基于这些信息量更高的基因。
# 【高变基因】选择细胞间变化较大的基因用于降维，nfeatures 控制数量；不等于删掉其他基因的原始表达。
combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = n_variable_features
)

# ScaleData：对高变基因做中心化/标准化。第一轮只 scale 高变基因，可以明显降低内存占用。
# 【缩放】对用于分析的基因做中心化/标准化，使不同基因更便于进入 PCA；具体基因由 features 指定。
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)

# RunPCA：把高维表达矩阵压缩成主成分，便于后续建邻居图和可视化。
# 【PCA】把大量基因的变化压缩成主成分；npcs 是计算数量，后续 dims 决定实际使用哪些。
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = n_pcs, verbose = FALSE)

# FindNeighbors + FindClusters：根据 PCA 空间中的相似性构建细胞图，并进行聚类。
# 【邻居图】在选定主成分空间中寻找表达相似的细胞，为图聚类提供连接关系。
combined <- FindNeighbors(combined, dims = use_dims)
# 【聚类】利用邻居图划分细胞群；resolution 越大通常分得越细，编号本身没有生物学身份。
combined <- FindClusters(combined, resolution = cluster_resolution)

# RunUMAP：把细胞降到二维，方便观察 cluster 和样本分布。
# 【UMAP】将表达相似性压缩到低维便于展示；二维距离不能解释为组织空间距离。
combined <- RunUMAP(combined, dims = use_dims)

# 按样本来源上色的 UMAP。单样本时主要用于确认图正常生成。
pdf(file.path(out_dir, "umap_by_sample.pdf"), width = 8, height = 6)
# 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
# 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
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
# 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
saveRDS(combined, file.path(out_dir, "gallbladder_cancer_basic_seurat.rds"))

message("Done. Outputs written to: ", out_dir)
