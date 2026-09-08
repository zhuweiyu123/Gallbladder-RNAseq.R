# ========================================================================
# 【中文阅读指南】四样本合并分析与 STK31/NK 主流程
# 输入：tissue1、tissue2、tissue4、tissue5 的 10x 数据；路径在 base_dir 中。
# 流程：分别质控 → 合并 → 标准化和聚类 → 细胞注释 → STK31 分组 → 差异/GO/通讯及后续图。
# 输出：results/merged_basic_seurat 和 results/merged_stk31_nk_analysis 及其子目录。
# 这是包含多段历史后续分析的长脚本；中途会重新定义 out_dir 和部分函数，按从上到下顺序阅读。
# merge 合并细胞数据，并不自动完成批次效应校正；UMAP 的样本混合情况需要另行判断。
# 此处 Epithelial 是标志基因注释；名称中的 tumor 不代表已完成 CNV 恶性鉴定。
# high 通常由上皮细胞 STK31 的 75% 分位数确定；阈值为零时改用表达 > 0，详见分组函数。
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
# base_dir 是 10x 原始矩阵所在目录；out_dir 是四个样本合并分析的结果输出目录。
base_dir <- "/home/zhuweiyu/codex-r/gallbladder_cancer"
out_dir <- "/home/zhuweiyu/codex-r/results/merged_basic_seurat"
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
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

# 按样本来源上色的 UMAP，用来观察不同样本是否混合、是否有明显批次差异。
pdf(file.path(out_dir, "umap_by_sample.pdf"), width = 8, height = 6)
# 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
# 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
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

if (packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  combined <- JoinLayers(combined, assay = "RNA")
}

# 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
saveRDS(combined, file.path(out_dir, "gallbladder_cancer_merged_basic_seurat.rds"))

message("Basic analysis done. Outputs written to: ", out_dir)

# ============================================================
# STK31-NK 完整分析（在合并数据上）
# 用途：胆囊癌合并数据（tissue1/2/4/5），聚焦 STK31
#   [1] 定位 STK31 高表达肿瘤上皮细胞群
#   [2] 鉴定 NK 细胞候选群（NK marker score）
#   [3] 差异分析（肿瘤上皮细胞内 STK31 high vs low、STK31-high 肿瘤上皮 vs NK/vs all，NK vs all）
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
tumor_epithelial_celltype_label <- "Epithelial"
# 【参数 tumor_epithelial_stk31_high_quantile】上皮内 high 的分位数设定，0.75 表示第 75 百分位；实际阈值和并列值会影响 high 数量。
tumor_epithelial_stk31_high_quantile <- 0.75
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

# 【函数：available_genes】把想看的基因与对象实际行名取交集，避免访问不存在的基因。
# 返回基因名向量；返回长度为零表示这些基因不在当前矩阵内。
available_genes <- function(object, genes) {
  # 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
  intersect(genes, rownames(object))
}

# 【函数：fetch_gene_matrix】从 Seurat 对象提取指定基因的表达矩阵，并兼容不同版本的 layer/slot 接口。
# 返回矩阵通常是基因×细胞；data 是标准化层，counts 是原始计数层。
fetch_gene_matrix <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NULL)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

# 【函数：mean_nonzero】计算 x>0 的比例：TRUE 按 1、FALSE 按 0 求平均。
# 虽然函数名含 mean，但这里返回的是检测比例，并非阳性细胞的平均表达量。
mean_nonzero <- function(x) mean(x > 0)

# 【函数：summarize_gene_by_cluster】对指定基因逐聚类统计表达水平和检出比例，返回聚类摘要表。
# 用于定位 STK31 信号集中在哪些群；群均值不能替代每个细胞的表达状态。
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

# 【函数：score_nk_clusters】用可检测到的 NK 标志基因给每个聚类计算平均表达分数。
# 返回群级评分用于候选筛选；T 细胞也可能表达部分细胞毒标志，因此需进一步验证。
score_nk_clusters <- function(object, markers) {
  markers <- available_genes(object, markers)
  if (length(markers) == 0) stop("None of the NK marker genes were found.")
  expr <- fetch_gene_matrix(object, markers)
  # 【按细胞求均值】在基因×细胞矩阵中，colMeans 对每列跨基因求平均，常用来构建逐细胞模块分数。
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

# 【函数：score_broad_celltypes】按各类标志基因集合计算聚类分数，再选择得分最高的类型作为初步注释。
# 这是依赖人工标志列表的启发式分类；相近分数或混合标志需要核查。
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

# 【函数：choose_stk31_clusters】根据聚类级 STK31 平均表达分位数等条件选出候选聚类。
# 返回群编号；这一步的群筛选与后续上皮细胞内 high/low 分组是两个不同层次。
choose_stk31_clusters <- function(summary_table) {
  cutoff <- as.numeric(stats::quantile(summary_table$avg_expr, top_stk31_quantile, na.rm = TRUE))
  selected <- summary_table$cluster[
    summary_table$avg_expr >= cutoff & summary_table$pct_expr >= min_pct_for_stk31_cluster
  ]
  if (length(selected) == 0) selected <- summary_table$cluster[1]
  selected
}

# 【函数：choose_nk_clusters】依据 NK 群评分的高分位阈值等条件筛选 NK 候选群。
# 返回候选聚类编号；自动筛选不等于已确认所有入选细胞都是 NK。
choose_nk_clusters <- function(score_table) {
  if (length(manual_nk_clusters) > 0) return(manual_nk_clusters)
  cutoff <- as.numeric(stats::quantile(score_table$nk_score, 0.90, na.rm = TRUE))
  selected <- score_table$cluster[score_table$nk_score >= cutoff & score_table$nk_marker_pct > 0]
  if (length(selected) == 0) selected <- score_table$cluster[1]
  selected
}

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 【函数：run_marker_test】输入对象和两组细胞，调用差异表达检验并把基因表写入 output_file。
# 重点看比较方向、avg_log2FC 和调整后 P 值；正 logFC 表示第一组相对第二组更高。
run_marker_test <- function(object, cells.1, cells.2, output_file) {
  if (length(cells.1) < 3 || length(cells.2) < 3) {
    warning("Too few cells for marker test: ", output_file)
    return(data.frame())
  }
  if (run_exact_wilcox) {
    object$.comparison_group <- "unused"
    object$.comparison_group[colnames(object) %in% cells.1] <- "group_1"
    object$.comparison_group[colnames(object) %in% cells.2] <- "group_2"
    # 【差异表达】比较 ident.1 与 ident.2；avg_log2FC>0 表示第一组更高，p_val_adj 是多重检验调整后 P 值。
    # min.pct/logfc.threshold 是进入检验的筛选条件；细胞级检验不自动控制患者内相关性。
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
    # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
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

# 【函数：average_expression_for_genes】对指定 cells 中每个目标基因计算表达摘要，供配体受体候选评分使用。
# 返回按基因组织的结果；表达矩阵里缺失的基因应与真实零表达区分。
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

# 【函数：infer_lr_links】沿着候选表逐个读取配体与受体，结合发送细胞和接收细胞的表达摘要打分。
# 这是人工候选配对的描述性筛选；不能把表达乘积当成已证实的分子结合或通讯概率。
infer_lr_links <- function(object, tumor_epithelial_high_cells, nk_cells, lr_table) {
  lr_genes <- unique(c(lr_table$ligand, lr_table$receptor))
  tumor_epithelial_high_expr <- average_expression_for_genes(object, tumor_epithelial_high_cells, lr_genes)
  nk_expr <- average_expression_for_genes(object, nk_cells, lr_genes)
  names(tumor_epithelial_high_expr)[names(tumor_epithelial_high_expr) != "gene"] <- paste0("tumor_epithelial_stk31_high_", names(tumor_epithelial_high_expr)[names(tumor_epithelial_high_expr) != "gene"])
  names(nk_expr)[names(nk_expr) != "gene"] <- paste0("nk_group_", names(nk_expr)[names(nk_expr) != "gene"])
  forward <- merge(lr_table, tumor_epithelial_high_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  forward <- merge(forward, nk_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  forward$direction <- "STK31_high_tumor_epithelial_to_NK_cell"
  reverse <- merge(lr_table, nk_expr, by.x = "ligand", by.y = "gene", all.x = TRUE)
  reverse <- merge(reverse, tumor_epithelial_high_expr, by.x = "receptor", by.y = "gene", all.x = TRUE)
  reverse$direction <- "NK_cell_to_STK31_high_tumor_epithelial"
  reverse <- reverse[, names(forward)]
  links <- rbind(forward, reverse)
  links[is.na(links)] <- 0
  links$interaction_score <- sqrt(links$tumor_epithelial_stk31_high_avg_expr * links$nk_group_avg_expr)
  links <- links[order(-links$interaction_score, -links$tumor_epithelial_stk31_high_pct_expr, -links$nk_group_pct_expr), ]
  links
}

# 【函数：summarize_gene_by_group】按 group_vec 给每个细胞分组，再汇总所选基因的平均表达和检测情况。
# group_vec 顺序应与矩阵细胞列一致；返回长表，便于画图和写 CSV。
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

# 【函数：assign_tumor_epithelial_stk31_group】只在指定上皮细胞内计算 STK31 分位数阈值，随后生成 high/low 标签。
# 阈值>0 时使用 >= 阈值；阈值为零时使用 >0；其余细胞保持 Other。
# 返回 list，含每个细胞的标签、实际阈值与各组数量；先看这些数量再解读下游差异。
assign_tumor_epithelial_stk31_group <- function(object, stk31_expr, celltype_label, high_quantile = 0.75) {
  tumor_epithelial_idx <- as.character(object$manual_celltype) == celltype_label
  if (sum(tumor_epithelial_idx) == 0) {
    stop("No tumor epithelial cells found with manual_celltype == ", celltype_label)
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
    group = group,
    cutoff = cutoff,
    tumor_epithelial_cells = sum(tumor_epithelial_idx),
    high_cells = sum(high_idx),
    low_cells = sum(low_idx)
  )
}

# 【函数：run_go_if_available】检查 GO 所需包及可用差异基因，再尝试执行富集并保存结果。
# 不足以分析时按本函数分支跳过；应区分未运行、无显著富集和成功检出。
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
  # 【基因 ID】将常见基因符号转换为数据库识别的 ENTREZID；有些符号无法匹配或一对多，需看转换后数量。
  converted <- clusterProfiler::bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  if (nrow(converted) < 10) {
    message("Fewer than 10 converted genes; skip GO enrichment for ", output_prefix)
    return(invisible(NULL))
  }
  # 【GO 过度代表分析】检查候选基因是否集中于某些 GO 条目；BP 是生物过程，BH 用于多重检验调整。
  # 如果没有显式传入 universe，就使用数据库可用背景；富集不直接说明通路被激活。
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

legacy_stk31_high_clusters_qc_only <- choose_stk31_clusters(stk31_cluster_summary)

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
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
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
stk31_expr <- as.numeric(fetch_gene_matrix(combined, target_gene)[target_gene, ])
tumor_epithelial_assignment <- assign_tumor_epithelial_stk31_group(
  combined,
  stk31_expr,
  celltype_label = tumor_epithelial_celltype_label,
  high_quantile = tumor_epithelial_stk31_high_quantile
)
combined$tumor_epithelial_stk31_group <- tumor_epithelial_assignment$group
combined$stk31_group <- combined$tumor_epithelial_stk31_group
combined$nk_group <- ifelse(as.character(combined$seurat_clusters) %in% nk_clusters, "NK_cell", "Other")
combined$relationship_group <- ifelse(
  combined$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial",
  ifelse(
    combined$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial", "STK31_low_tumor_epithelial",
    ifelse(combined$nk_group == "NK_cell", "NK_cell", "Other")
  )
)
combined$relationship_group <- factor(
  combined$relationship_group,
  levels = c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell", "Other")
)
combined$manual_celltype_stk31_tumor_epithelial_split <- as.character(combined$manual_celltype)
combined$manual_celltype_stk31_tumor_epithelial_split[
  combined$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"
] <- "STK31_high_tumor_epithelial"
combined$manual_celltype_stk31_tumor_epithelial_split[
  combined$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"
] <- "STK31_low_tumor_epithelial"
combined$manual_celltype_stk31_tumor_epithelial_split <- factor(combined$manual_celltype_stk31_tumor_epithelial_split)

tumor_epithelial_stk31_high_cells <- colnames(combined)[combined$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"]
tumor_epithelial_stk31_low_cells <- colnames(combined)[combined$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"]
nk_cells <- colnames(combined)[combined$nk_group == "NK_cell"]
if (length(tumor_epithelial_stk31_high_cells) == 0) stop("No STK31-high tumor epithelial cells found.")
if (length(tumor_epithelial_stk31_low_cells) == 0) stop("No STK31-low tumor epithelial cells found.")
if (length(nk_cells) == 0) stop("No NK cells found.")
other_cells <- colnames(combined)[combined$relationship_group == "Other"]

# --- STK31 阳性细胞分布 ---
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

# 输出分析总览。
summary_table <- data.frame(
  item = c(
    "target_gene", "tumor_epithelial_celltype_label", "tumor_epithelial_stk31_high_cutoff",
    "tumor_epithelial_cells", "stk31_high_tumor_epithelial_cells", "stk31_low_tumor_epithelial_cells",
    "legacy_stk31_high_clusters_qc_only", "nk_candidate_clusters", "nk_candidate_cells", "other_cells",
    "detected_nk_markers", "samples_analyzed"
  ),
  value = c(
    target_gene,
    tumor_epithelial_celltype_label,
    tumor_epithelial_assignment$cutoff,
    tumor_epithelial_assignment$tumor_epithelial_cells,
    length(tumor_epithelial_stk31_high_cells),
    length(tumor_epithelial_stk31_low_cells),
    paste(legacy_stk31_high_clusters_qc_only, collapse = ";"),
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
  file.path(stk31_out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  DimPlot(combined, reduction = "umap", group.by = "manual_celltype_stk31_tumor_epithelial_split", label = TRUE, repel = TRUE) +
    ggtitle("Manual annotation with STK31 high/low tumor epithelial split (merged)")
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
    ggtitle("STK31 high/low tumor epithelial cells and NK cells (merged)")
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

# --- UMAP 中心距离：STK31-high tumor epithelial vs NK_cell ---
if ("umap" %in% names(combined@reductions)) {
  umap <- Embeddings(combined, "umap")
  centroids <- aggregate(umap, by = list(group = combined$relationship_group), FUN = mean)
  write.csv(centroids, file.path(stk31_out_dir, "relationship_group_umap_centroids.csv"), row.names = FALSE)
  if (all(c("STK31_high_tumor_epithelial", "NK_cell") %in% centroids$group)) {
    # 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
    umap_cols <- setdiff(colnames(centroids), "group")[1:2]
    stk31_center <- as.numeric(centroids[centroids$group == "STK31_high_tumor_epithelial", umap_cols])
    nk_center <- as.numeric(centroids[centroids$group == "NK_cell", umap_cols])
    distance_table <- data.frame(
      comparison = "STK31_high_tumor_epithelial_vs_NK_cell",
      umap_centroid_distance = sqrt(sum((stk31_center - nk_center)^2))
    )
    write.csv(distance_table, file.path(stk31_out_dir, "stk31_nk_umap_centroid_distance.csv"), row.names = FALSE)
  }
}

# ---------- 差异分析 ----------
message("Running differential expression tests")
tumor_epithelial_high_vs_low <- run_marker_test(
  combined, cells.1 = tumor_epithelial_stk31_high_cells, cells.2 = tumor_epithelial_stk31_low_cells,
  output_file = file.path(stk31_out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv")
)
tumor_epithelial_high_vs_nk <- run_marker_test(
  combined, cells.1 = tumor_epithelial_stk31_high_cells, cells.2 = nk_cells,
  output_file = file.path(stk31_out_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv")
)
tumor_epithelial_high_vs_other <- run_marker_test(
  combined, cells.1 = tumor_epithelial_stk31_high_cells, cells.2 = c(nk_cells, tumor_epithelial_stk31_low_cells, other_cells),
  output_file = file.path(stk31_out_dir, "markers_stk31_high_tumor_epithelial_vs_all_other.csv")
)
nk_vs_other <- run_marker_test(
  combined, cells.1 = nk_cells, cells.2 = c(tumor_epithelial_stk31_high_cells, tumor_epithelial_stk31_low_cells, other_cells),
  output_file = file.path(stk31_out_dir, "markers_nk_vs_all_other.csv")
)

# ---------- 候选 ligand-receptor 互作 ----------
message("Inferring candidate ligand-receptor pathways")
lr_links <- infer_lr_links(combined, tumor_epithelial_stk31_high_cells, nk_cells, lr_reference)
write.csv(lr_links, file.path(stk31_out_dir, "candidate_ligand_receptor_pathways.csv"), row.names = FALSE)

pathway_summary <- aggregate(
  interaction_score ~ direction + pathway,
  data = lr_links, FUN = max
)
pathway_summary <- pathway_summary[order(-pathway_summary$interaction_score), ]
write.csv(pathway_summary, file.path(stk31_out_dir, "candidate_pathway_summary.csv"), row.names = FALSE)

# --- GO 富集（可选）---
message("Trying optional GO enrichment")
run_go_if_available(tumor_epithelial_high_vs_low, file.path(stk31_out_dir, "stk31_high_vs_low_tumor_epithelial"))
run_go_if_available(tumor_epithelial_high_vs_nk, file.path(stk31_out_dir, "stk31_high_tumor_epithelial_vs_nk"))
run_go_if_available(tumor_epithelial_high_vs_other, file.path(stk31_out_dir, "stk31_high_tumor_epithelial_vs_all_other"))
run_go_if_available(nk_vs_other, file.path(stk31_out_dir, "nk_vs_all_other"))

# --- 每个 cluster 的 marker 基因（仅用于细胞身份 QC，不作为 STK31-high 分组依据）---
message("Finding all cluster markers (this may take a while on merged data)")
# 【聚类标志】逐群与其他细胞比较，寻找解释细胞身份的 marker；仍需要检查标志组合及样本组成。
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
message("- ", file.path(stk31_out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv"))
message("- ", file.path(stk31_out_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv"))
message("- ", file.path(stk31_out_dir, "candidate_ligand_receptor_pathways.csv"))

# ============================================================
# 【后续区段 1】重新读取基础对象和注释，按细胞类型拆分 STK31 信号，并汇总每个样本中的 NK 状态。
# Integrated from 05_merged_stk31_nk_followup.R
# This section used to live in a separate follow-up script.
# ============================================================
# ============================================================
# merged STK31-NK follow-up analysis
# ============================================================
# Purpose:
#   [1] Split STK31 signal by cell type after merged analysis
#   [2] Re-test STK31-positive vs STK31-negative inside epithelial cells
#   [3] Summarize STK31-NK relationship within each sample
#   [4] Score NK functional states by sample
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260630)

basic_result_dir <- "/home/zhuweiyu/codex-r/results/merged_basic_seurat"
main_result_dir <- "/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis"
input_rds <- file.path(basic_result_dir, "gallbladder_cancer_merged_basic_seurat.rds")
out_dir <- file.path(main_result_dir, "followup")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
# 【参数 min_cells_for_de】差异分析两组需达到的最低细胞数；过小的群会停止或跳过检验。
min_cells_for_de <- 3
run_exact_wilcox <- TRUE

nk_function_sets <- list(
  NK_cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH"),
  NK_activation = c("IFNG", "TNF", "CD69", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5"),
  NK_checkpoint = c("TIGIT", "HAVCR2", "LAG3", "PDCD1", "CTLA4", "TOX"),
  NK_migration = c("CXCR3", "CCR5", "XCR1", "CX3CR1", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2"),
  TGFbeta_response = c("TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "TAGLN", "ACTA2")
)

identity_marker_sets <- list(
  Epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"),
  NK_cell = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"),
  Myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCN1", "TYROBP", "CST3"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "CD68"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "ESAM")
)

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
}

# 【函数：available_genes】把想看的基因与对象实际行名取交集，避免访问不存在的基因。
# 返回基因名向量；返回长度为零表示这些基因不在当前矩阵内。
available_genes <- function(object, genes) {
  intersect(genes, rownames(object))
}

# 【函数：fetch_gene_matrix】从 Seurat 对象提取指定基因的表达矩阵，并兼容不同版本的 layer/slot 接口。
# 返回矩阵通常是基因×细胞；data 是标准化层，counts 是原始计数层。
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

# 【函数：read_summary_item】从 item/value 格式的旧汇总表读取指定项目，并把分号分隔的值拆开。
# 文件、列或项目缺失时按函数分支返回 fallback。
read_summary_item <- function(path, item_name, fallback = character(0)) {
  if (!file.exists(path)) {
    return(fallback)
  }
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  summary_df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("item", "value") %in% colnames(summary_df))) {
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

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 【函数：summarize_gene_by_group】按 group_vec 给每个细胞分组，再汇总所选基因的平均表达和检测情况。
# group_vec 顺序应与矩阵细胞列一致；返回长表，便于画图和写 CSV。
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

# 【函数：run_marker_test】输入对象和两组细胞，调用差异表达检验并把基因表写入 output_file。
# 重点看比较方向、avg_log2FC 和调整后 P 值；正 logFC 表示第一组相对第二组更高。
run_marker_test <- function(object, cells.1, cells.2, output_file) {
  if (length(cells.1) < min_cells_for_de || length(cells.2) < min_cells_for_de) {
    warning("Too few cells for marker test: ", output_file)
    empty <- data.frame()
    write.csv(empty, output_file, row.names = FALSE)
    return(empty)
  }

  if (run_exact_wilcox) {
    object$.comparison_group <- "unused"
    object$.comparison_group[colnames(object) %in% cells.1] <- "group_1"
    object$.comparison_group[colnames(object) %in% cells.2] <- "group_2"
    markers <- FindMarkers(
      object,
      ident.1 = "group_1",
      ident.2 = "group_2",
      group.by = ".comparison_group",
      logfc.threshold = 0,
      min.pct = 0.05,
      test.use = "wilcox"
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
      avg_expr_group1 = avg1,
      avg_expr_group2 = avg2,
      pct.1 = pct1,
      pct.2 = pct2,
      pct_diff = pct1 - pct2,
      score = abs(log2((avg1 + 1e-9) / (avg2 + 1e-9))) * pmax(pct1, pct2),
      stringsAsFactors = FALSE
    )
    markers <- markers[pmax(markers$pct.1, markers$pct.2) >= 0.05, ]
    markers <- markers[order(-markers$score, -abs(markers$avg_log2FC)), ]
  }

  write.csv(markers, output_file, row.names = FALSE)
  markers
}

# 【函数：module_score】对可用基因的标准化表达逐细胞取平均，得到基因集合的简易表达分数。
# 本函数没有 AddModuleScore 的匹配背景扣除；分数不能直接视为功能活性的实测值。
module_score <- function(object, genes) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(object)))
  }
  expr <- fetch_gene_matrix(object, genes)
  as.numeric(Matrix::colMeans(expr))
}

# 【函数：summarize_numeric_by_group】把数值向量 values 按 group_vec 汇总，value_name 用于标记指标。
# 用于把逐细胞分数压缩成样本或聚类级摘要。
summarize_numeric_by_group <- function(values, group_vec, value_name, group_name = "group") {
  group_vec <- as.character(group_vec)
  groups <- sort(unique(group_vec))
  out <- do.call(rbind, lapply(groups, function(grp) {
    idx <- group_vec == grp & !is.na(values)
    data.frame(
      group = grp,
      cells = sum(group_vec == grp),
      mean_score = ifelse(any(idx), mean(values[idx]), NA_real_),
      median_score = ifelse(any(idx), median(values[idx]), NA_real_),
      pct_score_positive = ifelse(any(idx), mean(values[idx] > 0), NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  names(out)[1] <- group_name
  out$score_name <- value_name
  out
}

stop_if_missing(input_rds)
stop_if_missing(file.path(main_result_dir, "analysis_summary.csv"))
stop_if_missing(file.path(main_result_dir, "cluster_celltype_annotation.csv"))

message("Reading merged Seurat object: ", input_rds)
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj, assay = "RNA")
}

if (!target_gene %in% rownames(obj)) {
  stop(target_gene, " was not found in the Seurat object.")
}
if (!"sample" %in% colnames(obj@meta.data)) {
  stop("The merged object must contain a sample metadata column.")
}

analysis_summary_file <- file.path(main_result_dir, "analysis_summary.csv")
annotation_file <- file.path(main_result_dir, "cluster_celltype_annotation.csv")
nk_clusters <- read_summary_item(analysis_summary_file, "nk_candidate_clusters")

cluster_annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)
cluster_to_manual <- setNames(cluster_annotation$manual_celltype, as.character(cluster_annotation$cluster))
cluster_to_broad <- setNames(cluster_annotation$broad_celltype, as.character(cluster_annotation$cluster))

obj$manual_celltype <- unname(cluster_to_manual[as.character(obj$seurat_clusters)])
obj$broad_celltype <- unname(cluster_to_broad[as.character(obj$seurat_clusters)])
obj$manual_celltype[is.na(obj$manual_celltype) | !nzchar(obj$manual_celltype)] <- "Unassigned"
obj$broad_celltype[is.na(obj$broad_celltype) | !nzchar(obj$broad_celltype)] <- "Unassigned"
obj$manual_celltype <- factor(obj$manual_celltype)
obj$broad_celltype <- factor(obj$broad_celltype)

stk31_expr <- as.numeric(fetch_gene_matrix(obj, target_gene)[target_gene, ])
obj$stk31_positive <- stk31_expr > 0
tumor_epithelial_assignment <- assign_tumor_epithelial_stk31_group(
  obj,
  stk31_expr,
  celltype_label = tumor_epithelial_celltype_label,
  high_quantile = tumor_epithelial_stk31_high_quantile
)
obj$tumor_epithelial_stk31_group <- tumor_epithelial_assignment$group
obj$stk31_group <- obj$tumor_epithelial_stk31_group
obj$nk_group <- ifelse(as.character(obj$seurat_clusters) %in% nk_clusters, "NK_cell", "Other")
obj$relationship_group <- ifelse(
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial",
  ifelse(
    obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial", "STK31_low_tumor_epithelial",
    ifelse(obj$nk_group == "NK_cell", "NK_cell", "Other")
  )
)

message("Summarizing STK31-positive cells by cell type and sample")
celltype_total <- as.data.frame(table(obj$manual_celltype), stringsAsFactors = FALSE)
colnames(celltype_total) <- c("manual_celltype", "total_cells")
celltype_pos <- as.data.frame(table(obj$manual_celltype[obj$stk31_positive]), stringsAsFactors = FALSE)
colnames(celltype_pos) <- c("manual_celltype", "stk31_positive_cells")
stk31_by_celltype <- merge(celltype_total, celltype_pos, by = "manual_celltype", all.x = TRUE)
stk31_by_celltype$stk31_positive_cells[is.na(stk31_by_celltype$stk31_positive_cells)] <- 0
stk31_by_celltype$pct_stk31_positive_within_celltype <- stk31_by_celltype$stk31_positive_cells / stk31_by_celltype$total_cells
stk31_by_celltype$pct_of_all_stk31_positive <- stk31_by_celltype$stk31_positive_cells / sum(stk31_by_celltype$stk31_positive_cells)
stk31_expr_by_celltype <- summarize_gene_by_group(obj, target_gene, obj$manual_celltype, "manual_celltype")
stk31_by_celltype <- merge(stk31_by_celltype, stk31_expr_by_celltype, by = "manual_celltype", all.x = TRUE)
stk31_by_celltype <- stk31_by_celltype[order(-stk31_by_celltype$stk31_positive_cells), ]
write.csv(stk31_by_celltype, file.path(out_dir, "stk31_positive_by_celltype.csv"), row.names = FALSE)

sample_celltype_total <- as.data.frame(table(obj$sample, obj$manual_celltype), stringsAsFactors = FALSE)
colnames(sample_celltype_total) <- c("sample", "manual_celltype", "total_cells")
sample_celltype_pos <- as.data.frame(table(obj$sample[obj$stk31_positive], obj$manual_celltype[obj$stk31_positive]), stringsAsFactors = FALSE)
colnames(sample_celltype_pos) <- c("sample", "manual_celltype", "stk31_positive_cells")
stk31_by_sample_celltype <- merge(sample_celltype_total, sample_celltype_pos, by = c("sample", "manual_celltype"), all.x = TRUE)
stk31_by_sample_celltype$stk31_positive_cells[is.na(stk31_by_sample_celltype$stk31_positive_cells)] <- 0
stk31_by_sample_celltype$pct_stk31_positive_within_sample_celltype <- stk31_by_sample_celltype$stk31_positive_cells / stk31_by_sample_celltype$total_cells
stk31_by_sample_celltype <- stk31_by_sample_celltype[order(stk31_by_sample_celltype$sample, -stk31_by_sample_celltype$stk31_positive_cells), ]
write.csv(stk31_by_sample_celltype, file.path(out_dir, "stk31_positive_by_sample_celltype.csv"), row.names = FALSE)

message("Running tumor-epithelial-only STK31 high vs low comparison")
epithelial_cells <- colnames(obj)[obj$manual_celltype == tumor_epithelial_celltype_label]
epithelial <- subset(obj, cells = epithelial_cells)
# 【按名称对齐】match(x,y) 返回 x 各元素在 y 中的位置，没找到返回 NA；用来保证标签和表达对应同一细胞。
epithelial$tumor_epithelial_stk31_group <- obj$tumor_epithelial_stk31_group[match(colnames(epithelial), colnames(obj))]
epithelial$epithelial_stk31_group <- epithelial$tumor_epithelial_stk31_group

epithelial_summary <- data.frame(
  metric = c("tumor_epithelial_cells", "stk31_high_tumor_epithelial_cells", "stk31_low_tumor_epithelial_cells", "tumor_epithelial_stk31_high_cutoff"),
  value = c(
    ncol(epithelial),
    sum(epithelial$epithelial_stk31_group == "STK31_high_tumor_epithelial"),
    sum(epithelial$epithelial_stk31_group == "STK31_low_tumor_epithelial"),
    tumor_epithelial_assignment$cutoff
  ),
  stringsAsFactors = FALSE
)
write.csv(epithelial_summary, file.path(out_dir, "epithelial_stk31_group_summary.csv"), row.names = FALSE)

epithelial_stk31_by_sample <- as.data.frame(table(epithelial$sample, epithelial$epithelial_stk31_group), stringsAsFactors = FALSE)
colnames(epithelial_stk31_by_sample) <- c("sample", "epithelial_stk31_group", "cells")
write.csv(epithelial_stk31_by_sample, file.path(out_dir, "epithelial_stk31_group_by_sample.csv"), row.names = FALSE)

epi_pos_cells <- colnames(epithelial)[epithelial$epithelial_stk31_group == "STK31_high_tumor_epithelial"]
epi_neg_cells <- colnames(epithelial)[epithelial$epithelial_stk31_group == "STK31_low_tumor_epithelial"]
epithelial_markers <- run_marker_test(
  epithelial,
  cells.1 = epi_pos_cells,
  cells.2 = epi_neg_cells,
  output_file = file.path(out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv")
)
write.csv(head(epithelial_markers, 50), file.path(out_dir, "top_markers_stk31_high_vs_low_tumor_epithelial.csv"), row.names = FALSE)

write_plot(
  file.path(out_dir, "umap_epithelial_stk31_positive_negative.pdf"),
  DimPlot(epithelial, reduction = "umap", group.by = "epithelial_stk31_group") +
    ggtitle("STK31 high vs low tumor epithelial cells"),
  width = 8,
  height = 6
)
write_plot(
  file.path(out_dir, "umap_epithelial_stk31_by_sample.pdf"),
  DimPlot(epithelial, reduction = "umap", group.by = "sample") + ggtitle("Epithelial cells by sample"),
  width = 8,
  height = 6
)
write_plot(
  file.path(out_dir, "epithelial_identity_marker_dotplot.pdf"),
  DotPlot(epithelial, features = available_genes(epithelial, unique(unlist(identity_marker_sets))), group.by = "epithelial_stk31_group") +
    RotatedAxis() + ggtitle("Identity markers in tumor epithelial STK31 high/low groups"),
  width = 14,
  height = 5
)

message("Summarizing per-sample STK31-NK relationship")
samples <- sort(unique(as.character(obj$sample)))
sample_relationship <- do.call(rbind, lapply(samples, function(sample_id) {
  sample_idx <- as.character(obj$sample) == sample_id
  sample_cells <- colnames(obj)[sample_idx]
  nk_idx <- sample_idx & obj$nk_group == "NK_cell"
  tumor_epi_high_idx <- sample_idx & obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"
  tumor_epi_low_idx <- sample_idx & obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"
  epithelial_idx <- sample_idx & obj$manual_celltype == tumor_epithelial_celltype_label

  distance <- NA_real_
  if ("umap" %in% names(obj@reductions) && sum(nk_idx) > 0 && sum(tumor_epi_high_idx) > 0) {
    umap <- Embeddings(obj, "umap")
    nk_center <- colMeans(umap[nk_idx, 1:2, drop = FALSE])
    stk31_center <- colMeans(umap[tumor_epi_high_idx, 1:2, drop = FALSE])
    distance <- sqrt(sum((nk_center - stk31_center)^2))
  }

  data.frame(
    sample = sample_id,
    total_cells = length(sample_cells),
    stk31_positive_cells = sum(sample_idx & obj$stk31_positive),
    pct_stk31_positive = sum(sample_idx & obj$stk31_positive) / length(sample_cells),
    epithelial_cells = sum(epithelial_idx),
    epithelial_stk31_positive_cells = sum(epithelial_idx & obj$stk31_positive),
    pct_epithelial_stk31_positive = ifelse(sum(epithelial_idx) > 0, sum(epithelial_idx & obj$stk31_positive) / sum(epithelial_idx), NA_real_),
    nk_cells = sum(nk_idx),
    pct_nk_cells = sum(nk_idx) / length(sample_cells),
    stk31_high_tumor_epithelial_cells = sum(tumor_epi_high_idx),
    pct_stk31_high_tumor_epithelial_cells = sum(tumor_epi_high_idx) / length(sample_cells),
    stk31_low_tumor_epithelial_cells = sum(tumor_epi_low_idx),
    pct_stk31_low_tumor_epithelial_cells = sum(tumor_epi_low_idx) / length(sample_cells),
    stk31_high_tumor_epithelial_vs_nk_umap_distance = distance,
    stringsAsFactors = FALSE
  )
}))
write.csv(sample_relationship, file.path(out_dir, "sample_level_stk31_nk_relationship.csv"), row.names = FALSE)

message("Scoring NK functional states")
score_columns <- character(0)
for (score_name in names(nk_function_sets)) {
  col_name <- paste0(score_name, "_score")
  obj[[col_name]] <- module_score(obj, nk_function_sets[[score_name]])
  score_columns <- c(score_columns, col_name)
}

nk_cells_all <- colnames(obj)[obj$nk_group == "NK_cell"]
nk_obj <- subset(obj, cells = nk_cells_all)
nk_scores_by_sample <- do.call(rbind, lapply(score_columns, function(col_name) {
  summarize_numeric_by_group(
    values = nk_obj@meta.data[[col_name]],
    group_vec = nk_obj$sample,
    value_name = sub("_score$", "", col_name),
    group_name = "sample"
  )
}))
nk_scores_by_sample <- nk_scores_by_sample[, c("sample", "score_name", "cells", "mean_score", "median_score", "pct_score_positive")]
write.csv(nk_scores_by_sample, file.path(out_dir, "nk_function_scores_by_sample.csv"), row.names = FALSE)

nk_scores_by_cluster <- do.call(rbind, lapply(score_columns, function(col_name) {
  summarize_numeric_by_group(
    values = obj@meta.data[[col_name]],
    group_vec = obj$seurat_clusters,
    value_name = sub("_score$", "", col_name),
    group_name = "cluster"
  )
}))
nk_scores_by_cluster <- nk_scores_by_cluster[, c("cluster", "score_name", "cells", "mean_score", "median_score", "pct_score_positive")]
write.csv(nk_scores_by_cluster, file.path(out_dir, "nk_function_scores_by_cluster.csv"), row.names = FALSE)

score_plot_df <- nk_scores_by_sample
write_plot(
  file.path(out_dir, "nk_function_scores_by_sample.pdf"),
  # 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
  ggplot(score_plot_df, aes(x = sample, y = mean_score, fill = score_name)) +
    geom_col(position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "Mean module score", fill = "Score", title = "NK functional scores by sample"),
  width = 10,
  height = 5
)

write_plot(
  file.path(out_dir, "nk_function_score_violin_by_sample.pdf"),
  VlnPlot(nk_obj, features = score_columns, group.by = "sample", pt.size = 0.01, ncol = 2) +
    ggtitle("NK functional score distributions by sample"),
  width = 12,
  height = 8
)

followup_summary <- data.frame(
  item = c(
    "target_gene",
    "tumor_epithelial_celltype_label",
    "tumor_epithelial_stk31_high_cutoff",
    "nk_clusters_from_main_analysis",
    "samples_analyzed",
    "tumor_epithelial_cells",
    "stk31_high_tumor_epithelial_cells",
    "stk31_low_tumor_epithelial_cells",
    "nk_cells"
  ),
  value = c(
    target_gene,
    tumor_epithelial_celltype_label,
    tumor_epithelial_assignment$cutoff,
    paste(nk_clusters, collapse = ";"),
    paste(samples, collapse = ";"),
    ncol(epithelial),
    sum(epithelial$epithelial_stk31_group == "STK31_high_tumor_epithelial"),
    sum(epithelial$epithelial_stk31_group == "STK31_low_tumor_epithelial"),
    ncol(nk_obj)
  ),
  stringsAsFactors = FALSE
)
write.csv(followup_summary, file.path(out_dir, "followup_summary.csv"), row.names = FALSE)

message("Follow-up analysis done. Key outputs:")
message("- ", file.path(out_dir, "stk31_positive_by_celltype.csv"))
message("- ", file.path(out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv"))
message("- ", file.path(out_dir, "sample_level_stk31_nk_relationship.csv"))
message("- ", file.path(out_dir, "nk_function_scores_by_sample.csv"))

# ============================================================
# 【后续区段 2】进入 CellChat 推断和 GO/候选通路作图；注意这一段重新设置输入输出路径。
# Integrated from 06_merged_cellchat_go_plots.R
# This section used to live in a separate follow-up script.
# ============================================================
# ============================================================
# merged CellChat and GO plot follow-up
# ============================================================
# Purpose:
#   [1] Run CellChat on merged STK31-high epithelial/NK annotated cells
#   [2] Plot cell-cell interaction circle/heatmap/bubble figures
#   [3] Plot GO enrichment heatmaps from existing GO CSV outputs
#   [4] Plot candidate ligand-receptor pathway heatmap/bubble figures
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  if (requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
})

set.seed(20260703)

basic_result_dir <- "/home/zhuweiyu/codex-r/results/merged_basic_seurat"
main_result_dir <- "/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis"
followup_dir <- "/home/zhuweiyu/codex-r/results/merged_stk31_nk_followup"
input_rds <- file.path(basic_result_dir, "gallbladder_cancer_merged_basic_seurat.rds")
out_dir <- file.path(main_result_dir, "cellchat_go_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
min_cells_per_group <- 10
# 【参数 max_cells_per_group】每种细胞群的抽样上限；限制计算规模，并非要求各样本真实细胞比例相同。
max_cells_per_group <- 1500

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)
}

go_files <- c(
  STK31_high_vs_low_tumor_epithelial = file.path(main_result_dir, "stk31_high_vs_low_tumor_epithelial_go_bp.csv"),
  STK31_high_tumor_epithelial_vs_NK = file.path(main_result_dir, "stk31_high_tumor_epithelial_vs_nk_go_bp.csv"),
  STK31_high_tumor_epithelial_vs_all_other = file.path(main_result_dir, "stk31_high_tumor_epithelial_vs_all_other_go_bp.csv"),
  NK_vs_all_other = file.path(main_result_dir, "nk_vs_all_other_go_bp.csv")
)

candidate_lr_file <- file.path(main_result_dir, "candidate_ligand_receptor_pathways.csv")
candidate_pathway_file <- file.path(main_result_dir, "candidate_pathway_summary.csv")
annotation_file <- file.path(main_result_dir, "cluster_celltype_annotation.csv")
analysis_summary_file <- file.path(main_result_dir, "analysis_summary.csv")

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
}

# 【函数：available_genes】把想看的基因与对象实际行名取交集，避免访问不存在的基因。
# 返回基因名向量；返回长度为零表示这些基因不在当前矩阵内。
available_genes <- function(object, genes) {
  intersect(genes, rownames(object))
}

# 【函数：fetch_gene_matrix】从 Seurat 对象提取指定基因的表达矩阵，并兼容不同版本的 layer/slot 接口。
# 返回矩阵通常是基因×细胞；data 是标准化层，counts 是原始计数层。
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

# 【函数：read_summary_item】从 item/value 格式的旧汇总表读取指定项目，并把分号分隔的值拆开。
# 文件、列或项目缺失时按函数分支返回 fallback。
read_summary_item <- function(path, item_name, fallback = character(0)) {
  if (!file.exists(path)) {
    return(fallback)
  }
  summary_df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("item", "value") %in% colnames(summary_df))) {
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

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 【函数：write_base_pdf】打开 PDF 绘图设备，执行 expr 中的绘图表达式，再关闭设备。
# 基础 R 图与 ggplot 的保存方式不同；关闭设备后文件才完整写出。
write_base_pdf <- function(path, expr, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  force(expr)
  dev.off()
}

# 【函数：cap_labels】限制长标签的字符数，避免 GO 名称等文字把图的布局挤坏。
cap_labels <- function(x, max_chars = 60) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

# 【函数：safe_read_csv】在读取 CSV 前检查可用性；缺少文件时按本函数约定返回空结果。
safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# 【函数：plot_go_heatmap】读入多组 GO 富集表，选出代表性条目并拼接热图。
# 颜色编码的是富集统计量，不能据此推断该过程上调或下调。
plot_go_heatmap <- function(go_paths, out_prefix, top_n = 15) {
  go_tables <- lapply(names(go_paths), function(name) {
    x <- safe_read_csv(go_paths[[name]])
    if (is.null(x) || nrow(x) == 0 || !all(c("Description", "p.adjust") %in% colnames(x))) {
      return(NULL)
    }
    if (!"pvalue" %in% colnames(x)) {
      x$pvalue <- x$p.adjust
    }
    if (!"Count" %in% colnames(x)) {
      # 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
      x$Count <- vapply(strsplit(x$GeneRatio, "/"), function(z) as.numeric(z[1]), numeric(1))
    }
    x$comparison <- name
    x <- x[order(x$p.adjust, x$pvalue), ]
    head(x, top_n)
  })
  go_all <- do.call(rbind, go_tables)
  if (is.null(go_all) || nrow(go_all) == 0) {
    warning("No GO rows available for heatmap.")
    return(invisible(NULL))
  }

  go_all$score <- -log10(pmax(go_all$p.adjust, 1e-300))
  write.csv(go_all, paste0(out_prefix, "_go_heatmap_input.csv"), row.names = FALSE)

  terms <- unique(go_all$Description)
  comparisons <- names(go_paths)
  mat <- matrix(0, nrow = length(terms), ncol = length(comparisons), dimnames = list(terms, comparisons))
  for (i in seq_len(nrow(go_all))) {
    mat[go_all$Description[i], go_all$comparison[i]] <- max(mat[go_all$Description[i], go_all$comparison[i]], go_all$score[i])
  }
  mat <- mat[order(rowMeans(mat), decreasing = TRUE), , drop = FALSE]
  rownames(mat) <- cap_labels(rownames(mat), 70)

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(paste0(out_prefix, "_go_enrichment_heatmap.pdf"), width = 10, height = max(6, 0.25 * nrow(mat) + 2))
    pheatmap::pheatmap(
      mat,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      color = colorRampPalette(c("white", "#fee08b", "#f46d43", "#7f0000"))(100),
      main = "GO BP enrichment (-log10 adjusted P)"
    )
    dev.off()
  } else {
    message("Skipping GO heatmap because pheatmap is not installed.")
  }

  go_dot <- go_all[order(go_all$comparison, go_all$p.adjust), ]
  go_dot$Description_short <- factor(cap_labels(go_dot$Description, 70), levels = rev(unique(cap_labels(go_dot$Description, 70))))
  write_plot(
    paste0(out_prefix, "_go_enrichment_dotplot.pdf"),
    ggplot(go_dot, aes(x = comparison, y = Description_short, size = Count, color = score)) +
      geom_point() +
      scale_color_gradient(low = "#4575b4", high = "#d73027") +
      theme_bw() +
      labs(x = "Comparison", y = "GO BP term", color = "-log10(FDR)", size = "Genes", title = "GO BP enrichment") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)),
    width = 10,
    height = max(6, 0.25 * length(unique(go_dot$Description_short)) + 2)
  )
}

# 【函数：plot_candidate_lr】读入人工候选配体受体及通路摘要，绘制气泡图/热图等。
# 这些图的评分来自候选面板，应与 CellChatDB 的正式通讯结果区分。
plot_candidate_lr <- function(lr_file, pathway_file, out_prefix) {
  lr <- safe_read_csv(lr_file)
  pathway <- safe_read_csv(pathway_file)
  if (is.null(lr) || nrow(lr) == 0) {
    warning("No candidate ligand-receptor table found.")
    return(invisible(NULL))
  }
  if (is.null(pathway) || nrow(pathway) == 0) {
    pathway <- aggregate(interaction_score ~ direction + pathway, data = lr, FUN = max)
  }

  pathway <- pathway[order(-pathway$interaction_score), ]
  write.csv(pathway, paste0(out_prefix, "_candidate_pathway_plot_input.csv"), row.names = FALSE)
  top_pathways <- unique(head(pathway$pathway, 20))
  heat_df <- pathway[pathway$pathway %in% top_pathways, ]
  mat <- xtabs(interaction_score ~ pathway + direction, data = heat_df)
  mat <- mat[order(rowMeans(mat), decreasing = TRUE), , drop = FALSE]

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(paste0(out_prefix, "_candidate_pathway_heatmap.pdf"), width = 8, height = max(5, 0.35 * nrow(mat) + 2))
    pheatmap::pheatmap(
      mat,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      color = colorRampPalette(c("white", "#abd9e9", "#2c7bb6", "#08306b"))(100),
      main = "Candidate STK31/NK interaction pathways"
    )
    dev.off()
  } else {
    message("Skipping candidate pathway heatmap because pheatmap is not installed.")
  }

  lr$pair <- paste(lr$ligand, lr$receptor, sep = " -> ")
  lr_top <- head(lr[order(-lr$interaction_score), ], 30)
  lr_top$pair <- factor(lr_top$pair, levels = rev(unique(lr_top$pair)))
  write_plot(
    paste0(out_prefix, "_candidate_lr_bubble.pdf"),
    ggplot(lr_top, aes(x = direction, y = pair, size = interaction_score, color = pathway)) +
      geom_point(alpha = 0.85) +
      theme_bw() +
      labs(x = "Direction", y = "Ligand -> receptor", size = "Score", color = "Pathway", title = "Top candidate ligand-receptor pairs") +
      theme(axis.text.x = element_text(angle = 25, hjust = 1)),
    width = 11,
    height = 8
  )
}

# 【函数：run_cellchat】建立 CellChat 对象，使用人类配体受体库推断细胞群间通讯，并汇总通路和输出。
# 输入细胞标签必须与表达矩阵列一一对应；少量细胞群和表达预处理会影响推断。
# 返回内容依本脚本定义，主流程会继续提取通讯表或保存图；这是计算推断结果。
run_cellchat <- function(object, out_prefix) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    warning("CellChat is not installed; skipping CellChat analysis.")
    return(invisible(NULL))
  }

  data_input <- fetch_gene_matrix(object, rownames(object), slot = "data")
  meta <- data.frame(
    labels = object$cellchat_group,
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )

  # 【通讯输入】表达矩阵列名与 meta 行名必须对应；group.by 指定用哪个标签定义发送和接收细胞群。
  cellchat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  # 【通讯候选筛选】先找群中高表达的基因，再通过配体受体数据库筛选候选相互作用。
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  # 【通讯推断】根据群表达与配体受体库估计通讯分值；raw.use=TRUE 指未投影的表达，不等于原始 counts。
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  # 【通讯过滤】去掉不满足最少细胞数要求的群的通讯，min.cells 是本次阈值。
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_per_group)
  # 【通路汇总】把配体受体层面的通讯聚合到通路层面，再由 aggregateNet 汇总群间网络。
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  saveRDS(cellchat, paste0(out_prefix, "_cellchat_object.rds"))

  # 【通讯表】从 CellChat 对象提取 source、target、ligand、receptor、prob 等字段，便于后续筛选与作图。
  communication <- CellChat::subsetCommunication(cellchat)
  write.csv(communication, paste0(out_prefix, "_cellchat_communications.csv"), row.names = FALSE)

  if (!is.null(cellchat@net$count)) {
    write.csv(cellchat@net$count, paste0(out_prefix, "_cellchat_interaction_count_matrix.csv"))
  }
  if (!is.null(cellchat@net$weight)) {
    write.csv(cellchat@net$weight, paste0(out_prefix, "_cellchat_interaction_weight_matrix.csv"))
  }

  group_size <- as.numeric(table(cellchat@idents))
  names(group_size) <- names(table(cellchat@idents))

  write_base_pdf(
    paste0(out_prefix, "_cellchat_circle_count.pdf"),
    CellChat::netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions"),
    width = 8,
    height = 8
  )
  write_base_pdf(
    paste0(out_prefix, "_cellchat_circle_weight.pdf"),
    CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights"),
    width = 8,
    height = 8
  )
  write_base_pdf(
    paste0(out_prefix, "_cellchat_heatmap_count.pdf"),
    print(CellChat::netVisual_heatmap(cellchat, measure = "count")),
    width = 8,
    height = 7
  )
  write_base_pdf(
    paste0(out_prefix, "_cellchat_heatmap_weight.pdf"),
    print(CellChat::netVisual_heatmap(cellchat, measure = "weight")),
    width = 8,
    height = 7
  )

  focus_groups <- intersect(c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell", "Macrophage", "Mast_cell", "Myeloid"), levels(cellchat@idents))
  if (length(focus_groups) >= 2) {
    write_base_pdf(
      paste0(out_prefix, "_cellchat_bubble_focus_groups.pdf"),
      print(CellChat::netVisual_bubble(cellchat, sources.use = focus_groups, targets.use = focus_groups, remove.isolate = FALSE)),
      width = 12,
      height = 8
    )
  }

  pathway_summary <- data.frame()
  if (length(cellchat@netP$pathways) > 0) {
    pathway_summary <- do.call(rbind, lapply(cellchat@netP$pathways, function(pathway_name) {
      prob <- cellchat@netP$prob[, , pathway_name]
      data.frame(
        pathway = pathway_name,
        total_probability = sum(prob),
        max_probability = max(prob),
        stringsAsFactors = FALSE
      )
    }))
    pathway_summary <- pathway_summary[order(-pathway_summary$total_probability), ]
    write.csv(pathway_summary, paste0(out_prefix, "_cellchat_pathway_summary.csv"), row.names = FALSE)
  }

  invisible(cellchat)
}

message("Loading merged Seurat object")
stop_if_missing(input_rds)
stop_if_missing(annotation_file)
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj, assay = "RNA")
}

cluster_annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)
cluster_to_manual <- setNames(cluster_annotation$manual_celltype, as.character(cluster_annotation$cluster))
obj$manual_celltype <- unname(cluster_to_manual[as.character(obj$seurat_clusters)])
obj$manual_celltype[is.na(obj$manual_celltype) | !nzchar(obj$manual_celltype)] <- "Unassigned"

nk_clusters <- read_summary_item(analysis_summary_file, "nk_candidate_clusters")
stk31_expr <- as.numeric(fetch_gene_matrix(obj, target_gene)[target_gene, ])
tumor_epithelial_assignment <- assign_tumor_epithelial_stk31_group(
  obj,
  stk31_expr,
  celltype_label = tumor_epithelial_celltype_label,
  high_quantile = tumor_epithelial_stk31_high_quantile
)
obj$tumor_epithelial_stk31_group <- tumor_epithelial_assignment$group
obj$nk_group <- ifelse(as.character(obj$seurat_clusters) %in% nk_clusters, "NK_cell", "Other")
obj$cellchat_group <- as.character(obj$manual_celltype)
obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"] <- "STK31_high_tumor_epithelial"
obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"] <- "STK31_low_tumor_epithelial"
obj$cellchat_group[obj$nk_group == "NK_cell"] <- "NK_cell"

group_counts <- sort(table(obj$cellchat_group), decreasing = TRUE)
keep_groups <- names(group_counts[group_counts >= min_cells_per_group])
obj <- subset(obj, cells = colnames(obj)[obj$cellchat_group %in% keep_groups])
obj$cellchat_group <- factor(obj$cellchat_group)
group_counts_before_sampling <- as.data.frame(table(obj$cellchat_group), stringsAsFactors = FALSE)
colnames(group_counts_before_sampling) <- c("cellchat_group", "cells_before_sampling")

cells_by_group <- split(colnames(obj), obj$cellchat_group)
keep_cells <- unlist(lapply(cells_by_group, function(cells) {
  if (length(cells) > max_cells_per_group) {
    # 【抽样】从候选集合抽取元素；replace=TRUE 是有放回抽样，同一个细胞可能重复出现。
    sample(cells, max_cells_per_group)
  } else {
    cells
  }
}), use.names = FALSE)
obj <- subset(obj, cells = keep_cells)
obj$cellchat_group <- droplevels(factor(obj$cellchat_group))

group_counts_after_sampling <- as.data.frame(table(obj$cellchat_group), stringsAsFactors = FALSE)
colnames(group_counts_after_sampling) <- c("cellchat_group", "cells_after_sampling")
group_counts_output <- merge(group_counts_before_sampling, group_counts_after_sampling, by = "cellchat_group", all.x = TRUE)
group_counts_output$cells_after_sampling[is.na(group_counts_output$cells_after_sampling)] <- 0
write.csv(group_counts_output, file.path(out_dir, "cellchat_group_counts.csv"), row.names = FALSE)

message("Plotting GO enrichment heatmaps")
plot_go_heatmap(go_files, file.path(out_dir, "merged"))

message("Plotting candidate ligand-receptor figures")
plot_candidate_lr(candidate_lr_file, candidate_pathway_file, file.path(out_dir, "merged"))

message("Running CellChat")
run_cellchat(obj, file.path(out_dir, "merged"))

message("CellChat/GO plotting done. Outputs written to: ", out_dir)

# ============================================================
# 【后续区段 3】从完整 CellChat 表筛选 STK31-high 上皮↔NK，聚焦 HLA、TGFb、IFNG、TIGIT 等机制轴。
# Integrated from 07_focused_cellchat_stk31_nk_figures.R
# This section used to live in a separate follow-up script.
# ============================================================
# ============================================================
# Focused STK31-high tumor epithelial <-> NK CellChat figures
# ============================================================
# Purpose:
#   [1] Extract CellChat communications between STK31-high tumor epithelial cells and NK cells
#   [2] Export focused tables for MHC-I, TGFb, IFNG/IFN-II, and TIGIT-related axes
#   [3] Generate manuscript-style circle/bubble/heatmap/bar figures for these interactions
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  if (requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
})

cellchat_dir <- file.path("/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis", "cellchat_go_plots")
out_dir <- file.path("/home/zhuweiyu/codex-r/results/merged_stk31_nk_analysis", "cellchat_focused_stk31_nk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comm_file <- file.path(cellchat_dir, "merged_cellchat_communications.csv")
count_file <- file.path(cellchat_dir, "merged_cellchat_interaction_count_matrix.csv")
weight_file <- file.path(cellchat_dir, "merged_cellchat_interaction_weight_matrix.csv")

source_group <- "STK31_high_tumor_epithelial"
target_group <- "NK_cell"

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
}

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 【函数：draw_two_node_circle】把两个细胞群画成节点，用有方向的连线展示二者通讯摘要。
# value_col 指定线条使用的指标；先看图例确认是强度、数量还是其他量。
draw_two_node_circle <- function(path, edge_df, value_col, title, legend_title, color_col = NULL, width = 7, height = 6) {
  edge_df <- edge_df[is.finite(edge_df[[value_col]]) & edge_df[[value_col]] > 0, , drop = FALSE]
  pdf(path, width = width, height = height)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(1, 1, 3, 4), xpd = NA)
  plot.new()
  plot.window(xlim = c(-1.35, 1.85), ylim = c(-1.05, 1.05), asp = 1)

  node_pos <- data.frame(
    group = c("STK31_high_tumor_epithelial", "NK_cell"),
    x = c(-0.65, 0.65),
    y = c(0, 0),
    stringsAsFactors = FALSE
  )
  node_colors <- c(STK31_high_tumor_epithelial = "#d73027", NK_cell = "#4575b4")
  edge_colors <- c(
    "MHC-I / HLA axis" = "#d73027",
    "TGFb axis" = "#1a9850",
    "IFNG / IFN axis" = "#4575b4",
    "TIGIT / NECTIN axis" = "#984ea3",
    "Overall" = "#4d4d4d"
  )

  symbols(node_pos$x, node_pos$y, circles = c(0.18, 0.18), inches = FALSE, add = TRUE, bg = node_colors[node_pos$group], fg = "white", lwd = 2)
  text(node_pos$x, node_pos$y - 0.33, node_pos$group, cex = 0.9)
  title(main = title, cex.main = 1.05)

  if (nrow(edge_df) > 0) {
    max_value <- max(edge_df[[value_col]], na.rm = TRUE)
    for (i in seq_len(nrow(edge_df))) {
      direction <- edge_df$direction[i]
      source <- ifelse(direction == "STK31_high_tumor_epithelial_to_NK_cell", "STK31_high_tumor_epithelial", "NK_cell")
      target <- ifelse(direction == "STK31_high_tumor_epithelial_to_NK_cell", "NK_cell", "STK31_high_tumor_epithelial")
      x0 <- node_pos$x[node_pos$group == source]
      y0 <- node_pos$y[node_pos$group == source]
      x1 <- node_pos$x[node_pos$group == target]
      y1 <- node_pos$y[node_pos$group == target]
      same_direction <- which(edge_df$direction == direction)
      rank_in_direction <- match(i, same_direction)
      curve_span <- seq(0.35, 0.75, length.out = length(same_direction))
      curve_height <- ifelse(source == "STK31_high_tumor_epithelial", curve_span[rank_in_direction], -curve_span[rank_in_direction])
      xs <- seq(x0, x1, length.out = 80)
      ys <- curve_height * sin(seq(0, pi, length.out = 80))
      edge_color <- if (is.null(color_col)) "#4d4d4d" else edge_colors[edge_df[[color_col]][i]]
      if (is.na(edge_color)) edge_color <- "#4d4d4d"
      line_width <- 1.5 + 6 * edge_df[[value_col]][i] / max_value
      lines(xs, ys, col = edge_color, lwd = line_width)
      arrows(xs[72], ys[72], xs[79], ys[79], length = 0.11, angle = 25, col = edge_color, lwd = line_width)
      label_x <- mean(c(x0, x1))
      label_y <- curve_height + ifelse(curve_height > 0, 0.12, -0.12)
      text(label_x, label_y, signif(edge_df[[value_col]][i], 3), col = edge_color, cex = 0.8)
    }
  } else {
    text(0, 0.55, "No interactions detected", cex = 0.9)
  }

  if (!is.null(color_col) && nrow(edge_df) > 0) {
    present_axes <- unique(edge_df[[color_col]])
    present_axes <- present_axes[present_axes %in% names(edge_colors)]
    legend(
      x = 1.08,
      y = 0.55,
      legend = present_axes,
      col = edge_colors[present_axes],
      lwd = 4,
      bty = "n",
      horiz = FALSE,
      cex = 0.72,
      y.intersp = 1.1,
      title = legend_title
    )
  }
}

# 【函数：axis_label】根据通路、配体和受体名称，将相互作用归入预先指定的机制轴标签。
axis_label <- function(pathway, ligand, receptor) {
  if (pathway == "MHC-I" || grepl("^HLA-", ligand)) {
    return("MHC-I / HLA axis")
  }
  if (pathway == "TGFb" || grepl("^TGFB|TGF", ligand) || grepl("TGF", receptor)) {
    return("TGFb axis")
  }
  if (pathway %in% c("IFN-II", "IFN-I") || grepl("^IFNG$|^IFN", ligand) || grepl("IFNGR", receptor)) {
    return("IFNG / IFN axis")
  }
  if (pathway %in% c("TIGIT", "NECTIN") || grepl("TIGIT|NECTIN|PVR", ligand) || grepl("TIGIT|NECTIN|PVR", receptor)) {
    return("TIGIT / NECTIN axis")
  }
  "Other"
}

# 【函数：direction_label】把 source/target 细胞群转换成易读的发送→接收方向名称。
direction_label <- function(source, target) {
  ifelse(
    source == source_group & target == target_group,
    "STK31_high_tumor_epithelial_to_NK_cell",
    ifelse(
      source == target_group & target == source_group,
      "NK_cell_to_STK31_high_tumor_epithelial",
      paste(source, target, sep = "_to_")
    )
  )
}

if (!file.exists(comm_file)) {
  message("Skipping focused CellChat STK31/NK figures because CellChat communications were not generated: ", comm_file)
} else {
comm <- read.csv(comm_file, stringsAsFactors = FALSE, check.names = FALSE)
comm$prob <- as.numeric(comm$prob)
comm$pval <- as.numeric(comm$pval)
comm$direction <- direction_label(comm$source, comm$target)
comm$lr_pair <- paste(comm$ligand, comm$receptor, sep = " -> ")

focused <- comm[
  (comm$source == source_group & comm$target == target_group) |
    (comm$source == target_group & comm$target == source_group),
]
focused <- focused[order(focused$direction, -focused$prob), ]
write.csv(focused, file.path(out_dir, "stk31_high_tumor_epithelial_nk_all_cellchat_communications.csv"), row.names = FALSE)

focused$mechanism_axis <- mapply(axis_label, focused$pathway_name, focused$ligand, focused$receptor)
mechanism <- focused[focused$mechanism_axis != "Other", ]
mechanism <- mechanism[order(mechanism$mechanism_axis, mechanism$direction, -mechanism$prob), ]
write.csv(mechanism, file.path(out_dir, "stk31_high_tumor_epithelial_nk_mechanism_axes_cellchat.csv"), row.names = FALSE)

axis_summary <- aggregate(
  prob ~ mechanism_axis + direction,
  data = mechanism,
  FUN = sum
)
colnames(axis_summary)[colnames(axis_summary) == "prob"] <- "total_probability"
axis_count <- aggregate(
  lr_pair ~ mechanism_axis + direction,
  data = mechanism,
  FUN = length
)
colnames(axis_count)[colnames(axis_count) == "lr_pair"] <- "n_interactions"
axis_summary <- merge(axis_summary, axis_count, by = c("mechanism_axis", "direction"), all = TRUE)
axis_summary <- axis_summary[order(axis_summary$mechanism_axis, axis_summary$direction), ]
write.csv(axis_summary, file.path(out_dir, "stk31_high_tumor_epithelial_nk_mechanism_axis_summary.csv"), row.names = FALSE)

top_mechanism <- mechanism[order(-mechanism$prob), ]
top_mechanism <- top_mechanism[seq_len(min(30, nrow(top_mechanism))), ]
top_mechanism$lr_pair <- factor(top_mechanism$lr_pair, levels = rev(unique(top_mechanism$lr_pair)))
top_mechanism$direction <- factor(top_mechanism$direction, levels = c("STK31_high_tumor_epithelial_to_NK_cell", "NK_cell_to_STK31_high_tumor_epithelial"))
top_mechanism$mechanism_axis <- factor(
  top_mechanism$mechanism_axis,
  levels = c("MHC-I / HLA axis", "TGFb axis", "IFNG / IFN axis", "TIGIT / NECTIN axis")
)

write_plot(
  file.path(out_dir, "stk31_nk_mechanism_axes_bubble.pdf"),
  ggplot(top_mechanism, aes(x = direction, y = lr_pair, size = prob, color = mechanism_axis)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2, 8)) +
    theme_bw() +
    labs(
      x = "Direction",
      y = "Ligand -> receptor",
      color = "Mechanism axis",
      size = "CellChat probability",
      title = "Focused CellChat interactions between STK31-high tumor epithelial cells and NK cells"
    ) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)),
  width = 11,
  height = 7.5
)

if (nrow(axis_summary) > 0 && requireNamespace("pheatmap", quietly = TRUE)) {
  axis_mat <- xtabs(total_probability ~ mechanism_axis + direction, data = axis_summary)
  pdf(file.path(out_dir, "stk31_nk_mechanism_axis_heatmap.pdf"), width = 8, height = 5.5)
  pheatmap::pheatmap(
    axis_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color = colorRampPalette(c("white", "#fee08b", "#f46d43", "#7f0000"))(100),
    main = "STK31-high tumor epithelial / NK mechanism axes"
  )
  dev.off()
} else if (!requireNamespace("pheatmap", quietly = TRUE)) {
  message("Skipping focused mechanism heatmap because pheatmap is not installed.")
}

write_plot(
  file.path(out_dir, "stk31_nk_mechanism_axis_barplot.pdf"),
  ggplot(axis_summary, aes(x = mechanism_axis, y = total_probability, fill = direction)) +
    geom_col(position = "dodge", width = 0.75) +
    theme_bw() +
    labs(
      x = "Mechanism axis",
      y = "Total CellChat probability",
      fill = "Direction",
      title = "Directionality of focused mechanism axes"
    ) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)),
  width = 9,
  height = 5.5
)

draw_two_node_circle(
  file.path(out_dir, "stk31_nk_mechanism_axis_circle.pdf"),
  axis_summary,
  value_col = "total_probability",
  title = "Mechanism-axis circle: STK31-high tumor epithelial cells and NK cells",
  legend_title = "Mechanism axis",
  color_col = "mechanism_axis",
  width = 7.5,
  height = 6.5
)

if (file.exists(count_file) && file.exists(weight_file)) {
  count_matrix <- read.csv(count_file, row.names = 1, check.names = FALSE)
  weight_matrix <- read.csv(weight_file, row.names = 1, check.names = FALSE)
  pair_summary <- data.frame(
    metric = c("interaction_count", "interaction_count", "interaction_weight", "interaction_weight"),
    direction = c("STK31_high_tumor_epithelial_to_NK_cell", "NK_cell_to_STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial_to_NK_cell", "NK_cell_to_STK31_high_tumor_epithelial"),
    value = c(
      count_matrix[source_group, target_group],
      count_matrix[target_group, source_group],
      weight_matrix[source_group, target_group],
      weight_matrix[target_group, source_group]
    ),
    stringsAsFactors = FALSE
  )
  write.csv(pair_summary, file.path(out_dir, "stk31_high_tumor_epithelial_nk_pair_level_cellchat_summary.csv"), row.names = FALSE)

  count_edges <- pair_summary[pair_summary$metric == "interaction_count", , drop = FALSE]
  count_edges$overall_axis <- "Overall"
  draw_two_node_circle(
    file.path(out_dir, "stk31_nk_pair_interaction_count_circle.pdf"),
    count_edges,
    value_col = "value",
    title = "Interaction-count circle: STK31-high tumor epithelial cells and NK cells",
    legend_title = "",
    color_col = "overall_axis",
    width = 7,
    height = 6
  )

  weight_edges <- pair_summary[pair_summary$metric == "interaction_weight", , drop = FALSE]
  weight_edges$overall_axis <- "Overall"
  draw_two_node_circle(
    file.path(out_dir, "stk31_nk_pair_interaction_weight_circle.pdf"),
    weight_edges,
    value_col = "value",
    title = "Interaction-weight circle: STK31-high tumor epithelial cells and NK cells",
    legend_title = "",
    color_col = "overall_axis",
    width = 7,
    height = 6
  )
}

message("Focused CellChat STK31/NK outputs written to: ", out_dir)
}
