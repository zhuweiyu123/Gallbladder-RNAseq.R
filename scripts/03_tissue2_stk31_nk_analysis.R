# ========================================================================
# 【中文阅读指南】tissue2 的 STK31/NK 分析与身份验证
# 输入：脚本 01 保存的 gallbladder_cancer_basic_seurat.rds 和相应结果目录。
# 流程：读取对象 → 标志基因评分/注释 → 上皮内 STK31 分组 → 差异基因 → 候选通讯 → 验证图。
# 输出：results/tissue2_stk31_nk_analysis；身份和功能验证另存于 validation 子目录。
# high 的分位数只在上皮细胞中计算；阈值为零时使用表达 > 0，并不保证恰好四分之一细胞入组。
# 本脚本只有一个样本；细胞层面的差异结果不能当成多个患者之间的重复验证。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# ============================================================
# tissue2 STK31-NK 完整分析脚本
# ============================================================
# 用途：胆囊癌 tissue2 单细胞数据，聚焦 STK31
#   [1] Locate STK31-high tumor epithelial cells
#   [2] 鉴定 NK 细胞候选群（NK marker score）
#   [3] DE: epithelial STK31 high vs low, STK31-high tumor epithelial vs NK/all, NK vs all
#   [4] 候选 ligand-receptor 互作初筛
#   [5] 细胞身份验证（身份 marker DotPlot / UMAP）
#   [6] NK 功能模块评分（细胞毒/活化/耗竭/迁移/TGF-beta）
#   [7] 通路关键基因证据图
#   [8] GO BP 富集（需 clusterProfiler，没装则跳过）
# ============================================================
# 输入：gallbladder_cancer_basic_seurat.rds
# 输出：主分析 -> tissue2_stk31_nk_analysis/
#       验证分析 -> tissue2_stk31_nk_analysis/validation/
# ============================================================
# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# ---------- 配置参数 ----------

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260624)

# 基础分析结果路径
base_result_dir <- "/home/zhuweiyu/codex-r/results/tissue2_basic_seurat"
input_rds <- file.path(base_result_dir, "gallbladder_cancer_basic_seurat.rds")
out_dir <- "/home/zhuweiyu/codex-r/results/tissue2_stk31_nk_analysis"
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 关注基因和判定参数
target_gene <- "STK31"
min_pct_for_stk31_cluster <- 0.01
top_stk31_quantile <- 0.75
tumor_epithelial_celltype_label <- "Epithelial"
# 【参数 tumor_epithelial_stk31_high_quantile】上皮内 high 的分位数设定，0.75 表示第 75 百分位；实际阈值和并列值会影响 high 数量。
tumor_epithelial_stk31_high_quantile <- 0.75
run_exact_wilcox <- TRUE

# 手动指定 NK cluster（留空则自动按 NK marker score 筛选）
manual_nk_clusters <- character(0)

# NK 标志基因列表
nk_marker_genes <- c(
  "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH", "KLRD1", "KLRF1",
  "FCGR3A", "NCAM1", "TYROBP", "CTSW", "CST7", "SPON2", "XCL1", "XCL2"
)

# 粗略细胞类型标志基因：给 cluster 做自动注释提示
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

# 【人工候选表】列出待考察的配体、受体及通路名称；后面用两端表达打分，不等同于自动检验全部通讯。
# Candidate ligand-receptor table for STK31-high tumor epithelial <-> NK screening
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

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, "\nPlease run scripts/gallbladder_cancer_basic_analysis.R first.")
  }
}

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
  if (length(genes) == 0) {
    return(NULL)
  }
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]

  }
}

# 【函数：mean_nonzero】计算 x>0 的比例：TRUE 按 1、FALSE 按 0 求平均。
# 虽然函数名含 mean，但这里返回的是检测比例，并非阳性细胞的平均表达量。
mean_nonzero <- function(x) {
  mean(x > 0)
}

# 统计基因按 cluster 的表达概况

# 【函数：summarize_gene_by_cluster】对指定基因逐聚类统计表达水平和检出比例，返回聚类摘要表。
# 用于定位 STK31 信号集中在哪些群；群均值不能替代每个细胞的表达状态。
summarize_gene_by_cluster <- function(object, gene) {

  expr <- as.numeric(fetch_gene_matrix(object, gene)[gene, ])
  clusters <- as.character(object$seurat_clusters)

  # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
  summary <- do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
    idx <- clusters == cluster_id
    data.frame(
      cluster = cluster_id,
      cells = sum(idx),
      avg_expr = mean(expr[idx]),
      median_expr = median(expr[idx]),
      pct_expr = mean_nonzero(expr[idx]),
      stringsAsFactors = FALSE
    )
  }))
  summary[order(-summary$avg_expr, -summary$pct_expr), ]
}

# 对每个 cluster 计算 NK marker 平均得分

# 【函数：score_nk_clusters】用可检测到的 NK 标志基因给每个聚类计算平均表达分数。
# 返回群级评分用于候选筛选；T 细胞也可能表达部分细胞毒标志，因此需进一步验证。
score_nk_clusters <- function(object, markers) {
  markers <- available_genes(object, markers)
  if (length(markers) == 0) {
    stop("None of the NK marker genes were found in the object.")
  }

  expr <- fetch_gene_matrix(object, markers)
  # 【按细胞求均值】在基因×细胞矩阵中，colMeans 对每列跨基因求平均，常用来构建逐细胞模块分数。
  cell_score <- Matrix::colMeans(expr)
  clusters <- as.character(object$seurat_clusters)

  score_table <- do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
    idx <- clusters == cluster_id
    data.frame(

      cluster = cluster_id,
      cells = sum(idx),
      nk_score = mean(cell_score[idx]),
      nk_marker_pct = mean(cell_score[idx] > 0),
      detected_markers = paste(markers, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  score_table[order(-score_table$nk_score, -score_table$nk_marker_pct), ]
}

# 对每个 cluster 计算大类细胞类型 marker 得分

# 【函数：score_broad_celltypes】按各类标志基因集合计算聚类分数，再选择得分最高的类型作为初步注释。
# 这是依赖人工标志列表的启发式分类；相近分数或混合标志需要核查。
score_broad_celltypes <- function(object, marker_list) {
  clusters <- as.character(object$seurat_clusters)

  scores <- do.call(rbind, lapply(names(marker_list), function(cell_type) {
    markers <- available_genes(object, marker_list[[cell_type]])
    if (length(markers) == 0) {
      return(NULL)
    }

    expr <- fetch_gene_matrix(object, markers)
    cell_score <- Matrix::colMeans(expr)
    do.call(rbind, lapply(sort(unique(clusters)), function(cluster_id) {
      idx <- clusters == cluster_id

      data.frame(
        cluster = cluster_id,
        cell_type = cell_type,
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

# 【旧群级筛选】此处候选聚类用于表达质控，不能替代后续按单细胞表达建立的 STK31 分组。
# QC only: legacy cluster-level STK31 expression screening

# 【函数：choose_stk31_clusters】根据聚类级 STK31 平均表达分位数等条件选出候选聚类。
# 返回群编号；这一步的群筛选与后续上皮细胞内 high/low 分组是两个不同层次。
choose_stk31_clusters <- function(summary_table) {
  cutoff <- as.numeric(stats::quantile(summary_table$avg_expr, top_stk31_quantile, na.rm = TRUE))
  selected <- summary_table$cluster[
    summary_table$avg_expr >= cutoff & summary_table$pct_expr >= min_pct_for_stk31_cluster
  ]
  if (length(selected) == 0) {
    selected <- summary_table$cluster[1]
  }
  selected
}

# 按分位数自动筛选 NK 候选 cluster

# 【函数：choose_nk_clusters】依据 NK 群评分的高分位阈值等条件筛选 NK 候选群。
# 返回候选聚类编号；自动筛选不等于已确认所有入选细胞都是 NK。
choose_nk_clusters <- function(score_table) {
  if (length(manual_nk_clusters) > 0) {
    return(manual_nk_clusters)
  }

  cutoff <- as.numeric(stats::quantile(score_table$nk_score, 0.90, na.rm = TRUE))
  selected <- score_table$cluster[score_table$nk_score >= cutoff & score_table$nk_marker_pct > 0]
  if (length(selected) == 0) {
    selected <- score_table$cluster[1]
  }
  selected
}

# PDF 单图输出

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 差异分析：支持 Wilcoxon 和手工 logFC 两种模式

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
    # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
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

  # 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
  write.csv(markers, output_file, row.names = FALSE)
  markers
}

# 计算给定基因集在指定细胞中的表达

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

# 基于表达量推断候选 ligand-receptor 互作

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

# 【函数：infer_lr_links】沿着候选表逐个读取配体与受体，结合发送细胞和接收细胞的表达摘要打分。
# 这是人工候选配对的描述性筛选；不能把表达乘积当成已证实的分子结合或通讯概率。
infer_lr_links <- function(object, tumor_epithelial_stk31_high_cells, nk_cells, lr_table) {
  lr_genes <- unique(c(lr_table$ligand, lr_table$receptor))
# --- STK31 阳性细胞分布 ---
  tumor_epithelial_high_expr <- average_expression_for_genes(object, tumor_epithelial_stk31_high_cells, lr_genes)
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

# GO BP 富集（需 clusterProfiler，没装则跳过）

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
  converted <- clusterProfiler::bitr(

    sig,

    fromType = "SYMBOL",

    toType = "ENTREZID",

    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  if (nrow(converted) < 10) {
    message("Fewer than 10 converted genes; skip GO enrichment for ", output_prefix)
    return(invisible(NULL))
  }

  # 【GO 过度代表分析】检查候选基因是否集中于某些 GO 条目；BP 是生物过程，BH 用于多重检验调整。
  # 如果没有显式传入 universe，就使用数据库可用背景；富集不直接说明通路被激活。
  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,

    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  write.csv(as.data.frame(ego), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)

}

stop_if_missing(input_rds)

message("Reading Seurat object: ", input_rds)
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

if (!target_gene %in% rownames(obj)) {

  stop(target_gene, " was not found in the Seurat object. Please check gene symbols in features.tsv.gz.")
}

# ============================================================

message("Summarizing STK31 expression by cluster")
# ---------- 主分析开始 ----------
stk31_cluster_summary <- summarize_gene_by_cluster(obj, target_gene)
write.csv(stk31_cluster_summary, file.path(out_dir, "stk31_expression_by_cluster.csv"), row.names = FALSE)
legacy_stk31_high_clusters_qc_only <- choose_stk31_clusters(stk31_cluster_summary)

# ---------- 验证扩展配置 ----------

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

# NK 功能模块评分基因集

nk_function_sets <- list(
  NK_cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "GZMH"),
  NK_activation = c("IFNG", "TNF", "CD69", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5"),
  NK_checkpoint = c("TIGIT", "HAVCR2", "LAG3", "PDCD1", "CTLA4", "TOX"),
  NK_migration = c("CXCR3", "CCR5", "XCR1", "CX3CR1", "CCL3", "CCL4", "CCL5", "XCL1", "XCL2"),
  TGFbeta_response = c("TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "TAGLN", "ACTA2")
)

# 候选通路关键基因：在验证图中单独展示

pathway_panel_genes <- unique(c(
  "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
  "NECTIN2", "PVR", "TIGIT", "HAVCR2", "LAG3", "PDCD1", "CD274",
  "KLRK1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "KIR2DL1", "KIR3DL1",
  "IL18", "IL18R1", "IL15", "IL2RB", "IFNG", "IFNGR1",
  "TNFSF10", "TNFRSF10B", "FASLG", "FAS", "ICAM1", "ITGAL",
  "CXCL9", "CXCL10", "CXCL11", "CCL5", "XCL1", "XCL2"
))

# 全部验证特征基因（上面三组取并集）

validation_feature_panel <- unique(c(
  unlist(identity_marker_sets, use.names = FALSE),
  unlist(nk_function_sets, use.names = FALSE),
# 候选通路关键基因：在验证图中单独展示
  pathway_panel_genes
))

# ---------- 验证扩展工具函数 ----------

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

# 多页 PDF 输出

# 【函数：write_multi_page_pdf】将 plots 中的多个图依次写入同一 PDF，便于连续查看验证结果。
# 每次 print 对应一页，最后关闭设备完成保存。
write_multi_page_pdf <- function(path, plots, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  for (plt in plots) print(plt)
  dev.off()
}

# 取差异分析 top N 基因

# 【函数：top_de_genes】从差异基因结果中选择靠前的候选基因，用于展示；n 控制数量。
# 展示名单只是完整差异表的子集，显著性仍应查看调整后 P 值。
top_de_genes <- function(markers, n = 25) {
  if (is.null(markers) || nrow(markers) == 0) return(data.frame())
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC), markers$gene), ]
  head(markers, n)
}

# GO 富集（可选，没装 clusterProfiler 则跳过）

# 【函数：run_optional_go】可选 GO 分析封装：从差异表选基因、转换 ID，并将富集结果写到指定目录。
# 阅读时留意包检查和基因数检查；跳过不等于生物学上没有相关过程。
run_optional_go <- function(markers, prefix, go_out_dir) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Skipping GO enrichment for ", prefix, ": missing packages")
    return(invisible(NULL))
  }
  if (!all(c("gene", "avg_log2FC") %in% colnames(markers))) return(invisible(NULL))
  if ("p_val_adj" %in% colnames(markers)) {
    sig <- markers$gene[markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25]
  } else {
    sig <- markers$gene[abs(markers$avg_log2FC) > 0.25]
  }
  sig <- unique(sig[!is.na(sig)])
  if (length(sig) < 10) { message("Skipping GO for ", prefix, ": <10 sig genes"); return(invisible(NULL)) }
  converted <- clusterProfiler::bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  if (nrow(converted) < 10) { message("Skipping GO for ", prefix, ": <10 converted"); return(invisible(NULL)) }
  ego <- clusterProfiler::enrichGO(gene = unique(converted$ENTREZID), OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)
  write.csv(as.data.frame(ego), file.path(go_out_dir, paste0(prefix, "_go_bp.csv")), row.names = FALSE)
}

# --- NK cluster 自动筛选 ---

message("Scoring NK-like clusters")
nk_cluster_scores <- score_nk_clusters(obj, nk_marker_genes)
write.csv(nk_cluster_scores, file.path(out_dir, "nk_marker_score_by_cluster.csv"), row.names = FALSE)
nk_clusters <- choose_nk_clusters(nk_cluster_scores)

# --- 自动细胞类型注释 ---

message("Scoring broad cell types for cluster annotation hints")
broad_scores <- score_broad_celltypes(obj, broad_celltype_markers)
write.csv(broad_scores$long, file.path(out_dir, "broad_celltype_scores_by_cluster.csv"), row.names = FALSE)
write.csv(broad_scores$best, file.path(out_dir, "broad_celltype_best_guess_by_cluster.csv"), row.names = FALSE)

cluster_to_celltype <- setNames(as.character(broad_scores$best$cell_type), as.character(broad_scores$best$cluster))
obj$broad_celltype <- unname(cluster_to_celltype[as.character(obj$seurat_clusters)])
obj$broad_celltype[is.na(obj$broad_celltype)] <- "Unassigned"
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
obj$broad_celltype <- as.factor(obj$broad_celltype)

# 将自动注释和手动 NK 注释写回对象

obj$manual_celltype <- as.character(obj$broad_celltype)
obj$manual_celltype[as.character(obj$seurat_clusters) %in% nk_clusters] <- "NK_cell"
# 将自动注释和手动 NK 注释写回对象
obj$manual_celltype <- as.factor(obj$manual_celltype)
obj$cluster_celltype_label <- as.character(obj$manual_celltype)

cluster_celltype_annotation <- unique(data.frame(
  cluster = as.character(obj$seurat_clusters),

  broad_celltype = as.character(obj$broad_celltype),
  manual_celltype = as.character(obj$manual_celltype),
  nk_candidate_by_marker_score = as.character(obj$seurat_clusters) %in% nk_clusters,
  cluster_celltype_label = obj$cluster_celltype_label,
  stringsAsFactors = FALSE
))
cluster_celltype_annotation <- cluster_celltype_annotation[order(as.numeric(cluster_celltype_annotation$cluster)), ]
write.csv(cluster_celltype_annotation, file.path(out_dir, "cluster_celltype_annotation.csv"), row.names = FALSE)
write.csv(cluster_celltype_annotation, file.path(out_dir, "manual_celltype_annotation.csv"), row.names = FALSE)

# 【分析分组】在上皮内区分 STK31-high/low，再单列 NK 与 Other；群名及颜色在后续图表中沿用。
# --- Define analysis groups: STK31-high/low tumor epithelial, NK_cell, Other ---

stk31_expr <- as.numeric(fetch_gene_matrix(obj, target_gene)[target_gene, ])
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
# 【分析分组】在上皮内区分 STK31-high/low，再单列 NK 与 Other；群名及颜色在后续图表中沿用。
# --- Define analysis groups: STK31-high/low tumor epithelial, NK_cell, Other ---
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial",
  ifelse(
    obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial", "STK31_low_tumor_epithelial",
    ifelse(obj$nk_group == "NK_cell", "NK_cell", "Other")
  )
)
obj$relationship_group <- factor(
  obj$relationship_group,
  levels = c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell", "Other")
)
obj$manual_celltype_stk31_tumor_epithelial_split <- as.character(obj$manual_celltype)
obj$manual_celltype_stk31_tumor_epithelial_split[
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"
] <- "STK31_high_tumor_epithelial"
obj$manual_celltype_stk31_tumor_epithelial_split[
  obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"
] <- "STK31_low_tumor_epithelial"
obj$manual_celltype_stk31_tumor_epithelial_split <- factor(obj$manual_celltype_stk31_tumor_epithelial_split)

tumor_epithelial_stk31_high_cells <- colnames(obj)[obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"]
tumor_epithelial_stk31_low_cells <- colnames(obj)[obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"]
nk_cells <- colnames(obj)[obj$nk_group == "NK_cell"]
if (length(tumor_epithelial_stk31_high_cells) == 0) stop("No STK31-high tumor epithelial cells found.")
if (length(tumor_epithelial_stk31_low_cells) == 0) stop("No STK31-low tumor epithelial cells found.")
if (length(nk_cells) == 0) stop("No NK cells found. Check nk_clusters.")
other_cells <- colnames(obj)[obj$relationship_group == "Other"]

# --- STK31 阳性细胞分布 ---

stk31_positive_distribution <- as.data.frame(table(obj$seurat_clusters[stk31_expr > 0]))
colnames(stk31_positive_distribution) <- c("cluster", "stk31_positive_cells")
stk31_positive_distribution$pct_of_all_stk31_positive <- stk31_positive_distribution$stk31_positive_cells / sum(stk31_positive_distribution$stk31_positive_cells)
stk31_positive_distribution <- merge(
  stk31_positive_distribution,
# ---------- 主分析开始 ----------
  stk31_cluster_summary[, c("cluster", "cells", "avg_expr", "pct_expr")],
  by = "cluster",

  all.x = TRUE
)
stk31_positive_distribution <- stk31_positive_distribution[order(-stk31_positive_distribution$stk31_positive_cells), ]
write.csv(stk31_positive_distribution, file.path(out_dir, "stk31_positive_cell_distribution.csv"), row.names = FALSE)

# 输出分析总览表

summary_table <- data.frame(
  item = c(
    "target_gene", "tumor_epithelial_celltype_label", "tumor_epithelial_stk31_high_cutoff",
    "tumor_epithelial_cells", "stk31_high_tumor_epithelial_cells", "stk31_low_tumor_epithelial_cells",
    "legacy_stk31_high_clusters_qc_only", "nk_candidate_clusters", "nk_candidate_cells", "other_cells",
    "detected_nk_markers"
  ),

  value = c(

# 关注基因和判定参数

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
    paste(available_genes(obj, nk_marker_genes), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(out_dir, "analysis_summary.csv"), row.names = FALSE)

# ---------- 生成可视化 ----------

message("Writing STK31 and NK visualization")

write_plot(

  file.path(out_dir, "umap_broad_celltype_annotation.pdf"),

  # 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
  # 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
  DimPlot(obj, reduction = "umap", group.by = "broad_celltype", label = TRUE, repel = TRUE) +
    ggtitle("Marker-based broad cell type annotation")
)
write_plot(
  file.path(out_dir, "umap_manual_celltype_annotation.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "manual_celltype", label = TRUE, repel = TRUE) +
    ggtitle("Manual annotation with NK marker score candidates")
)
write_plot(
  file.path(out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "manual_celltype_stk31_tumor_epithelial_split", label = TRUE, repel = TRUE) +
    ggtitle("Manual annotation with STK31 high/low tumor epithelial split")
)
write_plot(
  file.path(out_dir, "umap_cluster_celltype_labels.pdf"),

  DimPlot(obj, reduction = "umap", group.by = "cluster_celltype_label", label = TRUE, repel = TRUE) +
    ggtitle("Cluster labels with manual cell type annotation"),
  width = 10,
  height = 7

)
write_plot(
  file.path(out_dir, "umap_stk31_expression.pdf"),
  FeaturePlot(obj, features = target_gene, reduction = "umap", order = TRUE) + ggtitle("STK31 expression")
)
write_plot(

  file.path(out_dir, "umap_stk31_high_and_nk_groups.pdf"),

  DimPlot(obj, reduction = "umap", group.by = "relationship_group", label = FALSE) +

    ggtitle("STK31-high tumor epithelial cells and NK candidate cells")

)
write_plot(
  file.path(out_dir, "stk31_violin_by_cluster.pdf"),
  VlnPlot(obj, features = target_gene, group.by = "seurat_clusters", pt.size = 0.01) +
    ggtitle("STK31 by cluster"),
  width = 10,

  height = 5
)
write_plot(
  file.path(out_dir, "nk_marker_dotplot.pdf"),
  DotPlot(obj, features = available_genes(obj, nk_marker_genes), group.by = "seurat_clusters") +
    RotatedAxis() + ggtitle("NK marker expression by cluster"),
  width = 12,
  height = 6
)
write_plot(

  file.path(out_dir, "broad_celltype_marker_dotplot.pdf"),
  DotPlot(obj, features = unique(unlist(lapply(broad_celltype_markers, function(x) available_genes(obj, x)))), group.by = "seurat_clusters") +
    RotatedAxis() + ggtitle("Broad cell type markers by cluster"),
  width = 16,
  height = 7
)

if (!"umap" %in% names(obj@reductions)) {
  stop("UMAP reduction not found. Please run RunUMAP first.")
}
# 【二维中心距离】比较两组 UMAP 坐标均值的位置；这是可视化描述，不证明细胞真实相邻或发生通讯。
# --- UMAP centroid distance: STK31-high tumor epithelial vs NK_cell ---
umap <- Embeddings(obj, "umap")

centroids <- aggregate(
  umap,
  by = list(group = obj$relationship_group),
  FUN = mean
)
write.csv(centroids, file.path(out_dir, "relationship_group_umap_centroids.csv"), row.names = FALSE)
if (all(c("STK31_high_tumor_epithelial", "NK_cell") %in% centroids$group)) {
  # 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
  umap_cols <- setdiff(colnames(centroids), "group")[1:2]
  stk31_center <- as.numeric(centroids[centroids$group == "STK31_high_tumor_epithelial", umap_cols])
  nk_center <- as.numeric(centroids[centroids$group == "NK_cell", umap_cols])
  distance_table <- data.frame(
    comparison = "STK31_high_tumor_epithelial_vs_NK_cell",
    umap_centroid_distance = sqrt(sum((stk31_center - nk_center)^2))

  )
  write.csv(distance_table, file.path(out_dir, "stk31_nk_umap_centroid_distance.csv"), row.names = FALSE)
}

# ============================================================

message("Running differential expression tests")
# ---------- 差异分析（Wilcoxon 秩和检验）----------
tumor_epithelial_high_vs_low <- run_marker_test(
  obj,
  cells.1 = tumor_epithelial_stk31_high_cells,
  cells.2 = tumor_epithelial_stk31_low_cells,
  output_file = file.path(out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv")
)
tumor_epithelial_high_vs_nk <- run_marker_test(
  obj,
  cells.1 = tumor_epithelial_stk31_high_cells,
  cells.2 = nk_cells,
  output_file = file.path(out_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv")
)
tumor_epithelial_high_vs_other <- run_marker_test(
  obj,

  cells.1 = tumor_epithelial_stk31_high_cells,
  cells.2 = c(nk_cells, tumor_epithelial_stk31_low_cells, other_cells),
  output_file = file.path(out_dir, "markers_stk31_high_tumor_epithelial_vs_all_other.csv")
)
nk_vs_other <- run_marker_test(
  obj,
  cells.1 = nk_cells,
  cells.2 = c(tumor_epithelial_stk31_high_cells, tumor_epithelial_stk31_low_cells, other_cells),
  output_file = file.path(out_dir, "markers_nk_vs_all_other.csv")
)

# ---------- 候选 ligand-receptor 互作 ----------

message("Inferring candidate ligand-receptor pathways")
lr_links <- infer_lr_links(obj, tumor_epithelial_stk31_high_cells, nk_cells, lr_reference)
write.csv(lr_links, file.path(out_dir, "candidate_ligand_receptor_pathways.csv"), row.names = FALSE)

pathway_summary <- aggregate(

  interaction_score ~ direction + pathway,
  data = lr_links,
  FUN = max
)
pathway_summary <- pathway_summary[order(-pathway_summary$interaction_score), ]
write.csv(pathway_summary, file.path(out_dir, "candidate_pathway_summary.csv"), row.names = FALSE)

# --- GO 富集（可选）---

message("Trying optional GO enrichment")
# GO BP 富集（需 clusterProfiler，没装则跳过）
run_go_if_available(tumor_epithelial_high_vs_low, file.path(out_dir, "stk31_high_vs_low_tumor_epithelial"))
# GO BP 富集（需 clusterProfiler，没装则跳过）
run_go_if_available(tumor_epithelial_high_vs_nk, file.path(out_dir, "stk31_high_tumor_epithelial_vs_nk"))
run_go_if_available(tumor_epithelial_high_vs_other, file.path(out_dir, "stk31_high_tumor_epithelial_vs_all_other"))
# GO BP 富集（需 clusterProfiler，没装则跳过）
run_go_if_available(nk_vs_other, file.path(out_dir, "nk_vs_all_other"))

# 保存主分析 Seurat 对象

# 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
saveRDS(obj, file.path(out_dir, "tissue2_stk31_nk_annotated_seurat.rds"))

# ============================================================

# 【身份核查】通过 marker 检查聚类标签是否合理；marker 列表不用于重新定义 STK31-high。
# --- Cluster markers for identity QC only, not for STK31-high grouping ---
message("Finding all cluster markers (this may take a while)")
# 【聚类标志】逐群与其他细胞比较，寻找解释细胞身份的 marker；仍需要检查标志组合及样本组成。
all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  all_markers,
  file.path(out_dir, "all_cluster_markers.csv"),
  row.names = FALSE
)

message("Done. Key outputs:")
message("- ", file.path(out_dir, "all_cluster_markers.csv"))
message("- ", file.path(out_dir, "stk31_expression_by_cluster.csv"))
message("- ", file.path(out_dir, "nk_marker_score_by_cluster.csv"))
message("- ", file.path(out_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv"))
message("- ", file.path(out_dir, "candidate_ligand_receptor_pathways.csv"))

# 结果输出到主目录下的 validation/ 子目录

validation_dir <- file.path(out_dir, "validation")
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 验证扩展开始 ----------

message("Starting validation analysis")

# 定义验证细胞类型：保留手动注释，NK cluster 标为 NK_cell

obj$validation_celltype <- as.character(obj$manual_celltype)
obj$validation_celltype[is.na(obj$validation_celltype) | !nzchar(obj$validation_celltype)] <- "Unassigned"
obj$validation_celltype[as.character(obj$seurat_clusters) %in% nk_clusters] <- "NK_cell"
obj$validation_celltype <- factor(obj$validation_celltype)
obj$validation_celltype_stk31_tumor_epithelial_split <- as.character(obj$validation_celltype)
obj$validation_celltype_stk31_tumor_epithelial_split[
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"
] <- "STK31_high_tumor_epithelial"
obj$validation_celltype_stk31_tumor_epithelial_split[
  obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"
] <- "STK31_low_tumor_epithelial"
obj$validation_celltype_stk31_tumor_epithelial_split <- factor(obj$validation_celltype_stk31_tumor_epithelial_split)

obj$validation_group <- ifelse(
# 【验证分组】把前面确定的分组用于身份与功能验证图；先核对每组细胞数及标志表达。
# --- Define validation groups: STK31-high/low tumor epithelial, NK_cell, Other ---
  obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial", "STK31_high_tumor_epithelial",
  ifelse(
    obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial", "STK31_low_tumor_epithelial",
    ifelse(obj$nk_group == "NK_cell", "NK_cell", "Other")
  )
)
obj$validation_group <- factor(
  obj$validation_group,
  levels = c("STK31_high_tumor_epithelial", "STK31_low_tumor_epithelial", "NK_cell", "Other")
)

# 输出 cluster 细胞注释表

cluster_annotation <- unique(data.frame(
  cluster = as.character(obj$seurat_clusters),
  validation_celltype = as.character(obj$validation_celltype),
  tumor_epithelial_stk31_group = as.character(obj$tumor_epithelial_stk31_group),
  nk_group = as.character(obj$nk_group),
  stringsAsFactors = FALSE
))
cluster_annotation <- cluster_annotation[order(as.numeric(cluster_annotation$cluster)), ]
write.csv(cluster_annotation, file.path(validation_dir, "cluster_validation_annotation.csv"), row.names = FALSE)

celltype_counts <- as.data.frame(table(obj$validation_celltype), stringsAsFactors = FALSE)
colnames(celltype_counts) <- c("validation_celltype", "cell_count")
celltype_counts$pct_cells <- celltype_counts$cell_count / sum(celltype_counts$cell_count)

write.csv(celltype_counts, file.path(validation_dir, "validation_celltype_counts.csv"), row.names = FALSE)

# --- 验证 UMAP ---

message("Writing validation UMAPs")
write_plot(
  file.path(validation_dir, "umap_validation_celltypes.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "validation_celltype", label = TRUE, repel = TRUE) +

    ggtitle("Validation cell types")
)
write_plot(
  file.path(validation_dir, "umap_validation_celltypes_with_stk31_tumor_epithelial_split.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "validation_celltype_stk31_tumor_epithelial_split", label = TRUE, repel = TRUE) +
    ggtitle("Validation cell types with STK31 high/low tumor epithelial split")
)
write_plot(
  file.path(validation_dir, "umap_validation_groups.pdf"),
  DimPlot(obj, reduction = "umap", group.by = "validation_group", label = TRUE, repel = TRUE) +
    ggtitle("STK31-high tumor epithelial and NK_cell groups")
)

# --- 批量验证 FeaturePlot ---

validation_feature_genes <- available_genes(obj, validation_feature_panel)
if (length(validation_feature_genes) > 0) {
  feature_plots <- lapply(validation_feature_genes, function(gene) {
    FeaturePlot(obj, features = gene, reduction = "umap", order = TRUE) + ggtitle(gene)
  })
# 多页 PDF 输出
  write_multi_page_pdf(
    file.path(validation_dir, "validation_feature_umaps.pdf"),
    feature_plots, width = 7, height = 6
  )
}

# --- 身份 marker DotPlot ---

identity_genes <- unique(unlist(identity_marker_sets, use.names = FALSE))
# --- 身份 marker DotPlot ---
identity_genes <- available_genes(obj, identity_genes)
if (length(identity_genes) > 0) {
  write_plot(
    file.path(validation_dir, "identity_marker_dotplot_by_celltype.pdf"),

    DotPlot(obj, features = identity_genes, group.by = "validation_celltype") +
      RotatedAxis() + ggtitle("Identity markers by validation cell type"),
    width = 16, height = 7
  )
  identity_table <- summarize_gene_by_group(obj, identity_genes, obj$validation_celltype, group_name = "validation_celltype")
  if (nrow(identity_table) > 0) {
    write.csv(identity_table, file.path(validation_dir, "identity_marker_expression_by_celltype.csv"), row.names = FALSE)
  }
}

# --- NK cluster 自动筛选 ---

message("Scoring NK function signatures")
score_names <- character(0)
for (score_name in names(nk_function_sets)) {
  genes <- available_genes(obj, nk_function_sets[[score_name]])
  if (length(genes) < 2) next
  obj <- AddModuleScore(obj, features = list(genes), name = score_name, assay = DefaultAssay(obj))
  score_names <- c(score_names, paste0(score_name, "1"))
}
if (length(score_names) > 0) {

  score_rename <- setNames(names(nk_function_sets)[seq_along(score_names)], score_names)
  for (old_name in names(score_rename)) {
    obj[[score_rename[[old_name]]]] <- obj[[old_name]]
  }
}

# --- NK 评分可视化 ---

nk_score_columns <- intersect(names(nk_function_sets), colnames(obj@meta.data))
if (length(nk_score_columns) > 0) {
  score_feature_plots <- lapply(nk_score_columns, function(sn) {
    FeaturePlot(obj, features = sn, reduction = "umap", order = TRUE) + ggtitle(sn)
  })
# 多页 PDF 输出
  write_multi_page_pdf(
    file.path(validation_dir, "nk_signature_umaps.pdf"),
    score_feature_plots, width = 7, height = 6
  )
  score_violin_plots <- lapply(nk_score_columns, function(sn) {
    VlnPlot(obj, features = sn, group.by = "validation_celltype", pt.size = 0.01) +
      RotatedAxis() + ggtitle(sn)
  })
# 多页 PDF 输出
  write_multi_page_pdf(
    file.path(validation_dir, "nk_signature_violin_by_celltype.pdf"),
    score_violin_plots, width = 10, height = 5
  )
  nk_score_summary <- do.call(rbind, lapply(nk_score_columns, function(sn) {
    data.frame(
      validation_celltype = as.character(obj$validation_celltype),
      score_name = sn,
      score_value = obj@meta.data[[sn]],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(nk_score_summary, file.path(validation_dir, "nk_signature_scores_by_cell.csv"), row.names = FALSE)
  nk_by_type <- aggregate(score_value ~ validation_celltype + score_name, data = nk_score_summary, FUN = mean)
  write.csv(nk_by_type, file.path(validation_dir, "nk_signature_scores_by_celltype.csv"), row.names = FALSE)
}

# --- 通路基因 DotPlot ---

message("Writing pathway evidence plots")
pathway_genes <- available_genes(obj, pathway_panel_genes)
if (length(pathway_genes) > 0) {
  write_plot(
    file.path(validation_dir, "pathway_gene_dotplot_by_celltype.pdf"),
    DotPlot(obj, features = pathway_genes, group.by = "validation_celltype") +
      RotatedAxis() + ggtitle("Candidate pathway genes by validation cell type"),
    width = 18, height = 8
  )
  pw_expr <- summarize_gene_by_group(obj, pathway_genes, obj$validation_celltype, group_name = "validation_celltype")
  write.csv(pw_expr, file.path(validation_dir, "pathway_gene_expression_by_celltype.csv"), row.names = FALSE)
}

# 【函数：copy_if_exists】只有源文件存在时才复制到目标位置，便于汇集已成功生成的结果。
copy_if_exists <- function(src, dst) {
  # 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
  if (file.exists(src)) tryCatch(file.copy(src, dst, overwrite = TRUE), error = function(e) NULL)
}
copy_if_exists(file.path(out_dir, "candidate_pathway_summary.csv"), file.path(validation_dir, "candidate_pathway_summary.csv"))
copy_if_exists(file.path(out_dir, "candidate_ligand_receptor_pathways.csv"), file.path(validation_dir, "candidate_ligand_receptor_pathways.csv"))

lr_file <- file.path(out_dir, "candidate_ligand_receptor_pathways.csv")
if (file.exists(lr_file)) {
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  candidate_lr <- read.csv(lr_file, stringsAsFactors = FALSE)
  top_lr <- candidate_lr[order(-candidate_lr$interaction_score), ]
  top_lr <- head(top_lr, min(20, nrow(top_lr)))
  if (nrow(top_lr) > 0) {
    write.csv(top_lr, file.path(validation_dir, "top_candidate_interactions.csv"), row.names = FALSE)
  }
}

# --- 差异基因 top 40 + GO 富集 ---

message("Exporting top DE genes and running optional GO enrichment")
if (!is.null(tumor_epithelial_high_vs_low) && nrow(tumor_epithelial_high_vs_low) > 0) {
  write.csv(top_de_genes(tumor_epithelial_high_vs_low, 40), file.path(validation_dir, "top_stk31_high_vs_low_tumor_epithelial_genes.csv"), row.names = FALSE)
# GO 富集（可选，没装 clusterProfiler 则跳过）
  run_optional_go(tumor_epithelial_high_vs_low, "stk31_high_vs_low_tumor_epithelial", validation_dir)
}
if (!is.null(tumor_epithelial_high_vs_nk) && nrow(tumor_epithelial_high_vs_nk) > 0) {
  write.csv(top_de_genes(tumor_epithelial_high_vs_nk, 40), file.path(validation_dir, "top_stk31_high_tumor_epithelial_vs_nk_genes.csv"), row.names = FALSE)
# GO 富集（可选，没装 clusterProfiler 则跳过）
  run_optional_go(tumor_epithelial_high_vs_nk, "stk31_high_tumor_epithelial_vs_nk", validation_dir)
}
if (!is.null(tumor_epithelial_high_vs_other) && nrow(tumor_epithelial_high_vs_other) > 0) {
  write.csv(top_de_genes(tumor_epithelial_high_vs_other, 40), file.path(validation_dir, "top_stk31_high_tumor_epithelial_vs_all_genes.csv"), row.names = FALSE)
  run_optional_go(tumor_epithelial_high_vs_other, "stk31_high_tumor_epithelial_vs_all_other", validation_dir)
}
if (!is.null(nk_vs_other) && nrow(nk_vs_other) > 0) {
  write.csv(top_de_genes(nk_vs_other, 40), file.path(validation_dir, "top_nk_vs_all_genes.csv"), row.names = FALSE)
# GO 富集（可选，没装 clusterProfiler 则跳过）
  run_optional_go(nk_vs_other, "nk_vs_all_other", validation_dir)
}

# --- 验证分析总览 ---

val_summary <- data.frame(
  item = c("target_gene", "tumor_epithelial_celltype_label", "tumor_epithelial_stk31_high_cutoff",
           "tumor_epithelial_cells", "stk31_high_tumor_epithelial_cells", "stk31_low_tumor_epithelial_cells",
           "legacy_stk31_high_clusters_qc_only", "nk_candidate_clusters", "nk_candidate_cells",
           "validation_groups", "validation_celltypes"),
  value = c(
# 关注基因和判定参数
    target_gene,
    tumor_epithelial_celltype_label,
    tumor_epithelial_assignment$cutoff,
    tumor_epithelial_assignment$tumor_epithelial_cells,
    sum(obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"),
    sum(obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"),
    paste(legacy_stk31_high_clusters_qc_only, collapse = ";"),
    paste(nk_clusters, collapse = ";"),
    sum(obj$nk_group == "NK_cell"),
    paste(levels(obj$validation_group), collapse = ";"),
    paste(levels(obj$validation_celltype), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(val_summary, file.path(validation_dir, "validation_summary.csv"), row.names = FALSE)

# 保存验证 Seurat 对象

saveRDS(obj, file.path(validation_dir, "tissue2_stk31_nk_validated_seurat.rds"))

message("Validation done. Outputs in: ", validation_dir)
