# ========================================================================
# 【中文阅读指南】细化 T/NK 注释，再重跑 STK31 与通讯分析
# 输入：合并基础对象、旧细胞类型注释表和旧分析汇总表。
# 阶段 1：提取 T/NK 候选、重新聚类，以 T、NK、排除标志基因细化注释并质控。
# 阶段 2：冻结新标签，重建 STK31 分组，计算差异基因、候选配体受体和 NK 功能分数。
# 阶段 3：CellChat 比较 STK31-high/low 上皮与高置信 NK 之间的双向通讯。
# 输出分别在 merged_tnk_refined_annotation、merged_stk31_nk_refined_analysis、merged_cellchat_refined_high_low_nk。
# STK31_STAGE_START 默认为 1；选择 2 或 3 需要前面阶段的对象和检查结果已存在。
# NK_cell 是最终高置信分组；T_NK_ambiguous 保留身份不确定的细胞，不能直接并入 NK。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
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

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
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
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(stage1_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stage2_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stage3_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
tumor_epithelial_celltype_label <- "Epithelial"
# 【参数 tumor_epithelial_stk31_high_quantile】上皮内 high 的分位数设定，0.75 表示第 75 百分位；实际阈值和并列值会影响 high 数量。
tumor_epithelial_stk31_high_quantile <- 0.75
# 【参数 min_cells_for_de】差异分析两组需达到的最低细胞数；过小的群会停止或跳过检验。
min_cells_for_de <- 10
# 【参数 min_cells_per_cellchat_group】每个通讯群需达到的最低细胞数；不是患者样本量。
min_cells_per_cellchat_group <- 10
# 【参数 max_cells_per_cellchat_group】通讯分析每群的抽样上限，用于控制内存与耗时；实际群大小另有导出表。
max_cells_per_cellchat_group <- 1500
run_exact_wilcox <- TRUE
# Speed control for large DE (without presto). Descriptive only.
# 【参数 max_cells_per_ident_de】每组用于差异检验的最大细胞数；抽样会影响精确数值，需保留随机种子。
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
max_cells_per_ident_de <- as.integer(Sys.getenv("STK31_MAX_CELLS_DE", "500"))
# Resume control: "1" full run; "2" skip stage1 if QC passed; "3" stage3 only
# 【断点续跑】1 从注释开始，2/3 从后续阶段开始；后续阶段仍需要前面保存的对象和质控依据。
stage_start <- as.integer(Sys.getenv("STK31_STAGE_START", "1"))

# Subclustering
# 【参数 tnk_resolution】T/NK 子集图聚类分辨率，越大通常得到更多子群。
tnk_resolution <- 0.4
# 【参数 tnk_n_variable】T/NK 子集选择的高变基因数，用于后续降维。
tnk_n_variable <- 2000
# 【参数 tnk_n_pcs】在 T/NK 子集中计算的 PCA 主成分数量。
tnk_n_pcs <- 30
# 【参数 tnk_use_dims】从已计算的 PCA 中实际用于邻居图/UMAP 的维度范围。
tnk_use_dims <- 1:20
# 【参数 nk_score_cluster_quantile】按群 NK 分数分位数扩充候选的阈值；值越大筛选通常越严格。
nk_score_cluster_quantile <- 0.85
# 【参数 score_delta】T/NK 分数差达到这一量级才认为有明确偏好；这是注释规则参数。
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
# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

# 【函数：available_genes】把想看的基因与对象实际行名取交集，避免访问不存在的基因。
# 返回基因名向量；返回长度为零表示这些基因不在当前矩阵内。
# 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
available_genes <- function(object, genes) intersect(genes, rownames(object))

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

# 【函数：module_score】对可用基因的标准化表达逐细胞取平均，得到基因集合的简易表达分数。
# 本函数没有 AddModuleScore 的匹配背景扣除；分数不能直接视为功能活性的实测值。
module_score <- function(object, genes) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(rep(NA_real_, ncol(object)))
  # 【按细胞求均值】在基因×细胞矩阵中，colMeans 对每列跨基因求平均，常用来构建逐细胞模块分数。
  as.numeric(Matrix::colMeans(fetch_gene_matrix(object, genes)))
}

# 【函数：gene_mean】对指定基因集合、指定细胞的表达取平均，得到一个汇总数值。
# 没有可用基因或细胞时返回 NA，表示无法计算。
gene_mean <- function(object, genes, cells = NULL) {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NA_real_)
  mat <- fetch_gene_matrix(object, genes)
  if (!is.null(cells)) mat <- mat[, intersect(cells, colnames(mat)), drop = FALSE]
  if (ncol(mat) == 0) return(NA_real_)
  mean(as.numeric(Matrix::colMeans(mat)))
}

# 【函数：read_summary_item】从 item/value 格式的旧汇总表读取指定项目，并把分号分隔的值拆开。
# 文件、列或项目缺失时按函数分支返回 fallback。
read_summary_item <- function(path, item_name, fallback = character(0)) {
  if (!file.exists(path)) return(fallback)
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  summary_df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("item", "value") %in% colnames(summary_df))) return(fallback)
  idx <- which(summary_df$item == item_name)
  if (length(idx) == 0) return(fallback)
  value <- summary_df$value[idx[1]]
  if (is.na(value) || !nzchar(as.character(value))) return(fallback)
  unlist(strsplit(as.character(value), ";", fixed = TRUE))
}

# 【函数：assign_tumor_epithelial_stk31_group】只在指定上皮细胞内计算 STK31 分位数阈值，随后生成 high/low 标签。
# 阈值>0 时使用 >= 阈值；阈值为零时使用 >0；其余细胞保持 Other。
# 返回 list，含每个细胞的标签、实际阈值与各组数量；先看这些数量再解读下游差异。
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

# 【函数：run_marker_test】输入对象和两组细胞，调用差异表达检验并把基因表写入 output_file。
# 重点看比较方向、avg_log2FC 和调整后 P 值；正 logFC 表示第一组相对第二组更高。
run_marker_test <- function(object, cells.1, cells.2, output_file) {
  if (length(cells.1) < min_cells_for_de || length(cells.2) < min_cells_for_de) {
    warning("Too few cells for marker test: ", output_file)
    empty <- data.frame()
    # 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
    write.csv(empty, output_file, row.names = FALSE)
    return(empty)
  }
  object$.comparison_group <- "unused"
  object$.comparison_group[colnames(object) %in% cells.1] <- "group_1"
  object$.comparison_group[colnames(object) %in% cells.2] <- "group_2"
  # 【差异表达】比较 ident.1 与 ident.2；avg_log2FC>0 表示第一组更高，p_val_adj 是多重检验调整后 P 值。
  # min.pct/logfc.threshold 是进入检验的筛选条件；细胞级检验不自动控制患者内相关性。
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
    # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
    avg_expr = as.numeric(Matrix::rowMeans(expr)),
    pct_expr = as.numeric(Matrix::rowMeans(expr > 0)),
    stringsAsFactors = FALSE
  )
}

# 【函数：infer_lr_links】沿着候选表逐个读取配体与受体，结合发送细胞和接收细胞的表达摘要打分。
# 这是人工候选配对的描述性筛选；不能把表达乘积当成已证实的分子结合或通讯概率。
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

# 【函数：label_evidence】比较 high/low 两个通讯值，按最小强度和比值阈值标记更高、相近或未检出。
# 默认比值阈值 1.25 是规则阈值，并不等于显著性检验。
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

# 【函数：write_session_info】保存 R/包版本及本次关键配置，方便追溯结果与重现运行环境。
write_session_info <- function(path) {
  sink(path)
  cat("timestamp:", as.character(Sys.time()), "\n")
  cat("script: 06_refine_t_nk_and_rerun_stk31_nk.R\n")
  cat("project_dir:", project_dir, "\n")
  cat("input_rds:", input_rds, "\n")
  cat("tnk_resolution:", tnk_resolution, "\n")
  cat("score_delta:", score_delta, "\n")
  cat("tumor_epithelial_stk31_high_quantile:", tumor_epithelial_stk31_high_quantile, "\n")
  # 【运行记录】记录 R 与已加载包的版本，帮助以后解释同一代码为何可能得到不同结果。
  print(sessionInfo())
  sink()
}

# 【函数：list_output_files】列出指定结果目录的文件，形成运行产物清单。
list_output_files <- function(dirs) {
  # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
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
  # 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
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
    # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
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
# 【阶段 1】在 T/NK 候选子集中重新做表达处理和聚类，先把细胞身份分清，再进入下游比较。
# STAGE 1: strict T/NK re-annotation
# ============================================================
message("=== STAGE 1: T/NK refinement ===")

# 【扩充候选】检查所有原聚类的 NK 标志分数，避免仅局限于原先手工选的 0/1/5 群。
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

# 【取细胞子集】按 cells 或条件保留需要的细胞；这一步改变本次分析范围，要同步核对后面的分母。
tnk <- subset(obj, cells = candidate_cells)
# 【子集重分析】对子集重新标准化、选择高变基因和降维，使 T/NK 内部差异更易分辨。
# Reprocess subset
message("Reclustering candidate T/NK cells: ", ncol(tnk))
# 【标准化】默认 LogNormalize 将每个细胞按总计数缩放再取 log1p，减少测序深度差异的影响。
tnk <- NormalizeData(tnk, verbose = FALSE)
# 【高变基因】选择细胞间变化较大的基因用于降维，nfeatures 控制数量；不等于删掉其他基因的原始表达。
tnk <- FindVariableFeatures(tnk, nfeatures = tnk_n_variable, verbose = FALSE)
# 【缩放】对用于分析的基因做中心化/标准化，使不同基因更便于进入 PCA；具体基因由 features 指定。
tnk <- ScaleData(tnk, verbose = FALSE)
# 【PCA】把大量基因的变化压缩成主成分；npcs 是计算数量，后续 dims 决定实际使用哪些。
tnk <- RunPCA(tnk, npcs = tnk_n_pcs, verbose = FALSE)
# 【邻居图】在选定主成分空间中寻找表达相似的细胞，为图聚类提供连接关系。
tnk <- FindNeighbors(tnk, dims = tnk_use_dims, verbose = FALSE)
# 【聚类】利用邻居图划分细胞群；resolution 越大通常分得越细，编号本身没有生物学身份。
tnk <- FindClusters(tnk, resolution = tnk_resolution, verbose = FALSE)
# 【UMAP】将表达相似性压缩到低维便于展示；二维距离不能解释为组织空间距离。
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

# 【参考阈值】用当前子聚类评分的中位数/分位数设置规则；这些是本批数据的相对阈值。
# Global thresholds from subcluster means (conservative)
med_nk <- stats::median(sub_scores$mean_nk_score, na.rm = TRUE)
med_t <- stats::median(sub_scores$mean_t_score, na.rm = TRUE)
med_ex <- stats::median(sub_scores$mean_exclude, na.rm = TRUE)
q75_ex <- as.numeric(stats::quantile(sub_scores$mean_exclude, 0.75, na.rm = TRUE))

# 【函数：assign_subcluster_label】综合 T、NK 和排除标志的群均值，给子聚类分配初步精细标签。
# 规则依次判断污染、明确 NK、明确 T 和不确定群；后面还有逐细胞修正。
# 【注释规则】先检查污染标志，再比较 T 与 NK 核心标志；不是看到 NKG7 就直接认定 NK。
assign_subcluster_label <- function(row) {
  # 【污染候选】排除标志较高且 T/NK 核心标志弱时，优先标记为其他/污染，保留后续复核空间。
  # Contaminant if exclusion markers dominate and lymphoid scores weak
  if (is.finite(row$mean_exclude) && row$mean_exclude >= q75_ex &&
      row$mean_core_nk < med_nk && row$mean_core_t < med_t) {
    return("Other_or_contaminant")
  }
  nk_high <- row$mean_nk_score >= med_nk && row$mean_core_nk >= med_nk * 0.9
  t_high <- row$mean_t_score >= med_t && row$mean_core_t >= med_t * 0.9
  nk_pref <- row$mean_nk_minus_t >= score_delta
  t_pref <- row$mean_t_minus_nk >= score_delta
  # 【高置信 NK】同时要求 NK 特征高、相对 T 特征占优，且 CD3/TCR 相关特征不过高。
  # High-confidence NK: NK markers high, TCR/CD3 axis not co-high
  if (nk_high && nk_pref && row$mean_core_t < row$mean_core_nk &&
      row$mean_core_t <= med_t) {
    return("NK_cell")
  }
  # 【明确 T】T 标志与相对偏好满足条件时归 T；不凭单个细胞毒基因判断。
  # Clear T
  if (t_high && t_pref && row$mean_core_t > row$mean_core_nk) {
    return("T_cell")
  }
  # 【保留不确定性】T 与 NK 同时较高或无法清晰区分时，优先保留 T_NK_ambiguous 标签。
  # Both high or unclear
  if (nk_high && t_high) return("T_NK_ambiguous")
  if (nk_high && !t_pref) return("T_NK_ambiguous")
  if (t_high) return("T_cell")
  "T_NK_ambiguous"
}

# 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
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
# 【逐细胞修正】从子聚类标签出发，用单细胞的 T/NK 核心分数修正不确定标签。
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
# 【最终 NK 约束】还要满足逐细胞核心 NK 高于核心 T，且 NK 分数达到设定门槛。
# Final high-confidence NK: require core NK > core T and NK score above median
final_nk_ok <- refined == "NK_cell" & cell_core_nk > cell_core_t & cell_nk >= med_cell_nk
refined[refined == "NK_cell" & !final_nk_ok] <- "T_NK_ambiguous"
tnk$refined_celltype <- refined

# analysis_celltype inside subset: only high-confidence NK as NK_cell
tnk$analysis_celltype <- tnk$refined_celltype
tnk$analysis_celltype[tnk$refined_celltype != "NK_cell" & tnk$refined_celltype != "T_cell"] <-
  ifelse(tnk$refined_celltype[tnk$refined_celltype != "NK_cell" & tnk$refined_celltype != "T_cell"] ==
           "Other_or_contaminant", "Other_or_contaminant", "T_NK_ambiguous")

# 【回填原对象】把子集的新标签按细胞名写回完整对象；未进入候选集的细胞另行保留原标签。
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
# 【旧标签处理】没有通过新 NK 规则的细胞，不能仅因为旧聚类曾叫 NK 就继续保留 NK_cell。
# Never keep legacy NK label for cells that failed high-confidence NK
legacy_nk_mask <- as.character(obj$seurat_clusters) %in% old_nk_clusters
not_new_nk <- obj$analysis_celltype != "NK_cell"
# If was legacy "NK" via cluster but not high-conf, set to refined or T_cell/ambiguous
obj$analysis_celltype[legacy_nk_mask & not_new_nk & obj$refined_celltype == "Not_in_TNK_candidate"] <- "T_NK_ambiguous"
# Epithelial and others from legacy stay
obj$analysis_celltype[obj$legacy_manual_celltype == tumor_epithelial_celltype_label &
                        obj$refined_celltype == "Not_in_TNK_candidate"] <- tumor_epithelial_celltype_label

# 【逐细胞导出】保留旧标签、新标签与条形码，方便追溯每个细胞的注释来源。
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

# 【旧群到新类型】汇总各原聚类被拆分到新类型的数量，观察是否存在混合群。
# By original cluster summary
by_cluster <- do.call(rbind, lapply(cluster_ids, function(cl) {
  idx <- as.character(obj$seurat_clusters) == cl
  # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
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

# 【检查样本来源】统计各子聚类来自哪些样本，判断是否主要由单一样本贡献。
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
  # 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
  # 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
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

# 【阶段验收】综合 T/NK 标志、数量与样本分布，记录是否满足进入下游的条件。
# QC summary / acceptance
high_nk_cells <- colnames(obj)[obj$analysis_celltype == "NK_cell"]
high_t_cells <- colnames(obj)[obj$analysis_celltype == "T_cell"]
# For T comparison use refined T in candidates + legacy T outside
t_compare_cells <- colnames(obj)[obj$refined_celltype == "T_cell" |
                                   (obj$analysis_celltype == "T_cell" & obj$refined_celltype == "Not_in_TNK_candidate")]
if (length(t_compare_cells) < 20) t_compare_cells <- high_t_cells

# 【函数：split_cluster】统计旧聚类 cl 被新注释分成了多少 NK/T/不确定/污染细胞。
# 用于检查新旧标签差异，不会在这里重新聚类。
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

# 【身份质控】要求最终 NK 的 CD3/TCR 相关表达低于 T，NK 核心标志高于对照；后面汇总是否通过。
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

# 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
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
# 【阶段 2】沿用通过质控的新身份标签，重新进行 STK31 分组、差异与候选通讯分析。
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

# 【敏感性分析】换用不同 STK31 阳性/高表达定义，查看细胞集合和后续结果是否改变。
# Sensitivity definitions
epi_idx <- obj$analysis_celltype == tumor_epithelial_celltype_label
epi_cells <- colnames(obj)[epi_idx]
# 1) counts > 0 on data layer already >0 approx; also try counts layer if available
# 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
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
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
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

# 【新旧 NK 重叠】对比旧候选与高置信 NK 的细胞集合，说明细化注释保留或排除了哪些细胞。
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

# 【样本摘要】每个样本分别计算上皮 STK31 与 NK 相关指标，避免把细胞数量误当样本数量。
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

# 【NK 功能摘要】只对高置信 NK 按样本汇总基因集平均表达；分数是表达代理指标。
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
    # 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
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
# 【阶段 3】在冻结的细胞身份上推断通讯，比较 high/low 上皮与 NK 的两个方向。
# STAGE 3: CellChat high vs low differential
# ============================================================
message("=== STAGE 3: CellChat high-vs-low ===")
# 【续跑检查】若跳过前两个阶段，需要从已有对象恢复本阶段所需的 STK31 标签。
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

  # 【通讯分组】将 STK31-high、low 上皮和 NK 等定义为 CellChat 的发送/接收群。
  # Build cellchat groups
  obj$cellchat_group <- as.character(obj$analysis_celltype)
  obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_high_tumor_epithelial"] <-
    "STK31_high_tumor_epithelial"
  obj$cellchat_group[obj$tumor_epithelial_stk31_group == "STK31_low_tumor_epithelial"] <-
    "STK31_low_tumor_epithelial"
  obj$cellchat_group[obj$analysis_celltype == "NK_cell"] <- "NK_cell"
  # 【合并稀少类别】将部分不确定类型汇总，避免大量小群干扰通讯输入；以代码中的映射为准。
  # Collapse rare ambiguous labels
  obj$cellchat_group[obj$cellchat_group %in% c("T_NK_ambiguous", "Other_or_contaminant", "Unassigned")] <-
    "Other_immune_or_unassigned"

  group_counts <- sort(table(obj$cellchat_group), decreasing = TRUE)
  keep_groups <- names(group_counts[group_counts >= min_cells_per_cellchat_group])
  # 【关键群纳入】核心比较群也必须满足最少细胞数要求，不能仅为了画图而忽略阈值。
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
    # 【抽样】从候选集合抽取元素；replace=TRUE 是有放回抽样，同一个细胞可能重复出现。
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
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_per_cellchat_group)
  # 【通路汇总】把配体受体层面的通讯聚合到通路层面，再由 aggregateNet 汇总群间网络。
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(stage3_dir, "merged_cellchat_object.rds"))

  # 【通讯表】从 CellChat 对象提取 source、target、ligand、receptor、prob 等字段，便于后续筛选与作图。
  communication <- CellChat::subsetCommunication(cellchat)
  write.csv(communication, file.path(stage3_dir, "merged_cellchat_communications.csv"), row.names = FALSE)

  # Helper to extract pair probs
  # 【函数：get_pair_prob】按 source 和 target 从通讯表取出指定方向的全部相互作用行。
  # 这里不重新计算概率，只筛选已有 prob 数据；没有匹配时返回空表。
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
  # 【函数：make_diff_table】将 high 与 low 通讯按通路/配体/受体键对齐，计算差值、比值并标记方向。
  # 返回的是两种上下文的描述性对照；不能把差值本身当作差异通讯的统计 P 值。
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
    # 【函数：add_key】统一通路字段名称，再把通路、配体、受体拼成键，用来匹配 high/low 两张表。
    # 缺少字段时使用当前代码的占位名称；空表返回空键。
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

  # 【机制轴对照】将相关配体受体归为候选轴，比较 high/low 上下文中的汇总通讯值。
  # Mechanism axis comparison
  axis_rows <- list()
  for (axis_name in names(mechanism_axes)) {
    genes <- mechanism_axes[[axis_name]]
    # 【函数：pick_axis】按机制基因及通路名称匹配相关通讯，分别加总传入的两列概率。
    # 这是基于名称规则的候选轴归类，需留意宽泛关键词可能包含的相互作用。
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

  # 【HLA 专项核查】专门检查 MHC-I/HLA 相关相互作用，避免仅依靠通路总分解释机制。
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

  # 【高低组热图】展示重点配体受体在两种上下文的分数；缺失和零值的处理见构表代码。
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

  # 【圆形网络】节点是细胞群，连线编码通讯数量或强度；先看调用的网络矩阵和图例。
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
