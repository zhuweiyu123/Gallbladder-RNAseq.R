#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】补齐 CNV 分组下的表达、GO/GSEA 和 CellChat 图
# 输入：脚本 11 的整合对象及已有差异表；CODEX_R_PROJECT_DIR 可指定项目目录。
# 流程：回读并核对分组 → 免疫靶点图/热图/模块分数 → 已有差异表做 GO/GSEA → 重算 CellChat。
# 输出：CNV 整合结果下 figures_stk31_de_go_cellchat 和 tables_stk31_de_go_cellchat。
# STK31-high/low 仍按恶性上皮中的原始计数 >0/==0；不要与旧分位数定义混用。
# ORA 用筛选后的差异基因；GSEA 用读入差异表中可排序的基因，不代表自动恢复了全基因组列表。
# 所需 GO 包缺失时会写跳过说明；占位说明不是成功完成富集的证据。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================

# Complete the missing CNV-binary STK31 figures that require the full Seurat
# object and R/CellChat/clusterProfiler environment.

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260902)

# 【函数：pick_project_dir】依次尝试环境变量、当前目录和服务器默认目录，选择含目标 RDS 的项目目录。
# 若都找不到就停止，避免在错误目录继续运行。
pick_project_dir <- function() {
  candidates <- unique(c(
    # 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
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
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 【参数 padj_cut】多重检验调整后 P 值的筛选上限；不能用原始 P 值直接替代。
padj_cut <- 0.05
# 【参数 logfc_cut】筛选差异基因的绝对 log2 倍数变化下限；0.25 约对应 1.19 倍变化。
logfc_cut <- 0.25
# 【参数 min_cells_cellchat】CNV 通讯分析允许进入网络的最少群细胞数。
min_cells_cellchat <- as.integer(Sys.getenv("CELLCHAT_MIN_CELLS", "10"))
# 【参数 max_cells_cellchat】CNV 通讯分析每群抽样上限；与最终导出的实际用量一起核对。
max_cells_cellchat <- as.integer(Sys.getenv("CELLCHAT_MAX_CELLS_PER_GROUP", "500"))

# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

# 【函数：save_both】把同一张 ggplot 图导出为 PNG 与 PDF；stem 决定基础文件名。
# 尺寸参数影响字体和图形布局，同名运行会更新对应输出。
save_both <- function(plot, stem, width, height) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# 【函数：fetch_layer】按指定层、基因和细胞取矩阵；未指定基因或细胞时保留全部。
# 先核对行列名称，counts 与 data 的数值尺度不同，不能混用。
fetch_layer <- function(object, layer_name = "data", genes = NULL, cells = NULL) {
  mat <- if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = "RNA", layer = layer_name)
  } else {
    GetAssayData(object, assay = "RNA", slot = layer_name)
  }
  # 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
  if (!is.null(genes)) genes <- intersect(genes, rownames(mat))
  if (is.null(genes)) genes <- rownames(mat)
  if (is.null(cells)) cells <- colnames(mat)
  mat[genes, cells, drop = FALSE]
}

# 【函数：zscore_rows】先转置、按列 scale、再转回，等价于对原矩阵每一行做 z-score。
# 颜色表示某基因相对自身均值的高低；常数行标准差为零时可能产生非有限值。
zscore_rows <- function(mat) {
  t(scale(t(as.matrix(mat))))
}

message("Reading CNV-binary annotated Seurat object: ", object_path)
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  obj <- JoinLayers(obj, assay = "RNA")
}
meta <- obj[[]]
need_cols <- c("patient_id", "epithelial_binary_label", "analysis_celltype_cnv_binary")
stop_if_not(all(need_cols %in% colnames(meta)),
            # 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
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
# 【分组】仅在 malignant 细胞内按 counts>0 取 high；未检出 counts==0 的恶性细胞为 low。
stk31_high <- malignant & stk31_counts > 0
stk31_low <- malignant & stk31_counts == 0
nk_cells <- if ("analysis_celltype_previous_or_new" %in% colnames(meta)) {
  !is.na(meta$analysis_celltype_previous_or_new) &
    meta$analysis_celltype_previous_or_new == "NK_cell"
} else {
  !is.na(meta$analysis_celltype_cnv_binary) &
    meta$analysis_celltype_cnv_binary == "NK_cell"
}

# 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
# := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
validation <- data.table(
  check = c("malignant_epithelial", "normal_epithelial", "stk31_high_malignant",
            "stk31_low_malignant", "nk_cells", "p2_likely_in_malignant_rule"),
  observed = c(sum(malignant, na.rm = TRUE), sum(normal, na.rm = TRUE),
               sum(stk31_high, na.rm = TRUE), sum(stk31_low, na.rm = TRUE),
               sum(nk_cells, na.rm = TRUE), "malignant_likely encoded as malignant"),
  expected = c(1047, 3716, 117, 930, 1479, "malignant_likely encoded as malignant")
)
validation[, pass := observed == expected]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(validation, file.path(table_dir, "R_object_missing_figure_validation.tsv"), sep = "\t")
stop_if_not(all(validation$pass), "Object validation failed")

group_label <- rep(NA_character_, ncol(obj))
group_label[stk31_high] <- "STK31_high_malignant_epithelial"
group_label[stk31_low] <- "STK31_low_malignant_epithelial"
group_label[normal] <- "Normal_epithelial"
group_label[nk_cells] <- "NK_cell"
obj$stk31_cnv_binary_plot_group <- group_label

# ------------------------------------------------------------
# 【表达概览】比较各分组中免疫靶基因的平均表达与检出比例，绘制点图。
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
# 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
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
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
dot_rows[, group := factor(group, levels = plot_groups)]
dot_rows[, gene := factor(gene, levels = rev(immune_genes))]
# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
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
# 【表达热图】按患者及 high/low 汇总恶性上皮表达，再选择基因绘图；读颜色前先看是否做了逐行标准化。
# Top expression heatmap for STK31-high vs low malignant cells
# ------------------------------------------------------------
# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
de_hl <- fread(file.path(de_dir, "markers_stk31_high_vs_low_malignant_epithelial.csv"))
top_up <- head(de_hl[avg_log2FC > 0][order(p_val_adj, -abs(avg_log2FC)), gene], 15)
top_down <- head(de_hl[avg_log2FC < 0][order(p_val_adj, -abs(avg_log2FC)), gene], 15)
heat_genes <- intersect(unique(c(top_up, top_down)), rownames(data_mat))
avg_rows <- rbindlist(lapply(unique(meta$patient_id[malignant]), function(patient) {
  rbindlist(lapply(c("STK31_high", "STK31_low"), function(status) {
    status_mask <- if (status == "STK31_high") stk31_high else stk31_low
    idx <- which(meta$patient_id == patient & status_mask)
    if (length(idx) == 0) return(NULL)
    # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
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
# 【模块评分】把人工定义的免疫基因集合概括为每个细胞的表达分数，再按分组比较。
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
  # 【按细胞求均值】在基因×细胞矩阵中，colMeans 对每列跨基因求平均，常用来构建逐细胞模块分数。
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
# 【富集分析】使用已经计算好的差异表；并不是在这个区段重新运行差异检验。
# GO / GSEA
# ------------------------------------------------------------
# 【函数：write_go_skip】把 GO/GSEA 未执行的原因保存到结果目录，使缺少结果可追溯。
write_go_skip <- function(prefix, note) {
  fwrite(data.table(note = note), file.path(table_dir, paste0(prefix, "_GO_GSEA_skipped.csv")))
}

# 【函数：run_go_gsea】从已有差异表读入 fold change/P 值/调整后 P 值，通过参数适配 edgeR 或 Seurat 列名。
# 排序分数为 sign(logFC)×[-log10(P)]；GSEA 用排序列表，ORA 用通过筛选的基因。
# 二者回答的问题不同；这里只使用差异表提供的基因，并输出 ID 转换与跳过记录。
run_go_gsea <- function(de_path, prefix, fc_col, p_col, padj_col) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    write_go_skip(prefix, "clusterProfiler/org.Hs.eg.db not installed")
    return(invisible(NULL))
  }
  de <- fread(de_path)
  de <- de[is.finite(get(fc_col)) & is.finite(get(p_col))]
  # 【GSEA 排序】正负由 fold change 决定，绝对值由 -log10(P) 决定；pmax 防止 P=0 导致无穷大。
  de[, rank_score := sign(get(fc_col)) * -log10(pmax(get(p_col), 1e-300))]
  de <- de[order(-abs(rank_score))]
  de <- de[!duplicated(gene)]
  # 【基因 ID】将常见基因符号转换为数据库识别的 ENTREZID；有些符号无法匹配或一对多，需看转换后数量。
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

  # 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
  gsea <- tryCatch(
    # 【GSEA】使用有方向的基因排序列表；NES 符号表示该基因集偏向排序顶部还是底部，显著性另看调整后 P 值。
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
  # 【GO 过度代表分析】检查候选基因是否集中于某些 GO 条目；BP 是生物过程，BH 用于多重检验调整。
  # 如果没有显式传入 universe，就使用数据库可用背景；富集不直接说明通路被激活。
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

# 【函数：summarize_cellchat_axes】按预定义目标细胞对和机制基因集合汇总 CellChat 通讯表。
# 轴级总分是已推断相互作用的摘要，不能独立证明因果调控。
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
  # 【按名称对齐】match(x,y) 返回 x 各元素在 y 中的位置，没找到返回 NA；用来保证标签和表达对应同一细胞。
  split_cells <- split(keep_cells, obj$cellchat_cnv_binary_group[match(keep_cells, cells)])
  keep_cells <- unlist(lapply(split_cells, function(x) {
    # 【抽样】从候选集合抽取元素；replace=TRUE 是有放回抽样，同一个细胞可能重复出现。
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
  # 【通讯输入】表达矩阵列名与 meta 行名必须对应；group.by 指定用哪个标签定义发送和接收细胞群。
  cellchat <- CellChat::createCellChat(object = cc_data, meta = cc_meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  # 【通讯候选筛选】先找群中高表达的基因，再通过配体受体数据库筛选候选相互作用。
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  # 【通讯推断】根据群表达与配体受体库估计通讯分值；raw.use=TRUE 指未投影的表达，不等于原始 counts。
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  # 【通讯过滤】去掉不满足最少细胞数要求的群的通讯，min.cells 是本次阈值。
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells_cellchat)
  # 【通路汇总】把配体受体层面的通讯聚合到通路层面，再由 aggregateNet 汇总群间网络。
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  # 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
  saveRDS(cellchat, file.path(table_dir, "CellChat_cnv_binary_object.rds"))
  # 【通讯表】从 CellChat 对象提取 source、target、ligand、receptor、prob 等字段，便于后续筛选与作图。
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
