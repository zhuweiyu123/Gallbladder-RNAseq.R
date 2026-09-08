# ========================================================================
# 【中文阅读指南】补充 STK31-high 上皮差异分析与 GO
# 输入：脚本 06 的 refined_merged_object.rds，沿用该对象中的高/低组和 NK 注释。
# 流程：确定比较细胞集合 → Wilcoxon 差异分析 → 导出高低方向基因 → GO BP → 绘制摘要。
# 输出：results/merged_stk31_nk_refined_analysis 下的差异表、GO 表和相关图。
# STK31_MAX_CELLS_DE 默认每组最多抽取 500 个细胞做差异检验，以控制计算时间。
# high 与 low 的比较较聚焦于上皮状态；high 与所有其他细胞的比较还包含细胞类型组成的差异。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# ============================================================
# 10_refined_stk31_high_vs_low_and_all_go.R
# On frozen refined object:
#   1) STK31-high epithelial vs STK31-low epithelial DE + GO BP
#   2) STK31-high epithelial vs all other cells DE + GO BP
# Epithelial from analysis_celltype == "Epithelial"
# STK31-high: epithelial counts > 0 (matches freeze when q-cutoff==0)
# ============================================================

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
input_rds <- file.path(
  project_dir, "results/merged_tnk_refined_annotation/refined_merged_object.rds"
)
out_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
epi_label <- "Epithelial"
# 【参数 padj_cut】多重检验调整后 P 值的筛选上限；不能用原始 P 值直接替代。
padj_cut <- 0.05
# 【参数 logfc_cut】筛选差异基因的绝对 log2 倍数变化下限；0.25 约对应 1.19 倍变化。
logfc_cut <- 0.25
# 【参数 min_genes_go】进入 GO 的最少候选基因数，过少时按代码跳过并写说明。
min_genes_go <- 10
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
max_cells_per_ident <- as.integer(Sys.getenv("STK31_MAX_CELLS_DE", "500"))

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing: ", path)
}

# 【函数：available_genes】把想看的基因与对象实际行名取交集，避免访问不存在的基因。
# 返回基因名向量；返回长度为零表示这些基因不在当前矩阵内。
# 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
available_genes <- function(object, genes) intersect(genes, rownames(object))

# 【函数：fetch_mat】提取当前 assay 中指定基因的矩阵，并兼容 Seurat 的 layer/slot 差异。
# drop=FALSE 保留二维结构，防止只取一个基因时自动变成向量。
fetch_mat <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NULL)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

# 【函数：run_de】为两组细胞建立临时标签，运行 FindMarkers，排序并导出差异基因表。
# 第一组/第二组的传入顺序决定 logFC 正负；每组抽样上限见 max.cells.per.ident。
run_de <- function(object, cells.1, cells.2, out_file) {
  if (length(cells.1) < 10 || length(cells.2) < 10) {
    stop("Too few cells: n1=", length(cells.1), " n2=", length(cells.2))
  }
  object$.cmp <- "unused"
  object$.cmp[colnames(object) %in% cells.1] <- "g1"
  object$.cmp[colnames(object) %in% cells.2] <- "g2"
  message("DE: n1=", length(cells.1), " n2=", length(cells.2),
          " max.cells.per.ident=", max_cells_per_ident)
  # 【差异表达】比较 ident.1 与 ident.2；avg_log2FC>0 表示第一组更高，p_val_adj 是多重检验调整后 P 值。
  # min.pct/logfc.threshold 是进入检验的筛选条件；细胞级检验不自动控制患者内相关性。
  markers <- FindMarkers(
    object,
    ident.1 = "g1",
    ident.2 = "g2",
    group.by = ".cmp",
    logfc.threshold = 0.1,
    min.pct = 0.1,
    test.use = "wilcox",
    max.cells.per.ident = max_cells_per_ident
  )
  markers$gene <- rownames(markers)
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
  # 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
  write.csv(markers, out_file, row.names = FALSE)
  markers
}

# 【函数：run_go_bp】筛选差异基因、转换 SYMBOL 到 ENTREZID，再进行 GO 生物过程（BP）富集。
# 主列表包含满足阈值的上调和下调基因；上/下调数量另外记录。
# 富集解释依赖候选列表与背景，具体背景是否显式设置以 enrichGO 调用为准。
run_go_bp <- function(markers, output_prefix, comparison_name) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("clusterProfiler/org.Hs.eg.db missing")
  }
  sig <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      !is.na(markers$avg_log2FC) & abs(markers$avg_log2FC) > logfc_cut
  ]
  sig_up <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      markers$avg_log2FC > logfc_cut
  ]
  sig_down <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      markers$avg_log2FC < -logfc_cut
  ]
  sig <- unique(sig[!is.na(sig)])
  write.csv(
    data.frame(
      comparison = comparison_name,
      n_sig = length(sig),
      n_sig_up = length(unique(sig_up[!is.na(sig_up)])),
      n_sig_down = length(unique(sig_down[!is.na(sig_down)])),
      padj_cut = padj_cut,
      logfc_cut = logfc_cut,
      stringsAsFactors = FALSE
    ),
    paste0(output_prefix, "_go_input_summary.csv"),
    row.names = FALSE
  )
  message(comparison_name, ": n_sig=", length(sig),
          " up=", length(unique(sig_up[!is.na(sig_up)])),
          " down=", length(unique(sig_down[!is.na(sig_down)])))

  if (length(sig) < min_genes_go) {
    write.csv(data.frame(note = "skipped_few_genes"), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  suppressMessages({
    # 【基因 ID】将常见基因符号转换为数据库识别的 ENTREZID；有些符号无法匹配或一对多，需看转换后数量。
    converted <- clusterProfiler::bitr(
      sig, fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  })
  if (nrow(converted) < min_genes_go) {
    write.csv(data.frame(note = "skipped_few_entrez"), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  # 【GO 过度代表分析】检查候选基因是否集中于某些 GO 条目；BP 是生物过程，BH 用于多重检验调整。
  # 如果没有显式传入 universe，就使用数据库可用背景；富集不直接说明通路被激活。
  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
  message(comparison_name, ": n_GO=", nrow(ego_df))

  if (nrow(ego_df) > 0) {
    top <- head(ego_df[order(ego_df$p.adjust), ], 20)
    # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
    top$Description <- factor(top$Description, levels = rev(top$Description))
    if (!"Count" %in% colnames(top)) {
      # 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
      top$Count <- vapply(strsplit(as.character(top$GeneRatio), "/"), function(z) as.numeric(z[1]), numeric(1))
    }
    top$score <- -log10(pmax(top$p.adjust, 1e-300))
    # 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
    p <- ggplot(top, aes(x = score, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#d73027", high = "#4575b4", trans = "log10") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0("GO BP: ", comparison_name),
        subtitle = "Refined object; epithelial=analysis_celltype; STK31-high=counts>0",
        x = "-log10(FDR)", y = NULL
      )
    # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
    ggsave(paste0(output_prefix, "_go_bp_dotplot.pdf"), p,
           width = 10, height = max(5, 0.28 * nrow(top) + 2))
  }
  invisible(ego_df)
}

# 【函数：top_genes_table】从差异表整理前 n 个展示基因和相关统计值。
# 前几名的排序不等于生物学重要性的最终排序。
top_genes_table <- function(markers, n = 15) {
  up <- head(markers[markers$avg_log2FC > 0, ], n)
  down <- head(markers[markers$avg_log2FC < 0, ], n)
  list(up = up, down = down)
}

message("=== 10 STK31 high vs low / vs all GO ===")
message("Start: ", Sys.time())
stop_if_missing(input_rds)

# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  try(obj <- JoinLayers(obj, assay = "RNA"), silent = TRUE)
}
if (!"analysis_celltype" %in% colnames(obj@meta.data)) stop("analysis_celltype missing")

# 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
stk31_counts <- tryCatch(
  as.numeric(fetch_mat(obj, target_gene, "counts")[target_gene, ]),
  error = function(e) as.numeric(fetch_mat(obj, target_gene, "data")[target_gene, ])
)
is_epi <- obj$analysis_celltype == epi_label
if (sum(is_epi) == 0 && "legacy_manual_celltype" %in% colnames(obj@meta.data)) {
  is_epi <- obj$legacy_manual_celltype == epi_label
}
high_cells <- colnames(obj)[is_epi & stk31_counts > 0]
low_cells <- colnames(obj)[is_epi & !(stk31_counts > 0)]
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
other_cells <- setdiff(colnames(obj), high_cells)

summary_df <- data.frame(
  item = c(
    "epithelial_cells", "stk31_high_epithelial", "stk31_low_epithelial",
    "all_other_vs_high", "stk31_definition", "object"
  ),
  value = c(
    sum(is_epi), length(high_cells), length(low_cells),
    length(other_cells), "epithelial & counts>0",
    basename(input_rds)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "stk31_high_low_all_de_summary.csv"), row.names = FALSE)
print(summary_df)

message("1) DE high vs low epithelial")
m_hl <- run_de(
  obj, high_cells, low_cells,
  file.path(out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv")
)
message("2) GO high vs low epithelial")
run_go_bp(
  m_hl,
  file.path(out_dir, "stk31_high_vs_low_tumor_epithelial"),
  "STK31_high_vs_low_epithelial"
)

message("3) DE high epithelial vs all other")
m_ha <- run_de(
  obj, high_cells, other_cells,
  file.path(out_dir, "markers_stk31_high_tumor_epithelial_vs_all_other.csv")
)
message("4) GO high epithelial vs all other")
run_go_bp(
  m_ha,
  file.path(out_dir, "stk31_high_tumor_epithelial_vs_all_other"),
  "STK31_high_epithelial_vs_all_other"
)

# top gene exports for quick reading
write.csv(head(m_hl[m_hl$avg_log2FC > 0, ], 50),
          file.path(out_dir, "top50_up_stk31_high_vs_low_epithelial.csv"), row.names = FALSE)
write.csv(head(m_hl[m_hl$avg_log2FC < 0, ], 50),
          file.path(out_dir, "top50_down_stk31_high_vs_low_epithelial.csv"), row.names = FALSE)
write.csv(head(m_ha[m_ha$avg_log2FC > 0, ], 50),
          file.path(out_dir, "top50_up_stk31_high_vs_all_other.csv"), row.names = FALSE)
write.csv(head(m_ha[m_ha$avg_log2FC < 0, ], 50),
          file.path(out_dir, "top50_down_stk31_high_vs_all_other.csv"), row.names = FALSE)

# brief text summary of top GO
# 【函数：summarize_go_top】读取已生成的 GO 表并整理排名靠前的条目，便于总览多个比较。
summarize_go_top <- function(path, n = 10) {
  if (!file.exists(path)) return(character(0))
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  g <- read.csv(path, stringsAsFactors = FALSE)
  if (!"Description" %in% names(g) || nrow(g) == 0) return(character(0))
  head(g$Description[order(g$p.adjust)], n)
}

hl_go <- summarize_go_top(file.path(out_dir, "stk31_high_vs_low_tumor_epithelial_go_bp.csv"))
ha_go <- summarize_go_top(file.path(out_dir, "stk31_high_tumor_epithelial_vs_all_other_go_bp.csv"))

writeLines(c(
  "# STK31-high epithelial DE/GO (refined object)",
  paste0("Date: ", Sys.time()),
  paste0("High cells: ", length(high_cells), "; Low epithelial: ", length(low_cells),
         "; All other: ", length(other_cells)),
  "",
  "## Meaning",
  "- high vs low: within-epithelial STK31-associated transcriptome (primary for STK31 biology).",
  "- high vs all other: high-epi vs rest of atlas (mixes cell-type identity + STK31); interpret carefully.",
  "",
  "## Top GO high vs low",
  paste0("- ", hl_go),
  "",
  "## Top GO high vs all other",
  paste0("- ", ha_go)
), file.path(out_dir, "stk31_high_low_all_de_go_README.md"))

message("Done: ", out_dir)
message("End: ", Sys.time())
