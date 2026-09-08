# ========================================================================
# 【中文阅读指南】对细化注释后的 NK 差异基因做 GO 富集
# 输入：results/merged_stk31_nk_refined_analysis 中已有的 NK 比较差异基因 CSV。
# 流程：筛选调整后 P 值和倍数变化 → 统计上/下调数量 → 转换基因 ID → GO BP 富集 → 导出表图。
# 输出：与差异基因表同目录的 GO 输入摘要、GO BP 结果和图。
# padj_cut=0.05、logfc_cut=0.25；min_genes=10 是此脚本设定的最低基因数，不是普遍生物学定律。
# 主富集列表合并了满足条件的上调与下调基因；富集出现某过程，不等于证明该过程被激活。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# ============================================================
# 09_refined_nk_go_enrichment.R
# GO BP for refined high-confidence NK DE comparisons.
# Inputs (already from script 06 Stage2):
#   markers_nk_vs_all_other.csv
#   markers_stk31_high_tumor_epithelial_vs_nk.csv
# NK definition: analysis_celltype == "NK_cell" (n=1551)
# ============================================================

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
in_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
out_dir <- in_dir
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 【参数 padj_cut】多重检验调整后 P 值的筛选上限；不能用原始 P 值直接替代。
padj_cut <- 0.05
# 【参数 logfc_cut】筛选差异基因的绝对 log2 倍数变化下限；0.25 约对应 1.19 倍变化。
logfc_cut <- 0.25
min_genes <- 10

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

# 【函数：run_go_bp】筛选差异基因、转换 SYMBOL 到 ENTREZID，再进行 GO 生物过程（BP）富集。
# 主列表包含满足阈值的上调和下调基因；上/下调数量另外记录。
# 富集解释依赖候选列表与背景，具体背景是否显式设置以 enrichGO 调用为准。
run_go_bp <- function(markers, output_prefix, comparison_name) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("clusterProfiler/org.Hs.eg.db not installed")
  }

  if (!"gene" %in% colnames(markers)) {
    if (!is.null(rownames(markers)) && !all(rownames(markers) %in% c(as.character(seq_len(nrow(markers)))))) {
      markers$gene <- rownames(markers)
    } else {
      stop("markers table has no gene column: ", comparison_name)
    }
  }

  if ("p_val_adj" %in% colnames(markers)) {
    sig <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & abs(markers$avg_log2FC) > logfc_cut
    ]
    # also report up/down separately for optional files
    sig_up <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & markers$avg_log2FC > logfc_cut
    ]
    sig_down <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & markers$avg_log2FC < -logfc_cut
    ]
  } else {
    stop("Expected p_val_adj column in markers for ", comparison_name)
  }

  sig <- unique(sig[!is.na(sig) & nzchar(sig)])
  message(comparison_name, ": n_sig_genes=", length(sig),
          " (up=", length(unique(sig_up)), ", down=", length(unique(sig_down)), ")")

  # 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
  write.csv(
    data.frame(
      comparison = comparison_name,
      n_sig = length(sig),
      n_sig_up = length(unique(sig_up[!is.na(sig_up)])),
      n_sig_down = length(unique(sig_down[!is.na(sig_down)])),
      padj_cut = padj_cut,
      logfc_cut = logfc_cut,
      nk_definition = "analysis_celltype==NK_cell refined high-confidence",
      stringsAsFactors = FALSE
    ),
    paste0(output_prefix, "_go_input_summary.csv"),
    row.names = FALSE
  )

  if (length(sig) < min_genes) {
    warning("Fewer than ", min_genes, " significant genes; skip GO for ", comparison_name)
    write.csv(
      data.frame(note = paste0("skipped: n_sig=", length(sig))),
      paste0(output_prefix, "_go_bp.csv"),
      row.names = FALSE
    )
    return(invisible(NULL))
  }

  suppressMessages({
    # 【基因 ID】将常见基因符号转换为数据库识别的 ENTREZID；有些符号无法匹配或一对多，需看转换后数量。
    converted <- clusterProfiler::bitr(
      sig,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  })
  message(comparison_name, ": n_entrez=", length(unique(converted$ENTREZID)))
  if (nrow(converted) < min_genes) {
    warning("Fewer than ", min_genes, " converted genes; skip GO for ", comparison_name)
    write.csv(
      data.frame(note = paste0("skipped: n_converted=", nrow(converted))),
      paste0(output_prefix, "_go_bp.csv"),
      row.names = FALSE
    )
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
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
  message(comparison_name, ": n_GO_terms=", nrow(ego_df))

  # simple dotplot of top terms
  if (nrow(ego_df) > 0) {
    top <- head(ego_df[order(ego_df$p.adjust), ], 20)
    # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
    top$Description <- factor(top$Description, levels = rev(top$Description))
    if (!"Count" %in% colnames(top)) {
      # 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
      top$Count <- vapply(strsplit(top$GeneRatio, "/"), function(z) as.numeric(z[1]), numeric(1))
    }
    top$score <- -log10(pmax(top$p.adjust, 1e-300))
    # 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
    p <- ggplot(top, aes(x = score, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#d73027", high = "#4575b4", trans = "log10") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0("GO BP: ", comparison_name),
        subtitle = "Refined high-confidence NK (analysis_celltype==NK_cell)",
        x = "-log10(adjusted P)",
        y = NULL,
        size = "Genes",
        color = "FDR"
      )
    # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
    ggsave(paste0(output_prefix, "_go_bp_dotplot.pdf"), p, width = 10, height = max(5, 0.28 * nrow(top) + 2))
  }

  invisible(ego_df)
}

message("=== 09 refined NK GO enrichment ===")
message("Start: ", Sys.time())

nk_file <- file.path(in_dir, "markers_nk_vs_all_other.csv")
stk31_nk_file <- file.path(in_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv")
stop_if_missing(nk_file)
stop_if_missing(stk31_nk_file)

# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
nk_markers <- read.csv(nk_file, stringsAsFactors = FALSE, check.names = FALSE)
stk31_nk_markers <- read.csv(stk31_nk_file, stringsAsFactors = FALSE, check.names = FALSE)

message("Running GO: refined NK vs all other")
run_go_bp(
  nk_markers,
  file.path(out_dir, "nk_vs_all_other"),
  "refined_NK_vs_all_other"
)

message("Running GO: STK31-high epithelial vs refined NK")
run_go_bp(
  stk31_nk_markers,
  file.path(out_dir, "stk31_high_tumor_epithelial_vs_nk"),
  "STK31_high_epithelial_vs_refined_NK"
)

# inventory note
note <- c(
  "# Refined NK GO outputs",
  paste0("Date: ", Sys.time()),
  "NK definition: analysis_celltype == NK_cell (refined high-confidence, n=1551)",
  "DE inputs:",
  paste0("  - ", nk_file),
  paste0("  - ", stk31_nk_file),
  paste0("Thresholds: p_val_adj < ", padj_cut, ", |avg_log2FC| > ", logfc_cut),
  "Outputs:",
  "  - nk_vs_all_other_go_bp.csv / _go_bp_dotplot.pdf / _go_input_summary.csv",
  "  - stk31_high_tumor_epithelial_vs_nk_go_bp.csv / _go_bp_dotplot.pdf / _go_input_summary.csv",
  "Legacy GO under merged_stk31_nk_analysis/ used old NK=7419 and must not be used for primary inference."
)
writeLines(note, file.path(out_dir, "refined_nk_go_README.md"))

message("Done. Outputs in: ", out_dir)
message("End: ", Sys.time())
