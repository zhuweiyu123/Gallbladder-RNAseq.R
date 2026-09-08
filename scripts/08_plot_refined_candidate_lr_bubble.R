# ========================================================================
# 【中文阅读指南】将人工整理的候选配体受体画成气泡图
# 输入：脚本 06 结果中的 candidate_ligand_receptor_pathways.csv，另可读取通路汇总表。
# 流程：检查列名 → 删除非有限评分 → 按 interaction_score 排序 → 选择前 top_n 个候选 → 绘图。
# 输出：results/merged_stk31_nk_refined_analysis 中的候选气泡图及相关通路图。
# 默认 top_n=30，可用 STK31_LR_BUBBLE_TOP_N 环境变量覆盖；这只改变展示数量。
# 气泡大小表示候选评分、颜色表示通路；输入来自人工候选面板，并不是 CellChatDB 的正式推断结果。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# Plot custom-panel candidate L-R bubble using refined NK results.
# Input: results/merged_stk31_nk_refined_analysis/candidate_ligand_receptor_pathways.csv
# Output: same directory + note that this is user-curated panel, not CellChatDB.

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
in_file <- file.path(
  project_dir,
  "results/merged_stk31_nk_refined_analysis/candidate_ligand_receptor_pathways.csv"
)
out_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
top_n <- as.integer(Sys.getenv("STK31_LR_BUBBLE_TOP_N", "30"))

if (!file.exists(in_file)) stop("Missing: ", in_file)

# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
lr <- read.csv(in_file, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ligand", "receptor", "pathway", "direction", "interaction_score")
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
miss <- setdiff(need, names(lr))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

# 【清理与排序】先去掉 NA/Inf，再按分数降序排列；排名仅限这份人工候选列表。
lr <- lr[is.finite(lr$interaction_score), , drop = FALSE]
lr <- lr[order(-lr$interaction_score), , drop = FALSE]
lr_top <- head(lr, top_n)
lr_top$pair <- paste(lr_top$ligand, lr_top$receptor, sep = " -> ")
# keep order by score for y axis
# 【控制显示顺序】把字符转为有固定 levels 的因子，防止纵轴默认按字母排序。
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
lr_top$pair <- factor(lr_top$pair, levels = rev(unique(lr_top$pair)))

# shorter direction labels for axis
lr_top$direction_short <- lr_top$direction
lr_top$direction_short <- gsub(
  "STK31_high_tumor_epithelial_to_NK_cell",
  "STK31-high epi -> refined NK",
  lr_top$direction_short
)
lr_top$direction_short <- gsub(
  "NK_cell_to_STK31_high_tumor_epithelial",
  "refined NK -> STK31-high epi",
  lr_top$direction_short
)

# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p <- ggplot(
  lr_top,
  aes(x = direction_short, y = pair, size = interaction_score, color = pathway)
) +
  geom_point(alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = "Top candidate ligand-receptor pairs (user-curated panel)",
    subtitle = paste0(
      "Refined NK (analysis_celltype==NK_cell, n~1551); ",
      "NOT CellChatDB. Top ", top_n, " by interaction_score."
    ),
    x = "Direction",
    y = "Ligand -> Receptor",
    size = "Score",
    color = "Pathway"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

out_pdf <- file.path(out_dir, "refined_candidate_lr_bubble.pdf")
# 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
ggsave(out_pdf, p, width = 11, height = 8)
message("Wrote: ", out_pdf)

# also pathway-level bubble/bar helper from summary if present
sum_file <- file.path(out_dir, "candidate_pathway_summary.csv")
if (file.exists(sum_file)) {
  ps <- read.csv(sum_file, stringsAsFactors = FALSE)
  ps <- ps[order(-ps$interaction_score), , drop = FALSE]
  ps_top <- head(ps, 20)
  ps_top$pathway <- factor(ps_top$pathway, levels = rev(unique(ps_top$pathway)))
  ps_top$direction_short <- gsub(
    "STK31_high_tumor_epithelial_to_NK_cell",
    "STK31-high epi -> refined NK",
    ps_top$direction
  )
  ps_top$direction_short <- gsub(
    "NK_cell_to_STK31_high_tumor_epithelial",
    "refined NK -> STK31-high epi",
    ps_top$direction_short
  )
  p2 <- ggplot(
    ps_top,
    aes(x = direction_short, y = pathway, size = interaction_score, color = pathway)
  ) +
    geom_point(alpha = 0.9) +
    theme_bw(base_size = 12) +
    labs(
      title = "Candidate pathway summary (user-curated panel, refined NK)",
      subtitle = "Max interaction_score per pathway x direction",
      x = "Direction", y = "Pathway", size = "Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )
  out2 <- file.path(out_dir, "refined_candidate_pathway_bubble.pdf")
  ggsave(out2, p2, width = 10, height = 7)
  message("Wrote: ", out2)
}

writeLines(
  c(
    "# refined_candidate_lr_bubble",
    "",
    paste0("Date: ", Sys.time()),
    paste0("Input: ", in_file),
    paste0("Top N pairs: ", top_n),
    "",
    "This figure uses the **user-curated** ligand-receptor panel",
    "(scripts/06 lr_reference), scored on STK31-high epithelial vs **refined NK**.",
    "It is NOT generated from CellChatDB.",
    "For CellChat high-vs-low bubble see:",
    "results/merged_cellchat_refined_high_low_nk/cellchat_high_low_nk_bubble.pdf"
  ),
  file.path(out_dir, "refined_candidate_lr_bubble_README.md")
)

message("Done.")
