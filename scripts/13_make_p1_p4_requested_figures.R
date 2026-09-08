#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】针对 P1/P4 生成上皮细胞相关图
# 输入：脚本 11 的 merged_cnv_binary_annotated.rds。
# 流程：取指定患者上皮 → 按 CNV 二分类汇总目标基因的检测比例与表达 → 绘图并保存数据。
# 输出：results/merged_cnv_binary_annotation 的 figures 和 tables 子目录。
# 比例的分母取决于当前筛选的患者与细胞群；不要将局部比例当作全数据比例。
# 先看 group_stats 和 groups/group_levels 的实际定义，再理解气泡大小及颜色。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

root <- "/home/zhuweiyu/codex-r/results/merged_cnv_binary_annotation"
object_path <- file.path(root, "merged_cnv_binary_annotated.rds")
figure_dir <- file.path(root, "figures")
table_dir <- file.path(root, "tables")

# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
# 【函数：safe_fraction】计算 numerator/denominator；分母不大于零时返回 NA，避免无意义比例。
safe_fraction <- function(numerator, denominator) ifelse(denominator > 0, numerator / denominator, NA_real_)
# 【函数：fmt_pct】将 0–1 的比例乘 100 后格式化为百分数字符串；digits 控制小数位数。
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
# 【函数：save_figure】将图写入相应结果目录，通常同时生成 PNG 和 PDF。
# stem 是不含扩展名的文件名；绘图数据在主流程中另外导出。
save_figure <- function(plot, stem, width, height) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  obj <- JoinLayers(obj, assay = "RNA")
}
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "RNA counts must be sparse")
genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M")
stop_if_not(all(genes %in% rownames(counts)), "A requested gene is absent")

# 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
# := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
meta <- as.data.table(obj[[]], keep.rownames = "cell_id")
# 【患者子集】只保留 P1/P4 的上皮二分类；P4 对应旧命名 tissue5，而不是 tissue4。
meta <- meta[patient_id %in% c("P1", "P4") &
  epithelial_binary_label %in% c("Normal epithelial cells", "Malignant epithelial cells")]
meta[, display_group := paste(patient_id, sub(" epithelial cells$", "", epithelial_binary_label))]
group_levels <- c("P1 Normal", "P1 Malignant", "P4 Normal", "P4 Malignant")
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
meta[, display_group := factor(display_group, levels = group_levels)]
expected_n <- c("P1 Normal" = 509L, "P1 Malignant" = 251L,
                "P4 Normal" = 47L, "P4 Malignant" = 99L)
observed_n <- table(meta$display_group)
stop_if_not(identical(as.integer(observed_n), as.integer(expected_n)),
            "Unexpected P1/P4 epithelial group sizes")

# 【函数：group_stats】对当前设定的患者/细胞群汇总 gene 的表达值和阳性比例。
# 返回可用于绘图的分组统计表；比例的分母是各组实际细胞数。
group_stats <- function(gene) {
  # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
  rbindlist(lapply(group_levels, function(group) {
    cells <- meta[display_group == group, cell_id]
    values <- as.numeric(counts[gene, cells, drop = TRUE])
    library_umi <- sum(counts[, cells, drop = FALSE])
    data.table(
      gene = gene,
      display_group = group,
      patient_id = sub(" .*", "", group),
      epithelial_group = sub("^[^ ]+ ", "", group),
      cell_n = length(cells),
      RNA_detected_n = sum(values > 0),
      # 【检测比例】mean(values>0) 等于表达该基因的细胞数/当前组细胞总数；零计数不一定意味着真实不表达。
      RNA_detected_fraction = mean(values > 0),
      total_raw_UMI = sum(values),
      group_library_UMI = library_umi,
      # 【伪批量 CPM】该基因在组内总 UMI ÷ 该组所有基因总 UMI × 100万。
      # 这是按组聚合后的表达量；不同于先算每个细胞 CPM 再取平均，也未做 TMM 校正。
      pseudobulk_CPM = sum(values) / library_umi * 1e6,
      mean_raw_UMI_per_cell = mean(values)
    )
  }))
}

stk <- group_stats("STK31")
stk[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
stk[, display_group := factor(display_group, levels = group_levels)]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(stk, file.path(table_dir, "Figure_P1_P4_STK31_epithelial_stats.csv"))

# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p_stk_a <- ggplot(stk, aes(x = display_group, y = "STK31")) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(3, 13),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "STK31 in P1/P4 epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

stk_malignant <- stk[epithelial_group == "Malignant"]
p_stk_b <- ggplot(stk_malignant, aes(x = patient_id, y = RNA_detected_fraction)) +
  geom_col(width = 0.58, fill = "#8c2d04") +
  geom_text(
    aes(label = paste0(RNA_detected_n, " / ", cell_n, " cells\n",
                       fmt_pct(RNA_detected_fraction, 2))),
    vjust = -0.35, size = 3.4
  ) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "STK31 RNA detection fraction",
       title = "STK31 in malignant epithelial cells") +
  theme_classic(base_size = 11)

p_stk <- p_stk_a + p_stk_b + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(p_stk, "Figure_P1_P4_STK31_epithelial_localization", 11, 4.8)

hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
hla <- rbindlist(lapply(hla_genes, group_stats))
hla <- hla[epithelial_group == "Malignant"]
hla[, log1p_pseudobulk_CPM := log1p(pseudobulk_CPM)]
hla[, gene := factor(gene, levels = rev(hla_genes))]
hla[, display_group := factor(display_group, levels = c("P1 Malignant", "P4 Malignant"))]
fwrite(hla, file.path(table_dir, "Figure_P1_P4_malignant_HLAI_stats.csv"))

p_hla <- ggplot(hla, aes(x = display_group, y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL, title = "HLA-I-related genes in P1/P4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(p_hla, "Figure_P1_P4_malignant_epithelial_HLAI", 7.2, 4.8)

required <- c(
  file.path(figure_dir, paste0("Figure_P1_P4_STK31_epithelial_localization.", c("png", "pdf"))),
  file.path(figure_dir, paste0("Figure_P1_P4_malignant_epithelial_HLAI.", c("png", "pdf"))),
  file.path(table_dir, "Figure_P1_P4_STK31_epithelial_stats.csv"),
  file.path(table_dir, "Figure_P1_P4_malignant_HLAI_stats.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "A requested P1/P4 output is absent or empty")
message("P1_P4_REQUESTED_FIGURES_PASS")
print(stk)
print(hla)
