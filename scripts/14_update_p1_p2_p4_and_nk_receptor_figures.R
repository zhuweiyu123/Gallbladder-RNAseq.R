#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】更新 P1/P2/P4 上皮与本地 NK 受体图
# 输入：CNV 整合对象；上皮用最终二分类，NK 使用保留的旧高置信注释。
# 流程：筛选 P1/P2/P4 → 汇总 STK31/HLA 指标 → 绘制上皮图 → 汇总 NK 抑制性受体表达。
# 输出：CNV 整合结果目录中的图和绘图统计表，具体名称见 save_figure 与 fwrite 调用。
# P2 中 malignant_likely 按已有二分类进入恶性组；本脚本不重新估计 CNV。
# 检测比例和平均表达回答不同问题：前者是有多少细胞表达，后者是该群体平均表达水平。
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

project <- "/home/zhuweiyu/codex-r"
local_root <- file.path(project, "results", "merged_cnv_binary_annotation")
object_path <- file.path(local_root, "merged_cnv_binary_annotated.rds")
local_figure_dir <- file.path(local_root, "figures")
local_table_dir <- file.path(local_root, "tables")

# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
# 【函数：fmt_pct】将 0–1 的比例乘 100 后格式化为百分数字符串；digits 控制小数位数。
fmt_pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
# 【函数：save_figure】将图写入相应结果目录，通常同时生成 PNG 和 PDF。
# stem 是不含扩展名的文件名；绘图数据在主流程中另外导出。
save_figure <- function(plot, directory, stem, width, height) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(file.path(directory, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# 【上皮作图】限定三位患者，用最终 CNV 工作标签汇总 STK31/HLA，既有分类规则保持不变。
# P1/P2/P4 epithelial figures from the final integrated object.
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  obj <- JoinLayers(obj, assay = "RNA")
}
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "RNA counts must be sparse")
target_genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
stop_if_not(all(target_genes %in% rownames(counts)), "A requested epithelial target gene is absent")

patients <- c("P1", "P2", "P4")
# 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
# := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
meta_all <- as.data.table(obj[[]], keep.rownames = "cell_id")
meta <- meta_all[
  patient_id %in% patients &
    epithelial_binary_label %in% c("Normal epithelial cells", "Malignant epithelial cells")
]
meta[, display_group := paste(patient_id, sub(" epithelial cells$", "", epithelial_binary_label))]
group_levels <- as.vector(rbind(
  paste(patients, "Normal"), paste(patients, "Malignant")
))
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
meta[, display_group := factor(display_group, levels = group_levels)]
expected_n <- c(
  "P1 Normal" = 509L, "P1 Malignant" = 251L,
  "P2 Normal" = 1485L, "P2 Malignant" = 697L,
  "P4 Normal" = 47L, "P4 Malignant" = 99L
)
stop_if_not(
  identical(as.integer(table(meta$display_group)), as.integer(expected_n[group_levels])),
  "Unexpected P1/P2/P4 epithelial group sizes"
)

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
      RNA_detected_fraction = mean(values > 0),
      total_raw_UMI = sum(values),
      group_library_UMI = library_umi,
      # 【表达尺度】组内该基因总 UMI/组内全部 UMI×100万是伪批量 CPM；后续 log1p 用于压缩绘图动态范围。
      pseudobulk_CPM = sum(values) / library_umi * 1e6,
      mean_raw_UMI_per_cell = mean(values)
    )
  }))
}

stk <- group_stats("STK31")
stk[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  display_group = factor(display_group, levels = group_levels)
)]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(stk, file.path(local_table_dir, "Figure_P1_P2_P4_STK31_epithelial_stats.csv"))

# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p_stk_a <- ggplot(stk, aes(x = display_group, y = "STK31")) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(3, 13),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#dceefb", high = "#084594") +
  labs(x = NULL, y = NULL, title = "STK31 in P1/P2/P4 epithelial cells") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

stk_malignant <- stk[epithelial_group == "Malignant"]
stk_malignant[, patient_id := factor(patient_id, levels = patients)]
p_stk_b <- ggplot(stk_malignant, aes(x = patient_id, y = RNA_detected_fraction)) +
  geom_col(width = 0.58, fill = "#8c2d04") +
  geom_text(
    aes(label = paste0(RNA_detected_n, " / ", cell_n, " cells\n",
                       fmt_pct(RNA_detected_fraction, 2))),
    vjust = -0.35, size = 3.2
  ) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "STK31 RNA detection fraction",
       title = "STK31 in malignant epithelial cells") +
  theme_classic(base_size = 11)

p_stk <- p_stk_a + p_stk_b + plot_layout(widths = c(1.75, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(
  p_stk, local_figure_dir, "Figure_P1_P2_P4_STK31_epithelial_localization",
  12.5, 5.0
)

hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
hla <- rbindlist(lapply(hla_genes, group_stats))[epithelial_group == "Malignant"]
hla[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  gene = factor(gene, levels = rev(hla_genes)),
  display_group = factor(display_group, levels = paste(patients, "Malignant"))
)]
fwrite(hla, file.path(local_table_dir, "Figure_P1_P2_P4_malignant_HLAI_plus_HLAE_stats.csv"))

p_hla <- ggplot(hla, aes(x = display_group, y = gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL,
       title = "HLA-I-related genes in P1/P2/P4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(
  p_hla, local_figure_dir, "Figure_P1_P2_P4_malignant_epithelial_HLAI_plus_HLAE",
  8.4, 5.4
)

# 【本地 NK 作图】切换到冻结 NK 标签与指定子聚类，统计抑制性受体；其分母不同于前面的上皮图。
# Local P1/P2/P4 NK inhibitory-receptor figure using the frozen local NK annotation.
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRD1")
# 【缺失基因】先确认哪些目标不在矩阵中；本批数据预期缺少 KIR2DL2，因此本地图只用可用受体。
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
missing_nk_receptors <- setdiff(nk_receptors, rownames(counts))
stop_if_not(identical(missing_nk_receptors, "KIR2DL2"),
            paste("Unexpected missing local NK receptors:",
                  paste(missing_nk_receptors, collapse = ", ")))
available_nk_receptors <- setdiff(nk_receptors, missing_nk_receptors)
# 【NK 纳入条件】除 P1/P2/P4 与 NK_cell 外，还要求旧子聚类为 4、7、12；并非把所有 T/NK 候选都画进来。
nk_meta <- meta_all[
  patient_id %in% patients &
    analysis_celltype_previous_or_new == "NK_cell" &
    as.character(previous_tnk_subcluster) %in% c("4", "7", "12")
]
stop_if_not(nrow(nk_meta) == 1129L, "Unexpected strict P1/P2/P4 local NK total")
nk_counts <- counts[, nk_meta$cell_id, drop = FALSE]
nk_meta[, `:=`(
  nk_subtype = paste0("NK_C", as.character(previous_tnk_subcluster)),
  all_gene_raw_umi = as.numeric(Matrix::colSums(nk_counts))
)]
stop_if_not(all(nk_meta$all_gene_raw_umi > 0), "A local NK library size is invalid")

subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^NK_C", "", nk_subtype))]
setorder(subtype_sizes, cluster_number, nk_subtype)
nk_dot <- rbindlist(lapply(available_nk_receptors, function(gene) {
  working <- nk_meta[, .(nk_subtype, all_gene_raw_umi)]
  working[, raw_UMI := as.numeric(nk_counts[gene, ])]
  working[, .(
    NK_cell_n = .N,
    RNA_detected_n = sum(raw_UMI > 0),
    RNA_detected_fraction = mean(raw_UMI > 0),
    total_raw_UMI = sum(raw_UMI),
    mean_logCP10K = mean(log1p(raw_UMI / all_gene_raw_umi * 1e4))
  ), by = nk_subtype][, gene := gene]
}))
nk_dot <- merge(
  nk_dot, subtype_sizes[, .(nk_subtype, NK_cell_n, cluster_number)],
  by = c("nk_subtype", "NK_cell_n"), all.x = TRUE, sort = FALSE
)
nk_unavailable <- copy(subtype_sizes)
nk_unavailable[, `:=`(
  RNA_detected_n = NA_integer_,
  RNA_detected_fraction = NA_real_,
  total_raw_UMI = NA_real_,
  mean_logCP10K = NA_real_,
  gene = "KIR2DL2"
)]
nk_dot <- rbindlist(list(nk_dot, nk_unavailable), use.names = TRUE, fill = TRUE)
nk_dot[, `:=`(
  gene = factor(gene, levels = rev(nk_receptors)),
  subtype_display = paste0(nk_subtype, "\n(n=", NK_cell_n, ")")
)]
subtype_levels <- nk_dot[gene == "KLRC1"][order(cluster_number, nk_subtype), subtype_display]
nk_dot[, subtype_display := factor(subtype_display, levels = subtype_levels)]
fwrite(
  nk_dot,
  file.path(local_table_dir, "Figure_P1_P2_P4_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)

p_nk <- ggplot(nk_dot, aes(subtype_display, gene)) +
  geom_point(
    data = nk_dot[!is.na(RNA_detected_fraction)],
    aes(size = RNA_detected_fraction, color = mean_logCP10K)
  ) +
  geom_point(
    data = nk_dot[is.na(RNA_detected_fraction)],
    shape = 4, color = "#777777", size = 3, stroke = 1
  ) +
  scale_size_continuous(
    name = "% RNA detected", range = c(1.4, 10),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "mean log1p(CP10K)", low = "#e0f3db", high = "#08589e") +
  labs(
    x = "Strict NK subcluster (P1/P2/P4 cell number shown)", y = NULL,
    title = "NK inhibitory receptors in local P1/P2/P4 strict NK cells",
    caption = "KIR2DL2 is absent from the local count feature matrix; × denotes unavailable, not zero."
  ) +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        axis.ticks = element_blank())
fig3_width <- max(10, 0.85 * uniqueN(nk_dot$nk_subtype))
save_figure(
  p_nk, local_figure_dir, "Figure_P1_P2_P4_NK_inhibitory_receptor_subtypes_plus_KIR2DL2",
  fig3_width, 6.5
)

required <- c(
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_STK31_epithelial_localization.", c("png", "pdf")
  )),
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_malignant_epithelial_HLAI_plus_HLAE.", c("png", "pdf")
  )),
  file.path(local_table_dir, "Figure_P1_P2_P4_STK31_epithelial_stats.csv"),
  file.path(local_table_dir, "Figure_P1_P2_P4_malignant_HLAI_plus_HLAE_stats.csv"),
  file.path(local_figure_dir, paste0(
    "Figure_P1_P2_P4_NK_inhibitory_receptor_subtypes_plus_KIR2DL2.", c("png", "pdf")
  )),
  file.path(local_table_dir, "Figure_P1_P2_P4_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "An updated figure or table is absent or empty")
message("LOCAL_P1_P2_P4_STRICT_NK_HLAE_AND_KIR2DL2_FIGURES_PASS")
print(stk)
print(hla)
print(nk_dot[gene == "KIR2DL2"])
