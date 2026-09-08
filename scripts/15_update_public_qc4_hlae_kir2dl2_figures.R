#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】为公共 QC4 图补充 HLA-E 与 KIR2DL2
# 输入：GBC_STK31_public_scRNA 中的目标计数、细胞/样本信息以及 NK.RDS。
# 流程：按 QC4 样本选择细胞 → 计算恶性上皮 HLA 指标 → 计算 NK/NKT 亚型受体指标 → 更新图表。
# 输出：11_figures/current_paper_minimal/paired_QC4 与对应 12_tables 目录。
# HLA-E 作为额外目标基因，KIR2DL2 作为额外 NK 受体；可用基因集合先与实际矩阵取交集。
# 基因缺失与基因存在但计数为零需分开理解；图只描述所选样本与亚型的表达。
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
  library(ggplot2)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
figure_dir <- file.path(project, "11_figures", "current_paper_minimal", "paired_QC4")
table_dir <- file.path(project, "12_tables", "current_paper_minimal", "paired_QC4")
# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
# 【函数：save_figure】将图写入相应结果目录，通常同时生成 PNG 和 PDF。
# stem 是不含扩展名的文件名；绘图数据在主流程中另外导出。
save_figure <- function(plot, stem, width, height) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

qc4 <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")

# 【上皮面板】在选定 QC4 样本的恶性上皮中，统计 HLA-I 相关基因并补充 HLA-E。
# QC4 malignant epithelial HLA-I figure with HLA-E.
hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M")
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
target_counts <- readRDS(file.path(
  project, "02_processed_data", "malignant_epithelial", "target_counts",
  "malignant_target_raw_counts.rds"
))
# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
mal_meta <- fread(file.path(
  project, "02_processed_data", "malignant_epithelial",
  "malignant_cell_metadata.csv.gz"
), showProgress = FALSE)
mal_library <- fread(file.path(
  project, "02_processed_data", "malignant_epithelial",
  "malignant_cell_library_sizes.tsv"
), showProgress = FALSE)
stop_if_not(identical(colnames(target_counts), mal_meta$barcode),
            "Malignant target matrix and metadata order differ")
stop_if_not(
  identical(as.integer(mal_library$malignant_col), as.integer(mal_meta$malignant_col)),
  "Malignant library rows differ from metadata"
)
stop_if_not(all(hla_genes %in% rownames(target_counts)),
            "A requested QC4 HLA-I target is absent")
mal_idx <- which(mal_meta$sample_name %in% qc4)
stop_if_not(length(mal_idx) == 13368L,
            "QC4 malignant epithelial cell total must be 13,368")
qc_target <- target_counts[, mal_idx, drop = FALSE]
qc_library <- mal_library$all_gene_raw_umi[mal_idx]
# 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
hla_stats <- rbindlist(lapply(hla_genes, function(gene) {
  values <- as.numeric(qc_target[gene, ])
  # 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
  # := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
  data.table(
    gene = gene,
    malignant_cell_n = length(values),
    RNA_detected_n = sum(values > 0),
    RNA_detected_fraction = mean(values > 0),
    total_raw_UMI = sum(values),
    pseudobulk_CPM = sum(values) / sum(qc_library) * 1e6
  )
}))
hla_stats[, `:=`(
  log1p_pseudobulk_CPM = log1p(pseudobulk_CPM),
  # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
  gene = factor(gene, levels = rev(hla_genes))
)]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(
  hla_stats,
  file.path(table_dir, "QC4_Figure2_malignant_HLAI_stats_plus_HLAE.csv")
)

# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p_hla <- ggplot(hla_stats, aes("Malignant epithelial cells", gene)) +
  geom_point(aes(size = RNA_detected_fraction, color = log1p_pseudobulk_CPM)) +
  scale_size_continuous(
    name = "% RNA detected", range = c(5, 16),
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  scale_color_gradient(name = "log1p pseudobulk CPM", low = "#fee8c8", high = "#b30000") +
  labs(x = NULL, y = NULL,
       title = "HLA-I-related genes in paired-QC4 malignant epithelial cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0), axis.ticks = element_blank())
save_figure(
  p_hla, "QC4_Figure2_malignant_epithelial_HLAI_plus_HLAE",
  6.8, 5.4
)

rm(target_counts, qc_target)
gc()

# 【NK/NKT 面板】改用 NK 对象及其亚型，补充 KIR2DL2；这一面板的细胞分母与上皮面板不同。
# QC4 public NK/NKT subtype figure with KIR2DL2.
nk_receptors <- c("KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRD1")
nk <- readRDS(file.path(project, "01_raw_data", "NK.RDS"))
stop_if_not(inherits(nk, "Seurat") && DefaultAssay(nk) == "RNA",
            "NK.RDS is not the expected RNA Seurat object")
nk_counts_all <- LayerData(nk[["RNA"]], layer = "counts")
nk_meta_all <- as.data.table(nk[[]], keep.rownames = "barcode")
stop_if_not(identical(colnames(nk_counts_all), nk_meta_all$barcode),
            "NK counts and metadata order differ")
nk_idx <- which(nk_meta_all$orig.ident %in% qc4)
stop_if_not(length(nk_idx) == 1135L, "QC4 NK cell total must be 1,135")
nk_counts <- nk_counts_all[, nk_idx, drop = FALSE]
nk_meta <- nk_meta_all[nk_idx]
nk_meta[, `:=`(
  nk_subtype = as.character(celltype),
  all_gene_raw_umi = as.numeric(Matrix::colSums(nk_counts))
)]
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
missing_nk_receptors <- setdiff(nk_receptors, rownames(nk_counts))
stop_if_not(identical(missing_nk_receptors, "KIR2DL2"),
            paste("Unexpected missing QC4 NK receptors:",
                  paste(missing_nk_receptors, collapse = ", ")))
available_nk_receptors <- setdiff(nk_receptors, missing_nk_receptors)
stop_if_not(all(nk_meta$all_gene_raw_umi > 0), "QC4 NK library sizes are invalid")

subtype_sizes <- nk_meta[, .(NK_cell_n = .N), by = nk_subtype]
subtype_sizes[, cluster_number := as.integer(sub("^.*_C([0-9]+)_.*$", "\\1", nk_subtype))]
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
  file.path(table_dir, "QC4_Figure3_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
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
    x = "Author NK/NKT subtype (QC4 cell number shown)", y = NULL,
    title = "NK inhibitory receptors in paired-QC4 primary samples",
    caption = "KIR2DL2 is absent from the QC4 NK count feature matrix; × denotes unavailable, not zero."
  ) +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        axis.ticks = element_blank())
fig3_width <- max(10, 0.85 * uniqueN(nk_dot$nk_subtype))
save_figure(
  p_nk, "QC4_Figure3_NK_inhibitory_receptor_subtypes_plus_KIR2DL2",
  fig3_width, 6.5
)

required <- c(
  file.path(figure_dir, paste0(
    "QC4_Figure2_malignant_epithelial_HLAI_plus_HLAE.", c("png", "pdf")
  )),
  file.path(figure_dir, paste0(
    "QC4_Figure3_NK_inhibitory_receptor_subtypes_plus_KIR2DL2.", c("png", "pdf")
  )),
  file.path(table_dir, "QC4_Figure2_malignant_HLAI_stats_plus_HLAE.csv"),
  file.path(table_dir, "QC4_Figure3_NK_receptor_subtype_stats_plus_KIR2DL2.csv")
)
stop_if_not(all(file.exists(required) & file.info(required)$size > 0),
            "An updated public QC4 output is absent or empty")
message("PUBLIC_QC4_HLAE_KIR2DL2_FIGURES_PASS")
print(hla_stats)
print(nk_dot[gene == "KIR2DL2"])
