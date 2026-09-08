# ========================================================================
# 【中文阅读指南】沿用已有结果绘制 STK31 分组 UMAP
# 输入：合并样本与 tissue2 的基础 Seurat 对象、细胞注释表和分析汇总表。
# 流程：读取已有 STK31 阈值 → 在上皮细胞内部拆分高/低组 → 沿用已有 UMAP 坐标绘图。
# 输出：两套 STK31/NK 结果目录中的分组 UMAP PDF；分组对象在内存中使用，本脚本不另存 RDS。
# 阅读重点：add_split_celltype 负责改标签，plot_split_umap 负责画图；这里没有重新计算 UMAP。
# 图中两个群靠近仅表示降维后的表达相似性，不能据此判断组织空间距离或实际接触。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

base_dir <- "/home/zhuweiyu/codex-r"
target_gene <- "STK31"
tumor_epithelial_celltype_label <- "Epithelial"
high_label <- "STK31_high_tumor_epithelial"
low_label <- "STK31_low_tumor_epithelial"

# 【函数：read_summary_value】从已有分析汇总表中读取指定 item 对应的值，供后续恢复分组阈值使用。
read_summary_value <- function(path, item) {
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  summary <- read.csv(path, stringsAsFactors = FALSE)
  value <- summary$value[summary$item == item]
  if (length(value) == 0) stop("Missing summary item: ", item, " in ", path)
  value[[1]]
}

# 【函数：get_gene_expr】取指定基因的表达向量，兼容不同 SeuratObject 版本。
# 向量顺序与对象细胞顺序一致，后面据此生成分组。
get_gene_expr <- function(object, gene) {
  assay <- DefaultAssay(object)
  # 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
  mat <- tryCatch(
    GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) GetAssayData(object, assay = assay, slot = "data")
  )
  if (!gene %in% rownames(mat)) stop(gene, " was not found in assay ", assay)
  as.numeric(mat[gene, ])
}

# 【函数：add_split_celltype】在原有细胞类型标签基础上，仅把上皮拆成 STK31-high 与 low。
# cutoff 来自已有汇总，不在本函数中重新拟合聚类或降维。
add_split_celltype <- function(object, cutoff) {
  if (!"manual_celltype" %in% colnames(object@meta.data)) {
    stop("manual_celltype is missing from object metadata")
  }
  stk31_expr <- get_gene_expr(object, target_gene)
  tumor_epithelial_idx <- as.character(object$manual_celltype) == tumor_epithelial_celltype_label
  if (sum(tumor_epithelial_idx) == 0) {
    stop("No tumor epithelial cells found with manual_celltype == ", tumor_epithelial_celltype_label)
  }

  group <- rep("Other", length(stk31_expr))
  # 【严格不等号】本绘图脚本以表达 > cutoff 为 high、<=cutoff 为 low；与旧主分析阈值>0时的 >= 写法不同。
  # 恰好等于阈值的细胞可能分组不同，这是现有代码行为，此次仅作注释说明。
  group[tumor_epithelial_idx & stk31_expr > cutoff] <- high_label
  group[tumor_epithelial_idx & stk31_expr <= cutoff] <- low_label
  names(group) <- colnames(object)

  object$tumor_epithelial_stk31_group <- group
  object$manual_celltype_stk31_tumor_epithelial_split <- as.character(object$manual_celltype)
  object$manual_celltype_stk31_tumor_epithelial_split[group == high_label] <- high_label
  object$manual_celltype_stk31_tumor_epithelial_split[group == low_label] <- low_label
  # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
  object$manual_celltype_stk31_tumor_epithelial_split <- factor(object$manual_celltype_stk31_tumor_epithelial_split)
  object
}

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(path, plot, width = 8, height = 6) {
  pdf(path, width = width, height = height)
  print(plot)
  dev.off()
}

# 【函数：plot_split_umap】使用对象已有 UMAP 坐标，按拆分后的细胞类型上色并输出图。
# 颜色表示类别；坐标轴不是基因表达量，也不是组织中的物理位置。
plot_split_umap <- function(object, path, title) {
  write_plot(
    path,
    # 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
    # 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
    DimPlot(
      object,
      reduction = "umap",
      group.by = "manual_celltype_stk31_tumor_epithelial_split",
      label = TRUE,
      repel = TRUE
    ) + ggtitle(title)
  )
  message("Wrote ", path)
}

message("Creating merged manual celltype UMAP with STK31 tumor epithelial split")
merged_out_dir <- file.path(base_dir, "results", "merged_stk31_nk_analysis")
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
merged_obj <- readRDS(file.path(base_dir, "results", "merged_basic_seurat", "gallbladder_cancer_merged_basic_seurat.rds"))
merged_annotation <- read.csv(file.path(merged_out_dir, "cluster_celltype_annotation.csv"), stringsAsFactors = FALSE)
cluster_to_manual <- setNames(merged_annotation$manual_celltype, as.character(merged_annotation$cluster))
merged_obj$manual_celltype <- unname(cluster_to_manual[as.character(merged_obj$seurat_clusters)])
merged_obj$manual_celltype[is.na(merged_obj$manual_celltype) | !nzchar(merged_obj$manual_celltype)] <- "Unassigned"
merged_obj$manual_celltype <- factor(merged_obj$manual_celltype)
merged_cutoff <- as.numeric(read_summary_value(file.path(merged_out_dir, "analysis_summary.csv"), "tumor_epithelial_stk31_high_cutoff"))
merged_obj <- add_split_celltype(merged_obj, merged_cutoff)
plot_split_umap(
  merged_obj,
  file.path(merged_out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  "Manual annotation with STK31 high/low tumor epithelial split (merged)"
)

message("Creating tissue2 manual celltype UMAP with STK31 tumor epithelial split")
tissue2_out_dir <- file.path(base_dir, "results", "tissue2_stk31_nk_analysis")
tissue2_obj <- readRDS(file.path(tissue2_out_dir, "tissue2_stk31_nk_annotated_seurat.rds"))
tissue2_cutoff <- as.numeric(read_summary_value(file.path(tissue2_out_dir, "analysis_summary.csv"), "tumor_epithelial_stk31_high_cutoff"))
tissue2_obj <- add_split_celltype(tissue2_obj, tissue2_cutoff)
plot_split_umap(
  tissue2_obj,
  file.path(tissue2_out_dir, "umap_manual_celltype_with_stk31_tumor_epithelial_split.pdf"),
  "Manual annotation with STK31 high/low tumor epithelial split"
)

validation_path <- file.path(tissue2_out_dir, "validation", "tissue2_stk31_nk_validated_seurat.rds")
if (file.exists(validation_path)) {
  message("Creating validation UMAP with STK31 tumor epithelial split")
  validation_obj <- readRDS(validation_path)
  validation_obj <- add_split_celltype(validation_obj, tissue2_cutoff)
  validation_obj$validation_celltype_stk31_tumor_epithelial_split <- validation_obj$manual_celltype_stk31_tumor_epithelial_split
  write_plot(
    file.path(tissue2_out_dir, "validation", "umap_validation_celltypes_with_stk31_tumor_epithelial_split.pdf"),
    DimPlot(
      validation_obj,
      reduction = "umap",
      group.by = "validation_celltype_stk31_tumor_epithelial_split",
      label = TRUE,
      repel = TRUE
    ) + ggtitle("Validation cell types with STK31 high/low tumor epithelial split")
  )
}

message("Done")
