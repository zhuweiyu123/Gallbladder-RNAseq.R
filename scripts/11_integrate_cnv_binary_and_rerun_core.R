#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】整合 CNV 二分类注释并重跑核心分析
# 输入：新的 all_local_annotated.rds、旧逐细胞冻结注释、local_epithelial_cnv_metadata.csv.gz。
# 流程：转换样本前缀并按条形码匹配 → 转移旧注释 → 合并 CNV 状态 → 恶性/正常二分类 → 差异与伪批量。
# 输出：results/merged_cnv_binary_annotation 下的整合对象、注释表、差异表、图和运行说明。
# 注意映射：tissue1→P1、tissue2→P2、tissue4→P3、tissue5→P4；不能把 tissue4 直接当作 P4。
# 本项目规则将 malignant_high_confidence/likely 归恶性，其余指定状态和缺失 CNV 的上皮归 Normal。
# Normal 在这里是工作分组名称，包含 uncertain、低质量和缺少 CNV 的细胞，并非全都已证实正常。
# 恶性上皮内 STK31 high=counts>0，low=counts==0；这里不再使用旧流程的 75% 分位数定义。
# 同时输出细胞级 Wilcoxon 和患者配对 edgeR 伪批量结果；后者按患者×组别汇总原始计数。
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
  library(edgeR)
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260821)

project <- "/home/zhuweiyu/codex-r"
cnv_project <- file.path(project, "GBC_STK31_public_scRNA", "07_local_gbc_cnv")
new_object_path <- file.path(
  cnv_project, "03_celltype_annotation", "objects", "all_local_annotated.rds"
)
old_annotation_path <- file.path(
  project, "results", "merged_tnk_refined_annotation", "refined_celltype_annotation_by_cell.csv"
)
cnv_metadata_path <- file.path(
  cnv_project, "07_consensus_annotation", "local_epithelial_cnv_metadata.csv.gz"
)
out_dir <- file.path(project, "results", "merged_cnv_binary_annotation")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
de_dir <- file.path(out_dir, "differential_expression")
for (path in c(out_dir, figure_dir, table_dir, de_dir)) {
  # 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
# 【函数：save_plot】同一绘图对象同时保存为 PNG 和 PDF，文件名前缀由 stem 决定。
# PNG 便于预览，PDF 便于排版；width/height 控制输出尺寸。
save_plot <- function(plot, stem, width, height) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

message("Reading new QC/annotation object")
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(new_object_path)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  obj <- JoinLayers(obj, assay = "RNA")
}
stop_if_not(ncol(obj) == 27108L, "Unexpected new-object cell count")
stop_if_not(!anyDuplicated(colnames(obj)), "New-object cell names are not unique")
stop_if_not("celltype" %in% colnames(obj[[]]), "New object lacks celltype")
stop_if_not("patient_id" %in% colnames(obj[[]]), "New object lacks patient_id")

message("Reading frozen previous annotation and CNV consensus")
# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
old <- fread(old_annotation_path)
cnv <- fread(cnv_metadata_path)
stop_if_not(nrow(old) == 28751L && !anyDuplicated(old$cell), "Unexpected old annotation")
stop_if_not(nrow(cnv) == 4304L && !anyDuplicated(cnv$barcode), "Unexpected CNV metadata")

# 【样本名对齐】旧条形码带 tissue 前缀，新对象用 P 前缀；先转换名称，再用 match 按名字匹配。
prefix_map <- c(tissue1 = "P1", tissue2 = "P2", tissue4 = "P3", tissue5 = "P4")
old[, old_sample_prefix := sub("_.*$", "", cell)]
old[, barcode_raw := sub("^[^_]+_", "", cell)]
old[, new_barcode := paste0(prefix_map[old_sample_prefix], "_", barcode_raw)]
stop_if_not(!anyNA(old$new_barcode), "Old-to-new prefix conversion failed")
stop_if_not(!anyDuplicated(old$new_barcode), "Converted old barcodes are not unique")

new_cells <- colnames(obj)
# 【关键对齐】match 返回新细胞在旧表的位置；未匹配是 NA，不能直接依靠两张表的行号拼接。
# 【按名称对齐】match(x,y) 返回 x 各元素在 y 中的位置，没找到返回 NA；用来保证标签和表达对应同一细胞。
old_index <- match(new_cells, old$new_barcode)
cnv_index <- match(new_cells, cnv$barcode)
matched_old <- !is.na(old_index)
matched_cnv <- !is.na(cnv_index)
stop_if_not(sum(matched_old) == 25821L, "Unexpected old/new overlap")
stop_if_not(sum(matched_cnv) == 4304L, "All CNV epithelial cells must map to new object")

previous_columns <- c(
  "legacy_broad_celltype", "legacy_manual_celltype", "tnk_subcluster",
  "refined_celltype", "analysis_celltype"
)
for (column in previous_columns) {
  values <- rep(NA_character_, length(new_cells))
  values[matched_old] <- as.character(old[[column]][old_index[matched_old]])
  obj[[paste0("previous_", column)]] <- values
}
obj$previous_annotation_matched <- matched_old

fallback_map <- c(
  "Epithelial" = "Epithelial",
  "T/NK" = "T_NK_unresolved_new",
  "B/Plasma" = "B_Plasma_unresolved_new",
  "Myeloid" = "Myeloid",
  "Endothelial" = "Endothelial",
  "Fibroblast/Stromal" = "Fibroblast_Stromal",
  "Mast" = "Mast_cell",
  "Ambiguous" = "Ambiguous"
)
base_celltype <- as.character(obj$previous_analysis_celltype)
fallback <- unname(fallback_map[as.character(obj$celltype)])
fallback[is.na(fallback)] <- paste0("New_", as.character(obj$celltype[is.na(fallback)]))
base_celltype[!matched_old] <- fallback[!matched_old]
# Resolve old unsplit epithelial labels using the new object as the annotation base.
old_epithelial_conflict <- matched_old &
  base_celltype == "Epithelial" & as.character(obj$celltype) != "Epithelial"
base_celltype[old_epithelial_conflict] <- fallback[old_epithelial_conflict]
stop_if_not(!anyNA(base_celltype), "Base celltype contains missing values")
obj$analysis_celltype_previous_or_new <- base_celltype

cnv_columns <- c(
  "cnv_consensus_status", "external_reference_cnv_status",
  "internal_reference_cnv_status", "cnv_leiden_external", "cnv_score_external",
  "cnv_burden_cell_external", "cnv_leiden_internal", "cnv_score_internal",
  "cnv_burden_cell_internal", "infercnvpy_version", "cnv_parameter_set"
)
for (column in cnv_columns) {
  if (is.numeric(cnv[[column]]) || is.integer(cnv[[column]])) {
    values <- rep(NA_real_, length(new_cells))
  } else {
    values <- rep(NA_character_, length(new_cells))
  }
  values[matched_cnv] <- cnv[[column]][cnv_index[matched_cnv]]
  obj[[column]] <- values
}
obj$cnv_metadata_matched <- matched_cnv

is_new_epithelial <- as.character(obj$celltype) == "Epithelial"
stop_if_not(sum(is_new_epithelial) == 4763L, "Unexpected new epithelial cell count")
stop_if_not(sum(is_new_epithelial & matched_cnv) == 4304L,
            "Unexpected CNV-assessed epithelial count")

# 【项目二分类规则】high_confidence/likely 归恶性；Normal 还含 uncertain、低质量、缺少 CNV 的上皮。
# 这里复现已有分组约定，Normal 标签不表示每个细胞都已经证实非恶性。
malignant_statuses <- c("malignant_high_confidence", "malignant_likely")
normal_statuses <- c("cnv_normal_like", "uncertain", "unresolved_low_quality")
binary_label <- rep(NA_character_, length(new_cells))
binary_label[is_new_epithelial & obj$cnv_consensus_status %in% malignant_statuses] <-
  "Malignant epithelial cells"
binary_label[is_new_epithelial & obj$cnv_consensus_status %in% normal_statuses] <-
  "Normal epithelial cells"
# User-specified fallback: epithelial cells without a CNV result are assigned normal.
binary_label[is_new_epithelial & is.na(obj$cnv_consensus_status)] <-
  "Normal epithelial cells"
stop_if_not(!anyNA(binary_label[is_new_epithelial]), "Epithelial binary labelling failed")
stop_if_not(sum(binary_label == "Malignant epithelial cells", na.rm = TRUE) == 1047L,
            "Unexpected malignant epithelial count")
stop_if_not(sum(binary_label == "Normal epithelial cells", na.rm = TRUE) == 3716L,
            "Unexpected normal epithelial count")
obj$epithelial_binary_label <- binary_label

final_celltype <- base_celltype
final_celltype[is_new_epithelial] <- binary_label[is_new_epithelial]
obj$analysis_celltype_cnv_binary <- final_celltype

label_source <- ifelse(
  matched_old, "previous_frozen_analysis_celltype", "new_annotation_fallback"
)
label_source[old_epithelial_conflict] <-
  "new_annotation_resolved_old_epithelial_conflict"
label_source[is_new_epithelial & matched_cnv] <- "dual_reference_cnv_user_binary_rule"
label_source[is_new_epithelial & !matched_cnv] <-
  "new_epithelial_without_cnv_assigned_normal_by_user_rule"
obj$analysis_celltype_cnv_binary_source <- label_source

# 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
# := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
mapping <- as.data.table(obj[[]], keep.rownames = "cell_id")
mapping[, old_cell := ifelse(matched_old, old$cell[old_index], NA_character_)]
mapping[, old_annotation_available := matched_old]
mapping[, cnv_annotation_available := matched_cnv]
mapping[, patient_id := as.character(patient_id)]
mapping[, new_celltype := as.character(celltype)]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(mapping, file.path(table_dir, "cell_level_annotation_with_cnv_binary.csv.gz"))

count_patient_binary <- mapping[
  new_celltype == "Epithelial",
  .N,
  by = .(patient_id, cnv_consensus_status, epithelial_binary_label)
][order(patient_id, epithelial_binary_label, cnv_consensus_status)]
fwrite(count_patient_binary, file.path(table_dir, "epithelial_binary_counts_by_patient_status.tsv"), sep = "\t")

count_final_celltype <- mapping[, .N, by = .(patient_id, analysis_celltype_cnv_binary)][
  order(patient_id, analysis_celltype_cnv_binary)
]
fwrite(count_final_celltype, file.path(table_dir, "final_celltype_counts_by_patient.tsv"), sep = "\t")

transfer_summary <- data.table(
  metric = c(
    "new_object_cells", "old_annotation_cells", "old_new_overlap",
    "new_cells_without_old_annotation", "old_cells_not_in_new_object",
    "new_epithelial_cells", "cnv_assessed_epithelial_cells",
    "epithelial_without_cnv_assigned_normal", "malignant_epithelial_cells",
    "normal_epithelial_cells", "old_epithelial_conflicts_resolved_from_new_annotation"
  ),
  value = c(
    ncol(obj), nrow(old), sum(matched_old), sum(!matched_old),
    sum(!old$new_barcode %in% new_cells), sum(is_new_epithelial),
    sum(is_new_epithelial & matched_cnv), sum(is_new_epithelial & !matched_cnv),
    sum(binary_label == "Malignant epithelial cells", na.rm = TRUE),
    sum(binary_label == "Normal epithelial cells", na.rm = TRUE),
    sum(old_epithelial_conflict)
  )
)
fwrite(transfer_summary, file.path(table_dir, "annotation_transfer_summary.tsv"), sep = "\t")

output_rds <- file.path(out_dir, "merged_cnv_binary_annotated.rds")

message("Generating annotation figures")
if ("umap" %in% names(obj@reductions)) {
  # 【读图】DimPlot 按类别上色，FeaturePlot 按连续表达上色；DotPlot 的点大小通常是检出比例、颜色是平均表达。
  # 若使用了缩放，颜色表示相对值；本次图的具体分组以 group.by/Idents 为准。
  p_all <- DimPlot(
    obj, reduction = "umap", group.by = "analysis_celltype_cnv_binary",
    label = TRUE, repel = TRUE, raster = TRUE
  ) + ggtitle("Final cell types with CNV-defined epithelial binary labels") +
    theme_classic(base_size = 11)
  save_plot(p_all, "UMAP_final_celltype_cnv_binary", 11, 8)

  # 【取细胞子集】按 cells 或条件保留需要的细胞；这一步改变本次分析范围，要同步核对后面的分母。
  epi_obj <- subset(obj, cells = new_cells[is_new_epithelial])
  p_epi <- DimPlot(
    epi_obj, reduction = "umap", group.by = "epithelial_binary_label",
    split.by = "patient_id", ncol = 2, raster = TRUE,
    cols = c("Malignant epithelial cells" = "#b2182b", "Normal epithelial cells" = "#2166ac")
  ) + ggtitle("CNV-derived epithelial binary label by patient")
  save_plot(p_epi, "UMAP_epithelial_binary_by_patient", 12, 9)
}

message("Preparing expression groups")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "Joined RNA counts are not sparse")
stop_if_not("STK31" %in% rownames(counts), "STK31 missing")
stk31_counts <- as.numeric(counts["STK31", ])
malignant_cells <- new_cells[obj$epithelial_binary_label %in% "Malignant epithelial cells"]
normal_cells <- new_cells[obj$epithelial_binary_label %in% "Normal epithelial cells"]
# 【本版 STK31 定义】恶性上皮内 counts>0 叫 high、counts==0 叫 low；不是旧脚本中的分位数分组。
stk31_high_malignant <- new_cells[
  obj$epithelial_binary_label %in% "Malignant epithelial cells" & stk31_counts > 0
]
stk31_low_malignant <- new_cells[
  obj$epithelial_binary_label %in% "Malignant epithelial cells" & stk31_counts == 0
]
nk_cells <- new_cells[obj$analysis_celltype_previous_or_new == "NK_cell"]
# 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
all_other_than_high <- setdiff(new_cells, stk31_high_malignant)

de_group_summary <- data.table(
  group = c(
    "malignant_epithelial", "normal_epithelial", "stk31_high_malignant",
    "stk31_low_malignant", "previous_high_confidence_NK", "all_other_than_stk31_high"
  ),
  cells = c(
    length(malignant_cells), length(normal_cells), length(stk31_high_malignant),
    length(stk31_low_malignant), length(nk_cells), length(all_other_than_high)
  )
)
fwrite(de_group_summary, file.path(de_dir, "de_group_counts.tsv"), sep = "\t")

# 【函数：run_de】为两组细胞建立临时标签，运行 FindMarkers，排序并导出差异基因表。
# 第一组/第二组的传入顺序决定 logFC 正负；每组抽样上限见 max.cells.per.ident。
run_de <- function(object, cells_1, cells_2, comparison) {
  stop_if_not(length(cells_1) >= 10L && length(cells_2) >= 10L,
              paste("Too few cells for", comparison))
  object$.de_group <- "unused"
  object$.de_group[colnames(object) %in% cells_1] <- "group1"
  object$.de_group[colnames(object) %in% cells_2] <- "group2"
  message("DE ", comparison, ": ", length(cells_1), " vs ", length(cells_2))
  # 【差异表达】比较 ident.1 与 ident.2；avg_log2FC>0 表示第一组更高，p_val_adj 是多重检验调整后 P 值。
  # min.pct/logfc.threshold 是进入检验的筛选条件；细胞级检验不自动控制患者内相关性。
  result <- FindMarkers(
    object, ident.1 = "group1", ident.2 = "group2", group.by = ".de_group",
    logfc.threshold = 0.1, min.pct = 0.1, test.use = "wilcox",
    max.cells.per.ident = 500, random.seed = 20260821, verbose = FALSE
  )
  result$gene <- rownames(result)
  result <- result[order(result$p_val_adj, -abs(result$avg_log2FC)), ]
  fwrite(as.data.table(result), file.path(de_dir, paste0(comparison, ".csv")))
  fwrite(as.data.table(head(result[result$avg_log2FC > 0, ], 50)),
         file.path(de_dir, paste0(comparison, "_top50_group1_up.csv")))
  fwrite(as.data.table(head(result[result$avg_log2FC < 0, ], 50)),
         file.path(de_dir, paste0(comparison, "_top50_group1_down.csv")))
  result
}

de_malignant_normal <- run_de(
  obj, malignant_cells, normal_cells,
  "markers_malignant_vs_normal_epithelial_celllevel_wilcox"
)
de_high_low <- run_de(
  obj, stk31_high_malignant, stk31_low_malignant,
  "markers_stk31_high_vs_low_malignant_epithelial"
)
de_high_nk <- run_de(
  obj, stk31_high_malignant, nk_cells,
  "markers_stk31_high_malignant_epithelial_vs_refined_nk"
)
de_high_all <- run_de(
  obj, stk31_high_malignant, all_other_than_high,
  "markers_stk31_high_malignant_epithelial_vs_all_other"
)

# 【患者级伪批量】将同一患者同一组的原始 counts 相加，再比较配对的恶性/Normal 组。
# 同一患者贡献一对汇总样本，能够在设计矩阵中控制患者间基线差异。
message("Running patient-paired epithelial pseudobulk edgeR")
epi_index <- which(is_new_epithelial)
epi_patient <- as.character(obj$patient_id[epi_index])
epi_binary_group <- ifelse(
  obj$epithelial_binary_label[epi_index] == "Malignant epithelial cells",
  "Malignant", "Normal"
)
pair_availability <- as.data.table(table(
  patient_id = epi_patient,
  epithelial_group = epi_binary_group
))
fwrite(pair_availability, file.path(de_dir, "pseudobulk_pair_availability.tsv"), sep = "\t")
availability_wide <- dcast(pair_availability, patient_id ~ epithelial_group, value.var = "N", fill = 0)
# 【选择可配对患者】必须同时有恶性与 Normal 细胞；本批次对应 P1、P2、P4。
paired_patients <- availability_wide[
  Malignant > 0 & Normal > 0, as.character(patient_id)
]
stop_if_not(identical(paired_patients, c("P1", "P2", "P4")),
            "Unexpected set of paired patients")
paired_epi_keep <- epi_patient %in% paired_patients
epi_index_paired <- epi_index[paired_epi_keep]
pb_group <- paste0(
  epi_patient[paired_epi_keep], "__", epi_binary_group[paired_epi_keep]
)
pb_levels <- as.vector(sapply(
  paired_patients,
  function(patient) paste0(patient, "__", c("Normal", "Malignant"))
))
# 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
pb_group <- factor(pb_group, levels = pb_levels)
# 【计数聚合】稀疏指示矩阵每列代表患者×组别；counts %*% aggregation 就是按组求和。
# 【设计矩阵】把患者/分组等变量编码为模型列；配对分析需同时考虑患者效应和目标组别效应。
aggregation <- sparse.model.matrix(~ 0 + pb_group)
colnames(aggregation) <- sub("^pb_group", "", colnames(aggregation))
pb_counts <- counts[, epi_index_paired, drop = FALSE] %*% aggregation
stop_if_not(ncol(pb_counts) == 6L, "Expected six paired patient-by-binary pseudobulks")
pb_meta <- data.table(pseudobulk_id = colnames(pb_counts))
pb_meta[, c("patient_id", "epithelial_group") := tstrsplit(pseudobulk_id, "__", fixed = TRUE)]
pb_meta[, epithelial_group := factor(epithelial_group, levels = c("Normal", "Malignant"))]
pb_meta[, patient_id := factor(patient_id, levels = paired_patients)]
stop_if_not(all(table(pb_meta$patient_id, pb_meta$epithelial_group) == 1L),
            "Pseudobulk pairing is incomplete")

# 【edgeR 输入】把按样本汇总的原始 counts 包装成差异分析对象；这里应使用计数而非 log 标准化矩阵。
dge <- DGEList(counts = pb_counts, samples = as.data.frame(pb_meta))
design <- model.matrix(~ patient_id + epithelial_group, data = pb_meta)
# 【低表达过滤】结合样本库大小和设计矩阵剔除信息不足的基因，减少无效检验。
keep_genes <- filterByExpr(dge, design = design)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
# 【库大小校正】计算组成偏差校正因子；常用 TMM 用于使不同伪批量样本的表达更可比。
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
# 【组别检验】在拟合的准似然模型中检验指定系数/对比；logFC 方向由设计矩阵和 coef 决定。
test <- glmQLFTest(fit, coef = "epithelial_groupMalignant")
pb_de <- topTags(test, n = Inf, sort.by = "PValue")$table
pb_de$gene <- rownames(pb_de)
pb_de <- pb_de[, c("gene", setdiff(colnames(pb_de), "gene"))]
fwrite(as.data.table(pb_de), file.path(de_dir, "pseudobulk_edgeR_paired_malignant_vs_normal.csv"))
fwrite(pb_meta, file.path(de_dir, "pseudobulk_sample_metadata.tsv"), sep = "\t")
fwrite(as.data.table(as.matrix(pb_counts), keep.rownames = "gene"),
       file.path(de_dir, "pseudobulk_raw_counts.tsv.gz"), sep = "\t")

pb_cpm <- cpm(dge, log = TRUE, prior.count = 1)
# 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
target_genes <- intersect(c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M"), rownames(pb_cpm))
fwrite(as.data.table(pb_cpm[target_genes, , drop = FALSE], keep.rownames = "gene"),
       file.path(de_dir, "pseudobulk_target_genes_log2CPM.tsv"), sep = "\t")

volcano <- as.data.table(pb_de)
volcano[, significant := FDR < 0.05 & abs(logFC) > 0.25]
volcano[, neg_log10_fdr := -log10(pmax(FDR, 1e-300))]
# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p_volcano <- ggplot(volcano, aes(x = logFC, y = neg_log10_fdr, color = significant)) +
  geom_point(alpha = 0.55, size = 1) +
  scale_color_manual(values = c(`FALSE` = "#bdbdbd", `TRUE` = "#b2182b")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, color = "#666666") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#666666") +
  labs(
    title = "Patient-paired pseudobulk: malignant vs normal epithelial",
    subtitle = "edgeR QL; design ~ patient + epithelial group; positive logFC = malignant higher",
    x = "log2 fold change", y = "-log10 FDR", color = "FDR<0.05 & |logFC|>0.25"
  ) + theme_classic(base_size = 11)
save_plot(p_volcano, "Volcano_pseudobulk_malignant_vs_normal_epithelial", 8, 6)

if ("umap" %in% names(obj@reductions)) {
  stk31_high_flag <- rep("Other", ncol(obj))
  stk31_high_flag[colnames(obj) %in% stk31_low_malignant] <- "STK31-low malignant epithelial"
  stk31_high_flag[colnames(obj) %in% stk31_high_malignant] <- "STK31-high malignant epithelial"
  obj$stk31_malignant_group <- factor(
    stk31_high_flag,
    levels = c("Other", "STK31-low malignant epithelial", "STK31-high malignant epithelial")
  )
  p_relationship <- DimPlot(
    obj, reduction = "umap", group.by = "stk31_malignant_group", raster = TRUE,
    cols = c("#d9d9d9", "#2166ac", "#b2182b")
  ) + ggtitle("STK31-detected and undetected malignant epithelial cells")
  save_plot(p_relationship, "UMAP_STK31_high_low_malignant_epithelial", 9, 7)
}

message("Saving derived annotated object")
# 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
saveRDS(obj, output_rds, compress = TRUE)

validation <- data.table(
  check = c(
    "new_object_cells", "old_new_overlap", "all_cnv_cells_mapped",
    "new_epithelial_cells", "malignant_binary_cells", "normal_binary_cells",
    "cnv_assessed_epithelial", "no_cnv_epithelial_user_normal_fallback",
    "non_epithelial_matched_previous_labels_preserved", "pseudobulk_columns",
    "paired_patient_group_completeness", "no_unsplit_epithelial_final_label",
    "output_rds_exists"
  ),
  observed = as.character(c(
    ncol(obj), sum(matched_old), sum(matched_cnv), sum(is_new_epithelial),
    length(malignant_cells), length(normal_cells), sum(is_new_epithelial & matched_cnv),
    sum(is_new_epithelial & !matched_cnv),
    all(final_celltype[matched_old & !is_new_epithelial] == base_celltype[matched_old & !is_new_epithelial]),
    ncol(pb_counts), all(table(pb_meta$patient_id, pb_meta$epithelial_group) == 1L),
    !any(final_celltype == "Epithelial"),
    file.exists(output_rds)
  )),
  expected = as.character(c(
    27108, 25821, 4304, 4763, 1047, 3716, 4304, 459,
    TRUE, 6, TRUE, TRUE, TRUE
  ))
)
validation[, pass := observed == expected]
fwrite(validation, file.path(table_dir, "validation_checks.tsv"), sep = "\t")
stop_if_not(all(validation$pass), "Final validation failed")

summary_lines <- c(
  "# CNV binary epithelial integration and core rerun",
  "",
  "## User-requested binary rule",
  "",
  "- malignant_high_confidence + malignant_likely -> Malignant epithelial cells",
  "- cnv_normal_like + uncertain + unresolved_low_quality -> Normal epithelial cells",
  "- 459 new-object epithelial cells without a CNV result -> Normal epithelial cells by explicit user fallback rule",
  "- Original CNV status and previous frozen annotation are retained in separate fields.",
  "",
  "## Object and mapping",
  "",
  paste0("- New object cells: ", ncol(obj)),
  paste0("- Previous frozen annotation matched: ", sum(matched_old)),
  paste0("- Old unsplit epithelial conflicts resolved from new annotation: ", sum(old_epithelial_conflict)),
  paste0("- CNV epithelial mapped: ", sum(matched_cnv), "/4304"),
  paste0("- Malignant epithelial cells: ", length(malignant_cells)),
  paste0("- Normal epithelial cells: ", length(normal_cells)),
  paste0("- Previous refined NK retained in new cell universe: ", length(nk_cells)),
  "",
  "## Differential expression",
  "",
  "- Patient-paired pseudobulk malignant vs normal epithelial: edgeR QL, design ~ patient + epithelial group; P1/P2/P4 only because P3 has no malignant cells under the user rule.",
  "- P2 malignant cells are all malignant_likely and parameter-sensitive; the n=3 paired pseudobulk result is exploratory.",
  "- Cell-level Wilcoxon outputs are provided for continuity with the old workflow; max 500 cells per group and must not be treated as patient-level inference.",
  "- STK31-high is raw RNA count > 0 within malignant epithelial cells, matching the previous workflow.",
  "- No GO, CellChat, ligand-receptor or causal analysis was run in this integration step.",
  "",
  "## Main output",
  "",
  paste0("- Seurat object: `", output_rds, "`")
)
writeLines(summary_lines, file.path(out_dir, "README.md"), useBytes = TRUE)

# 【运行记录】记录 R 与已加载包的版本，帮助以后解释同一代码为何可能得到不同结果。
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
message("PASS: CNV binary integration and core rerun complete")
print(transfer_summary)
print(de_group_summary)
print(validation)
