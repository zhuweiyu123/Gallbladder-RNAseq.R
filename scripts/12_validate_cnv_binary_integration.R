#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】独立回读检查 CNV 整合结果
# 输入：脚本 11 保存的对象、结果文件和 pseudobulk_sample_metadata.tsv。
# 流程：重新读盘 → 检查细胞数/唯一条形码/注释映射/必需文件 → 检查配对患者与伪批量数。
# 输出：tables/independent_readback_validation.tsv；全部通过时打印 INDEPENDENT_READBACK_PASS。
# expected 中的 27108、1047、117 等值是这一批数据的验收基准；换数据时应重新核对基准。
# 检查通过说明符合既定数据规则，不代表 CNV 分类或生物学假设已经得到实验验证。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

root <- "/home/zhuweiyu/codex-r/results/merged_cnv_binary_annotation"
rds_path <- file.path(root, "merged_cnv_binary_annotated.rds")
# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(rds_path)
# 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
# := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
meta <- as.data.table(obj[[]], keep.rownames = "barcode")
malignant_statuses <- c("malignant_high_confidence", "malignant_likely")
normal_statuses <- c("cnv_normal_like", "uncertain", "unresolved_low_quality")
is_epi <- as.character(meta$celltype) == "Epithelial"
is_malignant <- meta$epithelial_binary_label %in% "Malignant epithelial cells"
is_normal <- meta$epithelial_binary_label %in% "Normal epithelial cells"

# 【文件完整性】列出脚本 11 应生成的文件，后面同时检查存在且大小>0。
required_files <- c(
  "tables/cell_level_annotation_with_cnv_binary.csv.gz",
  "tables/epithelial_binary_counts_by_patient_status.tsv",
  "tables/final_celltype_counts_by_patient.tsv",
  "tables/annotation_transfer_summary.tsv",
  "differential_expression/markers_malignant_vs_normal_epithelial_celllevel_wilcox.csv",
  "differential_expression/markers_stk31_high_vs_low_malignant_epithelial.csv",
  "differential_expression/markers_stk31_high_malignant_epithelial_vs_refined_nk.csv",
  "differential_expression/markers_stk31_high_malignant_epithelial_vs_all_other.csv",
  "differential_expression/pseudobulk_edgeR_paired_malignant_vs_normal.csv",
  "differential_expression/pseudobulk_sample_metadata.tsv",
  "figures/UMAP_final_celltype_cnv_binary.png",
  "figures/UMAP_epithelial_binary_by_patient.png",
  "figures/UMAP_STK31_high_low_malignant_epithelial.png",
  "figures/Volcano_pseudobulk_malignant_vs_normal_epithelial.png",
  "README.md", "sessionInfo.txt"
)

# 【回读验收表】observed 是读盘得到的实际值，expected 是当前数据的已确认基准。
# 每项 pass 都为 TRUE 才继续；不要为了让检查通过而直接改 expected。
checks <- data.table(
  check = c(
    "rds_cells", "barcode_unique", "umap_present", "old_annotation_matched",
    "cnv_matched", "epithelial_total", "malignant_total", "normal_total",
    "malignant_status_mapping_exact", "normal_status_mapping_exact",
    "missing_cnv_epithelial_normal", "missing_cnv_source_explicit",
    "non_epithelial_not_binary_epithelial", "no_unsplit_epithelial_final_label",
    "old_epithelial_conflicts_resolved", "stk31_high_malignant",
    "stk31_low_malignant", "all_required_files_present_nonempty"
  ),
  observed = as.character(c(
    ncol(obj), !anyDuplicated(meta$barcode), "umap" %in% names(obj@reductions),
    sum(meta$previous_annotation_matched), sum(meta$cnv_metadata_matched),
    sum(is_epi), sum(is_malignant, na.rm = TRUE), sum(is_normal, na.rm = TRUE),
    all(meta$cnv_consensus_status[is_malignant] %in% malignant_statuses),
    all(meta$cnv_consensus_status[is_normal & !is.na(meta$cnv_consensus_status)] %in% normal_statuses),
    sum(is_epi & is.na(meta$cnv_consensus_status) & is_normal, na.rm = TRUE),
    all(meta$analysis_celltype_cnv_binary_source[
      is_epi & is.na(meta$cnv_consensus_status)
    ] == "new_epithelial_without_cnv_assigned_normal_by_user_rule"),
    all(!meta$analysis_celltype_cnv_binary[!is_epi] %in%
          c("Malignant epithelial cells", "Normal epithelial cells")),
    !any(meta$analysis_celltype_cnv_binary == "Epithelial"),
    sum(meta$analysis_celltype_cnv_binary_source ==
          "new_annotation_resolved_old_epithelial_conflict"),
    sum(meta$stk31_malignant_group == "STK31-high malignant epithelial", na.rm = TRUE),
    sum(meta$stk31_malignant_group == "STK31-low malignant epithelial", na.rm = TRUE),
    all(file.exists(file.path(root, required_files)) & file.info(file.path(root, required_files))$size > 0)
  )),
  expected = as.character(c(
    27108, TRUE, TRUE, 25821, 4304, 4763, 1047, 3716,
    TRUE, TRUE, 459, TRUE, TRUE, TRUE, 89, 117, 930, TRUE
  ))
)
checks[, pass := observed == expected]
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
fwrite(checks, file.path(root, "tables", "independent_readback_validation.tsv"), sep = "\t")
print(checks)
stop_if_not(all(checks$pass), "Independent read-back validation failed")

# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
pb_meta <- fread(file.path(root, "differential_expression", "pseudobulk_sample_metadata.tsv"))
stop_if_not(
  identical(sort(unique(as.character(pb_meta$patient_id))), c("P1", "P2", "P4")),
  "Unexpected pseudobulk paired patients"
)
stop_if_not(nrow(pb_meta) == 6L, "Unexpected pseudobulk sample count")

message("INDEPENDENT_READBACK_PASS")
