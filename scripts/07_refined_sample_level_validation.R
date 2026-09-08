# ========================================================================
# 【中文阅读指南】以样本为单位检查候选通讯的稳定性
# 输入：脚本 06 的 refined_merged_object.rds；沿用已经确定的细胞身份，不再次修改注释。
# 流程：每个样本分别统计高/低上皮和 NK → 按发送/接收方向评分 → 留一验证 → 细胞重采样 → 分组敏感性。
# 输出：默认 results/refined_sample_level_validation_v2，含各样本候选轴汇总及验证表。
# 这里的独立生物学单位是样本（本项目 n=4）；几千个细胞不等于几千个独立患者。
# 样本内 bootstrap 反映所采细胞的波动，不增加样本数；LOSO 检查去掉一个样本后结论是否改变。
# 基因在矩阵行名中存在，不代表在当前细胞群检测到；代码分别记录 genes_present 与 genes_detected。
# 此处的模块乘积是描述性候选分数，不能作为 CellChat 通讯概率解释。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# ============================================================
# 07_refined_sample_level_validation.R  (method-revised v2)
# Sample-level exploratory validation using FROZEN refined labels.
# Statistical unit = sample (n=4). Exploratory only.
# NK: analysis_celltype == "NK_cell"
#
# PI method fixes:
#  1) Axis-specific sender/receiver directions
#  2) Detection != presence in matrix
#  3) Bootstrap fields renamed; not used to upgrade evidence
#  4) LOSO stability separate from direction consistency
#  5) STK31 sensitivity reports set identity / non-informative
# Default output: results/refined_sample_level_validation_v2/
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
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
out_dir <- Sys.getenv(
  "STK31_SAMPLE_VAL_OUT",
  file.path(project_dir, "results/refined_sample_level_validation_v2")
)
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
epi_label <- "Epithelial"
# 【参数 n_boot】样本内细胞重采样的重复次数，控制稳定性估计的计算量。
n_boot <- as.integer(Sys.getenv("STK31_SAMPLE_BOOT", "200"))
# 【参数 min_high_cells】每样本 high 组的最低细胞数，用于判定该样本是否可评价。
min_high_cells <- 3
# 【参数 min_nk_cells】每样本 NK 的最低细胞数；数量不足时不应强行给出比较方向。
min_nk_cells <- 10
# 【参数 min_low_cells】每样本 low 组的最低细胞数，与 high/NK 门槛共同决定可评价性。
min_low_cells <- 20
# Detection thresholds (explicit)
# 【参数 min_mean_expr】判定基因实际检出时使用的平均表达下限，需与比例条件一起看。
min_mean_expr <- 0
# 【参数 min_pct_expr】实际检出所需的最低阳性细胞比例；0.01 表示 1%。
min_pct_expr <- 0.01
# 【参数 similar_abs】high-low 绝对差小于此阈值标记相近；这是描述性门槛，不是 P 值。
similar_abs <- 0.02

# Direction-aware axes (sender / receiver)
# 【候选轴及方向】ligands 在发送细胞评分，receptors 在接收细胞评分；同一路径反向需另算。
axes <- list(
  "MHC-I/HLA_epi_to_NK" = list(
    sender = "epithelial", receiver = "NK",
    ligands = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M"),
    receptors = c("KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR2DL4")
  ),
  "NECTIN2/PVR_epi_to_NK" = list(
    sender = "epithelial", receiver = "NK",
    ligands = c("NECTIN2", "PVR"),
    receptors = c("TIGIT", "CD96")
  ),
  "MICA/B_epi_to_NK" = list(
    sender = "epithelial", receiver = "NK",
    ligands = c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3"),
    receptors = c("KLRK1")
  ),
  "IFNG_NK_to_epi" = list(
    sender = "NK", receiver = "epithelial",
    ligands = c("IFNG"),
    receptors = c("IFNGR1", "IFNGR2")
  ),
  "TGFb_NK_to_epi" = list(
    sender = "NK", receiver = "epithelial",
    ligands = c("TGFB1", "TGFB2", "TGFB3"),
    receptors = c("TGFBR1", "TGFBR2")
  ),
  "TGFb_epi_to_NK" = list(
    sender = "epithelial", receiver = "NK",
    ligands = c("TGFB1", "TGFB2", "TGFB3"),
    receptors = c("TGFBR1", "TGFBR2")
  )
)

# 【函数：stop_if_missing】先检查输入文件是否存在；缺失时立刻 stop，避免下游产生误导性结果。
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
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

# 【函数：module_score_cells】在指定细胞中计算模块均值、阳性比例及每个基因的检测情况。
# genes_present 只表示矩阵有这一行；genes_detected 还需要满足表达和比例阈值。
module_score_cells <- function(object, genes, cells) {
  genes_present <- available_genes(object, genes)
  if (length(genes_present) == 0 || length(cells) == 0) {
    return(list(
      mean_expr = NA_real_, pct_pos = NA_real_, n_cells = length(cells),
      genes_present = character(0), genes_detected = character(0)
    ))
  }
  mat <- fetch_mat(object, genes_present)
  mat <- mat[, intersect(cells, colnames(mat)), drop = FALSE]
  if (ncol(mat) == 0) {
    return(list(
      mean_expr = NA_real_, pct_pos = NA_real_, n_cells = 0,
      genes_present = genes_present, genes_detected = character(0)
    ))
  }
  # per-gene detection in these cells
  # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
  gene_mean <- as.numeric(Matrix::rowMeans(mat))
  gene_pct <- as.numeric(Matrix::rowMeans(mat > 0))
  names(gene_mean) <- genes_present
  names(gene_pct) <- genes_present
  detected <- genes_present[
    (gene_mean > min_mean_expr) & (gene_pct >= min_pct_expr)
  ]
  # 【按细胞求均值】在基因×细胞矩阵中，colMeans 对每列跨基因求平均，常用来构建逐细胞模块分数。
  cell_score <- as.numeric(Matrix::colMeans(mat))
  list(
    mean_expr = mean(cell_score),
    pct_pos = mean(cell_score > 0),
    n_cells = length(cell_score),
    genes_present = genes_present,
    genes_detected = detected,
    gene_mean = gene_mean,
    gene_pct = gene_pct
  )
}

# 【函数：is_detected_module】综合可用基因、实际检出基因和模块数值，判断该模块能否作为检测到的信号。
# 请同时查看基因级阈值；这里包含或条件，不能简单理解成所有指标都超过同一阈值。
is_detected_module <- function(mod) {
  if (is.null(mod) || length(mod$genes_present) == 0) return(FALSE)
  if (!is.finite(mod$mean_expr) || !is.finite(mod$pct_pos)) return(FALSE)
  # module-level: at least one gene detected AND module mean/pct above floor
  length(mod$genes_detected) > 0 &&
    (mod$mean_expr > min_mean_expr) &&
    (mod$pct_pos >= min_pct_expr || mod$mean_expr > 0)
}

# 【函数：label_direction】根据 high-low 差值标记方向；绝对差小于 similar_abs 记为 similar。
# 非有限值记 insufficient，表示不能评价；方向标签本身不提供统计显著性。
label_direction <- function(delta) {
  if (!is.finite(delta)) return("insufficient")
  if (abs(delta) < similar_abs) return("similar")
  if (delta > 0) return("higher_in_high")
  "higher_in_low"
}

# 【函数：write_session】将会话版本与验证配置写入文件，供之后核对计算环境。
write_session <- function(path) {
  sink(path)
  cat("timestamp:", as.character(Sys.time()), "\n")
  cat("script: 07_refined_sample_level_validation.R (v2 method revision)\n")
  cat("input_rds:", input_rds, "\n")
  cat("out_dir:", out_dir, "\n")
  cat("n_boot:", n_boot, "\n")
  cat("min_mean_expr:", min_mean_expr, "\n")
  cat("min_pct_expr:", min_pct_expr, "\n")
  cat("similar_abs:", similar_abs, "\n")
  cat("NOTE: statistical unit is sample (n=4).\n")
  cat("Bootstrap unit: cells within each sample (NOT independent samples).\n")
  cat("evidence_level must NOT be upgraded by within-sample bootstrap.\n")
  # 【运行记录】记录 R 与已加载包的版本，帮助以后解释同一代码为何可能得到不同结果。
  print(sessionInfo())
  sink()
}

message("=== 07 refined sample-level validation v2 ===")
message("Start: ", Sys.time())
write_session(file.path(out_dir, "sessionInfo.txt"))
stop_if_missing(input_rds)

message("Loading frozen refined object")
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  try(obj <- JoinLayers(obj, assay = "RNA"), silent = TRUE)
}

if (!"analysis_celltype" %in% colnames(obj@meta.data)) {
  stop("analysis_celltype missing")
}
if (!"sample" %in% colnames(obj@meta.data)) stop("sample missing")

obj$is_refined_nk <- obj$analysis_celltype == "NK_cell"
n_nk <- sum(obj$is_refined_nk)
message("Refined NK: ", n_nk)

stk31_data <- as.numeric(fetch_mat(obj, target_gene, "data")[target_gene, ])
# 【错误分支】尝试运行代码，失败时进入 error 函数；应阅读返回值，区分正常结果与跳过/失败说明。
stk31_counts <- tryCatch(
  as.numeric(fetch_mat(obj, target_gene, "counts")[target_gene, ]),
  error = function(e) stk31_data
)
names(stk31_data) <- colnames(obj)
names(stk31_counts) <- colnames(obj)

is_epi <- obj$analysis_celltype == epi_label
if (sum(is_epi) == 0 && "legacy_manual_celltype" %in% colnames(obj@meta.data)) {
  is_epi <- obj$legacy_manual_celltype == epi_label
}
if (sum(is_epi) == 0) stop("No epithelial cells")

epi_high_primary <- is_epi & stk31_counts > 0
epi_low_primary <- is_epi & !epi_high_primary
samples <- sort(unique(as.character(obj$sample)))
message("Samples: ", paste(samples, collapse = ", "))

# 【分组数量】逐样本统计 high、low 和 NK 数，先判断样本能否满足后续比较的最低数量。
# ---------- group counts ----------
# 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
group_counts <- do.call(rbind, lapply(samples, function(s) {
  idx <- obj$sample == s
  n_high <- sum(idx & epi_high_primary)
  n_low <- sum(idx & epi_low_primary)
  n_nk_s <- sum(idx & obj$is_refined_nk)
  data.frame(
    sample = s,
    total_cells = sum(idx),
    epithelial_cells = sum(idx & is_epi),
    stk31_high_epithelial = n_high,
    stk31_low_epithelial = n_low,
    refined_nk = n_nk_s,
    pct_nk = n_nk_s / sum(idx),
    high_cells_insufficient = n_high < min_high_cells,
    nk_cells_insufficient = n_nk_s < min_nk_cells,
    low_cells_insufficient = n_low < min_low_cells,
    evaluable_for_high_low = (n_high >= min_high_cells) &
      (n_low >= min_low_cells) & (n_nk_s >= min_nk_cells),
    stringsAsFactors = FALSE
  )
}))
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
write.csv(group_counts, file.path(out_dir, "sample_level_group_counts.csv"), row.names = FALSE)

# 【样本内计算】每个样本单独取 high/low 上皮和 NK，避免先混合全部细胞后丢掉样本差异。
# ---------- axis scoring with direction ----------
axis_rows <- list()
effect_rows <- list()

for (axis_name in names(axes)) {
  ax <- axes[[axis_name]]
  lig <- ax$ligands
  rec <- ax$receptors
  sender <- ax$sender
  receiver <- ax$receiver

  for (s in samples) {
    idx <- obj$sample == s
    high_cells <- colnames(obj)[idx & epi_high_primary]
    low_cells <- colnames(obj)[idx & epi_low_primary]
    nk_cells <- colnames(obj)[idx & obj$is_refined_nk]
    all_epi_s <- colnames(obj)[idx & is_epi]

    genes_present_lig <- available_genes(obj, lig)
    genes_present_rec <- available_genes(obj, rec)

    n_high <- length(high_cells)
    n_low <- length(low_cells)
    n_nk_s <- length(nk_cells)
    cell_ok <- group_counts$evaluable_for_high_low[group_counts$sample == s]

    # 【按方向取表达】上皮→NK 与 NK→上皮使用不同端的配体/受体，避免把发送接收颠倒。
    # Compute sender/receiver modules by direction
    if (sender == "epithelial" && receiver == "NK") {
      # 【上皮发送】分别计算 high/low 上皮配体分数，并与同一样本的 NK 受体分数组合。
      # compare high vs low epithelial ligand; * NK receptor
      high_lig <- module_score_cells(obj, lig, high_cells)
      low_lig <- module_score_cells(obj, lig, low_cells)
      nk_rec <- module_score_cells(obj, rec, nk_cells)

      genes_det_sender <- unique(c(high_lig$genes_detected, low_lig$genes_detected))
      genes_det_receiver <- nk_rec$genes_detected

      # detection for direction classification: need ligand in high OR low epi,
      # and receptor in NK
      lig_det_high <- is_detected_module(high_lig)
      lig_det_low <- is_detected_module(low_lig)
      rec_det <- is_detected_module(nk_rec)
      ligand_detected <- lig_det_high || lig_det_low
      receptor_detected <- rec_det

      sender_lig_mean_high <- high_lig$mean_expr
      sender_lig_mean_low <- low_lig$mean_expr
      sender_lig_pct_high <- high_lig$pct_pos
      sender_lig_pct_low <- low_lig$pct_pos
      # report high as primary sender_ligand_pct for table
      sender_ligand_pct <- high_lig$pct_pos
      receiver_receptor_mean <- nk_rec$mean_expr
      receiver_receptor_pct <- nk_rec$pct_pos

      score_high <- as.numeric(high_lig$mean_expr * nk_rec$mean_expr)
      score_low <- as.numeric(low_lig$mean_expr * nk_rec$mean_expr)
      delta <- score_high - score_low

      not_detected_reason <- NA_character_
      if (!cell_ok) {
        direction <- "insufficient"
        not_detected_reason <- "insufficient_cell_counts"
      } else if (length(genes_present_lig) == 0 || length(genes_present_rec) == 0) {
        direction <- "not_detected"
        not_detected_reason <- "genes_absent_from_matrix"
      } else if (!ligand_detected && !receptor_detected) {
        direction <- "not_detected"
        not_detected_reason <- "ligand_and_receptor_not_detected_in_target_cells"
      } else if (!ligand_detected) {
        direction <- "not_detected"
        not_detected_reason <- "ligand_not_detected_in_sender_epithelial"
      } else if (!receptor_detected) {
        direction <- "not_detected"
        not_detected_reason <- "receptor_not_detected_in_receiver_NK"
      } else {
        direction <- label_direction(delta)
        not_detected_reason <- ""
      }

    } else if (sender == "NK" && receiver == "epithelial") {
      # 【NK 发送】NK 配体分数共用，比较 high/low 上皮中的受体分数。
      # NK ligand common; compare high vs low epithelial receptor
      nk_lig <- module_score_cells(obj, lig, nk_cells)
      high_rec <- module_score_cells(obj, rec, high_cells)
      low_rec <- module_score_cells(obj, rec, low_cells)

      genes_det_sender <- nk_lig$genes_detected
      genes_det_receiver <- unique(c(high_rec$genes_detected, low_rec$genes_detected))

      ligand_detected <- is_detected_module(nk_lig)
      receptor_detected <- is_detected_module(high_rec) || is_detected_module(low_rec)

      sender_lig_mean_high <- nk_lig$mean_expr  # same NK sender
      sender_lig_mean_low <- nk_lig$mean_expr
      sender_lig_pct_high <- nk_lig$pct_pos
      sender_lig_pct_low <- nk_lig$pct_pos
      sender_ligand_pct <- nk_lig$pct_pos
      receiver_receptor_mean <- high_rec$mean_expr
      receiver_receptor_pct <- high_rec$pct_pos

      score_high <- as.numeric(nk_lig$mean_expr * high_rec$mean_expr)
      score_low <- as.numeric(nk_lig$mean_expr * low_rec$mean_expr)
      delta <- score_high - score_low

      not_detected_reason <- NA_character_
      if (!cell_ok) {
        direction <- "insufficient"
        not_detected_reason <- "insufficient_cell_counts"
      } else if (length(genes_present_lig) == 0 || length(genes_present_rec) == 0) {
        direction <- "not_detected"
        not_detected_reason <- "genes_absent_from_matrix"
      } else if (!ligand_detected && !receptor_detected) {
        direction <- "not_detected"
        not_detected_reason <- "ligand_and_receptor_not_detected_in_target_cells"
      } else if (!ligand_detected) {
        direction <- "not_detected"
        not_detected_reason <- "ligand_not_detected_in_sender_NK"
      } else if (!receptor_detected) {
        direction <- "not_detected"
        not_detected_reason <- "receptor_not_detected_in_receiver_epithelial"
      } else {
        direction <- label_direction(delta)
        not_detected_reason <- ""
      }

    } else {
      stop("Unsupported sender/receiver: ", sender, " / ", receiver)
    }

    axis_rows[[length(axis_rows) + 1]] <- data.frame(
      sample = s,
      axis = axis_name,
      sender = sender,
      receiver = receiver,
      direction_biology = paste0(sender, "_to_", receiver),
      n_high = n_high,
      n_low = n_low,
      n_nk = n_nk_s,
      genes_present_in_matrix_ligand = paste(genes_present_lig, collapse = ";"),
      genes_present_in_matrix_receptor = paste(genes_present_rec, collapse = ";"),
      genes_detected_in_sender = paste(genes_det_sender, collapse = ";"),
      genes_detected_in_receiver = paste(genes_det_receiver, collapse = ";"),
      ligand_detected = ligand_detected,
      receptor_detected = receptor_detected,
      sender_ligand_mean_for_high_context = sender_lig_mean_high,
      sender_ligand_mean_for_low_context = sender_lig_mean_low,
      sender_ligand_pct = sender_ligand_pct,
      sender_ligand_pct_high_context = sender_lig_pct_high,
      sender_ligand_pct_low_context = sender_lig_pct_low,
      receiver_receptor_mean_high_context = receiver_receptor_mean,
      receiver_receptor_pct = receiver_receptor_pct,
      interaction_score_high = score_high,
      interaction_score_low = score_low,
      delta_interaction_high_minus_low = delta,
      direction = direction,
      not_detected_reason = not_detected_reason,
      evaluable_cell_counts = cell_ok,
      min_mean_expr = min_mean_expr,
      min_pct_expr = min_pct_expr,
      stringsAsFactors = FALSE
    )

    effect_rows[[length(effect_rows) + 1]] <- data.frame(
      sample = s,
      axis = axis_name,
      sender = sender,
      receiver = receiver,
      delta_interaction_high_minus_low = delta,
      direction = direction,
      ligand_detected = ligand_detected,
      receptor_detected = receptor_detected,
      not_detected_reason = not_detected_reason,
      evaluable_cell_counts = cell_ok,
      stringsAsFactors = FALSE
    )
  }
}

axis_df <- do.call(rbind, axis_rows)
effect_df <- do.call(rbind, effect_rows)
write.csv(axis_df, file.path(out_dir, "sample_level_candidate_axis_summary.csv"), row.names = FALSE)
write.csv(effect_df, file.path(out_dir, "sample_level_high_low_effect_directions.csv"), row.names = FALSE)

# 【留一法】每次排除一个样本，重新看其余样本的多数方向；检验结论是否被某一例主导。
# ---------- LOSO ----------
# Full-sample majority among directional labels excluding insufficient/not_detected for majority of "comparable" dirs
loso_rows <- list()
for (axis_name in names(axes)) {
  full <- effect_df[effect_df$axis == axis_name, ]
  full_dirs <- full$direction[full$direction %in% c("higher_in_high", "similar", "higher_in_low")]
  if (length(full_dirs) > 0) {
    full_tab <- table(full_dirs)
    full_majority <- names(full_tab)[which.max(full_tab)]
  } else {
    full_majority <- if (all(full$direction == "not_detected")) "not_detected" else "insufficient"
  }

  for (left_out in samples) {
    sub <- effect_df[effect_df$axis == axis_name & effect_df$sample != left_out, ]
    dirs <- sub$direction[sub$direction %in% c("higher_in_high", "similar", "higher_in_low")]
    n_dir <- length(dirs)
    if (n_dir == 0) {
      maj <- if (all(sub$direction == "not_detected")) "not_detected" else "insufficient"
      maj_frac <- NA_real_
      counts_str <- paste(names(table(sub$direction)), as.integer(table(sub$direction)), sep = "=", collapse = ";")
      stable <- FALSE
    } else {
      tab <- table(dirs)
      maj <- names(tab)[which.max(tab)]
      maj_frac <- as.numeric(max(tab) / n_dir)
      counts_str <- paste(names(tab), as.integer(tab), sep = "=", collapse = ";")
      # majority pattern stable if majority fraction >= 2/3 among remaining directional samples
      stable <- n_dir >= 2 && maj_frac >= (2 / 3)
    }
    changed <- !identical(as.character(maj), as.character(full_majority))
    loso_rows[[length(loso_rows) + 1]] <- data.frame(
      axis = axis_name,
      left_out_sample = left_out,
      n_remaining_directional = n_dir,
      majority_direction = maj,
      majority_fraction = maj_frac,
      remaining_direction_counts = counts_str,
      full_data_majority_direction = full_majority,
      majority_direction_changed_vs_full = changed,
      loso_majority_pattern_stable = stable,
      mean_delta_remaining = ifelse(
        n_dir > 0,
        mean(sub$delta_interaction_high_minus_low[sub$direction %in% c("higher_in_high", "similar", "higher_in_low")], na.rm = TRUE),
        NA_real_
      ),
      note = "loso_majority_pattern_stable does NOT imply direction_consistent_across_samples",
      stringsAsFactors = FALSE
    )
  }
}
loso_df <- do.call(rbind, loso_rows)
write.csv(loso_df, file.path(out_dir, "leave_one_sample_out_summary.csv"), row.names = FALSE)

# 【细胞重采样】在每个样本内部有放回抽细胞，观察估计波动；这不会产生新的独立患者。
# ---------- within-sample cell bootstrap (fast: pre-extract cell scores) ----------
message("Running descriptive within-sample cell bootstrap (n=", n_boot, ")")

# 【缓存分数】先计算各基因集合的逐细胞均值，bootstrap 时直接抽分数，减少重复读大矩阵。
# Precompute per-cell module scores for all unique gene sets
all_gene_sets <- unique(unlist(lapply(axes, function(a) list(
  paste(sort(a$ligands), collapse = "|"),
  paste(sort(a$receptors), collapse = "|")
))))
# map gene-set key -> named numeric vector of cell scores
cell_score_cache <- list()
# 【函数：get_cell_scores】从预先缓存的基因集合分数取逐细胞向量，减少 bootstrap 中重复提取矩阵的成本。
get_cell_scores <- function(genes) {
  key <- paste(sort(available_genes(obj, genes)), collapse = "|")
  if (!nzchar(key)) return(setNames(rep(NA_real_, ncol(obj)), colnames(obj)))
  if (!is.null(cell_score_cache[[key]])) return(cell_score_cache[[key]])
  mat <- fetch_mat(obj, available_genes(obj, genes))
  sc <- as.numeric(Matrix::colMeans(mat))
  names(sc) <- colnames(obj)
  cell_score_cache[[key]] <<- sc
  sc
}

# 【函数：detect_from_scores】依据当前抽样后的分数判断是否达到检测条件，供重采样循环使用。
detect_from_scores <- function(scores) {
  if (length(scores) == 0 || all(!is.finite(scores))) return(FALSE)
  mean_s <- mean(scores, na.rm = TRUE)
  pct_s <- mean(scores > 0, na.rm = TRUE)
  is.finite(mean_s) && mean_s > min_mean_expr && (pct_s >= min_pct_expr || mean_s > 0)
}

boot_rows <- list()
for (axis_name in names(axes)) {
  ax <- axes[[axis_name]]
  lig <- ax$ligands
  rec <- ax$receptors
  sender <- ax$sender
  receiver <- ax$receiver
  lig_scores_all <- get_cell_scores(lig)
  rec_scores_all <- get_cell_scores(rec)

  for (s in samples) {
    idx <- obj$sample == s
    high_cells <- colnames(obj)[idx & epi_high_primary]
    low_cells <- colnames(obj)[idx & epi_low_primary]
    nk_cells <- colnames(obj)[idx & obj$is_refined_nk]

    base_row <- effect_df[effect_df$axis == axis_name & effect_df$sample == s, ]
    if (nrow(base_row) == 0 || base_row$direction[1] %in% c("insufficient", "not_detected")) {
      boot_rows[[length(boot_rows) + 1]] <- data.frame(
        sample = s, axis = axis_name, sender = sender, receiver = receiver,
        n_boot = n_boot,
        bootstrap_unit = "cells_within_sample",
        mean_delta = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        prop_higher_in_high = NA_real_, prop_higher_in_low = NA_real_,
        prop_similar = NA_real_, prop_not_detected = NA_real_,
        within_sample_cell_bootstrap_stable = FALSE,
        note = paste0(
          "skipped_bootstrap; base_direction=",
          ifelse(nrow(base_row), base_row$direction[1], "missing"),
          "; DESCRIPTIVE only; NOT sample-level inference"
        ),
        stringsAsFactors = FALSE
      )
      next
    }

    if (length(high_cells) < min_high_cells || length(low_cells) < min_low_cells ||
        length(nk_cells) < min_nk_cells ||
        all(is.na(lig_scores_all)) || all(is.na(rec_scores_all))) {
      boot_rows[[length(boot_rows) + 1]] <- data.frame(
        sample = s, axis = axis_name, sender = sender, receiver = receiver,
        n_boot = n_boot, bootstrap_unit = "cells_within_sample",
        mean_delta = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        prop_higher_in_high = NA_real_, prop_higher_in_low = NA_real_,
        prop_similar = NA_real_, prop_not_detected = NA_real_,
        within_sample_cell_bootstrap_stable = FALSE,
        note = "insufficient_for_bootstrap; DESCRIPTIVE only",
        stringsAsFactors = FALSE
      )
      next
    }

    h_sc_lig <- lig_scores_all[high_cells]
    l_sc_lig <- lig_scores_all[low_cells]
    n_sc_lig <- lig_scores_all[nk_cells]
    h_sc_rec <- rec_scores_all[high_cells]
    l_sc_rec <- rec_scores_all[low_cells]
    n_sc_rec <- rec_scores_all[nk_cells]

    nh <- length(h_sc_lig)
    nl <- length(l_sc_lig)
    nn <- length(n_sc_lig)
    deltas <- numeric(n_boot)
    dirs_b <- character(n_boot)

    for (b in seq_len(n_boot)) {
      if (sender == "epithelial" && receiver == "NK") {
        h_s <- h_sc_lig[sample.int(nh, nh, replace = TRUE)]
        l_s <- l_sc_lig[sample.int(nl, nl, replace = TRUE)]
        n_s <- n_sc_rec[sample.int(nn, nn, replace = TRUE)]
        lig_ok <- detect_from_scores(h_s) || detect_from_scores(l_s)
        rec_ok <- detect_from_scores(n_s)
        if (!lig_ok || !rec_ok) {
          deltas[b] <- NA_real_
          dirs_b[b] <- "not_detected"
        } else {
          deltas[b] <- mean(h_s) * mean(n_s) - mean(l_s) * mean(n_s)
          dirs_b[b] <- label_direction(deltas[b])
        }
      } else {
        n_s <- n_sc_lig[sample.int(nn, nn, replace = TRUE)]
        h_s <- h_sc_rec[sample.int(nh, nh, replace = TRUE)]
        l_s <- l_sc_rec[sample.int(nl, nl, replace = TRUE)]
        lig_ok <- detect_from_scores(n_s)
        rec_ok <- detect_from_scores(h_s) || detect_from_scores(l_s)
        if (!lig_ok || !rec_ok) {
          deltas[b] <- NA_real_
          dirs_b[b] <- "not_detected"
        } else {
          deltas[b] <- mean(n_s) * mean(h_s) - mean(n_s) * mean(l_s)
          dirs_b[b] <- label_direction(deltas[b])
        }
      }
    }

    prop_h <- mean(dirs_b == "higher_in_high")
    prop_l <- mean(dirs_b == "higher_in_low")
    prop_s <- mean(dirs_b == "similar")
    prop_nd <- mean(dirs_b == "not_detected")
    within_stable <- max(prop_h, prop_l, prop_s, prop_nd) >= 0.6
    fin <- deltas[is.finite(deltas)]
    boot_rows[[length(boot_rows) + 1]] <- data.frame(
      sample = s, axis = axis_name, sender = sender, receiver = receiver,
      n_boot = n_boot,
      bootstrap_unit = "cells_within_sample",
      mean_delta = ifelse(length(fin), mean(fin), NA_real_),
      ci_low = ifelse(length(fin), as.numeric(stats::quantile(fin, 0.025)), NA_real_),
      ci_high = ifelse(length(fin), as.numeric(stats::quantile(fin, 0.975)), NA_real_),
      prop_higher_in_high = prop_h,
      prop_higher_in_low = prop_l,
      prop_similar = prop_s,
      prop_not_detected = prop_nd,
      within_sample_cell_bootstrap_stable = within_stable,
      note = "DESCRIPTIVE only; cells resampled within sample; NOT sample-level independent replicates; must NOT upgrade evidence_level",
      stringsAsFactors = FALSE
    )
  }
}
boot_df <- do.call(rbind, boot_rows)
write.csv(
  boot_df,
  file.path(out_dir, "descriptive_within_sample_cell_bootstrap_summary.csv"),
  row.names = FALSE
)

# 【定义敏感性】比较 counts>0、全体上皮 q75、每样本上皮 q75 三种集合；集合完全相同不构成独立重复验证。
# ---------- STK31 definition sensitivity with set identity ----------
# 【三种定义】原始 counts>0、全体上皮 q75、逐样本上皮 q75；后面用集合重叠检验是否真有不同。
# three definitions of high cells
def1_cells <- colnames(obj)[is_epi & stk31_counts > 0]
q_global <- as.numeric(stats::quantile(stk31_data[is_epi], 0.75, na.rm = TRUE))
if (!is.finite(q_global) || q_global == 0) {
  def2_cells <- colnames(obj)[is_epi & stk31_data > 0]
} else {
  def2_cells <- colnames(obj)[is_epi & stk31_data >= q_global]
}

def3_cells <- character(0)
for (s in samples) {
  epi_s <- is_epi & obj$sample == s
  if (sum(epi_s) < 5) {
    def3_cells <- c(def3_cells, colnames(obj)[epi_s & stk31_data > 0])
  } else {
    q_s <- as.numeric(stats::quantile(stk31_data[epi_s], 0.75, na.rm = TRUE))
    if (!is.finite(q_s) || q_s == 0) {
      def3_cells <- c(def3_cells, colnames(obj)[epi_s & stk31_data > 0])
    } else {
      def3_cells <- c(def3_cells, colnames(obj)[epi_s & stk31_data >= q_s])
    }
  }
}
def3_cells <- unique(def3_cells)

# 【函数：jaccard】计算两个细胞集合的交集大小/并集大小，衡量分组定义的重合程度。
# 1 表示集合相同，0 表示没有交集；两个集合都为空时返回 NA。
jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

identical_12 <- setequal(def1_cells, def2_cells)
identical_13 <- setequal(def1_cells, def3_cells)
identical_23 <- setequal(def2_cells, def3_cells)
all_identical <- identical_12 && identical_13 && identical_23

sens_global <- data.frame(
  definition = c("counts_gt_0", "epithelial_global_q75", "epithelial_per_sample_q75"),
  n_high_cells = c(length(def1_cells), length(def2_cells), length(def3_cells)),
  global_q75_cutoff = c(NA, q_global, NA),
  stringsAsFactors = FALSE
)

sens_sets <- data.frame(
  comparison = c("counts_vs_global_q75", "counts_vs_sample_q75", "global_q75_vs_sample_q75"),
  jaccard = c(jaccard(def1_cells, def2_cells), jaccard(def1_cells, def3_cells), jaccard(def2_cells, def3_cells)),
  identical_cell_set = c(identical_12, identical_13, identical_23),
  stringsAsFactors = FALSE
)

sens_sample <- do.call(rbind, lapply(samples, function(s) {
  epi_s <- colnames(obj)[is_epi & obj$sample == s]
  data.frame(
    sample = s,
    n_epithelial = length(epi_s),
    n_high_counts_gt0 = sum(epi_s %in% def1_cells),
    n_high_global_q75 = sum(epi_s %in% def2_cells),
    n_high_sample_q75 = sum(epi_s %in% def3_cells),
    stringsAsFactors = FALSE
  )
}))

sens_conclusion <- if (all_identical) {
  "Sensitivity analysis is non-informative because all tested rules collapse to STK31 detected versus undetected under the current zero-inflated expression distribution."
} else {
  "STK31 high definitions are not fully identical; interpret axis sensitivity with caution."
}

# write combined sensitivity file: sample counts + set metrics + conclusion rows
sens_out <- rbind(
  data.frame(
    record_type = "sample_counts", sample = sens_sample$sample,
    n_epithelial = sens_sample$n_epithelial,
    n_high_counts_gt0 = sens_sample$n_high_counts_gt0,
    n_high_global_q75 = sens_sample$n_high_global_q75,
    n_high_sample_q75 = sens_sample$n_high_sample_q75,
    jaccard = NA_real_, identical_cell_set = NA,
    conclusion = "",
    stringsAsFactors = FALSE
  ),
  data.frame(
    record_type = "set_comparison", sample = sens_sets$comparison,
    n_epithelial = NA_integer_,
    n_high_counts_gt0 = NA_integer_,
    n_high_global_q75 = NA_integer_,
    n_high_sample_q75 = NA_integer_,
    jaccard = sens_sets$jaccard,
    identical_cell_set = sens_sets$identical_cell_set,
    conclusion = "",
    stringsAsFactors = FALSE
  ),
  data.frame(
    record_type = "conclusion", sample = "all",
    n_epithelial = sum(is_epi),
    n_high_counts_gt0 = length(def1_cells),
    n_high_global_q75 = length(def2_cells),
    n_high_sample_q75 = length(def3_cells),
    jaccard = if (all_identical) 1 else mean(sens_sets$jaccard),
    identical_cell_set = all_identical,
    conclusion = sens_conclusion,
    stringsAsFactors = FALSE
  )
)
write.csv(sens_out, file.path(out_dir, "stk31_definition_sample_level_sensitivity.csv"), row.names = FALSE)

# 【证据汇总】把检出、可评价样本数、方向一致性和留一稳定性放在一起解读。
# ---------- validation summary per axis ----------
summary_rows <- list()
for (axis_name in names(axes)) {
  ax <- axes[[axis_name]]
  sub <- effect_df[effect_df$axis == axis_name, ]
  n_high <- sum(sub$direction == "higher_in_high")
  n_sim <- sum(sub$direction == "similar")
  n_low <- sum(sub$direction == "higher_in_low")
  n_ins <- sum(sub$direction == "insufficient")
  n_nd <- sum(sub$direction == "not_detected")
  n_eval_dir <- n_high + n_sim + n_low

  dirs <- sub$direction[sub$direction %in% c("higher_in_high", "similar", "higher_in_low")]
  direction_consistent <- length(dirs) > 0 && length(unique(dirs)) == 1

  loso_sub <- loso_df[loso_df$axis == axis_name, ]
  # axis-level LOSO majority pattern stable if ALL folds stable AND majority not changing vs full in most folds
  if (nrow(loso_sub) == 0) {
    loso_stable <- FALSE
  } else {
    loso_stable <- all(loso_sub$loso_majority_pattern_stable) &&
      mean(loso_sub$majority_direction_changed_vs_full) <= 0.25
  }

  boot_sub <- boot_df[boot_df$axis == axis_name, ]
  # only count samples where bootstrap was actually run (not skipped)
  ran <- !grepl("^skipped_bootstrap|^insufficient", boot_sub$note)
  if (!any(ran)) {
    boot_cons <- "not_evaluable"
  } else if (all(boot_sub$within_sample_cell_bootstrap_stable[ran])) {
    boot_cons <- "all_evaluable_samples_stable"
  } else if (any(boot_sub$within_sample_cell_bootstrap_stable[ran])) {
    boot_cons <- "some_samples_stable"
  } else {
    boot_cons <- "no_samples_stable"
  }

  # 【证据等级】bootstrap 的稳定并不能弥补只有少量独立样本的限制，不能据此升级机制证据。
  # evidence_level: MUST NOT upgrade based on bootstrap
  if (n_nd == nrow(sub) || (n_eval_dir == 0 && n_nd > 0 && n_ins == 0)) {
    evidence <- "not_detected"
    interp <- "Sender ligand and/or receiver receptor not detected under min_mean_expr/min_pct_expr in target cells."
  } else if (n_eval_dir == 0) {
    evidence <- "insufficient"
    interp <- "Insufficient cells or no directional labels after detection filters."
  } else if (direction_consistent && unique(dirs) == "similar") {
    evidence <- "exploratory_consistent"
    interp <- "Evaluable samples agree on similar (exploratory; n_sample=4). Within-sample bootstrap does not upgrade this."
  } else if (direction_consistent) {
    evidence <- "exploratory_consistent"
    interp <- paste0(
      "Evaluable samples agree on ", unique(dirs),
      " (exploratory; n_sample=4). Within-sample bootstrap does not upgrade this."
    )
  } else {
    evidence <- "exploratory_mixed"
    interp <- paste0(
      "Directions are mixed across the four samples. ",
      "The majority category may remain stable in LOSO summaries, ",
      "but this does not establish cross-sample directional consistency. ",
      "Within-sample cell bootstrap is descriptive only."
    )
  }

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    axis = axis_name,
    sender = ax$sender,
    receiver = ax$receiver,
    direction_biology = paste0(ax$sender, "_to_", ax$receiver),
    samples_evaluable_directional = n_eval_dir,
    samples_higher_in_high = n_high,
    samples_similar = n_sim,
    samples_higher_in_low = n_low,
    samples_insufficient = n_ins,
    samples_not_detected = n_nd,
    direction_consistent_across_samples = direction_consistent,
    loso_majority_pattern_stable = loso_stable,
    within_sample_bootstrap_consistency = boot_cons,
    evidence_level = evidence,
    interpretation = interp,
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "sample_level_validation_summary.csv"), row.names = FALSE)

# ---------- plots ----------
pdf(file.path(out_dir, "sample_level_candidate_axis_plot.pdf"), width = 12, height = 7)
plot_df <- axis_df
plot_df$delta_plot <- plot_df$delta_interaction_high_minus_low
plot_df$delta_plot[plot_df$direction %in% c("not_detected", "insufficient")] <- NA
# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
p1 <- ggplot(plot_df, aes(x = sample, y = delta_plot, fill = direction)) +
  geom_col(na.rm = TRUE) +
  facet_wrap(~axis, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Sample-level high-vs-low delta (direction-aware axes, v2)",
    subtitle = "NA bars omitted for not_detected/insufficient; unit=sample; exploratory",
    y = "delta (high context - low context)", x = NULL
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p1)
dev.off()

pdf(file.path(out_dir, "leave_one_sample_out_plot.pdf"), width = 11, height = 6)
p2 <- ggplot(
  loso_df,
  aes(x = left_out_sample, y = mean_delta_remaining, fill = majority_direction)
) +
  geom_col(na.rm = TRUE) +
  facet_wrap(~axis, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "LOSO mean delta among remaining samples (v2)",
    subtitle = "loso_majority_pattern_stable ≠ direction_consistent_across_samples",
    x = "Left-out sample", y = "Mean delta (remaining)"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p2)
dev.off()

# 【解读边界】将样本量、表达检测及观察性推断等限制随结果保存，便于写报告时准确引用。
# ---------- limitations ----------
lim <- c(
  "# Sample-level validation limitations (v2 method revision)",
  "",
  paste0("Date: ", Sys.time()),
  paste0("Output: ", out_dir),
  "",
  "## Method revisions vs v1",
  "- Axes have explicit sender/receiver (epi->NK vs NK->epi).",
  "- Detection requires expression in target cells (min_mean_expr / min_pct_expr), not only matrix presence.",
  "- Bootstrap fields: within_sample_cell_bootstrap_stable / within_sample_bootstrap_consistency.",
  "- LOSO majority pattern stability is separate from direction_consistent_across_samples.",
  "- STK31 definition sensitivity reports set identity; non-informative if all rules collapse.",
  "",
  "## Design limits",
  "- Only 4 samples. All results exploratory.",
  "- Statistical unit = sample, not cell.",
  "- Interaction score is a proxy, not CellChat probability, not causal.",
  "",
  "## Bootstrap",
  paste0("- n_boot = ", n_boot),
  "- Unit: cells within each sample.",
  "- Does NOT create independent sample replicates.",
  "- Must NOT upgrade evidence_level.",
  "",
  "## STK31 sensitivity",
  sens_conclusion,
  paste0("- identical_cell_set across three definitions: ", all_identical),
  paste0("- n_high counts_gt0 / global_q75 / sample_q75: ",
         length(def1_cells), " / ", length(def2_cells), " / ", length(def3_cells)),
  "",
  "## Conclusion boundary",
  "v2 sample-level results remain exploratory and do not upgrade PI main conclusions:",
  "CellChat refined: MHC-I/TIGIT/TGFb/NKG2D high-vs-low similar candidates; IFNG not_detected;",
  "no STK31-high-specific enhancement; no causality.",
  "",
  "Legacy v1 directory results/refined_sample_level_validation/ is superseded_method_issue."
)
writeLines(lim, file.path(out_dir, "analysis_limitations.md"))

inv <- list.files(out_dir, full.names = TRUE)
write.csv(
  data.frame(file = inv, stringsAsFactors = FALSE),
  file.path(out_dir, "output_file_inventory.csv"),
  row.names = FALSE
)

message("=== v2 summary ===")
print(summary_df[, c(
  "axis", "sender", "receiver", "samples_higher_in_high", "samples_similar",
  "samples_higher_in_low", "samples_not_detected",
  "direction_consistent_across_samples", "loso_majority_pattern_stable",
  "within_sample_bootstrap_consistency", "evidence_level"
)])
message("STK31 sensitivity: ", sens_conclusion)
message("Outputs: ", out_dir)
message("End: ", Sys.time())
