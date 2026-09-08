#!/usr/bin/env Rscript
# ========================================================================
# 【中文阅读指南】公共 QC4 和四个本地样本的 CellChat 通路分析
# 输入：公共全细胞 10x 矩阵/元数据/流式提取程序，以及脚本 11 的本地 CNV 整合对象。
# 流程：统一标签 → 每群限量抽样 → 标准化 → 各分析单位运行 CellChat → 汇总通路及 WNT。
# 公共 QC4 的四个样本合为一个分析单位；本地 tissue1/2/4/5 分别分析，总计五套通讯结果。
# 输出：results/qc4_local_sample_cellchat_pathways，各数据集子目录和 00_combined 汇总目录。
# CELLCHAT_MAX_CELLS_PER_GROUP 默认 500；少于 min_cells=10 的群不进入通讯估计。
# WNT 未检出会记录状态和占位图；含义是在当前数据/阈值下未检出，不能直接认定通路不存在。
# 通讯概率依赖表达、分群和算法配置；跨数据集比较用于探索，不能只凭概率大小判断真实信号强弱。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================

# CellChat pathway analysis for the public paired-QC4 cohort and four local
# gallbladder samples. The public cohort is treated as one analysis unit;
# local samples are analysed separately. Cell labels and per-group sampling
# are harmonized before inference. WNT is exported as a dedicated result set.

# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(CellChat)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260827)

project_dir <- "/home/zhuweiyu/codex-r"
public_dir <- file.path(project_dir, "GBC_STK31_public_scRNA")
local_rds <- file.path(
  project_dir, "results/merged_cnv_binary_annotation/merged_cnv_binary_annotated.rds"
)
out_dir <- file.path(project_dir, "results/qc4_local_sample_cellchat_pathways")
combined_dir <- file.path(out_dir, "00_combined")
qc4_input_dir <- file.path(out_dir, "00_qc4_downsampled_10x")
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc4_input_dir, recursive = TRUE, showWarnings = FALSE)

# 【分析单位】QC4 四个公共样本在本脚本中合并运行；本地四个样本分别运行，结果不能误读成八个独立网络。
qc4_samples <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")
local_samples <- c("tissue1", "tissue2", "tissue4", "tissue5")
dataset_order <- c("QC4_public", local_samples)
min_cells <- 10L
# 【参数 max_cells_per_group】每种细胞群的抽样上限；限制计算规模，并非要求各样本真实细胞比例相同。
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
max_cells_per_group <- as.integer(Sys.getenv("CELLCHAT_MAX_CELLS_PER_GROUP", "500"))
stopifnot(max_cells_per_group >= min_cells)

group_levels <- c(
  "Malignant epithelial cells", "Normal epithelial cells", "T cells", "NK cells",
  "B cells", "Plasma cells", "Myeloid cells", "Mast cells",
  "Fibroblast/stromal cells", "Endothelial cells"
)

group_colors <- c(
  "Malignant epithelial cells" = "#B2182B",
  "Normal epithelial cells" = "#EF8A62",
  "T cells" = "#2166AC",
  "NK cells" = "#053061",
  "B cells" = "#67A9CF",
  "Plasma cells" = "#92C5DE",
  "Myeloid cells" = "#4D9221",
  "Mast cells" = "#A6D96A",
  "Fibroblast/stromal cells" = "#9970AB",
  "Endothelial cells" = "#F1B6DA"
)

# 【函数：stop_if_not】把关键数据约束写成检查：只有 ok 明确为 TRUE 才继续，否则报告 message 并停止。
stop_if_not <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

# 【函数：write_note】将说明文字写入文件，用于记录未检出、跳过或结果阅读说明。
write_note <- function(path, text) {
  writeLines(text, path, useBytes = TRUE)
}

# 【函数：write_placeholder_pdf】无法绘制实际网络时生成带原因的占位 PDF。
# 占位图可让输出完整，但其中的说明不代表分析已检出信号。
write_placeholder_pdf <- function(path, title, subtitle = NULL) {
  # 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = title, size = 5) +
    annotate("text", x = 0, y = -0.15, label = subtitle %||% "", size = 3.5) +
    xlim(-1, 1) + ylim(-0.5, 0.5) + theme_void()
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(path, p, width = 7, height = 4, bg = "white")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# 【函数：safe_pdf】封装 PDF 设备的打开和关闭，on.exit 确保退出函数时关闭图形设备。
safe_pdf <- function(path, expr, width = 8, height = 8) {
  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# 【函数：save_both】把同一张 ggplot 图导出为 PNG 与 PDF；stem 决定基础文件名。
# 尺寸参数影响字体和图形布局，同名运行会更新对应输出。
save_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 320, bg = "white")
}

# 【函数：harmonize_public】把公共数据的 celltype/subtype 映射到统一大类，便于与本地结果比较。
# 无法映射的标签保留 NA，并由后续筛选排除。
harmonize_public <- function(celltype, subtype) {
  output <- rep(NA_character_, length(celltype))
  output[subtype == "Malignant epithelial cells"] <- "Malignant epithelial cells"
  output[subtype == "Normal epithelial cells"] <- "Normal epithelial cells"
  output[celltype %in% c("CD8+ T cell", "CD4+ T cells")] <- "T cells"
  output[celltype == "NK cells"] <- "NK cells"
  output[celltype == "B cells"] <- "B cells"
  output[celltype == "Plasma cells"] <- "Plasma cells"
  output[celltype %in% c(
    "Monocytes & Macrophages", "Neutrophils", "Dendritic cells"
  )] <- "Myeloid cells"
  output[celltype == "Mast cells"] <- "Mast cells"
  output[celltype == "Mesenchymal cells"] <- "Fibroblast/stromal cells"
  output[celltype == "Endothelial cells"] <- "Endothelial cells"
  output
}

# 【函数：harmonize_local】把本地精细标签映射成与公共数据一致的名称，例如 NK_cell→NK cells。
# 这一步改命名和汇总口径，不重新进行细胞身份鉴定。
harmonize_local <- function(label) {
  mapping <- c(
    "Malignant epithelial cells" = "Malignant epithelial cells",
    "Normal epithelial cells" = "Normal epithelial cells",
    "T_cell" = "T cells",
    "NK_cell" = "NK cells",
    "B_cell" = "B cells",
    "Plasma_cell" = "Plasma cells",
    "Macrophage" = "Myeloid cells",
    "Myeloid" = "Myeloid cells",
    "Mast_cell" = "Mast cells",
    "Fibroblast" = "Fibroblast/stromal cells",
    "Fibroblast_Stromal" = "Fibroblast/stromal cells",
    "Endothelial" = "Endothelial cells"
  )
  unname(mapping[as.character(label)])
}

# 【函数：sample_up_to】候选不超过 n 时全取，超过时随机取 n 个；默认无放回。
# 返回原索引的子集，配合固定随机种子提高重复抽样的一致性。
sample_up_to <- function(index, n) {
  # 【抽样】从候选集合抽取元素；replace=TRUE 是有放回抽样，同一个细胞可能重复出现。
  if (length(index) <= n) index else sample(index, n)
}

# 【函数：balanced_qc4_indices】每种细胞先给各 QC4 患者分配基础名额，再从剩余细胞补足群上限。
# 目标是减少单一患者主导；某患者细胞不足时不保证四者最终等量。
balanced_qc4_indices <- function(meta, max_per_group) {
  selected <- integer(0)
  for (group in group_levels) {
    group_index <- which(meta$harmonized_celltype == group)
    if (length(group_index) == 0L) next
    by_patient <- split(group_index, meta$sample_name[group_index])
    by_patient <- by_patient[qc4_samples[qc4_samples %in% names(by_patient)]]
    base_quota <- max(1L, floor(max_per_group / length(qc4_samples)))
    first_pass <- unlist(
      # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
      lapply(by_patient, sample_up_to, n = base_quota), use.names = FALSE
    )
    remaining_quota <- max_per_group - length(first_pass)
    # 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
    remaining <- setdiff(group_index, first_pass)
    fill <- if (remaining_quota > 0L) sample_up_to(remaining, remaining_quota) else integer(0)
    selected <- c(selected, first_pass, fill)
  }
  sort(unique(selected))
}

# 【函数：prepare_qc4_input】核对公共元数据与条形码顺序，选择 QC4 细胞，再调用外部流式程序提取矩阵列。
# 输出一个较小的 10x 输入目录及抽样记录，避免在 R 中一次装入整个大矩阵。
# 已有非空矩阵会被复用；若改变抽样参数，应核对缓存矩阵是否仍对应本次元数据。
prepare_qc4_input <- function() {
  matrix_path <- file.path(
    public_dir, "02_processed_data/All_single_cells/counts/10X_counts/matrix.mtx.gz"
  )
  features_path <- file.path(
    public_dir, "02_processed_data/All_single_cells/counts/10X_counts/features.tsv.gz"
  )
  barcodes_path <- file.path(
    public_dir, "02_processed_data/All_single_cells/counts/10X_counts/barcodes.tsv.gz"
  )
  metadata_path <- file.path(
    public_dir, "02_processed_data/All_single_cells/GBC_Metadata.txt"
  )
  extractor <- file.path(public_dir, "00_tools/bin/stage3_stream_extract_malignant")
  required <- c(matrix_path, features_path, barcodes_path, metadata_path, extractor)
  stop_if_not(all(file.exists(required)), paste("Missing QC4 input:", paste(required[!file.exists(required)], collapse = ", ")))

  message("Reading QC4 metadata and selecting balanced cells")
  # 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
  meta <- fread(metadata_path, showProgress = FALSE)[name != "type"]
  barcodes <- fread(cmd = paste("gzip -cd", shQuote(barcodes_path)), header = FALSE,
                    col.names = "name", showProgress = FALSE)
  features <- fread(cmd = paste("gzip -cd", shQuote(features_path)), header = FALSE,
                    showProgress = FALSE)
  stop_if_not(nrow(meta) == nrow(barcodes) && identical(meta$name, barcodes$name),
              "Public metadata and barcode order differ")
  meta[, harmonized_celltype := harmonize_public(celltype, subtype)]
  eligible <- meta$sample_name %in% qc4_samples & !is.na(meta$harmonized_celltype)
  selected_relative <- balanced_qc4_indices(meta[eligible], max_cells_per_group)
  selected <- which(eligible)[selected_relative]
  qc_meta <- copy(meta[selected])
  qc_meta[, source_dataset := "QC4_public"]
  stop_if_not(nrow(qc_meta) > 0L, "No QC4 cells selected")
  stop_if_not(all(table(qc_meta$harmonized_celltype) >= min_cells),
              "A selected QC4 cell group is below min_cells")

  cell_counts <- qc_meta[, .N, by = .(sample_name, harmonized_celltype)]
  # 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
  fwrite(cell_counts, file.path(qc4_input_dir, "selected_cell_counts_by_patient_group.csv"))
  fwrite(qc_meta, file.path(qc4_input_dir, "metadata.csv.gz"))

  col_map <- integer(nrow(meta))
  col_map[selected] <- seq_along(selected)
  row_map <- integer(nrow(features))
  row_map[1L] <- 1L
  col_map_path <- file.path(qc4_input_dir, "full_to_qc4_downsampled_col.int32.bin")
  row_map_path <- file.path(qc4_input_dir, "dummy_target_row.int32.bin")
  writeBin(as.integer(col_map), col_map_path, size = 4L, endian = .Platform$endian)
  writeBin(as.integer(row_map), row_map_path, size = 4L, endian = .Platform$endian)

  out_matrix <- file.path(qc4_input_dir, "matrix.mtx.gz")
  dummy_matrix <- file.path(qc4_input_dir, "dummy_target_matrix.mtx.gz")
  library_path <- file.path(qc4_input_dir, "cell_library_sizes.tsv")
  audit_path <- file.path(qc4_input_dir, "stream_extract_audit.txt")
  extract_tmp_dir <- file.path(qc4_input_dir, "stream_extract_tmp")
  dir.create(extract_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(out_matrix) || file.info(out_matrix)$size == 0L) {
    message("Streaming the selected QC4 columns from the 1.1-million-cell matrix")
    args <- c(
      matrix_path, col_map_path, row_map_path,
      as.character(nrow(features)), as.character(nrow(meta)), as.character(nrow(qc_meta)), "1",
      out_matrix, dummy_matrix, library_path, audit_path, extract_tmp_dir
    )
    status <- system2(extractor, args = args)
    stop_if_not(status == 0L, paste("QC4 stream extractor failed with status", status))
  }
  stop_if_not(file.exists(out_matrix) && file.info(out_matrix)$size > 0L,
              "QC4 subset matrix was not created")
  file.copy(features_path, file.path(qc4_input_dir, "features.tsv.gz"), overwrite = TRUE)
  barcode_connection <- gzfile(file.path(qc4_input_dir, "barcodes.tsv.gz"), open = "wt")
  writeLines(qc_meta$name, barcode_connection, useBytes = TRUE)
  close(barcode_connection)
  # 【快速表格】data.table 是高效表格结构；DT[条件,计算,by=分组] 表示先筛行，再按组汇总。
  # := 在表内更新列，.N 是当前组的行数；对逐细胞表而言通常就是细胞数。
  fwrite(data.table(name = qc_meta$name, labels = qc_meta$harmonized_celltype),
         file.path(qc4_input_dir, "cellchat_metadata.csv.gz"))
  rm(meta, barcodes, features, col_map, row_map)
  gc()
  invisible(qc_meta)
}

# 【函数：read_qc4_counts】读取前面生成的 QC4 小型 10x 矩阵与分组信息，供 CellChat 使用。
read_qc4_counts <- function() {
  # 【10x 数据】读取 matrix、features 和 barcodes；通常返回基因×细胞的稀疏原始计数矩阵。
  counts <- Read10X(data.dir = qc4_input_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[["Gene Expression"]] %||% counts[[1L]]
  meta <- fread(file.path(qc4_input_dir, "cellchat_metadata.csv.gz"))
  stop_if_not(identical(colnames(counts), meta$name), "QC4 matrix and metadata order differ")
  list(counts = counts, labels = meta$labels, cells = meta$name)
}

# 【函数：downsample_local】只取指定本地样本，按统一细胞类型限定抽样数量。
# 返回所选细胞及标签，后面用这些名字取 counts 对应列。
downsample_local <- function(object, sample_name) {
  index <- which(object$sample_id == sample_name & !is.na(object$harmonized_celltype))
  selected <- unlist(lapply(group_levels, function(group) {
    group_index <- index[object$harmonized_celltype[index] == group]
    sample_up_to(group_index, max_cells_per_group)
  }), use.names = FALSE)
  selected <- sort(unique(selected))
  list(cells = colnames(object)[selected], labels = object$harmonized_celltype[selected])
}

# 【函数：extract_pathway_summary】把 CellChat 每条通路的网络和通讯表整理为数据集级摘要。
# 总概率等指标用于描述当前推断网络，不是患者级显著性检验。
extract_pathway_summary <- function(cellchat, communication, dataset) {
  pathways <- cellchat@netP$pathways
  if (length(pathways) == 0L) {
    return(data.frame(
      dataset = character(), pathway = character(), total_probability = numeric(),
      active_source_target_edges = integer(), ligand_receptor_rows = integer()
    ))
  }
  do.call(rbind, lapply(seq_along(pathways), function(i) {
    pathway <- pathways[i]
    matrix <- cellchat@netP$prob[, , i, drop = TRUE]
    pathway_column <- if ("pathway_name" %in% names(communication)) {
      communication$pathway_name
    } else if ("pathway" %in% names(communication)) {
      communication$pathway
    } else {
      rep(NA_character_, nrow(communication))
    }
    data.frame(
      dataset = dataset,
      pathway = pathway,
      total_probability = sum(matrix, na.rm = TRUE),
      active_source_target_edges = sum(matrix > 0, na.rm = TRUE),
      ligand_receptor_rows = sum(pathway_column == pathway, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

# 【函数：extract_wnt_network】提取 WNT 通路中发送群→接收群的网络数值，并附上数据集名称。
# 未检出时按函数中的空结果分支处理。
extract_wnt_network <- function(cellchat, dataset) {
  match_index <- which(toupper(cellchat@netP$pathways) == "WNT")
  if (length(match_index) == 0L) {
    return(data.frame(
      dataset = character(), source = character(), target = character(), probability = numeric()
    ))
  }
  matrix <- cellchat@netP$prob[, , match_index[1L], drop = TRUE]
  grid <- expand.grid(
    source = rownames(matrix), target = colnames(matrix), stringsAsFactors = FALSE
  )
  grid$probability <- as.numeric(matrix)
  grid$dataset <- dataset
  grid[, c("dataset", "source", "target", "probability")]
}

# 【函数：plot_cellchat_outputs】绘制当前数据集的总体通讯网络和重点通路图，并保存到 dataset_dir。
plot_cellchat_outputs <- function(cellchat, dataset, dataset_dir, pathway_summary) {
  group_size <- as.numeric(table(cellchat@idents))
  names(group_size) <- names(table(cellchat@idents))
  colors <- unname(group_colors[names(group_size)])

  safe_pdf(file.path(dataset_dir, "cellchat_global_circle_count_weight.pdf"), {
    old_par <- par(mfrow = c(1, 2), xpd = TRUE)
    on.exit(par(old_par), add = TRUE)
    CellChat::netVisual_circle(
      cellchat@net$count, vertex.weight = group_size, color.use = colors,
      weight.scale = TRUE, label.edge = FALSE,
      title.name = paste(dataset, "interaction count")
    )
    CellChat::netVisual_circle(
      cellchat@net$weight, vertex.weight = group_size, color.use = colors,
      weight.scale = TRUE, label.edge = FALSE,
      title.name = paste(dataset, "interaction weight")
    )
  }, width = 14, height = 7)

  safe_pdf(file.path(dataset_dir, "cellchat_global_heatmap_count_weight.pdf"), {
    print(CellChat::netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues"))
    print(CellChat::netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds"))
  }, width = 9, height = 8)

  top_pathways <- head(pathway_summary$pathway[order(-pathway_summary$total_probability)], 10L)
  if (length(top_pathways) > 0L) {
    safe_pdf(file.path(dataset_dir, "cellchat_top10_pathway_networks.pdf"), {
      for (pathway in top_pathways) {
        try(CellChat::netVisual_aggregate(
          cellchat, signaling = pathway, layout = "circle",
          color.use = colors, vertex.weight = group_size,
          weight.scale = TRUE, label.edge = FALSE
        ), silent = TRUE)
      }
    }, width = 9, height = 9)
  }
}

# 【函数：plot_wnt_outputs】专门输出 WNT 网络、配体受体和检测状态；没有 WNT 时写明原因。
plot_wnt_outputs <- function(cellchat, dataset, dataset_dir) {
  wnt_dir <- file.path(dataset_dir, "WNT")
  dir.create(wnt_dir, recursive = TRUE, showWarnings = FALSE)
  detected <- any(toupper(cellchat@netP$pathways) == "WNT")
  status <- data.frame(
    dataset = dataset,
    WNT_detected = detected,
    # 【按名称对齐】match(x,y) 返回 x 各元素在 y 中的位置，没找到返回 NA；用来保证标签和表达对应同一细胞。
    matched_pathway = if (detected) cellchat@netP$pathways[match("WNT", toupper(cellchat@netP$pathways))] else NA_character_,
    stringsAsFactors = FALSE
  )
  write.csv(status, file.path(wnt_dir, "WNT_status.csv"), row.names = FALSE)
  if (!detected) {
    write_placeholder_pdf(
      file.path(wnt_dir, "WNT_not_detected.pdf"),
      paste(dataset, "— WNT not detected"),
      "No WNT signaling pathway passed the CellChat inference filters."
    )
    write.csv(data.frame(note = "WNT not detected"),
              file.path(wnt_dir, "WNT_ligand_receptor_interactions.csv"), row.names = FALSE)
    return(status)
  }

  # 【通讯表】从 CellChat 对象提取 source、target、ligand、receptor、prob 等字段，便于后续筛选与作图。
  wnt_communication <- CellChat::subsetCommunication(cellchat, signaling = "WNT")
  write.csv(wnt_communication,
            file.path(wnt_dir, "WNT_ligand_receptor_interactions.csv"), row.names = FALSE)
  group_size <- as.numeric(table(cellchat@idents))
  names(group_size) <- names(table(cellchat@idents))
  colors <- unname(group_colors[names(group_size)])

  safe_pdf(file.path(wnt_dir, "WNT_circle_network.pdf"), {
    CellChat::netVisual_aggregate(
      cellchat, signaling = "WNT", layout = "circle", color.use = colors,
      vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE
    )
  }, width = 9, height = 9)

  try(safe_pdf(file.path(wnt_dir, "WNT_chord_network.pdf"), {
    CellChat::netVisual_aggregate(
      cellchat, signaling = "WNT", layout = "chord", color.use = colors
    )
  }, width = 10, height = 10), silent = TRUE)

  try(safe_pdf(file.path(wnt_dir, "WNT_ligand_receptor_bubble.pdf"), {
    print(CellChat::netVisual_bubble(
      cellchat, signaling = "WNT", remove.isolate = FALSE,
      title.name = paste(dataset, "WNT ligand-receptor interactions")
    ))
  }, width = 12, height = 8), silent = TRUE)

  try(safe_pdf(file.path(wnt_dir, "WNT_ligand_receptor_contribution.pdf"), {
    print(CellChat::netAnalysis_contribution(cellchat, signaling = "WNT"))
  }, width = 9, height = 6), silent = TRUE)

  status
}

# 【函数：run_cellchat】建立 CellChat 对象，使用人类配体受体库推断细胞群间通讯，并汇总通路和输出。
# 输入细胞标签必须与表达矩阵列一一对应；少量细胞群和表达预处理会影响推断。
# 返回内容依本脚本定义，主流程会继续提取通讯表或保存图；这是计算推断结果。
run_cellchat <- function(counts, cells, labels, dataset) {
  dataset_dir <- file.path(out_dir, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
  labels <- factor(as.character(labels), levels = group_levels)
  labels <- droplevels(labels)
  keep <- !is.na(labels)
  counts <- counts[, keep, drop = FALSE]
  cells <- cells[keep]
  labels <- droplevels(labels[keep])
  stop_if_not(identical(colnames(counts), cells), paste(dataset, "count columns and cells differ"))
  valid_groups <- names(which(table(labels) >= min_cells))
  keep <- labels %in% valid_groups
  counts <- counts[, keep, drop = FALSE]
  cells <- cells[keep]
  labels <- droplevels(labels[keep])
  stop_if_not(nlevels(labels) >= 2L, paste(dataset, "has fewer than two valid cell groups"))

  count_table <- data.frame(
    dataset = dataset,
    cell_group = names(table(labels)),
    cells_used = as.integer(table(labels)),
    stringsAsFactors = FALSE
  )
  write.csv(count_table, file.path(dataset_dir, "cell_counts_used.csv"), row.names = FALSE)

  # 【创建对象】把 counts 与细胞信息装进 Seurat；对象同时保存表达矩阵、元数据及后续降维结果。
  seurat <- CreateSeuratObject(
    counts = counts,
    meta.data = data.frame(labels = labels, row.names = cells)
  )
  # 【标准化】默认 LogNormalize 将每个细胞按总计数缩放再取 log1p，减少测序深度差异的影响。
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
  data_input <- LayerData(seurat[["RNA"]], layer = "data")
  meta <- data.frame(labels = labels, row.names = cells, stringsAsFactors = FALSE)
  # 【通讯输入】表达矩阵列名与 meta 行名必须对应；group.by 指定用哪个标签定义发送和接收细胞群。
  cellchat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  # 【通讯候选筛选】先找群中高表达的基因，再通过配体受体数据库筛选候选相互作用。
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  # 【通讯推断】根据群表达与配体受体库估计通讯分值；raw.use=TRUE 指未投影的表达，不等于原始 counts。
  cellchat <- CellChat::computeCommunProb(
    # 【CellChat 参数】triMean 汇总群表达；raw.use=TRUE 指未做投影的表达数据，不是直接使用原始 UMI counts。
    # population.size=FALSE 表示此处不把群体细胞比例纳入通讯计算的该项权重。
    cellchat, type = "triMean", raw.use = TRUE, population.size = FALSE
  )
  # 【通讯过滤】去掉不满足最少细胞数要求的群的通讯，min.cells 是本次阈值。
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  # 【通路汇总】把配体受体层面的通讯聚合到通路层面，再由 aggregateNet 汇总群间网络。
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  # 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
  saveRDS(cellchat, file.path(dataset_dir, "cellchat_object.rds"), compress = FALSE)

  communication <- CellChat::subsetCommunication(cellchat)
  write.csv(communication, file.path(dataset_dir, "cellchat_all_interactions.csv"), row.names = FALSE)
  pathway_summary <- extract_pathway_summary(cellchat, communication, dataset)
  pathway_summary <- pathway_summary[order(-pathway_summary$total_probability), , drop = FALSE]
  write.csv(pathway_summary, file.path(dataset_dir, "cellchat_pathway_summary.csv"), row.names = FALSE)
  plot_cellchat_outputs(cellchat, dataset, dataset_dir, pathway_summary)
  wnt_status <- plot_wnt_outputs(cellchat, dataset, dataset_dir)
  wnt_network <- extract_wnt_network(cellchat, dataset)

  rm(seurat, data_input, meta, counts)
  gc()
  list(
    cellchat = cellchat,
    pathway_summary = pathway_summary,
    wnt_status = wnt_status,
    wnt_network = wnt_network,
    cell_counts = count_table
  )
}

# 【函数：make_combined_plots】把各数据集摘要放在共同的图表中，对照通路及 WNT 的分布。
# 如逐行 z-score，则颜色表示同一行的相对差异，不是原始通讯强度。
make_combined_plots <- function(results) {
  pathway <- rbindlist(lapply(results, `[[`, "pathway_summary"), fill = TRUE)
  wnt_status <- rbindlist(lapply(results, `[[`, "wnt_status"), fill = TRUE)
  wnt_network <- rbindlist(lapply(results, `[[`, "wnt_network"), fill = TRUE)
  cell_counts <- rbindlist(lapply(results, `[[`, "cell_counts"), fill = TRUE)
  fwrite(pathway, file.path(combined_dir, "all_dataset_pathway_summary.csv"))
  fwrite(wnt_status, file.path(combined_dir, "WNT_detection_summary.csv"))
  fwrite(wnt_network, file.path(combined_dir, "WNT_source_target_probabilities.csv"))
  fwrite(cell_counts, file.path(combined_dir, "all_dataset_cell_counts_used.csv"))

  if (nrow(pathway) > 0L) {
    pathway[, dataset := factor(dataset, levels = dataset_order)]
    top <- pathway[, .(max_probability = max(total_probability)), by = pathway][
      order(-max_probability)
    ][1:min(.N, 30L), pathway]
    plot_data <- pathway[pathway %in% top]
    pathway_order <- pathway[, .(maximum = max(total_probability)), by = pathway][
      order(maximum), pathway
    ]
    plot_data[, pathway := factor(pathway, levels = pathway_order)]
    p_dot <- ggplot(
      plot_data,
      aes(dataset, pathway, size = active_source_target_edges, color = total_probability)
    ) +
      geom_point(alpha = 0.9) +
      scale_color_viridis_c(option = "C", trans = "sqrt") +
      scale_size_continuous(range = c(1.5, 9)) +
      labs(
        x = NULL, y = NULL, color = "Total CellChat\nprobability",
        size = "Active source-target\nedges",
        title = "CellChat signaling pathways: public QC4 and local samples",
        subtitle = "Top 30 pathways by maximum total probability; identical inference settings"
      ) +
      theme_bw(base_size = 10.5) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_both(p_dot, file.path(combined_dir, "pathway_comparison_dotplot"),
              10, max(8, 0.28 * length(top) + 2.5))

    wide <- dcast(pathway, pathway ~ dataset, value.var = "total_probability", fill = 0)
    matrix <- as.matrix(wide[, -1])
    rownames(matrix) <- wide$pathway
    row_sd <- apply(matrix, 1, sd)
    # 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
    z <- t(vapply(seq_len(nrow(matrix)), function(i) {
      if (is.finite(row_sd[i]) && row_sd[i] > 0) {
        as.numeric(scale(matrix[i, ]))
      } else {
        rep(0, ncol(matrix))
      }
    }, numeric(ncol(matrix))))
    colnames(z) <- colnames(matrix)
    rownames(z) <- rownames(matrix)
    # 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
    z <- z[intersect(rev(pathway_order), rownames(z)), , drop = FALSE]
    z <- head(z, 40L)
    heat <- as.data.table(as.table(z))
    setnames(heat, c("pathway", "dataset", "z_score"))
    heat[, dataset := factor(dataset, levels = dataset_order)]
    heat[, pathway := factor(pathway, levels = rev(rownames(z)))]
    p_heat <- ggplot(heat, aes(dataset, pathway, fill = z_score)) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
      labs(
        x = NULL, y = NULL, fill = "Pathway-wise\nz-score",
        title = "Relative pathway activity across datasets",
        subtitle = "Top 40 pathways; z-scored within pathway for pattern comparison"
      ) +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
    save_both(p_heat, file.path(combined_dir, "pathway_comparison_zscore_heatmap"),
              9, max(9, 0.24 * nrow(z) + 2.5))
  }

  if (nrow(wnt_network) > 0L) {
    wnt_network[, dataset := factor(dataset, levels = dataset_order)]
    wnt_network[, source := factor(source, levels = group_levels)]
    wnt_network[, target := factor(target, levels = rev(group_levels))]
    wnt_not_detected <- wnt_status[WNT_detected == FALSE]
    wnt_not_detected[, `:=`(
      dataset = factor(dataset, levels = dataset_order),
      source = factor("B cells", levels = group_levels),
      target = factor("B cells", levels = rev(group_levels)),
      label = "WNT not detected"
    )]
    p_wnt <- ggplot(wnt_network, aes(source, target, fill = probability)) +
      geom_tile(color = "grey92", linewidth = 0.15) +
      geom_text(
        data = wnt_not_detected,
        aes(source, target, label = label),
        inherit.aes = FALSE, color = "grey35", size = 4.5
      ) +
      facet_wrap(~dataset, ncol = 2, drop = FALSE) +
      scale_fill_viridis_c(option = "B", trans = "sqrt") +
      labs(
        x = "WNT sender", y = "WNT receiver", fill = "CellChat\nprobability",
        title = "WNT signaling source-target networks",
        subtitle = "Shared cell-group order and color scale"
      ) +
      theme_bw(base_size = 8.5) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1),
        strip.background = element_rect(fill = "grey95")
      )
    save_both(p_wnt, file.path(combined_dir, "WNT_source_target_heatmap"), 16, 13)
  } else {
    write_placeholder_pdf(
      file.path(combined_dir, "WNT_not_detected_in_any_dataset.pdf"),
      "WNT not detected in any dataset"
    )
  }
}

message("=== QC4 + local sample CellChat pathway analysis ===")
message("Output: ", out_dir)
message("Max cells per cell group: ", max_cells_per_group)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
  options(future.globals.maxSize = 12 * 1024^3)
}

prepare_qc4_input()
qc4 <- read_qc4_counts()
results <- list()
message("Running CellChat: QC4_public")
results[["QC4_public"]] <- run_cellchat(qc4$counts, qc4$cells, qc4$labels, "QC4_public")
rm(qc4)
gc()

message("Reading latest local CNV-binary/refined annotation object")
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
local_object <- readRDS(local_rds)
stop_if_not(all(c("sample_id", "analysis_celltype_cnv_binary") %in% colnames(local_object[[]])),
            "Required local annotation columns are absent")
local_object$harmonized_celltype <- harmonize_local(local_object$analysis_celltype_cnv_binary)
if (inherits(local_object[["RNA"]], "Assay5")) {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  local_object <- JoinLayers(local_object, assay = "RNA")
}
local_counts <- LayerData(local_object[["RNA"]], layer = "counts")

for (sample_name in local_samples) {
  message("Running CellChat: ", sample_name)
  selection <- downsample_local(local_object, sample_name)
  results[[sample_name]] <- run_cellchat(
    local_counts[, selection$cells, drop = FALSE],
    selection$cells, selection$labels, sample_name
  )
}

make_combined_plots(results)

manifest <- rbindlist(lapply(dataset_order, function(dataset) {
  data.table(
    dataset = dataset,
    object = file.path(out_dir, dataset, "cellchat_object.rds"),
    all_interactions = file.path(out_dir, dataset, "cellchat_all_interactions.csv"),
    pathway_summary = file.path(out_dir, dataset, "cellchat_pathway_summary.csv"),
    WNT_status = file.path(out_dir, dataset, "WNT", "WNT_status.csv"),
    WNT_interactions = file.path(out_dir, dataset, "WNT", "WNT_ligand_receptor_interactions.csv")
  )
}))
manifest[, files_complete :=
  file.exists(object) & file.exists(all_interactions) & file.exists(pathway_summary) &
  file.exists(WNT_status) & file.exists(WNT_interactions)]
fwrite(manifest, file.path(out_dir, "run_manifest.csv"))
stop_if_not(all(manifest$files_complete), "One or more dataset result sets are incomplete")

# 【运行记录】记录 R 与已加载包的版本，帮助以后解释同一代码为何可能得到不同结果。
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
write_note(
  file.path(out_dir, "README_results.txt"),
  c(
    "Analysis units: public QC4 combined, tissue1, tissue2, tissue4, tissue5.",
    paste0("Cell groups were harmonized and downsampled to <=", max_cells_per_group, " cells/group."),
    "CellChat settings: human database, triMean, population.size=FALSE, min.cells=10.",
    "The 00_combined directory contains pathway comparison figures and the WNT source-target comparison.",
    "Each dataset directory contains global networks, top pathway networks, tables, and a dedicated WNT directory.",
    "A missing WNT pathway is reported explicitly and is not interpreted as biological absence."
  )
)
message("QC4_LOCAL_SAMPLE_CELLCHAT_PATHWAYS_PASS")
