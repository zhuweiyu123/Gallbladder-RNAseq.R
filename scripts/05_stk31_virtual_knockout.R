# ========================================================================
# 【中文阅读指南】STK31 虚拟敲除与候选方向探索
# 输入：四样本合并对象和 cluster_celltype_annotation.csv；只取注释为 Epithelial 的细胞。
# 流程：筛选网络基因 → scTenifoldKnk 虚拟移除 STK31 出边 → 匹配阳性/阴性细胞估计方向 → 排序。
# 输出：默认 results/merged_stk31_virtual_knockout；包含网络对象、扰动表、方向预测和稳定性图表。
# STK31-positive 在这里指原始 counts > 0；STK31_KO_* 环境变量可控制规模、重复次数和输出目录。
# 网络 distance 表示扰动幅度，不自带上调/下调方向；方向来自另一个观察性匹配模型。
# 预测变化取 -(阳性组 - 匹配阴性组)，是用于提出实验假设的近似，并非真实敲除后的测量值。
# 若 STK31 在推断网络中没有出边，脚本会标记网络结果不可解释；不能把零扰动当作无生物学作用。
# 阅读顺序：文件开头的路径/参数 → 工具函数 → 主流程；函数定义本身不会执行分析。
# R 入门：<- 是赋值；$ 取一列/一个成员；[行,列] 取子集；c() 建向量；list() 装不同类型对象。
# NA 表示缺失，不等于 0；counts 是原始计数，data 通常是 log 标准化表达。
# 运行环境：本项目默认在 Linux 服务器 /home/zhuweiyu/codex-r 下用 Rscript 运行。
# 本次中文注释用于解释现有实现；原有计算语句、参数、输出名称保持不变。
# ========================================================================
# 【加载依赖】library 加载本脚本用到的包；外层只隐藏启动提示，不会安装缺失的包。
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# 【可重复性】固定随机数起点，使同一环境下的抽样/随机算法更易复现；不同包版本仍可能产生差异。
set.seed(20260714)

# 固定远程服务器项目路径；本脚本默认在 /home/zhuweiyu/codex-r 上运行。
project_dir <- "/home/zhuweiyu/codex-r"
input_rds <- file.path(
  project_dir,
  "results/merged_basic_seurat/gallbladder_cancer_merged_basic_seurat.rds"
)
annotation_file <- file.path(
  project_dir,
  "results/merged_stk31_nk_analysis/cluster_celltype_annotation.csv"
)
# 【可调参数】Sys.getenv 先读环境变量，未设置时使用代码中的默认值；as.integer/as.numeric 把文本转为数值。
out_dir <- Sys.getenv(
  "STK31_KO_OUT_DIR",
  file.path(project_dir, "results/merged_stk31_virtual_knockout")
)
# 【输出目录】recursive=TRUE 可连同父目录一起建立；路径变量决定结果实际写到哪里。
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
epithelial_label <- "Epithelial"
# 【参数 max_genes】虚拟敲除网络允许的基因规模目标；规模越大通常越耗内存和时间。
max_genes <- as.integer(Sys.getenv("STK31_KO_MAX_GENES", "1500"))
# 【参数 n_networks】构建的网络重复数量，用于网络估计；增加会延长运行时间。
n_networks <- as.integer(Sys.getenv("STK31_KO_NETWORKS", "10"))
# 【参数 n_network_cells】每次网络构建所抽取的细胞数量，上限还会受可用细胞数限制。
n_network_cells <- as.integer(Sys.getenv("STK31_KO_NETWORK_CELLS", "1000"))
# 【参数 n_bootstrap】方向模型的 bootstrap 重复次数；增加可改善重采样精度，但不增加独立样本数。
n_bootstrap <- as.integer(Sys.getenv("STK31_KO_BOOTSTRAP", "100"))
# 【参数 n_cores】并行核心数，应与服务器资源配置相匹配。
n_cores <- as.integer(Sys.getenv("STK31_KO_CORES", "4"))
# 【参数 min_gene_pct】网络候选基因的最低检测细胞比例；0.01 即 1%，另结合最低细胞数规则。
min_gene_pct <- as.numeric(Sys.getenv("STK31_KO_MIN_GENE_PCT", "0.01"))
# 【参数 reuse_network】是否复用已有网络对象；复用前要确认旧对象来自相同输入和网络配置。
reuse_network <- tolower(Sys.getenv("STK31_KO_REUSE_NETWORK", "false")) == "true"

required_packages <- c("scTenifoldKnk")
missing_packages <- required_packages[
  # 【固定返回类型】vapply 与 lapply 类似，但需要指定每次返回的类型和长度，便于尽早发现不一致。
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with install.packages('scTenifoldKnk')."
  )
}

# 本分析围绕 STK31-HLA-I-NK 机制假设；这些基因会被优先保留用于方向汇总。
mechanism_sets <- list(
  HLA_I = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M"),
  Antigen_processing = c(
    "TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9", "PSMB10",
    "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3"
  ),
  IFN_response = c(
    "IFNGR1", "IFNGR2", "STAT1", "STAT2", "IRF1", "IRF7", "JAK1",
    "JAK2", "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"
  ),
  TIGIT_NECTIN = c("NECTIN2", "NECTIN3", "PVR", "TIGIT", "CD96"),
  NKG2D_ligands = c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "KLRK1"),
  TGF_beta = c(
    "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
    "ACVR1", "ACVR1B", "GDF15"
  ),
  NK_interaction = c(
    "TNFSF10", "TNFRSF10A", "TNFRSF10B", "FAS", "FASLG", "LGALS9",
    "HAVCR2", "ICAM1", "ITGAL"
  )
)

# 兼容 Seurat v4/v5：v5 使用 layer，旧版本使用 slot。
# 【函数：get_assay_matrix】按 SeuratObject 版本选择 layer 或 slot 参数，从 RNA assay 取指定矩阵。
# 网络输入使用 counts；画图或评分可使用标准化表达，具体由调用处决定。
get_assay_matrix <- function(object, layer) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = "RNA", layer = layer)
  } else {
    GetAssayData(object, assay = "RNA", slot = layer)
  }
}

# 【函数：write_plot】把传入的绘图对象写到指定文件；width、height 控制版面大小。
# 将画图与保存封装起来，可以让同一套输出保持一致尺寸。
write_plot <- function(filename, plot, width = 8, height = 6) {
  # 【保存图形】输出格式由扩展名决定；width/height 默认按英寸，dpi 主要影响位图清晰度。
  ggsave(
    filename = file.path(out_dir, filename), plot = plot,
    width = width, height = height, units = "in"
  )
}

# 【函数：scale_numeric】把数值转为标准化分数，便于不同量纲的指标组合。
# 没有变化或有效值不足时返回零，避免标准差为零引起异常。
scale_numeric <- function(x) {
  x <- as.numeric(x)
  if (length(unique(x[is.finite(x)])) < 2) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# 将原始 counts 转为 log-normalized 表达矩阵；这里不做表达插补。
# 【函数：make_log_normalized】将每个细胞的 counts 除以该细胞总计数，再乘 10000 并做 log1p。
# log1p(x)=log(1+x)，可保留零值；这里没有为未检出的基因补值。
make_log_normalized <- function(counts) {
  library_size <- Matrix::colSums(counts)
  if (any(library_size <= 0)) stop("Zero-library epithelial cells remain.")
  normalized <- counts %*% Matrix::Diagonal(x = 10000 / library_size)
  dimnames(normalized) <- dimnames(counts)
  normalized@x <- log1p(normalized@x)
  normalized
}

# 网络基因 = 检出率足够的高变基因 + 非全零的机制候选基因。
# 【函数：select_network_genes】从足够多细胞检测到的基因中按表达方差挑选，并优先保留可检测的机制基因。
# maximum 控制网络规模；稠密网络的内存随基因数增加而快速增长。
select_network_genes <- function(counts, log_normalized, maximum, required_genes) {
  detected <- Matrix::rowSums(counts > 0)
  minimum_cells <- max(3L, ceiling(ncol(counts) * min_gene_pct))
  eligible <- names(detected)[detected >= minimum_cells]
  # 【集合交集】intersect 只保留两份名单共有的元素，常用于筛选当前数据真正有的基因/细胞。
  required <- intersect(required_genes, names(detected)[detected > 0])
  eligible <- union(eligible, required)

  x <- log_normalized[eligible, , drop = FALSE]
  # 【按基因求均值】在基因×细胞矩阵中，rowMeans 对每一行跨细胞求平均；若输入是 >0 的逻辑矩阵，得到检出比例。
  gene_mean <- Matrix::rowMeans(x)
  gene_var <- Matrix::rowMeans(x ^ 2) - gene_mean ^ 2
  gene_var[!is.finite(gene_var)] <- 0
  ranked <- names(sort(gene_var, decreasing = TRUE))

  variable_slots <- max(0L, maximum - length(required))
  # 【集合差集】setdiff(a,b) 返回 a 中不属于 b 的元素，用于找缺失基因或定义比较的另一组。
  selected <- unique(c(required, head(setdiff(ranked, required), variable_slots)))
  selected[selected %in% rownames(counts)]
}

# 在同一样本内，为 STK31 阳性细胞匹配 QC 水平最接近的 STK31 阴性细胞。
# 【函数：nearest_negative_pool】根据 log1p(UMI数) 和 log1p(基因数) 标准化后的距离，为阳性细胞找最近的阴性候选。
# 返回每个阳性细胞的候选索引；调用处先限制同一样本，k 控制候选数量。
nearest_negative_pool <- function(metadata, positive_indices, negative_indices, k = 10L) {
  covariates <- cbind(
    log1p(metadata$nCount_RNA),
    log1p(metadata$nFeature_RNA)
  )
  covariates <- apply(covariates, 2, function(x) {
    spread <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / spread
  })

  # 【批量处理】lapply 对向量/列表的每个元素执行一次函数，返回列表；rbind/do.call 可再把结果按行拼表。
  lapply(positive_indices, function(i) {
    distance <- rowSums((covariates[negative_indices, , drop = FALSE] -
      matrix(covariates[i, ], nrow = length(negative_indices), ncol = 2, byrow = TRUE)) ^ 2)
    negative_indices[head(order(distance), min(k, length(distance)))]
  })
}

# 方向预测是观察性匹配分析，不是 scTenifoldKnk 网络敲除结果。
# 计算方式：预测敲除变化 = -(STK31 阳性细胞 - 匹配阴性细胞)。
# 只有 bootstrap CI 不跨 0 且 >=75% 样本方向一致时，才给出升高/降低方向。
# 【函数：estimate_direction】各样本内匹配 STK31 阳性/阴性细胞，以 -(阳性均值-阴性均值) 近似预测敲除方向。
# 对细胞有放回抽样形成区间，再要求至少 75% 样本方向一致才标记增加/降低。
# 返回总体预测、逐样本结果和留一结果；这一步独立于网络扰动分析，仍属于观察性探索。
estimate_direction <- function(expression, metadata, positive, samples, bootstraps) {
  if (!is.null(names(positive))) {
    positive <- positive[colnames(expression)]
  }
  positive <- as.logical(positive)
  if (length(positive) != ncol(expression) || anyNA(positive)) {
    stop("STK31 detection status is not aligned to expression columns.")
  }

  sample_names <- sort(unique(as.character(samples)))
  sample_effects <- list()
  bootstrap_effects <- vector("list", bootstraps)
  pools <- list()

  for (sample_name in sample_names) {
    sample_idx <- which(samples == sample_name)
    pos_idx <- sample_idx[positive[sample_idx]]
    neg_idx <- sample_idx[!positive[sample_idx]]
    message(
      "Direction model ", sample_name, ": ", length(pos_idx),
      " STK31-positive / ", length(neg_idx), " negative epithelial cells"
    )
    if (length(pos_idx) < 3 || length(neg_idx) < 3) next

    pools[[sample_name]] <- list(
      positive = pos_idx,
      nearest_negative = nearest_negative_pool(metadata, pos_idx, neg_idx)
    )
    matched_neg <- vapply(pools[[sample_name]]$nearest_negative, `[`, integer(1), 1L)
    sample_effects[[sample_name]] <- -(
      Matrix::rowMeans(expression[, pos_idx, drop = FALSE]) -
        Matrix::rowMeans(expression[, matched_neg, drop = FALSE])
    )
  }

  if (length(sample_effects) < 2) {
    stop("Fewer than two samples contain enough STK31-positive epithelial cells.")
  }

  sample_matrix <- do.call(cbind, sample_effects)
  rownames(sample_matrix) <- rownames(expression)
  colnames(sample_matrix) <- names(sample_effects)

  for (iteration in seq_len(bootstraps)) {
    iteration_effects <- lapply(pools, function(pool) {
      # 【抽样】从候选集合抽取元素；replace=TRUE 是有放回抽样，同一个细胞可能重复出现。
      selected <- sample(seq_along(pool$positive), length(pool$positive), replace = TRUE)
      pos_idx <- pool$positive[selected]
      neg_idx <- vapply(selected, function(j) {
        sample(pool$nearest_negative[[j]], 1L)
      }, integer(1))
      -(Matrix::rowMeans(expression[, pos_idx, drop = FALSE]) -
        Matrix::rowMeans(expression[, neg_idx, drop = FALSE]))
    })
    bootstrap_effects[[iteration]] <- rowMeans(do.call(cbind, iteration_effects))
  }

  bootstrap_matrix <- do.call(cbind, bootstrap_effects)
  rownames(bootstrap_matrix) <- rownames(expression)
  combined <- rowMeans(sample_matrix)
  ci_low <- apply(bootstrap_matrix, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(bootstrap_matrix, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  expected_sign <- sign(combined)
  consistency <- rowMeans(sign(sample_matrix) == expected_sign)

  direction <- ifelse(
    ci_low > 0 & consistency >= 0.75, "predicted_increase",
    ifelse(ci_high < 0 & consistency >= 0.75, "predicted_decrease", "uncertain")
  )

  result <- data.frame(
    gene = rownames(expression),
    predicted_ko_change = combined,
    bootstrap_ci_low = ci_low,
    bootstrap_ci_high = ci_high,
    sample_direction_consistency = consistency,
    samples_evaluated = ncol(sample_matrix),
    predicted_direction = direction,
    stringsAsFactors = FALSE
  )

  stability <- data.frame(
    gene = rep(rownames(sample_matrix), times = ncol(sample_matrix)),
    sample = rep(colnames(sample_matrix), each = nrow(sample_matrix)),
    predicted_ko_change = as.vector(sample_matrix),
    stringsAsFactors = FALSE
  )

  leave_one_out <- do.call(rbind, lapply(seq_len(ncol(sample_matrix)), function(i) {
    retained <- sample_matrix[, -i, drop = FALSE]
    effect <- rowMeans(retained)
    data.frame(
      gene = rownames(sample_matrix),
      held_out_sample = colnames(sample_matrix)[i],
      predicted_ko_change = effect,
      predicted_direction = ifelse(
        effect > 0, "predicted_increase",
        ifelse(effect < 0, "predicted_decrease", "no_change")
      ),
      stringsAsFactors = FALSE
    )
  }))

  list(result = result, stability = stability, leave_one_out = leave_one_out)
}

# 给基因标注属于哪个机制集合，方便后续解释 HLA-I、IFN、NK 互作等方向。
# 【函数：annotate_mechanisms】查询每个基因属于哪些人工机制集合，添加便于阅读的分类标签。
# 一个基因可能属于多个集合；不在集合中不代表没有生物学功能。
annotate_mechanisms <- function(genes) {
  vapply(genes, function(gene) {
    labels <- names(mechanism_sets)[vapply(mechanism_sets, function(x) gene %in% x, logical(1))]
    if (length(labels) == 0) "Other" else paste(labels, collapse = ";")
  }, character(1))
}

# 当 STK31 网络敲除不可解释时，不做通路富集伪显著性，只汇总机制基因的观察性方向。
# 【函数：summarize_curated_directions】按人工机制集合汇总预测增加、降低和不确定的基因数量及名称。
# 汇总依赖前面的方向模型，不是重新进行一轮实验或通路检验。
summarize_curated_directions <- function(direction_table) {
  output <- lapply(names(mechanism_sets), function(pathway) {
    members <- intersect(mechanism_sets[[pathway]], direction_table$gene)
    # 【按名称对齐】match(x,y) 返回 x 各元素在 y 中的位置，没找到返回 NA；用来保证标签和表达对应同一细胞。
    selected <- direction_table[match(members, direction_table$gene), , drop = FALSE]
    data.frame(
      pathway = pathway,
      genes_evaluated = length(members),
      predicted_increase = sum(selected$predicted_direction == "predicted_increase"),
      predicted_decrease = sum(selected$predicted_direction == "predicted_decrease"),
      uncertain = sum(selected$predicted_direction == "uncertain"),
      increase_genes = paste(selected$gene[
        selected$predicted_direction == "predicted_increase"
      ], collapse = ";"),
      decrease_genes = paste(selected$gene[
        selected$predicted_direction == "predicted_decrease"
      ], collapse = ";"),
      mean_predicted_ko_change = if (length(members) == 0) {
        NA_real_
      } else {
        mean(selected$predicted_ko_change, na.rm = TRUE)
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

# 读取合并后的 Seurat 对象，并接入前面人工整理的 cluster-celltype 注释。
message("Reading merged Seurat object: ", input_rds)
# 【读取中间对象】readRDS 恢复之前保存的 R 对象；检查文件路径和对象来自哪一版注释。
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  # 【Seurat v5 数据层】将拆分的数据层合并以供下游读取；这不是批次校正，也不是重新聚类。
  obj <- JoinLayers(obj, assay = "RNA")
}

if (!file.exists(annotation_file)) stop("Missing annotation file: ", annotation_file)
# 【读入表格】把已有 CSV/TSV 读入内存；后续列名检查用于确认文件格式符合预期。
annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)
cluster_to_celltype <- setNames(
  annotation$manual_celltype,
  as.character(annotation$cluster)
)
obj$manual_celltype <- unname(cluster_to_celltype[as.character(obj$seurat_clusters)])
obj$manual_celltype[is.na(obj$manual_celltype)] <- "Unassigned"

epithelial_cells <- colnames(obj)[obj$manual_celltype == epithelial_label]
if (length(epithelial_cells) == 0) stop("No epithelial cells found.")
# 【取细胞子集】按 cells 或条件保留需要的细胞；这一步改变本次分析范围，要同步核对后面的分母。
epithelial <- subset(obj, cells = epithelial_cells)
counts <- get_assay_matrix(epithelial, "counts")
if (!target_gene %in% rownames(counts)) stop(target_gene, " is absent from RNA counts.")

# 这里的 STK31-positive 指 counts > 0，不是高表达分组；低表达数据要避免过度解释。
stk31_counts <- as.numeric(counts[target_gene, ])
stk31_positive <- stk31_counts > 0
names(stk31_positive) <- colnames(counts)
if (sum(stk31_positive) < 20) {
  stop("Too few STK31-positive epithelial cells for exploratory knockout analysis.")
}

log_normalized <- make_log_normalized(counts)
required_genes <- unique(c(target_gene, unlist(mechanism_sets, use.names = FALSE)))
network_genes <- select_network_genes(
  counts, log_normalized, max_genes, required_genes
)
if (!target_gene %in% network_genes) stop("STK31 was removed during gene selection.")

counts_network <- as.matrix(counts[network_genes, , drop = FALSE])
storage.mode(counts_network) <- "numeric"
n_network_cells <- min(n_network_cells, ncol(counts_network))

# scTenifoldKnk 使用原始 counts 推断网络，并通过清零 STK31 的出边模拟虚拟敲除。
message(
  "Running scTenifoldKnk with ", nrow(counts_network), " genes and ",
  ncol(counts_network), " epithelial cells"
)
ko_file <- file.path(out_dir, "stk31_virtual_knockout_object.rds")
# 重画图或改下游表格时可复用已完成的网络对象，避免重复跑耗时网络构建。
if (reuse_network && file.exists(ko_file)) {
  message("Reusing existing scTenifoldKnk object: ", ko_file)
  ko <- readRDS(ko_file)
} else {
  ko <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = counts_network,
    gKO = target_gene,
    qc = FALSE,
    nc_nNet = n_networks,
    nc_nCells = n_network_cells,
    nc_nComp = 3,
    nc_q = 0.9,
    td_K = 3,
    ma_nDim = 2,
    nCores = n_cores
  )
  # 【保存中间对象】保存完整 R 对象供后续继续分析；RDS 需用 readRDS 读取，不能当 CSV 打开。
  saveRDS(ko, ko_file)
}

network_results <- ko$diffRegulation
network_results <- network_results[network_results$gene != target_gene, , drop = FALSE]
network_results <- network_results[order(network_results$p.adj, -network_results$distance), ]
wt_network <- as.matrix(ko$tensorNetworks$WT)
# 如果 STK31 在 WT 网络中没有出边，网络敲除层不能作为调控证据。
stk31_outdegree <- sum(abs(wt_network[target_gene, ]), na.rm = TRUE)
stk31_connected_genes <- sum(abs(wt_network[target_gene, ]) > 0, na.rm = TRUE)
network_interpretable <- stk31_outdegree > 0
network_results$mechanism_category <- annotate_mechanisms(network_results$gene)
network_results$network_evidence_status <- if (network_interpretable) {
  "interpretable"
} else {
  "not_interpretable_zero_STK31_connectivity"
}
# 【导出表格】将当前统计或注释写入文件，方便用 Excel 查看；文件名与目录见本次调用。
write.csv(
  network_results,
  file.path(out_dir, "stk31_network_perturbation.csv"),
  row.names = FALSE
)
if (!network_interpretable) {
  warning("STK31 has zero outdegree in the inferred network; KO results are not interpretable.")
}

metadata <- epithelial@meta.data
metadata <- metadata[colnames(log_normalized), , drop = FALSE]
# 观察性方向模型单独运行，给湿实验提供候选方向，不等同于因果证明。
direction <- estimate_direction(
  expression = log_normalized[network_genes, , drop = FALSE],
  metadata = metadata,
  positive = stk31_positive,
  samples = as.character(metadata$sample),
  bootstraps = n_bootstrap
)
direction$result <- direction$result[direction$result$gene != target_gene, ]
direction$stability <- direction$stability[direction$stability$gene != target_gene, ]
direction$leave_one_out <- direction$leave_one_out[
  direction$leave_one_out$gene != target_gene,
]
write.csv(
  direction$result,
  file.path(out_dir, "stk31_ko_direction_predictions.csv"),
  row.names = FALSE
)
write.csv(
  direction$stability,
  file.path(out_dir, "stk31_ko_sample_stability.csv"),
  row.names = FALSE
)
write.csv(
  direction$leave_one_out,
  file.path(out_dir, "stk31_ko_leave_one_sample_out.csv"),
  row.names = FALSE
)

priority <- merge(network_results, direction$result, by = "gene", all.x = TRUE)
priority$mechanism_category <- annotate_mechanisms(priority$gene)
priority$direction_score <- scale_numeric(abs(priority$predicted_ko_change))
# 只有网络可解释时才把网络分数纳入优先级；否则只按观察性方向和样本一致性排序。
if (network_interpretable) {
  priority$network_score <- scale_numeric(log1p(priority$distance))
  priority$priority_score <- priority$network_score + priority$direction_score +
    priority$sample_direction_consistency
  network_top_cutoff <- stats::quantile(priority$distance, 0.90, na.rm = TRUE)
  priority$confidence <- ifelse(
    !is.na(priority$p.adj) & priority$p.adj < 0.05 &
      priority$predicted_direction != "uncertain" &
      priority$sample_direction_consistency >= 0.75,
    "high",
    ifelse(
      priority$distance >= network_top_cutoff &
        priority$predicted_direction != "uncertain",
      "moderate", "exploratory"
    )
  )
} else {
  priority$network_score <- NA_real_
  priority$priority_score <- priority$direction_score +
    priority$sample_direction_consistency
  priority$confidence <- ifelse(
    priority$predicted_direction != "uncertain" &
      priority$sample_direction_consistency >= 0.75,
    "observational_supported", "observational_exploratory"
  )
}
priority <- priority[order(-priority$priority_score), ]
write.csv(
  priority,
  file.path(out_dir, "stk31_ko_priority_candidates.csv"),
  row.names = FALSE
)

pathways <- summarize_curated_directions(direction$result)
write.csv(
  pathways,
  file.path(out_dir, "stk31_ko_mechanism_direction_summary.csv"),
  row.names = FALSE
)

sample_summary <- do.call(rbind, lapply(split(seq_len(ncol(epithelial)), metadata$sample), function(idx) {
  data.frame(
    sample = as.character(metadata$sample[idx[1]]),
    epithelial_cells = length(idx),
    stk31_positive_cells = sum(stk31_positive[idx]),
    pct_stk31_positive = 100 * mean(stk31_positive[idx]),
    stringsAsFactors = FALSE
  )
}))
write.csv(sample_summary, file.path(out_dir, "stk31_expression_basis.csv"), row.names = FALSE)

analysis_summary <- data.frame(
  item = c(
    "target_gene", "epithelial_cells", "stk31_positive_cells",
    "pct_stk31_positive", "samples", "network_genes", "network_count",
    "network_cells_per_run", "bootstrap_iterations", "stk31_network_outdegree",
    "stk31_connected_genes", "network_interpretation",
    "direction_interpretation"
  ),
  value = c(
    target_gene, ncol(epithelial), sum(stk31_positive), 100 * mean(stk31_positive),
    length(unique(metadata$sample)), length(network_genes), n_networks,
    n_network_cells, n_bootstrap, stk31_outdegree, stk31_connected_genes,
    if (network_interpretable) "interpretable" else "not_interpretable",
    "Association-based prediction; experimental validation required"
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_summary, file.path(out_dir, "analysis_summary.csv"), row.names = FALSE)

# 图 1：展示每个样本里 STK31-positive 上皮细胞比例，判断数据基础是否足够。
# 【ggplot 图层】aes 把表格列映射到坐标/颜色/大小，后面的 + 逐层加入点、线、主题和标签。
expression_plot <- ggplot(
  sample_summary,
  aes(x = sample, y = pct_stk31_positive, fill = sample)
) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(stk31_positive_cells, "/", epithelial_cells)),
    vjust = -0.4, size = 3.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "STK31 detection in merged epithelial cells",
    x = NULL, y = "STK31-positive epithelial cells (%)"
  ) +
  theme_classic(base_size = 12)
write_plot("01_stk31_expression_basis.pdf", expression_plot, 7, 5)

# 图 2：网络层结果；若 STK31 没有出边，则输出诊断图而不是误导性的扰动排名。
if (network_interpretable) {
  top_network <- head(network_results[order(-network_results$distance), ], 25)
  # 【类别顺序】factor 把文本变成类别；levels 控制图表顺序，也可能影响模型参考组。
  top_network$gene <- factor(top_network$gene, levels = rev(top_network$gene))
  network_plot <- ggplot(
    top_network,
    aes(x = distance, y = gene, color = p.adj < 0.05)
  ) +
    geom_segment(aes(x = 0, xend = distance, yend = gene), color = "grey80") +
    geom_point(size = 2.8) +
    scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
    labs(
      title = "Top genes perturbed by virtual STK31 knockout",
      x = "Network manifold displacement", y = NULL, color = "FDR < 0.05"
    ) +
    theme_classic(base_size = 11)
} else {
  network_plot <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.58,
      label = "STK31 has no outgoing edges in the inferred network",
      size = 5, fontface = "bold"
    ) +
    annotate(
      "text", x = 0.5, y = 0.42,
      label = paste0(
        "Network virtual knockout is not interpretable\n",
        "Use the sample-matched direction results only as exploratory hypotheses"
      ),
      size = 4, color = "#4D4D4D"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = "STK31 network knockout diagnostic") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
}
write_plot("02_stk31_network_perturbation.pdf", network_plot, 8, 7)

# 图 3：观察性方向预测，横轴为预测敲除后的表达变化，线段为 bootstrap 置信区间。
plot_candidates <- priority[
  priority$confidence %in% c("high", "moderate", "observational_supported") |
    priority$mechanism_category != "Other",
]
plot_candidates <- head(plot_candidates[order(-plot_candidates$priority_score), ], 30)
plot_candidates$gene <- factor(plot_candidates$gene, levels = rev(plot_candidates$gene))
direction_colors <- c(
  predicted_decrease = "#2166AC",
  uncertain = "#7F7F7F",
  predicted_increase = "#B2182B"
)
direction_plot <- ggplot(
  plot_candidates,
  aes(x = predicted_ko_change, y = gene, color = predicted_direction)
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_segment(
    data = plot_candidates,
    aes(
      x = bootstrap_ci_low, xend = bootstrap_ci_high,
      y = gene, yend = gene,
      color = predicted_direction
    ),
    inherit.aes = FALSE, linewidth = 0.5, alpha = 0.6
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(values = direction_colors) +
  labs(
    title = "Exploratory expression direction after STK31 loss",
    subtitle = "Reverse of sample-matched STK31-positive vs negative expression; not network KO evidence",
    x = "Predicted KO change in mean log-normalized expression",
    y = NULL, color = "Prediction"
  ) +
  theme_classic(base_size = 11)
write_plot("03_stk31_ko_direction_predictions.pdf", direction_plot, 9, 8)

# 图 4：逐样本热图，查看候选基因方向是否由单一样本驱动。
stability_genes <- as.character(plot_candidates$gene)
stability_plot_data <- direction$stability[
  direction$stability$gene %in% stability_genes,
]
stability_plot_data$gene <- factor(
  stability_plot_data$gene,
  levels = rev(stability_genes)
)
stability_plot <- ggplot(
  stability_plot_data,
  aes(x = sample, y = gene, fill = predicted_ko_change)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  labs(
    title = "Cross-sample stability of predicted STK31-KO direction",
    x = NULL, y = NULL, fill = "Predicted KO\nchange"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
write_plot("04_stk31_ko_sample_stability.pdf", stability_plot, 7, 8)

# 图 5：leave-one-sample-out 稳定性；每次去掉一个样本后重新看方向是否保持。
leave_one_out_plot_data <- direction$leave_one_out[
  direction$leave_one_out$gene %in% stability_genes,
]
leave_one_out_plot_data$gene <- factor(
  leave_one_out_plot_data$gene,
  levels = rev(stability_genes)
)
leave_one_out_plot <- ggplot(
  leave_one_out_plot_data,
  aes(x = held_out_sample, y = gene, fill = predicted_ko_change)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  labs(
    title = "Leave-one-sample-out robustness of STK31-KO predictions",
    x = "Held-out sample", y = NULL, fill = "Predicted KO\nchange"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
write_plot("05_stk31_ko_leave_one_sample_out.pdf", leave_one_out_plot, 7, 8)

# 图 6：机制集合方向汇总，只统计观察性支持的升高/降低基因数，不做通路富集显著性。
pathway_plot_data <- rbind(
  data.frame(
    pathway = pathways$pathway,
    direction = "Predicted increase",
    genes = pathways$predicted_increase
  ),
  data.frame(
    pathway = pathways$pathway,
    direction = "Predicted decrease",
    genes = pathways$predicted_decrease
  )
)
pathway_plot_data$pathway <- factor(
  pathway_plot_data$pathway,
  levels = rev(pathways$pathway)
)
pathway_plot <- ggplot(
  pathway_plot_data,
  aes(x = genes, y = pathway, fill = direction)
) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(
    `Predicted increase` = "#B2182B",
    `Predicted decrease` = "#2166AC"
  )) +
  labs(
    title = "Exploratory direction summary for curated mechanisms",
    subtitle = "Counts require wet-lab validation; this is not pathway enrichment",
    x = "Genes with supported predicted direction", y = NULL, fill = NULL
  ) +
  theme_classic(base_size = 11)
write_plot("06_stk31_ko_mechanism_direction_summary.pdf", pathway_plot, 8, 5.5)

unlink(file.path(out_dir, c(
  "stk31_ko_pathway_enrichment.csv",
  "06_stk31_ko_pathway_enrichment.pdf"
)))

# 【运行记录】记录 R 与已加载包的版本，帮助以后解释同一代码为何可能得到不同结果。
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("STK31 virtual knockout analysis completed: ", out_dir)
