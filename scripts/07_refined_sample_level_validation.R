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

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
input_rds <- file.path(
  project_dir, "results/merged_tnk_refined_annotation/refined_merged_object.rds"
)
out_dir <- Sys.getenv(
  "STK31_SAMPLE_VAL_OUT",
  file.path(project_dir, "results/refined_sample_level_validation_v2")
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
epi_label <- "Epithelial"
n_boot <- as.integer(Sys.getenv("STK31_SAMPLE_BOOT", "200"))
min_high_cells <- 3
min_nk_cells <- 10
min_low_cells <- 20
# Detection thresholds (explicit)
min_mean_expr <- 0
min_pct_expr <- 0.01
similar_abs <- 0.02

# Direction-aware axes (sender / receiver)
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

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

available_genes <- function(object, genes) intersect(genes, rownames(object))

fetch_mat <- function(object, genes, slot = "data") {
  genes <- available_genes(object, genes)
  if (length(genes) == 0) return(NULL)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(object, assay = DefaultAssay(object), layer = slot)[genes, , drop = FALSE]
  } else {
    GetAssayData(object, assay = DefaultAssay(object), slot = slot)[genes, , drop = FALSE]
  }
}

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
  gene_mean <- as.numeric(Matrix::rowMeans(mat))
  gene_pct <- as.numeric(Matrix::rowMeans(mat > 0))
  names(gene_mean) <- genes_present
  names(gene_pct) <- genes_present
  detected <- genes_present[
    (gene_mean > min_mean_expr) & (gene_pct >= min_pct_expr)
  ]
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

is_detected_module <- function(mod) {
  if (is.null(mod) || length(mod$genes_present) == 0) return(FALSE)
  if (!is.finite(mod$mean_expr) || !is.finite(mod$pct_pos)) return(FALSE)
  # module-level: at least one gene detected AND module mean/pct above floor
  length(mod$genes_detected) > 0 &&
    (mod$mean_expr > min_mean_expr) &&
    (mod$pct_pos >= min_pct_expr || mod$mean_expr > 0)
}

label_direction <- function(delta) {
  if (!is.finite(delta)) return("insufficient")
  if (abs(delta) < similar_abs) return("similar")
  if (delta > 0) return("higher_in_high")
  "higher_in_low"
}

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
  print(sessionInfo())
  sink()
}

message("=== 07 refined sample-level validation v2 ===")
message("Start: ", Sys.time())
write_session(file.path(out_dir, "sessionInfo.txt"))
stop_if_missing(input_rds)

message("Loading frozen refined object")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
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

# ---------- group counts ----------
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
write.csv(group_counts, file.path(out_dir, "sample_level_group_counts.csv"), row.names = FALSE)

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

    # Compute sender/receiver modules by direction
    if (sender == "epithelial" && receiver == "NK") {
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

# ---------- within-sample cell bootstrap (fast: pre-extract cell scores) ----------
message("Running descriptive within-sample cell bootstrap (n=", n_boot, ")")

# Precompute per-cell module scores for all unique gene sets
all_gene_sets <- unique(unlist(lapply(axes, function(a) list(
  paste(sort(a$ligands), collapse = "|"),
  paste(sort(a$receptors), collapse = "|")
))))
# map gene-set key -> named numeric vector of cell scores
cell_score_cache <- list()
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

# ---------- STK31 definition sensitivity with set identity ----------
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
