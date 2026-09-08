#!/usr/bin/env Rscript

# CellChat pathway analysis for the public paired-QC4 cohort and four local
# gallbladder samples. The public cohort is treated as one analysis unit;
# local samples are analysed separately. Cell labels and per-group sampling
# are harmonized before inference. WNT is exported as a dedicated result set.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(CellChat)
})

set.seed(20260827)

project_dir <- "/home/zhuweiyu/codex-r"
public_dir <- file.path(project_dir, "GBC_STK31_public_scRNA")
local_rds <- file.path(
  project_dir, "results/merged_cnv_binary_annotation/merged_cnv_binary_annotated.rds"
)
out_dir <- file.path(project_dir, "results/qc4_local_sample_cellchat_pathways")
combined_dir <- file.path(out_dir, "00_combined")
qc4_input_dir <- file.path(out_dir, "00_qc4_downsampled_10x")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc4_input_dir, recursive = TRUE, showWarnings = FALSE)

qc4_samples <- c("GBC_033_P", "GBC_047_P", "GBC_056_P", "GBC_073_P")
local_samples <- c("tissue1", "tissue2", "tissue4", "tissue5")
dataset_order <- c("QC4_public", local_samples)
min_cells <- 10L
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

stop_if_not <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

write_note <- function(path, text) {
  writeLines(text, path, useBytes = TRUE)
}

write_placeholder_pdf <- function(path, title, subtitle = NULL) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = title, size = 5) +
    annotate("text", x = 0, y = -0.15, label = subtitle %||% "", size = 3.5) +
    xlim(-1, 1) + ylim(-0.5, 0.5) + theme_void()
  ggsave(path, p, width = 7, height = 4, bg = "white")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

safe_pdf <- function(path, expr, width = 8, height = 8) {
  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

save_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 320, bg = "white")
}

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

sample_up_to <- function(index, n) {
  if (length(index) <= n) index else sample(index, n)
}

balanced_qc4_indices <- function(meta, max_per_group) {
  selected <- integer(0)
  for (group in group_levels) {
    group_index <- which(meta$harmonized_celltype == group)
    if (length(group_index) == 0L) next
    by_patient <- split(group_index, meta$sample_name[group_index])
    by_patient <- by_patient[qc4_samples[qc4_samples %in% names(by_patient)]]
    base_quota <- max(1L, floor(max_per_group / length(qc4_samples)))
    first_pass <- unlist(
      lapply(by_patient, sample_up_to, n = base_quota), use.names = FALSE
    )
    remaining_quota <- max_per_group - length(first_pass)
    remaining <- setdiff(group_index, first_pass)
    fill <- if (remaining_quota > 0L) sample_up_to(remaining, remaining_quota) else integer(0)
    selected <- c(selected, first_pass, fill)
  }
  sort(unique(selected))
}

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
  fwrite(data.table(name = qc_meta$name, labels = qc_meta$harmonized_celltype),
         file.path(qc4_input_dir, "cellchat_metadata.csv.gz"))
  rm(meta, barcodes, features, col_map, row_map)
  gc()
  invisible(qc_meta)
}

read_qc4_counts <- function() {
  counts <- Read10X(data.dir = qc4_input_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[["Gene Expression"]] %||% counts[[1L]]
  meta <- fread(file.path(qc4_input_dir, "cellchat_metadata.csv.gz"))
  stop_if_not(identical(colnames(counts), meta$name), "QC4 matrix and metadata order differ")
  list(counts = counts, labels = meta$labels, cells = meta$name)
}

downsample_local <- function(object, sample_name) {
  index <- which(object$sample_id == sample_name & !is.na(object$harmonized_celltype))
  selected <- unlist(lapply(group_levels, function(group) {
    group_index <- index[object$harmonized_celltype[index] == group]
    sample_up_to(group_index, max_cells_per_group)
  }), use.names = FALSE)
  selected <- sort(unique(selected))
  list(cells = colnames(object)[selected], labels = object$harmonized_celltype[selected])
}

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

plot_wnt_outputs <- function(cellchat, dataset, dataset_dir) {
  wnt_dir <- file.path(dataset_dir, "WNT")
  dir.create(wnt_dir, recursive = TRUE, showWarnings = FALSE)
  detected <- any(toupper(cellchat@netP$pathways) == "WNT")
  status <- data.frame(
    dataset = dataset,
    WNT_detected = detected,
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

run_cellchat <- function(counts, cells, labels, dataset) {
  dataset_dir <- file.path(out_dir, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
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

  seurat <- CreateSeuratObject(
    counts = counts,
    meta.data = data.frame(labels = labels, row.names = cells)
  )
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
  data_input <- LayerData(seurat[["RNA"]], layer = "data")
  meta <- data.frame(labels = labels, row.names = cells, stringsAsFactors = FALSE)
  cellchat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat, type = "triMean", raw.use = TRUE, population.size = FALSE
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
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
    z <- t(vapply(seq_len(nrow(matrix)), function(i) {
      if (is.finite(row_sd[i]) && row_sd[i] > 0) {
        as.numeric(scale(matrix[i, ]))
      } else {
        rep(0, ncol(matrix))
      }
    }, numeric(ncol(matrix))))
    colnames(z) <- colnames(matrix)
    rownames(z) <- rownames(matrix)
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
local_object <- readRDS(local_rds)
stop_if_not(all(c("sample_id", "analysis_celltype_cnv_binary") %in% colnames(local_object[[]])),
            "Required local annotation columns are absent")
local_object$harmonized_celltype <- harmonize_local(local_object$analysis_celltype_cnv_binary)
if (inherits(local_object[["RNA"]], "Assay5")) {
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
