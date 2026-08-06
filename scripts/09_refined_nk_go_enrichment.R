# ============================================================
# 09_refined_nk_go_enrichment.R
# GO BP for refined high-confidence NK DE comparisons.
# Inputs (already from script 06 Stage2):
#   markers_nk_vs_all_other.csv
#   markers_stk31_high_tumor_epithelial_vs_nk.csv
# NK definition: analysis_celltype == "NK_cell" (n=1551)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
in_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
out_dir <- in_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_cut <- 0.05
logfc_cut <- 0.25
min_genes <- 10

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

run_go_bp <- function(markers, output_prefix, comparison_name) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("clusterProfiler/org.Hs.eg.db not installed")
  }

  if (!"gene" %in% colnames(markers)) {
    if (!is.null(rownames(markers)) && !all(rownames(markers) %in% c(as.character(seq_len(nrow(markers)))))) {
      markers$gene <- rownames(markers)
    } else {
      stop("markers table has no gene column: ", comparison_name)
    }
  }

  if ("p_val_adj" %in% colnames(markers)) {
    sig <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & abs(markers$avg_log2FC) > logfc_cut
    ]
    # also report up/down separately for optional files
    sig_up <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & markers$avg_log2FC > logfc_cut
    ]
    sig_down <- markers$gene[
      !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
        !is.na(markers$avg_log2FC) & markers$avg_log2FC < -logfc_cut
    ]
  } else {
    stop("Expected p_val_adj column in markers for ", comparison_name)
  }

  sig <- unique(sig[!is.na(sig) & nzchar(sig)])
  message(comparison_name, ": n_sig_genes=", length(sig),
          " (up=", length(unique(sig_up)), ", down=", length(unique(sig_down)), ")")

  write.csv(
    data.frame(
      comparison = comparison_name,
      n_sig = length(sig),
      n_sig_up = length(unique(sig_up[!is.na(sig_up)])),
      n_sig_down = length(unique(sig_down[!is.na(sig_down)])),
      padj_cut = padj_cut,
      logfc_cut = logfc_cut,
      nk_definition = "analysis_celltype==NK_cell refined high-confidence",
      stringsAsFactors = FALSE
    ),
    paste0(output_prefix, "_go_input_summary.csv"),
    row.names = FALSE
  )

  if (length(sig) < min_genes) {
    warning("Fewer than ", min_genes, " significant genes; skip GO for ", comparison_name)
    write.csv(
      data.frame(note = paste0("skipped: n_sig=", length(sig))),
      paste0(output_prefix, "_go_bp.csv"),
      row.names = FALSE
    )
    return(invisible(NULL))
  }

  suppressMessages({
    converted <- clusterProfiler::bitr(
      sig,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  })
  message(comparison_name, ": n_entrez=", length(unique(converted$ENTREZID)))
  if (nrow(converted) < min_genes) {
    warning("Fewer than ", min_genes, " converted genes; skip GO for ", comparison_name)
    write.csv(
      data.frame(note = paste0("skipped: n_converted=", nrow(converted))),
      paste0(output_prefix, "_go_bp.csv"),
      row.names = FALSE
    )
    return(invisible(NULL))
  }

  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
  message(comparison_name, ": n_GO_terms=", nrow(ego_df))

  # simple dotplot of top terms
  if (nrow(ego_df) > 0) {
    top <- head(ego_df[order(ego_df$p.adjust), ], 20)
    top$Description <- factor(top$Description, levels = rev(top$Description))
    if (!"Count" %in% colnames(top)) {
      top$Count <- vapply(strsplit(top$GeneRatio, "/"), function(z) as.numeric(z[1]), numeric(1))
    }
    top$score <- -log10(pmax(top$p.adjust, 1e-300))
    p <- ggplot(top, aes(x = score, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#d73027", high = "#4575b4", trans = "log10") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0("GO BP: ", comparison_name),
        subtitle = "Refined high-confidence NK (analysis_celltype==NK_cell)",
        x = "-log10(adjusted P)",
        y = NULL,
        size = "Genes",
        color = "FDR"
      )
    ggsave(paste0(output_prefix, "_go_bp_dotplot.pdf"), p, width = 10, height = max(5, 0.28 * nrow(top) + 2))
  }

  invisible(ego_df)
}

message("=== 09 refined NK GO enrichment ===")
message("Start: ", Sys.time())

nk_file <- file.path(in_dir, "markers_nk_vs_all_other.csv")
stk31_nk_file <- file.path(in_dir, "markers_stk31_high_tumor_epithelial_vs_nk.csv")
stop_if_missing(nk_file)
stop_if_missing(stk31_nk_file)

nk_markers <- read.csv(nk_file, stringsAsFactors = FALSE, check.names = FALSE)
stk31_nk_markers <- read.csv(stk31_nk_file, stringsAsFactors = FALSE, check.names = FALSE)

message("Running GO: refined NK vs all other")
run_go_bp(
  nk_markers,
  file.path(out_dir, "nk_vs_all_other"),
  "refined_NK_vs_all_other"
)

message("Running GO: STK31-high epithelial vs refined NK")
run_go_bp(
  stk31_nk_markers,
  file.path(out_dir, "stk31_high_tumor_epithelial_vs_nk"),
  "STK31_high_epithelial_vs_refined_NK"
)

# inventory note
note <- c(
  "# Refined NK GO outputs",
  paste0("Date: ", Sys.time()),
  "NK definition: analysis_celltype == NK_cell (refined high-confidence, n=1551)",
  "DE inputs:",
  paste0("  - ", nk_file),
  paste0("  - ", stk31_nk_file),
  paste0("Thresholds: p_val_adj < ", padj_cut, ", |avg_log2FC| > ", logfc_cut),
  "Outputs:",
  "  - nk_vs_all_other_go_bp.csv / _go_bp_dotplot.pdf / _go_input_summary.csv",
  "  - stk31_high_tumor_epithelial_vs_nk_go_bp.csv / _go_bp_dotplot.pdf / _go_input_summary.csv",
  "Legacy GO under merged_stk31_nk_analysis/ used old NK=7419 and must not be used for primary inference."
)
writeLines(note, file.path(out_dir, "refined_nk_go_README.md"))

message("Done. Outputs in: ", out_dir)
message("End: ", Sys.time())
