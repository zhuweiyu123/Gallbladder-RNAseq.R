# ============================================================
# 10_refined_stk31_high_vs_low_and_all_go.R
# On frozen refined object:
#   1) STK31-high epithelial vs STK31-low epithelial DE + GO BP
#   2) STK31-high epithelial vs all other cells DE + GO BP
# Epithelial from analysis_celltype == "Epithelial"
# STK31-high: epithelial counts > 0 (matches freeze when q-cutoff==0)
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
out_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "STK31"
epi_label <- "Epithelial"
padj_cut <- 0.05
logfc_cut <- 0.25
min_genes_go <- 10
max_cells_per_ident <- as.integer(Sys.getenv("STK31_MAX_CELLS_DE", "500"))

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing: ", path)
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

run_de <- function(object, cells.1, cells.2, out_file) {
  if (length(cells.1) < 10 || length(cells.2) < 10) {
    stop("Too few cells: n1=", length(cells.1), " n2=", length(cells.2))
  }
  object$.cmp <- "unused"
  object$.cmp[colnames(object) %in% cells.1] <- "g1"
  object$.cmp[colnames(object) %in% cells.2] <- "g2"
  message("DE: n1=", length(cells.1), " n2=", length(cells.2),
          " max.cells.per.ident=", max_cells_per_ident)
  markers <- FindMarkers(
    object,
    ident.1 = "g1",
    ident.2 = "g2",
    group.by = ".cmp",
    logfc.threshold = 0.1,
    min.pct = 0.1,
    test.use = "wilcox",
    max.cells.per.ident = max_cells_per_ident
  )
  markers$gene <- rownames(markers)
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
  write.csv(markers, out_file, row.names = FALSE)
  markers
}

run_go_bp <- function(markers, output_prefix, comparison_name) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("clusterProfiler/org.Hs.eg.db missing")
  }
  sig <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      !is.na(markers$avg_log2FC) & abs(markers$avg_log2FC) > logfc_cut
  ]
  sig_up <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      markers$avg_log2FC > logfc_cut
  ]
  sig_down <- markers$gene[
    !is.na(markers$p_val_adj) & markers$p_val_adj < padj_cut &
      markers$avg_log2FC < -logfc_cut
  ]
  sig <- unique(sig[!is.na(sig)])
  write.csv(
    data.frame(
      comparison = comparison_name,
      n_sig = length(sig),
      n_sig_up = length(unique(sig_up[!is.na(sig_up)])),
      n_sig_down = length(unique(sig_down[!is.na(sig_down)])),
      padj_cut = padj_cut,
      logfc_cut = logfc_cut,
      stringsAsFactors = FALSE
    ),
    paste0(output_prefix, "_go_input_summary.csv"),
    row.names = FALSE
  )
  message(comparison_name, ": n_sig=", length(sig),
          " up=", length(unique(sig_up[!is.na(sig_up)])),
          " down=", length(unique(sig_down[!is.na(sig_down)])))

  if (length(sig) < min_genes_go) {
    write.csv(data.frame(note = "skipped_few_genes"), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  suppressMessages({
    converted <- clusterProfiler::bitr(
      sig, fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  })
  if (nrow(converted) < min_genes_go) {
    write.csv(data.frame(note = "skipped_few_entrez"), paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  ego <- clusterProfiler::enrichGO(
    gene = unique(converted$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, paste0(output_prefix, "_go_bp.csv"), row.names = FALSE)
  message(comparison_name, ": n_GO=", nrow(ego_df))

  if (nrow(ego_df) > 0) {
    top <- head(ego_df[order(ego_df$p.adjust), ], 20)
    top$Description <- factor(top$Description, levels = rev(top$Description))
    if (!"Count" %in% colnames(top)) {
      top$Count <- vapply(strsplit(as.character(top$GeneRatio), "/"), function(z) as.numeric(z[1]), numeric(1))
    }
    top$score <- -log10(pmax(top$p.adjust, 1e-300))
    p <- ggplot(top, aes(x = score, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "#d73027", high = "#4575b4", trans = "log10") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0("GO BP: ", comparison_name),
        subtitle = "Refined object; epithelial=analysis_celltype; STK31-high=counts>0",
        x = "-log10(FDR)", y = NULL
      )
    ggsave(paste0(output_prefix, "_go_bp_dotplot.pdf"), p,
           width = 10, height = max(5, 0.28 * nrow(top) + 2))
  }
  invisible(ego_df)
}

top_genes_table <- function(markers, n = 15) {
  up <- head(markers[markers$avg_log2FC > 0, ], n)
  down <- head(markers[markers$avg_log2FC < 0, ], n)
  list(up = up, down = down)
}

message("=== 10 STK31 high vs low / vs all GO ===")
message("Start: ", Sys.time())
stop_if_missing(input_rds)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  try(obj <- JoinLayers(obj, assay = "RNA"), silent = TRUE)
}
if (!"analysis_celltype" %in% colnames(obj@meta.data)) stop("analysis_celltype missing")

stk31_counts <- tryCatch(
  as.numeric(fetch_mat(obj, target_gene, "counts")[target_gene, ]),
  error = function(e) as.numeric(fetch_mat(obj, target_gene, "data")[target_gene, ])
)
is_epi <- obj$analysis_celltype == epi_label
if (sum(is_epi) == 0 && "legacy_manual_celltype" %in% colnames(obj@meta.data)) {
  is_epi <- obj$legacy_manual_celltype == epi_label
}
high_cells <- colnames(obj)[is_epi & stk31_counts > 0]
low_cells <- colnames(obj)[is_epi & !(stk31_counts > 0)]
other_cells <- setdiff(colnames(obj), high_cells)

summary_df <- data.frame(
  item = c(
    "epithelial_cells", "stk31_high_epithelial", "stk31_low_epithelial",
    "all_other_vs_high", "stk31_definition", "object"
  ),
  value = c(
    sum(is_epi), length(high_cells), length(low_cells),
    length(other_cells), "epithelial & counts>0",
    basename(input_rds)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "stk31_high_low_all_de_summary.csv"), row.names = FALSE)
print(summary_df)

message("1) DE high vs low epithelial")
m_hl <- run_de(
  obj, high_cells, low_cells,
  file.path(out_dir, "markers_stk31_high_vs_low_tumor_epithelial.csv")
)
message("2) GO high vs low epithelial")
run_go_bp(
  m_hl,
  file.path(out_dir, "stk31_high_vs_low_tumor_epithelial"),
  "STK31_high_vs_low_epithelial"
)

message("3) DE high epithelial vs all other")
m_ha <- run_de(
  obj, high_cells, other_cells,
  file.path(out_dir, "markers_stk31_high_tumor_epithelial_vs_all_other.csv")
)
message("4) GO high epithelial vs all other")
run_go_bp(
  m_ha,
  file.path(out_dir, "stk31_high_tumor_epithelial_vs_all_other"),
  "STK31_high_epithelial_vs_all_other"
)

# top gene exports for quick reading
write.csv(head(m_hl[m_hl$avg_log2FC > 0, ], 50),
          file.path(out_dir, "top50_up_stk31_high_vs_low_epithelial.csv"), row.names = FALSE)
write.csv(head(m_hl[m_hl$avg_log2FC < 0, ], 50),
          file.path(out_dir, "top50_down_stk31_high_vs_low_epithelial.csv"), row.names = FALSE)
write.csv(head(m_ha[m_ha$avg_log2FC > 0, ], 50),
          file.path(out_dir, "top50_up_stk31_high_vs_all_other.csv"), row.names = FALSE)
write.csv(head(m_ha[m_ha$avg_log2FC < 0, ], 50),
          file.path(out_dir, "top50_down_stk31_high_vs_all_other.csv"), row.names = FALSE)

# brief text summary of top GO
summarize_go_top <- function(path, n = 10) {
  if (!file.exists(path)) return(character(0))
  g <- read.csv(path, stringsAsFactors = FALSE)
  if (!"Description" %in% names(g) || nrow(g) == 0) return(character(0))
  head(g$Description[order(g$p.adjust)], n)
}

hl_go <- summarize_go_top(file.path(out_dir, "stk31_high_vs_low_tumor_epithelial_go_bp.csv"))
ha_go <- summarize_go_top(file.path(out_dir, "stk31_high_tumor_epithelial_vs_all_other_go_bp.csv"))

writeLines(c(
  "# STK31-high epithelial DE/GO (refined object)",
  paste0("Date: ", Sys.time()),
  paste0("High cells: ", length(high_cells), "; Low epithelial: ", length(low_cells),
         "; All other: ", length(other_cells)),
  "",
  "## Meaning",
  "- high vs low: within-epithelial STK31-associated transcriptome (primary for STK31 biology).",
  "- high vs all other: high-epi vs rest of atlas (mixes cell-type identity + STK31); interpret carefully.",
  "",
  "## Top GO high vs low",
  paste0("- ", hl_go),
  "",
  "## Top GO high vs all other",
  paste0("- ", ha_go)
), file.path(out_dir, "stk31_high_low_all_de_go_README.md"))

message("Done: ", out_dir)
message("End: ", Sys.time())
