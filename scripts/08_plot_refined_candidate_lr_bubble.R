# Plot custom-panel candidate L-R bubble using refined NK results.
# Input: results/merged_stk31_nk_refined_analysis/candidate_ligand_receptor_pathways.csv
# Output: same directory + note that this is user-curated panel, not CellChatDB.

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20260717)

project_dir <- "/home/zhuweiyu/codex-r"
in_file <- file.path(
  project_dir,
  "results/merged_stk31_nk_refined_analysis/candidate_ligand_receptor_pathways.csv"
)
out_dir <- file.path(project_dir, "results/merged_stk31_nk_refined_analysis")
top_n <- as.integer(Sys.getenv("STK31_LR_BUBBLE_TOP_N", "30"))

if (!file.exists(in_file)) stop("Missing: ", in_file)

lr <- read.csv(in_file, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("ligand", "receptor", "pathway", "direction", "interaction_score")
miss <- setdiff(need, names(lr))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

lr <- lr[is.finite(lr$interaction_score), , drop = FALSE]
lr <- lr[order(-lr$interaction_score), , drop = FALSE]
lr_top <- head(lr, top_n)
lr_top$pair <- paste(lr_top$ligand, lr_top$receptor, sep = " -> ")
# keep order by score for y axis
lr_top$pair <- factor(lr_top$pair, levels = rev(unique(lr_top$pair)))

# shorter direction labels for axis
lr_top$direction_short <- lr_top$direction
lr_top$direction_short <- gsub(
  "STK31_high_tumor_epithelial_to_NK_cell",
  "STK31-high epi -> refined NK",
  lr_top$direction_short
)
lr_top$direction_short <- gsub(
  "NK_cell_to_STK31_high_tumor_epithelial",
  "refined NK -> STK31-high epi",
  lr_top$direction_short
)

p <- ggplot(
  lr_top,
  aes(x = direction_short, y = pair, size = interaction_score, color = pathway)
) +
  geom_point(alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = "Top candidate ligand-receptor pairs (user-curated panel)",
    subtitle = paste0(
      "Refined NK (analysis_celltype==NK_cell, n~1551); ",
      "NOT CellChatDB. Top ", top_n, " by interaction_score."
    ),
    x = "Direction",
    y = "Ligand -> Receptor",
    size = "Score",
    color = "Pathway"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

out_pdf <- file.path(out_dir, "refined_candidate_lr_bubble.pdf")
ggsave(out_pdf, p, width = 11, height = 8)
message("Wrote: ", out_pdf)

# also pathway-level bubble/bar helper from summary if present
sum_file <- file.path(out_dir, "candidate_pathway_summary.csv")
if (file.exists(sum_file)) {
  ps <- read.csv(sum_file, stringsAsFactors = FALSE)
  ps <- ps[order(-ps$interaction_score), , drop = FALSE]
  ps_top <- head(ps, 20)
  ps_top$pathway <- factor(ps_top$pathway, levels = rev(unique(ps_top$pathway)))
  ps_top$direction_short <- gsub(
    "STK31_high_tumor_epithelial_to_NK_cell",
    "STK31-high epi -> refined NK",
    ps_top$direction
  )
  ps_top$direction_short <- gsub(
    "NK_cell_to_STK31_high_tumor_epithelial",
    "refined NK -> STK31-high epi",
    ps_top$direction_short
  )
  p2 <- ggplot(
    ps_top,
    aes(x = direction_short, y = pathway, size = interaction_score, color = pathway)
  ) +
    geom_point(alpha = 0.9) +
    theme_bw(base_size = 12) +
    labs(
      title = "Candidate pathway summary (user-curated panel, refined NK)",
      subtitle = "Max interaction_score per pathway x direction",
      x = "Direction", y = "Pathway", size = "Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )
  out2 <- file.path(out_dir, "refined_candidate_pathway_bubble.pdf")
  ggsave(out2, p2, width = 10, height = 7)
  message("Wrote: ", out2)
}

writeLines(
  c(
    "# refined_candidate_lr_bubble",
    "",
    paste0("Date: ", Sys.time()),
    paste0("Input: ", in_file),
    paste0("Top N pairs: ", top_n),
    "",
    "This figure uses the **user-curated** ligand-receptor panel",
    "(scripts/06 lr_reference), scored on STK31-high epithelial vs **refined NK**.",
    "It is NOT generated from CellChatDB.",
    "For CellChat high-vs-low bubble see:",
    "results/merged_cellchat_refined_high_low_nk/cellchat_high_low_nk_bubble.pdf"
  ),
  file.path(out_dir, "refined_candidate_lr_bubble_README.md")
)

message("Done.")
