#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
processed <- file.path(project, "02_processed_data", "All_single_cells")
counts <- file.path(processed, "counts", "10X_counts")
out_dir <- file.path(project, "06_feasibility")

target_genes <- c(
  "STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2",
  "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRC1", "KLRD1",
  "NKG7", "GNLY", "PRF1", "GZMB", "FGFBP2", "FCGR3A", "XCL1", "XCL2"
)

meta <- fread(
  file.path(processed, "GBC_Metadata.txt"),
  sep = "\t",
  select = c("name", "celltype", "subtype", "sample_name"),
  showProgress = FALSE
)
schema <- meta$name == "type" & meta$celltype == "group" & meta$subtype == "group"
meta <- meta[!schema]

barcodes <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(counts, "barcodes.tsv.gz"))),
  header = FALSE,
  col.names = "barcode",
  showProgress = FALSE
)
features <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(counts, "features.tsv.gz"))),
  header = FALSE,
  col.names = c("feature_id", "feature_name", "feature_type"),
  showProgress = FALSE
)

nk <- readRDS(file.path(project, "01_raw_data", "NK.RDS"))
nk_cells <- colnames(nk)
nk_features <- rownames(nk)
full_nk_cells <- meta[celltype == "NK cells", name]

nk_meta <- as.data.table(nk[[]], keep.rownames = "cell_name")
sample_field_details <- character()
for (field in intersect(c("orig.ident", "sample.ident"), names(nk_meta))) {
  vals <- sort(unique(as.character(nk_meta[[field]])))
  full_samples <- sort(unique(meta$sample_name))
  sample_field_details <- c(
    sample_field_details,
    paste0("nk.", field, ".n_unique_nonmissing: ", length(vals)),
    paste0("nk.", field, ".n_missing: ", sum(is.na(nk_meta[[field]]))),
    paste0("nk.", field, ".not_in_full_metadata: ", paste(setdiff(vals, full_samples), collapse = " | ")),
    paste0("full_metadata.not_in_nk.", field, ": ", paste(setdiff(full_samples, vals), collapse = " | "))
  )
}

audit <- c(
  paste0("metadata_cells: ", nrow(meta)),
  paste0("barcode_rows: ", nrow(barcodes)),
  paste0("metadata_name_duplicates: ", sum(duplicated(meta$name))),
  paste0("barcode_duplicates: ", sum(duplicated(barcodes$barcode))),
  paste0("metadata_names_identical_to_barcodes_same_order: ", identical(meta$name, barcodes$barcode)),
  paste0("metadata_names_not_in_barcodes: ", sum(!meta$name %chin% barcodes$barcode)),
  paste0("barcodes_not_in_metadata_names: ", sum(!barcodes$barcode %chin% meta$name)),
  paste0("full_metadata_NK_cells: ", length(full_nk_cells)),
  paste0("NK_RDS_cells: ", length(nk_cells)),
  paste0("NK_RDS_cell_duplicates: ", sum(duplicated(nk_cells))),
  paste0("full_NK_names_identical_to_NK_RDS_same_order: ", identical(full_nk_cells, nk_cells)),
  paste0("full_NK_names_not_in_NK_RDS: ", sum(!full_nk_cells %chin% nk_cells)),
  paste0("NK_RDS_names_not_in_full_NK: ", sum(!nk_cells %chin% full_nk_cells)),
  paste0("feature_rows: ", nrow(features)),
  paste0("feature_name_duplicates: ", sum(duplicated(features$feature_name))),
  paste0("feature_types: ", paste(names(table(features$feature_type)), table(features$feature_type), sep = "=", collapse = " | ")),
  sample_field_details
)

writeLines(audit, file.path(out_dir, "identifier_mapping_audit.txt"), useBytes = TRUE)

presence <- data.table(
  gene = target_genes,
  present_in_full_counts = target_genes %chin% features$feature_name,
  present_in_NK_RDS = target_genes %chin% nk_features
)
fwrite(presence, file.path(out_dir, "target_gene_presence.csv"))

cat(paste(audit, collapse = "\n"), "\n")
cat("\nTARGET_GENE_PRESENCE\n")
print(presence)
