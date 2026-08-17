#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
nk_file <- file.path(project, "01_raw_data", "NK.RDS")
out_dir <- file.path(project, "06_feasibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c(
  "STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2",
  "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRC1", "KLRD1",
  "NKG7", "GNLY", "PRF1", "GZMB", "FGFBP2", "FCGR3A", "XCL1", "XCL2"
)

collapse_examples <- function(x, n = 10L) {
  values <- unique(as.character(x))
  values <- values[!is.na(values)]
  paste(head(values, n), collapse = " | ")
}

if (!file.exists(nk_file)) stop("Missing NK.RDS: ", nk_file)

nk <- readRDS(nk_file)
features <- rownames(nk)
cells <- colnames(nk)

lines <- c(
  paste0("file: ", nk_file),
  paste0("file_size_bytes: ", file.info(nk_file)$size),
  paste0("file_md5: ", unname(tools::md5sum(nk_file))),
  paste0("R_version: ", R.version.string),
  paste0("installed_Seurat_version: ", as.character(packageVersion("Seurat"))),
  paste0("installed_SeuratObject_version: ", as.character(packageVersion("SeuratObject"))),
  paste0("class: ", paste(class(nk), collapse = ", ")),
  paste0("dim_features_cells: ", paste(dim(nk), collapse = " x ")),
  paste0("unique_features: ", uniqueN(features)),
  paste0("unique_cells: ", uniqueN(cells)),
  paste0("duplicated_features: ", sum(duplicated(features))),
  paste0("duplicated_cells: ", sum(duplicated(cells)))
)

raw_counts_present <- FALSE
normalized_data_present <- FALSE
meta <- NULL

if (inherits(nk, "Seurat")) {
  assay_names <- names(nk@assays)
  lines <- c(
    lines,
    paste0("object_version_slot: ", as.character(nk@version)),
    paste0("assays: ", paste(assay_names, collapse = ", ")),
    paste0("DefaultAssay: ", DefaultAssay(nk))
  )
  for (assay_name in assay_names) {
    assay <- nk[[assay_name]]
    layer_names <- tryCatch(Layers(assay), error = function(e) character())
    slot_names <- slotNames(assay)
    lines <- c(
      lines,
      paste0("assay.", assay_name, ".class: ", paste(class(assay), collapse = ", ")),
      paste0("assay.", assay_name, ".dim: ", paste(dim(assay), collapse = " x ")),
      paste0("assay.", assay_name, ".layers: ", paste(layer_names, collapse = ", ")),
      paste0("assay.", assay_name, ".slots: ", paste(slot_names, collapse = ", "))
    )
    if ("counts" %in% layer_names) {
      counts_dim <- dim(LayerData(assay, layer = "counts"))
      raw_counts_present <- all(counts_dim > 0)
      lines <- c(lines, paste0("assay.", assay_name, ".counts_dim: ", paste(counts_dim, collapse = " x ")))
    }
    if ("data" %in% layer_names) {
      data_dim <- dim(LayerData(assay, layer = "data"))
      normalized_data_present <- all(data_dim > 0)
      lines <- c(lines, paste0("assay.", assay_name, ".data_dim: ", paste(data_dim, collapse = " x ")))
    }
  }
  reduction_names <- names(nk@reductions)
  lines <- c(lines, paste0("reductions: ", paste(reduction_names, collapse = ", ")))
  for (reduction_name in reduction_names) {
    emb_dim <- dim(Embeddings(nk, reduction = reduction_name))
    lines <- c(lines, paste0("reduction.", reduction_name, ".dim: ", paste(emb_dim, collapse = " x ")))
  }
  graph_names <- names(nk@graphs)
  neighbor_names <- names(nk@neighbors)
  lines <- c(
    lines,
    paste0("graphs: ", paste(graph_names, collapse = ", ")),
    paste0("neighbors: ", paste(neighbor_names, collapse = ", "))
  )
  meta <- as.data.table(nk[[]], keep.rownames = "cell_name")
} else {
  lines <- c(lines, "recognized_Seurat_object: FALSE")
}

lines <- c(
  lines,
  paste0("raw_counts_present: ", raw_counts_present),
  paste0("normalized_data_present: ", normalized_data_present)
)

if (!is.null(meta)) {
  fields <- setdiff(names(meta), "cell_name")
  lines <- c(
    lines,
    paste0("metadata_rows: ", nrow(meta)),
    paste0("metadata_columns: ", ncol(meta) - 1L),
    paste0("metadata_fields: ", paste(fields, collapse = ", "))
  )
  for (field in fields) {
    lines <- c(
      lines,
      paste0("metadata.", field, ".n_unique: ", uniqueN(meta[[field]])),
      paste0("metadata.", field, ".examples: ", collapse_examples(meta[[field]]))
    )
  }
  sample_fields <- grep("sample|orig.ident", fields, ignore.case = TRUE, value = TRUE)
  patient_fields <- grep("patient", fields, ignore.case = TRUE, value = TRUE)
  site_fields <- grep("site|tissue|organ|source", fields, ignore.case = TRUE, value = TRUE)
  subtype_fields <- grep("subtype|celltype|cell_type|cluster|ident", fields, ignore.case = TRUE, value = TRUE)
  lines <- c(
    lines,
    paste0("candidate_sample_fields: ", paste(sample_fields, collapse = ", ")),
    paste0("candidate_patient_fields: ", paste(patient_fields, collapse = ", ")),
    paste0("candidate_site_fields: ", paste(site_fields, collapse = ", ")),
    paste0("candidate_NK_subtype_fields: ", paste(subtype_fields, collapse = ", "))
  )
  fwrite(meta, file.path(out_dir, "nk_metadata.csv"))
}

writeLines(lines, file.path(out_dir, "nk_rds_structure.txt"), useBytes = TRUE)

presence <- data.table(
  gene = target_genes,
  present_in_full_counts = NA,
  present_in_NK_RDS = target_genes %chin% features
)
fwrite(presence, file.path(out_dir, "target_gene_presence_partial.csv"))

cat(paste(lines[grepl("^(class|dim_features_cells|assays|DefaultAssay|reductions|raw_counts_present|normalized_data_present|candidate_)", lines)], collapse = "\n"), "\n")
cat("NK_STRUCTURE_CHECK_COMPLETE\n")
