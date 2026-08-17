#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
processed <- file.path(project, "02_processed_data", "All_single_cells")
raw_dir <- file.path(project, "01_raw_data")
out_dir <- file.path(project, "06_feasibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(processed, "GBC_Metadata.txt")
sample_file <- file.path(processed, "sample_info.xlsx")
nk_file <- file.path(raw_dir, "NK.RDS")

required <- c(metadata_file, sample_file, nk_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required files: ", paste(missing, collapse = ", "))

target_genes <- c(
  "STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2",
  "KIR3DL1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KLRC1", "KLRD1",
  "NKG7", "GNLY", "PRF1", "GZMB", "FGFBP2", "FCGR3A", "XCL1", "XCL2"
)

collapse_examples <- function(x, n = 5L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  paste(head(x, n), collapse = " | ")
}

derive_patient_id <- function(x) {
  ans <- sub("^([[:alpha:]]+_[0-9]+).*$", "\\1", x)
  ans[!grepl("^[[:alpha:]]+_[0-9]+", x)] <- NA_character_
  ans
}

standardize_site <- function(site, histology) {
  site <- trimws(as.character(site))
  histology <- trimws(as.character(histology))
  out <- rep("Unknown", length(site))
  is_carcinoma <- grepl("carcinoma", histology, ignore.case = TRUE)
  out[site == "Primary" & is_carcinoma] <- "Primary tumor"
  out[site == "Primary" & !is_carcinoma] <- "Other"
  out[site == "Liver invasion"] <- "Liver invasion"
  out[site == "Liver metastatic lesion"] <- "Liver metastasis"
  out[site == "Lymphonodus metastatic lesion"] <- "Lymph node metastasis"
  out[site == "Omentum"] <- "Omental metastasis"
  out[site %in% c("PBMC", "Primary polyp")] <- "Other"
  out[grepl("adjacent|NAT", site, ignore.case = TRUE)] <- "Adjacent/NAT"
  out
}

message("Reading cell metadata with data.table::fread")
meta <- fread(metadata_file, sep = "\t", header = TRUE, showProgress = TRUE)
schema_rows <- meta[name == "type" & celltype == "group" & subtype == "group"]
if (nrow(schema_rows)) meta <- meta[!(name == "type" & celltype == "group" & subtype == "group")]
setnames(meta, trimws(names(meta)))

required_meta_cols <- c("name", "celltype", "subtype", "sample_name")
if (!all(required_meta_cols %in% names(meta))) {
  stop("Metadata columns differ from expected: ", paste(names(meta), collapse = ", "))
}

metadata_summary <- rbindlist(list(
  data.table(metric = "rows_excluding_schema_row", value = as.character(nrow(meta)), details = ""),
  data.table(metric = "columns", value = as.character(ncol(meta)), details = paste(names(meta), collapse = " | ")),
  data.table(metric = "schema_rows_removed", value = as.character(nrow(schema_rows)), details = ""),
  data.table(metric = "unique_sample_name", value = as.character(uniqueN(meta$sample_name)), details = ""),
  data.table(metric = "duplicated_name_rows", value = as.character(sum(duplicated(meta$name))), details = ""),
  data.table(
    metric = paste0("column_examples:", names(meta)),
    value = vapply(meta, collapse_examples, character(1)),
    details = "up to 5 unique values"
  )
), use.names = TRUE, fill = TRUE)

celltype_counts <- meta[, .(n_cells = .N), by = celltype][order(-n_cells, celltype)]
subtype_counts <- meta[, .(n_cells = .N), by = .(celltype, subtype)][order(celltype, -n_cells, subtype)]

sample_long <- meta[, .(n_cells = .N), by = .(sample_name, celltype)]
sample_wide <- dcast(sample_long, sample_name ~ celltype, value.var = "n_cells", fill = 0L)
sample_total <- meta[, .(total_cells = .N), by = sample_name]
sample_special <- meta[, .(
  malignant_epithelial_cells = sum(subtype == "Malignant epithelial cells"),
  nk_cells = sum(celltype == "NK cells")
), by = sample_name]
sample_counts <- Reduce(
  function(x, y) merge(x, y, by = "sample_name", all = TRUE),
  list(sample_total, sample_special, sample_wide)
)
setorder(sample_counts, sample_name)

metadata_summary <- rbind(
  metadata_summary,
  data.table(
    metric = c(
      "malignant_epithelial_cells", "nk_cells", "samples_with_malignant_epithelial",
      "samples_with_nk", "samples_with_both_malignant_and_nk"
    ),
    value = as.character(c(
      sum(meta$subtype == "Malignant epithelial cells"),
      sum(meta$celltype == "NK cells"),
      sum(sample_special$malignant_epithelial_cells > 0),
      sum(sample_special$nk_cells > 0),
      sum(sample_special$malignant_epithelial_cells > 0 & sample_special$nk_cells > 0)
    )),
    details = ""
  ),
  fill = TRUE
)

message("Reading sample_info.xlsx")
sample_info <- as.data.table(read_excel(sample_file, skip = 1, .name_repair = "unique"))
sample_info <- sample_info[!is.na(`Sample ID`) & trimws(as.character(`Sample ID`)) != ""]
sample_info[, `Sample ID` := as.character(`Sample ID`)]
sample_info[, sample_info_present := TRUE]
sample_info[, patient_id := derive_patient_id(`Sample ID`)]
sample_info[, standardized_site := standardize_site(Site, `Histological type`)]

mapping <- data.table(sample_name = sort(unique(meta$sample_name)))
mapping <- merge(mapping, sample_counts, by = "sample_name", all.x = TRUE)
mapping <- merge(
  mapping,
  sample_info,
  by.x = "sample_name",
  by.y = "Sample ID",
  all.x = TRUE,
  suffixes = c("", "_sample_info")
)
mapping[, sample_id := sample_name]
mapping[, mapping_status := fifelse(
  !is.na(sample_info_present) & sample_info_present,
  "Matched",
  "Unmatched_metadata_to_sample_info"
)]

sample_info_only <- sample_info[!`Sample ID` %in% unique(meta$sample_name)]

patient_sample_mapping <- mapping[, .(
  patient_id,
  sample_id,
  sample_name,
  site = Site,
  standardized_site,
  mapping_status,
  malignant_epithelial_cells,
  nk_cells
)]
setorder(patient_sample_mapping, patient_id, sample_id)

keep_cols <- intersect(c(
  "sample_name", "sample_id", "patient_id", "Site", "standardized_site",
  "TNM stage", "Histological type", "Liver metastasis", "Lymph node metastasis",
  "Metastatic group", "OS_month", "Event", "scRNA-seq", "mapping_status",
  "total_cells", "malignant_epithelial_cells", "nk_cells"
), names(mapping))
sample_site_mapping <- mapping[, ..keep_cols]

primary <- mapping[standardized_site == "Primary tumor"]
primary_patient <- primary[, .(
  malignant_epithelial_cells = sum(malignant_epithelial_cells, na.rm = TRUE),
  nk_cells = sum(nk_cells, na.rm = TRUE),
  n_primary_samples = .N
), by = patient_id]
primary_both_positive <- primary_patient[
  malignant_epithelial_cells > 0 & nk_cells > 0
]

# Descriptive sensitivity table only. No threshold is declared "sufficient".
threshold_sensitivity <- CJ(
  malignant_min_cells = c(1L, 10L, 20L, 50L, 100L),
  nk_min_cells = c(1L, 10L, 20L, 50L)
)
threshold_sensitivity[, n_primary_patients := mapply(
  function(malignant_min, nk_min) {
    uniqueN(primary_patient[
      malignant_epithelial_cells >= malignant_min & nk_cells >= nk_min,
      patient_id
    ])
  },
  malignant_min_cells,
  nk_min_cells
)]
threshold_sensitivity[, interpretation := paste0(
  "descriptive only; malignant >= ", malignant_min_cells,
  ", NK >= ", nk_min_cells,
  "; not a predefined sufficiency criterion"
)]

metadata_summary <- rbind(
  metadata_summary,
  data.table(
    metric = c(
      "sample_info_rows", "metadata_samples_matched_to_sample_info",
      "metadata_samples_unmatched_to_sample_info", "sample_info_rows_without_metadata_cells",
      "primary_tumor_patients", "primary_tumor_patients_with_both_celltypes_gt0",
      "primary_tumor_patients_with_sufficient_cells"
    ),
    value = as.character(c(
      nrow(sample_info), sum(mapping$mapping_status == "Matched"),
      sum(mapping$mapping_status != "Matched"), nrow(sample_info_only),
      uniqueN(primary$patient_id), uniqueN(primary_both_positive$patient_id), "TO_BE_CONFIRMED"
    )),
    details = c(
      "", "", "", "", "", "malignant >0 and NK >0",
      "minimum cell threshold was not specified; no threshold was guessed"
    )
  ),
  fill = TRUE
)

fwrite(metadata_summary, file.path(out_dir, "metadata_summary.csv"))
fwrite(celltype_counts, file.path(out_dir, "celltype_counts.csv"))
fwrite(subtype_counts, file.path(out_dir, "subtype_counts.csv"))
fwrite(sample_counts, file.path(out_dir, "sample_celltype_counts.csv"))
fwrite(patient_sample_mapping, file.path(out_dir, "patient_sample_mapping.csv"))
fwrite(sample_site_mapping, file.path(out_dir, "sample_site_mapping.csv"))
fwrite(sample_info_only, file.path(out_dir, "sample_info_without_metadata_cells.csv"))
fwrite(primary_patient, file.path(out_dir, "primary_patient_cell_counts.csv"))
fwrite(threshold_sensitivity, file.path(out_dir, "primary_patient_threshold_sensitivity.csv"))

message("Reading NK.RDS for structure inspection only")
nk <- readRDS(nk_file)
nk_class <- class(nk)
nk_features <- rownames(nk)
nk_cells <- colnames(nk)

structure_lines <- c(
  paste0("file: ", nk_file),
  paste0("file_size_bytes: ", file.info(nk_file)$size),
  paste0("R_version: ", R.version.string),
  paste0("installed_Seurat_version: ", if (requireNamespace("Seurat", quietly = TRUE)) as.character(packageVersion("Seurat")) else "NOT_INSTALLED"),
  paste0("installed_SeuratObject_version: ", if (requireNamespace("SeuratObject", quietly = TRUE)) as.character(packageVersion("SeuratObject")) else "NOT_INSTALLED"),
  paste0("class: ", paste(nk_class, collapse = ", ")),
  paste0("dim_features_cells: ", paste(dim(nk), collapse = " x ")),
  paste0("n_unique_features: ", uniqueN(nk_features)),
  paste0("n_unique_cells: ", uniqueN(nk_cells)),
  paste0("duplicated_features: ", sum(duplicated(nk_features))),
  paste0("duplicated_cells: ", sum(duplicated(nk_cells)))
)

nk_meta <- NULL
nk_assays <- character()
nk_default_assay <- NA_character_
raw_counts_present <- FALSE
normalized_data_present <- FALSE

if (inherits(nk, "Seurat")) {
  nk_assays <- names(nk@assays)
  nk_default_assay <- SeuratObject::DefaultAssay(nk)
  structure_lines <- c(
    structure_lines,
    paste0("object_version_slot: ", as.character(nk@version)),
    paste0("assays: ", paste(nk_assays, collapse = ", ")),
    paste0("DefaultAssay: ", nk_default_assay)
  )
  for (assay_name in nk_assays) {
    assay <- nk[[assay_name]]
    layers <- tryCatch(SeuratObject::Layers(assay), error = function(e) character())
    structure_lines <- c(
      structure_lines,
      paste0("assay.", assay_name, ".class: ", paste(class(assay), collapse = ", ")),
      paste0("assay.", assay_name, ".dim: ", paste(dim(assay), collapse = " x ")),
      paste0("assay.", assay_name, ".layers: ", paste(layers, collapse = ", "))
    )
    if ("counts" %in% layers) {
      counts_dim <- tryCatch(dim(SeuratObject::LayerData(assay, layer = "counts")), error = function(e) c(NA, NA))
      raw_counts_present <- raw_counts_present || all(!is.na(counts_dim)) && prod(counts_dim) > 0
      structure_lines <- c(structure_lines, paste0("assay.", assay_name, ".counts_dim: ", paste(counts_dim, collapse = " x ")))
    }
    if ("data" %in% layers) {
      data_dim <- tryCatch(dim(SeuratObject::LayerData(assay, layer = "data")), error = function(e) c(NA, NA))
      normalized_data_present <- normalized_data_present || all(!is.na(data_dim)) && prod(data_dim) > 0
      structure_lines <- c(structure_lines, paste0("assay.", assay_name, ".data_dim: ", paste(data_dim, collapse = " x ")))
    }
  }
  reductions <- names(nk@reductions)
  structure_lines <- c(structure_lines, paste0("reductions: ", paste(reductions, collapse = ", ")))
  for (reduction in reductions) {
    structure_lines <- c(
      structure_lines,
      paste0("reduction.", reduction, ".dim: ", paste(dim(SeuratObject::Embeddings(nk, reduction = reduction)), collapse = " x "))
    )
  }
  nk_meta <- nk[[]]
} else if (inherits(nk, "SingleCellExperiment")) {
  nk_assays <- SummarizedExperiment::assayNames(nk)
  raw_counts_present <- "counts" %in% nk_assays
  normalized_data_present <- any(c("logcounts", "data") %in% nk_assays)
  structure_lines <- c(structure_lines, paste0("assays: ", paste(nk_assays, collapse = ", ")))
  nk_meta <- as.data.frame(SummarizedExperiment::colData(nk))
} else {
  structure_lines <- c(structure_lines, "recognized_single_cell_object: FALSE")
}

structure_lines <- c(
  structure_lines,
  paste0("raw_counts_present: ", raw_counts_present),
  paste0("normalized_data_present: ", normalized_data_present)
)

if (!is.null(nk_meta)) {
  nk_meta_dt <- as.data.table(nk_meta, keep.rownames = "cell_name")
  meta_fields <- setdiff(names(nk_meta_dt), "cell_name")
  structure_lines <- c(
    structure_lines,
    paste0("metadata_n_rows: ", nrow(nk_meta_dt)),
    paste0("metadata_fields: ", paste(meta_fields, collapse = ", "))
  )
  for (field in meta_fields) {
    structure_lines <- c(
      structure_lines,
      paste0("metadata.", field, ".n_unique: ", uniqueN(nk_meta_dt[[field]])),
      paste0("metadata.", field, ".examples: ", collapse_examples(nk_meta_dt[[field]]))
    )
  }
  sample_fields <- grep("sample|orig.ident", meta_fields, ignore.case = TRUE, value = TRUE)
  patient_fields <- grep("patient", meta_fields, ignore.case = TRUE, value = TRUE)
  site_fields <- grep("site|tissue|organ|source", meta_fields, ignore.case = TRUE, value = TRUE)
  subtype_fields <- grep("subtype|celltype|cell_type|cluster|ident", meta_fields, ignore.case = TRUE, value = TRUE)
  structure_lines <- c(
    structure_lines,
    paste0("candidate_sample_fields: ", paste(sample_fields, collapse = ", ")),
    paste0("candidate_patient_fields: ", paste(patient_fields, collapse = ", ")),
    paste0("candidate_site_fields: ", paste(site_fields, collapse = ", ")),
    paste0("candidate_NK_subtype_fields: ", paste(subtype_fields, collapse = ", "))
  )
  for (field in sample_fields) {
    vals <- unique(as.character(nk_meta_dt[[field]]))
    structure_lines <- c(
      structure_lines,
      paste0("sample_field.", field, ".overlap_with_full_metadata: ", sum(vals %in% unique(meta$sample_name)), "/", length(vals)),
      paste0("sample_field.", field, ".overlap_with_sample_info: ", sum(vals %in% sample_info$`Sample ID`), "/", length(vals))
    )
  }
}

writeLines(structure_lines, file.path(out_dir, "nk_rds_structure.txt"), useBytes = TRUE)

target_presence <- data.table(
  gene = target_genes,
  present_in_full_counts = NA,
  present_in_NK_RDS = target_genes %chin% nk_features
)
fwrite(target_presence, file.path(out_dir, "target_gene_presence_partial.csv"))

cat("METADATA_N_ROWS=", nrow(meta), "\n", sep = "")
cat("METADATA_N_COLS=", ncol(meta), "\n", sep = "")
cat("MALIGNANT_EPITHELIAL_N=", sum(meta$subtype == "Malignant epithelial cells"), "\n", sep = "")
cat("NK_N=", sum(meta$celltype == "NK cells"), "\n", sep = "")
cat("SAMPLES_WITH_BOTH=", sum(sample_special$malignant_epithelial_cells > 0 & sample_special$nk_cells > 0), "\n", sep = "")
cat("PRIMARY_PATIENTS_WITH_BOTH_GT0=", uniqueN(primary_both_positive$patient_id), "\n", sep = "")
cat("NK_CLASS=", paste(nk_class, collapse = ","), "\n", sep = "")
cat("NK_RAW_COUNTS_PRESENT=", raw_counts_present, "\n", sep = "")
cat("STAGE2_METADATA_NK_INSPECTION_COMPLETE\n")
