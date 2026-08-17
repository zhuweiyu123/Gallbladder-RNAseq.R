#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(SeuratObject)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
nk_file <- file.path(project, "01_raw_data", "NK.RDS")
sample_file <- file.path(project, "02_processed_data", "All_single_cells", "sample_info.xlsx")
out <- file.path(project, "06_feasibility", "stage3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
safe_fraction <- function(numerator, denominator) fifelse(denominator > 0, numerator / denominator, NA_real_)
safe_cpm <- function(count, library_size) fifelse(library_size > 0, count / library_size * 1e6, NA_real_)

summarize_concentration <- function(data, compartment_name, unit_column, scope_name, metrics) {
  rbindlist(lapply(metrics, function(metric_name) {
    data[, {
      signal <- as.numeric(get(metric_name))
      unit_id <- as.character(get(unit_column))
      ordering <- order(-signal, unit_id)
      ordered_signal <- signal[ordering]
      ordered_unit <- unit_id[ordering]
      total_signal <- sum(ordered_signal)
      cumulative <- if (total_signal > 0) cumsum(ordered_signal) / total_signal else rep(NA_real_, length(ordered_signal))
      shares <- if (total_signal > 0) ordered_signal / total_signal else rep(NA_real_, length(ordered_signal))
      hhi <- if (total_signal > 0) sum(shares^2) else NA_real_
      .(
        compartment = compartment_name,
        scope = scope_name,
        unit_type = unit_column,
        signal_metric = metric_name,
        total_unit_n = .N,
        nonzero_unit_n = sum(signal > 0),
        nonzero_unit_fraction = mean(signal > 0),
        median_nonzero_signal = if (any(signal > 0)) median(signal[signal > 0]) else NA_real_,
        total_signal = total_signal,
        top1_share = if (total_signal > 0) sum(head(ordered_signal, 1L)) / total_signal else NA_real_,
        top3_share = if (total_signal > 0) sum(head(ordered_signal, 3L)) / total_signal else NA_real_,
        top5_share = if (total_signal > 0) sum(head(ordered_signal, 5L)) / total_signal else NA_real_,
        top10_share = if (total_signal > 0) sum(head(ordered_signal, 10L)) / total_signal else NA_real_,
        n_units_for_50pct = if (total_signal > 0) which(cumulative >= 0.5)[1] else NA_integer_,
        n_units_for_80pct = if (total_signal > 0) which(cumulative >= 0.8)[1] else NA_integer_,
        hhi = hhi,
        effective_unit_n = if (!is.na(hhi) && hhi > 0) 1 / hhi else NA_real_,
        max_unit_id = if (length(ordered_unit)) ordered_unit[1] else NA_character_,
        max_unit_signal = if (length(ordered_signal)) ordered_signal[1] else NA_real_,
        highly_concentrated_top1_gt50pct = if (total_signal > 0) ordered_signal[1] / total_signal > 0.5 else NA
      )
    }, by = gene]
  }), use.names = TRUE)
}

derive_patient_id <- function(x) {
  answer <- sub("^([[:alpha:]]+_[0-9]+).*$", "\\1", x)
  answer[!grepl("^[[:alpha:]]+_[0-9]+", x)] <- NA_character_
  answer
}

standardize_site <- function(site, histology) {
  site <- trimws(as.character(site))
  histology <- trimws(as.character(histology))
  answer <- rep("Unknown", length(site))
  carcinoma <- grepl("carcinoma", histology, ignore.case = TRUE)
  answer[site == "Primary" & carcinoma] <- "Primary tumor"
  answer[site == "Primary" & !carcinoma] <- "Other"
  answer[site == "Liver invasion"] <- "Liver invasion"
  answer[site == "Liver metastatic lesion"] <- "Liver metastasis"
  answer[site == "Lymphonodus metastatic lesion"] <- "Lymph node metastasis"
  answer[site == "Omentum"] <- "Omental metastasis"
  answer[site %in% c("PBMC", "Primary polyp")] <- "Other"
  answer[grepl("adjacent|NAT", site, ignore.case = TRUE)] <- "Adjacent/NAT"
  answer
}

wilson_interval <- function(successes, total, confidence = 0.95) {
  z <- qnorm(1 - (1 - confidence) / 2)
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  list(lower = pmax(0, center - half), upper = pmin(1, center + half))
}

target_genes <- c(
  "KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRC1", "KLRD1", "NKG7", "GNLY",
  "PRF1", "GZMB", "FGFBP2", "FCGR3A", "XCL1", "XCL2"
)

message("Reading NK.RDS")
nk <- readRDS(nk_file)
stop_if_not(inherits(nk, "Seurat"), "NK.RDS is not a Seurat object")
stop_if_not(DefaultAssay(nk) == "RNA", "NK DefaultAssay is not RNA")
counts <- LayerData(nk[["RNA"]], layer = "counts")
stop_if_not(inherits(counts, "sparseMatrix"), "NK raw counts are not sparse")
stop_if_not(!"KIR2DL2" %in% rownames(counts), "KIR2DL2 unexpectedly present; scope must be reviewed")
feature_status <- data.table(
  compartment = "NK",
  gene = c(target_genes, "KIR2DL2"),
  feature_present = c(target_genes, "KIR2DL2") %in% rownames(counts),
  analysis_status = fifelse(
    c(target_genes, "KIR2DL2") %in% rownames(counts),
    "feature_present_analyzed",
    "feature_absent_not_zero_filled"
  )
)
fwrite(feature_status, file.path(out, "nk_target_feature_status.csv"))
stop_if_not(all(feature_status[gene %in% target_genes, feature_present]), "An NK target gene is missing from the raw-count feature list")

nk_meta <- as.data.table(nk[[]], keep.rownames = "barcode")
stop_if_not(identical(colnames(counts), nk_meta$barcode), "NK counts columns differ from metadata row order")
stop_if_not(all(c("orig.ident", "celltype") %in% names(nk_meta)), "Required NK metadata fields missing")
stop_if_not(!anyNA(nk_meta$orig.ident), "Missing NK orig.ident")
stop_if_not(!anyNA(nk_meta$celltype), "Missing author NK subtype")
stop_if_not(uniqueN(nk_meta$orig.ident) == 133L, "Expected 133 NK-containing samples")
stop_if_not(uniqueN(nk_meta$celltype) == 14L, "Expected 14 author NK subtypes")

sample_info <- as.data.table(read_excel(sample_file, skip = 1, .name_repair = "unique"))
sample_info <- sample_info[!is.na(`Sample ID`) & trimws(as.character(`Sample ID`)) != ""]
sample_info[, `Sample ID` := as.character(`Sample ID`)]
sample_info[, patient_id := derive_patient_id(`Sample ID`)]
sample_info[, standardized_site := standardize_site(Site, `Histological type`)]
sample_index <- match(nk_meta$orig.ident, sample_info$`Sample ID`)
stop_if_not(!anyNA(sample_index), "An NK orig.ident does not map to sample_info")

cell_library <- as.numeric(Matrix::colSums(counts))
zero_library_cell_n <- sum(cell_library <= 0)
stop_if_not(zero_library_cell_n == 0L, "An NK cell has zero all-gene raw library size; structural sample validity must be reviewed")
cell_base <- data.table(
  barcode = nk_meta$barcode,
  sample_name = as.character(nk_meta$orig.ident),
  patient_id = sample_info$patient_id[sample_index],
  patient_id_source = "derived_from_sample_prefix",
  site = as.character(sample_info$Site[sample_index]),
  standardized_site = sample_info$standardized_site[sample_index],
  TNM_stage = as.character(sample_info$`TNM stage`[sample_index]),
  histological_type = as.character(sample_info$`Histological type`[sample_index]),
  metastatic_group = as.character(sample_info$`Metastatic group`[sample_index]),
  nk_subtype = as.character(nk_meta$celltype),
  all_gene_raw_umi = cell_library
)

sample_library <- cell_base[, .(
  NK_cell_n = .N,
  NK_pseudobulk_library_size = sum(all_gene_raw_umi)
), by = .(
  sample_name, patient_id, patient_id_source, site, standardized_site,
  TNM_stage, histological_type, metastatic_group
)]
stop_if_not(nrow(sample_library) == 133L, "Expected 133 NK sample summaries")

overall_list <- vector("list", length(target_genes))
sample_list <- vector("list", length(target_genes))
subtype_list <- vector("list", length(target_genes))
for (gene_index in seq_along(target_genes)) {
  gene <- target_genes[gene_index]
  values <- as.numeric(counts[gene, ])
  detected <- values > 0
  detected_values <- values[detected]
  overall_list[[gene_index]] <- data.table(
    gene = gene,
    feature_present = TRUE,
    total_NK_n = length(values),
    RNA_detected_n = sum(detected),
    RNA_detected_fraction = mean(detected),
    total_raw_UMI = sum(values),
    mean_raw_count_among_detected = if (length(detected_values)) mean(detected_values) else NA_real_,
    median_raw_count_among_detected = if (length(detected_values)) median(detected_values) else NA_real_
  )

  working <- cell_base[, .(sample_name, nk_subtype)]
  working[, `:=`(raw_count = values, RNA_detected = detected)]
  gene_sample <- working[, .(
    RNA_detected_n = sum(RNA_detected),
    pseudobulk_raw_count = sum(raw_count)
  ), by = sample_name]
  gene_sample <- merge(sample_library, gene_sample, by = "sample_name", all.x = TRUE, sort = FALSE)
  gene_sample[, `:=`(
    gene = gene,
    RNA_detected_fraction = safe_fraction(RNA_detected_n, NK_cell_n),
    pseudobulk_CPM = safe_cpm(pseudobulk_raw_count, NK_pseudobulk_library_size)
  )]
  setcolorder(gene_sample, c(
    "sample_name", "patient_id", "patient_id_source", "site", "standardized_site",
    "TNM_stage", "histological_type", "metastatic_group", "gene", "NK_cell_n",
    "RNA_detected_n", "RNA_detected_fraction", "pseudobulk_raw_count",
    "NK_pseudobulk_library_size", "pseudobulk_CPM"
  ))
  sample_list[[gene_index]] <- gene_sample

  gene_subtype <- working[, .(
    subtype_cell_n = .N,
    RNA_detected_n = sum(RNA_detected),
    pseudobulk_raw_count = sum(raw_count)
  ), by = nk_subtype]
  gene_subtype[, `:=`(
    gene = gene,
    RNA_detected_fraction = safe_fraction(RNA_detected_n, subtype_cell_n)
  )]
  interval <- wilson_interval(gene_subtype$RNA_detected_n, gene_subtype$subtype_cell_n)
  gene_subtype[, `:=`(
    wilson95_lower = interval$lower,
    wilson95_upper = interval$upper
  )]
  gene_subtype[, observed_detection_fraction_rank := frank(-RNA_detected_fraction, ties.method = "min")]
  gene_subtype[, observed_detected_n_rank := frank(-RNA_detected_n, ties.method = "min")]
  gene_subtype[, fraction_of_all_gene_RNA_detected_cells := {
    total_detected <- sum(RNA_detected_n)
    if (total_detected > 0) RNA_detected_n / total_detected else rep(NA_real_, .N)
  }]
  setcolorder(gene_subtype, c(
    "gene", "nk_subtype", "subtype_cell_n", "RNA_detected_n", "RNA_detected_fraction",
    "wilson95_lower", "wilson95_upper", "observed_detection_fraction_rank",
    "observed_detected_n_rank", "fraction_of_all_gene_RNA_detected_cells", "pseudobulk_raw_count"
  ))
  subtype_list[[gene_index]] <- gene_subtype
}

overall <- rbindlist(overall_list)
sample_long <- rbindlist(sample_list)
subtype_long <- rbindlist(subtype_list)
setorder(subtype_long, gene, observed_detection_fraction_rank, -subtype_cell_n, nk_subtype)
fwrite(overall, file.path(out, "nk_target_gene_detection.csv"))
fwrite(subtype_long, file.path(out, "nk_subtype_gene_detection.csv"))
fwrite(sample_long, file.path(out, "nk_sample_gene_detection.csv"))

kir3dl1 <- sample_long[gene == "KIR3DL1"]
kir3dl1[, is_primary_tumor := standardized_site == "Primary tumor"]
setorder(kir3dl1, standardized_site, -RNA_detected_n, sample_name)
fwrite(kir3dl1, file.path(out, "kir3dl1_sample_detection.csv"))

kir_thresholds <- rbindlist(lapply(c(1L, 3L, 5L, 10L, 20L), function(threshold) {
  rbindlist(list(
    data.table(
      scope = "all_NK_samples", threshold_RNA_detected_cell_n = threshold,
      sample_n = kir3dl1[RNA_detected_n >= threshold, .N],
      unique_patient_n = kir3dl1[RNA_detected_n >= threshold, uniqueN(patient_id)],
      patient_aggregated_threshold_n = kir3dl1[, .(patient_RNA_detected_n = sum(RNA_detected_n)), by = patient_id][
        patient_RNA_detected_n >= threshold, .N
      ],
      total_sample_n = nrow(kir3dl1), total_unique_patient_n = uniqueN(kir3dl1$patient_id)
    ),
    data.table(
      scope = "primary_tumor_NK_samples", threshold_RNA_detected_cell_n = threshold,
      sample_n = kir3dl1[standardized_site == "Primary tumor" & RNA_detected_n >= threshold, .N],
      unique_patient_n = kir3dl1[standardized_site == "Primary tumor" & RNA_detected_n >= threshold, uniqueN(patient_id)],
      patient_aggregated_threshold_n = kir3dl1[standardized_site == "Primary tumor", .(
        patient_RNA_detected_n = sum(RNA_detected_n)
      ), by = patient_id][patient_RNA_detected_n >= threshold, .N],
      total_sample_n = kir3dl1[standardized_site == "Primary tumor", .N],
      total_unique_patient_n = kir3dl1[standardized_site == "Primary tumor", uniqueN(patient_id)]
    )
  ))
}))
kir_thresholds[, sample_fraction := safe_fraction(sample_n, total_sample_n)]
fwrite(kir_thresholds, file.path(out, "kir3dl1_sample_threshold_summary.csv"))

concentration <- rbindlist(list(
  summarize_concentration(
    sample_long, "NK", "sample_name", "all_NK_samples",
    c("RNA_detected_n", "pseudobulk_raw_count")
  ),
  summarize_concentration(
    sample_long[standardized_site == "Primary tumor"], "NK", "sample_name",
    "primary_tumor_NK_samples", c("RNA_detected_n", "pseudobulk_raw_count")
  )
), use.names = TRUE)
setcolorder(concentration, c(
  "compartment", "gene", "scope", "unit_type", "signal_metric", "total_unit_n",
  "nonzero_unit_n", "nonzero_unit_fraction", "median_nonzero_signal", "total_signal",
  "top1_share", "top3_share", "top5_share", "top10_share", "n_units_for_50pct",
  "n_units_for_80pct", "hhi", "effective_unit_n", "max_unit_id", "max_unit_signal",
  "highly_concentrated_top1_gt50pct"
))
fwrite(concentration, file.path(out, "nk_signal_concentration.csv"))

subtype_concentration <- summarize_concentration(
  subtype_long, "NK", "nk_subtype", "all_author_NK_subtypes",
  c("RNA_detected_n", "pseudobulk_raw_count")
)
setcolorder(subtype_concentration, c(
  "compartment", "gene", "scope", "unit_type", "signal_metric", "total_unit_n",
  "nonzero_unit_n", "nonzero_unit_fraction", "median_nonzero_signal", "total_signal",
  "top1_share", "top3_share", "top5_share", "top10_share", "n_units_for_50pct",
  "n_units_for_80pct", "hhi", "effective_unit_n", "max_unit_id", "max_unit_signal",
  "highly_concentrated_top1_gt50pct"
))
fwrite(subtype_concentration, file.path(out, "nk_subtype_signal_concentration.csv"))

wide_measures <- dcast(
  sample_long,
  sample_name + patient_id + patient_id_source + site + standardized_site +
    TNM_stage + histological_type + metastatic_group + NK_cell_n +
    NK_pseudobulk_library_size ~ gene,
  value.var = c("RNA_detected_n", "RNA_detected_fraction", "pseudobulk_raw_count", "pseudobulk_CPM")
)
fwrite(wide_measures, file.path(out, "nk_sample_gene_detection_wide.csv"))

audit <- c(
  paste0("NK_matrix_dim=", paste(dim(counts), collapse = "x")),
  paste0("NK_matrix_nnz=", length(counts@x)),
  paste0("NK_cells=", ncol(counts)),
  paste0("NK_samples=", uniqueN(cell_base$sample_name)),
  paste0("NK_subtypes=", uniqueN(cell_base$nk_subtype)),
  paste0("all_cells_have_positive_library=", all(cell_library > 0)),
  paste0("zero_library_cell_n=", zero_library_cell_n),
  paste0("target_gene_sums=", paste(overall$gene, overall$total_raw_UMI, sep = ":", collapse = "|"))
)
writeLines(audit, file.path(out, "stage3_nk_detection_audit.txt"), useBytes = TRUE)
cat(paste(audit, collapse = "\n"), "\n")
