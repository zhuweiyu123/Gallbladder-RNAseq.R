#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
mal_dir <- file.path(project, "02_processed_data", "malignant_epithelial")
target_dir <- file.path(mal_dir, "target_counts")
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

message("Reading target-only malignant Matrix Market subset")
feature_manifest <- fread(file.path(out, "malignant_target_feature_manifest.csv"))
barcodes <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(target_dir, "barcodes.tsv.gz"))),
  header = FALSE, col.names = "barcode", showProgress = FALSE
)
target_counts <- readMM(gzfile(file.path(target_dir, "matrix.mtx.gz")))
target_counts <- as(target_counts, "CsparseMatrix")
rownames(target_counts) <- feature_manifest$gene
colnames(target_counts) <- barcodes$barcode

mal_meta <- fread(file.path(mal_dir, "malignant_cell_metadata.csv.gz"), showProgress = FALSE)
cell_library <- fread(file.path(mal_dir, "malignant_cell_library_sizes.tsv"))
stop_if_not(identical(dim(target_counts), c(nrow(feature_manifest), 96942L)), "Unexpected malignant target matrix dimension")
stop_if_not(identical(colnames(target_counts), mal_meta$barcode), "Malignant target columns differ from cell metadata")
stop_if_not(identical(cell_library$malignant_col, seq_len(nrow(mal_meta))), "Cell library index differs from manifest")
zero_library_cell_n <- sum(cell_library$all_gene_raw_umi <= 0)
stop_if_not(zero_library_cell_n == 0L, "A malignant cell has zero all-gene raw library size; structural sample validity must be reviewed")
stop_if_not(!anyDuplicated(rownames(target_counts)), "Duplicated target genes")
stop_if_not(!anyDuplicated(colnames(target_counts)), "Duplicated malignant barcodes")

saveRDS(target_counts, file.path(target_dir, "malignant_target_raw_counts.rds"), compress = "xz")

cell_base <- mal_meta[, .(
  malignant_col, barcode, sample_name, patient_id, patient_id_source,
  site, standardized_site, TNM_stage, histological_type, metastatic_group
)]
cell_base[, all_gene_raw_umi := cell_library$all_gene_raw_umi]

sample_library <- cell_base[, .(
  malignant_cell_n = .N,
  malignant_pseudobulk_library_size = sum(all_gene_raw_umi)
), by = .(
  sample_name, patient_id, patient_id_source, site, standardized_site,
  TNM_stage, histological_type, metastatic_group
)]
stop_if_not(nrow(sample_library) == 58L, "Expected 58 malignant-containing samples")

overall_list <- vector("list", nrow(feature_manifest))
sample_list <- vector("list", nrow(feature_manifest))
for (gene_index in seq_len(nrow(feature_manifest))) {
  gene <- feature_manifest$gene[gene_index]
  values <- as.numeric(target_counts[gene_index, ])
  positive <- values > 0
  positive_values <- values[positive]
  overall_list[[gene_index]] <- data.table(
    gene = gene,
    panel = feature_manifest$panel[gene_index],
    feature_present = TRUE,
    total_cell_n = length(values),
    RNA_detected_n = sum(positive),
    RNA_detected_fraction = mean(positive),
    total_raw_UMI = sum(values),
    mean_raw_count_among_detected = if (length(positive_values)) mean(positive_values) else NA_real_,
    median_raw_count_among_detected = if (length(positive_values)) median(positive_values) else NA_real_
  )
  working <- cell_base[, .(sample_name)]
  working[, `:=`(raw_count = values, RNA_detected = positive)]
  gene_sample <- working[, .(
    RNA_detected_n = sum(RNA_detected),
    pseudobulk_raw_count = sum(raw_count)
  ), by = sample_name]
  gene_sample <- merge(sample_library, gene_sample, by = "sample_name", all.x = TRUE, sort = FALSE)
  gene_sample[, `:=`(
    gene = gene,
    panel = feature_manifest$panel[gene_index],
    RNA_detected_fraction = safe_fraction(RNA_detected_n, malignant_cell_n),
    pseudobulk_CPM = safe_cpm(pseudobulk_raw_count, malignant_pseudobulk_library_size)
  )]
  setcolorder(gene_sample, c(
    "sample_name", "patient_id", "patient_id_source", "site", "standardized_site",
    "TNM_stage", "histological_type", "metastatic_group", "gene", "panel",
    "malignant_cell_n", "RNA_detected_n", "RNA_detected_fraction",
    "pseudobulk_raw_count", "malignant_pseudobulk_library_size", "pseudobulk_CPM"
  ))
  sample_list[[gene_index]] <- gene_sample
}

overall_all <- rbindlist(overall_list)
sample_long_all <- rbindlist(sample_list)
overall <- overall_all[panel == "prespecified"]
sample_long <- sample_long_all[panel == "prespecified"]
fwrite(overall, file.path(out, "malignant_target_gene_detection.csv"))
fwrite(sample_long, file.path(out, "malignant_sample_gene_detection.csv"))
fwrite(overall_all[panel != "prespecified"], file.path(out, "supplementary_hlae_malignant_detection.csv"))
fwrite(sample_long_all[panel != "prespecified"], file.path(out, "supplementary_hlae_malignant_sample_detection.csv"))

stk31 <- sample_long[gene == "STK31"]
stk31[, is_primary_tumor := standardized_site == "Primary tumor"]
setorder(stk31, standardized_site, -RNA_detected_n, sample_name)
fwrite(stk31, file.path(out, "stk31_sample_detection.csv"))

stk31_thresholds <- rbindlist(lapply(c(1L, 3L, 5L, 10L, 20L), function(threshold) {
  rbindlist(list(
    data.table(
      scope = "all_malignant_samples", threshold_RNA_detected_cell_n = threshold,
      sample_n = stk31[RNA_detected_n >= threshold, .N],
      unique_patient_n = stk31[RNA_detected_n >= threshold, uniqueN(patient_id)],
      total_sample_n = nrow(stk31), total_unique_patient_n = uniqueN(stk31$patient_id)
    ),
    data.table(
      scope = "primary_tumor_malignant_samples", threshold_RNA_detected_cell_n = threshold,
      sample_n = stk31[standardized_site == "Primary tumor" & RNA_detected_n >= threshold, .N],
      unique_patient_n = stk31[standardized_site == "Primary tumor" & RNA_detected_n >= threshold, uniqueN(patient_id)],
      total_sample_n = stk31[standardized_site == "Primary tumor", .N],
      total_unique_patient_n = stk31[standardized_site == "Primary tumor", uniqueN(patient_id)]
    )
  ))
}))
stk31_thresholds[, sample_fraction := safe_fraction(sample_n, total_sample_n)]
fwrite(stk31_thresholds, file.path(out, "stk31_sample_threshold_summary.csv"))

concentration <- rbindlist(list(
  summarize_concentration(
    sample_long_all, "malignant_epithelial", "sample_name", "all_malignant_samples",
    c("RNA_detected_n", "pseudobulk_raw_count")
  ),
  summarize_concentration(
    sample_long_all[standardized_site == "Primary tumor"], "malignant_epithelial", "sample_name",
    "primary_tumor_malignant_samples", c("RNA_detected_n", "pseudobulk_raw_count")
  )
), use.names = TRUE)
setcolorder(concentration, c(
  "compartment", "gene", "scope", "unit_type", "signal_metric", "total_unit_n",
  "nonzero_unit_n", "nonzero_unit_fraction", "median_nonzero_signal", "total_signal",
  "top1_share", "top3_share", "top5_share", "top10_share", "n_units_for_50pct",
  "n_units_for_80pct", "hhi", "effective_unit_n", "max_unit_id", "max_unit_signal",
  "highly_concentrated_top1_gt50pct"
))
fwrite(concentration, file.path(out, "malignant_signal_concentration.csv"))

audit <- c(
  paste0("target_matrix_dim=", paste(dim(target_counts), collapse = "x")),
  paste0("target_matrix_nnz=", length(target_counts@x)),
  paste0("target_matrix_total_raw_UMI=", sum(target_counts)),
  paste0("malignant_cells=", nrow(mal_meta)),
  paste0("malignant_samples=", uniqueN(mal_meta$sample_name)),
  paste0("all_cells_have_positive_library=", all(cell_library$all_gene_raw_umi > 0)),
  paste0("zero_library_cell_n=", zero_library_cell_n),
  paste0("target_gene_sums=", paste(overall_all$gene, overall_all$total_raw_UMI, sep = ":", collapse = "|"))
)
writeLines(audit, file.path(out, "stage3_malignant_detection_audit.txt"), useBytes = TRUE)
cat(paste(audit, collapse = "\n"), "\n")
