#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
sample_file <- file.path(project, "02_processed_data", "All_single_cells", "sample_info.xlsx")
out <- file.path(project, "06_feasibility", "stage3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

stop_if_not <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

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

make_wide <- function(long, compartment) {
  id_columns <- if (compartment == "malignant") {
    c("sample_name", "malignant_cell_n", "malignant_pseudobulk_library_size")
  } else {
    c("sample_name", "NK_cell_n", "NK_pseudobulk_library_size")
  }
  result <- dcast(
    long,
    as.formula(paste(paste(id_columns, collapse = " + "), "~ gene")),
    value.var = c("RNA_detected_n", "RNA_detected_fraction", "pseudobulk_raw_count", "pseudobulk_CPM")
  )
  metric_prefix <- c(
    RNA_detected_n = "RNA_detected_n",
    RNA_detected_fraction = "RNA_detected_fraction",
    pseudobulk_raw_count = "pseudobulk_raw_count",
    pseudobulk_CPM = "pseudobulk_CPM"
  )
  for (old_prefix in names(metric_prefix)) {
    columns <- grep(paste0("^", old_prefix, "_"), names(result), value = TRUE)
    genes <- sub(paste0("^", old_prefix, "_"), "", columns)
    setnames(result, columns, paste0(genes, "_", metric_prefix[[old_prefix]]))
  }
  result
}

mal_long <- rbindlist(list(
  fread(file.path(out, "malignant_sample_gene_detection.csv")),
  fread(file.path(out, "supplementary_hlae_malignant_sample_detection.csv"))
), use.names = TRUE, fill = TRUE)
nk_long <- fread(file.path(out, "nk_sample_gene_detection.csv"))
stop_if_not(uniqueN(mal_long$sample_name) == 58L, "Expected 58 malignant samples")
stop_if_not(uniqueN(nk_long$sample_name) == 133L, "Expected 133 NK samples")

sample_info <- as.data.table(read_excel(sample_file, skip = 1, .name_repair = "unique"))
sample_info <- sample_info[!is.na(`Sample ID`) & trimws(as.character(`Sample ID`)) != ""]
sample_info[, `Sample ID` := as.character(`Sample ID`)]
sample_info[, patient_id := derive_patient_id(`Sample ID`)]
sample_info[, standardized_site := standardize_site(Site, `Histological type`)]

primary <- sample_info[`scRNA-seq` == "Yes" & standardized_site == "Primary tumor", .(
  sample_name = `Sample ID`,
  patient_id,
  patient_id_source = "derived_from_sample_prefix",
  site = as.character(Site),
  standardized_site,
  TNM_stage = as.character(`TNM stage`),
  histological_type = as.character(`Histological type`),
  metastatic_group = as.character(`Metastatic group`)
)]
stop_if_not(nrow(primary) == 86L, "Expected 86 scRNA Primary tumor samples")
stop_if_not(!anyDuplicated(primary$sample_name), "Duplicated Primary sample_name")
stop_if_not(!anyNA(primary$patient_id), "A Primary sample has no derived patient identifier")

mal_wide <- make_wide(mal_long, "malignant")
nk_wide <- make_wide(nk_long, "NK")
combined <- merge(primary, mal_wide, by = "sample_name", all.x = TRUE, sort = FALSE)
combined <- merge(combined, nk_wide, by = "sample_name", all.x = TRUE, sort = FALSE)
combined[, `:=`(
  has_malignant_compartment = !is.na(malignant_cell_n),
  has_NK_compartment = !is.na(NK_cell_n)
)]
combined[, exact_compartment_matched := has_malignant_compartment & has_NK_compartment]
combined[is.na(malignant_cell_n), malignant_cell_n := 0L]
combined[is.na(NK_cell_n), NK_cell_n := 0L]
combined[, structurally_evaluable_reference :=
  exact_compartment_matched & malignant_cell_n >= 100L & NK_cell_n >= 100L &
    malignant_pseudobulk_library_size > 0 & NK_pseudobulk_library_size > 0]
combined[, patient_level_inference_reference := fifelse(
  structurally_evaluable_reference,
  "member_of_reference_structural_set",
  "not_member_of_reference_structural_set"
)]
combined[, match_scope := fifelse(exact_compartment_matched, "exact_same_primary_sample", "primary_sample_not_exactly_matched")]

required_gene_columns <- c(
  "STK31_RNA_detected_n", "STK31_RNA_detected_fraction", "STK31_pseudobulk_raw_count",
  "HLA-B_RNA_detected_fraction", "HLA-B_pseudobulk_raw_count",
  paste0(c("HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2"), "_pseudobulk_raw_count"),
  "KIR3DL1_RNA_detected_n", "KIR3DL1_RNA_detected_fraction", "KIR3DL1_pseudobulk_raw_count",
  paste0(c("KIR2DL1", "KIR2DL3", "KLRC1"), "_RNA_detected_n"),
  paste0(c("KIR2DL1", "KIR2DL3", "KLRC1"), "_RNA_detected_fraction"),
  paste0(c("KIR2DL1", "KIR2DL3", "KLRC1"), "_pseudobulk_raw_count")
)
stop_if_not(all(required_gene_columns %in% names(combined)), "Required matched-table gene columns missing")

setorder(combined, patient_id, sample_name)
fwrite(combined, file.path(out, "primary_all_samples_compartment_availability.csv"))
fwrite(combined[exact_compartment_matched == TRUE], file.path(out, "primary_matched_stk31_hla_nk.csv"))

structural_thresholds <- CJ(
  malignant_cell_threshold = c(20L, 50L, 100L, 200L),
  NK_cell_threshold = c(20L, 50L, 100L, 200L)
)
base_exact_sample_n <- combined[exact_compartment_matched == TRUE, .N]
base_exact_patient_n <- combined[exact_compartment_matched == TRUE, uniqueN(patient_id)]
structural_thresholds <- rbindlist(lapply(seq_len(nrow(structural_thresholds)), function(index) {
  malignant_threshold <- structural_thresholds$malignant_cell_threshold[index]
  nk_threshold <- structural_thresholds$NK_cell_threshold[index]
  eligible <- combined[
    exact_compartment_matched &
      malignant_cell_n >= malignant_threshold &
      NK_cell_n >= nk_threshold &
      malignant_pseudobulk_library_size > 0 &
      NK_pseudobulk_library_size > 0
  ]
  data.table(
    malignant_cell_threshold = malignant_threshold,
    NK_cell_threshold = nk_threshold,
    exact_sample_n = nrow(eligible),
    unique_patient_n = uniqueN(eligible$patient_id),
    eligible_sample_n = nrow(eligible),
    eligible_patient_n = uniqueN(eligible$patient_id),
    base_exact_primary_sample_n = base_exact_sample_n,
    base_exact_primary_patient_n = base_exact_patient_n,
    patient_count_rule = "any_exact_primary_sample_passes",
    sample_names = paste(sort(eligible$sample_name), collapse = "|"),
    patient_ids = paste(sort(unique(eligible$patient_id)), collapse = "|")
  )
}), use.names = TRUE)
structural_thresholds[, criterion_type := "structural_only_no_target_gene_filter"]
setcolorder(structural_thresholds, c(
  "criterion_type", "malignant_cell_threshold", "NK_cell_threshold",
  "exact_sample_n", "unique_patient_n", "eligible_sample_n", "eligible_patient_n",
  "base_exact_primary_sample_n",
  "base_exact_primary_patient_n", "patient_count_rule", "sample_names", "patient_ids"
))
fwrite(structural_thresholds, file.path(out, "primary_structural_sensitivity_stage3.csv"))

signal_thresholds <- CJ(
  malignant_cell_threshold = c(20L, 50L, 100L, 200L),
  NK_cell_threshold = c(20L, 50L, 100L, 200L),
  STK31_detected_cell_threshold = c(1L, 3L, 5L, 10L),
  KIR3DL1_detected_cell_threshold = c(1L, 3L, 5L, 10L)
)
signal_thresholds <- rbindlist(lapply(seq_len(nrow(signal_thresholds)), function(index) {
  malignant_threshold <- signal_thresholds$malignant_cell_threshold[index]
  nk_threshold <- signal_thresholds$NK_cell_threshold[index]
  stk31_threshold <- signal_thresholds$STK31_detected_cell_threshold[index]
  kir3dl1_threshold <- signal_thresholds$KIR3DL1_detected_cell_threshold[index]
  eligible <- combined[
    exact_compartment_matched &
      malignant_cell_n >= malignant_threshold &
      NK_cell_n >= nk_threshold &
      STK31_RNA_detected_n >= stk31_threshold &
      KIR3DL1_RNA_detected_n >= kir3dl1_threshold &
      malignant_pseudobulk_library_size > 0 &
      NK_pseudobulk_library_size > 0
  ]
  data.table(
    malignant_cell_threshold = malignant_threshold,
    NK_cell_threshold = nk_threshold,
    STK31_detected_cell_threshold = stk31_threshold,
    KIR3DL1_detected_cell_threshold = kir3dl1_threshold,
    exact_sample_n = nrow(eligible),
    unique_patient_n = uniqueN(eligible$patient_id),
    eligible_sample_n = nrow(eligible),
    eligible_patient_n = uniqueN(eligible$patient_id),
    base_exact_primary_sample_n = base_exact_sample_n,
    base_exact_primary_patient_n = base_exact_patient_n,
    patient_count_rule = "any_exact_primary_sample_passes",
    sample_names = paste(sort(eligible$sample_name), collapse = "|"),
    patient_ids = paste(sort(unique(eligible$patient_id)), collapse = "|")
  )
}), use.names = TRUE)
signal_thresholds[, `:=`(
  criterion_type = "signal_supported_descriptive_not_an_exclusion_rule",
  reference_structural_cell_thresholds = malignant_cell_threshold == 100L & NK_cell_threshold == 100L
)]
setcolorder(signal_thresholds, c(
  "criterion_type", "reference_structural_cell_thresholds",
  "malignant_cell_threshold", "NK_cell_threshold",
  "STK31_detected_cell_threshold", "KIR3DL1_detected_cell_threshold",
  "exact_sample_n", "unique_patient_n", "eligible_sample_n", "eligible_patient_n",
  "base_exact_primary_sample_n",
  "base_exact_primary_patient_n", "patient_count_rule", "sample_names", "patient_ids"
))
stop_if_not(nrow(signal_thresholds) == 256L, "Threshold grid is not 4x4x4x4")
fwrite(signal_thresholds, file.path(out, "primary_threshold_sensitivity_stage3.csv"))

reference <- combined[structurally_evaluable_reference == TRUE]
inference_label <- if (uniqueN(reference$patient_id) < 10L) {
  "not_suitable_for_patient_level_association"
} else if (uniqueN(reference$patient_id) < 20L) {
  "pilot_or_descriptive_only"
} else {
  "may_support_exploratory_association_pending_power_analysis"
}

coverage_genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2", "HLA-E", "KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRC1")
coverage <- rbindlist(lapply(coverage_genes, function(gene) {
  count_column <- paste0(gene, "_pseudobulk_raw_count")
  if (!count_column %in% names(reference)) return(NULL)
  values <- reference[[count_column]]
  fraction <- if (length(values)) mean(values > 0) else NA_real_
  label <- fifelse(
    is.na(fraction), "not_evaluable",
    fifelse(fraction >= 0.8, "broadly_detected", fifelse(fraction >= 0.5, "intermittently_detected", "sparse"))
  )
  ordered <- sort(values, decreasing = TRUE)
  total <- sum(ordered)
  data.table(
    reference_set = "exact_primary_malignant_ge100_NK_ge100",
    gene = gene,
    structurally_evaluable_sample_n = nrow(reference),
    structurally_evaluable_unique_patient_n = uniqueN(reference$patient_id),
    pseudobulk_nonzero_sample_n = sum(values > 0),
    pseudobulk_nonzero_fraction = fraction,
    technical_coverage_label = label,
    top1_share_of_total_counts = if (total > 0) ordered[1] / total else NA_real_,
    top3_share_of_total_counts = if (total > 0) sum(head(ordered, 3L)) / total else NA_real_,
    highly_concentrated_top1_gt50pct = if (total > 0) ordered[1] / total > 0.5 else NA
  )
}))
fwrite(coverage, file.path(out, "primary_reference_signal_coverage.csv"))

summary <- data.table(
  metric = c(
    "primary_sample_n", "primary_unique_patient_n", "exact_compartment_matched_sample_n",
    "exact_compartment_matched_unique_patient_n", "reference_structurally_evaluable_sample_n",
    "reference_structurally_evaluable_unique_patient_n", "reference_inference_label"
  ),
  value = c(
    nrow(combined), uniqueN(combined$patient_id),
    combined[exact_compartment_matched == TRUE, .N],
    combined[exact_compartment_matched == TRUE, uniqueN(patient_id)],
    nrow(reference), uniqueN(reference$patient_id), inference_label
  )
)
fwrite(summary, file.path(out, "primary_matching_summary.csv"))
cat(paste(summary$metric, summary$value, sep = "=", collapse = "\n"), "\n")
