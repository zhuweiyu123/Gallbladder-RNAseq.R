#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
out <- file.path(project, "06_feasibility")

mapping <- fread(file.path(out, "patient_sample_mapping.csv"))
sites <- fread(file.path(out, "sample_site_mapping.csv"))
sample_only <- fread(file.path(out, "sample_info_without_metadata_cells.csv"))
subtypes <- fread(file.path(out, "subtype_counts.csv"))
primary <- fread(file.path(out, "primary_patient_cell_counts.csv"))
thresholds <- fread(file.path(out, "primary_patient_threshold_sensitivity.csv"))

patient_multiplicity <- mapping[, .(
  n_samples = uniqueN(sample_id),
  sample_ids = paste(sort(unique(sample_id)), collapse = " | "),
  sites = paste(sort(unique(standardized_site)), collapse = " | ")
), by = patient_id][order(-n_samples, patient_id)]
fwrite(patient_multiplicity, file.path(out, "patient_sample_multiplicity.csv"))

patient_all_sites <- mapping[, .(
  malignant_epithelial_cells = sum(malignant_epithelial_cells, na.rm = TRUE),
  nk_cells = sum(nk_cells, na.rm = TRUE),
  n_samples = uniqueN(sample_id),
  site_categories = paste(sort(unique(standardized_site)), collapse = " | ")
), by = patient_id]
fwrite(patient_all_sites, file.path(out, "patient_all_sites_cell_counts.csv"))

site_counts <- sites[, .(n_samples = .N), by = standardized_site][order(-n_samples)]
fwrite(site_counts, file.path(out, "site_category_counts.csv"))

sample_only_summary <- sample_only[, .(n_samples = .N), by = .(
  scRNA_seq = `scRNA-seq`,
  standardized_site
)][order(scRNA_seq, standardized_site)]
fwrite(sample_only_summary, file.path(out, "sample_info_without_metadata_summary.csv"))

nk_subtypes <- subtypes[celltype == "NK cells"][order(-n_cells)]
fwrite(nk_subtypes, file.path(out, "nk_subtype_counts.csv"))

metrics <- c(
  paste0("metadata_mapped_patients: ", uniqueN(mapping$patient_id)),
  paste0("patients_with_multiple_metadata_samples: ", patient_multiplicity[n_samples > 1, .N]),
  paste0("maximum_samples_per_patient: ", max(patient_multiplicity$n_samples)),
  paste0("metadata_samples: ", uniqueN(mapping$sample_id)),
  paste0("primary_tumor_samples: ", sites[standardized_site == "Primary tumor", .N]),
  paste0("all_site_patients_malignant_gt0_and_nk_gt0: ", patient_all_sites[malignant_epithelial_cells > 0 & nk_cells > 0, uniqueN(patient_id)]),
  paste0("primary_tumor_patients: ", primary[, uniqueN(patient_id)]),
  paste0("primary_tumor_patients_malignant_gt0_and_nk_gt0: ", primary[malignant_epithelial_cells > 0 & nk_cells > 0, uniqueN(patient_id)]),
  paste0("sample_info_only_rows: ", nrow(sample_only)),
  paste0("sample_info_only_scRNA_yes: ", sample_only[`scRNA-seq` == "Yes", .N]),
  paste0("sample_info_only_scRNA_no: ", sample_only[`scRNA-seq` == "No", .N]),
  paste0("site_counts: ", paste(site_counts$standardized_site, site_counts$n_samples, sep = "=", collapse = " | ")),
  paste0("nk_subtypes: ", paste(nk_subtypes$subtype, nk_subtypes$n_cells, sep = "=", collapse = " | ")),
  "sufficient_cell_threshold: 待核实（研究方案未预设阈值）"
)

writeLines(metrics, file.path(out, "stage2_metrics.txt"), useBytes = TRUE)
cat(paste(metrics, collapse = "\n"), "\n\n")
print(thresholds)
