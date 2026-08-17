#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
out <- file.path(project, "06_feasibility", "stage3")
mal_dir <- file.path(project, "02_processed_data", "malignant_epithelial")

checks <- list()
add_check <- function(name, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, passed = isTRUE(passed), detail = as.character(detail)
  )
}
assert_check <- function(name, passed, detail = "") {
  add_check(name, passed, detail)
  if (!isTRUE(passed)) stop(sprintf("Validation failed: %s (%s)", name, detail), call. = FALSE)
}
near_equal <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && all(abs(as.numeric(x) - as.numeric(y)) <= tolerance, na.rm = TRUE) &&
    identical(is.na(x), is.na(y))
}

required_rows <- c(
  malignant_target_gene_detection.csv = 8L,
  malignant_sample_gene_detection.csv = 464L,
  stk31_sample_detection.csv = 58L,
  nk_target_gene_detection.csv = 13L,
  nk_subtype_gene_detection.csv = 182L,
  nk_sample_gene_detection.csv = 1729L,
  kir3dl1_sample_detection.csv = 133L,
  primary_threshold_sensitivity_stage3.csv = 256L
)
for (file_name in names(required_rows)) {
  path <- file.path(out, file_name)
  assert_check(paste0("file_exists_", file_name), file.exists(path) && file.info(path)$size > 0, path)
  observed <- nrow(fread(path))
  assert_check(paste0("row_count_", file_name), observed == required_rows[[file_name]],
               sprintf("observed=%d expected=%d", observed, required_rows[[file_name]]))
}

mal_genes <- c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2")
nk_genes <- c("KIR3DL1", "KIR2DL1", "KIR2DL3", "KLRC1", "KLRD1", "NKG7", "GNLY", "PRF1", "GZMB", "FGFBP2", "FCGR3A", "XCL1", "XCL2")
mal_overall <- fread(file.path(out, "malignant_target_gene_detection.csv"))
mal_sample <- fread(file.path(out, "malignant_sample_gene_detection.csv"))
nk_overall <- fread(file.path(out, "nk_target_gene_detection.csv"))
nk_sample <- fread(file.path(out, "nk_sample_gene_detection.csv"))
nk_subtype <- fread(file.path(out, "nk_subtype_gene_detection.csv"))
nk_feature_status <- fread(file.path(out, "nk_target_feature_status.csv"))

assert_check("malignant_gene_set", setequal(mal_overall$gene, mal_genes), paste(sort(mal_overall$gene), collapse = "|"))
assert_check("NK_gene_set", setequal(nk_overall$gene, nk_genes), paste(sort(nk_overall$gene), collapse = "|"))
assert_check("KIR2DL2_not_zero_filled", !"KIR2DL2" %in% c(mal_overall$gene, nk_overall$gene, mal_sample$gene, nk_sample$gene), "feature absent")
assert_check("KIR2DL2_feature_absent_recorded", nrow(nk_feature_status[gene == "KIR2DL2" & !feature_present &
               analysis_status == "feature_absent_not_zero_filled"]) == 1L)
assert_check("all_reported_features_present", all(mal_overall$feature_present) && all(nk_overall$feature_present), "TRUE required")
assert_check("malignant_sample_gene_cartesian_complete", uniqueN(mal_sample[, .(sample_name, gene)]) == 58L * 8L &&
               all(mal_sample[, .N, by = sample_name]$N == 8L), "58x8")
assert_check("NK_sample_gene_cartesian_complete", uniqueN(nk_sample[, .(sample_name, gene)]) == 133L * 13L &&
               all(nk_sample[, .N, by = sample_name]$N == 13L), "133x13")
assert_check("NK_subtype_gene_cartesian_complete", uniqueN(nk_subtype[, .(nk_subtype, gene)]) == 14L * 13L &&
               all(nk_subtype[, .N, by = nk_subtype]$N == 13L), "14x13")
assert_check("patient_id_source_malignant", identical(unique(mal_sample$patient_id_source), "derived_from_sample_prefix"))
assert_check("patient_id_source_NK", identical(unique(nk_sample$patient_id_source), "derived_from_sample_prefix"))
stk_threshold <- fread(file.path(out, "stk31_sample_threshold_summary.csv"))
kir_threshold <- fread(file.path(out, "kir3dl1_sample_threshold_summary.csv"))
assert_check("standalone_threshold_values", setequal(stk_threshold$threshold_RNA_detected_cell_n, c(1L,3L,5L,10L,20L)) &&
               setequal(kir_threshold$threshold_RNA_detected_cell_n, c(1L,3L,5L,10L,20L)))
assert_check("standalone_threshold_patient_count_semantics", all(stk_threshold$unique_patient_n <= stk_threshold$total_unique_patient_n) &&
               all(kir_threshold$unique_patient_n <= kir_threshold$sample_n) &&
               all(kir_threshold$patient_aggregated_threshold_n <= kir_threshold$total_unique_patient_n),
             "unique_patient_n counts passing samples; patient_aggregated_threshold_n sums samples first")

mal_by_gene <- mal_sample[, .(
  sample_raw = sum(as.numeric(pseudobulk_raw_count)),
  sample_detected = sum(as.numeric(RNA_detected_n))
), by = gene]
mal_compare <- merge(mal_overall, mal_by_gene, by = "gene")
assert_check("malignant_overall_raw_equals_sample_sum", all(mal_compare$total_raw_UMI == mal_compare$sample_raw))
assert_check("malignant_overall_detected_equals_sample_sum", all(mal_compare$RNA_detected_n == mal_compare$sample_detected))
assert_check("malignant_overall_fraction", near_equal(mal_overall$RNA_detected_fraction,
                                                       mal_overall$RNA_detected_n / mal_overall$total_cell_n))
assert_check("malignant_sample_fraction", near_equal(mal_sample$RNA_detected_fraction,
                                                      mal_sample$RNA_detected_n / mal_sample$malignant_cell_n))

nk_by_sample <- nk_sample[, .(
  sample_raw = sum(as.numeric(pseudobulk_raw_count)),
  sample_detected = sum(as.numeric(RNA_detected_n))
), by = gene]
nk_by_subtype <- nk_subtype[, .(
  subtype_raw = sum(as.numeric(pseudobulk_raw_count)),
  subtype_detected = sum(as.numeric(RNA_detected_n))
), by = gene]
nk_compare <- Reduce(function(x, y) merge(x, y, by = "gene"), list(nk_overall, nk_by_sample, nk_by_subtype))
assert_check("NK_overall_raw_equals_sample_sum", all(nk_compare$total_raw_UMI == nk_compare$sample_raw))
assert_check("NK_overall_detected_equals_sample_sum", all(nk_compare$RNA_detected_n == nk_compare$sample_detected))
assert_check("NK_overall_raw_equals_subtype_sum", all(nk_compare$total_raw_UMI == nk_compare$subtype_raw))
assert_check("NK_overall_detected_equals_subtype_sum", all(nk_compare$RNA_detected_n == nk_compare$subtype_detected))
assert_check("NK_overall_fraction", near_equal(nk_overall$RNA_detected_fraction,
                                                nk_overall$RNA_detected_n / nk_overall$total_NK_n))
assert_check("NK_sample_fraction", near_equal(nk_sample$RNA_detected_fraction,
                                               nk_sample$RNA_detected_n / nk_sample$NK_cell_n))
assert_check("NK_subtype_fraction", near_equal(nk_subtype$RNA_detected_fraction,
                                                nk_subtype$RNA_detected_n / nk_subtype$subtype_cell_n))
assert_check("detected_counts_within_denominators", all(mal_sample$RNA_detected_n <= mal_sample$malignant_cell_n) &&
               all(nk_sample$RNA_detected_n <= nk_sample$NK_cell_n) &&
               all(nk_subtype$RNA_detected_n <= nk_subtype$subtype_cell_n))

primary_exact <- fread(file.path(out, "primary_matched_stk31_hla_nk.csv"))
primary_all <- fread(file.path(out, "primary_all_samples_compartment_availability.csv"))
assert_check("primary_all_sample_count", nrow(primary_all) == 86L, nrow(primary_all))
assert_check("primary_exact_match_count", nrow(primary_exact) == 11L, nrow(primary_exact))
assert_check("primary_exact_only_exact_scope", all(primary_exact$standardized_site == "Primary tumor") &&
               all(primary_exact$exact_compartment_matched) &&
               all(primary_exact$match_scope == "exact_same_primary_sample"))
assert_check("primary_exact_unique_samples", !anyDuplicated(primary_exact$sample_name))
assert_check("primary_patient_id_source", identical(unique(primary_exact$patient_id_source), "derived_from_sample_prefix"))

grid <- fread(file.path(out, "primary_threshold_sensitivity_stage3.csv"))
grid_keys <- c("malignant_cell_threshold", "NK_cell_threshold", "STK31_detected_cell_threshold", "KIR3DL1_detected_cell_threshold")
assert_check("threshold_grid_unique_keys", !anyDuplicated(grid[, ..grid_keys]))
assert_check("threshold_grid_values", setequal(grid$malignant_cell_threshold, c(20L,50L,100L,200L)) &&
               setequal(grid$NK_cell_threshold, c(20L,50L,100L,200L)) &&
               setequal(grid$STK31_detected_cell_threshold, c(1L,3L,5L,10L)) &&
               setequal(grid$KIR3DL1_detected_cell_threshold, c(1L,3L,5L,10L)))
assert_check("threshold_grid_count_bounds", all(grid$eligible_patient_n <= grid$eligible_sample_n) &&
               all(grid$eligible_sample_n <= grid$base_exact_primary_sample_n) &&
               all(grid$base_exact_primary_sample_n == 11L))
assert_check("threshold_grid_patient_rule", identical(unique(grid$patient_count_rule), "any_exact_primary_sample_passes"))

is_monotone <- function(data, varied, fixed) {
  tests <- data[order(get(varied)), all(diff(eligible_sample_n) <= 0) && all(diff(eligible_patient_n) <= 0), by = fixed]
  all(tests$V1)
}
monotone_ok <- all(vapply(grid_keys, function(varied) is_monotone(grid, varied, setdiff(grid_keys, varied)), logical(1)))
assert_check("threshold_grid_monotone_nonincreasing", monotone_ok)

read_mtx_header <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection))
  banner <- readLines(connection, n = 1L)
  line <- readLines(connection, n = 1L)
  while (length(line) && startsWith(line, "%")) line <- readLines(connection, n = 1L)
  list(banner = banner, dims = scan(text = line, quiet = TRUE))
}
full_header <- read_mtx_header(file.path(mal_dir, "10X_counts", "matrix.mtx.gz"))
target_header <- read_mtx_header(file.path(mal_dir, "target_counts", "matrix.mtx.gz"))
assert_check("full_malignant_matrix_banner", identical(full_header$banner, "%%MatrixMarket matrix coordinate integer general"))
assert_check("full_malignant_matrix_header", identical(as.numeric(full_header$dims), c(32137, 96942, 203691975)), paste(full_header$dims, collapse = "x"))
assert_check("target_malignant_matrix_header", identical(as.numeric(target_header$dims), c(9, 96942, 384581)), paste(target_header$dims, collapse = "x"))

read_gz_lines <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection))
  readLines(connection, warn = FALSE)
}
full_features <- read_gz_lines(file.path(mal_dir, "10X_counts", "features.tsv.gz"))
full_barcodes <- read_gz_lines(file.path(mal_dir, "10X_counts", "barcodes.tsv.gz"))
target_features <- read_gz_lines(file.path(mal_dir, "target_counts", "features.tsv.gz"))
target_barcodes <- read_gz_lines(file.path(mal_dir, "target_counts", "barcodes.tsv.gz"))
assert_check("full_features_count_unique", length(full_features) == 32137L && !anyDuplicated(full_features), length(full_features))
assert_check("malignant_barcodes_count_unique", length(full_barcodes) == 96942L && !anyDuplicated(full_barcodes), length(full_barcodes))
assert_check("target_features_count_unique", length(target_features) == 9L && !anyDuplicated(target_features), length(target_features))
assert_check("full_target_barcode_identity", identical(full_barcodes, target_barcodes))

audit <- fread(file.path(out, "stage3_stream_extraction_audit.txt"), sep = "=", header = FALSE,
               col.names = c("key", "value"), fill = TRUE)
get_audit <- function(requested_key) as.numeric(audit[key == requested_key, value])
assert_check("stream_actual_equals_declared_nnz", get_audit("actual_nnz") == get_audit("declared_nnz") &&
               get_audit("actual_nnz") == 1457063953)
assert_check("stream_full_output_nnz", get_audit("malignant_output_nnz") == full_header$dims[3])
assert_check("stream_target_output_nnz", get_audit("target_output_nnz") == target_header$dims[3])
assert_check("stream_no_adjacent_duplicates", get_audit("adjacent_duplicate_coordinates") == 0)

result <- rbindlist(checks)
fwrite(result, file.path(out, "stage3_validation_checks.csv"))
summary_lines <- c(
  "STAGE3_VALIDATION_OK",
  paste0("checks_total=", nrow(result)),
  paste0("checks_passed=", sum(result$passed)),
  paste0("checks_failed=", sum(!result$passed)),
  paste0("validated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(summary_lines, file.path(out, "stage3_validation_summary.txt"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")
