#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project <- "/home/zhuweiyu/codex-r/GBC_STK31_public_scRNA"
processed <- file.path(project, "02_processed_data", "All_single_cells")
counts_dir <- file.path(processed, "counts", "10X_counts")
mal_dir <- file.path(project, "02_processed_data", "malignant_epithelial")
mal_10x <- file.path(mal_dir, "10X_counts")
target_dir <- file.path(mal_dir, "target_counts")
stage3_out <- file.path(project, "06_feasibility", "stage3")
dir.create(mal_10x, recursive = TRUE, showWarnings = FALSE)
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage3_out, recursive = TRUE, showWarnings = FALSE)

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

message("Reading and validating cell metadata")
meta <- fread(file.path(processed, "GBC_Metadata.txt"), sep = "\t", showProgress = TRUE)
schema_row <- meta$name == "type" & meta$celltype == "group" & meta$subtype == "group"
stop_if_not(sum(schema_row) == 1L, "Expected exactly one metadata schema row")
meta <- meta[!schema_row]
stop_if_not(identical(names(meta), c("name", "celltype", "subtype", "sample_name")),
            "Unexpected metadata columns")
meta[, full_col := .I]

message("Reading and validating 10x barcodes/features")
barcodes <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(counts_dir, "barcodes.tsv.gz"))),
  header = FALSE, col.names = "barcode", showProgress = FALSE
)
features <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(counts_dir, "features.tsv.gz"))),
  header = FALSE,
  col.names = c("feature_id", "feature_name", "feature_type"),
  showProgress = FALSE
)
stop_if_not(nrow(meta) == 1117245L, "Unexpected metadata cell count")
stop_if_not(nrow(barcodes) == nrow(meta), "Barcode/metadata size mismatch")
stop_if_not(!anyDuplicated(meta$name), "Duplicated metadata names")
stop_if_not(!anyDuplicated(barcodes$barcode), "Duplicated barcodes")
stop_if_not(identical(meta$name, barcodes$barcode), "Barcode order differs from metadata; positional extraction aborted")
stop_if_not(nrow(features) == 32137L, "Unexpected feature count")
stop_if_not(!anyDuplicated(features$feature_name), "Duplicated feature symbols")
stop_if_not(all(features$feature_type == "Gene Expression"), "Non-Gene-Expression feature found")

target_genes <- data.table(
  gene = c("STK31", "HLA-A", "HLA-B", "HLA-C", "B2M", "NLRC5", "TAP1", "TAP2", "HLA-E"),
  panel = c("prespecified", rep("prespecified", 7L), "supplementary_KLRC1_ligand")
)
target_genes[, full_row := match(gene, features$feature_name)]
stop_if_not(!anyNA(target_genes$full_row), "A malignant target gene is missing from features")
target_genes[, target_row := .I]

message("Reading sample_info.xlsx without modifying the workbook")
sample_info <- as.data.table(read_excel(
  file.path(processed, "sample_info.xlsx"), skip = 1, .name_repair = "unique"
))
sample_info <- sample_info[!is.na(`Sample ID`) & trimws(as.character(`Sample ID`)) != ""]
sample_info[, `Sample ID` := as.character(`Sample ID`)]
stop_if_not(!anyDuplicated(sample_info$`Sample ID`), "Duplicated Sample ID in sample_info")
sample_info[, patient_id := derive_patient_id(`Sample ID`)]
sample_info[, standardized_site := standardize_site(Site, `Histological type`)]

message("Building malignant epithelial manifest")
malignant <- meta[subtype == "Malignant epithelial cells"]
stop_if_not(nrow(malignant) == 96942L, "Unexpected malignant epithelial cell count")
malignant[, malignant_col := .I]
clinical_index <- match(malignant$sample_name, sample_info$`Sample ID`)
stop_if_not(!anyNA(clinical_index), "A malignant cell sample does not map to sample_info")
malignant[, `:=`(
  barcode = name,
  patient_id = sample_info$patient_id[clinical_index],
  patient_id_source = "derived_from_sample_prefix",
  site = as.character(sample_info$Site[clinical_index]),
  standardized_site = sample_info$standardized_site[clinical_index],
  TNM_stage = as.character(sample_info$`TNM stage`[clinical_index]),
  histological_type = as.character(sample_info$`Histological type`[clinical_index]),
  metastatic_group = as.character(sample_info$`Metastatic group`[clinical_index])
)]
malignant_manifest <- malignant[, .(
  malignant_col, full_col, barcode, sample_name, patient_id,
  patient_id_source, site, standardized_site, TNM_stage,
  histological_type, metastatic_group, celltype, subtype
)]
stop_if_not(identical(malignant_manifest$barcode, meta$name[malignant_manifest$full_col]),
            "Malignant manifest order validation failed")

fwrite(malignant_manifest, file.path(mal_dir, "malignant_cell_metadata.csv.gz"), compress = "gzip")
fwrite(malignant_manifest, file.path(stage3_out, "malignant_cell_manifest.csv.gz"), compress = "gzip")
fwrite(malignant_manifest[, .(barcode)], file.path(mal_10x, "barcodes.tsv.gz"),
       sep = "\t", col.names = FALSE, compress = "gzip")
fwrite(features, file.path(mal_10x, "features.tsv.gz"),
       sep = "\t", col.names = FALSE, compress = "gzip")
fwrite(malignant_manifest[, .(barcode)], file.path(target_dir, "barcodes.tsv.gz"),
       sep = "\t", col.names = FALSE, compress = "gzip")
fwrite(target_genes[, .(gene, gene, feature_type = "Gene Expression")],
       file.path(target_dir, "features.tsv.gz"), sep = "\t", col.names = FALSE, compress = "gzip")
fwrite(target_genes, file.path(stage3_out, "malignant_target_feature_manifest.csv"))

col_map <- integer(nrow(meta))
col_map[malignant_manifest$full_col] <- malignant_manifest$malignant_col
row_map <- integer(nrow(features))
row_map[target_genes$full_row] <- target_genes$target_row
col_con <- file(file.path(mal_dir, "full_to_malignant_col.int32.bin"), "wb")
writeBin(as.integer(col_map), col_con, size = 4L, endian = "little")
close(col_con)
row_con <- file(file.path(mal_dir, "full_to_target_row.int32.bin"), "wb")
writeBin(as.integer(row_map), row_con, size = 4L, endian = "little")
close(row_con)

audit <- c(
  paste0("metadata_cells=", nrow(meta)),
  paste0("features=", nrow(features)),
  paste0("malignant_cells=", nrow(malignant_manifest)),
  paste0("malignant_samples=", uniqueN(malignant_manifest$sample_name)),
  paste0("barcode_order_identical=", identical(meta$name, barcodes$barcode)),
  paste0("patient_id_source=derived_from_sample_prefix"),
  paste0("target_rows=", paste(target_genes$gene, target_genes$full_row, sep = ":", collapse = "|"))
)
writeLines(audit, file.path(stage3_out, "stage3_manifest_audit.txt"), useBytes = TRUE)
cat(paste(audit, collapse = "\n"), "\n")
