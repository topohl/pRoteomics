#!/usr/bin/env Rscript
# Global raw-derived joint compartment QC preprocessing. It is intentionally
# separate from the three dataset-specific ProTIGY/WGCNA contracts.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "joint_compartment_qc_utils.R"))

global_arg <- tolower(script_arg_value("--dataset", "all"))
if (!global_arg %in% c("all", "global")) stop("This is a global preprocessing script; use --dataset all or global.", call. = FALSE)
runtime <- list(script = "01_preprocessing/01_prepare_joint_protigy_input.r", stage = "joint_qc_preprocessing", dataset = "global", args = commandArgs(trailingOnly = TRUE), dry_run = is_dry_run(), started_at = Sys.time())
raw_file <- Sys.getenv("PROTEOMICS_JOINT_RAW_MATRIX_FILE", unset = path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"))
metadata_file <- Sys.getenv("PROTEOMICS_JOINT_METADATA_FILE", unset = path_metadata("TPE9_sample_metadata_males.xlsx"))
idmap_file <- Sys.getenv("PROTEOMICS_JOINT_IDMAPPING_FILE", unset = path_external("MOUSE_10090_idmapping.dat"))
output_root <- path_processed("01_preprocessing", "joint_compartment_qc", "global")
protigy_file <- path_processed("01_preprocessing", "protigy_input", "global", "joint_shared_core_log2_median_normalized_imputed.gct")
min_block <- joint_qc_env_number("PROTEOMICS_JOINT_MIN_DETECTION_PER_BLOCK", .70, 0, 1)
min_union <- joint_qc_env_number("PROTEOMICS_JOINT_UNION_MIN_DETECTION", .30, 0, 1)
max_imputed <- joint_qc_env_number("PROTEOMICS_JOINT_MAX_IMPUTED_FRACTION_WARN", .05, 0, 1)

if (runtime$dry_run) {
  cat("[DRY-RUN] joint raw matrix: ", normalizePath(raw_file, winslash = "/", mustWork = FALSE), " [", if (file.exists(raw_file)) "PASS" else "FAIL", "]\n", sep = "")
  cat("[DRY-RUN] joint metadata: ", normalizePath(metadata_file, winslash = "/", mustWork = FALSE), " [", if (file.exists(metadata_file)) "PASS" else "FAIL", "]\n", sep = "")
  cat("[DRY-RUN] mouse idmapping: ", normalizePath(idmap_file, winslash = "/", mustWork = FALSE), " [", if (file.exists(idmap_file)) "PASS" else "FAIL", "]\n", sep = "")
  cat("[DRY-RUN] output root: ", output_root, "\n", sep = "")
  cat("[DRY-RUN] balanced shared-core detection threshold: ", min_block, "; broad-union threshold: ", min_union, "\n", sep = "")
  quit(status = if (all(file.exists(c(raw_file, metadata_file, idmap_file)))) 0 else 1, save = "no")
}

required <- c("readxl", "dplyr", "AnnotationDbi", "org.Mm.eg.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
for (input in list(raw_matrix = raw_file, sample_metadata = metadata_file, mouse_idmapping = idmap_file)) {
  if (!file.exists(input[[1]])) stop("Required input does not exist: ", input[[1]], call. = FALSE)
}
record_input_resolution(script = runtime$script, dataset = "global", stage = runtime$stage, input_name = "joint_raw_matrix", expected_path = path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"), resolved_path = raw_file, resolution_mode = if (identical(normalizePath(raw_file, winslash = "/", mustWork = FALSE), normalizePath(path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"), winslash = "/", mustWork = FALSE))) "canonical" else "explicit_override", producer_script_or_artifact_id = "quicksearch.pg_matrix.tsv")
record_input_resolution(script = runtime$script, dataset = "global", stage = runtime$stage, input_name = "joint_sample_metadata", expected_path = path_metadata("TPE9_sample_metadata_males.xlsx"), resolved_path = metadata_file, resolution_mode = "canonical_or_explicit", producer_script_or_artifact_id = "data/metadata/TPE9_sample_metadata_males.xlsx")
record_input_resolution(script = runtime$script, dataset = "global", stage = runtime$stage, input_name = "joint_mouse_idmapping", expected_path = path_external("MOUSE_10090_idmapping.dat"), resolved_path = idmap_file, resolution_mode = "canonical_or_explicit", producer_script_or_artifact_id = "data/external/MOUSE_10090_idmapping.dat")

raw <- utils::read.delim(raw_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"), quote = "")
metadata0 <- as.data.frame(readxl::read_xlsx(metadata_file), check.names = FALSE)
annotation_columns <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
if (length(setdiff(annotation_columns, names(raw)))) stop("Raw unified matrix lacks required annotation columns: ", paste(setdiff(annotation_columns, names(raw)), collapse = ", "), call. = FALSE)
sample_columns <- setdiff(names(raw), annotation_columns)
if (length(sample_columns) < 3L) stop("Raw unified matrix has fewer than three quantitative sample columns.", call. = FALSE)
metadata <- joint_qc_prepare_metadata(metadata0, sample_columns)
if (any(!metadata$dataset %in% joint_qc_datasets())) stop("Canonical metadata does not assign every matched sample to neuron_neuropil, neuron_soma, or microglia.", call. = FALSE)
metadata_alignment_audit <- data.frame(n_matrix_sample_columns = length(sample_columns), n_metadata_rows = nrow(metadata0), n_exact_matched_samples = nrow(metadata), n_unmatched_matrix_samples = 0L, n_unmatched_metadata_rows = 0L, n_duplicated_matrix_sample_ids = 0L, n_duplicated_metadata_sample_ids = 0L, stringsAsFactors = FALSE)
sample_selection_audit <- metadata[, c("Sample", "dataset", "exclusion_status", "background_or_blank"), drop = FALSE]
keep_samples <- metadata$exclusion_status == "included" & !metadata$background_or_blank
if (!any(keep_samples)) stop("No non-excluded, non-background samples remain.", call. = FALSE)
metadata <- metadata[keep_samples, , drop = FALSE]
sample_columns <- sample_columns[keep_samples]

raw$.source_row_id <- seq_len(nrow(raw))
mapping <- qc_build_mapping_context(idmap_file = idmap_file)
canonical <- build_canonical_protein_group_tables(raw, dataset = "global", source_file = raw_file,
  entry_map = mapping$entry_map, gene_map = mapping$gene_map, accession_gene_map = mapping$accession_gene_map,
  reviewed_map = mapping$reviewed_map, manual_mapping = mapping$manual_mapping, gene_annotation_maps = mapping$gene_annotation_maps,
  manual_gene_annotation_overrides = mapping$manual_gene_annotation_overrides, uniprot_mapping_file_hash = mapping$idmap_sha256,
  strict = TRUE, identifier_col = "Protein.Group", feature_col = "Protein.Group")
feature_table <- canonical$wide
if (anyNA(feature_table$ProteinGroupID) || anyDuplicated(feature_table$ProteinGroupID)) stop("Canonical ProteinGroupID is missing or duplicated; refusing to repair identity.", call. = FALSE)
feature_table$FeatureDisplayLabel <- ifelse(is.na(feature_table$official_gene_symbol) | !nzchar(feature_table$official_gene_symbol), feature_table$original_identifier, feature_table$official_gene_symbol)
feature_table$joint_qc_eligible <- feature_table$protein_group_ambiguity_class != "mixed_species_or_contaminant"
feature_table$joint_qc_exclusion_reason <- ifelse(feature_table$joint_qc_eligible, NA_character_, "mixed_species_or_contaminant")

raw_numeric <- as.matrix(data.frame(lapply(raw[sample_columns], function(x) suppressWarnings(as.numeric(as.character(x)))), check.names = FALSE))
rownames(raw_numeric) <- feature_table$ProteinGroupID
colnames(raw_numeric) <- sample_columns
joint_validate_canonical_feature_matrix(feature_table, raw_numeric, "raw quantitative matrix")
positive <- joint_qc_log2_positive(raw_numeric)
joint_validate_canonical_feature_matrix(feature_table, positive$raw, "positive raw quantitative matrix")
joint_validate_canonical_feature_matrix(feature_table, positive$log2, "positive log2 quantitative matrix")
eligible_ids <- feature_table$ProteinGroupID[feature_table$joint_qc_eligible]
eligible_log2 <- joint_validate_feature_subset(eligible_ids, positive$log2, "eligible canonical feature subset")
technical <- joint_qc_select_technical_block(metadata)
if (!is.na(technical$warning)) warning(technical$warning, call. = FALSE)
metadata$compartment <- metadata$dataset
celltype_col <- joint_qc_first_col(metadata, c("celltype", "CellType"))
if (!is.na(celltype_col)) metadata$cell_type <- as.character(metadata[[celltype_col]])
region_col <- joint_qc_first_col(metadata, c("region", "Region")); layer_col <- joint_qc_first_col(metadata, c("layer", "Layer"))
if (!is.na(region_col)) metadata$region <- as.character(metadata[[region_col]])
if (!is.na(layer_col)) metadata$layer <- as.character(metadata[[layer_col]])
if (!is.na(region_col) && !is.na(layer_col)) metadata$region_layer <- paste(metadata$region, metadata$layer, sep = "_")
group_col <- joint_qc_first_col(metadata, c("experimental_group", "ExpGroup", "group", "Group")); if (!is.na(group_col)) metadata$experimental_group <- as.character(metadata[[group_col]])
animal_col <- joint_qc_first_col(metadata, c("animal_id", "AnimalID")); if (!is.na(animal_col)) metadata$animal_id <- as.character(metadata[[animal_col]])
run_col <- joint_qc_first_col(metadata, c("run", "Run", "run_order")); if (!is.na(run_col)) metadata$run <- as.character(metadata[[run_col]])
metadata$technical_block <- technical$block
sample_missingness <- data.frame(Sample = sample_columns, dataset = metadata$dataset, technical_block = technical$block,
  n_proteins = nrow(positive$log2), n_missing_raw = colSums(is.na(positive$log2)), fraction_missing_raw = colMeans(is.na(positive$log2)), stringsAsFactors = FALSE)
sample_counts <- as.data.frame(table(dataset = metadata$dataset, technical_block = technical$block), stringsAsFactors = FALSE)
names(sample_counts)[names(sample_counts) == "Freq"] <- "n_samples"
filters <- joint_qc_filter_universes(!is.na(eligible_log2), metadata, eligible_ids, technical$block, min_block, min_union)
feature_audit <- feature_table[, c("ProteinGroupID", "source_row_id", "original_identifier", "member_identifiers_original", "member_identifiers_canonical", "member_accessions", "protein_group_ambiguity_class", "joint_qc_eligible", "joint_qc_exclusion_reason"), drop = FALSE]
filter_idx <- match(feature_audit$ProteinGroupID, filters$audit$ProteinGroupID)
for (nm in setdiff(names(filters$audit), "ProteinGroupID")) feature_audit[[nm]] <- filters$audit[[nm]][filter_idx]
feature_audit$exclusion_reason[!feature_audit$joint_qc_eligible] <- feature_audit$joint_qc_exclusion_reason[!feature_audit$joint_qc_eligible]

primary_ids <- eligible_ids[filters$primary]
complete_ids <- eligible_ids[filters$complete]
broad_ids <- eligible_ids[filters$broad]
joint_qc_report_feature_subset_diagnostics(primary_ids, positive$log2, "primary shared-core before subset")
primary_log2 <- joint_validate_feature_subset(primary_ids, positive$log2, "primary shared-core subset")
stopifnot(identical(rownames(primary_log2), primary_ids))
normalization <- joint_qc_joint_median_normalize(primary_log2)
primary_normalized <- joint_validate_feature_subset(primary_ids, normalization$matrix, "primary normalized subset")
stopifnot(identical(rownames(primary_normalized), primary_ids))
imputation <- joint_qc_median_impute(primary_normalized, primary_ids)
imputation$matrix <- joint_validate_feature_subset(primary_ids, imputation$matrix, "primary imputed subset")
stopifnot(identical(rownames(imputation$matrix), primary_ids))
if (anyNA(imputation$matrix)) stop("Primary shared-core matrix still has missing values after deterministic median imputation.", call. = FALSE)
complete_log2 <- joint_validate_feature_subset(complete_ids, positive$log2, "complete-case subset")
stopifnot(identical(rownames(complete_log2), complete_ids))
if (anyNA(complete_log2)) stop("Complete-case sensitivity matrix contains missing values.", call. = FALSE)
broad_log2 <- joint_validate_feature_subset(broad_ids, positive$log2, "broad-union subset")
broad_binary <- (!is.na(broad_log2)) * 1L
rownames(broad_binary) <- broad_ids
colnames(broad_binary) <- colnames(broad_log2)
stopifnot(identical(rownames(broad_binary), broad_ids))

sample_norm <- data.frame(Sample = colnames(primary_log2), dataset = metadata$dataset, technical_block = technical$block,
  pre_normalization_median_log2 = normalization$pre_median, post_normalization_median_log2 = normalization$post_median,
  pre_normalization_q25_log2 = apply(primary_log2, 2, stats::quantile, probs = .25, na.rm = TRUE),
  pre_normalization_q75_log2 = apply(primary_log2, 2, stats::quantile, probs = .75, na.rm = TRUE),
  post_normalization_q25_log2 = apply(normalization$matrix, 2, stats::quantile, probs = .25, na.rm = TRUE),
  post_normalization_q75_log2 = apply(normalization$matrix, 2, stats::quantile, probs = .75, na.rm = TRUE), stringsAsFactors = FALSE)
impute_long <- which(imputation$missing, arr.ind = TRUE)
impute_audit <- if (nrow(impute_long)) data.frame(ProteinGroupID = rownames(imputation$matrix)[impute_long[, 1]], Sample = colnames(imputation$matrix)[impute_long[, 2]], dataset = metadata$dataset[impute_long[, 2]], technical_block = technical$block[impute_long[, 2]], original_value = NA_real_, imputed_value = imputation$matrix[imputation$missing], imputation_method = "global_joint_protein_wise_median", protein_median_used = imputation$protein_median[impute_long[, 1]], stringsAsFactors = FALSE) else data.frame(ProteinGroupID = character(), Sample = character(), dataset = character(), technical_block = character(), original_value = numeric(), imputed_value = numeric(), imputation_method = character(), protein_median_used = numeric())
sample_impute <- data.frame(Sample = colnames(imputation$matrix), dataset = metadata$dataset, technical_block = technical$block, n_imputed = colSums(imputation$missing), fraction_imputed = colMeans(imputation$missing), stringsAsFactors = FALSE)
block_key <- paste(metadata$dataset, technical$block, sep = "\r")
impute_by_block <- lapply(unique(block_key), function(key) {
  keep <- block_key == key
  data.frame(level = "dataset_x_technical_block", dataset = metadata$dataset[which(keep)[1]], technical_block = technical$block[which(keep)[1]], n_imputed = sum(imputation$missing[, keep, drop = FALSE]), n_values = length(imputation$missing[, keep, drop = FALSE]), fraction_imputed = mean(imputation$missing[, keep, drop = FALSE]))
})
impute_summary <- do.call(rbind, c(list(data.frame(level = "global", dataset = "all", technical_block = "all", n_imputed = sum(imputation$missing), n_values = length(imputation$missing), fraction_imputed = mean(imputation$missing))), impute_by_block))
if (mean(imputation$missing) > max_imputed) warning("Primary shared-core imputation fraction exceeds PROTEOMICS_JOINT_MAX_IMPUTED_FRACTION_WARN (", max_imputed, "): ", signif(mean(imputation$missing), 3), call. = FALSE)

dir_create(output_root); dir_create(dirname(protigy_file))
write_csv <- function(x, name) utils::write.csv(x, file.path(output_root, name), row.names = FALSE, na = "")
write_csv(feature_table, "canonical_feature_table.csv"); write_csv(canonical$bridge, "protein_group_member_bridge.csv")
write_csv(feature_audit, "feature_filter_audit.csv"); write_csv(filters$dataset_long, "observed_detection_by_dataset.csv"); write_csv(filters$block_long, "observed_detection_by_dataset_technical_block.csv")
write_csv(canonical$collision_audit, "protein_group_id_collision_audit.csv")
write_csv(feature_table[!feature_table$joint_qc_eligible, , drop = FALSE], "protein_group_exclusions.csv")
write_csv(feature_table[, c("ProteinGroupID", "source_file", "source_feature_id", "source_row_id", "original_identifier", "member_identifiers_original", "member_identifiers_canonical")], "source_row_provenance.csv")
write_csv(metadata_alignment_audit, "sample_metadata_alignment_audit.csv"); write_csv(sample_selection_audit, "sample_selection_exclusion_background_audit.csv"); write_csv(sample_counts, "compartment_technical_block_sample_counts.csv"); write_csv(sample_missingness, "sample_missingness_raw.csv")
write_csv(data.frame(annotation_column = annotation_columns, present = annotation_columns %in% names(raw), stringsAsFactors = FALSE), "raw_annotation_column_audit.csv"); write_csv(data.frame(quantitative_sample_column = names(raw_numeric), stringsAsFactors = FALSE), "raw_quantitative_sample_column_audit.csv")
write_csv(sample_norm, "sample_normalization_audit.csv"); write_csv(sample_impute, "imputation_footprint_by_sample.csv"); write_csv(imputation$audit, "imputation_footprint_by_protein.csv"); write_csv(impute_audit, "row_level_imputation_audit.csv"); write_csv(impute_summary, "imputation_summary_by_dataset_technical_block.csv")
write_csv(data.frame(input_name = c("raw_matrix", "metadata", "idmapping"), path = normalizePath(c(raw_file, metadata_file, idmap_file), winslash = "/", mustWork = TRUE), sha256 = vapply(c(raw_file, metadata_file, idmap_file), file_hash_sha256, character(1)), stringsAsFactors = FALSE), "input_manifest.csv")
write_csv(data.frame(metric = c("raw_rows", "raw_columns", "quantitative_samples", "retained_samples", "primary_shared_core_proteins", "complete_case_proteins", "broad_union_proteins", "primary_imputed_fraction", "zero_values", "nonfinite_values"), value = c(nrow(raw), ncol(raw), ncol(raw_numeric), nrow(metadata), length(primary_ids), length(complete_ids), length(broad_ids), mean(imputation$missing), sum(positive$zero), sum(positive$nonfinite)), stringsAsFactors = FALSE), "input_and_matrix_audit.csv")
joint_qc_write_matrix_tsv(imputation$matrix, file.path(output_root, "joint_shared_core_log2_median_normalized_imputed.tsv")); joint_qc_write_matrix_tsv(primary_log2, file.path(output_root, "joint_shared_core_log2_unnormalized.tsv")); joint_qc_write_matrix_tsv(complete_log2, file.path(output_root, "joint_complete_case_log2.tsv")); joint_qc_write_matrix_tsv(broad_binary, file.path(output_root, "joint_broad_union_detected_binary.tsv"))
joint_qc_write_gct_v13(imputation$matrix, metadata, feature_table, protigy_file)
saveRDS(list(contract_version = "joint_compartment_qc_v1", metadata = metadata, technical_block = technical, primary = list(matrix = imputation$matrix, unnormalized_log2 = primary_log2, feature_ids = primary_ids), complete_case = list(matrix = complete_log2, feature_ids = complete_ids), broad_union = list(detected_binary = broad_binary, feature_ids = broad_ids), feature_table = feature_table, feature_filter_audit = feature_audit), file.path(output_root, "joint_compartment_qc_matrices.rds"))
if (requireNamespace("writexl", quietly = TRUE)) writexl::write_xlsx(list(feature_filter = feature_audit, sample_normalization = sample_norm, sample_imputation = sample_impute, observed_detection = filters$dataset_long), file.path(output_root, "joint_compartment_qc_audit_bundle.xlsx"))
writeLines(c("# Joint compartment preprocessing", paste0("Input matrix: ", normalizePath(raw_file, winslash = "/", mustWork = TRUE)), paste0("Retained samples: ", nrow(metadata)), paste0("Technical block: ", technical$variable, " (", technical$mode, ")"), paste0("Primary/complete/union proteins: ", length(primary_ids), "/", length(complete_ids), "/", length(broad_ids)), paste0("Primary imputed fraction: ", signif(mean(imputation$missing), 4)), "Primary analysis is raw-positive log2, one joint sample-wise median normalization, then label-blind protein-wise median imputation. Microglia samples are microglia-enriched ROIs, not purified microglia."), file.path(output_root, "README.md"))
write_run_manifest(file.path(path_results("logs", "01_preprocessing", "01_prepare_joint_protigy_input", "global"), "run_manifest.yml"), inputs = list(raw_matrix = raw_file, metadata = metadata_file, idmapping = idmap_file), outputs = list(processed = output_root, protigy_gct = protigy_file), parameters = list(min_detection_per_dataset_batch = min_block, union_min_detection = min_union, normalization = "joint_sample_wise_median", imputation = "global_joint_protein_wise_median"), notes = "Global QC only; no combined DE/WGCNA or primary batch correction.")
