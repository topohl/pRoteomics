# ================================================================
# Script: 02_id_mapping/01b_animal_level_mapping_audit.r
# Scope: validation-only audit of isolated animal-level forward mapping
# Notes: Does not map proteins and does not write to historical mapping trees.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "mapping_branch_utils.R"))

datasets <- valid_datasets()
direction <- "forward"
corrected_extract_root <- path_processed("01_preprocessing", "gct_extractR_animal_level")
corrected_mapping_root <- path_processed("02_id_mapping_animal_level")
legacy_extract_root <- path_processed("01_preprocessing", "gct_extractR")
legacy_mapping_root <- path_processed("02_id_mapping")
audit_table_root <- path_results("tables", "02_id_mapping_animal_level", "mapping_branch_audit")
audit_report_root <- path_results("reports", "02_id_mapping_animal_level", "mapping_branch_audit")
audit_log_root <- path_results("logs", "02_id_mapping_animal_level", "mapping_branch_audit")
invisible(lapply(c(audit_table_root, audit_report_root, audit_log_root), dir_create))

read_contract_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    na.strings = ""
  )
}

is_missing_text <- function(x) is.na(x) | !nzchar(trimws(x)) | trimws(x) == "NA"
id_key <- function(x) ifelse(is.na(x), "<MISSING>", as.character(x))
same_text <- function(left, right) {
  left_missing <- is_missing_text(left)
  right_missing <- is_missing_text(right)
  (left_missing & right_missing) | (!left_missing & !right_missing & left == right)
}

contrast_files <- function(path, exclude_global = FALSE) {
  files <- sort(list.files(path, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE))
  if (exclude_global) files <- files[!grepl("^GLOBAL_", basename(files))]
  files
}

required_sus_res_files <- list(
  neuron_neuropil = c(
    "CA1slmsus_CA1slmres.csv", "CA1sosus_CA1sores.csv", "CA1srsus_CA1srres.csv",
    "CA2slmsus_CA2slmres.csv", "CA2sosus_CA2sores.csv", "CA2srsus_CA2srres.csv",
    "CA3sosus_CA3sores.csv", "CA3srsus_CA3srres.csv",
    "DGmosus_DGmores.csv", "DGposus_DGpores.csv"
  ),
  neuron_soma = c(
    "CA1spsus_CA1spres.csv", "CA2spsus_CA2spres.csv",
    "CA3spsus_CA3spres.csv", "DGsgsus_DGsgres.csv"
  ),
  microglia = c(
    "CA1microgliasus_CA1microgliares.csv", "CA2microgliasus_CA2microgliares.csv",
    "CA3microgliasus_CA3microgliares.csv", "DGmicrogliasus_DGmicrogliares.csv"
  )
)

file_validation_all <- list()
determinism_all <- list()
legacy_only_all <- list()
corrected_only_all <- list()
primary_all <- list()
manifest_all <- list()
dataset_summary_all <- list()

for (dataset in datasets) {
  expected_count <- unname(animal_level_mapping_expected_counts()[[dataset]])
  input_dir <- file.path(corrected_extract_root, dataset, direction)
  mapped_dir <- file.path(corrected_mapping_root, "mapped", dataset, direction, "per_file")
  unmapped_dir <- file.path(corrected_mapping_root, "unmapped", dataset, direction, "per_file")
  bridge_dir <- file.path(corrected_mapping_root, "member_bridge", dataset, direction, "per_file")
  legacy_input_dir <- file.path(legacy_extract_root, dataset, direction)
  legacy_mapped_dir <- file.path(legacy_mapping_root, "mapped", dataset, direction, "per_file")

  input_files <- contrast_files(input_dir)
  mapped_files <- contrast_files(mapped_dir)
  unmapped_files <- contrast_files(unmapped_dir, exclude_global = TRUE)
  bridge_files <- contrast_files(bridge_dir)
  names_in <- basename(input_files)
  names_mapped <- basename(mapped_files)
  names_unmapped <- basename(unmapped_files)
  names_bridge <- sub("_member_bridge\\.csv$", ".csv", basename(bridge_files))

  primary <- data.frame(
    dataset = dataset,
    expected_file = required_sus_res_files[[dataset]],
    input_present = required_sus_res_files[[dataset]] %in% names_in,
    mapped_present = required_sus_res_files[[dataset]] %in% names_mapped,
    stringsAsFactors = FALSE
  )
  primary_all[[dataset]] <- primary

  file_rows <- vector("list", length(input_files))
  discrepancy_rows <- list()
  legacy_only_rows <- list()
  corrected_only_rows <- list()

  for (i in seq_along(input_files)) {
    input_file <- input_files[[i]]
    filename <- basename(input_file)
    mapped_file <- file.path(mapped_dir, filename)
    unmapped_file <- file.path(unmapped_dir, filename)
    bridge_file <- file.path(bridge_dir, sub("\\.csv$", "_member_bridge.csv", filename))
    legacy_input_file <- file.path(legacy_input_dir, filename)
    legacy_mapped_file <- file.path(legacy_mapped_dir, filename)
    required_outputs_exist <- all(file.exists(c(mapped_file, unmapped_file, bridge_file)))
    failure <- character()

    if (!required_outputs_exist) {
      failure <- c(failure, "missing corrected mapped/unmapped/member-bridge output")
      file_rows[[i]] <- data.frame(
        dataset = dataset, source_file = filename, input_protein_rows = NA_integer_,
        mapped_output_rows = NA_integer_, unmapped_member_rows = NA_integer_,
        member_bridge_rows = NA_integer_, mapped_protein_count = NA_integer_,
        unmapped_protein_count = NA_integer_, mapping_percentage = NA_real_,
        duplicate_final_accession_rows = NA_integer_, duplicate_final_accession_values = NA_integer_,
        missing_SYMBOL_count = NA_integer_, missing_ENTREZ_count = NA_integer_,
        filename_preserved = FALSE, direction_preserved = FALSE,
        failed = TRUE, failure_reason = paste(failure, collapse = "; "), stringsAsFactors = FALSE
      )
      next
    }

    input <- read_contract_csv(input_file)
    mapped <- read_contract_csv(mapped_file)
    unmapped <- read_contract_csv(unmapped_file)
    bridge <- read_contract_csv(bridge_file)
    required_columns <- c(
      "original_identifier", "representative_accession", "official_gene_symbol",
      "official_entrez_id", "mapping_status", "source_file", "ProteinGroupID"
    )
    missing_columns <- setdiff(required_columns, names(mapped))
    if (length(missing_columns)) failure <- c(failure, paste0("missing mapped columns: ", paste(missing_columns, collapse = ",")))

    input_ids <- id_key(input[[1]])
    mapped_ids <- if ("original_identifier" %in% names(mapped)) id_key(mapped$original_identifier) else character()
    if (nrow(mapped) != nrow(input)) failure <- c(failure, "mapped row count differs from input")
    if (anyDuplicated(input_ids)) failure <- c(failure, "duplicate source protein IDs in corrected input")
    if (anyDuplicated(mapped_ids)) failure <- c(failure, "duplicate source protein IDs in corrected mapped output")
    if (!setequal(input_ids, mapped_ids)) failure <- c(failure, "source protein ID set changed during mapping")

    filename_preserved <- "source_file" %in% names(mapped) &&
      length(unique(mapped$source_file)) == 1L && identical(unique(mapped$source_file), filename) &&
      filename %in% names_mapped && filename %in% names_unmapped && filename %in% names_bridge
    direction_preserved <- grepl("(res_.+con|sus_.+res|sus_.+con)\\.csv$", filename) &&
      !grepl("(con_.+res|res_.+sus|con_.+sus)\\.csv$", filename)
    if (!filename_preserved) failure <- c(failure, "comparison filename/source_file not preserved")
    if (!direction_preserved) failure <- c(failure, "unexpected biological comparison direction")

    accession <- if ("representative_accession" %in% names(mapped)) mapped$representative_accession else rep(NA_character_, nrow(mapped))
    accession_nonmissing <- accession[!is_missing_text(accession)]
    mapped_status <- if ("mapping_status" %in% names(mapped)) mapped$mapping_status else rep(NA_character_, nrow(mapped))
    mapped_count <- sum(mapped_status == "mapped", na.rm = TRUE)
    unmapped_count <- nrow(mapped) - mapped_count
    missing_symbol <- if ("official_gene_symbol" %in% names(mapped)) sum(is_missing_text(mapped$official_gene_symbol)) else nrow(mapped)
    missing_entrez <- if ("official_entrez_id" %in% names(mapped)) sum(is_missing_text(mapped$official_entrez_id)) else nrow(mapped)

    file_rows[[i]] <- data.frame(
      dataset = dataset,
      source_file = filename,
      input_protein_rows = nrow(input),
      mapped_output_rows = nrow(mapped),
      unmapped_member_rows = nrow(unmapped),
      member_bridge_rows = nrow(bridge),
      mapped_protein_count = mapped_count,
      unmapped_protein_count = unmapped_count,
      mapping_percentage = if (nrow(mapped)) 100 * mapped_count / nrow(mapped) else NA_real_,
      duplicate_final_accession_rows = sum(duplicated(accession_nonmissing)),
      duplicate_final_accession_values = sum(table(accession_nonmissing) > 1L),
      missing_SYMBOL_count = missing_symbol,
      missing_ENTREZ_count = missing_entrez,
      filename_preserved = filename_preserved,
      direction_preserved = direction_preserved,
      failed = length(failure) > 0L,
      failure_reason = paste(failure, collapse = "; "),
      stringsAsFactors = FALSE
    )

    if (!file.exists(legacy_input_file) || !file.exists(legacy_mapped_file)) {
      discrepancy_rows[[filename]] <- data.frame(
        dataset = dataset, source_file = filename, source_protein_id = NA_character_,
        legacy_final_accession = NA_character_, corrected_final_accession = NA_character_,
        legacy_SYMBOL = NA_character_, corrected_SYMBOL = NA_character_,
        legacy_ENTREZ = NA_character_, corrected_ENTREZ = NA_character_,
        discrepancy = "missing paired legacy input or mapped file", stringsAsFactors = FALSE
      )
      next
    }

    legacy <- read_contract_csv(legacy_mapped_file)
    legacy_keys <- id_key(legacy$original_identifier)
    corrected_keys <- mapped_ids
    shared_keys <- intersect(legacy_keys, corrected_keys)
    legacy_index <- match(shared_keys, legacy_keys)
    corrected_index <- match(shared_keys, corrected_keys)
    legacy_shared <- legacy[legacy_index, , drop = FALSE]
    corrected_shared <- mapped[corrected_index, , drop = FALSE]
    same_accession <- same_text(legacy_shared$representative_accession, corrected_shared$representative_accession)
    same_symbol <- same_text(legacy_shared$official_gene_symbol, corrected_shared$official_gene_symbol)
    same_entrez <- same_text(legacy_shared$official_entrez_id, corrected_shared$official_entrez_id)
    discrepant <- !(same_accession & same_symbol & same_entrez)
    if (any(discrepant)) {
      discrepancy_rows[[filename]] <- data.frame(
        dataset = dataset,
        source_file = filename,
        source_protein_id = shared_keys[discrepant],
        legacy_final_accession = legacy_shared$representative_accession[discrepant],
        corrected_final_accession = corrected_shared$representative_accession[discrepant],
        legacy_SYMBOL = legacy_shared$official_gene_symbol[discrepant],
        corrected_SYMBOL = corrected_shared$official_gene_symbol[discrepant],
        legacy_ENTREZ = legacy_shared$official_entrez_id[discrepant],
        corrected_ENTREZ = corrected_shared$official_entrez_id[discrepant],
        discrepancy = "mapping contract differs for shared source ID",
        stringsAsFactors = FALSE
      )
    }
    legacy_only <- setdiff(legacy_keys, corrected_keys)
    corrected_only <- setdiff(corrected_keys, legacy_keys)
    if (length(legacy_only)) {
      idx <- match(legacy_only, legacy_keys)
      legacy_only_rows[[filename]] <- data.frame(
        dataset = dataset, source_file = filename, source_protein_id = legacy_only,
        final_accession = legacy$representative_accession[idx],
        SYMBOL = legacy$official_gene_symbol[idx], ENTREZ = legacy$official_entrez_id[idx],
        stringsAsFactors = FALSE
      )
    }
    if (length(corrected_only)) {
      idx <- match(corrected_only, corrected_keys)
      corrected_only_rows[[filename]] <- data.frame(
        dataset = dataset, source_file = filename, source_protein_id = corrected_only,
        final_accession = mapped$representative_accession[idx],
        SYMBOL = mapped$official_gene_symbol[idx], ENTREZ = mapped$official_entrez_id[idx],
        stringsAsFactors = FALSE
      )
    }

    for (entry in list(
      c(role = "corrected_input", path = input_file),
      c(role = "corrected_mapped", path = mapped_file),
      c(role = "legacy_mapped", path = legacy_mapped_file)
    )) {
      manifest_all[[length(manifest_all) + 1L]] <- data.frame(
        dataset = dataset, source_file = filename, role = entry[["role"]],
        path = normalizePath(entry[["path"]], winslash = "/", mustWork = TRUE),
        sha256 = file_hash_sha256(entry[["path"]]), stringsAsFactors = FALSE
      )
    }
  }

  file_validation <- do.call(rbind, file_rows)
  determinism <- if (length(discrepancy_rows)) do.call(rbind, discrepancy_rows) else data.frame(
    dataset = character(), source_file = character(), source_protein_id = character(),
    legacy_final_accession = character(), corrected_final_accession = character(),
    legacy_SYMBOL = character(), corrected_SYMBOL = character(), legacy_ENTREZ = character(),
    corrected_ENTREZ = character(), discrepancy = character(), stringsAsFactors = FALSE
  )
  legacy_only_df <- if (length(legacy_only_rows)) do.call(rbind, legacy_only_rows) else data.frame(
    dataset = character(), source_file = character(), source_protein_id = character(),
    final_accession = character(), SYMBOL = character(), ENTREZ = character(), stringsAsFactors = FALSE
  )
  corrected_only_df <- if (length(corrected_only_rows)) do.call(rbind, corrected_only_rows) else data.frame(
    dataset = character(), source_file = character(), source_protein_id = character(),
    final_accession = character(), SYMBOL = character(), ENTREZ = character(), stringsAsFactors = FALSE
  )

  dataset_dir <- file.path(audit_table_root, dataset)
  dir_create(dataset_dir)
  utils::write.csv(file_validation, file.path(dataset_dir, "file_level_validation.csv"), row.names = FALSE, na = "")
  utils::write.csv(determinism, file.path(dataset_dir, "mapping_determinism_discrepancies.csv"), row.names = FALSE, na = "")
  utils::write.csv(legacy_only_df, file.path(dataset_dir, "legacy_only_source_ids.csv"), row.names = FALSE, na = "")
  utils::write.csv(corrected_only_df, file.path(dataset_dir, "corrected_only_source_ids.csv"), row.names = FALSE, na = "")
  utils::write.csv(primary, file.path(dataset_dir, "required_sus_vs_res_files.csv"), row.names = FALSE, na = "")

  file_validation_all[[dataset]] <- file_validation
  determinism_all[[dataset]] <- determinism
  legacy_only_all[[dataset]] <- legacy_only_df
  corrected_only_all[[dataset]] <- corrected_only_df
  unique_legacy_only <- unique(legacy_only_df$source_protein_id)
  unique_corrected_only <- unique(corrected_only_df$source_protein_id)
  dataset_summary_all[[dataset]] <- data.frame(
    dataset = dataset,
    corrected_input_directory = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    corrected_mapped_directory = normalizePath(mapped_dir, winslash = "/", mustWork = TRUE),
    input_contrast_csv_count = length(input_files),
    mapped_output_csv_count = length(mapped_files),
    unmapped_output_csv_count = length(unmapped_files),
    member_bridge_output_csv_count = length(bridge_files),
    expected_contrast_count = expected_count,
    protein_rows_per_contrast_min = min(file_validation$input_protein_rows, na.rm = TRUE),
    protein_rows_per_contrast_max = max(file_validation$input_protein_rows, na.rm = TRUE),
    mapped_proteins_per_contrast_min = min(file_validation$mapped_protein_count, na.rm = TRUE),
    mapped_proteins_per_contrast_max = max(file_validation$mapped_protein_count, na.rm = TRUE),
    unmapped_proteins_per_contrast_min = min(file_validation$unmapped_protein_count, na.rm = TRUE),
    unmapped_proteins_per_contrast_max = max(file_validation$unmapped_protein_count, na.rm = TRUE),
    mapping_percentage_min = min(file_validation$mapping_percentage, na.rm = TRUE),
    mapping_percentage_max = max(file_validation$mapping_percentage, na.rm = TRUE),
    duplicate_final_accession_rows_max = max(file_validation$duplicate_final_accession_rows, na.rm = TRUE),
    missing_SYMBOL_count_max = max(file_validation$missing_SYMBOL_count, na.rm = TRUE),
    missing_ENTREZ_count_max = max(file_validation$missing_ENTREZ_count, na.rm = TRUE),
    failed_file_count = sum(file_validation$failed),
    required_sus_vs_res_count = nrow(primary),
    missing_required_sus_vs_res_count = sum(!primary$input_present | !primary$mapped_present),
    mapping_determinism_discrepancy_count = nrow(determinism),
    legacy_only_source_id_instances = nrow(legacy_only_df),
    legacy_only_unique_source_ids = length(unique_legacy_only),
    corrected_only_source_id_instances = nrow(corrected_only_df),
    corrected_only_unique_source_ids = length(unique_corrected_only),
    filenames_and_directions_preserved = all(file_validation$filename_preserved & file_validation$direction_preserved),
    stringsAsFactors = FALSE
  )
}

dataset_summary <- do.call(rbind, dataset_summary_all)
utils::write.csv(dataset_summary, file.path(audit_table_root, "animal_level_mapping_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(do.call(rbind, manifest_all), file.path(audit_log_root, "mapping_audit_source_manifest.csv"), row.names = FALSE, na = "")

report_lines <- c(
  "# Animal-level corrected ID-mapping audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "This report validates only the isolated corrected forward mapping branch. It does not run enrichment.",
  "",
  "| Dataset | Inputs | Mapped outputs | Protein rows | Mapped | Unmapped | Mapping % | Determinism differences | Legacy-only unique IDs | Failed files |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(dataset_summary)), function(i) {
    x <- dataset_summary[i, ]
    sprintf(
      "| %s | %d | %d | %d | %d | %d | %.4f | %d | %d | %d |",
      x$dataset, x$input_contrast_csv_count, x$mapped_output_csv_count,
      x$protein_rows_per_contrast_min, x$mapped_proteins_per_contrast_min,
      x$unmapped_proteins_per_contrast_min, x$mapping_percentage_min,
      x$mapping_determinism_discrepancy_count, x$legacy_only_unique_source_ids,
      x$failed_file_count
    )
  }, character(1)),
  "",
  "Mapping determinism is evaluated only for shared `original_identifier` values in same-named legacy and corrected contrasts. Legacy-only malformed/extra identifiers are reported separately."
)
writeLines(report_lines, file.path(audit_report_root, "animal_level_mapping_summary.md"))
write_session_info(file.path(audit_log_root, "sessionInfo.txt"))

fatal <- dataset_summary$input_contrast_csv_count != dataset_summary$expected_contrast_count |
  dataset_summary$mapped_output_csv_count != dataset_summary$expected_contrast_count |
  dataset_summary$unmapped_output_csv_count != dataset_summary$expected_contrast_count |
  dataset_summary$member_bridge_output_csv_count != dataset_summary$expected_contrast_count |
  dataset_summary$failed_file_count > 0L |
  dataset_summary$missing_required_sus_vs_res_count > 0L |
  dataset_summary$mapping_determinism_discrepancy_count > 0L |
  dataset_summary$corrected_only_unique_source_ids > 0L |
  !dataset_summary$filenames_and_directions_preserved
if (any(fatal)) {
  stop("Animal-level mapping audit failed; inspect ", audit_table_root, call. = FALSE)
}

print(dataset_summary)
cat("Animal-level mapping audit passed.\n")
