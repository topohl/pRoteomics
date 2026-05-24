# Validate the local PRIDE/ProteomeXchange package skeleton.
#
# The validation report is intentionally explicit. FAIL means the package is
# probably not ready for PRIDE upload. WARN means manual review is needed.

source("R/paths.R")
source("R/pride_helpers.R")

ensure_pride_dirs()
package_dir <- pride_package_dir()
metadata_dir <- repo_path("PRIDE_package", "00_metadata")
out_file <- file.path(metadata_dir, "validation_report.tsv")

add_check <- function(checks, check_name, status, detail) {
  rbind(checks, data.frame(
    check_name = check_name,
    status = status,
    detail = detail,
    stringsAsFactors = FALSE
  ))
}

checks <- data.frame(
  check_name = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

files <- list_pride_files(package_dir)
rel <- relative_to(files, package_dir)

sdrf_path <- file.path(metadata_dir, "sdrf_proteomics.tsv")
manifest_path <- file.path(metadata_dir, "pride_file_manifest.tsv")
checksum_path <- file.path(metadata_dir, "checksum.md5")

checks <- add_check(checks, "sdrf_present", ifelse(file.exists(sdrf_path), "PASS", "FAIL"), sdrf_path)
checks <- add_check(checks, "manifest_present", ifelse(file.exists(manifest_path), "PASS", "FAIL"), manifest_path)
checks <- add_check(checks, "checksum_present", ifelse(file.exists(checksum_path), "PASS", "FAIL"), checksum_path)
if (file.exists(checksum_path)) {
  checksum_lines <- readLines(checksum_path, warn = FALSE)
  checks <- add_check(
    checks,
    "md5_checksums_nonempty",
    ifelse(length(checksum_lines) > 0 && any(nzchar(checksum_lines)), "PASS", "FAIL"),
    paste(length(checksum_lines), "checksum line(s)")
  )
}

raw_files <- files[vapply(rel, classify_pride_file, character(1)) == "RAW"]
checks <- add_check(
  checks,
  "raw_files_present",
  ifelse(length(raw_files) > 0, "PASS", "FAIL"),
  paste(length(raw_files), "RAW-like files found in PRIDE_package")
)

search_or_result <- files[vapply(rel, classify_pride_file, character(1)) %in% c("SEARCH", "RESULT")]
search_files <- files[vapply(rel, classify_pride_file, character(1)) == "SEARCH"]
result_files <- files[vapply(rel, classify_pride_file, character(1)) == "RESULT"]
fasta_files <- files[vapply(rel, classify_pride_file, character(1)) == "FASTA"]
parameter_files <- files[vapply(rel, classify_pride_file, character(1)) == "PARAMETERS_FILE"]
checks <- add_check(
  checks,
  "search_or_result_files_present",
  ifelse(length(search_or_result) > 0, "PASS", "FAIL"),
  paste(length(search_or_result), "SEARCH/RESULT-like files found in PRIDE_package")
)
checks <- add_check(
  checks,
  "search_engine_outputs_present",
  ifelse(length(search_files) > 0, "PASS", "FAIL"),
  paste(length(search_files), "SEARCH-like files found in PRIDE_package")
)
checks <- add_check(
  checks,
  "processed_quantification_or_results_present",
  ifelse(length(result_files) > 0, "PASS", "FAIL"),
  paste(length(result_files), "RESULT-like files found in PRIDE_package")
)
checks <- add_check(
  checks,
  "fasta_database_present",
  ifelse(length(fasta_files) > 0, "PASS", "FAIL"),
  paste(length(fasta_files), "FASTA/database-like files found in PRIDE_package")
)
checks <- add_check(
  checks,
  "search_parameters_present",
  ifelse(length(parameter_files) > 0, "PASS", "FAIL"),
  paste(length(parameter_files), "PARAMETERS_FILE-like files found in PRIDE_package")
)

duplicated_names <- unique(basename(files)[duplicated(basename(files))])
checks <- add_check(
  checks,
  "duplicate_file_names",
  ifelse(length(duplicated_names) == 0, "PASS", "WARN"),
  ifelse(length(duplicated_names) == 0, "No duplicate basenames", paste(duplicated_names, collapse = "; "))
)

meta <- read_sample_metadata()
if (nrow(meta) == 0) {
  checks <- add_check(checks, "sample_metadata_found", "FAIL", "No sample_metadata_clean.tsv/csv found in expected locations")
} else {
  checks <- add_check(checks, "sample_metadata_found", "PASS", paste(nrow(meta), "rows"))
  required <- c("sample_id", "animal_id", "sex", "group", "region", "layer", "raw_file_name")
  present <- intersect(required, names(meta))
  missing <- setdiff(required, names(meta))
  checks <- add_check(
    checks,
    "sample_metadata_required_columns",
    ifelse(length(missing) == 0, "PASS", "WARN"),
    paste("present:", paste(present, collapse = ", "), "| missing:", paste(missing, collapse = ", "))
  )

  for (col in intersect(c("sex", "group", "region", "layer"), names(meta))) {
    n_missing <- sum(is.na(meta[[col]]) | meta[[col]] == "")
    checks <- add_check(
      checks,
      paste0("missing_metadata_", col),
      ifelse(n_missing == 0, "PASS", "WARN"),
      paste(n_missing, "missing values")
    )
  }

  raw_col <- intersect(c("raw_file_name", "raw_file", "data_file", "filename", "file_name"), names(meta))
  if (length(raw_col) == 0) {
    checks <- add_check(checks, "raw_file_mapping", "FAIL", "No raw_file_name/raw_file/data_file column found")
  } else {
    raw_values <- unique(meta[[raw_col[[1]]]])
    raw_values <- raw_values[!is.na(raw_values) & raw_values != ""]
    present_raw_basenames <- basename(raw_files)
    missing_raw <- setdiff(basename(raw_values), present_raw_basenames)
    checks <- add_check(
      checks,
      "raw_files_listed_in_metadata_present_in_package",
      ifelse(length(missing_raw) == 0 && length(raw_values) > 0, "PASS", "FAIL"),
      ifelse(length(missing_raw) == 0, paste(length(raw_values), "raw files mapped"), paste(missing_raw, collapse = "; "))
    )
  }
}

write_tsv(checks, out_file)
message("Wrote validation report: ", out_file)

if (any(checks$status == "FAIL")) {
  message("Validation contains FAIL entries. Review validation_report.tsv before PRIDE upload.")
}
