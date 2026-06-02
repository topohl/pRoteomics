#!/usr/bin/env Rscript
# Contract-aware validation for pg_matrix_onward PRIDE/journal staging.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
export_level <- cli$export_level
rules <- config$validation_rules[[export_level]]
if (is.null(rules)) {
  stop("No validation_rules for export_level: ", export_level, call. = FALSE)
}

checks <- data.frame(check = character(), status = character(), detail = character(), stringsAsFactors = FALSE)
add_check <- function(name, ok, detail, severity = "FAIL") {
  status <- if (isTRUE(ok)) "PASS" else severity
  checks <<- rbind(checks, data.frame(check = name, status = status, detail = detail, stringsAsFactors = FALSE))
}

severity_for <- function(rule_name, default = "FAIL") {
  rules[[rule_name]] %||% default
}

pg_present <- length(pg_matrix_input_paths(config)) > 0 && any(file.exists(pg_matrix_input_paths(config)))
add_check("pg_matrix_present", pg_present, paste(pg_matrix_input_paths(config), collapse = "; "), severity_for("pg_matrix_present"))

meta_files <- c(
  pride_submission_dir("metadata", "sample_metadata.tsv"),
  pride_submission_dir("metadata", "sdrf_like_metadata.tsv")
)
add_check(
  "sample_metadata_present",
  any(file.exists(meta_files)),
  paste(meta_files, collapse = "; "),
  severity_for("sample_metadata_present")
)
add_check(
  "sdrf_like_metadata_present",
  file.exists(meta_files[[2]]),
  meta_files[[2]],
  severity_for("sdrf_like_metadata_present")
)

processed_pkg <- pride_submission_dir("processed_data")
processed_files <- if (dir.exists(processed_pkg)) {
  list.files(processed_pkg, recursive = TRUE, full.names = TRUE)
} else character(0)
processed_files <- processed_files[file.exists(processed_files) & !file.info(processed_files)$isdir]
add_check(
  "processed_matrix_package_present",
  length(processed_files) > 0,
  paste0(length(processed_files), " file(s) under ", processed_pkg),
  severity_for("processed_matrix_package_present")
)

id_map_hits <- unlist(lapply(config$datasets, function(ds) {
  root <- file.path(repo_path(config$canonical_inputs$processed_id_mapping$root), ds, "forward", "per_file")
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = "\\.csv$", full.names = TRUE)
}))
add_check(
  "id_mapping_present",
  length(id_map_hits) > 0,
  paste0(length(id_map_hits), " mapped contrast CSV(s)"),
  severity_for("id_mapping_present")
)

raw_candidates <- c(
  if (dir.exists(path_raw())) list.files(path_raw(), recursive = TRUE, full.names = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("raw"))) list.files(pride_submission_dir("raw"), recursive = TRUE, full.names = TRUE) else character(0)
)
raw_candidates <- raw_candidates[file.exists(raw_candidates) & !file.info(raw_candidates)$isdir]
add_check(
  "raw_ms_vendor_present",
  length(raw_candidates) > 0,
  paste0(length(raw_candidates), " raw/vendor candidate file(s); external to repo entry point"),
  severity_for("raw_ms_vendor_present", "INFO")
)

search_root <- path_external("search_results")
search_candidates <- if (dir.exists(search_root)) list.files(search_root, recursive = TRUE, full.names = TRUE) else character(0)
add_check(
  "search_engine_outputs_present",
  length(search_candidates) > 0,
  paste0(length(search_candidates), " search-engine output candidate file(s); external"),
  severity_for("search_engine_outputs_present", "INFO")
)

fasta_candidates <- if (dir.exists(path_external())) {
  list.files(path_external(), pattern = "\\.(fasta|fa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
} else character(0)
add_check(
  "fasta_database_present",
  length(fasta_candidates) > 0,
  paste0(length(fasta_candidates), " FASTA candidate file(s); external"),
  severity_for("fasta_database_present", "WARN")
)

parameter_candidates <- if (dir.exists(path_external())) {
  list.files(path_external(), pattern = "(mqpar|parameter|params|settings)", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
} else character(0)
add_check(
  "search_parameters_present",
  length(parameter_candidates) > 0,
  paste0(length(parameter_candidates), " search-parameter candidate file(s); external"),
  severity_for("search_parameters_present", "WARN")
)

supp_dir <- pride_submission_dir("supplementary_tables")
supp_files <- if (dir.exists(supp_dir)) list.files(supp_dir, recursive = TRUE, full.names = TRUE) else character(0)
add_check(
  "supplementary_tables_present",
  length(supp_files) > 0,
  paste0(length(supp_files), " supplementary table(s)"),
  severity_for("supplementary_tables_present", "WARN")
)

manifest <- pride_submission_dir("manifests", "pride_file_manifest.tsv")
add_check("pride_manifest_present", file.exists(manifest), manifest, severity_for("pride_manifest_present", "WARN"))

methods_summary <- pride_submission_dir("methods", "methods_summary.md")
add_check("methods_summary_present", file.exists(methods_summary), methods_summary, severity_for("methods_summary_present", "WARN"))

claims <- path_results("tables", "biological_claims_table.csv")
add_check(
  "biological_claims_present",
  file.exists(claims),
  claims,
  severity_for("biological_claims_present", "INFO")
)

if (file.exists(manifest)) {
  mf <- utils::read.delim(manifest, stringsAsFactors = FALSE)
  add_check("manifest_has_sha256", "sha256" %in% names(mf) && all(nzchar(mf$sha256)), "sha256 column populated", "WARN")
}

report_tsv <- pride_submission_dir("validation", "validation_report.tsv")
summary_md <- pride_submission_dir("validation", "validation_summary.md")

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/05_validate_pride_submission.R")
  dry_run_line("Export level", export_level)
  apply(checks, 1, function(row) dry_run_line(row[["check"]], row[["detail"]], row[["status"]]))
  quit(status = 0, save = "no")
}

utils::write.table(checks, report_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write_validation_summary_md(checks, summary_md, export_level, config)
message("Wrote validation report: ", report_tsv)
message("Wrote validation summary: ", summary_md)
if (any(checks$status == "FAIL")) {
  stop("PRIDE export validation failed. See ", report_tsv, call. = FALSE)
}
