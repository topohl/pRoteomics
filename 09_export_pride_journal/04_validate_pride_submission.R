# Validates PRIDE/journal staging folders and manifest presence.
# Consumes pride_submission/; produces pride_submission/validation/validation_report.tsv.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

ensure_pride_dirs()
DRY_RUN <- is_dry_run()
checks <- data.frame(check = character(), status = character(), detail = character(), stringsAsFactors = FALSE)
add_check <- function(name, ok, detail, fail_status = "FAIL") {
  checks <<- rbind(checks, data.frame(check = name, status = if (ok) "PASS" else fail_status, detail = detail, stringsAsFactors = FALSE))
}

manifest <- pride_submission_dir("manifests", "pride_file_manifest.tsv")
add_check("manifest_exists", file.exists(manifest), manifest)
add_check("sample_metadata_exists", length(list.files(pride_submission_dir("metadata"), pattern = "sample_metadata|sdrf", full.names = TRUE)) > 0, pride_submission_dir("metadata"))
add_check("processed_data_dir_exists", dir.exists(pride_submission_dir("processed_data")), pride_submission_dir("processed_data"))
add_check("methods_dir_exists", dir.exists(pride_submission_dir("methods")), pride_submission_dir("methods"))
raw_candidates <- c(
  if (dir.exists(path_raw())) list.files(path_raw(), recursive = TRUE, full.names = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("raw"))) list.files(pride_submission_dir("raw"), recursive = TRUE, full.names = TRUE) else character(0)
)
raw_candidates <- raw_candidates[file.exists(raw_candidates) & !file.info(raw_candidates)$isdir]
search_candidates <- c(
  if (dir.exists(path_external("search_results"))) list.files(path_external("search_results"), recursive = TRUE, full.names = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("search_results"))) list.files(pride_submission_dir("search_results"), recursive = TRUE, full.names = TRUE) else character(0)
)
search_candidates <- search_candidates[file.exists(search_candidates) & !file.info(search_candidates)$isdir]
fasta_candidates <- c(
  if (dir.exists(path_external())) list.files(path_external(), pattern = "\\.(fasta|fa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("methods"))) list.files(pride_submission_dir("methods"), pattern = "\\.(fasta|fa)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character(0)
)
fasta_candidates <- fasta_candidates[file.exists(fasta_candidates) & !file.info(fasta_candidates)$isdir]
parameter_candidates <- c(
  if (dir.exists(path_external())) list.files(path_external(), pattern = "(mqpar|parameter|params|settings).*(xml|txt|json|yaml|yml|csv|tsv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("methods"))) list.files(pride_submission_dir("methods"), pattern = "(mqpar|parameter|params|settings).*(xml|txt|json|yaml|yml|csv|tsv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character(0)
)
parameter_candidates <- parameter_candidates[file.exists(parameter_candidates) & !file.info(parameter_candidates)$isdir]
processed_candidates <- c(
  if (dir.exists(path_processed("01_preprocessing"))) list.files(path_processed("01_preprocessing"), pattern = "\\.(csv|tsv|xlsx|rds|RDS)$", recursive = TRUE, full.names = TRUE) else character(0),
  if (dir.exists(pride_submission_dir("processed_data"))) list.files(pride_submission_dir("processed_data"), pattern = "\\.(csv|tsv|xlsx|rds|RDS)$", recursive = TRUE, full.names = TRUE) else character(0)
)
supplementary_candidates <- if (dir.exists(pride_submission_dir("supplementary_tables"))) {
  list.files(pride_submission_dir("supplementary_tables"), pattern = "\\.(csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE)
} else character(0)
methods_summary <- pride_submission_dir("methods", "methods_summary.md")
add_check("raw_or_vendor_files_present", length(raw_candidates) > 0, paste0(length(raw_candidates), " raw/vendor candidate file(s)"))
add_check("search_engine_outputs_present", length(search_candidates) > 0, paste0(length(search_candidates), " search-engine output candidate file(s)"))
add_check("fasta_database_present", length(fasta_candidates) > 0, paste0(length(fasta_candidates), " FASTA/database candidate file(s)"))
add_check("search_parameters_present", length(parameter_candidates) > 0, paste0(length(parameter_candidates), " search-parameter candidate file(s)"))
add_check("processed_matrices_present", length(processed_candidates) > 0, paste0(length(processed_candidates), " processed data candidate file(s)"))
add_check("supplementary_tables_present", length(supplementary_candidates) > 0, paste0(length(supplementary_candidates), " supplementary candidate table(s)"))
add_check("methods_summary_exists", file.exists(methods_summary), methods_summary)

if (file.exists(manifest)) {
  mf <- utils::read.delim(manifest, stringsAsFactors = FALSE)
  add_check("manifest_has_md5", "md5" %in% names(mf) && all(nzchar(mf$md5)), "md5 column populated")
  add_check("manifest_has_pride_flags", "intended_for_PRIDE" %in% names(mf), "intended_for_PRIDE column")
  add_check("manifest_lists_raw_files", "file_type" %in% names(mf) && any(mf$file_type %in% c("raw_ms", "RAW", "raw")), "raw file type rows in manifest", fail_status = "WARN")
  add_check("manifest_lists_processed_quantification", "file_type" %in% names(mf) && any(mf$file_type %in% c("processed_quantitative_matrix", "RESULT")), "processed quantitative rows in manifest")
}

metadata_files <- list.files(pride_submission_dir("metadata"), pattern = "\\.(tsv|csv)$", full.names = TRUE)
if (length(metadata_files) > 0) {
  meta <- tryCatch(utils::read.delim(metadata_files[[1]], stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  meta_names <- tolower(gsub("[^a-z0-9]+", "_", names(meta)))
  has_raw_file_name <- "raw_file_name" %in% meta_names || "comment_data_file" %in% meta_names
  add_check("sample_metadata_has_raw_file_name", has_raw_file_name, metadata_files[[1]])
  if (has_raw_file_name && length(raw_candidates) > 0) {
    raw_col <- names(meta)[match(TRUE, meta_names %in% c("raw_file_name", "comment_data_file"))]
    listed <- unique(basename(meta[[raw_col]]))
    listed <- listed[!is.na(listed) & nzchar(listed)]
    missing_raw <- setdiff(listed, basename(raw_candidates))
    add_check("sample_metadata_raw_files_traceable", length(listed) > 0 && length(missing_raw) == 0,
              if (length(missing_raw) == 0) paste(length(listed), "listed raw files traceable") else paste(missing_raw, collapse = "; "))
  }
}

if (isTRUE(DRY_RUN)) {
  dry_run_line("Script", "09_export_pride_journal/04_validate_pride_submission.R")
  apply(checks, 1, function(row) dry_run_line(row[["check"]], row[["detail"]], row[["status"]]))
  quit(status = if (any(checks$status == "FAIL")) 1 else 0, save = "no")
}

out_file <- pride_submission_dir("validation", "validation_report.tsv")
utils::write.table(checks, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote validation report: ", out_file)
