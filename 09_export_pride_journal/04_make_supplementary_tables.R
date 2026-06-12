#!/usr/bin/env Rscript
# ================================================================
# Script: 09_export_pride_journal/04_make_supplementary_tables.R
# Stage: export
# Scope: global
# Consumes: required results/tables/; optional results/source_data/.
# Produces: pride_submission/supplementary_tables/.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: PRIDE supplementary tables.
# ================================================================

# Journal/source-data supplementary tables (not PRIDE-required raw/search outputs).
# Copies only config-listed globs by default; no recursive results/ scan.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
datasets <- resolve_export_datasets(cli$dataset, config)
out_dir <- pride_submission_dir("supplementary_tables")
staging_manifest <- file.path(out_dir, "_supplementary_table_staging_manifest.tsv")

globs <- config$supplementary_table_globs
if (length(datasets) < length(config$datasets)) {
  globs <- globs[vapply(globs, function(g) any(vapply(datasets, function(ds) grepl(ds, g, fixed = TRUE), logical(1))) || !grepl("<dataset>", g, fixed = TRUE), logical(1))]
}
config$supplementary_table_globs <- globs
files <- supplementary_candidate_files(config, datasets, include_derived = TRUE)

stable_path_hash <- function(path, n = 8L) {
  codes <- utf8ToInt(gsub("\\\\", "/", relative_to(path, repo_root())))
  h <- 5381
  for (code in codes) {
    h <- (h * 33 + code) %% 4294967296
  }
  high <- floor(h / 65536)
  low <- h %% 65536
  substr(sprintf("%04x%04x", high, low), 1L, n)
}

trim_with_ext <- function(x, max_chars) {
  x <- safe_filename(x, max_chars = max_chars)
  ext <- tools::file_ext(x)
  stem <- tools::file_path_sans_ext(x)
  if (!nzchar(ext) || nchar(x, type = "chars") <= max_chars) return(substr(x, 1L, max_chars))
  ext_part <- paste0(".", ext)
  stem_limit <- max(1L, max_chars - nchar(ext_part, type = "chars"))
  paste0(substr(stem, 1L, stem_limit), ext_part)
}

supplementary_stage_names <- function(files, datasets, max_chars = 96L) {
  if (!length(files)) return(character())
  rel <- vapply(files, relative_to, character(1), root = repo_root())
  parts <- strsplit(gsub("\\\\", "/", rel), "/", fixed = TRUE)
  source_category <- vapply(parts, function(p) {
    if (length(p) >= 2L && identical(p[[1]], "results")) return(paste(p[1:2], collapse = "_"))
    p[[1]] %||% "source"
  }, character(1))
  analysis_folder <- vapply(parts, function(p) {
    hit <- p[grepl("^\\d{2}_", p)]
    if (length(hit)) hit[[1]] else "analysis"
  }, character(1))
  ds <- vapply(files, dataset_for_manifest_file, character(1), datasets = datasets)
  base <- basename(files)

  raw_names <- mapply(
    function(src, analysis, dataset, filename) paste(src, analysis, dataset, filename, sep = "__"),
    source_category,
    analysis_folder,
    ds,
    base,
    USE.NAMES = FALSE
  )
  staged <- vapply(raw_names, trim_with_ext, character(1), max_chars = max_chars)

  duplicate_names <- unique(staged[duplicated(staged) | duplicated(staged, fromLast = TRUE)])
  if (length(duplicate_names)) {
    for (name in duplicate_names) {
      idx <- which(staged == name)
      for (i in idx) {
        ext <- tools::file_ext(staged[[i]])
        stem <- tools::file_path_sans_ext(staged[[i]])
        suffix <- paste0("__", stable_path_hash(files[[i]]))
        ext_part <- if (nzchar(ext)) paste0(".", ext) else ""
        stem_limit <- max(1L, max_chars - nchar(suffix, type = "chars") - nchar(ext_part, type = "chars"))
        staged[[i]] <- paste0(substr(stem, 1L, stem_limit), suffix, ext_part)
      }
    }
  }

  staged
}

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/04_make_supplementary_tables.R")
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Supplementary candidate count", length(files), if (length(files) > 0) "PASS" else "WARN")
  dry_run_line("Supplementary target", out_dir)
  quit(status = 0, save = "no")
}

if (length(files) == 0) {
  warning(
    "No supplementary candidate tables matched config globs. ",
    "Run with --include-derived-results after generating results/tables outputs.",
    call. = FALSE
  )
  empty_status <- data.frame(
    source_path = character(),
    source_dataset = character(),
    staged_file = character(),
    target_path = character(),
    target_path_chars = integer(),
    copied = logical(),
    stringsAsFactors = FALSE
  )
  utils::write.table(empty_status, staging_manifest, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("Copied: 0; failed: 0")
  message("Staging manifest: ", staging_manifest)
} else {
  stage_names <- supplementary_stage_names(files, datasets)
  targets <- file.path(out_dir, stage_names)
  copy_status <- data.frame(
    source_path = vapply(files, relative_to, character(1), root = repo_root()),
    source_dataset = vapply(files, dataset_for_manifest_file, character(1), datasets = datasets),
    staged_file = stage_names,
    target_path = vapply(targets, relative_to, character(1), root = repo_root()),
    target_path_chars = nchar(normalizePath(targets, winslash = "/", mustWork = FALSE), type = "chars"),
    copied = FALSE,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(files)) {
    copy_status$copied[[i]] <- isTRUE(copy_export_file(files[[i]], targets[[i]], dry_run = FALSE))
  }
  utils::write.table(copy_status, staging_manifest, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  copied_count <- sum(copy_status$copied)
  failed_count <- sum(!copy_status$copied)
  message("Staged journal supplementary tables in: ", out_dir)
  message("Copied: ", copied_count, "; failed: ", failed_count)
  message("Staging manifest: ", staging_manifest)
  if (failed_count > 0L) {
    stop("Failed to stage ", failed_count, " supplementary table(s). See: ", staging_manifest, call. = FALSE)
  }
}
