#!/usr/bin/env Rscript
# Creates pride_submission/manifests/pride_file_manifest.tsv from canonical export paths only.
# Use --include-derived-results or --recursive for broader (non-default) scans.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
datasets <- resolve_export_datasets(cli$dataset, config)
export_level <- cli$export_level

collect_canonical_files <- function() {
  files <- character()
  for (ds in datasets) {
    files <- c(files, processed_files_for_dataset(ds, config, include_derived = cli$include_derived_results))
  }
  files <- c(files, pg_matrix_input_paths(config), metadata_input_paths(config))
  staged_meta <- c(
    pride_submission_dir("metadata", "sample_metadata.tsv"),
    pride_submission_dir("metadata", "sdrf_like_metadata.tsv")
  )
  files <- c(files, staged_meta[file.exists(staged_meta)])

  pride_root <- pride_submission_dir()
  if (dir.exists(pride_root)) {
  files <- c(files, list_files_shallow(pride_root, max_depth = 4L))
  }

  if (isTRUE(cli$recursive)) {
    extra_roots <- c(path_processed(), path_results("tables"), path_results("source_data"))
    extra_roots <- extra_roots[dir.exists(extra_roots)]
    files <- c(files, unlist(lapply(extra_roots, list.files, recursive = TRUE, full.names = TRUE, no.. = TRUE)))
  } else if (isTRUE(cli$include_derived_results)) {
    files <- c(files, supplementary_candidate_files(config, datasets, include_derived = TRUE))
  }

  files <- unique(normalizePath(files[file.exists(files)], winslash = "/", mustWork = FALSE))
  files[!file.info(files)$isdir]
}

files <- collect_canonical_files()
manifest_target <- pride_submission_dir("manifests", "pride_file_manifest.tsv")

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/01_make_pride_manifest.R")
  dry_run_line("Export level", export_level)
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Include derived results", cli$include_derived_results)
  dry_run_line("Recursive scan", cli$recursive)
  dry_run_line("Canonical file count", length(files))
  dry_run_line("Manifest target", manifest_target)
  quit(status = 0, save = "no")
}

manifest <- manifest_rows_for_files(files, datasets, config, hash_algo = config$manifest$hash_algorithm %||% "sha256")
manifest$export_level <- export_level
manifest$export_scope <- export_scope_label(export_level)
manifest$file_path <- manifest$relative_path

out_cols <- c(
  "file_path", "relative_path", "dataset", "export_category", "export_level", "export_scope",
  "size_bytes", "modified_time", "sha256",
  "intended_for_PRIDE", "intended_for_supplement"
)
manifest <- manifest[, intersect(out_cols, names(manifest)), drop = FALSE]

utils::write.table(manifest, manifest_target, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote PRIDE/journal manifest (canonical, non-recursive by default): ", manifest_target)
message("Files listed: ", nrow(manifest))
