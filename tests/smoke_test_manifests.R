#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "validation_utils.R"))

full_manifest_sweep <- tolower(Sys.getenv("PROTEOMICS_MANIFEST_SMOKE_FULL", unset = "")) %in% c("1", "true", "yes")

discover_manifest_files <- function() {
  roots <- if (isTRUE(full_manifest_sweep)) {
    c(path_processed(), path_results())
  } else {
    c(path_processed(), path_results("logs"), path_results("tables"))
  }
  rg <- Sys.which("rg")
  if (nzchar(rg)) {
    files <- system2(rg, c("--files", "--no-ignore", "-g", "*manifest*.csv", roots), stdout = TRUE, stderr = FALSE)
    files <- normalizePath(files, winslash = "/", mustWork = FALSE)
    files <- files[grepl("manifest.*\\.csv$", basename(files), ignore.case = TRUE)]
    files <- files[grepl("/(data/processed|results)/", files)]
    return(files)
  }
  unlist(lapply(roots, list.files, pattern = "manifest.*\\.csv$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
}

manifest_files <- discover_manifest_files()
manifest_files <- unique(manifest_files[file.exists(manifest_files)])
manifest_files <- manifest_files[!grepl("results/logs/pipeline/pipeline_manifest_", manifest_files, fixed = TRUE)]

manifest_limit <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_MANIFEST_SMOKE_LIMIT", unset = "25")))
if (!isTRUE(full_manifest_sweep) && length(manifest_files) > manifest_limit) {
  info <- file.info(manifest_files)
  manifest_files <- rownames(info)[order(info$mtime, decreasing = TRUE)][seq_len(manifest_limit)]
  message(
    "INFO smoke_test_manifests: validating newest ", manifest_limit,
    " manifest CSVs; set PROTEOMICS_MANIFEST_SMOKE_FULL=true for full historical sweep"
  )
}

fail <- character()
if (!length(manifest_files)) {
  message("SKIP smoke_test_manifests: no manifest CSV files found")
  quit(status = 0, save = "no")
}

for (f in manifest_files) {
  df <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    message("WARN could not read manifest: ", f)
    next
  }
  key_cols <- intersect(c("analysis_id", "dataset", "run_id", "ontology", "result_type", "comparison", "contrast", "output_table"), names(df))
  if (length(key_cols) >= 2) {
    dup <- duplicate_key_summary(df, key_cols)
    if (isTRUE(dup$n_duplicate_keys > 0)) fail <- c(fail, paste("Duplicate manifest keys in", f, ":", dup$n_duplicate_keys))
  }
  path_qc <- validate_manifest_paths(df, allow_missing = TRUE)
  missing_required <- path_qc[!path_qc$exists & nzchar(path_qc$path), , drop = FALSE]
  if (nrow(missing_required)) {
    message("WARN manifest has missing paths: ", f, " (", nrow(missing_required), ")")
  }
}

if (length(fail)) {
  message("FAIL smoke_test_manifests")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_manifests")
