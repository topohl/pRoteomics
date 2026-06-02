# Shared helpers for 09_export_pride_journal (manifest-driven, pg_matrix-onward).

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("validate_dataset", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}
if (!exists("read_sample_metadata", mode = "function")) {
  source(repo_path("R", "pride_helpers.R"))
}
if (!exists("resolve_dataset_inputs", mode = "function")) {
  source(repo_path("R", "dataset_inputs.R"))
}

export_config_path <- function() {
  repo_path("09_export_pride_journal", "config", "export_config.yml")
}

load_export_config <- function(path = export_config_path()) {
  if (!file.exists(path)) {
    stop("Export config not found: ", path, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read export config.", call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  cfg$config_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  cfg
}

export_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1]]
  }
  has <- function(flag) flag %in% args
  list(
    dataset = value_after("--dataset", default = Sys.getenv("PROTEOMICS_DATASET", unset = "all")),
    export_level = value_after("--export-level", default = "pg_matrix_onward"),
    include_derived_results = has("--include-derived-results"),
    recursive = has("--recursive"),
    dry_run = has("--dry-run") || is_dry_run(),
    skip_validation = has("--skip-validation"),
    skip_supplementary = has("--skip-supplementary"),
    skip_manuscript = has("--skip-manuscript"),
    skip_claims = has("--skip-claims")
  )
}

resolve_export_datasets <- function(dataset_arg, config) {
  if (is.null(dataset_arg) || !nzchar(dataset_arg) || identical(tolower(dataset_arg), "all")) {
    return(config$datasets)
  }
  validate_dataset(dataset_arg, source = "--dataset")
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(unname(openssl::sha256(file(path), algo = "sha256")))
  }
  unname(tools::md5sum(path))
}

list_files_nonrecursive <- function(dir, pattern = NULL) {
  if (!dir.exists(dir)) return(character())
  list.files(dir, pattern = pattern, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
}

list_files_shallow <- function(dir, pattern = NULL, max_depth = 2L) {
  if (!dir.exists(dir)) return(character())
  if (max_depth <= 0L) return(character())
  hits <- list_files_nonrecursive(dir, pattern)
  if (max_depth > 1L) {
    subdirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
    subdirs <- subdirs[basename(subdirs) != basename(dir)]
    deeper <- unlist(lapply(subdirs, list_files_shallow, pattern = pattern, max_depth = max_depth - 1L), use.names = FALSE)
    hits <- unique(c(hits, deeper))
  }
  hits[file.exists(hits) & !file.info(hits)$isdir]
}

glob_paths <- function(patterns, root = repo_root()) {
  patterns <- unique(as.character(patterns))
  patterns <- patterns[nzchar(patterns)]
  if (!length(patterns)) return(character())
  hits <- unlist(lapply(patterns, function(p) {
    full <- if (grepl("^([A-Za-z]:|/|~)", p)) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE) else repo_path(p)
    if (file.exists(full)) return(full)
    Sys.glob(full)
  }), use.names = FALSE)
  hits <- unique(normalizePath(hits[hits != "" & file.exists(hits)], winslash = "/", mustWork = FALSE))
  hits[!file.info(hits)$isdir]
}

pg_matrix_input_paths <- function(config) {
  glob_paths(config$canonical_inputs$pg_matrix$paths)
}

metadata_input_paths <- function(config) {
  glob_paths(config$canonical_inputs$sample_metadata$paths)
}

dataset_matrix_patterns <- function(dataset) {
  dataset <- validate_dataset(dataset)
  c(
    sprintf("^\\d{8}_pgmatrix_imputed_%s_[0-9]+samples_missing70pct(\\.xlsx|_with_metadata\\.xlsx)$", dataset),
    sprintf("^\\d{8}_pgmatrix_imputed_%s_[0-9]+samples_missing70pct_with_metadata\\.gct$", dataset)
  )
}

processed_files_for_dataset <- function(dataset, config, include_derived = FALSE) {
  dataset <- validate_dataset(dataset)
  files <- character()

  prep_root <- repo_path(config$canonical_inputs$processed_preprocessing$root)
  subdirs <- config$canonical_inputs$processed_preprocessing$per_dataset_subdirs
  for (sub in subdirs) {
    base <- file.path(prep_root, sub)
    if (sub %in% c("gct_extractR", "protigy_output")) {
      ds_dir <- file.path(base, dataset)
      if (dir.exists(ds_dir)) {
        files <- c(files, list_files_shallow(ds_dir, max_depth = 3L))
      }
      next
    }
    if (!dir.exists(base)) next
    if (sub %in% c("impute", "excel_convert", "morpheus")) {
      hits <- list.files(
        base,
        pattern = paste0("pgmatrix_imputed_", dataset),
        full.names = TRUE,
        recursive = FALSE,
        ignore.case = TRUE
      )
      files <- c(files, hits)
    }
  }

  map_root <- file.path(repo_path(config$canonical_inputs$processed_id_mapping$root), dataset)
  if (dir.exists(map_root)) {
    files <- c(files, list.files(map_root, pattern = "\\.(csv|tsv|xlsx|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
  }

  de_root <- repo_path(config$canonical_inputs$processed_de_enrichment$root)
  de_ds <- file.path(de_root, c("clusterProfiler", "compareGO"), dataset)
  for (d in de_ds) {
    if (dir.exists(d)) {
      files <- c(files, list.files(d, pattern = "\\.(csv|tsv|xlsx|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
    }
  }

  qc_root <- file.path(repo_path(config$canonical_inputs$qc_reports$root), dataset)
  if (dir.exists(qc_root)) {
    files <- c(files, list.files(qc_root, pattern = "\\.(csv|tsv|md|txt|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
  }

  if (isTRUE(include_derived)) {
    supp_globs <- config$supplementary_table_globs
    supp_globs <- supp_globs[vapply(supp_globs, function(g) grepl(dataset, g, fixed = TRUE) || grepl("\\*", g, fixed = TRUE), logical(1))]
    files <- c(files, glob_paths(supp_globs))
  }

  files <- unique(normalizePath(files[file.exists(files)], winslash = "/", mustWork = FALSE))
  files[!file.info(files)$isdir]
}

classify_export_category <- function(path, config) {
  p <- tolower(gsub("\\\\", "/", path))
  ext <- tolower(tools::file_ext(p))
  if (any(grepl(paste0("/", gsub("\\.", "\\\\.", basename(pg_matrix_input_paths(config)[1])), "$"), p, fixed = FALSE)) ||
      grepl("/pg_matrix/", p)) return("pg_matrix_input")
  if (grepl("/metadata/", p) || grepl("sample_metadata", p)) return("sample_metadata")
  if (grepl("/01_preprocessing/", p)) return("processed_preprocessing")
  if (grepl("/02_id_mapping/", p)) return("protein_id_mapping")
  if (grepl("/04_differential_expression_enrichment/", p)) return("differential_abundance_enrichment")
  if (grepl("/03_qc_exploration/", p) || grepl("/reports/03_qc", p)) return("qc_report")
  if (grepl("/results/tables/", p) || grepl("/results/source_data/", p)) return("derived_results")
  if (grepl("/pride_submission/", p)) return("pride_staging")
  if (grepl("/methods/", p)) return("provenance")
  if (ext %in% c("raw", "wiff", "mzml", "mzxml", "d", "tdf", "baf")) return("external_raw_ms")
  if (grepl("fasta|\\.fa$", p)) return("external_fasta")
  if (grepl("mqpar|search|parameter", p)) return("external_search_parameters")
  "other"
}

export_scope_label <- function(export_level = "pg_matrix_onward") {
  switch(
    export_level,
    pg_matrix_onward = "PRIDE/journal processed-data export from pg_matrix onward",
    .default = paste0("export_level=", export_level)
  )
}

manifest_rows_for_files <- function(files, datasets, config, hash_algo = "sha256") {
  if (!length(files)) {
    return(data.frame(
      relative_path = character(),
      dataset = character(),
      export_category = character(),
      size_bytes = numeric(),
      modified_time = character(),
      sha256 = character(),
      intended_for_PRIDE = logical(),
      intended_for_supplement = logical(),
      stringsAsFactors = FALSE
    ))
  }
  info <- file.info(files)
  rel <- vapply(files, relative_to, character(1), root = repo_root())
  ds <- vapply(files, function(f) {
    for (d in datasets) {
      if (grepl(d, f, fixed = TRUE)) return(d)
    }
  }, character(1))
  ds[is.na(ds)] <- "all"
  cats <- vapply(files, classify_export_category, character(1), config = config)
  hashes <- if (identical(hash_algo, "sha256")) vapply(files, sha256_file, character(1)) else unname(tools::md5sum(files))
  data.frame(
    relative_path = rel,
    dataset = ds,
    export_category = cats,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = hashes,
    intended_for_PRIDE = cats %in% c("pg_matrix_input", "sample_metadata", "processed_preprocessing", "protein_id_mapping", "pride_staging", "provenance"),
    intended_for_supplement = cats %in% c("differential_abundance_enrichment", "derived_results", "qc_report"),
    stringsAsFactors = FALSE
  )
}

copy_export_file <- function(src, dest, dry_run = FALSE) {
  dest_dir <- dirname(dest)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(dry_run)) {
    message("[DRY-RUN] copy ", src, " -> ", dest)
    return(invisible(FALSE))
  }
  src_norm <- normalizePath(src, winslash = "/", mustWork = FALSE)
  dest_norm <- normalizePath(dest, winslash = "/", mustWork = FALSE)
  if (identical(src_norm, dest_norm)) {
    return(invisible(TRUE))
  }
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) warning("Failed to copy: ", src, " -> ", dest, call. = FALSE)
  invisible(ok)
}

supplementary_candidate_files <- function(config, datasets, include_derived = TRUE) {
  globs <- config$supplementary_table_globs
  hits <- glob_paths(globs)
  if (!isTRUE(include_derived)) return(hits)
  hits
}

write_validation_summary_md <- function(checks, out_path, export_level, config) {
  n_fail <- sum(checks$status == "FAIL")
  n_warn <- sum(checks$status == "WARN")
  n_info <- sum(checks$status == "INFO")
  lines <- c(
    "# PRIDE export validation summary",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("Export level: ", export_level),
    paste0("Scope: ", export_scope_label(export_level)),
    "",
    config$scope_note,
    "",
    "## Summary",
    "",
    paste0("- FAIL: ", n_fail),
    paste0("- WARN: ", n_warn),
    paste0("- INFO: ", n_info),
    paste0("- PASS: ", sum(checks$status == "PASS")),
    "",
    "## Checks",
    ""
  )
  for (i in seq_len(nrow(checks))) {
    lines <- c(lines, paste0("- **", checks$check[[i]], "**: ", checks$status[[i]], " — ", checks$detail[[i]]))
  }
  lines <- c(
    lines,
    "",
    "## Reproducibility boundary",
    "",
    "The committed workflow is reproducible from the pg_matrix stage onward.",
    "Missing raw/vendor MS or search-engine files are expected in partial/processed PRIDE depositions.",
    ""
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, out_path)
  invisible(out_path)
}
