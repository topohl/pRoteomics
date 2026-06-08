# Helpers for the machine-readable pipeline registry.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("valid_datasets", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}

read_pipeline_registry <- function(path = repo_path("pipeline.yml")) {
  if (!file.exists(path)) {
    stop("Missing pipeline registry: ", path, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read pipeline.yml. Install it with install.packages('yaml') or renv::restore().", call. = FALSE)
  }
  registry <- yaml::read_yaml(path)
  validate_pipeline_registry(registry, path = path)
  registry
}

pipeline_stage_names <- function(registry) {
  names(registry$stages)
}

pipeline_steps <- function(registry, selected_stages, dataset = current_dataset(), include_unsupported = FALSE) {
  datasets <- if (identical(dataset, "all")) valid_datasets() else validate_dataset(dataset)
  out <- do.call(rbind, lapply(selected_stages, function(stage) {
    scripts <- registry$stages[[stage]]$scripts
    do.call(rbind, lapply(scripts, function(step) {
      supported <- unlist(step$datasets %||% registry$datasets, use.names = FALSE)
      scope <- as.character(step$scope %||% if ("global" %in% supported) "global" else "dataset_specific")
      row_datasets <- if (identical(scope, "global") || identical(supported, "global")) "global" else datasets
      rows <- lapply(row_datasets, function(ds) {
        is_supported <- identical(ds, "global") || ds %in% supported
        if (!include_unsupported && !is_supported) return(NULL)
        data.frame(
          dataset = ds,
          stage = stage,
          script = as.character(step$script),
          scope = scope,
          required = isTRUE(step$required),
          supported = is_supported,
          supported_datasets = paste(supported, collapse = ","),
          consumes_required = paste(unlist(step$consumes_required %||% character(), use.names = FALSE), collapse = "|"),
          consumes_optional = paste(unlist(step$consumes_optional %||% character(), use.names = FALSE), collapse = "|"),
          produces = paste(unlist(step$produces %||% character(), use.names = FALSE), collapse = "|"),
          recomputes_core_state = isTRUE(step$recomputes_core_state),
          safe_downstream_rerun = isTRUE(step$safe_downstream_rerun),
          notes = paste(unlist(step$notes %||% character(), use.names = FALSE), collapse = " "),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    }))
  }))
  rownames(out) <- NULL
  out
}

split_registry_paths <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  unlist(strsplit(as.character(x), "\\|", fixed = FALSE), use.names = FALSE)
}

resolve_registry_path <- function(path, dataset = current_dataset()) {
  path <- gsub("<dataset>", dataset, path, fixed = TRUE)
  path <- gsub("\\\\", "/", path)
  if (grepl("^([A-Za-z]:|/|\\\\\\\\)", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(repo_path(path), winslash = "/", mustWork = FALSE)
}

registry_path_exists <- function(path, dataset = current_dataset()) {
  if (!nzchar(path) || grepl(" or |configured|external ", path, ignore.case = TRUE)) return(NA)
  resolved <- resolve_registry_path(path, dataset = dataset)
  if (grepl("[*?\\[]", resolved)) return(length(Sys.glob(resolved)) > 0L)
  file.exists(resolved) || dir.exists(resolved)
}

registry_missing_paths <- function(paths, dataset = current_dataset()) {
  paths <- split_registry_paths(paths)
  exists <- vapply(paths, registry_path_exists, logical(1), dataset = dataset)
  paths[!is.na(exists) & !exists]
}

registry_present_paths <- function(paths, dataset = current_dataset()) {
  paths <- split_registry_paths(paths)
  exists <- vapply(paths, registry_path_exists, logical(1), dataset = dataset)
  paths[!is.na(exists) & exists]
}

validate_pipeline_registry <- function(registry, path = "pipeline.yml") {
  if (is.null(registry$stages) || !length(registry$stages)) {
    stop(path, " has no stages.", call. = FALSE)
  }
  datasets <- unlist(registry$datasets, use.names = FALSE)
  if (!setequal(datasets, valid_datasets())) {
    stop("pipeline.yml datasets must match valid_datasets(): ", paste(valid_datasets(), collapse = ", "), call. = FALSE)
  }
  allowed_datasets <- c(datasets, "global")

  for (stage in names(registry$stages)) {
    scripts <- registry$stages[[stage]]$scripts
    if (is.null(scripts) || !length(scripts)) {
      stop("Stage '", stage, "' has no scripts in pipeline.yml.", call. = FALSE)
    }
    for (step in scripts) {
      if (is.null(step$script) || !nzchar(step$script)) {
        stop("Stage '", stage, "' has an entry without script.", call. = FALSE)
      }
      if (is.null(step$required)) {
        stop("Registry entry missing required flag: ", step$script, call. = FALSE)
      }
      supported <- unlist(step$datasets %||% datasets, use.names = FALSE)
      unknown <- setdiff(supported, allowed_datasets)
      if (length(unknown)) {
        stop("Registry entry ", step$script, " has unsupported datasets: ", paste(unknown, collapse = ", "), call. = FALSE)
      }
      for (field in c("stage", "scope", "produces", "consumes_required", "consumes_optional", "recomputes_core_state", "safe_downstream_rerun", "notes")) {
        if (is.null(step[[field]])) {
          stop("Registry entry ", step$script, " is missing field: ", field, call. = FALSE)
        }
      }
    }
  }
  invisible(TRUE)
}

pipeline_audit_dir <- function() {
  dir_create(path_results("logs", "pipeline"))
}

active_analysis_scripts <- function() {
  roots <- c(
    "01_preprocessing", "02_id_mapping", "03_qc_exploration",
    "04_differential_expression_enrichment", "05_celltype_enrichment_EWCE",
    "06_modules_WGCNA", "07_spatial_networks", "08_behavior_physio_coupling",
    "09_export_pride_journal"
  )
  files <- unlist(lapply(roots, function(root) {
    if (!dir.exists(repo_path(root))) return(character())
    list.files(repo_path(root), pattern = "\\.[Rr]$", recursive = TRUE, full.names = FALSE)
  }), use.names = FALSE)
  files <- file.path(sub("/.*", "", files), files)
  # list.files(full.names = FALSE) above returns paths relative to each root; rebuild safely.
  files <- unlist(lapply(roots, function(root) {
    if (!dir.exists(repo_path(root))) return(character())
    file.path(root, list.files(repo_path(root), pattern = "\\.[Rr]$", recursive = TRUE, full.names = FALSE))
  }), use.names = FALSE)
  gsub("\\\\", "/", files)
}

guess_stage_for_script <- function(script) {
  if (grepl("^01_preprocessing|^02_id_mapping", script)) return("core")
  if (grepl("^03_qc_exploration/04b|^03_qc_exploration/05_empirical", script)) return("qc_global")
  if (grepl("^03_qc_exploration", script)) return("qc")
  if (grepl("^04_differential|^05_celltype", script)) return("enrichment")
  if (grepl("^06_modules_WGCNA/01_WGCNA", script)) return("modules_wgcna")
  if (grepl("^06_modules_WGCNA", script)) return("modules_downstream")
  if (grepl("^07_spatial", script)) return("networks")
  if (grepl("^08_behavior", script)) return("coupling")
  if (grepl("^09_export", script)) return("export")
  "unknown"
}

write_pipeline_validation_tables <- function(registry) {
  audit_dir <- pipeline_audit_dir()
  steps <- pipeline_steps(registry, pipeline_stage_names(registry), dataset = "all", include_unsupported = TRUE)
  missing <- unique(steps$script[!file.exists(repo_path(steps$script))])
  missing_tbl <- data.frame(script = missing, stringsAsFactors = FALSE)
  utils::write.csv(missing_tbl, file.path(audit_dir, "missing_registered_scripts.csv"), row.names = FALSE)

  registered <- unique(c(steps$script, vapply(registry$legacy %||% list(), function(x) as.character(x$script), character(1))))
  unregistered <- setdiff(active_analysis_scripts(), registered)
  unregistered <- unregistered[!grepl("/legacy/", unregistered, fixed = TRUE)]
  if (length(unregistered)) {
    unregistered_tbl <- data.frame(
      script = unregistered,
      guessed_stage = vapply(unregistered, guess_stage_for_script, character(1)),
      reason_not_registered = ifelse(grepl("^[0-9]{2}_.*/[0-9]{2}_|compat|legacy", unregistered), "compatibility wrapper or older numbering", "not present in pipeline.yml"),
      recommendation = "Audit and either add to pipeline.yml or mark as legacy with a replacement/status.",
      stringsAsFactors = FALSE
    )
  } else {
    unregistered_tbl <- data.frame(
      script = character(),
      guessed_stage = character(),
      reason_not_registered = character(),
      recommendation = character(),
      stringsAsFactors = FALSE
    )
  }
  utils::write.csv(unregistered_tbl, file.path(audit_dir, "unregistered_scripts.csv"), row.names = FALSE)
  invisible(list(missing = missing_tbl, unregistered = unregistered_tbl))
}

validate_pipeline_scripts_exist <- function(registry, fail = TRUE) {
  audit <- write_pipeline_validation_tables(registry)
  if (isTRUE(fail) && nrow(audit$missing)) {
    stop("Active script(s) listed in pipeline.yml do not exist:\n", paste(audit$missing$script, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

run_order_script_references <- function(path = repo_path("RUN_ORDER.md")) {
  if (!file.exists(path)) return(character())
  txt <- readLines(path, warn = FALSE)
  refs <- unlist(regmatches(txt, gregexpr("[0-9]{2}_[A-Za-z0-9_./ -]+\\.[Rr]", txt)), use.names = FALSE)
  refs <- gsub("\\\\", "/", trimws(refs))
  refs <- refs[!grepl("^(90_testing|99_deprecated)/", refs)]
  unique(refs[file.exists(repo_path(refs))])
}

validate_run_order_against_registry <- function(registry, run_order_path = repo_path("RUN_ORDER.md")) {
  active <- unique(pipeline_steps(registry, pipeline_stage_names(registry), dataset = "all", include_unsupported = TRUE)$script)
  legacy <- vapply(registry$legacy %||% list(), function(x) as.character(x$script), character(1))
  refs <- run_order_script_references(run_order_path)
  undocumented <- setdiff(refs, c(active, legacy))
  if (length(undocumented)) {
    stop(
      "RUN_ORDER.md references scripts that are neither active in pipeline.yml nor marked legacy:\n",
      paste(undocumented, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
