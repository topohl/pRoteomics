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

pipeline_steps <- function(registry, selected_stages, dataset = current_dataset()) {
  dataset <- validate_dataset(dataset)
  out <- do.call(rbind, lapply(selected_stages, function(stage) {
    scripts <- registry$stages[[stage]]$scripts
    do.call(rbind, lapply(scripts, function(step) {
      supported <- unlist(step$datasets %||% registry$datasets, use.names = FALSE)
      data.frame(
        stage = stage,
        script = as.character(step$script),
        required = isTRUE(step$required),
        supported = dataset %in% supported,
        supported_datasets = paste(supported, collapse = ","),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(out) <- NULL
  out
}

validate_pipeline_registry <- function(registry, path = "pipeline.yml") {
  if (is.null(registry$stages) || !length(registry$stages)) {
    stop(path, " has no stages.", call. = FALSE)
  }
  datasets <- unlist(registry$datasets, use.names = FALSE)
  if (!setequal(datasets, valid_datasets())) {
    stop("pipeline.yml datasets must match valid_datasets(): ", paste(valid_datasets(), collapse = ", "), call. = FALSE)
  }

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
      unknown <- setdiff(supported, datasets)
      if (length(unknown)) {
        stop("Registry entry ", step$script, " has unsupported datasets: ", paste(unknown, collapse = ", "), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

validate_pipeline_scripts_exist <- function(registry) {
  steps <- pipeline_steps(registry, pipeline_stage_names(registry), dataset = valid_datasets()[[1]])
  missing <- unique(steps$script[!file.exists(repo_path(steps$script))])
  if (length(missing)) {
    stop("Active script(s) listed in pipeline.yml do not exist:\n", paste(missing, collapse = "\n"), call. = FALSE)
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
  active <- unique(pipeline_steps(registry, pipeline_stage_names(registry), dataset = valid_datasets()[[1]])$script)
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
