# Shared path and input-contract helpers for isolated MapThatProt branches.

mapping_is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[\\\\/]|\\\\\\\\|/)", path)
}

resolve_mapping_root <- function(value, default_root) {
  value <- trimws(as.character(value))
  path <- if (!nzchar(value)) {
    default_root
  } else if (mapping_is_absolute_path(value)) {
    path.expand(value)
  } else {
    repo_path(value)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_mapthatprot_roots <- function(
    gct_extract_root = Sys.getenv("PROTEOMICS_GCT_EXTRACT_ROOT", unset = ""),
    mapping_output_root = Sys.getenv("PROTEOMICS_MAPPING_OUTPUT_ROOT", unset = "")) {
  list(
    gct_extract_root = resolve_mapping_root(
      gct_extract_root,
      path_processed("01_preprocessing", "gct_extractR")
    ),
    mapping_output_root = resolve_mapping_root(
      mapping_output_root,
      path_processed("02_id_mapping")
    )
  )
}

mapping_paths_equal <- function(left, right) {
  left <- normalizePath(left, winslash = "/", mustWork = FALSE)
  right <- normalizePath(right, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    left <- tolower(left)
    right <- tolower(right)
  }
  identical(left, right)
}

mapping_results_namespace <- function(mapping_output_root) {
  default_root <- path_processed("02_id_mapping")
  if (mapping_paths_equal(mapping_output_root, default_root)) {
    return("02_id_mapping")
  }
  namespace <- safe_filename(basename(normalizePath(
    mapping_output_root,
    winslash = "/",
    mustWork = FALSE
  )))
  if (!nzchar(namespace)) {
    stop("Could not derive a results namespace from PROTEOMICS_MAPPING_OUTPUT_ROOT.", call. = FALSE)
  }
  namespace
}

resolve_mapthatprot_paths <- function(dataset, direction = "forward", roots = resolve_mapthatprot_roots()) {
  dataset <- validate_dataset(dataset)
  if (!direction %in% c("forward", "reverse")) {
    stop("Mapping direction must be 'forward' or 'reverse'.", call. = FALSE)
  }
  namespace <- mapping_results_namespace(roots$mapping_output_root)
  result_root <- function(kind) path_results(kind, namespace, "MapThatProt_batch")
  list(
    dataset = dataset,
    direction = direction,
    analysis_namespace = namespace,
    gct_extract_root = roots$gct_extract_root,
    mapping_output_root = roots$mapping_output_root,
    raw_dir = file.path(roots$gct_extract_root, dataset, direction),
    mapped_dir = file.path(roots$mapping_output_root, "mapped", dataset, direction, "per_file"),
    unmapped_dir = file.path(roots$mapping_output_root, "unmapped", dataset, direction, "per_file"),
    member_bridge_dir = file.path(roots$mapping_output_root, "member_bridge", dataset, direction, "per_file"),
    tables_root = result_root("tables"),
    logs_root = result_root("logs"),
    reports_root = result_root("reports")
  )
}

animal_level_mapping_expected_counts <- function() {
  c(neuron_neuropil = 30L, neuron_soma = 12L, microglia = 12L)
}

is_animal_level_gct_extract_root <- function(path) {
  identical(
    tolower(basename(normalizePath(path, winslash = "/", mustWork = FALSE))),
    "gct_extractr_animal_level"
  )
}

list_mapthatprot_input_files <- function(raw_dir) {
  sort(list.files(raw_dir, pattern = ".*_.*\\.csv$", full.names = TRUE))
}

validate_mapthatprot_input_count <- function(dataset, direction, input_root, csv_files) {
  dataset <- validate_dataset(dataset)
  corrected <- is_animal_level_gct_extract_root(input_root)
  expected <- if (corrected && identical(direction, "forward")) {
    unname(animal_level_mapping_expected_counts()[[dataset]])
  } else {
    NA_integer_
  }
  observed <- length(csv_files)
  valid <- if (is.na(expected)) observed > 0L else identical(observed, expected)
  if (!valid) {
    expectation <- if (is.na(expected)) "at least one" else as.character(expected)
    stop(
      "MapThatProt input-count validation failed for ", dataset, "/", direction,
      ": expected ", expectation, " CSV file(s), found ", observed,
      " under ", normalizePath(file.path(input_root, dataset, direction), winslash = "/", mustWork = FALSE),
      ".",
      call. = FALSE
    )
  }
  list(corrected = corrected, expected = expected, observed = observed, valid = valid)
}

mapthatprot_resolution_report <- function(paths, input_count_validation) {
  expected <- if (is.na(input_count_validation$expected)) "at least 1" else as.character(input_count_validation$expected)
  c(
    "Resolved dataset" = paths$dataset,
    "Mapping direction" = paths$direction,
    "GCT extract input root" = paths$gct_extract_root,
    "Mapping output root" = paths$mapping_output_root,
    "Results namespace" = paths$analysis_namespace,
    "Raw contrast directory" = paths$raw_dir,
    "Expected raw contrast CSV count" = expected,
    "Raw contrast CSV count" = as.character(input_count_validation$observed),
    "Mapped output directory" = paths$mapped_dir,
    "Unmapped output directory" = paths$unmapped_dir,
    "Member bridge output directory" = paths$member_bridge_dir
  )
}
