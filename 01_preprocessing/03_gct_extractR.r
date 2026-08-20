# ================================================================
# Script: 01_preprocessing/03_gct_extractR.r
# Stage: core
# Scope: dataset_specific
# Consumes: PROTEOMICS_GCT_INPUT_ROOT/<dataset>/*.gct (default data/processed/01_preprocessing/protigy_output).
# Produces: PROTEOMICS_GCT_OUTPUT_ROOT/<dataset>/{forward,reverse}/*.csv and indexComparisons.csv
# Notes: Splits physical ProTigy statistical-result GCT fields; statistical fields may be row descriptors.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "protigy_stat_gct_utils.R"))

MODULE_ID <- "01_preprocessing"
DEFAULT_SUBSTEP_ID <- "gct_extractR"
comparison_name <- current_dataset_from_cli()
use_label_map <- FALSE

truthy_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y")
}

has_cli_flag <- function(flags) {
  any(commandArgs(trailingOnly = TRUE) %in% flags) ||
    any(commandArgs(trailingOnly = FALSE) %in% flags)
}

force_rerun <- truthy_env("PROTEOMICS_FORCE_RERUN") ||
  truthy_env("PROTEOMICS_GCT_FORCE_RERUN") ||
  truthy_env("PROTEOMICS_RECOMPUTE") ||
  truthy_env("PROTEOMICS_GCT_RECOMPUTE") ||
  has_cli_flag(c("--force-rerun", "--recompute"))

roots <- resolve_gct_io_roots()
input_stem <- Sys.getenv("PROTEOMICS_GCT_INPUT_STEM", unset = "")
resolved <- resolve_single_protigy_gct(
  roots$input_root,
  comparison_name,
  input_stem = input_stem,
  use_manifest = TRUE
)
gct_path <- resolved$path
outdir <- file.path(roots$output_root, comparison_name)

default_output_root <- normalizePath(
  path_processed("01_preprocessing", "gct_extractR"),
  winslash = "/", mustWork = FALSE
)
analysis_namespace <- if (identical(roots$output_root, default_output_root)) {
  DEFAULT_SUBSTEP_ID
} else {
  safe_filename(basename(roots$output_root))
}
CANONICAL_PATHS <- module_paths(MODULE_ID, analysis_namespace)

animal_level_mode <- truthy_env("PROTEOMICS_GCT_STRICT_ANIMAL_LEVEL") ||
  identical(tolower(basename(roots$input_root)), "protigy_output_animal_level")
gct <- read_protigy_stat_gct(gct_path, comparison_name, strict_primary = animal_level_mode)
corrected_contract <- if (animal_level_mode) validate_corrected_protigy_gct_contract(gct) else NULL
diagnostic <- protigy_gct_diagnostic_row(gct)

if (is_dry_run()) {
  outdir_fwd <- file.path(outdir, "forward")
  outdir_rev <- file.path(outdir, "reverse")
  existing_tables <- c(
    list.files(outdir_fwd, pattern = "\\.csv$", full.names = TRUE),
    list.files(outdir_rev, pattern = "\\.csv$", full.names = TRUE)
  )
  dry_run_line("Script", "01_preprocessing/03_gct_extractR.r")
  dry_run_line("Dry-run mode", "no output folders or CSV files will be created")
  dry_run_line("Dataset", comparison_name)
  dry_run_line("Analysis/input root", roots$input_root, if (dir.exists(roots$input_root)) "PASS" else "FAIL")
  dry_run_line("Resolved dataset directory", resolved$dataset_dir, if (dir.exists(resolved$dataset_dir)) "PASS" else "FAIL")
  dry_run_line("Resolved GCT", gct_path, if (file.exists(gct_path)) "PASS" else "FAIL")
  dry_run_line(
    "GCT dimensions",
    paste0(
      "nrmat=", diagnostic$nrmat, ", ncmat=", diagnostic$ncmat,
      ", nrhd=", diagnostic$nrhd, ", nchd=", diagnostic$nchd,
      ", physical_fields=", diagnostic$observed_physical_fields
    )
  )
  dry_run_line("Detected comparison naming style", diagnostic$comparison_naming_style)
  dry_run_line("Detected comparisons", diagnostic$n_detected_comparisons)
  dry_run_line("Valid within-unit stress comparisons", diagnostic$n_valid_stress_comparisons)
  dry_run_line("Rejected comparisons", diagnostic$n_detected_comparisons - diagnostic$n_valid_stress_comparisons)
  dry_run_line("Available statistical metrics", diagnostic$metrics_found)
  dry_run_line("Strict animal-level mode", animal_level_mode)
  if (animal_level_mode) dry_run_line("Expected corrected comparisons", corrected_contract$expected_comparison_count, "PASS")
  dry_run_line("Analysis/output root", roots$output_root)
  dry_run_line("Forward output directory", outdir_fwd)
  dry_run_line("Reverse output directory", outdir_rev)
  dry_run_line("Index output", file.path(outdir, "indexComparisons.csv"))
  dry_run_line("Existing split contrast tables", length(existing_tables))
  dry_run_line("Recompute existing tables", force_rerun)
  quit(status = 0, save = "no")
}

required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

dir_create(CANONICAL_PATHS$logs)
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))

safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    stringr::str_replace_all("_+", "_")
}

swap_comparison <- function(comp_key, use_label_map = FALSE) {
  parts <- stringr::str_split(canonicalize_protigy_comparison(comp_key), "\\.over\\.", simplify = TRUE)
  if (ncol(parts) != 2L) return(NA_character_)
  rev <- paste0(parts[2], ".over.", parts[1])
  if (use_label_map) {
    rev <- stringr::str_replace_all(rev, "_1", "con")
    rev <- stringr::str_replace_all(rev, "_2", "res")
    rev <- stringr::str_replace_all(rev, "_3", "sus")
  }
  rev
}

parse_compkey <- function(key) {
  parsed <- parse_protigy_comparison(key, comparison_name)
  label_map <- c("1" = "con", "2" = "res", "3" = "sus")
  if (isTRUE(parsed$parsed[[1]])) {
    left <- paste0(gsub("_", "", parsed$left_unit[[1]], fixed = TRUE), label_map[[parsed$left_group[[1]]]])
    right <- paste0(gsub("_", "", parsed$right_unit[[1]], fixed = TRUE), label_map[[parsed$right_group[[1]]]])
    return(paste0(left, "_", right))
  }
  key2 <- stringr::str_replace_all(canonicalize_protigy_comparison(key), "\\.over\\.", "_")
  key2 <- stringr::str_replace_all(key2, "_1", "con")
  key2 <- stringr::str_replace_all(key2, "_2", "res")
  key2 <- stringr::str_replace_all(key2, "_3", "sus")
  stringr::str_replace_all(key2, "[^A-Za-z0-9_]", "")
}

parsed_fields <- gct$parsed_fields
parsed_fields <- parsed_fields[!duplicated(parsed_fields$metric_comparison_key), , drop = FALSE]
if (animal_level_mode) {
  accepted <- gct$comparison_validation$comparison[gct$comparison_validation$valid_stress_comparison]
  parsed_fields <- parsed_fields[parsed_fields$comparison %in% accepted, , drop = FALSE]
}
if (!nrow(parsed_fields)) stop("No valid Metric.Comparison fields detected.", call. = FALSE)
by_comparison <- split(parsed_fields, parsed_fields$comparison)

outdir_fwd <- file.path(outdir, "forward")
outdir_rev <- file.path(outdir, "reverse")
fs::dir_create(outdir_fwd)
fs::dir_create(outdir_rev)

recode_map <- c(
  "adj.P.Val" = "padj",
  "P.Value" = "pval",
  "logFC" = "log2fc",
  "RawlogFC" = "rawlog2fc",
  "Log.P.Value" = "logpval",
  "AveExpr" = "aveExpr",
  "RawAveExpr" = "rawAveExpr",
  "t" = "t"
)

written_index <- purrr::imap_dfr(by_comparison, function(field_rows, comp_key) {
  df_out <- data.frame(gene_symbol = as.character(gct$ids), stringsAsFactors = FALSE)
  metrics <- field_rows$metric
  output_names <- vapply(metrics, function(metric) {
    if (metric %in% names(recode_map)) recode_map[[metric]] else metric
  }, character(1))
  output_names <- make.unique(output_names)

  for (i in seq_len(nrow(field_rows))) {
    raw_values <- as.character(gct$data[[field_rows$column_internal[[i]]]])
    df_out[[output_names[[i]]]] <- if (identical(metrics[[i]], "significant") && animal_level_mode) {
      toupper(trimws(raw_values)) == "TRUE"
    } else {
      readr::parse_number(raw_values, na = c("", "NA", "NaN"))
    }
  }

  comp2 <- parse_compkey(comp_key)
  fwd_file <- file.path(outdir_fwd, paste0(safe_name(comp2), ".csv"))
  fwd_status <- "computed"
  if (file.exists(fwd_file) && !isTRUE(force_rerun)) {
    fwd_status <- "skipped_existing"
    message("Skipping existing table: ", fwd_file)
  } else {
    utils::write.csv(df_out, fwd_file, row.names = FALSE, quote = TRUE)
    message("Wrote: ", fwd_file)
  }

  df_rev <- reverse_protigy_metric_frame(df_out, stats::setNames(metrics, output_names))
  reversed_key <- swap_comparison(comp_key)
  rev_comp <- parse_compkey(reversed_key)
  rev_file <- file.path(outdir_rev, paste0(safe_name(rev_comp), ".csv"))
  rev_status <- "computed"
  if (file.exists(rev_file) && !isTRUE(force_rerun)) {
    rev_status <- "skipped_existing"
    message("Skipping existing reversed table: ", rev_file)
  } else {
    utils::write.csv(df_rev, rev_file, row.names = FALSE, quote = TRUE)
    message("Wrote reversed: ", rev_file)
  }

  tibble::tibble(
    comparison = comp_key,
    parsed_forward_comparison = comp2,
    parsed_reverse_comparison = rev_comp,
    n_columns = nrow(field_rows),
    columns_used = paste(field_rows$field, collapse = ";"),
    forward_file = fwd_file,
    reverse_file = rev_file,
    forward_status = fwd_status,
    reverse_status = rev_status
  )
})

readr::write_csv(written_index, file.path(outdir, "indexComparisons.csv"))

write_run_manifest(
  file.path(CANONICAL_PATHS$logs, comparison_name, "run_manifest.yml"),
  inputs = list(
    analysis_input_root = roots$input_root,
    resolved_dataset_directory = resolved$dataset_dir,
    gct_path = gct_path,
    gct_sha256 = gct$sha256
  ),
  outputs = list(
    analysis_output_root = roots$output_root,
    output_dir = outdir,
    forward_dir = outdir_fwd,
    reverse_dir = outdir_rev,
    index = file.path(outdir, "indexComparisons.csv")
  ),
  parameters = list(
    comparison_name = comparison_name,
    input_stem = resolved$input_stem,
    input_resolution = resolved$resolution,
    comparison_naming_style = gct$naming_style,
    strict_animal_level_mode = animal_level_mode,
    use_label_map = use_label_map,
    force_rerun = force_rerun,
    n_detected_comparisons = nrow(gct$comparison_validation),
    n_valid_stress_comparisons = sum(gct$comparison_validation$valid_stress_comparison),
    n_rejected_comparisons = sum(!gct$comparison_validation$valid_stress_comparison),
    n_comparisons_exported = nrow(written_index)
  )
)

message("Finished successfully.")
