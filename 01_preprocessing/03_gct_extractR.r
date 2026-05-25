# ===========================================================
# Consumes:
#   - ProTigy/limma GCT comparison output from data/processed/01_preprocessing/protigy_output/<comparison>/
# Produces:
#   - forward/reverse raw contrast CSVs for ID mapping under data/processed/01_preprocessing/gct_extractR/<comparison>/
# File contract:
#   - feeds docs/file_contracts.tsv object mapped_contrast_csv through 02_id_mapping/01_MapThatProt_batch.r
# ===========================================================
# GCT Comparison Splitter: Per-Comparison Forward and Reverse Output
# Robust version for Metric.Comparison formatted GCT files
# Handles duplicate column names by keeping first occurrence only
# Author: Tobias Pohl
# ===========================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "gct_extractR"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "tibble")

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "),
         ". Install them explicitly before running this script.", call. = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# -------------------------------
# Parameters
# -------------------------------

use_label_map <- FALSE   # TRUE = con/res/sus mapping
comparison_name <- Sys.getenv("PROTEOMICS_GCT_COMPARISON", unset = "neuron_neuropil")

# -------------------------------
# Input
# -------------------------------

input_stem <- Sys.getenv("PROTEOMICS_GCT_INPUT_STEM", unset = "")
if (!nzchar(input_stem)) {
  manifest_candidates <- c(
    path_processed("01_preprocessing", "protigy_output", comparison_name, "protigy_manifest.csv"),
    path_processed("01_preprocessing", "protigy_output", "protigy_manifest.csv")
  )
  manifest_file <- manifest_candidates[file.exists(manifest_candidates)][1]
  if (!is.na(manifest_file)) {
    manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
    path_col <- intersect(c("gct_path", "file", "path"), names(manifest))[1]
    if (!is.na(path_col)) {
      candidate <- manifest[[path_col]][1]
      input_stem <- tools::file_path_sans_ext(basename(candidate))
    }
  }
}
if (!nzchar(input_stem)) {
  stop(
    "Could not infer ProTigy GCT input stem. Set PROTEOMICS_GCT_INPUT_STEM ",
    "or provide a protigy_manifest.csv with a gct_path/file/path column under ",
    path_processed("01_preprocessing", "protigy_output", comparison_name),
    call. = FALSE
  )
}

gct_path <- path_processed("01_preprocessing", "protigy_output", comparison_name, paste0(input_stem, ".gct"))
outdir <- file.path(CANONICAL_PATHS$processed, comparison_name)

if (is_dry_run()) {
  outdir_fwd <- file.path(outdir, "forward")
  outdir_rev <- file.path(outdir, "reverse")
  dry_run_line("Script", "01_preprocessing/03_gct_extractR.r")
  dry_run_line("GCT input", gct_path, if (file.exists(gct_path)) "PASS" else "FAIL")
  dry_run_line("Forward output directory", outdir_fwd)
  dry_run_line("Reverse output directory", outdir_rev)
  dry_run_line("Index output", file.path(outdir, "indexComparisons.csv"))
  quit(status = if (file.exists(gct_path)) 0 else 1, save = "no")
}
if (!file.exists(gct_path)) stop("GCT input not found: ", gct_path, call. = FALSE)
load_required_packages(required_pkgs)
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))

# -------------------------------
# Helper Functions
# -------------------------------

safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    stringr::str_replace_all("_+", "_")
}

swap_comparison <- function(comp_key, use_label_map = FALSE) {
  parts <- stringr::str_split(comp_key, "\\.over\\.", simplify = TRUE)

  if (ncol(parts) != 2) return(NA_character_)

  rev <- paste0(parts[2], ".over.", parts[1])

  if (use_label_map) {
    rev <- stringr::str_replace_all(rev, "_1", "con")
    rev <- stringr::str_replace_all(rev, "_2", "res")
    rev <- stringr::str_replace_all(rev, "_3", "sus")
  }

  rev
}

split_col <- function(col) {
  m <- stringr::str_match(
    col,
    "^([A-Za-z0-9\\.]+)\\.([A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+)$"
  )

  if (is.na(m[1, 1])) {
    return(list(metric = NA_character_, comparison = NA_character_))
  }

  list(
    metric = m[1, 2],
    comparison = m[1, 3]
  )
}

parse_compkey <- function(key) {
  m <- stringr::str_match(
    key,
    "^([A-Za-z0-9]+)_([a-z]+)_([123])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([123])$"
  )

  label_map <- c("1" = "con", "2" = "res", "3" = "sus")

  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]
    g1 <- m[1, 3]
    l1 <- m[1, 4]

    r2 <- m[1, 5]
    g2 <- m[1, 6]
    l2 <- m[1, 7]

    left  <- paste0(r1, g1, label_map[[l1]])
    right <- paste0(r2, g2, label_map[[l2]])

    return(paste0(left, "_", right))
  }

  key2 <- stringr::str_replace_all(key, "\\.over\\.", "_")
  key2 <- stringr::str_replace_all(key2, "_1", "con")
  key2 <- stringr::str_replace_all(key2, "_2", "res")
  key2 <- stringr::str_replace_all(key2, "_3", "sus")
  key2 <- stringr::str_replace_all(key2, "[^A-Za-z0-9_]", "")

  key2
}

# -------------------------------
# Read GCT File
# -------------------------------

raw <- utils::read.delim(
  gct_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

raw_clean <- raw[-c(1:2), ]

header_row <- which(tolower(trimws(raw_clean[[1]])) == "id")[1]

if (is.na(header_row)) {
  stop("Could not find GCT header row. Expected first column to contain 'id'.")
}

col_names <- as.character(unlist(raw_clean[header_row, ]))

data <- raw_clean[(header_row + 1):nrow(raw_clean), , drop = FALSE]
colnames(data) <- col_names

# -------------------------------
# Remove Duplicate Columns
# Keep first occurrence only
# -------------------------------

dup_cols <- duplicated(names(data))

if (any(dup_cols)) {
  message(
    "Removing duplicate columns, keeping first occurrence: ",
    paste(unique(names(data)[dup_cols]), collapse = ", ")
  )

  data <- data[, !dup_cols, drop = FALSE]
}

if (!"id" %in% names(data)) {
  stop("No 'id' column found after assigning GCT header.")
}

if (anyDuplicated(names(data))) {
  stop(
    "Duplicate columns remain after cleanup: ",
    paste(unique(names(data)[duplicated(names(data))]), collapse = ", ")
  )
}

# -------------------------------
# Remove Annotation Rows
# -------------------------------

annotation_rows <- c(
  "AnimalID", "ReplicateGroup", "celltype", "celltype_layer", "layer", "region",
  "group", "group2", "ExpGroup", "celltype_ExpGroup", "region_ExpGroup",
  "celltype_layer_ExpGroup", "celltype_sublayer_ExpGroup", "plate",
  "sampleNumber", "shortname"
)

data <- data |>
  dplyr::filter(!id %in% annotation_rows, id != "na")

# -------------------------------
# Parse Feature Columns
# -------------------------------

feature_cols <- setdiff(names(data), "id")
split_info <- lapply(feature_cols, split_col)

comparison_keys <- setNames(
  vapply(split_info, `[[`, character(1), "comparison"),
  feature_cols
)

metric_keys <- setNames(
  vapply(split_info, `[[`, character(1), "metric"),
  feature_cols
)

valid_cols <- names(comparison_keys)[!is.na(comparison_keys)]

if (length(valid_cols) == 0) {
  stop("No valid Metric.Comparison columns detected.")
}

by_comparison <- split(valid_cols, comparison_keys[valid_cols])

# -------------------------------
# Convert Numeric Columns
# -------------------------------

data[valid_cols] <- lapply(data[valid_cols], readr::parse_number)

outdir_fwd <- file.path(outdir, "forward")
outdir_rev <- file.path(outdir, "reverse")

fs::dir_create(outdir_fwd)
fs::dir_create(outdir_rev)

# -------------------------------
# Metric Rename Map
# -------------------------------

recode_map <- c(
  "adj.P.Val"   = "padj",
  "P.Value"     = "pval",
  "logFC"       = "log2fc",
  "RawlogFC"    = "rawlog2fc",
  "Log.P.Value" = "logpval",
  "AveExpr"     = "aveExpr",
  "t"           = "t"
)

# -------------------------------
# Main Loop: Write Forward & Reverse CSVs
# -------------------------------

written_index <- purrr::imap_dfr(by_comparison, function(cols, comp_key) {

  df_out <- data |>
    dplyr::select(id, dplyr::all_of(cols)) |>
    dplyr::mutate(id = as.character(id))

  names(df_out)[1] <- "gene_symbol"

  metrics <- metric_keys[cols]

  new_names <- vapply(metrics, function(m) {
    if (m %in% names(recode_map)) {
      recode_map[[m]]
    } else {
      m
    }
  }, character(1))

  new_names <- make.unique(new_names)
  names(df_out)[-1] <- new_names

  comp2 <- parse_compkey(comp_key)

  fwd_file <- file.path(
    outdir_fwd,
    paste0(safe_name(comp2), ".csv")
  )

  utils::write.csv(
    df_out,
    fwd_file,
    row.names = FALSE,
    quote = TRUE
  )

  message("Wrote: ", fwd_file)

  rev_file <- NA_character_
  rev_comp <- NA_character_

  if (stringr::str_detect(comp_key, "\\.over\\.")) {

    df_rev <- df_out

    log_cols <- names(df_rev)[
      stringr::str_detect(
        names(df_rev),
        stringr::regex("log.*fc", ignore_case = TRUE)
      )
    ]

    for (col in log_cols) {
      df_rev[[col]] <- suppressWarnings(as.numeric(df_rev[[col]]) * -1)
    }

    m <- stringr::str_match(
      comp_key,
      "^([A-Za-z0-9]+)_([a-z]+)_([123])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([123])$"
    )

    label_map <- c("1" = "con", "2" = "res", "3" = "sus")

    if (!is.na(m[1, 1])) {
      r1 <- m[1, 2]
      g1 <- m[1, 3]
      l1 <- m[1, 4]

      r2 <- m[1, 5]
      g2 <- m[1, 6]
      l2 <- m[1, 7]

      left  <- paste0(r2, g2, label_map[[l2]])
      right <- paste0(r1, g1, label_map[[l1]])

      rev_comp <- paste0(left, "_", right)

    } else {
      rev_comp <- swap_comparison(comp_key, TRUE)
      rev_comp <- stringr::str_replace_all(rev_comp, "_1", "con")
      rev_comp <- stringr::str_replace_all(rev_comp, "_2", "res")
      rev_comp <- stringr::str_replace_all(rev_comp, "_3", "sus")
      rev_comp <- stringr::str_replace_all(rev_comp, "[^A-Za-z0-9_]", "")
    }

    rev_file <- file.path(
      outdir_rev,
      paste0(safe_name(rev_comp), ".csv")
    )

    utils::write.csv(
      df_rev,
      rev_file,
      row.names = FALSE,
      quote = TRUE
    )

    message("Wrote reversed: ", rev_file)
  }

  tibble::tibble(
    comparison = comp_key,
    parsed_forward_comparison = comp2,
    parsed_reverse_comparison = rev_comp,
    n_columns = length(cols),
    columns_used = paste(cols, collapse = ";"),
    forward_file = fwd_file,
    reverse_file = rev_file
  )
})

# -------------------------------
# Index File
# -------------------------------

readr::write_csv(
  written_index,
  file.path(outdir, "indexComparisons.csv")
)

write_run_manifest(
  file.path(CANONICAL_PATHS$logs, comparison_name, "run_manifest.yml"),
  inputs = list(gct_path = gct_path),
  outputs = list(
    output_dir = outdir,
    forward_dir = outdir_fwd,
    reverse_dir = outdir_rev,
    index = file.path(outdir, "indexComparisons.csv")
  ),
  parameters = list(
    comparison_name = comparison_name,
    input_stem = input_stem,
    use_label_map = use_label_map,
    n_comparisons_exported = nrow(written_index)
  )
)

message("Finished successfully.")
