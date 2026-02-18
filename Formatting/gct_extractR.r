# ===========================================================
# GCT Comparison Splitter: Per-Comparison Forward and Reverse Output
# Robust version for Metric.Comparison formatted GCT files
# Author: Tobias Pohl
# ===========================================================

# -------------------------------
# Library Setup
# -------------------------------
required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "tibble")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
pacman::p_load(char = required_pkgs)

# -------------------------------
# Parameters
# -------------------------------
use_label_map <- FALSE   # TRUE = con/res/sus mapping

# -------------------------------
# Input
# -------------------------------
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets")

input_file <- "20260218-pgmatrix-imputed-neuron-soma-71samples-missing70pct-with-metadata_Two-sample_mod_T_2026-02-18-transformed-p-val_n66x5538"
gct_path <- file.path(
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/protigy_output/neuron_soma",
  paste0(input_file, ".gct")
)

# -------------------------------
# Helper Functions
# -------------------------------

#' Clean up file names for safe output
safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    stringr::str_replace_all("_+", "_")
}

#' Swap comparison key for reverse output
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

#' Split column name into metric and comparison
split_col <- function(col) {
  m <- stringr::str_match(
    col,
    "^([A-Za-z0-9\\.]+)\\.([A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+)$"
  )
  if (is.na(m[1,1])) {
    return(list(metric = NA_character_, comparison = NA_character_))
  }
  list(
    metric = m[1,2],
    comparison = m[1,3]
  )
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
col_names <- as.character(unlist(raw_clean[header_row, ]))
data <- raw_clean[(header_row + 1):nrow(raw_clean), ]
colnames(data) <- col_names

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
by_comparison <- split(valid_cols, comparison_keys[valid_cols])

# -------------------------------
# Convert Numeric Columns
# -------------------------------

data[valid_cols] <- lapply(data[valid_cols], readr::parse_number)

# -------------------------------
# Output Folder Structure
# -------------------------------

outdir_base <- "raw"
fname <- basename(gct_path)
subfolder <- basename(fs::path_dir(gct_path))
if (is.na(subfolder)) subfolder <- "unknown-comparison"

outdir <- file.path(outdir_base, subfolder)
outdir_fwd <- file.path(outdir, "forward")
outdir_rev <- file.path(outdir, "reverse")

fs::dir_create(outdir_fwd)
fs::dir_create(outdir_rev)

# -------------------------------
# Metric Rename Map
# -------------------------------

recode_map <- c(
  "adj.P.Val" = "padj",
  "P.Value"   = "pval",
  "logFC"     = "log2fc",
  "RawlogFC"  = "rawlog2fc",
  "Log.P.Value" = "logpval",
  "AveExpr"   = "aveExpr",
  "t"         = "t"
)

# -------------------------------
# Main Loop: Write Forward & Reverse CSVs
# -------------------------------

purrr::iwalk(by_comparison, function(cols, comp_key) {
  df_out <- data |>
    dplyr::select(id, dplyr::all_of(cols)) |>
    dplyr::mutate(id = as.character(id))
  names(df_out)[1] <- "gene_symbol"

  # Rename metric columns
  metrics <- metric_keys[cols]
  new_names <- vapply(metrics, function(m) {
    if (use_label_map && m %in% names(recode_map)) recode_map[[m]] else m
  }, character(1))
  new_names <- make.unique(new_names)
  names(df_out)[-1] <- new_names

  # Improved: Parse comp_key for output file naming like CA1spres_CA1spcon
  # Accepts e.g. CA1_sp_2.over.CA1_sp_1 or similar
  parse_compkey <- function(key) {
    # Try to match region_group pattern: e.g. CA1_sp_2.over.CA1_sp_1
    m <- stringr::str_match(key, "^([A-Za-z0-9]+)_([a-z]+)_([123])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([123])$")
    label_map <- c("1" = "con", "2" = "res", "3" = "sus")
    if (!is.na(m[1,1])) {
      r1 <- m[1,2]; g1 <- m[1,3]; l1 <- m[1,4];
      r2 <- m[1,5]; g2 <- m[1,6]; l2 <- m[1,7];
      left <- paste0(r1, g1, label_map[[l1]])
      right <- paste0(r2, g2, label_map[[l2]])
      return(paste0(left, "_", right))
    }
    # fallback: just replace .over. with _ and _1/2/3 with con/res/sus
    key2 <- stringr::str_replace_all(key, "\\.over\\.", "_")
    key2 <- stringr::str_replace_all(key2, "_1", "con")
    key2 <- stringr::str_replace_all(key2, "_2", "res")
    key2 <- stringr::str_replace_all(key2, "_3", "sus")
    key2 <- stringr::str_replace_all(key2, "[^A-Za-z0-9_]", "")
    return(key2)
  }

  comp2 <- parse_compkey(comp_key)
  fwd_file <- file.path(outdir_fwd, paste0(safe_name(comp2), ".csv"))
  utils::write.csv(df_out, fwd_file, row.names = FALSE, quote = TRUE)
  message("Wrote: ", fwd_file)

  # Reverse output
  if (stringr::str_detect(comp_key, "\\.over\\.")) {
    df_rev <- df_out
    log_cols <- names(df_rev)[stringr::str_detect(
      names(df_rev),
      stringr::regex("log.*fc", ignore_case = TRUE)
    )]
    for (col in log_cols) {
      df_rev[[col]] <- suppressWarnings(as.numeric(df_rev[[col]]) * -1)
    }
    # Swap left/right for reverse
    m <- stringr::str_match(comp_key, "^([A-Za-z0-9]+)_([a-z]+)_([123])\\.over\\.([A-Za-z0-9]+)_([a-z]+)_([123])$")
    label_map <- c("1" = "con", "2" = "res", "3" = "sus")
    if (!is.na(m[1,1])) {
      r1 <- m[1,2]; g1 <- m[1,3]; l1 <- m[1,4];
      r2 <- m[1,5]; g2 <- m[1,6]; l2 <- m[1,7];
      left <- paste0(r2, g2, label_map[[l2]])
      right <- paste0(r1, g1, label_map[[l1]])
      rev_comp <- paste0(left, "_", right)
    } else {
      rev_comp <- swap_comparison(comp_key, TRUE)
      rev_comp <- stringr::str_replace_all(rev_comp, "_1", "con")
      rev_comp <- stringr::str_replace_all(rev_comp, "_2", "res")
      rev_comp <- stringr::str_replace_all(rev_comp, "_3", "sus")
      rev_comp <- stringr::str_replace_all(rev_comp, "[^A-Za-z0-9_]", "")
    }
    rev_file <- file.path(outdir_rev, paste0(safe_name(rev_comp), ".csv"))
    utils::write.csv(df_rev, rev_file, row.names = FALSE, quote = TRUE)
    message("Wrote (reversed): ", rev_file)
  }
})

# -------------------------------
# Index Files
# -------------------------------

index_tbl <- tibble::tibble(
  comparison = names(by_comparison),
  n_columns = lengths(by_comparison),
  forward_file = file.path(
    outdir_fwd,
    paste0(
      safe_name(
        stringr::str_replace_all(names(by_comparison), "\\.over\\.", "_")
      ),
      ".csv"
    )
  ),
  reverse_file = file.path(
    outdir_rev,
    paste0(
      safe_name(
        vapply(names(by_comparison), swap_comparison, character(1), use_label_map)
      ),
      ".csv"
    )
  )
)

readr::write_csv(index_tbl, file.path(outdir, "indexComparisons.csv"))

message("Finished successfully.")