# -----------------------------------------------------------
# GCT Comparison Splitter: Per-Comparison Forward and Reverse Output
# -----------------------------------------------------------
# - Reads in GCT-format proteomics data
# - Extracts one CSV per group comparison
# - For each comparison, produces an additional "reverse" CSV
#   - In the reverse CSV, only log2fc/rawlog2fc columns are sign-inverted
#   - Other columns (padj, pval, etc.) are preserved
# - Supports optional label mapping ("con", "res", "sus") vs numeric grouping (_1/_2/_3)
# - Results are in two subfolders: "raw/<subfolder>/" (forward), "raw/<subfolder>/reverse/" (reverse)
# - All steps are commented for maintainability

# ------ Library setup ------
required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "CePa")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Using pacman to install/load missing packages: ", paste(missing_pkgs, collapse = ", "))
  tryCatch(
    pacman::p_load(missing_pkgs, character.only = TRUE),
    error = function(e) warning("pacman failed to install/load: ", conditionMessage(e))
  )
}
invisible(lapply(intersect(required_pkgs, rownames(installed.packages())), function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# ---- Parameters ----
use_label_map <- FALSE  # TRUE for "con"/"res"/"sus", FALSE for "_1"/"_2"/"_3"

# ---- Input ----
setwd("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets")
file <- "data/pg.matrix-neha-new-30oct_Two-sample_mod_T_2025-10-30-transformed-p-val_n120x5349"
gct_data <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct", paste0(file, ".gct"))

# ------- GCT reading and cleaning -------
raw <- tryCatch(
  utils::read.delim(gct_data, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = ""),
  error = function(e) stop("Failed to read GCT file: ", conditionMessage(e))
)
raw_clean <- raw[-c(1:2), ]

header_row <- which(tolower(trimws(as.character(raw_clean[[1]]))) == "id")
if (length(header_row) == 0) {
  header_row_any <- which(apply(raw_clean, 1, function(r) any(tolower(trimws(as.character(r))) == "id")))
  if (length(header_row_any) >= 1) {
    header_row <- header_row_any[1]
    warning("Header row containing 'id' was not found in first column; using first match found in any column at row: ", header_row)
  } else {
    sample_rows <- raw_clean[seq_len(min(10, nrow(raw_clean))), , drop = FALSE]
    stop("Failed to find a header row containing 'id' in the GCT file. First few rows of the cleaned data:\n", paste(capture.output(print(sample_rows)), collapse = "\n"))
  }
} else if (length(header_row) > 1) {
  warning("Multiple rows matching 'id' found in first column; using first occurrence.")
  header_row <- header_row[1]
}
stopifnot(length(header_row) == 1)
colnames <- as.character(unlist(raw_clean[header_row, ]))
data <- raw_clean[(header_row + 1):nrow(raw_clean), ]
colnames(data) <- colnames

# ---- Remove annotation lines ----
annotation_rows <- c(
  "AnimalID", "ReplicateGroup", "celltype", "celltype_layer", "layer", "region",
  "group", "group2", "ExpGroup", "celltype_ExpGroup", "region_ExpGroup",
  "celltype_layer_ExpGroup", "celltype_sublayer_ExpGroup", "plate",
  "sampleNumber",  "shortname"
)
data <- data %>% filter(!id %in% annotation_rows, id != "na")
num_cols <- setdiff(names(data), "id")
data[num_cols] <- lapply(data[num_cols], readr::parse_number)

# ------ Function: extract comparison key ------
extract_comparison <- function(col, use_label_map = TRUE) {
  comp <- stringr::str_extract(col, "[A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+")
  if (!is.na(comp)) return(comp)
  label_map <- c("1" = "con", "2" = "res", "3" = "sus")
  m <- stringr::str_match(col, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
    if (use_label_map) {
      n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
      n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
      return(paste0(r1, n1_label, "_", r2, n2_label))
    } else {
      return(paste0(r1, n1, "_", r2, n2))
    }
  }
  m2 <- stringr::str_match(col, "^([A-Za-z0-9]+)_([0-9]+)_([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m2[1, 1])) {
    r1 <- m2[1, 2]; n1 <- m2[1, 3]; r2 <- m2[1, 4]; n2 <- m2[1, 5]
    if (use_label_map) {
      n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
      n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
      return(paste0(r1, n1_label, "_", r2, n2_label))
    } else {
      return(paste0(r1, n1, "_", r2, n2))
    }
  }
  comp2 <- stringr::str_extract(col, "[A-Za-z0-9_]+_[A-Za-z0-9_]+")
  if (!is.na(comp2)) return(comp2)
  NA_character_
}

# ---- Create mapping of feature columns to comparisons ----
all_cols <- names(data)
feature_cols <- setdiff(all_cols, "id")
comparison_keys <- setNames(
  vapply(feature_cols, extract_comparison, character(1), use_label_map = use_label_map), 
  feature_cols)
valid_cols <- names(comparison_keys)[stringr::str_detect(comparison_keys, "\\.over\\.")]
valid_keys <- comparison_keys[valid_cols]
by_comparison <- split(valid_cols, valid_keys)

# ---- Output folder structure ----
outdir_base <- "raw"
fname <- basename(gct_data)
subfolder <- stringr::str_extract(fname, "(?<=E9-).*?(?=_Two-sample_mod_T)")
if (is.na(subfolder) || !nzchar(subfolder)) subfolder <- stringr::str_extract(fname, "(?<=E9-).*?(?=_Two-sample)")
if (is.na(subfolder) || !nzchar(subfolder)) subfolder <- "unknown-comparison"
outdir <- file.path(outdir_base, subfolder)
fs::dir_create(outdir)
outdir_rev <- file.path(outdir_base, subfolder, "reverse")
fs::dir_create(outdir_rev)

safe_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") %>%
    stringr::str_replace_all("_+", "_")
}

# ------- Main loop: Write files for each comparison -------
purrr::iwalk(by_comparison, function(cols, comp_key) {
  # -- Subset columns for this comparison --
  df_out <- data %>% mutate(id = as.character(id)) %>% select(id, all_of(cols))
  if (is.list(df_out$id)) {
    df_out$id <- vapply(df_out$id, function(x) paste(as.character(x), collapse = ";"), character(1))
  }
  names(df_out)[1] <- "gene_symbol" # always rename id
  # -- Metric column names --
  recode_map <- c(
    "adj.P.Val"   = "padj",
    "P.Value"     = "pval",
    "logFC"       = "log2fc",
    "Log.P.Value" = "logpval",
    "AveExpr"     = "aveExpr",
    "t"           = "t",
    "RawlogFC"    = "rawlog2fc"
  )
  metrics <- stringr::str_remove(cols, "\\.[A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+$")
  new_names <- vapply(metrics, function(m) {
    if (!is.na(m) && m %in% names(recode_map) && use_label_map) recode_map[[m]] else m
  }, character(1))
  new_names <- make.unique(new_names, sep = "_")
  names(df_out)[-1] <- new_names
  # -- Filename for original file --
  comp2 <- comp_key
  m <- stringr::str_match(comp2, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
    label_map <- c("1" = "con", "2" = "res", "3" = "sus")
    if (use_label_map) {
      n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
      n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
    } else {
      n1_label <- n1
      n2_label <- n2
    }
    comp2 <- paste0(r1, n1_label, "_", r2, n2_label)
  } else {
    comp2 <- stringr::str_replace_all(comp2, "\\.over\\.", "_")
    if (use_label_map) {
      comp2 <- stringr::str_replace_all(comp2, "_1", "con")
      comp2 <- stringr::str_replace_all(comp2, "_2", "res")
      comp2 <- stringr::str_replace_all(comp2, "_3", "sus")
    }
  }
  fname <- file.path(outdir, paste0(safe_name(comp2), ".csv"))
  utils::write.csv(df_out, fname, row.names = FALSE, quote = TRUE)
  message("Wrote: ", fname)

  # --- Reverse condition: log fold-changes only ---
  if (!is.na(m[1, 1])) {
    df_rev <- df_out
    # Only invert logFC/RawlogFC columns, keep others unchanged!
    for (col in names(df_rev)) {
      if (grepl("^(logFC|RawlogFC)(?:$|_)?", col, perl = TRUE)) {
        df_rev[[col]] <- suppressWarnings(as.numeric(as.character(df_rev[[col]]))) * -1
      }
    }
    r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
    label_map <- c("1" = "con", "2" = "res", "3" = "sus")
    if (use_label_map) {
      n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
      n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
    } else {
      n1_label <- n1
      n2_label <- n2
    }
    rev_comp2 <- paste0(r2, n2_label, "_", r1, n1_label)
    rev_fname <- file.path(outdir_rev, paste0(safe_name(rev_comp2), ".csv"))
    utils::write.csv(df_rev, rev_fname, row.names = FALSE, quote = TRUE)
    message("Wrote (reversed): ", rev_fname)
  }
})

# ---- Optional: index files (tracking) ----

# Forward index
index <- tibble::tibble(
  comparison = names(by_comparison),
  n_columns = lengths(by_comparison),
  file = vapply(names(by_comparison), function(comp_key) {
    comp2 <- comp_key
    m <- stringr::str_match(comp2, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
    if (!is.na(m[1, 1])) {
      r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
      label_map <- c("1" = "con", "2" = "res", "3" = "sus")
      if (use_label_map) {
        n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
        n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
      } else {
        n1_label <- n1
        n2_label <- n2
      }
      comp2 <- paste0(r1, n1_label, "_", r2, n2_label)
    } else {
      comp2 <- stringr::str_replace_all(comp2, "\\.over\\.", "_")
      if (use_label_map) {
        comp2 <- stringr::str_replace_all(comp2, "_1", "con")
        comp2 <- stringr::str_replace_all(comp2, "_2", "res")
        comp2 <- stringr::str_replace_all(comp2, "_3", "sus")
      }
    }
    file.path(outdir, paste0(safe_name(comp2), ".csv"))
  }, character(1))
)
readr::write_csv(index, file.path(outdir, "indexComparisons.csv"))
message("Wrote index: ", file.path(outdir, "indexComparisons.csv"))

# Reverse index
index_rev <- tibble::tibble(
  comparison = names(by_comparison),
  n_columns = lengths(by_comparison),
  file = vapply(names(by_comparison), function(comp_key) {
    m <- stringr::str_match(comp_key, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
    if (!is.na(m[1, 1])) {
      r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
      label_map <- c("1" = "con", "2" = "res", "3" = "sus")
      if (use_label_map) {
        n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
        n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
      } else {
        n1_label <- n1
        n2_label <- n2
      }
      rev_comp2 <- paste0(r2, n2_label, "_", r1, n1_label)
      file.path(outdir_rev, paste0(safe_name(rev_comp2), ".csv"))
    } else {
      NA_character_
    }
  }, character(1))
)
readr::write_csv(index_rev, file.path(outdir_rev, "indexComparisons_reverse.csv"))
message("Wrote reverse index: ", file.path(outdir_rev, "indexComparisons_reverse.csv"))
