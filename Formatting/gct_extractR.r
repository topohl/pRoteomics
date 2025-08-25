# Install pacman if needed, then use pacman to install/load packages (only if not already installed)
required_pkgs <- c("dplyr", "readr", "stringr", "purrr", "fs", "CePa")

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# find missing packages
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  message("Using pacman to install/load missing packages: ", paste(missing_pkgs, collapse = ", "))
  tryCatch(
    pacman::p_load(missing_pkgs, character.only = TRUE),
    error = function(e) warning("pacman failed to install/load: ", conditionMessage(e))
  )
}

# Ensure all requested packages are loaded (those already installed)
invisible(lapply(intersect(required_pkgs, rownames(installed.packages())), function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# --------- READ AND CLEAN GCT (same idea as before) ----------
# Adjust filename
# read in files from direction
# set wd
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets")

# read in .gct
gct_data <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/E9-region-expgroup_Two-sample_mod_T_2025-08-22_n323x5458.gct"

raw <- tryCatch(
  utils::read.delim(gct_data, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = ""),
  error = function(e) stop("Failed to read GCT file: ", conditionMessage(e))
)

# Remove initial header lines; adjust if needed
raw_clean <- raw[-c(1:2), ]

# Find header row (contains "id") using a robust, case-insensitive trimmed match.
# Prefer first-column match but fall back to searching any column and provide diagnostics if not found.
header_row <- which(tolower(trimws(as.character(raw_clean[[1]]))) == "id")

if (length(header_row) == 0) {
  # try to find 'id' anywhere in the row (any column)
  header_row_any <- which(apply(raw_clean, 1, function(r) any(tolower(trimws(as.character(r))) == "id")))
  if (length(header_row_any) >= 1) {
    header_row <- header_row_any[1]
    warning("Header row containing 'id' was not found in first column; using first match found in any column at row: ", header_row)
  } else {
    # diagnostic: show first few rows to help debug the file format
    sample_rows <- raw_clean[seq_len(min(10, nrow(raw_clean))), , drop = FALSE]
    stop("Failed to find a header row containing 'id' in the GCT file. First few rows of the cleaned data:\n",
         paste(capture.output(print(sample_rows)), collapse = "\n"))
  }
} else if (length(header_row) > 1) {
  warning("Multiple rows matching 'id' found in the first column; using the first occurrence.")
  header_row <- header_row[1]
}

stopifnot(length(header_row) == 1)

colnames <- as.character(unlist(raw_clean[header_row, ]))

# Extract data rows after header
data <- raw_clean[(header_row + 1):nrow(raw_clean), ]
colnames(data) <- colnames

# Remove annotation rows and "na"
annotation_rows <- c(
  "AnimalID", "ReplicateGroup", "celltype", "celltype_layer", "layer", "region",
  "group", "group2", "ExpGroup", "celltype_ExpGroup", "region_ExpGroup",
  "celltype_layer_ExpGroup", "celltype_sublayer_ExpGroup", "plate",
  "sampleNumber",  "shortname"
)

data <- data %>%
  filter(!id %in% annotation_rows, id != "na")

# Convert numerics (keep id as character)
num_cols <- setdiff(names(data), "id")
data[num_cols] <- lapply(data[num_cols], readr::parse_number)

# --------- PARSE COMPARISON KEYS AND SPLIT ----------
# We assume feature columns follow the pattern: <metric>.<comparison>
# Examples:
#   adj.P.Val.CA1_2.over.CA1_1
#   P.Value.CA1_3.over.CA1_1
#   logFC.DG_2.over.CA1_1
#   AveExpr.CA3_2.over.CA1_1
#   t.DG_1.over.CA1_1
#   RawlogFC.CA3_1.over.CA1_1
#
# We'll extract the comparison part as the substring after the first dot:
#   metric = adj.P.Val
#   comparison = CA1_2.over.CA1_1
#
# Note: columns named exactly "id" will be skipped from parsing.

all_cols <- names(data)
feature_cols <- setdiff(all_cols, "id")

# Helper: extract the comparison key from a feature column name.
# Behavior:
# - Prefer explicit ".over." tokens (e.g. "adj.P.Val.CA1_2.over.CA1_1" ->
#   "CA1_2.over.CA1_1")
# - If not present, try to parse common region_group patterns like "CA1_2_CA1_1"
# - Keep detection robust to metric names that contain dots (e.g. "adj.P.Val.*")
# - Return NA_character_ when no comparison-like pattern can be detected
# - Numeric group suffixes (1/2/3) will be mapped to labels (con/res/sus)
#   later when building filenames
extract_comparison <- function(col) {
  # Prefer the ".over." style first (e.g. "CA1_2.over.CA1_1")
  comp <- stringr::str_extract(col, "[A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+")
  if (!is.na(comp)) return(comp)
  # Map numeric group codes to labels: 1 -> con, 2 -> res, 3 -> sus
  label_map <- c("1" = "con", "2" = "res", "3" = "sus")

  # Try explicit ".over." style with explicit region_group tokens, e.g. "CA1_2.over.CA1_1"
  m <- stringr::str_match(col, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
    n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
    n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
    return(paste0(r1, n1_label, "_", r2, n2_label))
  }

  # Try underscore-separated tokens like "CA1_2_CA1_1"
  m2 <- stringr::str_match(col, "^([A-Za-z0-9]+)_([0-9]+)_([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m2[1, 1])) {
    r1 <- m2[1, 2]; n1 <- m2[1, 3]; r2 <- m2[1, 4]; n2 <- m2[1, 5]
    n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
    n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
    return(paste0(r1, n1_label, "_", r2, n2_label))
  }

  # Fallback to a simple underscore-capture (preserve previous behavior)
  comp2 <- stringr::str_extract(col, "[A-Za-z0-9_]+_[A-Za-z0-9_]+")
  if (!is.na(comp2)) return(comp2)

  NA_character_
}

comparison_keys <- setNames(vapply(feature_cols, extract_comparison, character(1)), feature_cols)

# Keep only those columns that actually have a comparison pattern with ".over."
# and avoid edge cases where metric names themselves might contain dots 
# (we already pick after first dot)
valid_cols <- names(comparison_keys)[str_detect(comparison_keys, "\\.over\\.")]
valid_keys <- comparison_keys[valid_cols]

# Build mapping: comparison -> columns
by_comparison <- split(valid_cols, valid_keys)

# --------- WRITE ONE CSV PER COMPARISON ----------
# Create an output directory
# base output dir
outdir_base <- "raw"

# extract substring between "E9-" and "_Two-sample_mod_T" from the gct file name
fname <- basename(gct_data)
subfolder <- stringr::str_extract(fname, "(?<=E9-).*?(?=_Two-sample_mod_T)")

# fallback to a more permissive pattern if the strict one fails
if (is.na(subfolder) || !nzchar(subfolder)) {
  subfolder <- stringr::str_extract(fname, "(?<=E9-).*?(?=_Two-sample)")
}

# final fallback name
if (is.na(subfolder) || !nzchar(subfolder)) {
  subfolder <- "unknown-comparison"
}

outdir <- file.path(outdir_base, subfolder)
dir_create(outdir)

# Safe file name helper
safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9._-]", "_") %>%
    str_replace_all("_+", "_")
}

# Iterate over each comparison group and write CSV with id + all columns in that group
# Ensure 'id' stays as a single character column (may contain ";" separators) and is not split.
iwalk(by_comparison, function(cols, comp_key) {
    df_out <- data %>%
        mutate(id = as.character(id)) %>%   # keep any ";"-separated identifiers in one column
        select(id, all_of(cols))

    # If id was stored as a list (edge case), collapse elements into a single ";"-separated string
    if (is.list(df_out$id)) {
        df_out$id <- vapply(df_out$id, function(x) paste(as.character(x), collapse = ";"), character(1))
    }

    # Rename 'id' to 'gene_symbol'
    names(df_out)[1] <- "gene_symbol"

    # Map metric prefixes to desired column names, after removing the comparison suffix
    recode_map <- c(
        "adj.P.Val"   = "padj",
        "P.Value"     = "pval",
        "logFC"       = "log2fc",
        "Log.P.Value" = "logpval",
        "AveExpr"     = "aveExpr",
        "t"           = "t",
        "RawlogFC"    = "rawlog2fc"
    )

  # remove the ".<comp>.over.<comp>" suffix to get the metric portion
  metrics <- stringr::str_remove(cols, "\\.[A-Za-z0-9_]+\\.over\\.[A-Za-z0-9_]+$")

  # produce new names, falling back to the metric itself if not in the recode_map
  new_names <- vapply(metrics, function(m) {
    if (!is.na(m) && m %in% names(recode_map)) recode_map[[m]] else m
  }, character(1))

  # ensure unique column names if needed
  new_names <- make.unique(new_names, sep = "_")

  # assign new names (keep 'id' as first column)
  names(df_out)[-1] <- new_names

  # Invert sign for log2fc and rawlog2fc columns (handle possible make.unique suffixes)
  invert_pattern <- grepl("^(log2fc|rawlog2fc)(?:$|_)", names(df_out), perl = TRUE)
  # skip id column if matched accidentally
  invert_idx <- which(invert_pattern & names(df_out) != "id")
  if (length(invert_idx) > 0) {
    for (j in invert_idx) {
      # coerce to numeric safely and multiply by -1
      df_out[[j]] <- suppressWarnings(as.numeric(as.character(df_out[[j]]))) * -1
    }
  }

  # derive comp2 (comparison string without ".over.") and build filename from it
  label_map <- c("1" = "con", "2" = "res", "3" = "sus")
  comp2 <- comp_key

  # Prefer exact parse of pattern like "REGION_# .over. REGION_#"
  m <- stringr::str_match(comp2, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
  if (!is.na(m[1, 1])) {
    r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
    n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
    n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
    comp2 <- paste0(r1, n1_label, "_", r2, n2_label)  # e.g. CA1res_CA1sus
  } else {
    # Fallback: replace ".over." with "_" then map numeric suffixes 1/2/3 to con/res/sus
    comp2 <- stringr::str_replace_all(comp2, "\\.over\\.", "_")
    comp2 <- stringr::str_replace_all(comp2, "_1", "con")
    comp2 <- stringr::str_replace_all(comp2, "_2", "res")
    comp2 <- stringr::str_replace_all(comp2, "_3", "sus")
  }

  fname <- file.path(outdir, paste0(safe_name(comp2), ".csv"))

  # Use base::write.csv with quote = TRUE so fields containing ";" remain a single quoted cell in the CSV
  utils::write.csv(df_out, fname, row.names = FALSE, quote = TRUE)
  message("Wrote: ", fname)
})

# --------- OPTIONAL: summary index of what was written ----------
index <- tibble(
  comparison = names(by_comparison),
  n_columns = lengths(by_comparison),
  file = vapply(names(by_comparison), function(comp_key) {
    comp2 <- comp_key
    m <- stringr::str_match(comp2, "^([A-Za-z0-9]+)_([0-9]+)\\.over\\.([A-Za-z0-9]+)_([0-9]+)$")
    if (!is.na(m[1, 1])) {
      r1 <- m[1, 2]; n1 <- m[1, 3]; r2 <- m[1, 4]; n2 <- m[1, 5]
      label_map <- c("1" = "con", "2" = "res", "3" = "sus")
      n1_label <- ifelse(n1 %in% names(label_map), label_map[[n1]], n1)
      n2_label <- ifelse(n2 %in% names(label_map), label_map[[n2]], n2)
      comp2 <- paste0(r1, n1_label, "_", r2, n2_label)
    } else {
      comp2 <- stringr::str_replace_all(comp2, "\\.over\\.", "_")
      comp2 <- stringr::str_replace_all(comp2, "_1", "con")
      comp2 <- stringr::str_replace_all(comp2, "_2", "res")
      comp2 <- stringr::str_replace_all(comp2, "_3", "sus")
    }
    file.path(outdir, paste0(safe_name(comp2), ".csv"))
  }, character(1))
)

write_csv(index, file.path(outdir, "indexComparisons.csv"))
message("Wrote index: ", file.path(outdir, "indexComparisons.csv"))