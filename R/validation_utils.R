# Lightweight validation and naming helpers shared across active scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

safe_name <- function(x, max_chars = 180) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- ifelse(nzchar(x), x, "unnamed")
  substr(x, 1, max_chars)
}

detect_column <- function(df, candidates, required = FALSE, context = "data frame") {
  nms <- names(df)
  nms_clean <- tolower(gsub("[^a-z0-9]", "", nms))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  idx <- match(cand_clean, nms_clean)
  idx <- idx[!is.na(idx)]
  if (length(idx)) return(nms[idx[[1]]])
  if (isTRUE(required)) {
    stop("Could not find required column in ", context, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- sub("\\.d$", "", x, ignore.case = TRUE)
  x <- trimws(tolower(x))
  x <- gsub("[[:space:]]+", "_", x)
  x
}

sample_overlap_summary <- function(matrix_samples, metadata_samples) {
  matrix_norm <- normalize_sample_id(matrix_samples)
  metadata_norm <- normalize_sample_id(metadata_samples)
  overlap <- intersect(matrix_norm, metadata_norm)
  data.frame(
    n_matrix_samples = length(unique(matrix_norm[!is.na(matrix_norm) & nzchar(matrix_norm)])),
    n_metadata_samples = length(unique(metadata_norm[!is.na(metadata_norm) & nzchar(metadata_norm)])),
    n_overlap = length(unique(overlap)),
    overlap_fraction_matrix = if (length(unique(matrix_norm))) length(unique(overlap)) / length(unique(matrix_norm)) else NA_real_,
    overlap_fraction_metadata = if (length(unique(metadata_norm))) length(unique(overlap)) / length(unique(metadata_norm)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

validate_manifest_paths <- function(manifest, path_cols = c("input_gene_file", "output_table", "output_plot"), allow_missing = TRUE) {
  if (is.null(manifest) || !nrow(manifest)) {
    return(data.frame(path_column = character(), path = character(), exists = logical(), stringsAsFactors = FALSE))
  }
  cols <- intersect(path_cols, names(manifest))
  out <- do.call(rbind, lapply(cols, function(col) {
    vals <- unique(as.character(manifest[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    data.frame(path_column = col, path = vals, exists = file.exists(vals), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) out <- data.frame(path_column = character(), path = character(), exists = logical(), stringsAsFactors = FALSE)
  if (!allow_missing && any(!out$exists)) {
    stop("Manifest contains missing paths:\n", paste(out$path[!out$exists], collapse = "\n"), call. = FALSE)
  }
  out
}

duplicate_key_summary <- function(df, keys) {
  keys <- intersect(keys, names(df))
  if (!length(keys) || is.null(df) || !nrow(df)) {
    return(data.frame(n_duplicate_keys = 0L, stringsAsFactors = FALSE))
  }
  key <- do.call(paste, c(df[keys], sep = "\r"))
  data.frame(n_duplicate_keys = sum(duplicated(key)), stringsAsFactors = FALSE)
}

interpretation_strength <- function(fdr = NA_real_, effect_size = NA_real_, n = NA_integer_, bootstrap_support = NA_real_) {
  fdr <- suppressWarnings(as.numeric(fdr))
  effect_size <- suppressWarnings(abs(as.numeric(effect_size)))
  n <- suppressWarnings(as.numeric(n))
  bootstrap_support <- suppressWarnings(as.numeric(bootstrap_support))
  if (!is.na(fdr) && fdr <= 0.05 && (is.na(effect_size) || effect_size >= 0.3) && (is.na(n) || n >= 6) && (is.na(bootstrap_support) || bootstrap_support >= 0.7)) {
    return("strong")
  }
  if (!is.na(fdr) && fdr <= 0.10 && (is.na(n) || n >= 6)) return("moderate")
  "exploratory"
}
