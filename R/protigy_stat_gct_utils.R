# Shared parsing, validation, and comparison helpers for ProTigy statistical-result GCTs.

protigy_supported_metrics <- function() {
  c(
    "signed.Log.P.Value", "Log.P.Value", "adj.P.Val", "P.Value",
    "RawAveExpr", "RawlogFC", "significant", "sign.logP", "AveExpr",
    "logFC", "t", "B"
  )
}

protigy_required_da_metrics <- function() {
  c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "significant", "Log.P.Value")
}

canonicalize_protigy_comparison <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x <- gsub("_over_", ".over.", x, fixed = TRUE)
  x
}

protigy_comparison_naming_style <- function(x) {
  x <- as.character(x)
  out <- rep("unrecognized", length(x))
  out[grepl(".over.", x, fixed = TRUE)] <- "historical_dot_over"
  out[grepl("_over_", x, fixed = TRUE)] <- "corrected_underscore_over"
  out
}

parse_protigy_stat_field <- function(field, metrics = protigy_supported_metrics()) {
  field <- as.character(field)
  metrics <- metrics[order(nchar(metrics), decreasing = TRUE, method = "radix")]
  metric <- comparison <- rep(NA_character_, length(field))
  naming_style <- rep("unrecognized", length(field))

  for (candidate in metrics) {
    prefix <- paste0(candidate, ".")
    hit <- is.na(metric) & startsWith(field, prefix)
    if (!any(hit)) next
    suffix <- substring(field[hit], nchar(prefix) + 1L)
    style <- protigy_comparison_naming_style(suffix)
    parsed <- style != "unrecognized"
    hit_index <- which(hit)
    metric[hit_index[parsed]] <- candidate
    comparison[hit_index[parsed]] <- canonicalize_protigy_comparison(suffix[parsed])
    naming_style[hit_index[parsed]] <- style[parsed]
  }

  data.frame(
    field = field,
    metric = metric,
    comparison = comparison,
    naming_style = naming_style,
    stringsAsFactors = FALSE
  )
}

canonicalize_protigy_unit_token <- function(unit, dataset) {
  dataset <- validate_dataset(dataset)
  unit <- trimws(as.character(unit))
  m <- regexec("^(CA[123]|DG)_(.+)$", unit, ignore.case = TRUE)
  p <- regmatches(unit, m)[[1]]
  if (!length(p)) return(NA_character_)
  region <- toupper(p[[2]])
  suffix <- tolower(p[[3]])

  valid <- switch(
    dataset,
    neuron_neuropil = paste(region, suffix, sep = "_") %in% c(
      "CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
      "CA3_so", "CA3_sr", "DG_mo", "DG_po"
    ),
    neuron_soma = paste(region, suffix, sep = "_") %in% c(
      "CA1_sp", "CA2_sp", "CA3_sp", "DG_sg"
    ),
    microglia = identical(suffix, "microglia"),
    FALSE
  )
  if (!isTRUE(valid)) return(NA_character_)
  paste(region, suffix, sep = "_")
}

protigy_spatial_unit_label <- function(unit, dataset) {
  unit <- canonicalize_protigy_unit_token(unit, dataset)
  if (is.na(unit)) return(NA_character_)
  if (identical(validate_dataset(dataset), "microglia")) {
    return(sub("_microglia$", "", unit))
  }
  unit
}

parse_protigy_comparison <- function(comparison, dataset) {
  canonical <- canonicalize_protigy_comparison(comparison)
  m <- regexec("^(.*)_([123])\\.over\\.(.*)_([123])$", canonical)
  p <- regmatches(canonical, m)[[1]]
  if (!length(p)) {
    return(data.frame(
      comparison = canonical, parsed = FALSE, left_unit = NA_character_,
      left_group = NA_character_, right_unit = NA_character_, right_group = NA_character_,
      left_canonical_unit = NA_character_, right_canonical_unit = NA_character_,
      spatial_unit = NA_character_, same_unit = FALSE, allowed_group_pair = FALSE,
      valid_stress_comparison = FALSE, biological_contrast = NA_character_,
      rejection_reason = "unparseable_comparison", stringsAsFactors = FALSE
    ))
  }

  left_unit <- p[[2]]
  left_group <- p[[3]]
  right_unit <- p[[4]]
  right_group <- p[[5]]
  left_canonical <- canonicalize_protigy_unit_token(left_unit, dataset)
  right_canonical <- canonicalize_protigy_unit_token(right_unit, dataset)
  units_valid <- !is.na(left_canonical) && !is.na(right_canonical)
  same_unit <- units_valid && identical(left_canonical, right_canonical)
  pair <- paste0(left_group, "/", right_group)
  allowed_pair <- pair %in% c("2/1", "3/2", "3/1")
  valid <- same_unit && allowed_pair
  biological <- c("2/1" = "RES_vs_CON", "3/2" = "SUS_vs_RES", "3/1" = "SUS_vs_CON")
  reason <- if (!units_valid) {
    "invalid_canonical_spatial_unit"
  } else if (!same_unit) {
    "cross_spatial_unit"
  } else if (!allowed_pair) {
    "unsupported_group_pair"
  } else {
    "accepted_primary_stress_contrast"
  }

  data.frame(
    comparison = canonical,
    parsed = TRUE,
    left_unit = left_unit,
    left_group = left_group,
    right_unit = right_unit,
    right_group = right_group,
    left_canonical_unit = left_canonical,
    right_canonical_unit = right_canonical,
    spatial_unit = if (valid) protigy_spatial_unit_label(left_canonical, dataset) else NA_character_,
    same_unit = same_unit,
    allowed_group_pair = allowed_pair,
    valid_stress_comparison = valid,
    biological_contrast = if (allowed_pair) unname(biological[[pair]]) else NA_character_,
    rejection_reason = reason,
    stringsAsFactors = FALSE
  )
}

validate_protigy_comparisons <- function(comparisons, dataset, strict_primary = FALSE) {
  comparisons <- as.character(comparisons)
  rows <- lapply(comparisons, parse_protigy_comparison, dataset = dataset)
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out$duplicate_comparison <- duplicated(out$comparison) | duplicated(out$comparison, fromLast = TRUE)
    rownames(out) <- NULL
  }
  if (isTRUE(strict_primary) && nrow(out) && any(!out$valid_stress_comparison)) {
    rejected <- unique(out$comparison[!out$valid_stress_comparison])
    stop(
      "Corrected animal-level GCT contains rejected comparison(s): ",
      paste(rejected, collapse = ", "),
      call. = FALSE
    )
  }
  out
}

protigy_metric_is_signed <- function(metric) {
  as.character(metric) %in% c("logFC", "RawlogFC", "t")
}

reverse_protigy_metric <- function(x, metric) {
  if (protigy_metric_is_signed(metric)) -suppressWarnings(as.numeric(x)) else x
}

reverse_protigy_metric_frame <- function(df, metric_by_column) {
  out <- df
  shared <- intersect(names(metric_by_column), names(out))
  for (column in shared) {
    if (protigy_metric_is_signed(metric_by_column[[column]])) {
      out[[column]] <- -suppressWarnings(as.numeric(out[[column]]))
    }
  }
  out
}

protigy_is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[\\\\/]|\\\\\\\\|/)", path)
}

resolve_protigy_root <- function(value, default_root) {
  value <- trimws(as.character(value))
  path <- if (!nzchar(value)) {
    default_root
  } else if (protigy_is_absolute_path(value)) {
    path.expand(value)
  } else {
    repo_path(value)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_gct_io_roots <- function(
    input_root = Sys.getenv("PROTEOMICS_GCT_INPUT_ROOT", unset = ""),
    output_root = Sys.getenv("PROTEOMICS_GCT_OUTPUT_ROOT", unset = "")) {
  list(
    input_root = resolve_protigy_root(
      input_root,
      path_processed("01_preprocessing", "protigy_output")
    ),
    output_root = resolve_protigy_root(
      output_root,
      path_processed("01_preprocessing", "gct_extractR")
    )
  )
}

resolve_single_protigy_gct <- function(root, dataset, input_stem = "", use_manifest = TRUE) {
  dataset <- validate_dataset(dataset)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  dataset_dir <- file.path(root, dataset)
  input_stem <- trimws(as.character(input_stem))

  if (!nzchar(input_stem) && isTRUE(use_manifest)) {
    manifest_candidates <- c(
      file.path(dataset_dir, "protigy_manifest.csv"),
      file.path(root, "protigy_manifest.csv")
    )
    manifest_file <- manifest_candidates[file.exists(manifest_candidates)][1]
    if (!is.na(manifest_file)) {
      manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
      path_col <- intersect(c("gct_path", "file", "path"), names(manifest))[1]
      if (!is.na(path_col) && nrow(manifest)) {
        candidate <- as.character(manifest[[path_col]][1])
        if (!is.na(candidate) && nzchar(candidate)) {
          candidate_path <- if (protigy_is_absolute_path(candidate)) {
            path.expand(candidate)
          } else {
            file.path(dataset_dir, basename(candidate))
          }
          if (file.exists(candidate_path)) {
            input_stem <- tools::file_path_sans_ext(basename(candidate_path))
          }
        }
      }
    }
  }

  if (nzchar(input_stem)) {
    path <- file.path(dataset_dir, paste0(input_stem, ".gct"))
    if (!file.exists(path)) stop("GCT input not found: ", path, call. = FALSE)
    return(list(
      root = root,
      dataset_dir = normalizePath(dataset_dir, winslash = "/", mustWork = FALSE),
      input_stem = input_stem,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      resolution = "input_stem_or_manifest"
    ))
  }

  candidates <- if (dir.exists(dataset_dir)) {
    sort(list.files(dataset_dir, pattern = "\\.gct$", full.names = TRUE, ignore.case = TRUE), method = "radix")
  } else {
    character()
  }
  if (length(candidates) != 1L) {
    detail <- if (!length(candidates)) "No .gct files were found." else paste0(
      "Multiple .gct files were found: ", paste(basename(candidates), collapse = ", ")
    )
    stop(
      "Could not resolve exactly one ProTigy GCT for dataset '", dataset,
      "' under ", dataset_dir, ". ", detail,
      " Set PROTEOMICS_GCT_INPUT_STEM when needed.",
      call. = FALSE
    )
  }
  list(
    root = root,
    dataset_dir = normalizePath(dataset_dir, winslash = "/", mustWork = FALSE),
    input_stem = tools::file_path_sans_ext(basename(candidates[[1]])),
    path = normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE),
    resolution = "single_gct_auto_detected"
  )
}

read_protigy_stat_gct <- function(path, dataset, strict_primary = FALSE) {
  dataset <- validate_dataset(dataset)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  first <- readLines(path, n = 2L, warn = FALSE)
  if (length(first) < 2L || trimws(strsplit(first[[1]], "\t", fixed = TRUE)[[1]][[1]]) != "#1.3") {
    stop("Expected a GCT v1.3 file: ", path, call. = FALSE)
  }
  dims <- suppressWarnings(as.integer(strsplit(first[[2]], "\t", fixed = TRUE)[[1]][1:4]))
  if (length(dims) != 4L || anyNA(dims)) stop("Invalid GCT v1.3 dimension line: ", path, call. = FALSE)
  names(dims) <- c("nrmat", "ncmat", "nrhd", "nchd")

  raw <- utils::read.delim(
    path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    check.names = FALSE, comment.char = "", quote = "", fill = TRUE,
    na.strings = NULL
  )
  header_candidates <- which(tolower(trimws(as.character(raw[[1]]))) == "id")
  if (!length(header_candidates)) stop("Could not find the GCT header row in: ", path, call. = FALSE)
  header_row <- header_candidates[[1]]
  original_names <- as.character(unlist(raw[header_row, , drop = TRUE], use.names = FALSE))
  data_start <- header_row + 1L + dims[["nchd"]]
  data_end <- data_start + dims[["nrmat"]] - 1L
  if (data_end > nrow(raw)) {
    stop("GCT contains fewer protein rows than declared: ", path, call. = FALSE)
  }
  data <- raw[data_start:data_end, , drop = FALSE]
  internal_names <- make.unique(ifelse(nzchar(original_names), original_names, "__blank__"), sep = "__duplicate_")
  names(data) <- internal_names
  id_column <- internal_names[[which(tolower(trimws(original_names)) == "id")[[1]]]]
  ids <- trimws(as.character(data[[id_column]]))

  fields <- parse_protigy_stat_field(original_names)
  fields$column_index <- seq_along(original_names)
  fields$column_internal <- internal_names
  parsed_fields <- fields[!is.na(fields$comparison), , drop = FALSE]
  parsed_fields$metric_comparison_key <- paste(parsed_fields$metric, parsed_fields$comparison, sep = "||")
  parsed_fields$duplicate_metric_comparison <- duplicated(parsed_fields$metric_comparison_key) |
    duplicated(parsed_fields$metric_comparison_key, fromLast = TRUE)

  comparisons <- unique(parsed_fields$comparison)
  comparison_validation <- validate_protigy_comparisons(comparisons, dataset, strict_primary = strict_primary)
  if (nrow(comparison_validation)) {
    style_map <- vapply(comparison_validation$comparison, function(comp) {
      styles <- unique(parsed_fields$naming_style[parsed_fields$comparison == comp])
      paste(sort(styles, method = "radix"), collapse = ";")
    }, character(1))
    comparison_validation$naming_style <- unname(style_map[comparison_validation$comparison])
  }

  description_index <- which(tolower(trimws(original_names)) == "description")
  description <- if (length(description_index)) {
    as.character(data[[internal_names[[description_index[[1]]]]]])
  } else {
    rep(NA_character_, length(ids))
  }

  required_metrics <- protigy_required_da_metrics()
  valid_comparisons <- comparison_validation$comparison[comparison_validation$valid_stress_comparison]
  required_grid <- if (length(valid_comparisons)) {
    expand.grid(metric = required_metrics, comparison = valid_comparisons, stringsAsFactors = FALSE)
  } else {
    data.frame(metric = character(), comparison = character(), stringsAsFactors = FALSE)
  }
  present_keys <- unique(parsed_fields$metric_comparison_key)
  required_grid$present <- paste(required_grid$metric, required_grid$comparison, sep = "||") %in% present_keys
  missing_required <- required_grid[!required_grid$present, c("metric", "comparison"), drop = FALSE]

  numeric_metrics <- setdiff(required_metrics, "significant")
  numeric_fields <- parsed_fields[
    parsed_fields$metric %in% numeric_metrics & parsed_fields$comparison %in% valid_comparisons,
    , drop = FALSE
  ]
  numeric_checks <- lapply(seq_len(nrow(numeric_fields)), function(i) {
    values <- trimws(as.character(data[[numeric_fields$column_internal[[i]]]]))
    missing <- is.na(values) | !nzchar(values) | toupper(values) %in% c("NA", "NAN")
    parsed <- suppressWarnings(as.numeric(values))
    data.frame(
      metric = numeric_fields$metric[[i]],
      comparison = numeric_fields$comparison[[i]],
      n_values = length(values),
      n_missing = sum(missing),
      n_non_numeric = sum(!missing & is.na(parsed)),
      stringsAsFactors = FALSE
    )
  })
  numeric_validation <- if (length(numeric_checks)) do.call(rbind, numeric_checks) else data.frame()

  logfc_rows <- parsed_fields[parsed_fields$metric == "logFC", , drop = FALSE]
  duplicate_comparisons <- sum(duplicated(logfc_rows$comparison))
  expected_physical_fields <- 1L + dims[["nrhd"]] + dims[["ncmat"]]
  observed_physical_fields <- length(original_names)
  styles <- sort(unique(parsed_fields$naming_style), method = "radix")

  list(
    dataset = dataset,
    path = path,
    sha256 = file_hash_sha256(path),
    dimensions = dims,
    expected_physical_fields = expected_physical_fields,
    observed_physical_fields = observed_physical_fields,
    physical_dimension_match = identical(as.integer(expected_physical_fields), as.integer(observed_physical_fields)),
    data = data,
    original_names = original_names,
    internal_names = internal_names,
    id_column = id_column,
    ids = ids,
    description = description,
    fields = fields,
    parsed_fields = parsed_fields,
    comparison_validation = comparison_validation,
    metrics_found = sort(unique(parsed_fields$metric), method = "radix"),
    missing_required_fields = missing_required,
    numeric_validation = numeric_validation,
    naming_style = if (!length(styles)) "none" else if (length(styles) == 1L) styles else paste0("mixed:", paste(styles, collapse = "+")),
    duplicate_comparisons = duplicate_comparisons,
    duplicate_metric_comparison_fields = sum(parsed_fields$duplicate_metric_comparison),
    protein_ids_unique = !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids),
    description_present = length(description_index) > 0L,
    description_aligned = length(description) == length(ids) && all(!is.na(description) & nzchar(trimws(description)))
  )
}

protigy_gct_diagnostic_row <- function(gct) {
  validation <- gct$comparison_validation
  valid <- validation$valid_stress_comparison
  cross <- validation$parsed & !validation$same_unit &
    !is.na(validation$left_canonical_unit) & !is.na(validation$right_canonical_unit)
  missing_metrics <- if (nrow(gct$missing_required_fields)) {
    paste(sort(unique(gct$missing_required_fields$metric), method = "radix"), collapse = ";")
  } else {
    ""
  }
  data.frame(
    dataset = gct$dataset,
    path = gct$path,
    sha256 = gct$sha256,
    nrmat = unname(gct$dimensions[["nrmat"]]),
    ncmat = unname(gct$dimensions[["ncmat"]]),
    nrhd = unname(gct$dimensions[["nrhd"]]),
    nchd = unname(gct$dimensions[["nchd"]]),
    observed_physical_fields = gct$observed_physical_fields,
    declared_physical_fields = gct$expected_physical_fields,
    physical_dimension_match = gct$physical_dimension_match,
    n_protein_ids = length(gct$ids),
    protein_ids_unique = gct$protein_ids_unique,
    description_present = gct$description_present,
    description_aligned = gct$description_aligned,
    comparison_naming_style = gct$naming_style,
    n_detected_comparisons = nrow(validation),
    n_valid_stress_comparisons = sum(valid),
    n_cross_unit_comparisons = sum(cross),
    n_duplicate_comparisons = gct$duplicate_comparisons,
    n_duplicate_metric_comparison_fields = gct$duplicate_metric_comparison_fields,
    metrics_found = paste(gct$metrics_found, collapse = ";"),
    metrics_missing = missing_metrics,
    n_missing_numeric_values = if (nrow(gct$numeric_validation)) sum(gct$numeric_validation$n_missing) else NA_integer_,
    n_non_numeric_values = if (nrow(gct$numeric_validation)) sum(gct$numeric_validation$n_non_numeric) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

validate_corrected_protigy_gct_contract <- function(gct) {
  validation <- gct$comparison_validation
  accepted <- validation[validation$valid_stress_comparison, , drop = FALSE]
  expected_units <- switch(
    gct$dataset,
    neuron_neuropil = c(
      "CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
      "CA3_so", "CA3_sr", "DG_mo", "DG_po"
    ),
    neuron_soma = c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg"),
    microglia = sort(unique(accepted$spatial_unit), method = "radix")
  )
  observed_units <- sort(unique(accepted$spatial_unit), method = "radix")
  if (!identical(sort(expected_units, method = "radix"), observed_units)) {
    stop(
      "Corrected GCT spatial-unit coverage mismatch for ", gct$dataset,
      ". Expected: ", paste(expected_units, collapse = ", "),
      "; observed: ", paste(observed_units, collapse = ", "),
      call. = FALSE
    )
  }
  expected_pairs <- c("2/1", "3/2", "3/1")
  coverage <- lapply(expected_units, function(unit) {
    x <- accepted[accepted$spatial_unit == unit, , drop = FALSE]
    pairs <- paste0(x$left_group, "/", x$right_group)
    data.frame(
      spatial_unit = unit,
      n_comparisons = nrow(x),
      observed_pairs = paste(sort(pairs, method = "radix"), collapse = ";"),
      complete = nrow(x) == 3L && setequal(pairs, expected_pairs) && !anyDuplicated(pairs),
      stringsAsFactors = FALSE
    )
  })
  coverage <- do.call(rbind, coverage)
  expected_count <- length(expected_units) * length(expected_pairs)
  failures <- character()
  if (nrow(validation) != expected_count) failures <- c(failures, "unexpected comparison count")
  if (any(!coverage$complete)) failures <- c(failures, "incomplete per-unit contrast coverage")
  if (!gct$protein_ids_unique) failures <- c(failures, "protein IDs are missing or duplicated")
  if (!gct$description_present || !gct$description_aligned) failures <- c(failures, "Description is missing or misaligned")
  if (nrow(gct$missing_required_fields)) failures <- c(failures, "required metric/comparison fields are missing")
  if (gct$duplicate_comparisons > 0L || gct$duplicate_metric_comparison_fields > 0L) {
    failures <- c(failures, "duplicate comparison fields are present")
  }
  if (nrow(gct$numeric_validation) && any(gct$numeric_validation$n_missing > 0L | gct$numeric_validation$n_non_numeric > 0L)) {
    failures <- c(failures, "required statistics contain missing or non-numeric values")
  }
  if (length(failures)) {
    stop(
      "Corrected ProTigy GCT contract failed for ", gct$dataset, ": ",
      paste(unique(failures), collapse = "; "),
      call. = FALSE
    )
  }
  list(
    expected_spatial_units = expected_units,
    expected_comparison_count = expected_count,
    coverage = coverage,
    status = "PASS"
  )
}

extract_protigy_metric_long <- function(gct, comparisons = NULL,
                                        metrics = c("logFC", "t", "P.Value", "adj.P.Val", "B", "AveExpr")) {
  if (is.null(comparisons)) {
    comparisons <- gct$comparison_validation$comparison[gct$comparison_validation$valid_stress_comparison]
  }
  comparisons <- canonicalize_protigy_comparison(comparisons)
  fields <- gct$parsed_fields[
    gct$parsed_fields$comparison %in% comparisons & gct$parsed_fields$metric %in% metrics,
    , drop = FALSE
  ]
  fields <- fields[!duplicated(fields$metric_comparison_key), , drop = FALSE]
  validation <- gct$comparison_validation
  rows <- lapply(comparisons, function(comp) {
    selected <- fields[fields$comparison == comp, , drop = FALSE]
    meta <- validation[validation$comparison == comp, , drop = FALSE]
    out <- data.frame(
      dataset = gct$dataset,
      spatial_unit = meta$spatial_unit[[1]],
      comparison = comp,
      biological_contrast = meta$biological_contrast[[1]],
      protein_id = gct$ids,
      description = gct$description,
      stringsAsFactors = FALSE
    )
    for (metric in metrics) {
      hit <- selected[selected$metric == metric, , drop = FALSE]
      out[[metric]] <- if (nrow(hit)) {
        suppressWarnings(as.numeric(as.character(gct$data[[hit$column_internal[[1]]]])))
      } else {
        rep(NA_real_, length(gct$ids))
      }
    }
    out
  })
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

match_protigy_proteins <- function(legacy_ids, animal_ids) {
  legacy_ids <- as.character(legacy_ids)
  animal_ids <- as.character(animal_ids)
  list(
    matched = sort(intersect(legacy_ids, animal_ids), method = "radix"),
    legacy_only = sort(setdiff(legacy_ids, animal_ids), method = "radix"),
    animal_only = sort(setdiff(animal_ids, legacy_ids), method = "radix")
  )
}

deterministic_signed_rank <- function(statistic, protein_id) {
  statistic <- suppressWarnings(as.numeric(statistic))
  protein_id <- as.character(protein_id)
  out <- rep(NA_integer_, length(statistic))
  finite <- is.finite(statistic) & !is.na(protein_id) & nzchar(protein_id)
  if (!any(finite)) return(out)
  idx <- which(finite)
  ordered <- idx[order(-statistic[idx], protein_id[idx], method = "radix")]
  out[ordered] <- seq_along(ordered)
  out
}
