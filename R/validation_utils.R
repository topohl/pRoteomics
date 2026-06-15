# Lightweight validation and naming helpers shared across active scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
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

known_pipeline_output_specs <- function() {
  group_effect_cols <- c(
    "dataset", "level", "endpoint_id", "endpoint_label", "contrast",
    "estimate", "SE", "p_value", "FDR_within_dataset_level", "FDR_global",
    "evidence_status", "n_samples", "n_animals", "model_type",
    "formula_used", "rank_deficient_model", "model_warning"
  )
  list(
    module_group_effects.csv = list(required_columns = group_effect_cols),
    supermodule_group_effects.csv = list(required_columns = group_effect_cols),
    WGCNA_module_biological_annotation.csv = list(
      required_columns = c(
        "dataset", "ModuleID", "ModuleColor", "microenvironment_class",
        "microenvironment_label", "classification_rationale", "interpretation_note"
      )
    ),
    WGCNA_supermodule_biological_annotation.csv = list(
      required_columns = c(
        "dataset", "SupermoduleID", "Supermodule_DisplayLabel",
        "dominant_microenvironment_class", "dominant_module_labels",
        "Supermodule_LabelRationale", "interpretation_note"
      ),
      recommended_columns = c(
        "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
        "Supermodule_CompositionShortLabel", "Supermodule_CompositionLabelSource",
        "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
        "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
        "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
        "n_member_modules_with_informative_labels",
        "fraction_member_modules_with_informative_labels"
      )
    ),
    WGCNA_module_group_effects_interpretable.csv = list(
      required_columns = c(
        "dataset", "level", "contrast", "estimate", "p_value", "FDR_global",
        "interpretation_sentence", "ModulePlotLabel", "Supermodule_PlotLabel",
        "Supermodule_FullAnnotationLabel"
      )
    ),
    WGCNA_supermodule_group_effects_interpretable.csv = list(
      required_columns = c(
        "dataset", "level", "contrast", "estimate", "p_value", "FDR_global",
        "interpretation_sentence", "Supermodule_PlotLabel",
        "Supermodule_FullAnnotationLabel", "dominant_microenvironment_class",
        "Supermodule_LabelRationale"
      )
    ),
    biological_claims.csv = list(
      required_columns = c(
        "claim_id", "dataset", "biological_program", "evidence_type",
        "claim_grade", "primary_evidence", "orthogonal_support",
        "major_limitation", "safe_interpretation", "unsafe_overinterpretation",
        "missingness_confounded", "plate_or_batch_confounded",
        "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
        "marker_contamination_risk", "qc_interpretation_flag"
      )
    )
  )
}

validate_known_pipeline_output <- function(path, dataset = NULL) {
  specs <- known_pipeline_output_specs()
  filename <- basename(path)
  if (!filename %in% names(specs)) {
    return(data.frame(path = path, validation_status = "not_applicable", validation_message = "", stringsAsFactors = FALSE))
  }
  if (!file.exists(path)) {
    return(data.frame(path = path, validation_status = "missing", validation_message = "Expected output file does not exist.", stringsAsFactors = FALSE))
  }

  df <- tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) e)
  if (inherits(df, "error")) {
    return(data.frame(path = path, validation_status = "error", validation_message = conditionMessage(df), stringsAsFactors = FALSE))
  }

  messages <- character()
  required <- specs[[filename]]$required_columns
  missing <- setdiff(required, names(df))
  if (length(missing)) messages <- c(messages, paste0("Missing required column(s): ", paste(missing, collapse = ", ")))
  recommended <- specs[[filename]]$recommended_columns %||% character()
  missing_recommended <- setdiff(recommended, names(df))
  if (length(missing_recommended)) messages <- c(messages, paste0("Missing recommended optional column(s): ", paste(missing_recommended, collapse = ", ")))

  if ("dataset" %in% names(df) && exists("valid_datasets", mode = "function")) {
    observed <- unique(as.character(df$dataset[!is.na(df$dataset) & nzchar(as.character(df$dataset))]))
    invalid <- setdiff(observed, valid_datasets())
    if (length(invalid)) messages <- c(messages, paste0("Invalid dataset value(s): ", paste(invalid, collapse = ", ")))
  }
  if (!is.null(dataset) && "dataset" %in% names(df)) {
    observed <- unique(as.character(df$dataset[!is.na(df$dataset) & nzchar(as.character(df$dataset))]))
    mismatch <- setdiff(observed, dataset)
    if (length(mismatch)) messages <- c(messages, paste0("Dataset values do not match selected dataset ", dataset, ": ", paste(mismatch, collapse = ", ")))
  }

  for (col in intersect(c("p_value", "FDR_within_dataset_level", "FDR_global"), names(df))) {
    values <- suppressWarnings(as.numeric(df[[col]]))
    bad <- !is.na(values) & (values < 0 | values > 1)
    if (any(bad)) messages <- c(messages, paste0(col, " has value(s) outside [0, 1]."))
  }
  for (col in intersect(c("n_samples", "n_animals"), names(df))) {
    values <- suppressWarnings(as.numeric(df[[col]]))
    bad <- !is.na(values) & values < 0
    if (any(bad)) messages <- c(messages, paste0(col, " has negative value(s)."))
  }

  if (identical(filename, "WGCNA_supermodule_biological_annotation.csv")) {
    if ("Supermodule_CompositionLabel" %in% names(df)) {
      comp <- trimws(as.character(df$Supermodule_CompositionLabel))
      comp_present <- !is.na(comp) & nzchar(comp)
      unresolved <- comp_present & grepl("mixed|unresolved|low-specificity", comp, ignore.case = TRUE)
      if (any(comp_present) && all(unresolved[comp_present])) {
        messages <- c(messages, "All supermodule composition labels are mixed/unresolved; member-module biological labels may not be propagating.")
      }
    } else {
      label_cols <- intersect(c("Supermodule_FinalLabel", "Supermodule_DisplayLabel", "dominant_module_labels"), names(df))
      if (length(label_cols)) {
        lab <- apply(df[label_cols], 1, paste, collapse = " ")
        if (length(lab) && all(grepl("mixed|unresolved|low-specificity", lab, ignore.case = TRUE))) {
          messages <- c(messages, "All supermodules are mixed/unresolved and no composition labels are available.")
        }
      }
    }
    if ("Supermodule_CompositionLabel" %in% names(df) &&
        "Supermodule_CompositionConfidence" %in% names(df)) {
      comp <- trimws(as.character(df$Supermodule_CompositionLabel))
      conf <- trimws(as.character(df$Supermodule_CompositionConfidence))
      if (any(!is.na(comp) & nzchar(comp)) && !any(!is.na(conf) & nzchar(conf))) {
        messages <- c(messages, "Composition labels are present but all confidence values are missing.")
      }
    }
    if (all(c("dataset", "dominant_microenvironment_class", "Supermodule_CompositionLabel") %in% names(df))) {
      micro_rows <- as.character(df$dataset) == "microglia"
      caution_dom <- as.character(df$dominant_microenvironment_class) %in% c(
        "shared_microenvironment", "neuropil_sensitive",
        "vascular_basement_membrane_ecm", "vascular_bbb_mural"
      )
      overclaim <- grepl("microglia activation|microglia state", as.character(df$Supermodule_CompositionLabel), ignore.case = TRUE)
      if (any(micro_rows & caution_dom & overclaim, na.rm = TRUE)) {
        messages <- c(messages, "Microglia composition labels overclaim microglia specificity despite shared/neuropil/vascular dominance.")
      }
    }
  }

  data.frame(
    path = path,
    validation_status = if (length(messages)) "warning" else "ok",
    validation_message = paste(messages, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

validate_pipeline_outputs <- function(paths, dataset = NULL) {
  paths <- unique(paths[file.exists(paths) & !dir.exists(paths)])
  if (!length(paths)) {
    return(data.frame(path = character(), validation_status = character(), validation_message = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(paths, validate_known_pipeline_output, dataset = dataset))
}
