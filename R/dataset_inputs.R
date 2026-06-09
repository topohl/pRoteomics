# Shared dataset-specific input contracts for WGCNA and module scoring.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("validate_dataset", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}

first_existing_path <- function(paths) {
  paths <- unique(normalizePath(paths[nzchar(paths)], winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

latest_matching_file <- function(root, pattern, recursive = TRUE) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)
  candidates <- list.files(root, pattern = pattern, full.names = TRUE, recursive = recursive)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

latest_matching_file_anywhere <- function(roots, pattern, recursive = TRUE) {
  candidates <- vapply(roots, latest_matching_file, character(1), pattern = pattern, recursive = recursive)
  candidates <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

dataset_filter_contract <- function(dataset = current_dataset()) {
  dataset <- validate_dataset(dataset)
  switch(
    dataset,
    neuron_neuropil = list(
      celltype_layer = c("neuron_neuropil", "neuropil"),
      celltype = c("neuron_neuropil", "neuropil", "neuron"),
      layers = c("SO", "SR", "SLM", "MO", "PO", "so", "sr", "slm", "mo", "po"),
      region_only = FALSE,
      description = "Neuron neuropil samples; prefer CellTypeLayer == neuron_neuropil with neuropil layers."
    ),
    neuron_soma = list(
      celltype_layer = c("neuron_soma", "soma"),
      celltype = c("neuron_soma", "soma", "neuron"),
      layers = c("SP", "SG", "sp", "sg"),
      region_only = FALSE,
      description = "Neuron soma samples; accept neuron_soma tokens or SP/SG soma-compatible layers."
    ),
    microglia = list(
      celltype_layer = c("microglia", "microglial"),
      celltype = c("microglia", "microglial"),
      layers = character(),
      region_only = TRUE,
      description = "Microglia samples; region-only spatial structure when no layer exists."
    )
  )
}

resolve_dataset_inputs <- function(dataset = current_dataset(), purpose = c("wgcna", "module_score")) {
  purpose <- match.arg(purpose)
  dataset <- validate_dataset(dataset)
  filter_contract <- dataset_filter_contract(dataset)
  idmap_file <- normalizePath(path_external("MOUSE_10090_idmapping.dat"), winslash = "/", mustWork = FALSE)

  if (identical(purpose, "wgcna")) {
    explicit_expression <- Sys.getenv("PROTEOMICS_WGCNA_UPSTREAM_XLSX", unset = "")
    explicit_metadata <- Sys.getenv("PROTEOMICS_WGCNA_UPSTREAM_META_XLSX", unset = "")
    expression_file <- if (nzchar(explicit_expression)) {
      normalizePath(explicit_expression, winslash = "/", mustWork = FALSE)
    } else {
      latest_matching_file(
        path_processed("01_preprocessing", "impute"),
        paste0("^\\d{8}_pgmatrix_imputed_", dataset, "_[0-9]+samples_missing70pct\\.xlsx$"),
        recursive = FALSE
      )
    }
    metadata_file <- if (nzchar(explicit_metadata)) {
      normalizePath(explicit_metadata, winslash = "/", mustWork = FALSE)
    } else {
      path_metadata("TPE9_sample_metadata_males.xlsx")
    }
    matrix_format <- "imputed_expression_wide"
  } else {
    explicit_expression <- Sys.getenv("PROTEOMICS_MODULE_SCORE_PROTEIN_FILE", unset = "")
    explicit_metadata <- Sys.getenv("PROTEOMICS_MODULE_SCORE_METADATA_FILE", unset = "")
    expression_file <- if (nzchar(explicit_expression)) {
      normalizePath(explicit_expression, winslash = "/", mustWork = FALSE)
    } else {
      latest_matching_file_anywhere(
        c(
          path_processed("01_preprocessing"),
          path_processed("01_preprocessing", "excel_convert"),
          path_processed("morpheus")
        ),
        paste0("^\\d{8}_pgmatrix_imputed_", dataset, "_[0-9]+samples_missing70pct_with_metadata\\.xlsx$")
      )
    }
    canonical_metadata <- path_processed(
      "01_preprocessing",
      "06_merged_metadata_module_score",
      dataset,
      "sample_metadata_merged_clean_for_module_scores.xlsx"
    )
    metadata_file <- first_existing_path(c(canonical_metadata))
    if (is.na(metadata_file) && nzchar(explicit_metadata)) {
      metadata_file <- normalizePath(explicit_metadata, winslash = "/", mustWork = FALSE)
    }
    legacy_dataset_candidates <- c(
      path_results("module_scores", dataset, "sample_metadata_merged_clean_for_module_scores.xlsx"),
      path_processed("01_preprocessing", dataset, "sample_metadata_merged_clean_for_module_scores.xlsx")
    )
    if (is.na(metadata_file)) {
      legacy_dataset_metadata <- first_existing_path(legacy_dataset_candidates)
      if (!is.na(legacy_dataset_metadata)) {
        warning(
          "Using legacy dataset-scoped module-score metadata fallback: ", legacy_dataset_metadata,
          ". Regenerate canonical metadata at data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/.",
          call. = FALSE
        )
        metadata_file <- legacy_dataset_metadata
      }
    }
    allow_global_metadata_fallback <- tolower(Sys.getenv("PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA", unset = "")) %in% c("1", "true", "yes")
    if (is.na(metadata_file) && isTRUE(allow_global_metadata_fallback)) {
      global_metadata <- normalizePath(path_results("module_scores", "sample_metadata_merged_clean_for_module_scores.xlsx"), winslash = "/", mustWork = FALSE)
      if (file.exists(global_metadata)) {
        warning(
          "Using legacy global module-score metadata fallback: ", global_metadata,
          " for dataset '", dataset, "'. Regenerate canonical dataset-scoped metadata at data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/",
          call. = FALSE
        )
        metadata_file <- global_metadata
      }
    }
    matrix_format <- "morpheus_with_metadata_rows"
  }

  diagnostics <- c(
    paste0("dataset=", dataset),
    paste0("purpose=", purpose),
    paste0("matrix_format=", matrix_format),
    paste0("expression_file_exists=", file.exists(expression_file)),
    paste0("metadata_file_exists=", file.exists(metadata_file)),
    paste0("idmap_file_exists=", file.exists(idmap_file)),
    filter_contract$description
  )

  list(
    dataset = dataset,
    purpose = purpose,
    expression_file = expression_file,
    metadata_file = metadata_file,
    idmap_file = idmap_file,
    matrix_format = matrix_format,
    protein_id_col_candidates = c("gene_symbol", "T: Protein.Names", "Genes", "Protein.Group", "ProteinID", "UniProt"),
    sample_id_col_candidates = c("Sample", "sample_id", "SampleID", "SampleColumn", "row.names", "id"),
    expected_filter = filter_contract,
    diagnostics = diagnostics
  )
}

metadata_matches_dataset <- function(metadata, dataset = current_dataset()) {
  dataset <- validate_dataset(dataset)
  contract <- dataset_filter_contract(dataset)
  norm <- function(x) {
    x <- tolower(trimws(as.character(x)))
    gsub("[[:space:]-]+", "_", x)
  }
  keep <- rep(TRUE, nrow(metadata))
  if ("Exclude" %in% names(metadata)) {
    keep <- keep & (is.na(metadata$Exclude) | metadata$Exclude == FALSE)
  } else if ("exclude" %in% names(metadata)) {
    keep <- keep & (is.na(metadata$exclude) | metadata$exclude == FALSE)
  }

  token_hit <- rep(FALSE, nrow(metadata))
  if ("CellTypeLayer" %in% names(metadata)) {
    token_hit <- token_hit | norm(metadata$CellTypeLayer) %in% norm(contract$celltype_layer)
  }
  if ("celltype_layer" %in% names(metadata)) {
    token_hit <- token_hit | norm(metadata$celltype_layer) %in% norm(contract$celltype_layer)
  }
  if ("CellType" %in% names(metadata)) {
    token_hit <- token_hit | norm(metadata$CellType) %in% norm(contract$celltype)
  }
  if ("celltype" %in% names(metadata)) {
    token_hit <- token_hit | norm(metadata$celltype) %in% norm(contract$celltype)
  }
  if (length(contract$layers)) {
    layer_hit <- rep(FALSE, nrow(metadata))
    if ("Layer" %in% names(metadata)) layer_hit <- layer_hit | norm(metadata$Layer) %in% norm(contract$layers)
    if ("layer" %in% names(metadata)) layer_hit <- layer_hit | norm(metadata$layer) %in% norm(contract$layers)
    token_hit <- token_hit | layer_hit
  }
  if (isTRUE(contract$region_only) && !any(token_hit, na.rm = TRUE)) {
    token_hit <- rep(TRUE, nrow(metadata))
  }
  keep & token_hit
}
