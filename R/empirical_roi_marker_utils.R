# Utilities for animal-level empirical ROI marker discovery.

empirical_roi_marker_contract_version <- function() {
  "empirical_roi_marker_v2_animal_paired_limma"
}

empirical_roi_dataset_levels <- function() {
  c("neuron_neuropil", "neuron_soma", "microglia")
}

empirical_roi_required_metadata_fields <- function() {
  c(
    "sample_id", "AnimalID", "ReplicateGroup", "region", "layer",
    "celltype_layer", "ExpGroup"
  )
}

empirical_roi_nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

empirical_roi_collapse_unique <- function(x) {
  x <- sort(unique(as.character(x[empirical_roi_nonempty(x)])), method = "radix")
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

validate_known_empirical_metadata_corrections <- function(metadata, require_all = TRUE) {
  required <- empirical_roi_required_metadata_fields()
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop(
      "Metadata correction audit is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  expected <- data.frame(
    sample_number = c("S129", "S130", "S131", "S132"),
    AnimalID_expected = "A0003",
    ReplicateGroup_expected = c("Left", "Right", "Left", "Right"),
    region_expected = "CA2",
    layer_expected = c("sr", "sr", "sp", "sp"),
    celltype_layer_expected = c(
      "neuron_neuropil", "neuron_neuropil", "neuron_soma", "neuron_soma"
    ),
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(expected)), function(i) {
    token <- expected$sample_number[[i]]
    hit <- grepl(paste0("(^|[_/\\\\])", token, "([_.]|$)"), as.character(metadata$sample_id), perl = TRUE)
    n_matches <- sum(hit)
    if (isTRUE(require_all) && n_matches != 1L) {
      stop("Expected exactly one authoritative metadata row for ", token, "; found ", n_matches, ".", call. = FALSE)
    }
    if (n_matches == 0L) {
      return(data.frame(
        sample_number = token,
        n_matches = 0L,
        correction_valid = NA,
        sample_id = NA_character_,
        AnimalID = NA_character_,
        ReplicateGroup = NA_character_,
        region = NA_character_,
        layer = NA_character_,
        celltype_layer = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    if (n_matches > 1L) {
      stop("Metadata sample token is not unique: ", token, call. = FALSE)
    }
    observed <- metadata[which(hit), , drop = FALSE]
    valid <- identical(as.character(observed$AnimalID[[1]]), expected$AnimalID_expected[[i]]) &&
      identical(as.character(observed$ReplicateGroup[[1]]), expected$ReplicateGroup_expected[[i]]) &&
      identical(as.character(observed$region[[1]]), expected$region_expected[[i]]) &&
      identical(tolower(as.character(observed$layer[[1]])), expected$layer_expected[[i]]) &&
      identical(as.character(observed$celltype_layer[[1]]), expected$celltype_layer_expected[[i]])
    data.frame(
      sample_number = token,
      n_matches = 1L,
      correction_valid = valid,
      sample_id = as.character(observed$sample_id[[1]]),
      AnimalID = as.character(observed$AnimalID[[1]]),
      ReplicateGroup = as.character(observed$ReplicateGroup[[1]]),
      region = as.character(observed$region[[1]]),
      layer = tolower(as.character(observed$layer[[1]])),
      celltype_layer = as.character(observed$celltype_layer[[1]]),
      stringsAsFactors = FALSE
    )
  })
  audit <- dplyr::bind_rows(rows)
  if (isTRUE(require_all) && !all(audit$correction_valid %in% TRUE)) {
    stop("Authoritative S129/S130/S131/S132 metadata corrections are not valid.", call. = FALSE)
  }
  audit
}

validate_empirical_roi_metadata <- function(mat, metadata, dataset) {
  dataset <- match.arg(dataset, empirical_roi_dataset_levels())
  if (!is.matrix(mat) || is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Empirical ROI expression input must be a named ProteinGroupID x sample matrix.", call. = FALSE)
  }
  if (ncol(mat) != nrow(metadata)) {
    stop("Empirical ROI matrix columns and metadata rows have different lengths.", call. = FALSE)
  }
  if (!identical(rownames(metadata), colnames(mat))) {
    stop("Empirical ROI metadata must be aligned exactly to matrix columns.", call. = FALSE)
  }
  required <- empirical_roi_required_metadata_fields()
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop("Empirical ROI metadata is missing required fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(rownames(mat)) || any(!empirical_roi_nonempty(rownames(mat)))) {
    stop("ProteinGroupID row names must be unique and nonmissing.", call. = FALSE)
  }
  if (anyDuplicated(as.character(metadata$sample_id)) ||
      !identical(as.character(metadata$sample_id), colnames(mat))) {
    stop("Structured metadata sample_id must uniquely and exactly match matrix column names.", call. = FALSE)
  }
  for (column in required) {
    if (!all(empirical_roi_nonempty(metadata[[column]]))) {
      stop("Structured metadata field has missing/blank values: ", column, call. = FALSE)
    }
  }
  if (!all(as.character(metadata$celltype_layer) == dataset)) {
    stop("Structured celltype_layer disagrees with the requested dataset: ", dataset, call. = FALSE)
  }
  if (!all(as.character(metadata$ReplicateGroup) %in% c("Left", "Right"))) {
    stop("ReplicateGroup must contain only exact Left/Right values.", call. = FALSE)
  }
  if (!all(as.character(metadata$region) %in% c("CA1", "CA2", "CA3", "DG"))) {
    stop("Structured region must be one of CA1, CA2, CA3, or DG.", call. = FALSE)
  }
  allowed_layers <- list(
    neuron_neuropil = c("slm", "so", "sr", "mo", "po"),
    neuron_soma = c("sp", "sg"),
    microglia = "microglia"
  )
  if (!all(tolower(as.character(metadata$layer)) %in% allowed_layers[[dataset]])) {
    stop("Structured layer contains values outside the dataset contract for ", dataset, ".", call. = FALSE)
  }
  animal_groups <- metadata |>
    dplyr::transmute(AnimalID = as.character(.data$AnimalID), ExpGroup = as.character(.data$ExpGroup)) |>
    dplyr::distinct() |>
    dplyr::count(.data$AnimalID, name = "n_groups")
  if (any(animal_groups$n_groups != 1L)) {
    stop("Each AnimalID must map to exactly one ExpGroup; ExpGroup is provenance only.", call. = FALSE)
  }
  invisible(TRUE)
}

aggregate_empirical_matrix_equal_weight <- function(mat, metadata, group_cols) {
  missing <- setdiff(group_cols, names(metadata))
  if (length(missing)) stop("Aggregation metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (ncol(mat) != nrow(metadata)) stop("Aggregation matrix/metadata dimensions disagree.", call. = FALSE)
  group_frame <- data.frame(lapply(metadata[group_cols], as.character), check.names = FALSE)
  if (any(!vapply(group_frame, function(x) all(empirical_roi_nonempty(x)), logical(1)))) {
    stop("Aggregation group fields must be complete.", call. = FALSE)
  }
  key <- do.call(paste, c(group_frame, sep = "\u001f"))
  group_levels <- sort(unique(key), method = "radix")
  index <- split(seq_len(ncol(mat)), factor(key, levels = group_levels), drop = TRUE)
  values <- lapply(index, function(columns) {
    out <- rowMeans(mat[, columns, drop = FALSE], na.rm = TRUE)
    out[is.nan(out)] <- NA_real_
    out
  })
  aggregated <- do.call(cbind, values)
  if (is.null(dim(aggregated))) aggregated <- matrix(aggregated, nrow = nrow(mat))
  rownames(aggregated) <- rownames(mat)
  colnames(aggregated) <- group_levels
  first_rows <- vapply(index, function(columns) columns[[1]], integer(1))
  group_metadata <- group_frame[first_rows, , drop = FALSE]
  rownames(group_metadata) <- group_levels
  group_metadata$n_input_units <- as.integer(lengths(index))
  list(mat = aggregated, meta = group_metadata)
}

empirical_roi_spatial_unit_qc <- function(metadata, dataset) {
  dataset_name <- dataset
  spatial_cols <- if (dataset == "neuron_neuropil") {
    c("AnimalID", "region", "layer")
  } else {
    c("AnimalID", "region")
  }
  metadata |>
    dplyr::mutate(
      AnimalID = as.character(.data$AnimalID),
      region = as.character(.data$region),
      layer = tolower(as.character(.data$layer)),
      ReplicateGroup = as.character(.data$ReplicateGroup),
      sample_id = as.character(.data$sample_id)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(spatial_cols))) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_left = sum(.data$ReplicateGroup == "Left"),
      n_right = sum(.data$ReplicateGroup == "Right"),
      sample_ids = paste(sort(.data$sample_id, method = "radix"), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dataset = dataset_name,
      side_status = dplyr::case_when(
        .data$n_samples == 2L & .data$n_left == 1L & .data$n_right == 1L ~ "bilateral",
        .data$n_samples == 1L & (.data$n_left + .data$n_right) == 1L ~ "one_sided",
        TRUE ~ "invalid"
      ),
      spatial_unit = if (dataset_name == "neuron_neuropil") {
        paste(.data$AnimalID, .data$region, .data$layer, sep = "::")
      } else {
        paste(.data$AnimalID, .data$region, sep = "::")
      }
    )
}

aggregate_empirical_roi_dataset <- function(mat, metadata, dataset) {
  dataset <- match.arg(dataset, empirical_roi_dataset_levels())
  validate_empirical_roi_metadata(mat, metadata, dataset)
  metadata <- metadata |>
    dplyr::mutate(
      AnimalID = as.character(.data$AnimalID),
      ReplicateGroup = as.character(.data$ReplicateGroup),
      region = as.character(.data$region),
      layer = tolower(as.character(.data$layer)),
      ExpGroup = as.character(.data$ExpGroup)
    )
  spatial_qc <- empirical_roi_spatial_unit_qc(metadata, dataset)
  if (any(spatial_qc$side_status == "invalid")) {
    stop("Invalid hemisphere structure in ", dataset, ": expected bilateral units or one valid side.", call. = FALSE)
  }

  spatial_cols <- if (dataset == "neuron_neuropil") {
    c("AnimalID", "region", "layer")
  } else {
    c("AnimalID", "region")
  }
  hemisphere <- aggregate_empirical_matrix_equal_weight(mat, metadata, spatial_cols)
  region <- if (dataset == "neuron_neuropil") {
    aggregate_empirical_matrix_equal_weight(hemisphere$mat, hemisphere$meta, c("AnimalID", "region"))
  } else {
    hemisphere
  }

  region_coverage <- region$meta |>
    dplyr::count(.data$AnimalID, name = "n_regions")
  if (any(region_coverage$n_regions != 4L)) {
    stop("Equal region weighting requires CA1, CA2, CA3, and DG for every AnimalID in ", dataset, ".", call. = FALSE)
  }
  animal <- aggregate_empirical_matrix_equal_weight(region$mat, region$meta, "AnimalID")
  animal_group <- metadata |>
    dplyr::distinct(.data$AnimalID, .data$ExpGroup)
  final_meta <- animal$meta |>
    dplyr::select("AnimalID") |>
    dplyr::left_join(animal_group, by = "AnimalID") |>
    dplyr::mutate(dataset = dataset)
  final_meta$observation_id <- paste(final_meta$AnimalID, final_meta$dataset, sep = "::")
  colnames(animal$mat) <- final_meta$observation_id
  rownames(final_meta) <- final_meta$observation_id

  one_sided <- spatial_qc$spatial_unit[spatial_qc$side_status == "one_sided"]
  aggregation_qc <- data.frame(
    dataset = dataset,
    n_samples = nrow(metadata),
    n_animals = dplyr::n_distinct(metadata$AnimalID),
    n_spatial_units = nrow(spatial_qc),
    n_bilateral_spatial_units = sum(spatial_qc$side_status == "bilateral"),
    n_one_sided_spatial_units = sum(spatial_qc$side_status == "one_sided"),
    n_invalid_spatial_units = sum(spatial_qc$side_status == "invalid"),
    one_sided_spatial_units = if (length(one_sided)) paste(sort(one_sided, method = "radix"), collapse = ";") else "",
    n_animal_region_units = ncol(region$mat),
    n_final_animal_dataset_rows = ncol(animal$mat),
    neuropil_layer_before_region_weighting = dataset == "neuron_neuropil",
    biological_replicate_unit = "AnimalID",
    hemisphere_role = "ReplicateGroup_Left_Right_technical_within_animal",
    stringsAsFactors = FALSE
  )
  list(
    mat = animal$mat,
    meta = final_meta,
    spatial_unit_qc = spatial_qc,
    aggregation_qc = aggregation_qc,
    hemisphere_mat = hemisphere$mat,
    hemisphere_meta = hemisphere$meta,
    region_mat = region$mat,
    region_meta = region$meta
  )
}

empirical_roi_normalize_member_set <- function(x) {
  vapply(as.character(x), function(value) {
    members <- trimws(unlist(strsplit(value, ";", fixed = TRUE), use.names = FALSE))
    members <- sort(unique(members[empirical_roi_nonempty(members)]), method = "radix")
    if (length(members)) paste(members, collapse = ";") else NA_character_
  }, character(1))
}

build_empirical_roi_protein_index <- function(feature_tables) {
  datasets <- empirical_roi_dataset_levels()
  if (!identical(sort(names(feature_tables)), sort(datasets))) {
    stop("Feature tables must contain exactly the three empirical ROI datasets.", call. = FALSE)
  }
  required <- c("ProteinGroupID", "member_accessions", "mapping_status")
  rows <- dplyr::bind_rows(lapply(datasets, function(dataset) {
    feature_table <- feature_tables[[dataset]]
    missing <- setdiff(required, names(feature_table))
    if (length(missing)) {
      stop("Canonical feature table for ", dataset, " is missing protein identity fields: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    accessions <- empirical_roi_normalize_member_set(feature_table$member_accessions)
    fully_mapped <- as.character(feature_table$mapping_status) == "mapped" & empirical_roi_nonempty(accessions)
    dplyr::tibble(
      dataset = dataset,
      canonical_ProteinGroupID = as.character(feature_table$ProteinGroupID),
      exact_member_accessions = accessions,
      protein_identity_alignment_eligible = fully_mapped,
      empirical_protein_group_key = dplyr::if_else(
        fully_mapped,
        paste0("EPG:exact_accession_set:", accessions),
        paste0("EPG:dataset_only:", dataset, ":", as.character(feature_table$ProteinGroupID))
      )
    )
  }))
  if (anyDuplicated(rows[c("dataset", "canonical_ProteinGroupID")])) {
    stop("Canonical ProteinGroupID is duplicated within a dataset feature table.", call. = FALSE)
  }
  duplicate_keys <- rows |>
    dplyr::count(.data$dataset, .data$empirical_protein_group_key, name = "n_groups_for_key")
  rows <- rows |>
    dplyr::left_join(duplicate_keys, by = c("dataset", "empirical_protein_group_key")) |>
    dplyr::group_by(.data$empirical_protein_group_key) |>
    dplyr::mutate(n_datasets_with_exact_group = dplyr::n_distinct(.data$dataset)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      cross_dataset_model_eligible = .data$protein_identity_alignment_eligible &
        .data$n_groups_for_key == 1L & .data$n_datasets_with_exact_group == length(datasets),
      protein_identity_level = "exact_fully_mapped_canonical_accession_set_no_gene_collapse"
    ) |>
    dplyr::arrange(.data$empirical_protein_group_key, .data$dataset, .data$canonical_ProteinGroupID)
  rows
}

align_empirical_roi_aggregated <- function(aggregated, protein_index) {
  datasets <- empirical_roi_dataset_levels()
  eligible_keys <- sort(unique(protein_index$empirical_protein_group_key[
    protein_index$cross_dataset_model_eligible %in% TRUE
  ]), method = "radix")
  if (!length(eligible_keys)) {
    stop("No exact canonical protein groups are uniquely present in all three compartments.", call. = FALSE)
  }
  setNames(lapply(datasets, function(dataset) {
    dataset_name <- dataset
    mapping <- protein_index |>
      dplyr::filter(.data$dataset == dataset_name, .data$cross_dataset_model_eligible) |>
      dplyr::arrange(match(.data$empirical_protein_group_key, eligible_keys))
    source <- aggregated[[dataset]]
    row_index <- match(mapping$canonical_ProteinGroupID, rownames(source$mat))
    if (anyNA(row_index) || !identical(mapping$empirical_protein_group_key, eligible_keys)) {
      stop("Protein-group alignment failed for ", dataset, ".", call. = FALSE)
    }
    source$mat <- source$mat[row_index, , drop = FALSE]
    rownames(source$mat) <- mapping$empirical_protein_group_key
    source
  }), datasets)
}

rekey_empirical_roi_table <- function(table, dataset, protein_index) {
  dataset_name <- dataset
  mapping <- protein_index |>
    dplyr::filter(.data$dataset == dataset_name, .data$cross_dataset_model_eligible) |>
    dplyr::select("canonical_ProteinGroupID", "empirical_protein_group_key")
  table |>
    dplyr::rename(canonical_ProteinGroupID = "ProteinGroupID") |>
    dplyr::inner_join(mapping, by = "canonical_ProteinGroupID") |>
    dplyr::mutate(ProteinGroupID = .data$empirical_protein_group_key) |>
    dplyr::select(-"empirical_protein_group_key")
}

map_empirical_roi_raw_detection <- function(protein_index, global_features, raw_detection) {
  required_features <- c("ProteinGroupID", "member_accessions", "mapping_status")
  required_detection <- c("ProteinGroupID", "dataset", "observed_detection_rate_raw")
  missing <- c(
    setdiff(required_features, names(global_features)),
    setdiff(required_detection, names(raw_detection))
  )
  if (length(missing)) {
    stop("Global raw-detection provenance is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  global_index <- global_features |>
    dplyr::transmute(
      global_ProteinGroupID = as.character(.data$ProteinGroupID),
      exact_member_accessions = empirical_roi_normalize_member_set(.data$member_accessions),
      global_mapping_status = as.character(.data$mapping_status)
    ) |>
    dplyr::filter(
      .data$global_mapping_status == "mapped",
      empirical_roi_nonempty(.data$exact_member_accessions)
    )
  if (anyDuplicated(global_index$exact_member_accessions)) {
    stop("Global raw-detection feature table has duplicated exact accession sets.", call. = FALSE)
  }
  keys <- protein_index |>
    dplyr::filter(.data$cross_dataset_model_eligible) |>
    dplyr::distinct(.data$empirical_protein_group_key, .data$exact_member_accessions)
  if (anyDuplicated(keys$empirical_protein_group_key) || anyDuplicated(keys$exact_member_accessions)) {
    stop("Eligible empirical protein identity is not one-to-one with exact accession sets.", call. = FALSE)
  }
  tidyr::crossing(keys, dataset = empirical_roi_dataset_levels()) |>
    dplyr::left_join(global_index, by = "exact_member_accessions") |>
    dplyr::left_join(
      dplyr::transmute(
        raw_detection,
        global_ProteinGroupID = as.character(.data$ProteinGroupID),
        dataset = as.character(.data$dataset),
        observed_detection_rate_raw = suppressWarnings(as.numeric(.data$observed_detection_rate_raw))
      ),
      by = c("global_ProteinGroupID", "dataset")
    ) |>
    dplyr::mutate(
      ProteinGroupID = .data$empirical_protein_group_key,
      detection_rate = .data$observed_detection_rate_raw,
      detectability_source = dplyr::if_else(
        is.finite(.data$observed_detection_rate_raw),
        "raw_unified_pre_imputation_exact_accession_bridge",
        "raw_unified_exact_accession_not_available"
      )
    ) |>
    dplyr::select(
      "ProteinGroupID", "dataset", "observed_detection_rate_raw", "detection_rate",
      "detectability_source", "global_ProteinGroupID", "exact_member_accessions"
    ) |>
    dplyr::arrange(.data$ProteinGroupID, .data$dataset)
}

read_empirical_roi_raw_detection <- function(
    protein_index,
    root = file.path("data", "processed", "01_preprocessing", "joint_compartment_qc", "global")) {
  feature_file <- file.path(root, "canonical_feature_table.csv")
  detection_file <- file.path(root, "observed_detection_by_dataset.csv")
  if (!file.exists(feature_file) || !file.exists(detection_file)) {
    stop("Required global raw-detection provenance files are missing under: ", root, call. = FALSE)
  }
  global_features <- readr::read_csv(
    feature_file,
    col_select = dplyr::all_of(c("ProteinGroupID", "member_accessions", "mapping_status")),
    show_col_types = FALSE, progress = FALSE
  )
  raw_detection <- readr::read_csv(
    detection_file,
    col_select = dplyr::all_of(c("ProteinGroupID", "dataset", "observed_detection_rate_raw")),
    show_col_types = FALSE, progress = FALSE
  )
  map_empirical_roi_raw_detection(protein_index, global_features, raw_detection)
}

build_empirical_roi_model_input <- function(aggregated) {
  datasets <- empirical_roi_dataset_levels()
  missing <- setdiff(datasets, names(aggregated))
  if (length(missing)) stop("Aggregated model input is missing datasets: ", paste(missing, collapse = ", "), call. = FALSE)
  protein_ids <- sort(unique(unlist(lapply(aggregated[datasets], function(x) rownames(x$mat)))), method = "radix")
  observation_meta <- dplyr::bind_rows(lapply(datasets, function(dataset) aggregated[[dataset]]$meta)) |>
    dplyr::mutate(
      AnimalID = as.character(.data$AnimalID),
      dataset = factor(as.character(.data$dataset), levels = datasets)
    ) |>
    dplyr::arrange(.data$AnimalID, .data$dataset)
  if (anyDuplicated(observation_meta$observation_id)) {
    stop("Final AnimalID x dataset observations are not unique.", call. = FALSE)
  }
  expression <- matrix(
    NA_real_, nrow = length(protein_ids), ncol = nrow(observation_meta),
    dimnames = list(protein_ids, observation_meta$observation_id)
  )
  for (dataset in datasets) {
    source <- aggregated[[dataset]]$mat
    expression[rownames(source), colnames(source)] <- source
  }
  observation_meta$AnimalID <- factor(
    observation_meta$AnimalID,
    levels = sort(unique(as.character(observation_meta$AnimalID)), method = "radix")
  )
  design <- stats::model.matrix(~ AnimalID + dataset, data = observation_meta)
  if (qr(design)$rank != ncol(design)) {
    stop("AnimalID + dataset design is not full rank.", call. = FALSE)
  }
  list(expression = expression, observation_meta = observation_meta, design = design)
}

build_empirical_roi_limma_inference <- function(model_input) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for empirical ROI marker inference.", call. = FALSE)
  }
  expression <- model_input$expression
  observation_meta <- model_input$observation_meta
  design <- model_input$design
  model_warnings <- character()
  capture_warning <- function(w) {
    model_warnings <<- unique(c(model_warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
  fit <- withCallingHandlers(limma::lmFit(expression, design), warning = capture_warning)
  contrast_matrix <- limma::makeContrasts(
    microglia_vs_neuron_neuropil = datasetmicroglia,
    microglia_vs_neuron_soma = datasetmicroglia - datasetneuron_soma,
    neuron_neuropil_vs_microglia = -datasetmicroglia,
    neuron_soma_vs_microglia = datasetneuron_soma - datasetmicroglia,
    levels = design
  )
  contrasted <- withCallingHandlers(limma::contrasts.fit(fit, contrast_matrix), warning = capture_warning)
  moderated <- withCallingHandlers(limma::eBayes(contrasted), warning = capture_warning)
  coefficients <- moderated$coefficients
  p_values <- moderated$p.value
  fdr <- apply(p_values, 2L, stats::p.adjust, method = "BH")
  if (is.null(dim(fdr))) fdr <- matrix(fdr, ncol = ncol(p_values), dimnames = dimnames(p_values))

  finite <- is.finite(expression) & !is.na(expression)
  datasets <- empirical_roi_dataset_levels()
  observation_counts <- vapply(datasets, function(dataset) {
    rowSums(finite[, as.character(observation_meta$dataset) == dataset, drop = FALSE])
  }, numeric(nrow(expression)))
  colnames(observation_counts) <- datasets
  animal_ids <- as.character(observation_meta$AnimalID)
  n_animals <- apply(finite, 1L, function(ok) length(unique(animal_ids[ok])))
  finite_coefficients <- rowSums(is.finite(coefficients)) == ncol(coefficients)
  finite_p_values <- rowSums(is.finite(p_values)) == ncol(p_values)
  finite_fdr <- rowSums(is.finite(fdr)) == ncol(fdr)
  testable_all <- rowSums(observation_counts > 0L) == length(datasets) & fit$df.residual > 0L

  stats <- data.frame(
    ProteinGroupID = rownames(expression),
    n_observations = rowSums(finite),
    n_animals = as.integer(n_animals),
    n_observations_neuron_neuropil = as.integer(observation_counts[, "neuron_neuropil"]),
    n_observations_neuron_soma = as.integer(observation_counts[, "neuron_soma"]),
    n_observations_microglia = as.integer(observation_counts[, "microglia"]),
    residual_df = as.numeric(fit$df.residual),
    nominal_design_rank = qr(design)$rank,
    nominal_residual_df = nrow(design) - qr(design)$rank,
    testable_in_all_three_compartments = testable_all,
    finite_coefficients_all_contrasts = finite_coefficients,
    finite_p_values_all_contrasts = finite_p_values,
    finite_FDR_all_contrasts = finite_fdr,
    model_fit_valid = testable_all & finite_coefficients & finite_p_values & finite_fdr,
    logFC_microglia_vs_neuropil = coefficients[, "microglia_vs_neuron_neuropil"],
    logFC_microglia_vs_soma = coefficients[, "microglia_vs_neuron_soma"],
    logFC_neuropil_vs_microglia = coefficients[, "neuron_neuropil_vs_microglia"],
    logFC_soma_vs_microglia = coefficients[, "neuron_soma_vs_microglia"],
    p_value_microglia_vs_neuropil = p_values[, "microglia_vs_neuron_neuropil"],
    p_value_microglia_vs_soma = p_values[, "microglia_vs_neuron_soma"],
    p_value_neuropil_vs_microglia = p_values[, "neuron_neuropil_vs_microglia"],
    p_value_soma_vs_microglia = p_values[, "neuron_soma_vs_microglia"],
    FDR_microglia_vs_neuropil = fdr[, "microglia_vs_neuron_neuropil"],
    FDR_microglia_vs_soma = fdr[, "microglia_vs_neuron_soma"],
    FDR_neuropil_vs_microglia = fdr[, "neuron_neuropil_vs_microglia"],
    FDR_soma_vs_microglia = fdr[, "neuron_soma_vs_microglia"],
    biological_replicate_unit = "AnimalID",
    model_description = "limma eBayes paired fixed-effect design: ~ AnimalID + dataset after hemisphere/layer/region aggregation",
    marker_contract_version = empirical_roi_marker_contract_version(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    stats = stats,
    design = design,
    contrasts = contrast_matrix,
    model_warnings = model_warnings
  )
}

empirical_roi_mapping_provenance <- function(feature_tables, protein_index) {
  if (is.null(names(feature_tables)) || any(!nzchar(names(feature_tables)))) {
    stop("Feature-table list must be named by dataset.", call. = FALSE)
  }
  required <- c(
    "ProteinGroupID", "original_identifier", "member_accessions", "official_gene_symbol",
    "protein_group_gene_annotation_status", "gene_level_claim_allowed",
    "gene_annotation_contract_version", "uniprot_mapping_file_hash", "orgdb_package_version"
  )
  rows <- dplyr::bind_rows(lapply(names(feature_tables), function(dataset) {
    dataset_name <- dataset
    feature_table <- feature_tables[[dataset]]
    missing <- setdiff(required, names(feature_table))
    if (length(missing)) {
      stop("Canonical feature table for ", dataset, " is missing mapping provenance: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    index_rows <- protein_index |>
      dplyr::filter(.data$dataset == dataset_name, .data$cross_dataset_model_eligible) |>
      dplyr::select("canonical_ProteinGroupID", "empirical_protein_group_key",
                    "protein_identity_level", "exact_member_accessions")
    feature_table |>
      dplyr::transmute(
        dataset = dataset,
        canonical_ProteinGroupID = as.character(.data$ProteinGroupID),
        original_identifier = as.character(.data$original_identifier),
        member_accessions = as.character(.data$member_accessions),
        official_gene_symbol = as.character(.data$official_gene_symbol),
        protein_group_gene_annotation_status = as.character(.data$protein_group_gene_annotation_status),
        gene_level_claim_allowed = .data$gene_level_claim_allowed %in% TRUE,
        gene_annotation_contract_version = as.character(.data$gene_annotation_contract_version),
        uniprot_mapping_file_hash = as.character(.data$uniprot_mapping_file_hash),
        orgdb_package_version = as.character(.data$orgdb_package_version)
      ) |>
      dplyr::inner_join(index_rows, by = "canonical_ProteinGroupID")
  }))
  rows |>
    dplyr::group_by(.data$empirical_protein_group_key) |>
    dplyr::summarise(
      ProteinGroupID = dplyr::first(.data$empirical_protein_group_key),
      ProteinGroupID_neuron_neuropil = empirical_roi_collapse_unique(.data$canonical_ProteinGroupID[.data$dataset == "neuron_neuropil"]),
      ProteinGroupID_neuron_soma = empirical_roi_collapse_unique(.data$canonical_ProteinGroupID[.data$dataset == "neuron_soma"]),
      ProteinGroupID_microglia = empirical_roi_collapse_unique(.data$canonical_ProteinGroupID[.data$dataset == "microglia"]),
      canonical_ProteinGroupIDs = empirical_roi_collapse_unique(.data$canonical_ProteinGroupID),
      protein_identity_level = dplyr::first(.data$protein_identity_level),
      exact_member_accessions = dplyr::first(.data$exact_member_accessions),
      mapping_source_datasets = empirical_roi_collapse_unique(.data$dataset),
      original_identifiers = empirical_roi_collapse_unique(.data$original_identifier),
      member_accessions = empirical_roi_collapse_unique(.data$member_accessions),
      official_gene_symbols = empirical_roi_collapse_unique(.data$official_gene_symbol),
      n_official_gene_symbols = dplyr::n_distinct(.data$official_gene_symbol[empirical_roi_nonempty(.data$official_gene_symbol)]),
      protein_group_gene_annotation_status = empirical_roi_collapse_unique(.data$protein_group_gene_annotation_status),
      mapping_claim_allowed_all_sources = dplyr::n_distinct(.data$dataset) == 3L && all(.data$gene_level_claim_allowed),
      gene_annotation_contract_version = empirical_roi_collapse_unique(.data$gene_annotation_contract_version),
      uniprot_mapping_file_hash = empirical_roi_collapse_unique(.data$uniprot_mapping_file_hash),
      orgdb_package_version = empirical_roi_collapse_unique(.data$orgdb_package_version),
      .groups = "drop"
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      official_gene_symbol = dplyr::if_else(
        .data$n_official_gene_symbols == 1L,
        .data$official_gene_symbols,
        NA_character_
      ),
      gene_mapping_claim_allowed = .data$mapping_claim_allowed_all_sources &
        .data$n_official_gene_symbols == 1L &
        .data$protein_group_gene_annotation_status == "concordant_official_gene",
      mapping_provenance_status = dplyr::case_when(
        .data$gene_mapping_claim_allowed ~ "canonical_gene_claim_allowed",
        .data$n_official_gene_symbols > 1L ~ "conflicting_official_gene_symbols",
        .data$n_official_gene_symbols == 0L ~ "no_official_gene_symbol",
        TRUE ~ "canonical_gene_claim_not_allowed"
      )
    )
}

empirical_roi_fdr_pass <- function(x, threshold) {
  is.finite(x) & !is.na(x) & x <= threshold
}

select_empirical_roi_marker_sets <- function(stats, min_detection = 0.30,
                                             min_abs_logfc = 0.50,
                                             fdr_threshold = 0.10) {
  required <- c(
    "ProteinGroupID", "official_gene_symbol", "gene_mapping_claim_allowed",
    "detection_rate_microglia", "detection_rate_neuropil", "detection_rate_soma",
    "logFC_microglia_vs_neuropil", "logFC_microglia_vs_soma",
    "logFC_neuropil_vs_microglia", "logFC_soma_vs_microglia",
    "p_value_microglia_vs_neuropil", "p_value_microglia_vs_soma",
    "p_value_neuropil_vs_microglia", "p_value_soma_vs_microglia",
    "FDR_microglia_vs_neuropil", "FDR_microglia_vs_soma",
    "FDR_neuropil_vs_microglia", "FDR_soma_vs_microglia",
    "residual_df"
  )
  missing <- setdiff(required, names(stats))
  if (length(missing)) stop("Marker selection input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  stats <- stats |>
    dplyr::mutate(
      dataset_enriched_in = dplyr::case_when(
        .data$logFC_microglia_vs_neuropil >= min_abs_logfc &
          .data$logFC_microglia_vs_soma >= min_abs_logfc ~ "microglia",
        .data$logFC_neuropil_vs_microglia >= min_abs_logfc ~ "neuron_neuropil",
        .data$logFC_soma_vs_microglia >= min_abs_logfc ~ "neuron_soma",
        TRUE ~ "shared_or_ambiguous"
      )
    )

  strong_detection <- min(0.80, min_detection + 0.20)
  strong_logfc <- min_abs_logfc + 0.50
  valid_micro_neuropil <- is.finite(stats$logFC_microglia_vs_neuropil) &
    empirical_roi_fdr_pass(stats$FDR_microglia_vs_neuropil, fdr_threshold) &
    is.finite(stats$p_value_microglia_vs_neuropil) & stats$residual_df > 0
  valid_micro_soma <- is.finite(stats$logFC_microglia_vs_soma) &
    empirical_roi_fdr_pass(stats$FDR_microglia_vs_soma, fdr_threshold) &
    is.finite(stats$p_value_microglia_vs_soma) & stats$residual_df > 0
  valid_neuropil_micro <- is.finite(stats$logFC_neuropil_vs_microglia) &
    empirical_roi_fdr_pass(stats$FDR_neuropil_vs_microglia, fdr_threshold) &
    is.finite(stats$p_value_neuropil_vs_microglia) & stats$residual_df > 0
  valid_soma_micro <- is.finite(stats$logFC_soma_vs_microglia) &
    empirical_roi_fdr_pass(stats$FDR_soma_vs_microglia, fdr_threshold) &
    is.finite(stats$p_value_soma_vs_microglia) & stats$residual_df > 0

  make_set <- function(marker_set, keep, comparison, confidence, notes,
                       evidence_type, inferential_valid, p_value, fdr, selection_rule) {
    tab <- stats[which(keep %in% TRUE), , drop = FALSE]
    valid <- as.logical(inferential_valid[which(keep %in% TRUE)])
    is_inferential <- identical(evidence_type, "inferential")
    claim_allowed <- is_inferential & valid & tab$gene_mapping_claim_allowed %in% TRUE
    tab |>
      dplyr::mutate(
        marker_set = marker_set,
        ProteinID = .data$ProteinGroupID,
        mapped_gene_symbol = .data$official_gene_symbol,
        GeneSymbol = dplyr::if_else(claim_allowed, .data$official_gene_symbol, NA_character_),
        comparison = comparison,
        p_value = as.numeric(p_value[which(keep %in% TRUE)]),
        FDR = as.numeric(fdr[which(keep %in% TRUE)]),
        FDR_pass = if (is_inferential) empirical_roi_fdr_pass(.data$FDR, fdr_threshold) else NA,
        marker_source = "empirical_roi_marker_sets",
        selection_rule = selection_rule,
        confidence = confidence,
        notes = notes,
        model_type = "limma_paired_fixed_effect",
        marker_evidence_type = evidence_type,
        inferential_test_valid = if (is_inferential) valid else FALSE,
        claim_allowed = claim_allowed,
        claim_allowed_reason = dplyr::case_when(
          evidence_type == "descriptive_only" ~ "descriptive_only_no_equivalence_test",
          !valid ~ "required_contrast_inference_not_valid",
          !.data$gene_mapping_claim_allowed ~ "canonical_gene_mapping_claim_not_allowed",
          TRUE ~ "inferential_and_canonical_mapping_supported"
        ),
        marker_contract_version = empirical_roi_marker_contract_version()
      ) |>
      dplyr::relocate(
        "marker_set", "ProteinGroupID", "ProteinID", "GeneSymbol",
        "mapped_gene_symbol", "marker_evidence_type", "inferential_test_valid",
        "claim_allowed", "claim_allowed_reason", "dataset_enriched_in",
        "comparison", "p_value", "FDR", "FDR_pass"
      )
  }

  micro_keep <- is.finite(stats$detection_rate_microglia) &
    stats$detection_rate_microglia >= min_detection &
    stats$logFC_microglia_vs_neuropil >= min_abs_logfc &
    stats$logFC_microglia_vs_soma >= min_abs_logfc &
    valid_micro_neuropil & valid_micro_soma
  high_micro_keep <- is.finite(stats$detection_rate_microglia) &
    stats$detection_rate_microglia >= strong_detection &
    stats$logFC_microglia_vs_neuropil >= strong_logfc &
    stats$logFC_microglia_vs_soma >= strong_logfc &
    valid_micro_neuropil & valid_micro_soma
  neuropil_keep <- is.finite(stats$detection_rate_neuropil) &
    stats$detection_rate_neuropil >= min_detection &
    stats$logFC_neuropil_vs_microglia >= min_abs_logfc & valid_neuropil_micro
  high_neuropil_keep <- is.finite(stats$detection_rate_neuropil) &
    stats$detection_rate_neuropil >= strong_detection &
    stats$logFC_neuropil_vs_microglia >= strong_logfc & valid_neuropil_micro
  soma_keep <- is.finite(stats$detection_rate_soma) &
    stats$detection_rate_soma >= min_detection &
    stats$logFC_soma_vs_microglia >= min_abs_logfc & valid_soma_micro
  shared_keep <- is.finite(stats$detection_rate_microglia) &
    is.finite(stats$detection_rate_neuropil) &
    stats$detection_rate_microglia >= min_detection &
    stats$detection_rate_neuropil >= min_detection &
    is.finite(stats$logFC_microglia_vs_neuropil) &
    abs(stats$logFC_microglia_vs_neuropil) < min_abs_logfc &
    stats$dataset_enriched_in == "shared_or_ambiguous"

  out <- dplyr::bind_rows(
    make_set(
      "empirical_microglia_roi_enriched", micro_keep,
      "microglia_vs_neuron_neuropil_and_soma", "empirical_roi_enriched",
      "ROI enrichment only; microglia ROI samples are not purified microglia.",
      "inferential", valid_micro_neuropil & valid_micro_soma,
      pmax(stats$p_value_microglia_vs_neuropil, stats$p_value_microglia_vs_soma),
      pmax(stats$FDR_microglia_vs_neuropil, stats$FDR_microglia_vs_soma),
      paste0("detection_microglia>=", min_detection, "; both logFC>=", min_abs_logfc, "; both BH FDR<=", fdr_threshold)
    ),
    make_set(
      "empirical_microglia_roi_high_confidence", high_micro_keep,
      "microglia_vs_neuron_neuropil_and_soma", "empirical_high_confidence",
      "Stronger empirical ROI enrichment; annotation only.",
      "inferential", valid_micro_neuropil & valid_micro_soma,
      pmax(stats$p_value_microglia_vs_neuropil, stats$p_value_microglia_vs_soma),
      pmax(stats$FDR_microglia_vs_neuropil, stats$FDR_microglia_vs_soma),
      paste0("detection_microglia>=", strong_detection, "; both logFC>=", strong_logfc, "; both BH FDR<=", fdr_threshold)
    ),
    make_set(
      "empirical_neuropil_enriched", neuropil_keep,
      "neuron_neuropil_vs_microglia", "empirical_roi_enriched",
      "Neuropil-sensitive marker evidence; do not subtract from WGCNA.",
      "inferential", valid_neuropil_micro,
      stats$p_value_neuropil_vs_microglia, stats$FDR_neuropil_vs_microglia,
      paste0("detection_neuropil>=", min_detection, "; logFC>=", min_abs_logfc, "; BH FDR<=", fdr_threshold)
    ),
    make_set(
      "empirical_neuropil_sensitive_high_confidence", high_neuropil_keep,
      "neuron_neuropil_vs_microglia", "empirical_high_confidence",
      "High-confidence neuropil-sensitive marker evidence.",
      "inferential", valid_neuropil_micro,
      stats$p_value_neuropil_vs_microglia, stats$FDR_neuropil_vs_microglia,
      paste0("detection_neuropil>=", strong_detection, "; logFC>=", strong_logfc, "; BH FDR<=", fdr_threshold)
    ),
    make_set(
      "empirical_soma_enriched", soma_keep,
      "neuron_soma_vs_microglia", "empirical_roi_enriched",
      "Soma-enriched marker evidence; annotation only.",
      "inferential", valid_soma_micro,
      stats$p_value_soma_vs_microglia, stats$FDR_soma_vs_microglia,
      paste0("detection_soma>=", min_detection, "; logFC>=", min_abs_logfc, "; BH FDR<=", fdr_threshold)
    ),
    make_set(
      "empirical_microglia_neuropil_shared", shared_keep,
      "microglia_neuropil_shared_descriptive", "empirical_shared_descriptive",
      "Descriptive co-detection/effect-size set only; no equivalence test was performed.",
      "descriptive_only", rep(FALSE, nrow(stats)),
      rep(NA_real_, nrow(stats)), rep(NA_real_, nrow(stats)),
      paste0("detection_microglia_and_neuropil>=", min_detection, "; abs(logFC_microglia_vs_neuropil)<", min_abs_logfc, "; descriptive only")
    )
  )
  if (!nrow(out)) return(out)
  out |>
    dplyr::arrange(.data$marker_set, .data$ProteinGroupID)
}
