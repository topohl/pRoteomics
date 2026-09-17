#!/usr/bin/env Rscript
# ================================================================
# Script: 03_qc_exploration/05_empirical_roi_marker_discovery.r
# Stage: qc_global
# Scope: global
# Consumes: required Stage 01 post-filter/imputed quantitative matrices for all datasets, mouse UniProt mapping, sample metadata; optional manual mappings.
# Produces: results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Global because it compares all three dataset matrices.
# ================================================================

# Discover experiment-specific ROI-enrichment marker sets across neuron
# neuropil, neuron soma, and microglia at the AnimalID biological-replicate
# level. Structured metadata defines anatomy; raw acquisition names do not.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "empirical_roi_marker_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()
proposed_only <- "--proposed-only" %in% args ||
  tolower(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_PROPOSED_ONLY", unset = "false")) == "true"
DATASETS <- empirical_roi_dataset_levels()
SUBSTEP_ID <- "05_empirical_roi_marker_discovery"
PATHS <- create_module_dirs("03_qc_exploration", SUBSTEP_ID)

min_detection <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_MIN_DETECTION", unset = "0.30")))
min_abs_logfc <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_MIN_ABS_LOGFC", unset = "0.50")))
fdr_threshold <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_FDR", unset = "0.10")))
if (!is.finite(min_detection)) min_detection <- 0.30
if (!is.finite(min_abs_logfc)) min_abs_logfc <- 0.50
if (!is.finite(fdr_threshold)) fdr_threshold <- 0.10

inputs <- setNames(lapply(DATASETS, function(dataset) {
  list(
    matrix = qc_resolve_matrix(
      dataset,
      env = paste0("PROTEOMICS_EMPIRICAL_MARKER_", toupper(dataset), "_MATRIX_FILE")
    ),
    metadata = qc_resolve_metadata(
      dataset,
      env = paste0("PROTEOMICS_EMPIRICAL_MARKER_", toupper(dataset), "_METADATA_FILE")
    )
  )
}), DATASETS)

canonical_marker_path <- file.path(PATHS$tables, "empirical_roi_marker_sets.csv")
proposed_marker_path <- file.path(PATHS$tables, "empirical_roi_marker_sets_proposed.csv")
selected_marker_path <- if (proposed_only) proposed_marker_path else canonical_marker_path

if (dry_run) {
  dry_run_line("Script", "03_qc_exploration/05_empirical_roi_marker_discovery.r")
  dry_run_line("Mode", if (proposed_only) "proposed-only (canonical output protected)" else "canonical")
  for (dataset in DATASETS) {
    dry_run_line(
      paste(dataset, "matrix"), inputs[[dataset]]$matrix,
      if (file.exists(inputs[[dataset]]$matrix)) "PASS" else "WARN"
    )
    dry_run_line(
      paste(dataset, "metadata"), inputs[[dataset]]$metadata,
      if (file.exists(inputs[[dataset]]$metadata)) "PASS" else "WARN"
    )
  }
  dry_run_line("Output", selected_marker_path, "INFO")
  quit(status = 0, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "readr", "limma")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

mapping_context <- qc_build_mapping_context()

read_dataset <- function(dataset) {
  if (!file.exists(inputs[[dataset]]$matrix)) {
    stop("Missing empirical marker matrix for ", dataset, ": ", inputs[[dataset]]$matrix, call. = FALSE)
  }
  expr <- qc_load_canonical_expression(
    inputs[[dataset]]$matrix,
    inputs[[dataset]]$metadata,
    dataset,
    mapping_context = mapping_context
  )
  if (!isTRUE(proposed_only)) {
    qc_write_canonical_feature_artifacts(expr, file.path(PATHS$tables, dataset))
  }
  aggregated <- aggregate_empirical_roi_dataset(expr$mat, expr$meta, dataset)
  feature_stats <- data.frame(
    dataset = dataset,
    ProteinGroupID = expr$feature_table$ProteinGroupID,
    availability_after_filtering = rowMeans(is.finite(expr$mat) & !is.na(expr$mat)),
    abundance_score_imputed = rowMeans(expr$mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  animal_means <- data.frame(
    ProteinGroupID = rownames(aggregated$mat),
    dataset = dataset,
    mean_abundance = rowMeans(aggregated$mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(
    expr = expr,
    aggregated = aggregated,
    feature_stats = feature_stats,
    animal_means = animal_means
  )
}

loaded <- setNames(lapply(DATASETS, read_dataset), DATASETS)
all_metadata <- dplyr::bind_rows(lapply(DATASETS, function(dataset) {
  loaded[[dataset]]$expr$meta |>
    dplyr::select(dplyr::all_of(empirical_roi_required_metadata_fields()))
})) |>
  dplyr::distinct(.data$sample_id, .keep_all = TRUE)
metadata_correction_audit <- validate_known_empirical_metadata_corrections(all_metadata, require_all = TRUE)

aggregated_unaligned <- setNames(lapply(DATASETS, function(dataset) loaded[[dataset]]$aggregated), DATASETS)
feature_tables <- setNames(lapply(DATASETS, function(dataset) loaded[[dataset]]$expr$feature_table), DATASETS)
protein_index <- build_empirical_roi_protein_index(feature_tables)
raw_detection_provenance <- read_empirical_roi_raw_detection(protein_index)
aggregated <- align_empirical_roi_aggregated(aggregated_unaligned, protein_index)
model_input <- build_empirical_roi_model_input(aggregated)
if (ncol(model_input$expression) != 27L) {
  stop("Expected exactly 27 AnimalID x dataset observations; found ", ncol(model_input$expression), ".", call. = FALSE)
}
model <- build_empirical_roi_limma_inference(model_input)

mapping_provenance <- empirical_roi_mapping_provenance(feature_tables, protein_index)

for (dataset in DATASETS) {
  loaded[[dataset]]$feature_stats <- rekey_empirical_roi_table(
    loaded[[dataset]]$feature_stats, dataset, protein_index
  )
  loaded[[dataset]]$animal_means <- rekey_empirical_roi_table(
    loaded[[dataset]]$animal_means, dataset, protein_index
  )
}

feature_stats_wide <- dplyr::bind_rows(lapply(loaded, `[[`, "feature_stats")) |>
  dplyr::select(
    "ProteinGroupID", "dataset", "availability_after_filtering", "abundance_score_imputed"
  ) |>
  tidyr::pivot_wider(
    id_cols = "ProteinGroupID",
    names_from = "dataset",
    values_from = c("availability_after_filtering", "abundance_score_imputed")
  )
raw_detection_wide <- raw_detection_provenance |>
  dplyr::select(
    "ProteinGroupID", "dataset", "detection_rate", "observed_detection_rate_raw", "detectability_source"
  ) |>
  tidyr::pivot_wider(
    id_cols = "ProteinGroupID", names_from = "dataset",
    values_from = c("detection_rate", "observed_detection_rate_raw", "detectability_source")
  )
animal_means_wide <- dplyr::bind_rows(lapply(loaded, `[[`, "animal_means")) |>
  dplyr::select("ProteinGroupID", "dataset", "mean_abundance") |>
  tidyr::pivot_wider(
    id_cols = "ProteinGroupID", names_from = "dataset", values_from = "mean_abundance",
    names_prefix = "mean_abundance_"
  )

protein_stats <- model$stats |>
  dplyr::left_join(mapping_provenance, by = "ProteinGroupID") |>
  dplyr::left_join(feature_stats_wide, by = "ProteinGroupID") |>
  dplyr::left_join(raw_detection_wide, by = "ProteinGroupID") |>
  dplyr::left_join(animal_means_wide, by = "ProteinGroupID") |>
  dplyr::mutate(
    detection_rate_microglia = .data$detection_rate_microglia,
    detection_rate_neuropil = .data$detection_rate_neuron_neuropil,
    detection_rate_soma = .data$detection_rate_neuron_soma,
    observed_detection_rate_raw_microglia = .data$observed_detection_rate_raw_microglia,
    observed_detection_rate_raw_neuropil = .data$observed_detection_rate_raw_neuron_neuropil,
    observed_detection_rate_raw_soma = .data$observed_detection_rate_raw_neuron_soma,
    availability_after_filtering_microglia = .data$availability_after_filtering_microglia,
    availability_after_filtering_neuropil = .data$availability_after_filtering_neuron_neuropil,
    availability_after_filtering_soma = .data$availability_after_filtering_neuron_soma,
    abundance_score_imputed_microglia = .data$abundance_score_imputed_microglia,
    abundance_score_imputed_neuropil = .data$abundance_score_imputed_neuron_neuropil,
    abundance_score_imputed_soma = .data$abundance_score_imputed_neuron_soma,
    mean_abundance_microglia = .data$mean_abundance_microglia,
    mean_abundance_neuropil = .data$mean_abundance_neuron_neuropil,
    mean_abundance_soma = .data$mean_abundance_neuron_soma
  )

out <- select_empirical_roi_marker_sets(
  protein_stats,
  min_detection = min_detection,
  min_abs_logfc = min_abs_logfc,
  fdr_threshold = fdr_threshold
)

aggregation_qc <- dplyr::bind_rows(lapply(aggregated, `[[`, "aggregation_qc"))
spatial_unit_qc <- dplyr::bind_rows(lapply(aggregated, `[[`, "spatial_unit_qc"))
residual_df_distribution <- model$stats |>
  dplyr::count(.data$residual_df, name = "n_ProteinGroupIDs") |>
  dplyr::arrange(.data$residual_df)
inference_validation <- data.frame(
  metric = c(
    "final_animal_x_dataset_observations",
    "design_rank",
    "nominal_complete_residual_df",
    "canonical_ProteinGroupID_rows_across_datasets",
    "exact_protein_groups_uniquely_present_all_three_compartments",
    "protein_identity_rows_not_cross_dataset_model_eligible",
    "ProteinGroupIDs_total",
    "ProteinGroupIDs_testable_all_three_compartments",
    "ProteinGroupIDs_finite_coefficients_all_contrasts",
    "ProteinGroupIDs_finite_p_values_all_contrasts",
    "ProteinGroupIDs_finite_BH_FDR_all_contrasts",
    "ProteinGroupIDs_complete_27_observations",
    "model_warning_count"
  ),
  value = as.character(c(
    ncol(model_input$expression),
    qr(model_input$design)$rank,
    nrow(model_input$design) - qr(model_input$design)$rank,
    nrow(protein_index),
    dplyr::n_distinct(protein_index$empirical_protein_group_key[protein_index$cross_dataset_model_eligible]),
    sum(!protein_index$cross_dataset_model_eligible),
    nrow(model$stats),
    sum(model$stats$testable_in_all_three_compartments),
    sum(model$stats$finite_coefficients_all_contrasts),
    sum(model$stats$finite_p_values_all_contrasts),
    sum(model$stats$finite_FDR_all_contrasts),
    sum(model$stats$n_observations == 27L),
    length(model$model_warnings)
  )),
  stringsAsFactors = FALSE
)
model_warning_table <- data.frame(
  warning = model$model_warnings,
  stringsAsFactors = FALSE
)
marker_counts <- out |>
  dplyr::count(.data$marker_set, .data$marker_evidence_type, .data$claim_allowed, name = "n_rows") |>
  dplyr::arrange(.data$marker_set, .data$marker_evidence_type, dplyr::desc(.data$claim_allowed))

if (isTRUE(proposed_only)) {
  proposed_outputs <- list(
    markers = proposed_marker_path,
    protein_inference = file.path(PATHS$tables, "empirical_roi_protein_inference_proposed.csv"),
    aggregation_qc = file.path(PATHS$tables, "empirical_roi_aggregation_qc_proposed.csv"),
    spatial_unit_qc = file.path(PATHS$tables, "empirical_roi_spatial_unit_qc_proposed.csv"),
    metadata_corrections = file.path(PATHS$tables, "empirical_roi_metadata_corrections_proposed.csv"),
    residual_df = file.path(PATHS$tables, "empirical_roi_residual_df_distribution_proposed.csv"),
    inference_validation = file.path(PATHS$tables, "empirical_roi_inference_validation_proposed.csv"),
    marker_counts = file.path(PATHS$tables, "empirical_roi_marker_counts_proposed.csv"),
    protein_identity_crosswalk = file.path(PATHS$tables, "empirical_roi_protein_identity_crosswalk_proposed.csv"),
    raw_detection_provenance = file.path(PATHS$tables, "empirical_roi_raw_detection_provenance_proposed.csv"),
    model_warnings = file.path(PATHS$tables, "empirical_roi_model_warnings_proposed.csv")
  )
  readr::write_csv(out, proposed_outputs$markers, na = "")
  readr::write_csv(protein_stats, proposed_outputs$protein_inference, na = "")
  readr::write_csv(aggregation_qc, proposed_outputs$aggregation_qc, na = "")
  readr::write_csv(spatial_unit_qc, proposed_outputs$spatial_unit_qc, na = "")
  readr::write_csv(metadata_correction_audit, proposed_outputs$metadata_corrections, na = "")
  readr::write_csv(residual_df_distribution, proposed_outputs$residual_df, na = "")
  readr::write_csv(inference_validation, proposed_outputs$inference_validation, na = "")
  readr::write_csv(marker_counts, proposed_outputs$marker_counts, na = "")
  readr::write_csv(protein_index, proposed_outputs$protein_identity_crosswalk, na = "")
  readr::write_csv(raw_detection_provenance, proposed_outputs$raw_detection_provenance, na = "")
  readr::write_csv(model_warning_table, proposed_outputs$model_warnings, na = "")
  output_manifest <- proposed_outputs
  manifest_path <- file.path(PATHS$logs, "run_manifest_proposed.yml")
} else {
  write_table_and_source(out, PATHS$tables, PATHS$source_data, "empirical_roi_marker_sets.csv")

  heat <- protein_stats |>
    dplyr::filter(.data$gene_mapping_claim_allowed, !is.na(.data$official_gene_symbol)) |>
    dplyr::select(
      GeneSymbol = .data$official_gene_symbol,
      .data$detection_rate_microglia, .data$detection_rate_neuropil, .data$detection_rate_soma
    ) |>
    dplyr::slice_max(
      order_by = pmax(.data$detection_rate_microglia, .data$detection_rate_neuropil, .data$detection_rate_soma, na.rm = TRUE),
      n = 75
    ) |>
    tidyr::pivot_longer(-"GeneSymbol", names_to = "dataset", values_to = "detection_rate")
  p1 <- ggplot2::ggplot(heat, ggplot2::aes(x = .data$dataset, y = .data$GeneSymbol, fill = .data$detection_rate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#2F6F73", na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Detection") +
    ggplot2::theme_classic(base_size = 7)
  ggplot2::ggsave(file.path(PATHS$figures, "empirical_marker_detection_heatmap.svg"), p1, width = 120, height = 150, units = "mm", device = svglite::svglite)

  lfc <- protein_stats |>
    dplyr::select(
      .data$ProteinGroupID, .data$logFC_microglia_vs_neuropil,
      .data$logFC_microglia_vs_soma, .data$logFC_neuropil_vs_microglia,
      .data$logFC_soma_vs_microglia
    ) |>
    tidyr::pivot_longer(-"ProteinGroupID", names_to = "comparison", values_to = "logFC")
  p2 <- ggplot2::ggplot(lfc, ggplot2::aes(x = .data$comparison, y = .data$logFC)) +
    ggplot2::geom_hline(yintercept = c(-min_abs_logfc, min_abs_logfc), linetype = "dashed", color = "grey55") +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92", color = "grey35") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Animal-level paired limma coefficient") +
    ggplot2::theme_classic(base_size = 8)
  ggplot2::ggsave(file.path(PATHS$figures, "empirical_marker_logFC_summary.svg"), p2, width = 140, height = 90, units = "mm", device = svglite::svglite)

  output_manifest <- list(
    markers = canonical_marker_path,
    figures = PATHS$figures
  )
  manifest_path <- file.path(PATHS$logs, "run_manifest.yml")
}

write_run_manifest(
  manifest_path,
  inputs = inputs,
  outputs = output_manifest,
  parameters = list(
    proposed_only = proposed_only,
    min_detection = min_detection,
    min_abs_logFC = min_abs_logfc,
    FDR = fdr_threshold,
    biological_replicate_unit = "AnimalID",
    model = "limma eBayes: ~ AnimalID + dataset",
    marker_contract_version = empirical_roi_marker_contract_version()
  ),
  notes = "Structured metadata defines anatomy. Hemispheres are averaged within AnimalID/spatial unit; neuropil layers then regions are equally weighted. ExpGroup is not a marker-defining covariate."
)

message(
  if (proposed_only) "Proposed empirical ROI marker validation complete: " else "Empirical ROI marker discovery complete: ",
  selected_marker_path
)
