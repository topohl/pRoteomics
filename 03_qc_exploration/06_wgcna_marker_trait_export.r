#!/usr/bin/env Rscript
# ================================================================
# Script: 03_qc_exploration/06_wgcna_marker_trait_export.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/; optional config/marker_panels/wgcna_reference_marker_sets.csv.
# Produces: results/tables/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Exports QC marker traits from existing WGCNA state; scheduled in QC but requires prior WGCNA outputs when run.
# ================================================================

#
# Export sample-level marker and compartment traits for downstream WGCNA annotation.
# These scores are annotation traits, not purity estimates and not default covariates.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "06_wgcna_marker_trait_export"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_WGCNA_MARKER_TRAIT_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_WGCNA_MARKER_TRAIT_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/06_wgcna_marker_trait_export.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      "Writes WGCNA marker traits by sample to tables and source_data.",
      "Uses reference marker registry and empirical ROI marker sets when present; legacy panels remain fallback.",
      "Marker scores are annotation traits only, not purity estimates and not default covariates."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Marker trait matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_WGCNA_MARKER_TRAIT_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- standardize_wgcna_metadata(expr$meta, DATASET)
marker_sets <- load_wgcna_marker_sets()
marker_source_metadata <- attr(marker_sets, "marker_source_metadata")
primary_panels <- names(marker_sets)

protein_key <- normalize_gene_token(rownames(mat))
z_mat <- t(scale(t(mat)))
z_mat[!is.finite(z_mat)] <- NA_real_
sample_scores <- lapply(names(marker_sets), function(panel) {
  idx <- which(protein_key %in% normalize_gene_token(marker_sets[[panel]]))
  raw_score <- if (length(idx)) colMeans(mat[idx, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(mat))
  z_score <- if (length(idx)) colMeans(z_mat[idx, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(mat))
  data.frame(
    Sample = colnames(mat),
    marker_panel = panel,
    n_markers_detected = length(idx),
    raw_marker_score = raw_score,
    z_marker_score = z_score,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

wide_score <- sample_scores |>
  dplyr::select("Sample", "marker_panel", "raw_marker_score", "z_marker_score") |>
  tidyr::pivot_wider(
    names_from = "marker_panel",
    values_from = c("raw_marker_score", "z_marker_score"),
    names_glue = "{sub('_marker_score$', '', .value)}_{marker_panel}_score"
  )
wide_n <- sample_scores |>
  dplyr::select("Sample", "marker_panel", "n_markers_detected") |>
  tidyr::pivot_wider(names_from = "marker_panel", values_from = "n_markers_detected", names_glue = "n_{marker_panel}_detected")

traits <- meta |>
  dplyr::select(dplyr::any_of(c("Sample", "AnimalID", "Region", "Layer", "RegionLayer", "SpatialUnit", "SpatialLabel", "StressGroup", "ExpGroup", "Sex", "Batch"))) |>
  dplyr::mutate(dataset = DATASET, .before = "Sample") |>
  dplyr::left_join(wide_score, by = "Sample") |>
  dplyr::left_join(wide_n, by = "Sample")

add_pair_scores <- function(df, micro_panel, neuro_panel, prefix = "") {
  raw_micro <- paste0("raw_", micro_panel, "_score")
  raw_neuro <- paste0("raw_", neuro_panel, "_score")
  z_micro <- paste0("z_", micro_panel, "_score")
  z_neuro <- paste0("z_", neuro_panel, "_score")
  out_prefix <- if (nzchar(prefix)) paste0(prefix, "_") else ""
  if (all(c(raw_micro, raw_neuro) %in% names(df))) {
    df[[paste0("raw_", out_prefix, "microglia_minus_neuropil_score")]] <- df[[raw_micro]] - df[[raw_neuro]]
    df[[paste0("raw_", out_prefix, "microglia_to_neuropil_ratio")]] <- ifelse(is.finite(df[[raw_neuro]]) & df[[raw_neuro]] != 0, df[[raw_micro]] / df[[raw_neuro]], NA_real_)
  }
  if (all(c(z_micro, z_neuro) %in% names(df))) {
    df[[paste0("z_", out_prefix, "microglia_minus_neuropil_score")]] <- df[[z_micro]] - df[[z_neuro]]
    df[[paste0("z_", out_prefix, "microglia_to_neuropil_ratio")]] <- ifelse(is.finite(df[[z_neuro]]) & df[[z_neuro]] != 0, df[[z_micro]] / df[[z_neuro]], NA_real_)
  }
  df
}

traits <- add_pair_scores(traits, "canonical_microglia_homeostatic", "canonical_neuronal_synaptic_neuropil")
traits <- add_pair_scores(traits, "microglia", "neuronal_synaptic_neuropil")
traits <- add_pair_scores(traits, "empirical_microglia_roi_enriched", "empirical_neuropil_enriched", "empirical")
if (!"microglia_minus_neuropil_score" %in% names(traits) && "raw_microglia_minus_neuropil_score" %in% names(traits)) traits$microglia_minus_neuropil_score <- traits$raw_microglia_minus_neuropil_score
if (!"microglia_to_neuropil_ratio" %in% names(traits) && "raw_microglia_to_neuropil_ratio" %in% names(traits)) traits$microglia_to_neuropil_ratio <- traits$raw_microglia_to_neuropil_ratio
traits$interpretation_note <- "Marker scores are WGCNA annotation traits only; not purity estimates and not default covariates."

out_file <- "wgcna_marker_traits_by_sample.csv"
write_table_and_source(traits, PATHS$tables, PATHS$source_data, out_file)

summary_tbl <- sample_scores |>
  dplyr::group_by(.data$marker_panel) |>
  dplyr::summarise(
    dataset = DATASET,
    n_samples = dplyr::n_distinct(.data$Sample),
    n_markers_detected = max(.data$n_markers_detected, na.rm = TRUE),
    mean_raw_score = mean(.data$raw_marker_score, na.rm = TRUE),
    median_raw_score = stats::median(.data$raw_marker_score, na.rm = TRUE),
    mean_z_score = mean(.data$z_marker_score, na.rm = TRUE),
    median_z_score = stats::median(.data$z_marker_score, na.rm = TRUE),
    .groups = "drop"
  )
write_table_and_source(summary_tbl, PATHS$tables, PATHS$source_data, "marker_trait_summary_by_dataset.csv")

p_summary <- sample_scores |>
  dplyr::mutate(marker_panel = factor(.data$marker_panel, levels = rev(primary_panels))) |>
  ggplot2::ggplot(ggplot2::aes(x = .data$marker_panel, y = .data$z_marker_score)) +
  ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.25, fill = "grey92", color = "grey35") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12), size = 1, alpha = 0.7, color = "#2F6F73") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean marker z-score") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "marker_trait_summary_by_dataset.svg"), p_summary,
                width = 150, height = 100, units = "mm", device = svglite::svglite)

if (DATASET == "microglia") {
  micro_x <- "z_canonical_neuronal_synaptic_neuropil_score"
  micro_y <- "z_canonical_microglia_homeostatic_score"
  if (all(c("z_empirical_neuropil_enriched_score", "z_empirical_microglia_roi_enriched_score") %in% names(traits))) {
    micro_x <- "z_empirical_neuropil_enriched_score"
    micro_y <- "z_empirical_microglia_roi_enriched_score"
  } else if (!all(c(micro_x, micro_y) %in% names(traits)) && all(c("z_neuronal_synaptic_neuropil_score", "z_microglia_score") %in% names(traits))) {
    micro_x <- "z_neuronal_synaptic_neuropil_score"
    micro_y <- "z_microglia_score"
  }
  p_micro <- traits |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[micro_x]], y = .data[[micro_y]], color = .data$SpatialLabel)) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::labs(x = "Neuronal/synaptic neuropil marker z-score", y = "Microglia marker z-score", color = "Spatial label") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_vs_neuropil_marker_score.svg"), p_micro,
                  width = 120, height = 95, units = "mm", device = svglite::svglite)
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file),
  outputs = list(
    traits = file.path(PATHS$tables, out_file),
    source_traits = file.path(PATHS$source_data, out_file),
    figures = PATHS$figures
  ),
  parameters = list(
    dataset = DATASET,
    marker_panels = primary_panels,
    marker_source_hierarchy = unique(marker_source_metadata$marker_source %||% "unknown"),
    marker_registry_version = attr(marker_sets, "marker_registry_version") %||% NA_character_,
    empirical_marker_set_version = attr(marker_sets, "empirical_marker_set_version") %||% NA_character_
  ),
  notes = "Marker scores are annotation traits for WGCNA reporting/interpretation only; they are not purity estimates and are not default model covariates."
)

message("WGCNA marker trait export complete for dataset: ", DATASET)
