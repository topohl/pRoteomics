#!/usr/bin/env Rscript
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
marker_sets <- wgcna_marker_sets()
primary_panels <- c(
  "microglia", "neuronal_synaptic_neuropil", "nuclear_soma", "astrocyte",
  "oligodendrocyte_myelin", "endothelial_pericyte_vascular",
  "mitochondrial_oxphos", "ribosomal_translation", "rnp_rna_processing"
)
marker_sets <- marker_sets[primary_panels]

protein_key <- normalize_gene_token(rownames(mat))
sample_scores <- lapply(names(marker_sets), function(panel) {
  idx <- which(protein_key %in% normalize_gene_token(marker_sets[[panel]]))
  score <- if (length(idx)) colMeans(mat[idx, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(mat))
  data.frame(
    Sample = colnames(mat),
    marker_panel = panel,
    n_markers_detected = length(idx),
    marker_score = score,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

wide_score <- sample_scores |>
  dplyr::select(.data$Sample, .data$marker_panel, .data$marker_score) |>
  tidyr::pivot_wider(names_from = .data$marker_panel, values_from = .data$marker_score, names_glue = "{marker_panel}_score")
wide_n <- sample_scores |>
  dplyr::select(.data$Sample, .data$marker_panel, .data$n_markers_detected) |>
  tidyr::pivot_wider(names_from = .data$marker_panel, values_from = .data$n_markers_detected, names_glue = "n_{marker_panel}_detected")

traits <- meta |>
  dplyr::select(dplyr::any_of(c("Sample", "Region", "Layer", "RegionLayer", "SpatialUnit", "SpatialLabel", "StressGroup", "ExpGroup", "Sex", "Batch"))) |>
  dplyr::mutate(dataset = DATASET, .before = .data$Sample) |>
  dplyr::left_join(wide_score, by = "Sample") |>
  dplyr::left_join(wide_n, by = "Sample") |>
  dplyr::mutate(
    microglia_minus_neuropil_score = .data$microglia_score - .data$neuronal_synaptic_neuropil_score,
    microglia_to_neuropil_ratio = dplyr::if_else(
      is.finite(.data$neuronal_synaptic_neuropil_score) & .data$neuronal_synaptic_neuropil_score != 0,
      .data$microglia_score / .data$neuronal_synaptic_neuropil_score,
      NA_real_
    ),
    interpretation_note = "Marker scores are WGCNA annotation traits only; not purity estimates and not default covariates."
  )

out_file <- "wgcna_marker_traits_by_sample.csv"
write_table_and_source(traits, PATHS$tables, PATHS$source_data, out_file)

summary_tbl <- sample_scores |>
  dplyr::group_by(.data$marker_panel) |>
  dplyr::summarise(
    dataset = DATASET,
    n_samples = dplyr::n_distinct(.data$Sample),
    n_markers_detected = max(.data$n_markers_detected, na.rm = TRUE),
    mean_score = mean(.data$marker_score, na.rm = TRUE),
    median_score = stats::median(.data$marker_score, na.rm = TRUE),
    .groups = "drop"
  )
write_table_and_source(summary_tbl, PATHS$tables, PATHS$source_data, "marker_trait_summary_by_dataset.csv")

p_summary <- sample_scores |>
  dplyr::mutate(marker_panel = factor(.data$marker_panel, levels = rev(primary_panels))) |>
  ggplot2::ggplot(ggplot2::aes(x = .data$marker_panel, y = .data$marker_score)) +
  ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.25, fill = "grey92", color = "grey35") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12), size = 1, alpha = 0.7, color = "#2F6F73") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean marker abundance") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "marker_trait_summary_by_dataset.svg"), p_summary,
                width = 150, height = 100, units = "mm", device = svglite::svglite)

if (DATASET == "microglia") {
  p_micro <- traits |>
    ggplot2::ggplot(ggplot2::aes(x = .data$neuronal_synaptic_neuropil_score, y = .data$microglia_score, color = .data$SpatialLabel)) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::labs(x = "Neuronal/synaptic neuropil marker score", y = "Microglia marker score", color = "Spatial label") +
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
  parameters = list(dataset = DATASET, marker_panels = primary_panels),
  notes = "Marker scores are annotation traits for WGCNA reporting/interpretation only; they are not purity estimates and are not default model covariates."
)

message("WGCNA marker trait export complete for dataset: ", DATASET)
