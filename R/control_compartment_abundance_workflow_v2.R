# Authoritative v2 workflow for control-only compartment-marker abundance.
# Sourced by 03_qc_exploration/04e_control_compartment_abundance_publication_figures.r.

source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "compartment_abundance_utils.R"))
source(repo_path("R", "plotting_nature.R"))
source(repo_path("R", "joint_compartment_qc_plotting.R"))
source(repo_path("R", "control_compartment_abundance_rendering.R"))

global_arg <- tolower(script_arg_value("--dataset", "global"))
if (!global_arg %in% c("all", "global")) {
  stop("This is a global script; use --dataset all or global.", call. = FALSE)
}
dry_run <- is_dry_run()
render_only <- "--render-only" %in% commandArgs(trailingOnly = TRUE)
compact_candidate_only <- "--compact-recognizable-candidate" %in%
  commandArgs(trailingOnly = TRUE)
compact_difference_only <- "--compact-enrichment-differences" %in%
  commandArgs(trailingOnly = TRUE)

as_integer_contract <- function(env, default, lower, upper) {
  value <- suppressWarnings(as.integer(Sys.getenv(env, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower || value > upper) {
    stop(env, " must be an integer from ", lower, " to ", upper, ".", call. = FALSE)
  }
  value
}

min_regions <- as_integer_contract(
  "PROTEOMICS_CONTROL_ABUNDANCE_MIN_REGIONS", 3L, 1L, 4L
)
min_con_animals <- as_integer_contract(
  "PROTEOMICS_CONTROL_ABUNDANCE_MIN_CON_ANIMALS", 2L, 1L, 3L
)
strict_con_animals <- as_integer_contract(
  "PROTEOMICS_CONTROL_ABUNDANCE_STRICT_CON_ANIMALS", 3L, 1L, 3L
)
min_score_proteins <- as_integer_contract(
  "PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_PROTEINS", 3L, 1L, 10000L
)
min_score_fraction <- suppressWarnings(as.numeric(Sys.getenv(
  "PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_FRACTION", unset = "0.5"
)))
if (!is.finite(min_score_fraction) || min_score_fraction < 0 || min_score_fraction > 1) {
  stop(
    "PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_FRACTION must be between 0 and 1.",
    call. = FALSE
  )
}

processed_root <- Sys.getenv(
  "PROTEOMICS_JOINT_QC_PROCESSED_DIR",
  unset = path_processed("01_preprocessing", "joint_compartment_qc", "global")
)
inputs <- c(
  joint_bundle = file.path(processed_root, "joint_compartment_qc_matrices.rds"),
  marker_registry = repo_path(
    "config", "marker_panels", "wgcna_reference_marker_sets.csv"
  )
)
substep <- "04e_control_compartment_abundance_publication_figures"
roots <- list(
  source_data = path_results(
    "source_data", "03_qc_exploration", substep, "global"
  ),
  figures = path_results(
    "figures", "03_qc_exploration", substep, "global"
  ),
  reports = path_results(
    "reports", "03_qc_exploration", substep, "global"
  ),
  logs = path_results(
    "logs", "03_qc_exploration", substep, "global"
  )
)
render_inputs <- c(
  dot_source = file.path(
    roots$source_data, "v2_17_marker_detection_dot_heatmap_source.csv"
  ),
  rank_source = file.path(
    roots$source_data, "v2_12_rank_abundance_data.csv"
  ),
  display_selection = file.path(
    roots$source_data, "v2_11_display_selection_provenance.csv"
  )
)
compact_candidate_inputs <- c(
  animal_abundance = file.path(
    roots$source_data, "v2_03_animal_level_abundance.csv"
  ),
  marker_mapping = file.path(
    roots$source_data, "v2_04_all_marker_canonical_mapping_provenance.csv"
  ),
  marker_detection = file.path(
    roots$source_data, "v2_05_marker_by_dataset_detection_abundance.csv"
  ),
  marker_direction = file.path(
    roots$source_data, "v2_06_protein_intended_direction.csv"
  ),
  eligible_alternatives = file.path(
    roots$source_data, "v2_11a_display_selection_eligible_alternatives.csv"
  ),
  rank_source = file.path(
    roots$source_data, "v2_12_rank_abundance_data.csv"
  )
)
compact_difference_inputs <- c(
  compact_source = file.path(
    roots$source_data, "figure2d_compact_recognizable_markers_v2_source_data.csv"
  )
)

if (dry_run) {
  required <- if (compact_difference_only) {
    compact_difference_inputs
  } else if (compact_candidate_only) {
    compact_candidate_inputs
  } else if (render_only) {
    render_inputs
  } else {
    inputs
  }
  for (nm in names(required)) {
    cat(
      "[DRY-RUN] ", nm, ": ",
      normalizePath(required[[nm]], winslash = "/", mustWork = FALSE),
      " [", if (file.exists(required[[nm]])) "PASS" else "FAIL", "]\n",
      sep = ""
    )
  }
  cat("[DRY-RUN] authoritative workflow: 04e v2\n")
  cat("[DRY-RUN] primary group: CON from ExpGroup\n")
  cat(
    "[DRY-RUN] marker abundance: raw-positive log2 + joint shared-core offsets; ",
    "observed non-imputed values only\n",
    sep = ""
  )
  cat(
    "[DRY-RUN] hierarchy: technical replicate -> layer -> region within hemisphere -> ",
    "equal hemispheres -> animal; >=", min_regions, " regions per valid hemisphere\n",
    sep = ""
  )
  cat(
    "[DRY-RUN] marker eligibility: external registry + canonical non-contaminant + ",
    "nonconflicting + intended detection >=", min_con_animals,
    "/3; shared core sensitivity only\n",
    sep = ""
  )
  cat(
    "[DRY-RUN] v2 outputs do not overwrite existing v1 analytical files: ",
    normalizePath(roots$source_data, winslash = "/", mustWork = FALSE), "\n",
    sep = ""
  )
  quit(status = if (all(file.exists(required))) 0 else 1, save = "no")
}

missing_inputs <- inputs[!file.exists(inputs)]
if (!render_only && !compact_candidate_only && !compact_difference_only &&
    length(missing_inputs)) {
  stop(
    "Missing required input(s): ",
    paste(names(missing_inputs), collapse = ", "),
    call. = FALSE
  )
}
missing_render_inputs <- render_inputs[!file.exists(render_inputs)]
if (render_only && length(missing_render_inputs)) {
  stop(
    "Missing completed v2 rendering input(s): ",
    paste(names(missing_render_inputs), collapse = ", "),
    call. = FALSE
  )
}
missing_compact_candidate_inputs <- compact_candidate_inputs[
  !file.exists(compact_candidate_inputs)
]
if (compact_candidate_only && length(missing_compact_candidate_inputs)) {
  stop(
    "Missing completed v2 compact-candidate input(s): ",
    paste(names(missing_compact_candidate_inputs), collapse = ", "),
    call. = FALSE
  )
}
missing_compact_difference_inputs <- compact_difference_inputs[
  !file.exists(compact_difference_inputs)
]
if (compact_difference_only && length(missing_compact_difference_inputs)) {
  stop(
    "Missing completed compact Figure 2d source: ",
    paste(names(missing_compact_difference_inputs), collapse = ", "),
    call. = FALSE
  )
}
allow_overwrite <- ca_as_logical(Sys.getenv(
  "PROTEOMICS_CONTROL_ABUNDANCE_V2_ALLOW_OVERWRITE",
  unset = "false"
))
if (!render_only && !compact_candidate_only && !compact_difference_only &&
    !allow_overwrite) {
  existing_v2_outputs <- c(
    list.files(
      roots$source_data, pattern = "^v2_", full.names = TRUE,
      recursive = FALSE
    ),
    list.files(
      roots$figures, pattern = "_v2_", full.names = TRUE,
      recursive = FALSE
    ),
    list.files(
      roots$reports, pattern = "^v2_", full.names = TRUE,
      recursive = FALSE
    ),
    list.files(
      roots$logs, pattern = "^v2_", full.names = TRUE,
      recursive = FALSE
    )
  )
  if (length(existing_v2_outputs)) {
    stop(
      "Refusing to overwrite existing v2 output(s): ",
      paste(
        normalizePath(
          existing_v2_outputs, winslash = "/", mustWork = FALSE
        ),
        collapse = ", "
      ),
      ". Set PROTEOMICS_CONTROL_ABUNDANCE_V2_ALLOW_OVERWRITE=true only for ",
      "an intentional analytical regeneration.",
      call. = FALSE
    )
  }
}
required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "scales", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

rendering_source_files <- c(
  "v2_18_rendering_marker_display.csv",
  "v2_18a_rendering_dot_display.csv",
  "v2_18b_rendering_rank_display.csv",
  "v2_18c_terminal_rank_drop_audit.csv"
)
rendering_figure_files <- c(
  "control_compartment_marker_detection_dot_heatmap_v2_120mm.svg",
  "control_compartment_marker_detection_dot_heatmap_v2_120mm.pdf",
  "control_rank_abundance_v2_extended_data_183mm.svg",
  "control_rank_abundance_v2_extended_data_183mm.pdf",
  "control_rank_abundance_v2_extended_data_183mm.png"
)

ca_require_candidate_outputs_writable <- function(
    paths, allow_overwrite, candidate_label, overwrite_env) {
  existing <- paths[file.exists(paths)]
  if (length(existing) && !allow_overwrite) {
    stop(
      "Refusing to overwrite existing ", candidate_label, " output(s): ",
      paste(normalizePath(existing, winslash = "/", mustWork = FALSE), collapse = ", "),
      ". Set ", overwrite_env,
      "=true only for intentional candidate rerendering.",
      call. = FALSE
    )
  }
  invisible(paths)
}

write_control_rendering_outputs_v2 <- function(rendering, roots, inputs_used,
                                               allow_overwrite = FALSE) {
  paths <- c(
    file.path(roots$source_data, rendering_source_files),
    file.path(roots$figures, rendering_figure_files),
    file.path(roots$reports, "v2_rendering_figure_legends.md"),
    file.path(roots$logs, "v2_rendering_manifest.yml")
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) && !allow_overwrite) {
    stop(
      "Refusing to overwrite existing v2 rendering output(s): ",
      paste(
        normalizePath(existing, winslash = "/", mustWork = FALSE),
        collapse = ", "
      ),
      ". Set PROTEOMICS_CONTROL_ABUNDANCE_V2_RENDER_ALLOW_OVERWRITE=true ",
      "only for an intentional v2 render regeneration.",
      call. = FALSE
    )
  }
  invisible(lapply(roots, dir_create))
  qc_write_csv(
    rendering$config,
    file.path(roots$source_data, rendering_source_files[[1]])
  )
  qc_write_csv(
    rendering$dot,
    file.path(roots$source_data, rendering_source_files[[2]])
  )
  qc_write_csv(
    rendering$rank_representatives,
    file.path(roots$source_data, rendering_source_files[[3]])
  )
  qc_write_csv(
    rendering$terminal_drop,
    file.path(roots$source_data, rendering_source_files[[4]])
  )

  dot_plot <- ca_build_control_marker_dot_heatmap_v2(rendering)
  rank_plot <- ca_build_control_rank_abundance_plot_v2(rendering)
  marker_tile_height_mm <- max(82, 22 + 3.7 * length(unique(rendering$dot$marker_label)))
  marker_tile_width_mm <- max(108, 54 + 18 * length(unique(rendering$dot$dataset_label)))
  ca_save_vector_figure(
    dot_plot, roots$figures,
    "control_compartment_marker_detection_dot_heatmap_v2_120mm",
    marker_tile_width_mm, marker_tile_height_mm
  )
  ca_save_vector_figure(
    rank_plot, roots$figures,
    "control_rank_abundance_v2_extended_data_183mm",
    183, 52
  )
  ggplot2::ggsave(
    file.path(
      roots$figures,
      "control_rank_abundance_v2_extended_data_183mm.png"
    ),
    rank_plot,
    width = 183, height = 52, units = "mm", dpi = 300,
    device = "png", limitsize = FALSE, bg = "white"
  )
  legends <- ca_control_rendering_legends_v2()
  writeLines(
    c(
      "# Main marker panel", "", legends$main, "",
      "# Extended Data rank-abundance panel", "", legends$extended
    ),
    file.path(roots$reports, "v2_rendering_figure_legends.md")
  )
  write_run_manifest(
    file.path(roots$logs, "v2_rendering_manifest.yml"),
    inputs = as.list(inputs_used),
    outputs = list(
      source_files = file.path(roots$source_data, rendering_source_files),
      figure_files = file.path(roots$figures, rendering_figure_files),
      legends = file.path(roots$reports, "v2_rendering_figure_legends.md")
    ),
    parameters = list(
      execution_mode = "render_only_from_completed_v2_source_tables",
      analytical_recalculation = FALSE,
      fill = "median_within_protein_centered_log2",
      mark_geometry = "uniform_tile; outline_colour retained as a non-area cue",
      unreliable_detection_symbol = "multiplication_sign",
      display_cap_log2 = 3
    ),
    notes = paste(
      "Rendering-only regeneration from completed v2 source tables.",
      "No inclusion, normalization, detection, aggregation, marker selection or",
      "analytical values were recalculated."
    )
  )
  invisible(list(dot = dot_plot, rank = rank_plot, outputs = paths))
}

write_compact_recognizable_candidate_v2 <- function(
    rendering, roots, inputs_used, allow_overwrite = FALSE) {
  source_paths <- c(
    audit = file.path(
      roots$source_data, "figure2d_requested_marker_eligibility_audit_v2.csv"
    ),
    source_data = file.path(
      roots$source_data, "figure2d_compact_recognizable_markers_v2_source_data.csv"
    ),
    provenance = file.path(
      roots$source_data, "figure2d_compact_recognizable_markers_v2_provenance.csv"
    )
  )
  figure_path <- file.path(
    roots$figures, "figure2d_compact_recognizable_markers_v2.svg"
  )
  manifest_path <- file.path(
    roots$logs, "figure2d_compact_recognizable_markers_v2_manifest.yml"
  )
  paths <- c(source_paths, figure = figure_path, manifest = manifest_path)
  ca_require_candidate_outputs_writable(
    paths, allow_overwrite, "compact Figure 2d candidate",
    "PROTEOMICS_CONTROL_ABUNDANCE_COMPACT_CANDIDATE_ALLOW_OVERWRITE"
  )
  invisible(lapply(roots, dir_create))
  qc_write_csv(rendering$audit, source_paths[["audit"]])
  qc_write_csv(rendering$dot, source_paths[["source_data"]])
  qc_write_csv(rendering$provenance, source_paths[["provenance"]])
  plot <- ca_build_compact_recognizable_markers_v2(rendering)
  ggplot2::ggsave(
    figure_path, plot, width = 132, height = 76, units = "mm",
    device = function(...) svglite::svglite(
      ..., pointsize = 7, fix_text_size = FALSE
    ),
    limitsize = FALSE, bg = "white"
  )
  write_run_manifest(
    manifest_path,
    inputs = as.list(inputs_used),
    outputs = list(
      marker_audit = source_paths[["audit"]],
      source_data = source_paths[["source_data"]],
      provenance = source_paths[["provenance"]],
      figure = figure_path
    ),
    parameters = list(
      execution_mode = "presentation_only_from_completed_authoritative_v2_tables",
      analytical_recalculation = FALSE,
      marker_eligibility_recalculation = FALSE,
      observed_cross_compartment_direction_used_for_selection = FALSE,
      observed_cross_compartment_effect_magnitude_used_for_selection = FALSE,
      marker_set = rendering$config$marker_gene,
      metric = "median_within_protein_centered_log2_across_valid_CON_animals",
      display_cap_log2 = 3,
      columns = c(
        "Neuron soma", "Neuron neuropil", "Microglia/PVM-enriched ROI"
      )
    ),
    notes = paste(
      "Compact manuscript candidate derived from completed authoritative v2 tables.",
      "Animal aggregation, normalization, detection, eligibility and direction values",
      "were not recomputed; the existing ca_center_animal_abundance presentation",
      "transform was applied to completed animal-level abundance rows.",
      "The recognizable marker set was frozen from external/configuration and",
      "intended-compartment eligibility only. Observed cross-compartment direction",
      "and effect magnitude were retained for post-selection audit, not selection."
    )
  )
  invisible(list(plot = plot, outputs = paths))
}

write_marker_enrichment_difference_candidate_v2 <- function(
    rendering, roots, inputs_used, allow_overwrite = FALSE) {
  source_path <- file.path(
    roots$source_data,
    "figure2d_compact_marker_enrichment_differences_v2_source_data.csv"
  )
  provenance_path <- file.path(
    roots$source_data,
    "figure2d_compact_marker_enrichment_differences_v2_provenance.csv"
  )
  figure_path <- file.path(
    roots$figures, "figure2d_compact_marker_enrichment_differences_v2.svg"
  )
  manifest_path <- file.path(
    roots$logs, "figure2d_compact_marker_enrichment_differences_v2_manifest.yml"
  )
  paths <- c(
    source_data = source_path, provenance = provenance_path,
    figure = figure_path, manifest = manifest_path
  )
  ca_require_candidate_outputs_writable(
    paths, allow_overwrite, "Figure 2d enrichment-difference candidate",
    "PROTEOMICS_CONTROL_ABUNDANCE_DIFFERENCE_CANDIDATE_ALLOW_OVERWRITE"
  )
  invisible(lapply(roots, dir_create))
  source_for_csv <- rendering$source
  exact_numeric_columns <- c(
    "intended_median_centered_log2", "comparator_median_centered_log2",
    "intended_minus_comparator_log2"
  )
  source_for_csv[exact_numeric_columns] <- lapply(
    source_for_csv[exact_numeric_columns],
    function(x) ifelse(is.finite(x), sprintf("%.17g", x), NA_character_)
  )
  qc_write_csv(source_for_csv, source_path)
  qc_write_csv(rendering$provenance, provenance_path)
  plot <- ca_build_marker_enrichment_differences_v2(rendering)
  ggplot2::ggsave(
    figure_path, plot, width = 132, height = 92, units = "mm",
    device = function(...) svglite::svglite(
      ..., pointsize = 7, fix_text_size = FALSE
    ),
    limitsize = FALSE, bg = "white"
  )
  write_run_manifest(
    manifest_path,
    inputs = as.list(inputs_used),
    outputs = list(
      source_data = source_path, provenance = provenance_path,
      figure = figure_path
    ),
    parameters = list(
      execution_mode = "presentation_only_from_completed_compact_candidate_source",
      analytical_recalculation = FALSE,
      marker_eligibility_recalculation = FALSE,
      marker_selection_recalculation = FALSE,
      observed_cross_compartment_direction_used_for_selection = FALSE,
      observed_cross_compartment_effect_magnitude_used_for_selection = FALSE,
      metric = paste(
        "intended preparation median centered log2 minus comparator preparation",
        "median centered log2"
      ),
      comparisons_per_marker = 2,
      rescaling = "none",
      inferential_statistics = "none"
    ),
    notes = paste(
      "Numerical values are exact arithmetic subtractions of two quantitatively",
      "reliable completed figure2d compact-source values. Comparisons lacking a",
      "reliable intended or comparator value remain explicitly not quantitatively",
      "estimable; their missing differences are not replaced.",
      "No model, p-value, aggregation, normalization, marker selection or",
      "eligibility calculation was run."
    )
  )
  invisible(list(plot = plot, outputs = paths))
}

if (compact_difference_only) {
  input_path <- compact_difference_inputs[["compact_source"]]
  rendering <- ca_prepare_marker_enrichment_differences_v2(
    compact_source = qc_read_table(input_path),
    source_path = normalizePath(input_path, winslash = "/", mustWork = TRUE),
    source_sha256 = file_hash_sha256(input_path)
  )
  allow_difference_overwrite <- ca_as_logical(Sys.getenv(
    "PROTEOMICS_CONTROL_ABUNDANCE_DIFFERENCE_CANDIDATE_ALLOW_OVERWRITE",
    unset = "false"
  ))
  write_marker_enrichment_difference_candidate_v2(
    rendering, roots, compact_difference_inputs, allow_difference_overwrite
  )
  message("04e compact marker-enrichment difference rendering complete.")
  quit(status = 0, save = "no")
}

if (compact_candidate_only) {
  rendering <- ca_prepare_compact_recognizable_markers_v2(
    animal_abundance = qc_read_table(compact_candidate_inputs[["animal_abundance"]]),
    marker_mapping = qc_read_table(compact_candidate_inputs[["marker_mapping"]]),
    marker_detection = qc_read_table(compact_candidate_inputs[["marker_detection"]]),
    marker_direction = qc_read_table(compact_candidate_inputs[["marker_direction"]]),
    eligible_alternatives = qc_read_table(
      compact_candidate_inputs[["eligible_alternatives"]]
    ),
    rank_table = qc_read_table(compact_candidate_inputs[["rank_source"]])
  )
  allow_compact_candidate_overwrite <- ca_as_logical(Sys.getenv(
    "PROTEOMICS_CONTROL_ABUNDANCE_COMPACT_CANDIDATE_ALLOW_OVERWRITE",
    unset = "false"
  ))
  write_compact_recognizable_candidate_v2(
    rendering, roots, compact_candidate_inputs,
    allow_compact_candidate_overwrite
  )
  message("04e compact recognizable-marker candidate rendering complete.")
  quit(status = 0, save = "no")
}

if (render_only) {
  rendering <- ca_prepare_control_rendering_sources_v2(
    dot_source = qc_read_table(render_inputs[["dot_source"]]),
    rank_table = qc_read_table(render_inputs[["rank_source"]]),
    display_selection = qc_read_table(render_inputs[["display_selection"]])
  )
  allow_render_overwrite <- ca_as_logical(Sys.getenv(
    "PROTEOMICS_CONTROL_ABUNDANCE_V2_RENDER_ALLOW_OVERWRITE",
    unset = "false"
  ))
  write_control_rendering_outputs_v2(
    rendering, roots, render_inputs, allow_render_overwrite
  )
  message("04e v2 render-only regeneration complete.")
  quit(status = 0, save = "no")
}

bundle <- readRDS(inputs[["joint_bundle"]])
reconstructed <- ca_validate_bundle_reconstruction(bundle)
feature_table <- reconstructed$feature_table
eligible_ids <- as.character(
  feature_table$ProteinGroupID[feature_table$joint_qc_eligible %in% TRUE]
)
shared_ids <- intersect(reconstructed$primary_feature_ids, eligible_ids)
broad_ids <- intersect(reconstructed$broad_feature_ids, eligible_ids)
analysis_mat <- reconstructed$normalized_nonimputed[
  eligible_ids, , drop = FALSE
]
meta <- reconstructed$metadata
if (!identical(
  colnames(analysis_mat),
  as.character(meta[[reconstructed$columns$sample]])
)) {
  stop("Reconstructed non-imputed values are not aligned to metadata.", call. = FALSE)
}

con <- ca_filter_expression_group(analysis_mat, meta, "CON")
if (!identical(con$resolved_group_column, "ExpGroup")) {
  stop(
    "The primary analysis must resolve experimental group from ExpGroup.",
    call. = FALSE
  )
}
animals <- sort(unique(as.character(con$meta[[con$columns$animal]])))
datasets <- c("neuron_soma", "neuron_neuropil", "microglia")
if (length(animals) != 3L) {
  stop("The v2 descriptive stability contract requires exactly three CON animals.", call. = FALSE)
}
if (!setequal(unique(as.character(con$meta[[con$columns$dataset]])), datasets)) {
  stop("CON metadata does not contain exactly the three required datasets.", call. = FALSE)
}

median_agg <- ca_aggregate_abundance_v2(
  con$mat, con$meta, "median", min_regions, require_both_hemispheres = FALSE
)
mean_agg <- ca_aggregate_abundance_v2(
  con$mat, con$meta, "mean", min_regions, require_both_hemispheres = FALSE
)
both_hemisphere_agg <- ca_aggregate_abundance_v2(
  con$mat, con$meta, "median", min_regions, require_both_hemispheres = TRUE
)
detection <- ca_detection_audit_v2(
  median_agg, min_con_animals, strict_con_animals
)
both_hemisphere_detection <- ca_detection_audit_v2(
  both_hemisphere_agg, min_con_animals, strict_con_animals
)
detection <- detection |>
  dplyr::mutate(
    in_joint_shared_core = .data$ProteinGroupID %in% shared_ids,
    in_broad_union = .data$ProteinGroupID %in% broad_ids,
    imputed_sample_count_in_values_used = 0L,
    imputed_fraction_in_values_used = 0,
    observed_detection_source =
      "raw_positive_joint_bundle_before_imputation",
    normalization_source =
      "joint_shared_core_sample_median_offsets_applied_to_raw_positive_log2"
  )

mean_estimates <- ca_dataset_estimates(
  mean_agg, "primary_region_balanced_mean_log2"
)
median_estimates <- ca_dataset_estimates(
  median_agg, "median_hierarchy_log2"
)
primary_rank <- ca_rank_table(
  mean_estimates, detection, "primary_eligible",
  "primary_region_balanced_mean_log2",
  "CON_primary_2_of_3", "CON_region_balanced_mean"
)
strict_rank <- ca_rank_table(
  mean_estimates, detection, "strict_eligible",
  "primary_region_balanced_mean_log2",
  "CON_strict_3_of_3", "CON_region_balanced_mean_strict"
)
median_rank <- ca_rank_table(
  median_estimates, detection, "primary_eligible",
  "median_hierarchy_log2",
  "CON_primary_2_of_3_median_hierarchy",
  "CON_median_hierarchy_sensitivity"
)
feature_annotation <- feature_table[, intersect(c(
  "ProteinGroupID", "FeatureDisplayLabel", "original_identifier",
  "member_accessions", "member_gene_symbols", "official_gene_symbol",
  "protein_group_ambiguity_class", "gene_level_claim_allowed",
  "joint_qc_eligible", "joint_qc_exclusion_reason"
), names(feature_table)), drop = FALSE]
rank_data <- dplyr::bind_rows(primary_rank, strict_rank, median_rank) |>
  dplyr::left_join(
    feature_annotation,
    by = "ProteinGroupID",
    relationship = "many-to-one"
  )

registry <- utils::read.csv(
  inputs[["marker_registry"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
registry <- ca_standardize_marker_registry_ranks(registry)
marker_matches <- ca_match_fidelity_markers(
  registry, feature_table, shared_ids
)
marker_eligibility <- ca_apply_marker_eligibility_v2(
  marker_matches, detection, shared_ids, broad_ids
)
marker_map_all <- ca_prepare_marker_map_v2(marker_eligibility)
marker_map_primary <- ca_prepare_marker_map_v2(
  marker_eligibility, "intended_primary_eligible"
)
marker_map_strict <- ca_prepare_marker_map_v2(
  marker_eligibility, "intended_strict_eligible"
)

animal_median <- ca_animal_abundance_table_v2(
  median_agg, "animal_log2_abundance"
)
animal_mean <- ca_animal_abundance_table_v2(
  mean_agg, "animal_log2_abundance"
)
animal_both <- ca_animal_abundance_table_v2(
  both_hemisphere_agg, "animal_log2_abundance"
)
hemisphere_abundance <- ca_hemisphere_abundance_table(
  median_agg, "hemisphere_log2_abundance"
)
centered <- ca_center_animal_abundance(animal_median)
centered_both <- ca_center_animal_abundance(animal_both)
marker_scores <- ca_score_markers_v2(
  centered, marker_map_primary, min_score_proteins, min_score_fraction
)
both_hemisphere_scores <- ca_score_markers_v2(
  centered_both, marker_map_primary, min_score_proteins, min_score_fraction
)

legacy_robust <- ca_robust_standardize(
  animal_mean[animal_mean$ProteinGroupID %in% shared_ids, , drop = FALSE]
)
legacy_marker_rows <- marker_eligibility |>
  dplyr::transmute(
    ProteinGroupID,
    marker_class = .data$fidelity_marker_class,
    marker_gene,
    marker_panel,
    fidelity_subpanel,
    primary_score_eligible = .data$shared_core_sensitivity_eligible
  )
legacy_scores <- ca_score_markers(
  legacy_robust, legacy_marker_rows,
  min_score_proteins, min_score_fraction
)

marker_dataset_audit <- ca_marker_dataset_audit_v2(
  marker_map_all, detection
)
protein_direction <- ca_protein_direction_v2(marker_dataset_audit)
display <- ca_select_display_markers_v2(marker_map_primary, 6L)

centered_dataset <- centered |>
  dplyr::group_by(.data$dataset, .data$ProteinGroupID) |>
  dplyr::summarise(
    median_centered_log2 = stats::median(.data$centered_log2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    median_centered_log2 = ifelse(
      is.finite(.data$median_centered_log2),
      .data$median_centered_log2,
      NA_real_
    )
  )
dot_source <- marker_dataset_audit |>
  dplyr::filter(.data$ProteinGroupID %in% display$selected$ProteinGroupID) |>
  dplyr::select(
    "dataset", "ProteinGroupID", "valid_animal_count",
    "valid_animal_fraction", "reliability_status", "intended_compartment"
  ) |>
  dplyr::left_join(
    centered_dataset,
    by = c("dataset", "ProteinGroupID"),
    relationship = "one-to-one"
  )
rendering <- ca_prepare_control_rendering_sources_v2(
  dot_source, rank_data, display$selected
)

hemisphere_markers <- hemisphere_abundance |>
  dplyr::filter(.data$ProteinGroupID %in% marker_map_primary$ProteinGroupID) |>
  dplyr::left_join(
    marker_map_primary |>
      dplyr::select(
        "ProteinGroupID", "marker_class", "marker_gene",
        "display_biological_subpanel"
      ),
    by = "ProteinGroupID",
    relationship = "many-to-many"
  )
left <- hemisphere_markers |>
  dplyr::filter(.data$hemisphere == "Left") |>
  dplyr::select(
    "dataset", "AnimalID", "marker_class", "ProteinGroupID",
    left_log2 = "hemisphere_log2_abundance",
    left_regions = "contributing_region_count"
  )
right <- hemisphere_markers |>
  dplyr::filter(.data$hemisphere == "Right") |>
  dplyr::select(
    "dataset", "AnimalID", "marker_class", "ProteinGroupID",
    right_log2 = "hemisphere_log2_abundance",
    right_regions = "contributing_region_count"
  )
left_right_pairs <- dplyr::full_join(
  left, right,
  by = c("dataset", "AnimalID", "marker_class", "ProteinGroupID"),
  relationship = "one-to-one"
) |>
  dplyr::mutate(
    both_hemispheres_valid =
      is.finite(.data$left_log2) & is.finite(.data$right_log2),
    right_minus_left_log2 = ifelse(
      .data$both_hemispheres_valid,
      .data$right_log2 - .data$left_log2,
      NA_real_
    ),
    absolute_left_right_difference_log2 = abs(.data$right_minus_left_log2)
  )
left_right_summary <- left_right_pairs |>
  dplyr::group_by(.data$dataset, .data$marker_class) |>
  dplyr::summarise(
    n_valid_pairs = sum(.data$both_hemispheres_valid),
    median_right_minus_left_log2 =
      stats::median(.data$right_minus_left_log2, na.rm = TRUE),
    median_absolute_left_right_difference_log2 =
      stats::median(.data$absolute_left_right_difference_log2, na.rm = TRUE),
    spearman_left_right = ifelse(
      sum(.data$both_hemispheres_valid) >= 3L,
      stats::cor(
        .data$left_log2[.data$both_hemispheres_valid],
        .data$right_log2[.data$both_hemispheres_valid],
        method = "spearman"
      ),
      NA_real_
    ),
    .groups = "drop"
  )

loao_cases <- ca_leave_one_animal_out_cases(animals)
loao_score_parts <- list()
loao_direction_parts <- list()
for (i in seq_len(nrow(loao_cases))) {
  omitted <- loao_cases$omitted_AnimalID[[i]]
  reduced <- animal_median[animal_median$AnimalID != omitted, , drop = FALSE]
  reduced_centered <- ca_center_animal_abundance(reduced)
  reduced_scores <- ca_score_markers_v2(
    reduced_centered, marker_map_primary,
    min_score_proteins, min_score_fraction
  )$scores |>
    dplyr::mutate(
      omission_id = loao_cases$omission_id[[i]],
      omitted_AnimalID = omitted
    )
  loao_score_parts[[i]] <- reduced_scores
  estimates <- reduced |>
    dplyr::group_by(.data$dataset, .data$ProteinGroupID) |>
    dplyr::summarise(
      omission_median_log2 = stats::median(
        .data$animal_log2_abundance, na.rm = TRUE
      ),
      .groups = "drop"
    )
  direction_rows <- lapply(seq_len(nrow(marker_map_primary)), function(j) {
    marker <- marker_map_primary[j, , drop = FALSE]
    x <- estimates[estimates$ProteinGroupID == marker$ProteinGroupID, , drop = FALSE]
    intended <- x$omission_median_log2[x$dataset == marker$intended_dataset]
    off <- x$omission_median_log2[x$dataset != marker$intended_dataset]
    off <- off[is.finite(off)]
    margin <- if (length(intended) && is.finite(intended) && length(off)) {
      intended - max(off)
    } else {
      NA_real_
    }
    data.frame(
      omission_id = loao_cases$omission_id[[i]],
      omitted_AnimalID = omitted,
      marker_class = marker$marker_class,
      ProteinGroupID = marker$ProteinGroupID,
      intended_margin_log2 = margin,
      intended_is_highest = if (is.finite(margin)) margin > 0 else NA,
      stringsAsFactors = FALSE
    )
  })
  loao_direction_parts[[i]] <- dplyr::bind_rows(direction_rows)
}
loao_scores <- dplyr::bind_rows(loao_score_parts)
loao_direction <- dplyr::bind_rows(loao_direction_parts)
loao_class_stability <- loao_scores |>
  dplyr::group_by(
    .data$omission_id, .data$omitted_AnimalID,
    .data$marker_class, .data$dataset
  ) |>
  dplyr::summarise(
    class_score = stats::median(
      .data$subpanel_balanced_median_centered_log2, na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::group_by(
    .data$omission_id, .data$omitted_AnimalID, .data$marker_class
  ) |>
  dplyr::mutate(
    intended_dataset = ca_intended_dataset(.data$marker_class),
    intended_class_score_highest =
      .data$class_score[.data$dataset == .data$intended_dataset] ==
      max(.data$class_score, na.rm = TRUE)
  ) |>
  dplyr::ungroup()

draws <- ca_exact_paired_bootstrap_draws(animals)
bootstrap_parts <- lapply(seq_len(nrow(draws)), function(i) {
  drawn <- as.character(unlist(
    draws[i, c("draw_1", "draw_2", "draw_3")],
    use.names = FALSE
  ))
  per_animal_scores <- marker_scores$scores
  index <- match(drawn, per_animal_scores$AnimalID)
  pieces <- lapply(seq_along(drawn), function(k) {
    x <- per_animal_scores[per_animal_scores$AnimalID == drawn[[k]], , drop = FALSE]
    x$draw_position <- k
    x
  })
  dplyr::bind_rows(pieces) |>
    dplyr::group_by(.data$dataset, .data$marker_class) |>
    dplyr::summarise(
      bootstrap_subpanel_balanced_median_centered_log2 =
        stats::median(
          .data$subpanel_balanced_median_centered_log2,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      bootstrap_draw_id = draws$bootstrap_draw_id[[i]],
      draw_1 = drawn[[1]], draw_2 = drawn[[2]], draw_3 = drawn[[3]],
      resampling_interpretation =
        "exact_27_ordered_descriptive_draws_not_inferential_CI"
    )
})
bootstrap_class_scores <- dplyr::bind_rows(bootstrap_parts)

requested_counts <- marker_eligibility |>
  dplyr::group_by(.data$fidelity_marker_class) |>
  dplyr::summarise(
    requested_external_genes = dplyr::n_distinct(.data$marker_gene),
    canonical_mapped_genes = dplyr::n_distinct(
      .data$marker_gene[
        .data$matched %in% TRUE &
          .data$canonical_marker_eligible %in% TRUE
      ]
    ),
    .groups = "drop"
  ) |>
  dplyr::rename(marker_class = "fidelity_marker_class")
eligible_counts <- marker_eligibility |>
  dplyr::group_by(.data$fidelity_marker_class) |>
  dplyr::summarise(
    intended_primary_eligible_proteins = dplyr::n_distinct(
      .data$ProteinGroupID[.data$intended_primary_eligible %in% TRUE]
    ),
    intended_strict_eligible_proteins = dplyr::n_distinct(
      .data$ProteinGroupID[.data$intended_strict_eligible %in% TRUE]
    ),
    shared_core_sensitivity_proteins = dplyr::n_distinct(
      .data$ProteinGroupID[
        .data$shared_core_sensitivity_eligible %in% TRUE
      ]
    ),
    .groups = "drop"
  ) |>
  dplyr::rename(marker_class = "fidelity_marker_class")
direction_summary <- protein_direction |>
  dplyr::group_by(.data$marker_class) |>
  dplyr::summarise(
    quantitatively_comparable_proteins = sum(
      .data$n_reliably_comparable_offtarget_datasets > 0L
    ),
    proteins_reliably_detected_only_in_intended_compartment = sum(
      .data$all_offtargets_not_reliably_detected
    ),
    fraction_intended_highest = mean(
      .data$intended_is_highest, na.rm = TRUE
    ),
    median_intended_margin_log2 = stats::median(
      .data$intended_margin_log2, na.rm = TRUE
    ),
    intended_margin_q1_log2 = stats::quantile(
      .data$intended_margin_log2, 0.25, na.rm = TRUE, names = FALSE
    ),
    intended_margin_q3_log2 = stats::quantile(
      .data$intended_margin_log2, 0.75, na.rm = TRUE, names = FALSE
    ),
    median_detection_advantage = stats::median(
      .data$detection_advantage, na.rm = TRUE
    ),
    .groups = "drop"
  )
subpanel_direction <- protein_direction |>
  dplyr::group_by(
    .data$marker_class, .data$display_biological_subpanel
  ) |>
  dplyr::summarise(
    subpanel_fraction_intended_highest =
      mean(.data$intended_is_highest, na.rm = TRUE),
    subpanel_median_intended_margin_log2 =
      stats::median(.data$intended_margin_log2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$marker_class) |>
  dplyr::summarise(
    subpanel_balanced_fraction_intended_highest =
      stats::median(.data$subpanel_fraction_intended_highest, na.rm = TRUE),
    subpanel_balanced_median_intended_margin_log2 =
      stats::median(.data$subpanel_median_intended_margin_log2, na.rm = TRUE),
    .groups = "drop"
  )
loao_summary <- loao_direction |>
  dplyr::group_by(.data$marker_class) |>
  dplyr::summarise(
    loao_fraction_protein_results_intended_highest =
      mean(.data$intended_is_highest, na.rm = TRUE),
    .groups = "drop"
  )
class_summary <- requested_counts |>
  dplyr::left_join(eligible_counts, by = "marker_class") |>
  dplyr::left_join(direction_summary, by = "marker_class") |>
  dplyr::left_join(subpanel_direction, by = "marker_class") |>
  dplyr::left_join(loao_summary, by = "marker_class") |>
  dplyr::mutate(
    analysis_type = "descriptive_validation_n_3_CON_animals",
    inferential_p_values = FALSE
  )

method_comparison <- marker_scores$scores |>
  dplyr::select(
    "dataset", "AnimalID", "marker_class",
    primary_median_centered_log2 = "median_centered_log2",
    primary_subpanel_balanced_median_centered_log2 =
      "subpanel_balanced_median_centered_log2"
  ) |>
  dplyr::left_join(
    legacy_scores$scores |>
      dplyr::select(
        "dataset", "AnimalID", "marker_class",
        legacy_shared_core_mean_robust_z = "mean_robust_z_score",
        legacy_shared_core_median_robust_z = "median_robust_z_score",
        legacy_shared_core_subpanel_balanced_mean_robust_z =
          "subpanel_balanced_mean_robust_z"
      ),
    by = c("dataset", "AnimalID", "marker_class"),
    relationship = "one-to-one"
  )

inclusion_audit <- con$meta |>
  dplyr::transmute(
    Sample = .data[[con$columns$sample]],
    dataset = .data[[con$columns$dataset]],
    AnimalID = .data[[con$columns$animal]],
    ExpGroup_raw = .data[[con$columns$group]],
    resolved_group = ca_normalize_group(.data[[con$columns$group]]),
    region = .data[[con$columns$region]],
    layer = .data[[con$columns$layer]],
    hemisphere = ca_normalize_hemisphere(.data[[con$columns$hemisphere]]),
    hemisphere_source_column = con$columns$hemisphere,
    plate = if (!is.na(con$columns$plate)) {
      .data[[con$columns$plate]]
    } else {
      NA_character_
    },
    included_primary_CON = TRUE
  )
normalization_audit <- data.frame(
  Sample = names(reconstructed$normalization_offsets),
  joint_normalization_offset_log2 =
    as.numeric(reconstructed$normalization_offsets),
  normalization_target_median_log2 =
    reconstructed$normalization_target_median,
  observed_primary_values_reproduced_to_tolerance = TRUE,
  validation_tolerance = 1e-12,
  stringsAsFactors = FALSE
)

contract_columns <- list(
  quantitative_value_status = "observed_nonimputed_or_not_applicable",
  normalization_source =
    "joint_shared_core_sample_median_offsets_applied_to_raw_positive_log2",
  feature_universe =
    "canonical_joint_qc_eligible_noncontaminant_external_fidelity_registry",
  animal_detection_threshold = paste0(">=", min_con_animals, "_of_3_CON_animals"),
  region_threshold = paste0(">=", min_regions, "_of_4_regions_per_valid_hemisphere"),
  hemisphere_handling = median_agg$hemisphere_handling,
  marker_source_policy =
    "external_SynGO_GO_MGI_Allen_with_documented_fallbacks",
  selection_policy =
    "external_rank_plus_intended_detection_plus_subpanel_balance_no_observed_effect"
)
add_contract <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(contract_columns)) {
    if (!nm %in% names(x)) x[[nm]] <- contract_columns[[nm]]
  }
  x
}

source_tables <- list(
  v2_00_joint_normalization_reconstruction_audit = add_contract(normalization_audit),
  v2_01_CON_sample_animal_hemisphere_inclusion = add_contract(inclusion_audit),
  v2_02_hemisphere_level_abundance_detection = add_contract(hemisphere_abundance),
  v2_03_animal_level_abundance = add_contract(dplyr::bind_rows(
    dplyr::mutate(animal_median, aggregation = "primary_median_hierarchy"),
    dplyr::mutate(animal_mean, aggregation = "mean_hierarchy_sensitivity")
  )),
  v2_04_all_marker_canonical_mapping_provenance = add_contract(marker_eligibility),
  v2_05_marker_by_dataset_detection_abundance = add_contract(marker_dataset_audit),
  v2_06_protein_intended_direction = add_contract(protein_direction),
  v2_07_animal_marker_class_summaries = add_contract(marker_scores$scores),
  v2_07a_animal_subpanel_summaries = add_contract(marker_scores$subpanel_scores),
  v2_08_marker_class_descriptive_summaries = add_contract(class_summary),
  v2_09_legacy_shared_core_robust_z_sensitivity = add_contract(legacy_scores$scores),
  v2_10_leave_one_animal_out_class_stability = add_contract(loao_class_stability),
  v2_10a_leave_one_animal_out_protein_direction = add_contract(loao_direction),
  v2_11_display_selection_provenance = add_contract(display$selected),
  v2_11a_display_selection_eligible_alternatives = add_contract(display$alternatives),
  v2_12_rank_abundance_data = add_contract(rank_data),
  v2_13_left_right_hemisphere_pairs = add_contract(left_right_pairs),
  v2_13a_left_right_hemisphere_concordance = add_contract(left_right_summary),
  v2_14_both_hemispheres_valid_sensitivity = add_contract(both_hemisphere_scores$scores),
  v2_15_exact_27_draw_descriptive_bootstrap = add_contract(bootstrap_class_scores),
  v2_16_primary_vs_legacy_method_comparison = add_contract(method_comparison),
  v2_17_marker_detection_dot_heatmap_source = add_contract(dot_source)
)
source_files <- paste0(names(source_tables), ".csv")
expected_outputs <- c(
  file.path(roots$source_data, source_files),
  file.path(roots$figures, rendering_figure_files),
  file.path(roots$reports, c(
    "v2_control_compartment_abundance_interpretation.md",
    "v2_rendering_figure_legends.md"
  )),
  file.path(roots$logs, c(
    "v2_run_manifest.yml", "v2_rendering_manifest.yml"
  ))
)
invisible(lapply(roots, dir_create))
for (i in seq_along(source_tables)) {
  qc_write_csv(
    source_tables[[i]],
    file.path(roots$source_data, source_files[[i]])
  )
}
write_control_rendering_outputs_v2(
  rendering,
  roots,
  inputs_used = c(
    render_inputs[["dot_source"]],
    render_inputs[["rank_source"]],
    render_inputs[["display_selection"]]
  ),
  allow_overwrite = allow_overwrite
)

writeLines(
  c(
    "# Control-only compartment abundance validation v2",
    "",
    "- Biological unit: animal (n = 3 CON animals), not ROI, layer, region or hemisphere.",
    "- Spatial preparations: neuronal soma, neuronal neuropil and microglia/PVM-enriched ROI; not purified cell types.",
    "- Quantitative matrix: raw-positive log2 values reconstructed from the joint bundle, normalized with the joint shared-core sample offsets, with observed missingness retained.",
    paste0(
      "- Hierarchy: technical replicates -> layers -> regions within hemisphere -> equal valid hemispheres -> animal. A hemisphere requires >=",
      min_regions,
      " regions; one valid hemisphere is allowed in the primary analysis and two are required in sensitivity."
    ),
    paste0(
      "- Primary marker eligibility: external fidelity registry, eligible canonical non-contaminant ProteinGroupID, no conflicting class, and >=",
      min_con_animals,
      "/3 valid animals in the intended compartment. Joint-shared-core membership is sensitivity metadata only."
    ),
    "- Primary class evidence: median within-protein centered log2 abundance and median of subpanel medians.",
    "- Off-target nondetection is retained as detection evidence and is never converted to zero or an imputed abundance contrast.",
    "- The exact 27 ordered animal-bootstrap draws and leave-one-animal-out results are descriptive stability analyses, not confidence intervals or inferential tests.",
    "- No ROI-level or sample-level p-values are reported."
  ),
  file.path(
    roots$reports,
    "v2_control_compartment_abundance_interpretation.md"
  )
)
write_run_manifest(
  file.path(roots$logs, "v2_run_manifest.yml"),
  inputs = as.list(inputs),
  outputs = list(
    source_files = file.path(roots$source_data, source_files),
    figure_files = file.path(roots$figures, rendering_figure_files),
    interpretation = file.path(
      roots$reports,
      "v2_control_compartment_abundance_interpretation.md"
    )
  ),
  parameters = list(
    contract_version = "control_compartment_abundance_v2",
    primary_group = "CON",
    resolved_group_column = con$resolved_group_column,
    included_samples = con$included_samples,
    included_animals = con$included_animals,
    hemisphere_column = median_agg$hemisphere_column,
    hemisphere_handling = median_agg$hemisphere_handling,
    min_regions_per_valid_hemisphere = min_regions,
    primary_min_CON_animals = min_con_animals,
    strict_min_CON_animals = strict_con_animals,
    primary_class_statistic = "median_within_protein_centered_log2",
    subpanel_balance_statistic = "median_of_subpanel_medians",
    legacy_sensitivity = "shared_core_robust_z_explicitly_nonprimary",
    exact_bootstrap = "all_27_ordered_descriptive_animal_draws",
    inferential_p_values = FALSE,
    normalized_nonimputed_reconstruction =
      "validated_against_joint_compartment_qc_v1_tolerance_1e-12"
  ),
  notes = paste(
    "Authoritative 04e control-only cross-compartment marker detection and abundance.",
    "Existing v1 outputs are retained and not overwritten.",
    "04d sample-level inferential claim path is deprecated."
  )
)

message("Control-only compartment abundance v2 workflow complete.")
