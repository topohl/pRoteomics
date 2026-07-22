#!/usr/bin/env Rscript
# Control-only compartment-marker fidelity and relative-abundance validation.
# This is a rendering-oriented global consumer of the completed joint-QC bundle.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "compartment_abundance_utils.R"))
source(repo_path("R", "joint_compartment_qc_plotting.R"))
source(repo_path("R", "control_compartment_abundance_rendering.R"))

global_arg <- tolower(script_arg_value("--dataset", "global"))
if (!global_arg %in% c("all", "global")) stop("This is a global script; use --dataset all or global.", call. = FALSE)
dry_run <- is_dry_run()
render_only <- "--render-only" %in% commandArgs(trailingOnly = TRUE)

as_integer_contract <- function(env, default, lower, upper) {
  value <- suppressWarnings(as.integer(Sys.getenv(env, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower || value > upper) {
    stop(env, " must be an integer from ", lower, " to ", upper, ".", call. = FALSE)
  }
  value
}

min_regions <- as_integer_contract("PROTEOMICS_CONTROL_ABUNDANCE_MIN_REGIONS", 3L, 1L, 4L)
min_con_animals <- as_integer_contract("PROTEOMICS_CONTROL_ABUNDANCE_MIN_CON_ANIMALS", 2L, 1L, 3L)
strict_con_animals <- as_integer_contract("PROTEOMICS_CONTROL_ABUNDANCE_STRICT_CON_ANIMALS", 3L, 1L, 3L)
min_score_proteins <- as_integer_contract("PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_PROTEINS", 3L, 1L, 10000L)
min_score_fraction <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_FRACTION", unset = "0.5")))
if (!is.finite(min_score_fraction) || min_score_fraction < 0 || min_score_fraction > 1) {
  stop("PROTEOMICS_CONTROL_ABUNDANCE_MIN_SCORE_FRACTION must be between 0 and 1.", call. = FALSE)
}

processed_root <- Sys.getenv(
  "PROTEOMICS_JOINT_QC_PROCESSED_DIR",
  unset = path_processed("01_preprocessing", "joint_compartment_qc", "global")
)
inputs <- c(
  joint_bundle = file.path(processed_root, "joint_compartment_qc_matrices.rds"),
  marker_registry = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"),
  publication_labels = repo_path("config", "marker_panels", "control_compartment_publication_labels.csv")
)
substep <- "04e_control_compartment_abundance_publication_figures"
roots <- list(
  source_data = path_results("source_data", "03_qc_exploration", substep, "global"),
  figures = path_results("figures", "03_qc_exploration", substep, "global"),
  reports = path_results("reports", "03_qc_exploration", substep, "global"),
  logs = path_results("logs", "03_qc_exploration", substep, "global")
)
render_inputs <- c(
  dot_source = file.path(roots$source_data, "07a_control_marker_dot_heatmap_source.csv"),
  rank_source = file.path(roots$source_data, "04_rank_abundance_with_percentile.csv"),
  marker_audit = file.path(roots$source_data, "05_marker_matches_and_exclusions.csv")
)

if (dry_run) {
  for (nm in names(inputs)) {
    cat("[DRY-RUN] ", nm, ": ", normalizePath(inputs[[nm]], winslash = "/", mustWork = FALSE),
        " [", if (file.exists(inputs[[nm]])) "PASS" else "FAIL", "]\n", sep = "")
  }
  cat("[DRY-RUN] group: CON (ExpGroup 1/CON only)\n")
  cat("[DRY-RUN] integer eligibility: >=", min_regions, " regions per animal; >=", min_con_animals,
      "/3 animals primary; ", strict_con_animals, "/3 strict\n", sep = "")
  cat("[DRY-RUN] source data: ", normalizePath(roots$source_data, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("[DRY-RUN] figures: ", normalizePath(roots$figures, winslash = "/", mustWork = FALSE), "\n", sep = "")
  if (render_only) {
    for (nm in names(render_inputs)) {
      cat("[DRY-RUN] render input ", nm, ": ", normalizePath(render_inputs[[nm]], winslash = "/", mustWork = FALSE),
          " [", if (file.exists(render_inputs[[nm]])) "PASS" else "FAIL", "]\n", sep = "")
    }
  }
  required_dry_inputs <- if (render_only) c(inputs[["publication_labels"]], render_inputs) else inputs
  quit(status = if (all(file.exists(required_dry_inputs))) 0 else 1, save = "no")
}

required_inputs <- if (render_only) inputs["publication_labels"] else inputs
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required input(s): ", paste(names(missing_inputs), collapse = ", "), call. = FALSE)
required_packages <- c("dplyr", "ggplot2", "ggrepel", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)

rendering_source_files <- c(
  "10_configured_marker_display.csv",
  "10a_marker_display_selection_provenance.csv",
  "10b_marker_display_eligible_alternatives.csv",
  "10c_marker_tile_display_source.csv",
  "10d_rank_representative_display_source.csv",
  "10e_terminal_rank_drop_audit.csv",
  "control_compartment_rendering_source_data.xlsx"
)
rendering_figure_files <- c(
  "control_compartment_marker_dot_heatmap_89mm.svg", "control_compartment_marker_dot_heatmap_89mm.pdf",
  "control_rank_abundance_extended_data_183mm.svg", "control_rank_abundance_extended_data_183mm.pdf",
  "control_rank_abundance_extended_data_183mm.png"
)

write_control_rendering_outputs <- function(rendering, roots, inputs_used, allow_overwrite = FALSE) {
  legends <- ca_control_rendering_legends()
  paths <- c(
    file.path(roots$source_data, rendering_source_files),
    file.path(roots$figures, rendering_figure_files),
    file.path(roots$reports, "rendering_figure_legends.md"),
    file.path(roots$logs, "rendering", "rendering_manifest.yml")
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) && !allow_overwrite) {
    stop(
      "Refusing to overwrite existing rendering output(s): ",
      paste(normalizePath(existing, winslash = "/", mustWork = FALSE), collapse = ", "),
      ". Set PROTEOMICS_CONTROL_ABUNDANCE_RENDER_ALLOW_OVERWRITE=true only for an intentional rendering regeneration.",
      call. = FALSE
    )
  }
  invisible(lapply(roots, dir_create))
  qc_write_csv(rendering$config, file.path(roots$source_data, rendering_source_files[[1]]))
  qc_write_csv(rendering$provenance, file.path(roots$source_data, rendering_source_files[[2]]))
  qc_write_csv(rendering$alternatives, file.path(roots$source_data, rendering_source_files[[3]]))
  qc_write_csv(rendering$dot, file.path(roots$source_data, rendering_source_files[[4]]))
  qc_write_csv(rendering$rank_representatives, file.path(roots$source_data, rendering_source_files[[5]]))
  qc_write_csv(rendering$terminal_drop, file.path(roots$source_data, rendering_source_files[[6]]))
  qc_write_xlsx(list(
    configured_markers = rendering$config,
    selection_provenance = rendering$provenance,
    eligible_alternatives = rendering$alternatives,
    tile_display = rendering$dot,
    rank_display = rendering$rank_representatives,
    terminal_drop = rendering$terminal_drop
  ), file.path(roots$source_data, rendering_source_files[[7]]))

  p_dot <- ca_build_control_marker_tile_heatmap(rendering)
  p_rank <- ca_build_control_rank_abundance_plot(rendering)
  ca_save_vector_figure(p_dot, roots$figures, "control_compartment_marker_dot_heatmap_89mm", 89, 92)
  ca_save_vector_figure(p_rank, roots$figures, "control_rank_abundance_extended_data_183mm", 120, 35)
  ggplot2::ggsave(
    file.path(roots$figures, "control_rank_abundance_extended_data_183mm.png"),
    p_rank, width = 120, height = 35, units = "mm", dpi = 300,
    device = "png", limitsize = FALSE, bg = "white"
  )

  writeLines(c(
    "# Main panel figure legend", "", legends$main, "",
    "# Extended Data figure legend", "", legends$extended
  ), file.path(roots$reports, "rendering_figure_legends.md"))
  write_run_manifest(
    file.path(roots$logs, "rendering", "rendering_manifest.yml"),
    inputs = as.list(inputs_used),
    outputs = list(
      source_files = file.path(roots$source_data, rendering_source_files),
      figure_files = file.path(roots$figures, rendering_figure_files),
      figure_legends = file.path(roots$reports, "rendering_figure_legends.md")
    ),
    parameters = list(
      execution_mode = "render_only_from_completed_04e_source_tables",
      configured_markers = nrow(rendering$config),
      deterministic_rank_labels = sum(rendering$rank_representatives$display_rank_label),
      marker_point_layer = "configured_representatives_only",
      tile_size_encoding = "none_all_configured_proteins_3_of_3_valid_animals",
      robust_z_display_cap = 3,
      robust_z_source_values = "uncapped",
      terminal_drop_definition = unique(rendering$terminal_drop$terminal_drop_definition)
    ),
    notes = paste(
      "Rendering-only regeneration from existing analytical source tables.",
      "No sample inclusion, normalization, detection, aggregation, rank, robust-z, or analytical source values were recalculated."
    )
  )
  invisible(list(dot = p_dot, rank = p_rank, outputs = paths))
}

if (render_only) {
  missing_render_inputs <- render_inputs[!file.exists(render_inputs)]
  if (length(missing_render_inputs)) {
    stop("Missing completed 04e rendering input(s): ", paste(names(missing_render_inputs), collapse = ", "), call. = FALSE)
  }
  display_config <- utils::read.csv(inputs[["publication_labels"]], check.names = FALSE, stringsAsFactors = FALSE)
  rendering <- ca_prepare_control_rendering_sources(
    dot_source = qc_read_table(render_inputs[["dot_source"]]),
    rank_table = qc_read_table(render_inputs[["rank_source"]]),
    marker_audit = qc_read_table(render_inputs[["marker_audit"]]),
    config = display_config
  )
  allow_render_overwrite <- ca_as_logical(Sys.getenv(
    "PROTEOMICS_CONTROL_ABUNDANCE_RENDER_ALLOW_OVERWRITE", unset = "false"
  ))
  write_control_rendering_outputs(
    rendering, roots,
    inputs_used = c(publication_labels = inputs[["publication_labels"]], render_inputs),
    allow_overwrite = allow_render_overwrite
  )
  terminal_summary <- rendering$terminal_drop |>
    dplyr::distinct(.data$dataset, .data$terminal_segment_n, .data$terminal_segment_start_rank,
                    .data$terminal_segment_initiating_drop_log2)
  message(
    "04e rendering-only regeneration complete; terminal segments: ",
    paste0(terminal_summary$dataset, "=", terminal_summary$terminal_segment_n, collapse = ", "), "."
  )
  quit(status = 0, save = "no")
}

bundle <- readRDS(inputs[["joint_bundle"]])
reconstructed <- ca_validate_bundle_reconstruction(bundle)
shared_ids <- reconstructed$primary_feature_ids
analysis_mat <- reconstructed$normalized_nonimputed
meta <- reconstructed$metadata
if (!identical(colnames(analysis_mat), as.character(meta[[reconstructed$columns$sample]]))) {
  stop("Reconstructed non-imputed values are not aligned to metadata.", call. = FALSE)
}

con <- ca_filter_expression_group(analysis_mat, meta, "CON")
if (!identical(con$resolved_group_column, "ExpGroup")) {
  stop("The primary analysis must resolve the experimental group from ExpGroup, not an anatomical group column.", call. = FALSE)
}
animals <- sort(unique(as.character(con$meta[[con$columns$animal]])))
if (length(animals) != 3L) stop("Exact stability analysis requires exactly three CON animals; found ", length(animals), ".", call. = FALSE)
datasets <- c("neuron_soma", "neuron_neuropil", "microglia")
if (!setequal(unique(as.character(con$meta[[con$columns$dataset]])), datasets)) {
  stop("CON metadata does not contain exactly the three required compartment datasets.", call. = FALSE)
}

primary_agg <- ca_aggregate_abundance(con$mat, con$meta, "mean", min_regions)
median_agg <- ca_aggregate_abundance(con$mat, con$meta, "median", min_regions)
all_agg <- ca_aggregate_abundance(analysis_mat, meta, "mean", min_regions)
all_animals <- length(unique(as.character(meta[[reconstructed$columns$animal]])))
all_min_animals <- as_integer_contract(
  "PROTEOMICS_CONTROL_ABUNDANCE_ALL_GROUP_MIN_ANIMALS",
  ceiling((2 / 3) * all_animals), 1L, all_animals
)

detection <- ca_detection_audit(primary_agg, min_con_animals, strict_con_animals)
detection$in_joint_shared_core <- detection$ProteinGroupID %in% shared_ids
detection$nonimputed_missing_sample_count <- detection$total_sample_count - detection$observed_sample_count
detection$nonimputed_missing_fraction <- detection$nonimputed_missing_sample_count / detection$total_sample_count
detection$imputed_sample_count_in_values_used <- 0L
detection$imputed_fraction <- 0
detection$upstream_primary_imputation_fraction <- ifelse(
  detection$in_joint_shared_core, detection$nonimputed_missing_fraction, NA_real_
)
detection$imputation_reason <- ifelse(
  detection$in_joint_shared_core,
  "shared_core_upstream_imputation_mask_audited_but_nonimputed_values_used",
  "outside_shared_core_no_imputed_quantitative_values_used"
)
detection$observed_detection_source <- "raw_positive_joint_bundle_before_imputation"
detection$valid_animal_contract <- paste0(">=", min_regions, "_of_4_regions")
all_detection <- ca_detection_audit(all_agg, all_min_animals, all_animals)

primary_estimates <- ca_dataset_estimates(primary_agg, "primary_region_balanced_mean_log2")
median_estimates <- ca_dataset_estimates(median_agg, "median_hierarchy_log2")
all_estimates <- ca_dataset_estimates(all_agg, "all_group_region_balanced_mean_log2")
primary_rank <- ca_rank_table(
  primary_estimates, detection, "primary_eligible", "primary_region_balanced_mean_log2",
  "CON_primary_2_of_3", "CON_region_balanced_mean"
)
median_rank <- ca_rank_table(
  median_estimates, detection, "primary_eligible", "median_hierarchy_log2",
  "CON_primary_2_of_3", "CON_median_hierarchy"
)
strict_rank <- ca_rank_table(
  primary_estimates, detection, "strict_eligible", "primary_region_balanced_mean_log2",
  "CON_strict_3_of_3", "CON_region_balanced_mean_strict"
)
all_rank <- ca_rank_table(
  all_estimates, all_detection, "primary_eligible", "all_group_region_balanced_mean_log2",
  paste0("all_groups_", all_min_animals, "_of_", all_animals), "all_groups_region_balanced_mean"
)

feature_annotation <- bundle$feature_table[, intersect(c(
  "ProteinGroupID", "FeatureDisplayLabel", "original_identifier", "member_accessions",
  "member_gene_symbols", "official_gene_symbol", "protein_group_ambiguity_class",
  "gene_level_claim_allowed"
), names(bundle$feature_table)), drop = FALSE]
rank_data <- dplyr::bind_rows(primary_rank, median_rank, strict_rank, all_rank) |>
  dplyr::left_join(feature_annotation, by = "ProteinGroupID", relationship = "many-to-one")

registry <- utils::read.csv(inputs[["marker_registry"]], check.names = FALSE, stringsAsFactors = FALSE)
marker_matches <- ca_match_fidelity_markers(registry, bundle$feature_table, shared_ids)
marker_map_raw <- marker_matches |>
  dplyr::filter(.data$primary_score_eligible %in% TRUE) |>
  dplyr::transmute(
    ProteinGroupID, marker_class = fidelity_marker_class, marker_gene,
    marker_panel, fidelity_subpanel, primary_score_eligible,
    source_type, source_name, source_reference, source_term_or_category,
    evidence_level, selection_rule, confidence, use_for, notes
  )
marker_map <- ca_prepare_marker_map(marker_map_raw)

primary_rank_annotated <- ca_join_marker_annotations(primary_rank, marker_map)
rank_data <- dplyr::bind_rows(
  primary_rank_annotated,
  rank_data |> dplyr::filter(.data$RankGroup != "CON_region_balanced_mean")
) |>
  dplyr::left_join(
    feature_annotation,
    by = "ProteinGroupID", relationship = "many-to-one", suffix = c("", ".annotation")
  )
for (nm in names(feature_annotation)) {
  alt <- paste0(nm, ".annotation")
  if (alt %in% names(rank_data)) rank_data[[alt]] <- NULL
}

animal_primary <- ca_animal_abundance_table(primary_agg, "animal_log2_abundance")
animal_median <- ca_animal_abundance_table(median_agg, "animal_log2_abundance")
animal_abundance <- dplyr::bind_rows(
  dplyr::mutate(animal_primary, aggregation = "primary_region_balanced_mean"),
  dplyr::mutate(animal_median, aggregation = "median_hierarchy_sensitivity")
)
robust <- ca_robust_standardize(animal_primary[animal_primary$ProteinGroupID %in% shared_ids, , drop = FALSE])
marker_scores <- ca_score_markers(robust, marker_map_raw, min_score_proteins, min_score_fraction)

display_config <- utils::read.csv(inputs[["publication_labels"]], check.names = FALSE, stringsAsFactors = FALSE)
labels <- display_config
names(labels)[names(labels) == "official_gene_symbol"] <- "publication_gene_symbol"
labels$marker_gene_key <- toupper(trimws(labels$publication_gene_symbol))
label_matches <- marker_matches |>
  dplyr::filter(.data$primary_score_eligible %in% TRUE, .data$gene_claim_eligible %in% TRUE) |>
  dplyr::inner_join(labels, by = c("marker_gene_key", "fidelity_marker_class" = "marker_class"), relationship = "many-to-many") |>
  dplyr::distinct(.data$fidelity_marker_class, .data$publication_gene_symbol, .data$ProteinGroupID, .keep_all = TRUE)

mean_median_concordance <- ca_rank_concordance(primary_rank, median_rank, "primary_mean", "median_hierarchy")
universe_concordance <- ca_rank_concordance(primary_rank, strict_rank, "2_of_3", "3_of_3")
all_group_concordance <- ca_rank_concordance(primary_rank, all_rank, "CON", "all_groups")
sensitivity_summary <- dplyr::bind_rows(
  mean_median_concordance$summary,
  universe_concordance$summary,
  all_group_concordance$summary
)
sensitivity_by_protein <- dplyr::bind_rows(
  dplyr::mutate(mean_median_concordance$by_protein, comparison = "primary_mean_vs_median_hierarchy"),
  dplyr::mutate(universe_concordance$by_protein, comparison = "2_of_3_vs_3_of_3"),
  dplyr::mutate(all_group_concordance$by_protein, comparison = "CON_vs_all_groups")
)
marker_ids <- sort(unique(marker_map$ProteinGroupID))
top_marker_stability <- all_group_concordance$by_protein |>
  dplyr::filter(.data$ProteinGroupID %in% marker_ids) |>
  dplyr::mutate(
    CON_top_decile = .data$RankPercentile_a >= 0.9,
    all_group_top_decile = .data$RankPercentile_b >= 0.9,
    top_decile_stable = .data$CON_top_decile == .data$all_group_top_decile,
    rank_percentile_change = .data$RankPercentile_b - .data$RankPercentile_a
  ) |>
  dplyr::left_join(marker_map, by = "ProteinGroupID", relationship = "many-to-many")
sensitivity_summary <- dplyr::bind_rows(
  sensitivity_summary,
  top_marker_stability |>
    dplyr::group_by(.data$dataset, .data$marker_class) |>
    dplyr::summarise(
      comparison = "CON_vs_all_groups_top_marker_stability",
      n_common_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      spearman_rank = NA_real_, spearman_rank_percentile = NA_real_,
      top_decile_stability_fraction = mean(.data$top_decile_stable),
      .groups = "drop"
    )
)

rank_vector <- function(values, protein_ids) {
  keep <- is.finite(values)
  out <- data.frame(ProteinGroupID = protein_ids, abundance = values, Rank = NA_integer_, RankPercentile = NA_real_)
  ord <- which(keep)[order(-values[keep], protein_ids[keep], method = "radix")]
  n <- length(ord)
  if (n) {
    out$Rank[ord] <- seq_len(n)
    out$RankPercentile[ord] <- if (n == 1L) 1 else 1 - (seq_len(n) - 1) / (n - 1)
  }
  out
}

primary_universe <- split(primary_rank$ProteinGroupID, primary_rank$dataset)
animal_matrix_for <- function(dataset) {
  x <- animal_primary[animal_primary$dataset == dataset, c("ProteinGroupID", "AnimalID", "animal_log2_abundance")]
  proteins <- primary_universe[[dataset]]
  x <- x[x$ProteinGroupID %in% proteins & x$AnimalID %in% animals, , drop = FALSE]
  ans <- matrix(NA_real_, length(proteins), length(animals), dimnames = list(proteins, animals))
  ans[cbind(match(x$ProteinGroupID, proteins), match(x$AnimalID, animals))] <- x$animal_log2_abundance
  ans
}
animal_matrices <- stats::setNames(lapply(datasets, animal_matrix_for), datasets)
draws <- ca_exact_paired_bootstrap_draws(animals)
bootstrap_parts <- list()
k <- 0L
for (i in seq_len(nrow(draws))) {
  drawn <- as.character(unlist(draws[i, c("draw_1", "draw_2", "draw_3")], use.names = FALSE))
  for (ds in datasets) {
    block <- animal_matrices[[ds]][, drawn, drop = FALSE]
    abundance <- rowMeans(block, na.rm = TRUE)
    abundance[rowSums(is.finite(block)) == 0L] <- NA_real_
    ranked <- rank_vector(abundance, rownames(block))
    ranked <- ranked[ranked$ProteinGroupID %in% marker_ids, , drop = FALSE]
    ranked$dataset <- ds
    ranked$bootstrap_draw_id <- draws$bootstrap_draw_id[[i]]
    ranked$draw_1 <- drawn[[1]]; ranked$draw_2 <- drawn[[2]]; ranked$draw_3 <- drawn[[3]]
    k <- k + 1L; bootstrap_parts[[k]] <- ranked
  }
}
bootstrap_exact <- dplyr::bind_rows(bootstrap_parts) |>
  dplyr::left_join(marker_map, by = "ProteinGroupID", relationship = "many-to-many")
bootstrap_summary <- bootstrap_exact |>
  dplyr::group_by(.data$dataset, .data$ProteinGroupID, .data$marker_class) |>
  dplyr::summarise(
    n_ordered_draws = dplyr::n_distinct(.data$bootstrap_draw_id),
    n_draws_ranked = sum(is.finite(.data$RankPercentile)),
    median_rank_percentile = stats::median(.data$RankPercentile, na.rm = TRUE),
    min_rank_percentile = min(.data$RankPercentile, na.rm = TRUE),
    max_rank_percentile = max(.data$RankPercentile, na.rm = TRUE),
    descriptive_range = .data$max_rank_percentile - .data$min_rank_percentile,
    interval_interpretation = "exact_descriptive_stability_range_not_inferential_CI",
    .groups = "drop"
  )

loao_cases <- ca_leave_one_animal_out_cases(animals)
loao_rank_parts <- list()
loao_score_parts <- list()
for (i in seq_len(nrow(loao_cases))) {
  omitted <- loao_cases$omitted_AnimalID[[i]]
  for (ds in datasets) {
    block <- animal_matrices[[ds]][, setdiff(animals, omitted), drop = FALSE]
    abundance <- rowMeans(block, na.rm = TRUE)
    abundance[rowSums(is.finite(block)) == 0L] <- NA_real_
    ranked <- rank_vector(abundance, rownames(block))
    ranked <- ranked[ranked$ProteinGroupID %in% marker_ids, , drop = FALSE]
    ranked$dataset <- ds; ranked$omission_id <- loao_cases$omission_id[[i]]; ranked$omitted_AnimalID <- omitted
    loao_rank_parts[[length(loao_rank_parts) + 1L]] <- ranked
  }
  reduced <- animal_primary[animal_primary$AnimalID != omitted, , drop = FALSE]
  reduced_robust <- ca_robust_standardize(reduced)
  reduced_scores <- ca_score_markers(reduced_robust, marker_map_raw, min_score_proteins, min_score_fraction)$scores
  reduced_scores$omission_id <- loao_cases$omission_id[[i]]
  reduced_scores$omitted_AnimalID <- omitted
  loao_score_parts[[i]] <- reduced_scores
}
loao_ranks <- dplyr::bind_rows(loao_rank_parts) |>
  dplyr::left_join(marker_map, by = "ProteinGroupID", relationship = "many-to-many")
loao_scores <- dplyr::bind_rows(loao_score_parts)

region_long <- ca_matrix_long(
  primary_agg$region$values, primary_agg$region$group_meta, "animal_region_log2"
)
region_influence <- region_long |>
  dplyr::filter(.data$ProteinGroupID %in% marker_ids) |>
  dplyr::group_by(.data$dataset, .data$ProteinGroupID, .data$region) |>
  dplyr::summarise(region_mean_log2 = mean(.data$animal_region_log2, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(
    primary_estimates |> dplyr::select("dataset", "ProteinGroupID", "primary_region_balanced_mean_log2"),
    by = c("dataset", "ProteinGroupID")
  ) |>
  dplyr::mutate(region_deviation_log2 = .data$region_mean_log2 - .data$primary_region_balanced_mean_log2) |>
  dplyr::left_join(marker_map, by = "ProteinGroupID", relationship = "many-to-many")

inclusion_audit <- con$meta |>
  dplyr::transmute(
    Sample = .data[[con$columns$sample]], dataset = .data[[con$columns$dataset]],
    AnimalID = .data[[con$columns$animal]], ExpGroup_raw = .data[[con$columns$group]],
    resolved_group = ca_normalize_group(.data[[con$columns$group]]),
    region = .data[[con$columns$region]], layer = .data[[con$columns$layer]],
    ReplicateGroup = if (!is.na(con$columns$replicate)) .data[[con$columns$replicate]] else NA_character_,
    plate = if (!is.na(con$columns$plate)) .data[[con$columns$plate]] else NA_character_,
    included_primary_CON = TRUE
  ) |>
  dplyr::group_by(.data$dataset) |>
  dplyr::mutate(dataset_included_samples = dplyr::n(), dataset_included_animals = dplyr::n_distinct(.data$AnimalID)) |>
  dplyr::ungroup()

dot_source <- robust |>
  dplyr::filter(.data$ProteinGroupID %in% label_matches$ProteinGroupID) |>
  dplyr::group_by(.data$dataset, .data$ProteinGroupID) |>
  dplyr::summarise(
    median_robust_standardized_abundance = stats::median(.data$robust_z, na.rm = TRUE),
    zero_MAD_excluded = all(.data$zero_MAD_excluded), .groups = "drop"
  ) |>
  dplyr::left_join(
    detection |> dplyr::select("dataset", "ProteinGroupID", raw_CON_valid_animal_detection_fraction = "valid_animal_fraction"),
    by = c("dataset", "ProteinGroupID")
  ) |>
  dplyr::left_join(
    label_matches |> dplyr::select("ProteinGroupID", marker_class = "fidelity_marker_class",
                                   "publication_gene_symbol", "label_order",
                                   rationale = "display_rationale") |> dplyr::distinct(),
    by = "ProteinGroupID", relationship = "many-to-many"
  ) |>
  dplyr::filter(is.finite(.data$median_robust_standardized_abundance))

rendering <- ca_prepare_control_rendering_sources(
  dot_source = dot_source,
  rank_table = rank_data,
  marker_audit = marker_matches,
  config = display_config
)

source_files <- c(
  "01_CON_sample_and_animal_inclusion_audit.csv",
  "02_animal_level_aggregated_abundance.csv",
  "03_CON_raw_detection_audit.csv",
  "04_rank_abundance_with_percentile.csv",
  "05_marker_matches_and_exclusions.csv",
  "06_proteinwise_robust_standardized_abundance.csv",
  "07_animal_level_marker_class_scores.csv",
  "08_rank_sensitivity_and_concordance.csv",
  "09_exact_bootstrap_rank_percentile_summary.csv",
  "09a_exact_bootstrap_rank_percentile_distribution.csv",
  "09b_leave_one_animal_out_rank_distribution.csv",
  "09c_leave_one_animal_out_marker_score_distribution.csv",
  "09d_region_influence_marker_sensitivity.csv",
  "control_compartment_abundance_source_data.xlsx"
)
dot_source_file <- "07a_control_marker_dot_heatmap_source.csv"
figure_files <- c(
  "control_compartment_marker_dot_heatmap_89mm.svg", "control_compartment_marker_dot_heatmap_89mm.pdf",
  "control_rank_abundance_extended_data_183mm.svg", "control_rank_abundance_extended_data_183mm.pdf"
)
expected_outputs <- c(
  file.path(roots$source_data, c(source_files, dot_source_file)),
  file.path(roots$figures, figure_files)
)
allow_overwrite <- ca_as_logical(Sys.getenv("PROTEOMICS_CONTROL_ABUNDANCE_ALLOW_OVERWRITE", unset = "false"))
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !allow_overwrite) {
  stop(
    "Refusing to overwrite existing control-abundance output(s): ",
    paste(normalizePath(existing_outputs, winslash = "/", mustWork = FALSE), collapse = ", "),
    ". Set PROTEOMICS_CONTROL_ABUNDANCE_ALLOW_OVERWRITE=true only for an intentional regeneration.",
    call. = FALSE
  )
}
invisible(lapply(roots, dir_create))

sensitivity_export <- dplyr::bind_rows(
  dplyr::mutate(sensitivity_summary, record_type = "summary"),
  dplyr::mutate(sensitivity_by_protein, record_type = "protein_rank_pair"),
  dplyr::mutate(top_marker_stability, record_type = "marker_top_decile_stability")
)
qc_write_csv(inclusion_audit, file.path(roots$source_data, source_files[[1]]))
qc_write_csv(animal_abundance, file.path(roots$source_data, source_files[[2]]))
qc_write_csv(detection, file.path(roots$source_data, source_files[[3]]))
qc_write_csv(rank_data, file.path(roots$source_data, source_files[[4]]))
qc_write_csv(marker_matches, file.path(roots$source_data, source_files[[5]]))
qc_write_csv(robust, file.path(roots$source_data, source_files[[6]]))
qc_write_csv(marker_scores$scores, file.path(roots$source_data, source_files[[7]]))
qc_write_csv(dot_source, file.path(roots$source_data, dot_source_file))
qc_write_csv(sensitivity_export, file.path(roots$source_data, source_files[[8]]))
qc_write_csv(bootstrap_summary, file.path(roots$source_data, source_files[[9]]))
qc_write_csv(bootstrap_exact, file.path(roots$source_data, source_files[[10]]))
qc_write_csv(loao_ranks, file.path(roots$source_data, source_files[[11]]))
qc_write_csv(loao_scores, file.path(roots$source_data, source_files[[12]]))
qc_write_csv(region_influence, file.path(roots$source_data, source_files[[13]]))
qc_write_xlsx(list(
  CON_inclusion = inclusion_audit,
  animal_abundance = animal_abundance,
  raw_detection = detection,
  rank_abundance = rank_data,
  marker_matches = marker_matches,
  robust_z = robust,
  marker_scores = marker_scores$scores,
  dot_heatmap_source = dot_source,
  rank_sensitivity = sensitivity_export,
  bootstrap_summary = bootstrap_summary,
  bootstrap_exact = bootstrap_exact,
  LOAO_ranks = loao_ranks,
  LOAO_scores = loao_scores,
  region_influence = region_influence
), file.path(roots$source_data, source_files[[14]]))

write_control_rendering_outputs(
  rendering, roots,
  inputs_used = c(
    publication_labels = inputs[["publication_labels"]],
    dot_source = file.path(roots$source_data, dot_source_file),
    rank_source = file.path(roots$source_data, source_files[[4]]),
    marker_audit = file.path(roots$source_data, source_files[[5]])
  ),
  allow_overwrite = allow_overwrite
)

report_lines <- c(
  "# Control-only compartment abundance validation",
  "",
  paste0("- Primary estimand: equal-weight arithmetic mean on the log2 scale across Left/Right, then layers, then CA1/CA2/CA3/DG, then CON animals; animals require >=", min_regions, " regions."),
  paste0("- Primary universe: >=", min_con_animals, "/3 valid CON animals; strict sensitivity: ", strict_con_animals, "/3."),
  "- Median Left/Right, layers, regions and animals is retained only as an aggregation sensitivity analysis.",
  "- Exact ordered paired-animal cluster bootstrap comprises all 27 draws and is descriptive stability, not an inferential confidence interval.",
  "- Compartment scores are descriptive; marker classes contain unequal numbers of proteins and do not have equal precision.",
  "- The three columns denote neuronal soma, neuronal neuropil and microglia/PVM-enriched ROI, not three cell types.",
  "- Abundance is relative jointly normalized log2 abundance, not cell fraction, purity, copy number, total hippocampal abundance or deconvolution."
)
writeLines(report_lines, file.path(roots$reports, "control_compartment_abundance_interpretation.md"))
write_run_manifest(
  file.path(roots$logs, "run_manifest.yml"),
  inputs = as.list(inputs),
  outputs = c(roots, list(source_files = c(source_files, dot_source_file), figure_files = figure_files)),
  parameters = list(
    primary_group = "CON", resolved_group_column = con$resolved_group_column,
    included_samples = con$included_samples, included_animals = con$included_animals,
    min_regions_per_animal = min_regions, primary_min_CON_animals = min_con_animals,
    strict_min_CON_animals = strict_con_animals, all_group_min_animals = all_min_animals,
    marker_score_min_proteins = min_score_proteins, marker_score_min_fraction = min_score_fraction,
    bootstrap = "all_27_ordered_paired_animal_draws", leave_one_out_cases = nrow(loao_cases),
    normalized_nonimputed_reconstruction = "exactly_validated_against_joint_compartment_qc_v1"
  ),
  notes = paste(
    "Control-only compartment-marker fidelity and region-balanced relative protein abundance.",
    "All-group ranks are a separately named sensitivity audit. Bootstrap ranges are descriptive only."
  )
)

message("Control-only compartment abundance publication consumer complete.")
