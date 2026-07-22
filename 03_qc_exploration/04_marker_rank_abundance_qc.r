# ================================================================
# Script: 03_qc_exploration/04_marker_rank_abundance_qc.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required Stage 01 post-filter/imputed quantitative matrix, mouse UniProt mapping, sample metadata; optional manual mappings and marker registry.
# Produces: results/tables/03_qc_exploration/04_marker_rank_abundance_qc/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Marker abundance QC; benefits from qc_global marker registries.
# ================================================================

# Dataset-aware rank-abundance and marker abundance sanity checks.
# Marker panels are abundance/compartment checks, not definitive purity estimates.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "compartment_abundance_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04_marker_rank_abundance_qc"
requested_group <- cli_arg_value(
  "--group",
  default = Sys.getenv("PROTEOMICS_RANK_ABUNDANCE_GROUP_FILTER", unset = ""),
  args = run$args
)
canonical_group <- if (nzchar(trimws(requested_group))) ca_normalize_group(requested_group)[[1]] else ""
PATHS <- ca_namespace_paths(qc_paths(SUBSTEP_ID, DATASET), canonical_group)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04_marker_rank_abundance_qc.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      paste0("Experimental-group filter: ", if (nzchar(canonical_group)) canonical_group else "all groups (legacy default)"),
      paste0("Output namespace: ", if (nzchar(canonical_group)) paste0("group_", canonical_group) else "existing all-group namespace"),
      "Writes rank tables, marker score tables, and SVG marker sanity-check plots."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggrepel", "scales", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Rank-abundance matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_load_canonical_expression(matrix_file, metadata_file, DATASET)
filtered <- ca_filter_expression_group(expr$mat, expr$meta, requested_group)
mat <- filtered$mat
meta <- filtered$meta
if (nzchar(canonical_group)) {
  invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))
}

marker_sets <- load_wgcna_marker_sets(include_empirical = FALSE)
marker_source_metadata <- attr(marker_sets, "marker_source_metadata")
marker_registry_version <- attr(marker_sets, "marker_registry_version") %||% NA_character_
if ("reference_microglia_pvm" %in% names(marker_sets)) {
  marker_sets$microglia <- marker_sets$reference_microglia_pvm
}
preferred_panels <- c(
  "microglia",
  "reference_microglia_pvm",
  "reference_cortical_excitatory_neuron",
  "reference_hippocampal_excitatory_neuron",
  "reference_inhibitory_interneuron",
  "reference_astrocyte",
  "reference_oligodendrocyte",
  "reference_vascular",
  "canonical_microglia_homeostatic",
  "canonical_microglia_phagolysosomal_state",
  "canonical_neuronal_synaptic_neuropil",
  "canonical_neuronal_soma_nuclear",
  "canonical_astrocyte",
  "canonical_oligodendrocyte_myelin",
  "canonical_endothelial_vascular",
  "mitochondrial_oxphos",
  "ribosomal_translation",
  "rnp_rna_processing"
)
preferred_panels <- preferred_panels[preferred_panels %in% names(marker_sets)]
marker_sets <- marker_sets[unique(c(preferred_panels, names(marker_sets)))]
attr(marker_sets, "marker_source_metadata") <- marker_source_metadata
attr(marker_sets, "marker_registry_version") <- marker_registry_version
if (!length(marker_sets)) stop("No marker panels available for rank-abundance QC.", call. = FALSE)

marker_registry <- dplyr::bind_rows(lapply(names(marker_sets), function(panel) {
  data.frame(marker_panel = panel, marker_gene = marker_sets[[panel]], stringsAsFactors = FALSE)
}))
marker_registry_file <- Sys.getenv(
  "PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE",
  unset = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
)
expr <- qc_add_input_manifest_paths(expr, c(marker_registry = marker_registry_file))
marker_matches <- qc_match_markers_to_protein_groups(
  marker_registry, expr$member_bridge, expr$feature_table,
  panel_col = "marker_panel", gene_col = "marker_gene"
)
qc_write_canonical_feature_artifacts(expr, PATHS$tables, marker_matches)

sample_scores <- lapply(names(marker_sets), function(panel) {
  ids <- marker_matches$matches$ProteinGroupID[marker_matches$matches$marker_panel == panel]
  idx <- match(ids, expr$feature_table$ProteinGroupID, nomatch = 0L)
  idx <- unique(idx[idx > 0L])
  if (!length(idx)) {
    return(data.frame(Sample = colnames(mat), marker_panel = panel, n_markers_detected = 0L,
                      marker_score = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(
    Sample = colnames(mat),
    marker_panel = panel,
    n_markers_detected = length(idx),
    marker_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
sample_scores <- dplyr::bind_rows(sample_scores) |>
  dplyr::left_join(meta, by = "Sample")

rank_group_cols <- intersect(c("Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "plate"), names(meta))
if (!length(rank_group_cols)) rank_group_cols <- "Sample"
rank_group_cols <- unique(c(intersect(c("Region", "region", "Layer", "layer"), rank_group_cols), rank_group_cols[1]))

long <- as.data.frame(mat, check.names = FALSE) |>
  tibble::rownames_to_column("ProteinGroupID") |>
  tidyr::pivot_longer(-ProteinGroupID, names_to = "Sample", values_to = "Log2Intensity") |>
  dplyr::left_join(meta, by = "Sample")

rank_data <- long |>
  dplyr::mutate(RankGroup = do.call(paste, c(dplyr::across(dplyr::all_of(rank_group_cols)), sep = " | "))) |>
  dplyr::filter(!is.na(Log2Intensity), !is.na(RankGroup), nzchar(RankGroup)) |>
  dplyr::group_by(RankGroup, ProteinGroupID) |>
  dplyr::summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(LinearValue = 2^MeanLog2) |>
  dplyr::group_by(RankGroup) |>
  dplyr::arrange(dplyr::desc(LinearValue), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup()

rank_data <- rank_data |>
  dplyr::left_join(
    expr$feature_table |> dplyr::select("ProteinGroupID", "FeatureDisplayLabel", "original_identifier", "member_accessions", "member_gene_symbols"),
    by = "ProteinGroupID"
  ) |>
  dplyr::left_join(
    marker_matches$matches |> dplyr::select("ProteinGroupID", "marker_panel", "matched_official_genes", "matched_member_accessions", "conflicting_marker_panels"),
    by = "ProteinGroupID",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(marker_panel = ifelse(is.na(.data$marker_panel), "none", .data$marker_panel))

qc_write_csv(rank_data, file.path(PATHS$tables, "rank_abundance_table.csv"))
qc_write_csv(sample_scores, file.path(PATHS$tables, "marker_scores_by_sample.csv"))
qc_write_xlsx(list(rank_abundance = rank_data, marker_scores = sample_scores),
              file.path(PATHS$tables, "rank_abundance_marker_qc.xlsx"))

plot_data <- rank_data |>
  dplyr::mutate(is_marker = marker_panel != "none")
label_data <- plot_data |>
  dplyr::filter(is_marker) |>
  dplyr::group_by(RankGroup, marker_panel) |>
  dplyr::slice_min(Rank, n = 3, with_ties = FALSE) |>
  dplyr::ungroup()

p_rank <- ggplot(plot_data, aes(Rank, LinearValue)) +
  geom_point(color = "grey78", alpha = 0.08, size = 0.2) +
  geom_point(data = label_data, aes(color = marker_panel), size = 1.2) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(label = FeatureDisplayLabel, fill = marker_panel),
    color = "white", size = 2, label.size = 0, max.overlaps = 100
  ) +
  scale_y_log10(labels = scales::label_number()) +
  facet_wrap(~RankGroup, scales = "free_x") +
  labs(x = "Protein rank", y = "Intensity (log10)", color = "Marker panel", fill = "Marker panel") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(PATHS$figures, "rank_abundance_marker_sanity_check.svg"), p_rank,
       width = 180, height = 140, units = "mm", device = svglite::svglite)

summary_vars <- intersect(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "plate"), names(sample_scores))
summary_vars <- summary_vars[!duplicated(tolower(summary_vars))]
if (length(summary_vars)) {
  score_summary <- sample_scores |>
    dplyr::mutate(dplyr::across(dplyr::all_of(summary_vars), as.character)) |>
    tidyr::pivot_longer(dplyr::all_of(summary_vars), names_to = "metadata_term", values_to = "metadata_value") |>
    dplyr::filter(!is.na(metadata_value), nzchar(as.character(metadata_value))) |>
    dplyr::group_by(marker_panel, metadata_term, metadata_value) |>
    dplyr::summarise(n = dplyr::n(), mean_score = mean(marker_score, na.rm = TRUE),
                     median_score = median(marker_score, na.rm = TRUE), .groups = "drop")
  qc_write_csv(score_summary, file.path(PATHS$tables, "marker_score_summary_by_metadata.csv"))

  p_score <- ggplot(sample_scores, aes(x = marker_panel, y = marker_score, color = .data[[summary_vars[[1]]]])) +
    geom_boxplot(outlier.shape = NA, color = "grey35") +
    geom_point(position = position_jitter(width = 0.16), alpha = 0.75, size = 1) +
    coord_flip() +
    labs(x = NULL, y = "Mean marker abundance", color = summary_vars[[1]]) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "marker_abundance_summary.svg"), p_score,
         width = 150, height = 100, units = "mm", device = svglite::svglite)
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(
    matrix = matrix_file,
    metadata = metadata_file,
    canonical_input_manifest = file.path(PATHS$tables, "input_path_hash_manifest.csv"),
    marker_registry = marker_registry_file
  ),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables),
  parameters = list(
    dataset = DATASET,
    marker_sets = names(marker_sets),
    marker_source_hierarchy = unique(marker_source_metadata$marker_source %||% "unknown"),
    marker_registry_version = attr(marker_sets, "marker_registry_version") %||% NA_character_,
    requested_group = if (nzchar(filtered$requested_group)) filtered$requested_group else "all",
    canonical_group = if (nzchar(filtered$canonical_group)) filtered$canonical_group else "all",
    resolved_group_column = filtered$resolved_group_column,
    included_samples = filtered$included_samples,
    included_animals = filtered$included_animals,
    allen_microglia_alias = if ("reference_microglia_pvm" %in% names(marker_sets)) "microglia uses reference_microglia_pvm" else NA_character_
  ),
  notes = "Marker abundance/compartment sanity checks only; not cell-type purity estimates. Allen microglia markers are labelled microglia_pvm in the source and aliased to microglia for this QC view when present."
)

message("Rank-abundance and marker QC complete for dataset: ", DATASET)
