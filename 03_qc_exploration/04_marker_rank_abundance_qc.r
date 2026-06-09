# ================================================================
# Script: 03_qc_exploration/04_marker_rank_abundance_qc.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional config/marker_panels/wgcna_reference_marker_sets.csv; results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv.
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

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04_marker_rank_abundance_qc"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04_marker_rank_abundance_qc.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c("Writes rank tables, marker score tables, and SVG marker sanity-check plots.")
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

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- expr$meta

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

gene_norm <- normalize_gene_token
protein_ids <- rownames(mat)
protein_key <- gene_norm(sub("_MOUSE$", "", protein_ids, ignore.case = TRUE))

sample_scores <- lapply(names(marker_sets), function(panel) {
  genes <- marker_sets[[panel]]
  idx <- which(protein_key %in% gene_norm(genes))
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
  tibble::rownames_to_column("Protein") |>
  tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Log2Intensity") |>
  dplyr::left_join(meta, by = "Sample")

rank_data <- long |>
  dplyr::mutate(RankGroup = do.call(paste, c(dplyr::across(dplyr::all_of(rank_group_cols)), sep = " | "))) |>
  dplyr::filter(!is.na(Log2Intensity), !is.na(RankGroup), nzchar(RankGroup)) |>
  dplyr::group_by(RankGroup, Protein) |>
  dplyr::summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(LinearValue = 2^MeanLog2) |>
  dplyr::group_by(RankGroup) |>
  dplyr::arrange(dplyr::desc(LinearValue), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup()

marker_lookup <- dplyr::bind_rows(lapply(names(marker_sets), function(panel) {
  data.frame(marker_panel = panel, marker = marker_sets[[panel]], marker_key = gene_norm(marker_sets[[panel]]))
})) |>
  dplyr::filter(nzchar(.data$marker_key)) |>
  dplyr::distinct(.data$marker_panel, .data$marker_key, .keep_all = TRUE)
rank_data <- rank_data |>
  dplyr::mutate(marker_key = gene_norm(sub("_MOUSE$", "", Protein, ignore.case = TRUE))) |>
  dplyr::left_join(marker_lookup, by = "marker_key", relationship = "many-to-many") |>
  dplyr::mutate(marker_panel = ifelse(is.na(marker_panel), "none", marker_panel))

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
    aes(label = Protein, fill = marker_panel),
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
    marker_registry = Sys.getenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE", unset = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"))
  ),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables),
  parameters = list(
    dataset = DATASET,
    marker_sets = names(marker_sets),
    marker_source_hierarchy = unique(marker_source_metadata$marker_source %||% "unknown"),
    marker_registry_version = attr(marker_sets, "marker_registry_version") %||% NA_character_,
    allen_microglia_alias = if ("reference_microglia_pvm" %in% names(marker_sets)) "microglia uses reference_microglia_pvm" else NA_character_
  ),
  notes = "Marker abundance/compartment sanity checks only; not cell-type purity estimates. Allen microglia markers are labelled microglia_pvm in the source and aliased to microglia for this QC view when present."
)

message("Rank-abundance and marker QC complete for dataset: ", DATASET)
