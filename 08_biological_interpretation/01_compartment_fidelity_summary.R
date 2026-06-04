# Cross-dataset compartment marker fidelity summary.
#
# Reads marker score tables from 03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r
# and produces compact soma/neuropil/microglia marker plots across the active
# compartment datasets. This is a biological interpretation / poster-summary layer,
# not formal cell-type purity estimation or deconvolution.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "01_compartment_fidelity_summary"
SOURCE_SUBSTEP_ID <- "04c_marker_detectability_and_wgcna_bridge"
PATHS <- create_module_dirs("08_biological_interpretation", file.path(SUBSTEP_ID, DATASET))

compartment_datasets <- strsplit(
  Sys.getenv("PROTEOMICS_COMPARTMENT_MARKER_DATASETS", unset = "neuron_soma,neuron_neuropil,microglia"),
  ","
)[[1]]
compartment_datasets <- trimws(compartment_datasets)
compartment_datasets <- compartment_datasets[nzchar(compartment_datasets)]

score_file_for_dataset <- function(dataset) {
  file.path(qc_paths(SOURCE_SUBSTEP_ID, dataset)$tables, "marker_detectability_by_sample.csv")
}

score_files <- stats::setNames(vapply(compartment_datasets, score_file_for_dataset, character(1)), compartment_datasets)

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "08_biological_interpretation/01_compartment_fidelity_summary.R",
    DATASET,
    matrix_file = NA_character_,
    metadata_file = NA_character_,
    paths = PATHS,
    extra = c(
      paste0("Input score table [", names(score_files), "]: ", score_files, ifelse(file.exists(score_files), " [PASS]", " [MISSING]")),
      "Writes cross-dataset compartment marker boxplots and heatmaps.",
      "Requires marker_detectability_by_sample.csv from 04c for at least two datasets; complete view requires neuron_soma, neuron_neuropil, and microglia."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

read_score_file <- function(dataset) {
  score_file <- score_file_for_dataset(dataset)
  if (!file.exists(score_file)) return(NULL)
  df <- tryCatch(qc_read_table(score_file), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  required_cols <- c("Sample", "marker_panel", "marker_score")
  if (!all(required_cols %in% names(df))) {
    warning("Skipping malformed marker score file for dataset '", dataset, "': ", score_file)
    return(NULL)
  }
  df$dataset <- dataset
  df
}

classify_compartment_marker_panel <- function(panel) {
  panel_l <- tolower(as.character(panel))
  dplyr::case_when(
    grepl("soma|nuclear|histone", panel_l) ~ "Soma markers",
    grepl("neuropil|synaptic|neuronal", panel_l) ~ "Neuropil markers",
    grepl("microglia", panel_l) ~ "Microglia markers",
    TRUE ~ NA_character_
  )
}

raw_scores <- dplyr::bind_rows(lapply(compartment_datasets, read_score_file))
if (is.null(raw_scores) || !nrow(raw_scores)) {
  stop("No 04c marker score tables found. Run 03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r for neuron_soma, neuron_neuropil, and microglia first.", call. = FALSE)
}

compartment_scores <- raw_scores |>
  dplyr::mutate(
    dataset = factor(as.character(dataset), levels = compartment_datasets),
    marker_class = classify_compartment_marker_panel(marker_panel),
    marker_score = as.numeric(marker_score)
  ) |>
  dplyr::filter(!is.na(dataset), !is.na(marker_class), !is.na(marker_score), is.finite(marker_score)) |>
  dplyr::group_by(dataset, Sample, marker_class) |>
  dplyr::summarise(
    marker_score = mean(marker_score, na.rm = TRUE),
    n_marker_panels = dplyr::n_distinct(marker_panel),
    marker_panels = collapse_unique(marker_panel),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    marker_class = factor(marker_class, levels = c("Soma markers", "Neuropil markers", "Microglia markers"))
  )

if (dplyr::n_distinct(compartment_scores$dataset) < 2) {
  stop("Compartment marker fidelity plots require at least two datasets. Found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", "), call. = FALSE)
}
if (dplyr::n_distinct(compartment_scores$marker_class) < 2) {
  stop("Compartment marker fidelity plots require at least two marker classes. Check marker panel names in the 04c output tables.", call. = FALSE)
}

qc_write_csv(compartment_scores, file.path(PATHS$tables, "compartment_marker_fidelity_scores.csv"))

summary_scores <- compartment_scores |>
  dplyr::group_by(marker_class, dataset) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    median_marker_score = median(marker_score, na.rm = TRUE),
    mean_marker_score = mean(marker_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score),
    mean_marker_score = ifelse(is.nan(mean_marker_score), NA_real_, mean_marker_score)
  )
qc_write_csv(summary_scores, file.path(PATHS$tables, "compartment_marker_fidelity_summary.csv"))

p_box <- ggplot(compartment_scores, aes(x = dataset, y = marker_score)) +
  geom_boxplot(outlier.shape = NA, color = "grey35") +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 0.8) +
  facet_wrap(~ marker_class, nrow = 1, scales = "free_y") +
  labs(x = "Dataset", y = "Mean marker abundance") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank())
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_boxplots.svg"), p_box, width = 180, height = 80, units = "mm", device = svglite::svglite)

p_heat <- ggplot(summary_scores, aes(x = dataset, y = marker_class, fill = median_marker_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(x = "Dataset", y = NULL, fill = "Median marker score") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_heatmap.svg"), p_heat, width = 105, height = 70, units = "mm", device = svglite::svglite)

z_scores <- compartment_scores |>
  dplyr::group_by(marker_class) |>
  dplyr::mutate(marker_score_z = as.numeric(scale(marker_score))) |>
  dplyr::ungroup() |>
  dplyr::group_by(marker_class, dataset) |>
  dplyr::summarise(median_marker_score_z = median(marker_score_z, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(median_marker_score_z = ifelse(is.nan(median_marker_score_z), NA_real_, median_marker_score_z))
qc_write_csv(z_scores, file.path(PATHS$tables, "compartment_marker_fidelity_zscore_summary.csv"))

p_z <- ggplot(z_scores, aes(x = dataset, y = marker_class, fill = median_marker_score_z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(x = "Dataset", y = NULL, fill = "Median z-score") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_zscore_heatmap.svg"), p_z, width = 105, height = 70, units = "mm", device = svglite::svglite)

notes <- c(
  "# Compartment marker fidelity summary", "",
  "This report compares soma, neuropil/synaptic, and microglia marker scores across compartment datasets.", "",
  "Interpretation constraints:",
  "- This is a compartment sanity check, not formal purity estimation or deconvolution.",
  "- The microglia marker class averages all available microglia marker panels, for example homeostatic and phagolysosomal panels.",
  "- Expected qualitative pattern: soma markers highest in neuron_soma, neuropil/synaptic markers highest in neuron_neuropil, and microglia markers highest in microglia.",
  "- Deviations can reflect true biological overlap, contamination, marker detectability limits, or missingness/imputation effects.", "",
  paste0("Datasets requested: ", paste(compartment_datasets, collapse = ", ")),
  paste0("Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", ")),
  paste0("Source substep: 03_qc_exploration/", SOURCE_SUBSTEP_ID)
)
writeLines(notes, file.path(PATHS$reports, "compartment_marker_fidelity_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = as.list(score_files),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables, reports = PATHS$reports),
  parameters = list(dataset = DATASET, compartment_datasets = compartment_datasets),
  notes = "Cross-dataset soma/neuropil/microglia marker fidelity summary from 04c marker_detectability_by_sample.csv outputs."
)

message("Compartment marker fidelity summary complete. Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", "))
