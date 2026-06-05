# Cross-dataset compartment marker fidelity QC.
#
# Reads explicit soma/neuropil/microglia marker-fidelity score tables from
# 04c_marker_detectability_and_wgcna_bridge.r and produces compact plots across
# the three active compartment datasets. This is QC/interpreter support, not
# formal purity estimation or deconvolution.
#
# Preferred input from 04c:
#   marker_fidelity_by_sample.csv
#
# Backward-compatible fallback input from older 04c versions:
#   marker_detectability_by_sample.csv
#
# The preferred input is intentionally narrow and uses only the explicitly
# curated fidelity marker classes:
#   Soma markers
#   Neuropil markers
#   Microglia/PVM markers
#
# Broad Allen cell-class reference panels remain useful in 04c broad QC outputs,
# but are not used here to define soma versus neuropil marker fidelity.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04d_compartment_marker_fidelity"
SOURCE_SUBSTEP_ID <- "04c_marker_detectability_and_wgcna_bridge"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)

compartment_datasets <- strsplit(
  Sys.getenv("PROTEOMICS_COMPARTMENT_MARKER_DATASETS", unset = "neuron_soma,neuron_neuropil,microglia"),
  ","
)[[1]]
compartment_datasets <- trimws(compartment_datasets)
compartment_datasets <- compartment_datasets[nzchar(compartment_datasets)]

fidelity_class_order <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")

fidelity_file_for_dataset <- function(dataset) {
  file.path(qc_paths(SOURCE_SUBSTEP_ID, dataset)$tables, "marker_fidelity_by_sample.csv")
}

legacy_score_file_for_dataset <- function(dataset) {
  file.path(qc_paths(SOURCE_SUBSTEP_ID, dataset)$tables, "marker_detectability_by_sample.csv")
}

fidelity_files <- stats::setNames(vapply(compartment_datasets, fidelity_file_for_dataset, character(1)), compartment_datasets)
legacy_score_files <- stats::setNames(vapply(compartment_datasets, legacy_score_file_for_dataset, character(1)), compartment_datasets)

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04d_compartment_marker_fidelity.r",
    DATASET,
    matrix_file = NA_character_,
    metadata_file = NA_character_,
    paths = PATHS,
    extra = c(
      paste0("Preferred fidelity table [", names(fidelity_files), "]: ", fidelity_files, ifelse(file.exists(fidelity_files), " [PASS]", " [MISSING]")),
      paste0("Fallback legacy score table [", names(legacy_score_files), "]: ", legacy_score_files, ifelse(file.exists(legacy_score_files), " [PASS]", " [MISSING]")),
      "Writes cross-dataset soma/neuropil/microglia marker fidelity boxplots and heatmaps.",
      "Preferred input is marker_fidelity_by_sample.csv from 04c; fallback uses marker_detectability_by_sample.csv with conservative panel-name matching."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

classify_legacy_compartment_marker_panel <- function(panel) {
  panel_l <- tolower(as.character(panel))
  dplyr::case_when(
    panel_l %in% c("canonical_neuronal_soma_nuclear", "nuclear_soma") ~ "Soma markers",
    panel_l %in% c("canonical_neuronal_synaptic_neuropil", "neuronal_synaptic_neuropil") ~ "Neuropil markers",
    panel_l %in% c(
      "canonical_microglia_homeostatic",
      "canonical_microglia_phagolysosomal_state",
      "reference_microglia_pvm",
      "microglia_homeostatic",
      "microglia_phagolysosomal"
    ) ~ "Microglia/PVM markers",
    TRUE ~ NA_character_
  )
}

read_fidelity_score_file <- function(dataset) {
  fidelity_file <- fidelity_file_for_dataset(dataset)
  if (!file.exists(fidelity_file)) return(NULL)

  df <- tryCatch(qc_read_table(fidelity_file), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)

  required_cols <- c("Sample", "fidelity_marker_class", "fidelity_marker_score")
  if (!all(required_cols %in% names(df))) {
    warning("Skipping malformed marker fidelity file for dataset '", dataset, "': ", fidelity_file)
    return(NULL)
  }

  out <- df |>
    dplyr::transmute(
      dataset = dataset,
      Sample = as.character(.data$Sample),
      marker_class = as.character(.data$fidelity_marker_class),
      marker_score = as.numeric(.data$fidelity_marker_score),
      n_detected_markers = if ("n_detected_markers" %in% names(df)) as.integer(.data$n_detected_markers) else NA_integer_,
      fraction_detected_markers = if ("fraction_detected_markers" %in% names(df)) as.numeric(.data$fraction_detected_markers) else NA_real_,
      source_table = "marker_fidelity_by_sample.csv"
    )

  out
}

read_legacy_score_file <- function(dataset) {
  score_file <- legacy_score_file_for_dataset(dataset)
  if (!file.exists(score_file)) return(NULL)

  df <- tryCatch(qc_read_table(score_file), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)

  required_cols <- c("Sample", "marker_panel", "marker_score")
  if (!all(required_cols %in% names(df))) {
    warning("Skipping malformed legacy marker score file for dataset '", dataset, "': ", score_file)
    return(NULL)
  }

  out <- df |>
    dplyr::mutate(
      dataset = dataset,
      marker_class = classify_legacy_compartment_marker_panel(.data$marker_panel),
      marker_score = as.numeric(.data$marker_score)
    ) |>
    dplyr::filter(!is.na(marker_class), is.finite(marker_score)) |>
    dplyr::group_by(dataset, Sample, marker_class) |>
    dplyr::summarise(
      marker_score = mean(marker_score, na.rm = TRUE),
      n_detected_markers = if ("n_detected_markers" %in% names(df)) sum(unique(n_detected_markers), na.rm = TRUE) else NA_integer_,
      fraction_detected_markers = if ("fraction_detected_markers" %in% names(df)) mean(fraction_detected_markers, na.rm = TRUE) else NA_real_,
      legacy_marker_panels = collapse_unique(marker_panel),
      .groups = "drop"
    ) |>
    dplyr::mutate(source_table = "marker_detectability_by_sample.csv")

  out
}

read_score_file <- function(dataset) {
  preferred <- read_fidelity_score_file(dataset)
  if (!is.null(preferred) && nrow(preferred)) return(preferred)

  fallback <- read_legacy_score_file(dataset)
  if (!is.null(fallback) && nrow(fallback)) {
    warning("Using legacy marker_detectability_by_sample.csv for dataset '", dataset,
            "'. Re-run 04c to generate marker_fidelity_by_sample.csv for the preferred explicit fidelity input.")
    return(fallback)
  }

  NULL
}

raw_scores <- dplyr::bind_rows(lapply(compartment_datasets, read_score_file))
if (is.null(raw_scores) || !nrow(raw_scores)) {
  stop("No 04c marker fidelity score tables found. Run 03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r for neuron_soma, neuron_neuropil, and microglia first.", call. = FALSE)
}

compartment_scores <- raw_scores |>
  dplyr::mutate(
    dataset = factor(as.character(dataset), levels = compartment_datasets),
    marker_class = factor(as.character(marker_class), levels = fidelity_class_order),
    marker_score = as.numeric(marker_score)
  ) |>
  dplyr::filter(!is.na(dataset), !is.na(marker_class), !is.na(marker_score), is.finite(marker_score)) |>
  dplyr::group_by(dataset, Sample, marker_class, source_table) |>
  dplyr::summarise(
    marker_score = mean(marker_score, na.rm = TRUE),
    n_detected_markers = if (all(is.na(n_detected_markers))) NA_integer_ else max(n_detected_markers, na.rm = TRUE),
    fraction_detected_markers = if (all(is.na(fraction_detected_markers))) NA_real_ else max(fraction_detected_markers, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(marker_class, dataset, Sample)

if (dplyr::n_distinct(compartment_scores$dataset) < 2) {
  stop("Compartment marker fidelity plots require at least two datasets. Found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", "), call. = FALSE)
}
if (dplyr::n_distinct(compartment_scores$marker_class) < 2) {
  stop("Compartment marker fidelity plots require at least two marker classes. Check marker_fidelity_by_sample.csv from 04c.", call. = FALSE)
}

qc_write_csv(compartment_scores, file.path(PATHS$tables, "compartment_marker_fidelity_scores.csv"))

summary_scores <- compartment_scores |>
  dplyr::group_by(marker_class, dataset) |>
  dplyr::summarise(
    n_samples = dplyr::n_distinct(Sample),
    median_marker_score = median(marker_score, na.rm = TRUE),
    mean_marker_score = mean(marker_score, na.rm = TRUE),
    median_detected_markers = median(n_detected_markers, na.rm = TRUE),
    median_fraction_detected_markers = median(fraction_detected_markers, na.rm = TRUE),
    source_tables = collapse_unique(source_table),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score),
    mean_marker_score = ifelse(is.nan(mean_marker_score), NA_real_, mean_marker_score),
    median_detected_markers = ifelse(is.nan(median_detected_markers), NA_real_, median_detected_markers),
    median_fraction_detected_markers = ifelse(is.nan(median_fraction_detected_markers), NA_real_, median_fraction_detected_markers)
  ) |>
  dplyr::arrange(marker_class, dataset)
qc_write_csv(summary_scores, file.path(PATHS$tables, "compartment_marker_fidelity_summary.csv"))

p_box <- ggplot(compartment_scores, aes(x = dataset, y = marker_score)) +
  geom_boxplot(outlier.shape = NA, color = "grey35") +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 0.8) +
  facet_wrap(~ marker_class, nrow = 1, scales = "free_y") +
  labs(x = "Dataset", y = "Mean fidelity marker abundance") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank())
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_boxplots.svg"), p_box, width = 180, height = 80, units = "mm", device = svglite::svglite)

p_heat <- ggplot(summary_scores, aes(x = dataset, y = marker_class, fill = median_marker_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(x = "Dataset", y = NULL, fill = "Median marker score") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_heatmap.svg"), p_heat, width = 115, height = 70, units = "mm", device = svglite::svglite)

z_scores <- compartment_scores |>
  dplyr::group_by(marker_class) |>
  dplyr::mutate(
    marker_score_z = {
      x <- marker_score
      sx <- stats::sd(x, na.rm = TRUE)
      mx <- mean(x, na.rm = TRUE)
      if (!is.finite(sx) || sx <= .Machine$double.eps) NA_real_ else (x - mx) / sx
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(marker_class, dataset) |>
  dplyr::summarise(median_marker_score_z = median(marker_score_z, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(median_marker_score_z = ifelse(is.nan(median_marker_score_z), NA_real_, median_marker_score_z)) |>
  dplyr::arrange(marker_class, dataset)
qc_write_csv(z_scores, file.path(PATHS$tables, "compartment_marker_fidelity_zscore_summary.csv"))

p_z <- ggplot(z_scores, aes(x = dataset, y = marker_class, fill = median_marker_score_z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(x = "Dataset", y = NULL, fill = "Median z-score") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "compartment_marker_fidelity_zscore_heatmap.svg"), p_z, width = 115, height = 70, units = "mm", device = svglite::svglite)

notes <- c(
  "# Compartment marker fidelity QC", "",
  "This report compares explicit soma, neuropil/synaptic, and microglia/PVM fidelity marker scores across compartment datasets.", "",
  "Interpretation constraints:",
  "- This is a compartment sanity check, not formal purity estimation or deconvolution.",
  "- Preferred input is marker_fidelity_by_sample.csv from 04c_marker_detectability_and_wgcna_bridge.r.",
  "- The fidelity classes are intentionally narrow: Soma markers, Neuropil markers, and Microglia/PVM markers.",
  "- Broad Allen neuronal cell-class panels are not used to define soma versus neuropil; they remain in the broader 04c marker-panel/marker-compartment QC outputs.",
  "- Expected qualitative pattern: soma markers highest in neuron_soma, neuropil/synaptic markers highest in neuron_neuropil, and microglia/PVM markers highest in microglia.",
  "- Deviations can reflect true biological overlap, contamination, marker detectability limits, or missingness/imputation effects.", "",
  paste0("Datasets requested: ", paste(compartment_datasets, collapse = ", ")),
  paste0("Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", ")),
  paste0("Marker classes found: ", paste(unique(as.character(compartment_scores$marker_class)), collapse = ", ")),
  paste0("Source table(s): ", paste(unique(compartment_scores$source_table), collapse = ", ")),
  paste0("Source substep: ", SOURCE_SUBSTEP_ID)
)
writeLines(notes, file.path(PATHS$reports, "compartment_marker_fidelity_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(preferred_fidelity_files = as.list(fidelity_files), fallback_legacy_score_files = as.list(legacy_score_files)),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables, reports = PATHS$reports),
  parameters = list(dataset = DATASET, compartment_datasets = compartment_datasets, marker_classes = fidelity_class_order),
  notes = "Cross-dataset soma/neuropil/microglia-PVM marker fidelity QC from 04c marker_fidelity_by_sample.csv outputs, with legacy fallback."
)

message("Compartment marker fidelity QC complete. Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", "))
