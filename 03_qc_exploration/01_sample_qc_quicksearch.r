# ================================================================
# Script: 03_qc_exploration/01_sample_qc_quicksearch.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/raw/pg_matrix/quicksearch.stats.annotated.xlsx; optional data/metadata/*.xlsx.
# Produces: results/figures/03_qc_exploration/01_sample_qc_quicksearch/<dataset>/; results/tables/03_qc_exploration/01_sample_qc_quicksearch/<dataset>/qc_summary_tables.xlsx.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Quicksearch sample-level QC.
# ================================================================

# ================================================================
# QC figures for spatial proteomics
# This script generates QC figures for the main manuscript and supplemental materials, as well as summary tables. # nolint
# Consumes an annotated QC stats file (e.g. from quicksearch) and produces consistent visualizations of key QC metrics across samples and cell types.
# The main figure focuses on core QC metrics like protein/precursor counts, signal, mass accuracy, and a composite QC score, while the supplemental figure includes additional metrics and missingness estimates.
# Outlier detection is performed using a robust z-score method, but outliers are not highlighted in the main figure to maintain visual clarity. Instead, a separate outlier summary figure and table are provided.
# The script is designed to be flexible to varying QC metrics and sample annotations, and includes error handling for missing data. It also saves all outputs in a specified results directory.
# Author: Tobias Pohl
# ================================================================

early_paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(early_paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))
early_run <- qc_args()
if (early_run$dry_run) {
  early_dataset <- early_run$dataset
  early_paths <- qc_paths("01_sample_qc_quicksearch", early_dataset)
  early_input <- Sys.getenv(
    "PROTEOMICS_QC_STATS_FILE",
    unset = path_raw("pg_matrix", "quicksearch.stats.annotated.xlsx")
  )
  early_out_dir <- Sys.getenv("PROTEOMICS_QC_PUBLICATION_DIR", unset = early_paths$figures)
  early_table_dir <- Sys.getenv("PROTEOMICS_QC_TABLE_DIR", unset = early_paths$tables)
  status <- qc_dry_run_contract(
    "03_qc_exploration/01_sample_qc_quicksearch.r",
    early_dataset,
    matrix_file = early_input,
    paths = list(figures = early_out_dir, tables = early_table_dir, logs = early_paths$logs),
    extra = c("Uses PROTEOMICS_QC_STATS_FILE for quicksearch stats override.")
  )
  quit(status = status, save = "no")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(scales)
  library(openxlsx)
  library(rlang)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
MODULE_ID <- "03_qc_exploration"
SUBSTEP_ID <- "01_sample_qc_quicksearch"
CANONICAL_PATHS <- qc_paths(SUBSTEP_ID, DATASET)

# ================================================================
# 1. Paths
# ================================================================

input_file <- Sys.getenv(
  "PROTEOMICS_QC_STATS_FILE",
  unset = path_raw("pg_matrix", "quicksearch.stats.annotated.xlsx")
)

out_dir <- Sys.getenv("PROTEOMICS_QC_PUBLICATION_DIR", unset = CANONICAL_PATHS$figures)
global_out_dir <- path_results("figures", MODULE_ID, SUBSTEP_ID, "global")
table_dir <- Sys.getenv("PROTEOMICS_QC_TABLE_DIR", unset = CANONICAL_PATHS$tables)
if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/01_sample_qc_quicksearch.r",
    DATASET,
    matrix_file = input_file,
    paths = list(figures = out_dir, tables = table_dir, logs = CANONICAL_PATHS$logs),
    extra = c("Uses PROTEOMICS_QC_STATS_FILE for quicksearch stats override.")
  )
  quit(status = status, save = "no")
}
if (!file.exists(input_file)) {
  stop(
    "QC stats input not found for dataset '", DATASET, "': ", input_file,
    ". Set PROTEOMICS_QC_STATS_FILE to an annotated quicksearch stats workbook.",
    call. = FALSE
  )
}
ensure_dir(out_dir)
ensure_dir(global_out_dir)
ensure_dir(table_dir)

# ================================================================
# 2. Load and clean data
# ================================================================

qc <- read_excel(input_file) %>%
  mutate(
    sample_id = as.character(sample_id),
    celltype_layer = as.character(celltype_layer)
  ) %>%
  filter(
    !grepl("background|blank|bg", sample_id, ignore.case = TRUE),
    !grepl("background|blank|bg", celltype_layer, ignore.case = TRUE)
  ) %>%
  mutate(celltype_layer = factor(celltype_layer))

add_qc_compartment <- function(df) {
  candidate_cols <- intersect(c("celltype_layer", "CellTypeLayer", "celltype", "CellType"), names(df))
  df$qc_compartment <- NA_character_

  for (col in candidate_cols) {
    vals <- trimws(as.character(df[[col]]))
    vals[!nzchar(vals)] <- NA_character_
    compartment_key <- tolower(vals)
    compartment_key <- gsub("[ -]+", "_", compartment_key)

    mapped <- case_when(
      compartment_key %in% c("microglia", "microglial") ~ "microglia",
      compartment_key %in% c("neuron_neuropil", "neuropil") ~ "neuropil",
      compartment_key %in% c("neuron_soma", "soma") ~ "soma",
      TRUE ~ NA_character_
    )

    fill <- is.na(df$qc_compartment) & !is.na(mapped)
    df$qc_compartment[fill] <- mapped[fill]
  }

  unmatched <- sort(unique(as.character(df$celltype_layer[is.na(df$qc_compartment)])))
  unmatched <- unmatched[!is.na(unmatched) & nzchar(unmatched)]
  if (length(unmatched) > 0) {
    warning(
      "Unmatched non-background celltype_layer values for global QC compartments: ",
      paste(unmatched, collapse = ", "),
      call. = FALSE
    )
  }

  df$qc_compartment <- factor(df$qc_compartment, levels = c("microglia", "neuropil", "soma"))
  df
}

qc_global <- add_qc_compartment(qc)

if ("celltype_layer" %in% names(qc)) {
  keep_dataset <- metadata_matches_dataset(qc, DATASET)
  if (any(keep_dataset)) {
    qc <- qc[keep_dataset, , drop = FALSE]
    qc$celltype_layer <- droplevels(factor(qc$celltype_layer))
  }
}

required_cols <- c("sample_id", "celltype_layer")
missing_required <- setdiff(required_cols, names(qc))

if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

# ================================================================
# 3. QC metric columns
# ================================================================

qc_metrics_all <- c(
  "Proteins.Identified",
  "Precursors.Identified",
  "MS1.Signal",
  "MS2.Signal",
  "FWHM.Scans",
  "FWHM.RT",
  "Median.Mass.Acc.MS1",
  "Median.Mass.Acc.MS2",
  "Median.Mass.Acc.MS1.Corrected",
  "Median.Mass.Acc.MS2.Corrected",
  "Normalisation.Instability",
  "Median.RT.Prediction.Acc",
  "Average.Peptide.Length",
  "Average.Peptide.Charge",
  "Average.Missed.Tryptic.Cleavages"
)

qc_metrics <- intersect(qc_metrics_all, names(qc))

qc <- qc %>%
  mutate(across(all_of(qc_metrics), as.numeric))

qc_global <- qc_global %>%
  mutate(across(all_of(intersect(qc_metrics_all, names(qc_global))), as.numeric))

has_cols <- function(x) all(x %in% names(qc))

# ================================================================
# 4. Publication-style palette
# ================================================================

celltype_cols <- c(
  "microglia" = "#4C78A8",
  "neuron_soma" = "#E45756",
  "neuron_neuropil" = "#72B7B2"
)

missing_levels <- setdiff(levels(qc$celltype_layer), names(celltype_cols))

if (length(missing_levels) > 0) {
  extra_cols <- hue_pal(l = 55, c = 70)(length(missing_levels))
  names(extra_cols) <- missing_levels
  celltype_cols <- c(celltype_cols, extra_cols)
}

neutral_cols <- c(
  dark = "#333333",
  mid = "#777777",
  light = "#BDBDBD",
  faint = "#E6E6E6"
)

# ================================================================
# 5. Theme and save helpers
# ================================================================

theme_publication_qc <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1, colour = "black"),
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(3, "mm"),
      legend.spacing.x = unit(2, "mm"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = neutral_cols["mid"], hjust = 0),
      plot.margin = margin(3, 3, 3, 3),
      panel.spacing = unit(2, "mm")
    )
}

save_svg <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    device = svglite
  )
}

save_global_svg <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(global_out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    device = svglite
  )
}

# ================================================================
# 6. Robust outlier detection
# Kept for tables, not drawn as rings in the main figure
# ================================================================

robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)

  if (is.na(mad_val) || mad_val == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - med) / mad_val
}

outlier_metrics <- intersect(
  c(
    "Proteins.Identified",
    "Precursors.Identified",
    "MS1.Signal",
    "MS2.Signal",
    "Normalisation.Instability",
    "Median.RT.Prediction.Acc",
    "Median.Mass.Acc.MS1",
    "Median.Mass.Acc.MS2"
  ),
  names(qc)
)

qc_outlier_z <- qc %>%
  group_by(celltype_layer) %>%
  mutate(
    across(
      all_of(outlier_metrics),
      robust_z,
      .names = "{.col}_robust_z"
    )
  ) %>%
  ungroup()

z_cols <- grep("_robust_z$", names(qc_outlier_z), value = TRUE)

qc_outlier_z <- qc_outlier_z %>%
  rowwise() %>%
  mutate(
    n_outlier_metrics = sum(abs(c_across(all_of(z_cols))) > 3, na.rm = TRUE),
    max_abs_robust_z = suppressWarnings(max(abs(c_across(all_of(z_cols))), na.rm = TRUE)),
    qc_outlier = n_outlier_metrics >= 2
  ) %>%
  ungroup() %>%
  mutate(
    max_abs_robust_z = ifelse(is.infinite(max_abs_robust_z), NA_real_, max_abs_robust_z)
  )

qc <- qc %>%
  left_join(
    qc_outlier_z %>%
      select(sample_id, n_outlier_metrics, max_abs_robust_z, qc_outlier),
    by = "sample_id"
  )

# ================================================================
# 7. Composite QC score
# Higher = better
# ================================================================

safe_scale <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

score_components <- list()

if ("Proteins.Identified" %in% names(qc)) {
  score_components$protein_depth <- safe_scale(qc$Proteins.Identified)
}

if ("Precursors.Identified" %in% names(qc)) {
  score_components$precursor_depth <- safe_scale(qc$Precursors.Identified)
}

if ("Normalisation.Instability" %in% names(qc)) {
  score_components$norm_stability <- -safe_scale(qc$Normalisation.Instability)
}

if ("Median.RT.Prediction.Acc" %in% names(qc)) {
  score_components$rt_accuracy <- -abs(safe_scale(qc$Median.RT.Prediction.Acc))
}

if ("Median.Mass.Acc.MS1" %in% names(qc)) {
  score_components$mass_ms1 <- -abs(safe_scale(qc$Median.Mass.Acc.MS1))
}

if ("Median.Mass.Acc.MS2" %in% names(qc)) {
  score_components$mass_ms2 <- -abs(safe_scale(qc$Median.Mass.Acc.MS2))
}

if (length(score_components) > 0) {
  score_mat <- as.data.frame(score_components)
  qc$QC.Score <- rowMeans(score_mat, na.rm = TRUE)
} else {
  qc$QC.Score <- NA_real_
}

# ================================================================
# 8. Plot helper functions
# No ring/outlier overlays
# ================================================================

plot_sample_bars <- function(df, yvar, ylab) {
  df %>%
    arrange(celltype_layer, .data[[yvar]]) %>%
    mutate(sample_order = factor(sample_id, levels = unique(sample_id))) %>%
    ggplot(aes(x = sample_order, y = .data[[yvar]], fill = celltype_layer)) +
    geom_col(width = 0.75, colour = NA) +
    facet_grid(. ~ celltype_layer, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = celltype_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = ylab) +
    theme_publication_qc() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.spacing.x = unit(1.5, "mm")
    )
}

global_compartment_cols <- c(
  "microglia" = "#a7d8d5",
  "neuropil" = "#156959",
  "soma" = "#757476"
)

plot_global_sample_bars <- function(df, yvar, ylab) {
  df %>%
    filter(!is.na(qc_compartment), !is.na(.data[[yvar]])) %>%
    arrange(qc_compartment, .data[[yvar]]) %>%
    mutate(
      sample_key = paste(qc_compartment, sample_id, row_number(), sep = "__"),
      sample_order = factor(sample_key, levels = sample_key)
    ) %>%
    ggplot(aes(x = sample_order, y = .data[[yvar]], fill = qc_compartment)) +
    geom_col(width = 1, colour = NA) +
    facet_grid(. ~ qc_compartment, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = global_compartment_cols) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = ylab) +
    theme_publication_qc() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.spacing.x = unit(1.5, "mm")
    )
}

plot_scatter_qc <- function(df, xvar, yvar, xlab, ylab, log_axes = FALSE) {
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]], colour = celltype_layer)) +
    geom_point(size = 1.4, alpha = 0.85) +
    scale_colour_manual(values = celltype_cols) +
    labs(x = xlab, y = ylab) +
    theme_publication_qc() +
    theme(aspect.ratio = 1)

  if (log_axes) {
    p <- p +
      scale_x_log10(labels = label_scientific()) +
      scale_y_log10(labels = label_scientific())
  }

  p
}

plot_box_jitter <- function(df, yvar, ylab) {
  ggplot(df, aes(x = celltype_layer, y = .data[[yvar]], fill = celltype_layer)) +
    geom_boxplot(
      width = 0.52,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.75,
      colour = "black"
    ) +
    geom_point(
      aes(colour = celltype_layer),
      position = position_jitter(width = 0.12, height = 0),
      size = 1,
      alpha = 0.65
    ) +
    scale_fill_manual(values = celltype_cols) +
    scale_colour_manual(values = celltype_cols, guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_publication_qc() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
}

# ================================================================
# 9. Main QC panels
# ================================================================

p_protein <- NULL
if ("Proteins.Identified" %in% names(qc)) {
  p_protein <- plot_sample_bars(qc, "Proteins.Identified", "Proteins identified")
}

p_precursor <- NULL
if ("Precursors.Identified" %in% names(qc)) {
  p_precursor <- plot_sample_bars(qc, "Precursors.Identified", "Precursors identified")
}

p_global_protein <- NULL
if ("Proteins.Identified" %in% names(qc_global)) {
  p_global_protein <- plot_global_sample_bars(qc_global, "Proteins.Identified", "Proteins identified")
}

p_global_precursor <- NULL
if ("Precursors.Identified" %in% names(qc_global)) {
  p_global_precursor <- plot_global_sample_bars(qc_global, "Precursors.Identified", "Precursors identified")
}

p_signal <- NULL
if (has_cols(c("MS1.Signal", "MS2.Signal"))) {
  p_signal <- plot_scatter_qc(
    qc,
    "MS1.Signal",
    "MS2.Signal",
    "MS1 signal",
    "MS2 signal",
    log_axes = TRUE
  )
}

p_mass <- NULL
if (has_cols(c("Median.Mass.Acc.MS1", "Median.Mass.Acc.MS2"))) {
  p_mass <- plot_scatter_qc(
    qc,
    "Median.Mass.Acc.MS1",
    "Median.Mass.Acc.MS2",
    "Median mass accuracy MS1 (ppm)",
    "Median mass accuracy MS2 (ppm)"
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = neutral_cols["mid"]) +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = neutral_cols["mid"])
}

p_norm <- NULL
if ("Normalisation.Instability" %in% names(qc)) {
  p_norm <- plot_box_jitter(qc, "Normalisation.Instability", "Normalisation instability")
}

p_rt <- NULL
if ("Median.RT.Prediction.Acc" %in% names(qc)) {
  p_rt <- plot_box_jitter(qc, "Median.RT.Prediction.Acc", "RT prediction accuracy")
}

p_qc_score <- NULL
if ("QC.Score" %in% names(qc)) {
  p_qc_score <- qc %>%
    arrange(QC.Score) %>%
    mutate(sample_order = factor(sample_id, levels = sample_id)) %>%
    ggplot(aes(x = sample_order, y = QC.Score, fill = celltype_layer)) +
    geom_col(width = 0.75, colour = NA) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = neutral_cols["mid"]) +
    scale_fill_manual(values = celltype_cols) +
    labs(x = NULL, y = "Composite QC score") +
    theme_publication_qc() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    )
}

# ================================================================
# 10. PCA of QC metrics
# ================================================================

p_pca <- NULL

pca_metrics <- intersect(
  c(
    "Proteins.Identified",
    "Precursors.Identified",
    "MS1.Signal",
    "MS2.Signal",
    "FWHM.Scans",
    "FWHM.RT",
    "Median.Mass.Acc.MS1",
    "Median.Mass.Acc.MS2",
    "Normalisation.Instability",
    "Median.RT.Prediction.Acc"
  ),
  names(qc)
)

if (length(pca_metrics) >= 3) {

  pca_df <- qc %>%
    select(sample_id, celltype_layer, qc_outlier, all_of(pca_metrics)) %>%
    drop_na(all_of(pca_metrics))

  if (nrow(pca_df) >= 4) {

    pca_mat <- pca_df %>%
      select(all_of(pca_metrics)) %>%
      scale()

    pca_res <- prcomp(pca_mat, center = TRUE, scale. = FALSE)

    pca_scores <- as.data.frame(pca_res$x[, 1:2]) %>%
      bind_cols(pca_df %>% select(sample_id, celltype_layer, qc_outlier))

    pca_centroids <- pca_scores %>%
      group_by(celltype_layer) %>%
      summarise(
        PC1 = mean(PC1, na.rm = TRUE),
        PC2 = mean(PC2, na.rm = TRUE),
        .groups = "drop"
      )

    pve <- 100 * summary(pca_res)$importance[2, 1:2]

    p_pca <- ggplot(pca_scores, aes(PC1, PC2, colour = celltype_layer)) +
      geom_point(size = 1.5, alpha = 0.9) +
      geom_point(
        data = pca_centroids,
        aes(PC1, PC2, fill = celltype_layer),
        shape = 23,
        size = 2.2,
        stroke = 0.3,
        colour = "black",
        inherit.aes = FALSE
      ) +
      scale_colour_manual(values = celltype_cols) +
      scale_fill_manual(values = celltype_cols, guide = "none") +
      labs(
        x = sprintf("PC1 (%.1f%%)", pve[1]),
        y = sprintf("PC2 (%.1f%%)", pve[2])
      ) +
      theme_publication_qc() +
      theme(aspect.ratio = 1)
  }
}

# ================================================================
# 11. Depth-quality relationship
# ================================================================

p_depth_quality <- NULL

if (has_cols(c("Proteins.Identified", "Normalisation.Instability"))) {
  p_depth_quality <- plot_scatter_qc(
    qc,
    "Proteins.Identified",
    "Normalisation.Instability",
    "Proteins identified",
    "Normalisation instability"
  ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.35,
      colour = neutral_cols["dark"],
      inherit.aes = FALSE,
      aes(x = Proteins.Identified, y = Normalisation.Instability)
    )
}

# ================================================================
# 12. Missingness estimate across QC metric columns
# ================================================================

p_missing <- NULL

if (length(qc_metrics) >= 3) {

  qc_missing <- qc %>%
    mutate(
      missing_fraction_qc_metrics = rowMeans(is.na(select(., all_of(qc_metrics))))
    )

  if ("MS1.Signal" %in% names(qc_missing)) {
    p_missing <- ggplot(
      qc_missing,
      aes(x = MS1.Signal, y = missing_fraction_qc_metrics, colour = celltype_layer)
    ) +
      geom_point(size = 1.4, alpha = 0.85) +
      scale_x_log10(labels = label_scientific()) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_colour_manual(values = celltype_cols) +
      labs(
        x = "MS1 signal",
        y = "Missing QC metric fraction"
      ) +
      theme_publication_qc() +
      theme(aspect.ratio = 1)
  } else {
    p_missing <- ggplot(
      qc_missing,
      aes(x = sample_id, y = missing_fraction_qc_metrics, fill = celltype_layer)
    ) +
      geom_col(width = 0.75, colour = NA) +
      scale_fill_manual(values = celltype_cols) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        x = NULL,
        y = "Missing QC metric fraction"
      ) +
      theme_publication_qc() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
}

# ================================================================
# 13. Coefficient of variation summary
# ================================================================

cv_metrics <- intersect(
  c(
    "Proteins.Identified",
    "Precursors.Identified",
    "MS1.Signal",
    "MS2.Signal"
  ),
  names(qc)
)

cv_summary <- NULL
p_cv <- NULL

if (length(cv_metrics) > 0) {

  cv_summary <- qc %>%
    group_by(celltype_layer) %>%
    summarise(
      across(
        all_of(cv_metrics),
        ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(cv_metrics),
      names_to = "metric",
      values_to = "cv"
    )

  p_cv <- ggplot(cv_summary, aes(x = metric, y = cv, fill = celltype_layer)) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = NA
    ) +
    scale_fill_manual(values = celltype_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Coefficient of variation") +
    theme_publication_qc() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

# ================================================================
# 14. Main manuscript QC figure
# ================================================================

main_plots <- list(
  p_protein,
  p_precursor,
  p_signal,
  p_mass,
  p_pca,
  p_norm,
  p_rt,
  p_depth_quality,
  p_qc_score
)

main_plots <- main_plots[!vapply(main_plots, is.null, logical(1))]

qc_main <- wrap_plots(main_plots, ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 9, face = "bold")
  )

save_svg(qc_main, "Fig_QC_main_publication_style_no_rings.svg", width = 18, height = 18)

global_qc_files <- character()

if (!is.null(p_global_protein)) {
  save_global_svg(p_global_protein, "Fig_QC_global_proteins_identified.svg", width = 12, height = 12)
  global_qc_files <- c(global_qc_files, "Fig_QC_global_proteins_identified.svg")
}

if (!is.null(p_global_precursor)) {
  save_global_svg(p_global_precursor, "Fig_QC_global_precursors_identified.svg", width = 12, height = 12)
  global_qc_files <- c(global_qc_files, "Fig_QC_global_precursors_identified.svg")
}

if (!is.null(p_global_protein) && !is.null(p_global_precursor)) {
  qc_global_depth_bars <- p_global_protein / p_global_precursor
  save_global_svg(qc_global_depth_bars, "Fig_QC_global_depth_bars.svg", width = 12, height = 12)
  global_qc_files <- c(global_qc_files, "Fig_QC_global_depth_bars.svg")
}

# ================================================================
# 15. Supplemental QC figure
# ================================================================

supp_plots <- list()

if (has_cols(c("FWHM.Scans", "FWHM.RT"))) {
  supp_plots$fwhm <- plot_scatter_qc(
    qc,
    "FWHM.Scans",
    "FWHM.RT",
    "FWHM scans",
    "FWHM RT"
  )
}

if (has_cols(c("Median.Mass.Acc.MS1.Corrected", "Median.Mass.Acc.MS2.Corrected"))) {
  supp_plots$mass_corrected <- plot_scatter_qc(
    qc,
    "Median.Mass.Acc.MS1.Corrected",
    "Median.Mass.Acc.MS2.Corrected",
    "Corrected MS1 mass accuracy (ppm)",
    "Corrected MS2 mass accuracy (ppm)"
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = neutral_cols["mid"]) +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = neutral_cols["mid"])
}

if ("Average.Peptide.Length" %in% names(qc)) {
  supp_plots$peptide_length <- plot_box_jitter(
    qc,
    "Average.Peptide.Length",
    "Average peptide length"
  )
}

if ("Average.Peptide.Charge" %in% names(qc)) {
  supp_plots$peptide_charge <- plot_box_jitter(
    qc,
    "Average.Peptide.Charge",
    "Average peptide charge"
  )
}

if ("Average.Missed.Tryptic.Cleavages" %in% names(qc)) {
  supp_plots$missed_cleavages <- plot_box_jitter(
    qc,
    "Average.Missed.Tryptic.Cleavages",
    "Missed tryptic cleavages"
  )
}

if (!is.null(p_missing)) {
  supp_plots$missingness <- p_missing
}

if (!is.null(p_cv)) {
  supp_plots$cv <- p_cv
}

if (length(supp_plots) > 0) {
  qc_supp <- wrap_plots(supp_plots, ncol = 3, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(
      legend.position = "bottom",
      plot.tag = element_text(size = 9, face = "bold")
    )

  save_svg(qc_supp, "Fig_QC_supplemental_publication_style_no_rings.svg", width = 18, height = 14)
}

# ================================================================
# 16. Dedicated outlier plot
# This keeps outlier information separate from the main manuscript figure
# ================================================================

p_outlier <- NULL

if ("max_abs_robust_z" %in% names(qc)) {
  p_outlier <- qc %>%
    arrange(max_abs_robust_z) %>%
    mutate(sample_order = factor(sample_id, levels = sample_id)) %>%
    ggplot(aes(x = sample_order, y = max_abs_robust_z, fill = celltype_layer)) +
    geom_col(width = 0.75, colour = NA) +
    geom_hline(
      yintercept = 3,
      linewidth = 0.3,
      linetype = "dashed",
      colour = neutral_cols["dark"]
    ) +
    scale_fill_manual(values = celltype_cols) +
    labs(
      x = NULL,
      y = "Maximum absolute robust z-score"
    ) +
    theme_publication_qc() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  save_svg(p_outlier, "Fig_QC_outlier_summary.svg", width = 12, height = 6)
}

# ================================================================
# 17. Correlation heatmap
# ================================================================

cor_plot <- NULL

if (length(pca_metrics) >= 3) {

  cor_mat <- qc %>%
    select(all_of(pca_metrics)) %>%
    cor(use = "pairwise.complete.obs")

  cor_df <- as.data.frame(as.table(cor_mat)) %>%
    rename(metric_x = Var1, metric_y = Var2, correlation = Freq)

  cor_plot <- ggplot(cor_df, aes(metric_x, metric_y, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#4C78A8",
      mid = "white",
      high = "#E45756",
      midpoint = 0,
      limits = c(-1, 1),
      name = "r"
    ) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_publication_qc() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  save_svg(cor_plot, "Fig_QC_correlation_heatmap.svg", width = 14, height = 12)
}

# ================================================================
# 18. Optional batch effect model
# ================================================================

batch_results <- NULL

if ("batch" %in% names(qc) && "Normalisation.Instability" %in% names(qc)) {

  batch_model <- lm(
    Normalisation.Instability ~ batch + celltype_layer,
    data = qc
  )

  batch_results <- as.data.frame(anova(batch_model))
  batch_results$term <- rownames(batch_results)
  rownames(batch_results) <- NULL

  write.csv(
    batch_results,
    file = file.path(table_dir, "qc_batch_model_normalisation_instability.csv"),
    row.names = FALSE
  )

  p_batch <- ggplot(qc, aes(x = factor(batch), y = Normalisation.Instability, fill = celltype_layer)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35, alpha = 0.75) +
    geom_point(
      aes(colour = celltype_layer),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
      size = 1,
      alpha = 0.65
    ) +
    scale_fill_manual(values = celltype_cols) +
    scale_colour_manual(values = celltype_cols, guide = "none") +
    labs(
      x = "Batch",
      y = "Normalisation instability"
    ) +
    theme_publication_qc()

  save_svg(p_batch, "Fig_QC_batch_normalisation_instability.svg", width = 10, height = 7)
}

# ================================================================
# 19. Summary tables
# ================================================================

qc_summary_by_celltype <- qc %>%
  group_by(celltype_layer) %>%
  summarise(
    n = n(),
    n_outliers = sum(qc_outlier, na.rm = TRUE),
    across(
      all_of(outlier_metrics),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        mad = ~ mad(.x, constant = 1.4826, na.rm = TRUE),
        cv = ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    QC.Score_mean = mean(QC.Score, na.rm = TRUE),
    QC.Score_sd = sd(QC.Score, na.rm = TRUE),
    .groups = "drop"
  )

qc_outliers <- qc_outlier_z %>%
  select(
    sample_id,
    celltype_layer,
    n_outlier_metrics,
    max_abs_robust_z,
    qc_outlier,
    all_of(z_cols)
  ) %>%
  arrange(desc(qc_outlier), desc(max_abs_robust_z))

write.csv(
  qc_summary_by_celltype,
  file = file.path(table_dir, "qc_summary_by_celltype_layer.csv"),
  row.names = FALSE
)

write.csv(
  qc_outliers,
  file = file.path(table_dir, "qc_robust_outlier_table.csv"),
  row.names = FALSE
)

if (!is.null(cv_summary)) {
  write.csv(
    cv_summary,
    file = file.path(table_dir, "qc_cv_summary.csv"),
    row.names = FALSE
  )
}

xlsx_list <- list(
  summary_by_celltype = qc_summary_by_celltype,
  robust_outliers = qc_outliers
)

if (!is.null(cv_summary)) {
  xlsx_list$cv_summary <- cv_summary
}

if (!is.null(batch_results)) {
  xlsx_list$batch_model <- batch_results
}

openxlsx::write.xlsx(
  xlsx_list,
  file = file.path(table_dir, "qc_summary_tables.xlsx"),
  overwrite = TRUE
)

write_run_manifest(
  file.path(CANONICAL_PATHS$logs, "run_manifest.yml"),
  inputs = list(qc_stats = input_file),
  outputs = list(figures = out_dir, tables = table_dir),
  parameters = list(dataset = DATASET),
  notes = "Dataset-aware sample-level quicksearch QC."
)

# ================================================================
# 20. Console output
# ================================================================

cat("\nQC manuscript figures saved to:\n")
cat(out_dir, "\n\n")

cat("Main figure:\n")
cat(" - Fig_QC_main_publication_style_no_rings.svg\n\n")

cat("Supplemental figure:\n")
cat(" - Fig_QC_supplemental_publication_style_no_rings.svg\n\n")

cat("Dedicated outlier figure:\n")
cat(" - Fig_QC_outlier_summary.svg\n\n")

cat("Global sample-depth figures:\n")
cat(global_out_dir, "\n")
if (length(global_qc_files) > 0) {
  cat(paste0(" - ", global_qc_files, collapse = "\n"), "\n\n")
} else {
  cat(" - none generated; required depth columns were not present.\n\n")
}

cat("Tables:\n")
cat(" - qc_summary_tables.xlsx\n")
cat(" - qc_summary_by_celltype_layer.csv\n")
cat(" - qc_robust_outlier_table.csv\n\n")

cat("Outlier rule:\n")
cat(" - Sample flagged if >=2 QC metrics have absolute robust z-score > 3 within celltype_layer.\n")
cat(" - Outliers are not overlaid on the main figure to avoid visual clutter.\n\n")

print(qc_main)
