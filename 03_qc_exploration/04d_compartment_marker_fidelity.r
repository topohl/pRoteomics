# ================================================================
# Script: 03_qc_exploration/04d_compartment_marker_fidelity.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional config/marker_panels/compartment_fidelity_marker_sets.csv; config/marker_panels/wgcna_reference_marker_sets.csv.
# Produces: results/tables/03_qc_exploration/04d_compartment_marker_fidelity/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Compartment marker fidelity QC.
# ================================================================

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

fidelity_protein_file_for_dataset <- function(dataset) {
  file.path(qc_paths(SOURCE_SUBSTEP_ID, dataset)$tables, "marker_fidelity_proteins_used_deduplicated.csv")
}

fidelity_files <- stats::setNames(vapply(compartment_datasets, fidelity_file_for_dataset, character(1)), compartment_datasets)
legacy_score_files <- stats::setNames(vapply(compartment_datasets, legacy_score_file_for_dataset, character(1)), compartment_datasets)
fidelity_protein_files <- stats::setNames(vapply(compartment_datasets, fidelity_protein_file_for_dataset, character(1)), compartment_datasets)

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
      paste0("Protein provenance table [", names(fidelity_protein_files), "]: ", fidelity_protein_files, ifelse(file.exists(fidelity_protein_files), " [PASS]", " [MISSING]")),
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


# Compact plotting defaults. Keep this generic so the codebase does not encode
# journal-specific naming while still producing clean figures.
#
# Use a generic graphics-device font by default. Explicit font names such as
# "Arial" can fail on Windows PDF/SVG devices with "invalid font type" when the
# font is installed for the OS but not registered for the active R graphics
# device. If a custom font is requested, use it only when the device can resolve
# it; otherwise fall back to the generic "sans" family.
safe_figure_font <- function(font = Sys.getenv("PROTEOMICS_FIGURE_FONT", unset = "sans")) {
  font <- trimws(as.character(font)[1])
  if (!nzchar(font)) return("sans")
  if (font %in% c("sans", "serif", "mono")) return(font)

  if (.Platform$OS.type == "windows") {
    known_windows_fonts <- tryCatch(names(grDevices::windowsFonts()), error = function(e) character())
    if (font %in% known_windows_fonts) return(font)

    warning(
      "Requested PROTEOMICS_FIGURE_FONT='", font,
      "' is not registered for the Windows R graphics device; using 'sans' instead.",
      call. = FALSE
    )
    return("sans")
  }

  font
}

figure_font <- safe_figure_font()

# Short labels avoid tilted x-axis text and make compact panels cleaner.
dataset_label_map <- c(
  neuron_soma = "SOMA",
  neuron_neuropil = "NEUROPIL",
  microglia = "MG/PVM"
)

# Muted, colorblind-safe qualitative palette. Used sparingly; most visual
# hierarchy is carried by geometry and spacing rather than saturated color.
dataset_palette <- c(
  neuron_soma = "#5B5A8C",
  neuron_neuropil = "#C9C5BF",
  microglia = "#D95F5F"
)

format_dataset_labels <- function(x) {
  x_chr <- as.character(x)
  out <- unname(dataset_label_map[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

compact_plot_theme <- function(base_size = 7, base_family = figure_font) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "grey12"),
      axis.title = ggplot2::element_text(size = base_size + 0.2, color = "grey10"),
      axis.text = ggplot2::element_text(size = base_size, color = "grey18"),
      axis.line = ggplot2::element_line(linewidth = 0.22, color = "grey18"),
      axis.ticks = ggplot2::element_line(linewidth = 0.22, color = "grey18"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = base_size + 0.2,
        face = "plain",
        color = "grey10",
        margin = ggplot2::margin(b = 2.5)
      ),
      legend.title = ggplot2::element_text(size = base_size - 0.2, color = "grey10"),
      legend.text = ggplot2::element_text(size = base_size - 0.4, color = "grey18"),
      legend.key.height = grid::unit(2.8, "mm"),
      legend.key.width = grid::unit(8.0, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title = ggplot2::element_text(size = base_size + 1.0, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.2, hjust = 0, color = "grey35"),
      plot.margin = ggplot2::margin(3, 4, 3, 3)
    )
}

save_compact_figure <- function(plot, filename, width_mm, height_mm) {
  ggplot2::ggsave(
    file.path(PATHS$figures, paste0(filename, ".svg")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = function(...) svglite::svglite(..., pointsize = 7, fix_text_size = FALSE),
    bg = "white"
  )

  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  pdf_args <- list(
    filename = file.path(PATHS$figures, paste0(filename, ".pdf")),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )
  if (identical(pdf_device, grDevices::pdf)) pdf_args$useDingbats <- FALSE
  do.call(ggplot2::ggsave, pdf_args)
}

pt_to_mm <- function(pt) {
  as.numeric(pt) * 25.4 / 72
}

finite_abs_limit <- function(x) {
  x <- abs(x[is.finite(x)])
  if (!length(x)) return(1)
  max(x, na.rm = TRUE)
}

format_p_value <- function(p, prefix = "adj. P") {
  out <- rep(paste0(prefix, " = NA"), length(p))
  finite <- is.finite(p)
  out[finite & p < 1e-4] <- paste0(prefix, " < 0.0001")
  out[finite & p >= 1e-4] <- paste0(prefix, " = ", formatC(p[finite & p >= 1e-4], format = "f", digits = 3))
  out
}

pairwise_marker_fidelity_tests <- function(scores, dataset_levels) {
  rows <- list()
  row_i <- 1L
  marker_classes <- fidelity_class_order[fidelity_class_order %in% as.character(unique(scores$marker_class))]

  for (cls in marker_classes) {
    cls_scores <- scores[as.character(scores$marker_class) == cls & is.finite(scores$marker_score), , drop = FALSE]
    datasets_present <- dataset_levels[dataset_levels %in% as.character(unique(cls_scores$dataset))]
    if (length(datasets_present) < 2) next

    comparisons <- utils::combn(datasets_present, 2, simplify = FALSE)
    for (comparison in comparisons) {
      x <- cls_scores$marker_score[as.character(cls_scores$dataset) == comparison[[1]]]
      y <- cls_scores$marker_score[as.character(cls_scores$dataset) == comparison[[2]]]
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]

      p_value <- NA_real_
      if (length(x) >= 2 && length(y) >= 2 && length(unique(c(x, y))) > 1) {
        p_value <- tryCatch(
          stats::wilcox.test(x, y, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }

      rows[[row_i]] <- data.frame(
        marker_class = cls,
        dataset_1 = comparison[[1]],
        dataset_2 = comparison[[2]],
        n_1 = length(x),
        n_2 = length(y),
        median_1 = if (length(x)) median(x, na.rm = TRUE) else NA_real_,
        median_2 = if (length(y)) median(y, na.rm = TRUE) else NA_real_,
        median_difference_1_minus_2 = if (length(x) && length(y)) median(x, na.rm = TRUE) - median(y, na.rm = TRUE) else NA_real_,
        p_value = p_value,
        test = "Wilcoxon rank-sum",
        p_adjust_method = "BH within marker_class",
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (!length(rows)) {
    return(data.frame(
      marker_class = character(),
      dataset_1 = character(),
      dataset_2 = character(),
      n_1 = integer(),
      n_2 = integer(),
      median_1 = numeric(),
      median_2 = numeric(),
      median_difference_1_minus_2 = numeric(),
      p_value = numeric(),
      p_adj_BH = numeric(),
      significant_BH_0_05 = logical(),
      p_label = character(),
      test = character(),
      p_adjust_method = character(),
      stringsAsFactors = FALSE
    ))
  }

  dplyr::bind_rows(rows) |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::mutate(
      p_adj_BH = stats::p.adjust(.data$p_value, method = "BH"),
      significant_BH_0_05 = is.finite(.data$p_adj_BH) & .data$p_adj_BH < 0.05,
      p_label = format_p_value(.data$p_adj_BH, prefix = "P")
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      factor(.data$marker_class, levels = fidelity_class_order, ordered = TRUE),
      .data$dataset_1,
      .data$dataset_2
    )
}

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

read_fidelity_protein_file <- function(dataset) {
  protein_file <- fidelity_protein_file_for_dataset(dataset)
  if (!file.exists(protein_file)) return(NULL)

  df <- tryCatch(qc_read_table(protein_file), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)

  optional_col <- function(name, default = NA_character_) {
    if (name %in% names(df)) df[[name]] else rep(default, nrow(df))
  }

  data.frame(
    dataset = dataset,
    fidelity_marker_class = as.character(optional_col("fidelity_marker_class")),
    marker_key = as.character(optional_col("marker_key")),
    requested_markers = as.character(optional_col("requested_markers")),
    matched_protein_id = as.character(optional_col("matched_protein_id")),
    used_in_fidelity_score = as.logical(optional_col("used_in_fidelity_score", FALSE)),
    contributing_marker_panels = as.character(optional_col("contributing_marker_panels")),
    contributing_fidelity_subpanels = as.character(optional_col("contributing_fidelity_subpanels")),
    source_types = as.character(optional_col("source_types")),
    source_names = as.character(optional_col("source_names")),
    source_references = as.character(optional_col("source_references")),
    source_terms_or_categories = as.character(optional_col("source_terms_or_categories")),
    evidence_levels = as.character(optional_col("evidence_levels")),
    confidences = as.character(optional_col("confidences")),
    include_in_fidelity_score_values = as.character(optional_col("include_in_fidelity_score_values")),
    n_marker_panels_containing_marker = suppressWarnings(as.integer(optional_col("n_marker_panels_containing_marker", NA_integer_))),
    n_matched_protein_ids = suppressWarnings(as.integer(optional_col("n_matched_protein_ids", NA_integer_))),
    fraction_nonmissing = suppressWarnings(as.numeric(optional_col("fraction_nonmissing", NA_real_))),
    mean_log2_abundance = suppressWarnings(as.numeric(optional_col("mean_log2_abundance", NA_real_))),
    stringsAsFactors = FALSE
  )
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

pairwise_tests <- pairwise_marker_fidelity_tests(compartment_scores, compartment_datasets)
qc_write_csv(pairwise_tests, file.path(PATHS$tables, "compartment_marker_fidelity_pairwise_tests.csv"))

message(
  "Compartment marker fidelity pairwise tests complete: ",
  sum(pairwise_tests$significant_BH_0_05 %in% TRUE, na.rm = TRUE),
  " BH-adjusted significant comparison(s) out of ",
  nrow(pairwise_tests),
  "."
)

protein_provenance <- dplyr::bind_rows(lapply(compartment_datasets, read_fidelity_protein_file))
if (!is.null(protein_provenance) && nrow(protein_provenance)) {
  protein_provenance <- protein_provenance |>
    dplyr::filter(!is.na(.data$fidelity_marker_class), nzchar(.data$fidelity_marker_class), !is.na(.data$marker_key), nzchar(.data$marker_key)) |>
    dplyr::arrange(.data$fidelity_marker_class, .data$dataset, .data$marker_key)
  qc_write_csv(protein_provenance, file.path(PATHS$tables, "compartment_marker_fidelity_protein_provenance.csv"))

  protein_source_summary <- protein_provenance |>
    dplyr::group_by(.data$dataset, .data$fidelity_marker_class, .data$source_names, .data$source_types, .data$evidence_levels, .data$confidences) |>
    dplyr::summarise(
      n_marker_keys = dplyr::n_distinct(.data$marker_key),
      n_used_marker_keys = dplyr::n_distinct(.data$marker_key[.data$used_in_fidelity_score %in% TRUE]),
      matched_protein_ids = collapse_unique(.data$matched_protein_id),
      contributing_marker_panels = collapse_unique(.data$contributing_marker_panels),
      source_terms_or_categories = collapse_unique(.data$source_terms_or_categories),
      source_references = collapse_unique(.data$source_references),
      .groups = "drop"
    ) |>
    dplyr::mutate(fraction_used_marker_keys = ifelse(.data$n_marker_keys > 0, .data$n_used_marker_keys / .data$n_marker_keys, NA_real_)) |>
    dplyr::arrange(.data$fidelity_marker_class, .data$dataset, .data$source_names)
  qc_write_csv(protein_source_summary, file.path(PATHS$tables, "compartment_marker_fidelity_source_summary.csv"))
} else {
  warning("No 04c marker_fidelity_proteins_used_deduplicated.csv files found; protein-level provenance summary was skipped.", call. = FALSE)
}

# Make compact display labels once and keep raw dataset keys for color mapping.
plot_scores <- compartment_scores |>
  dplyr::mutate(
    dataset_label = factor(format_dataset_labels(dataset), levels = format_dataset_labels(compartment_datasets)),
    dataset_x = as.numeric(.data$dataset_label),
    dataset_key = as.character(dataset),
    marker_class = factor(as.character(marker_class), levels = fidelity_class_order)
  )

dataset_x_breaks <- seq_along(compartment_datasets)
dataset_x_labels <- format_dataset_labels(compartment_datasets)

stat_annotations <- pairwise_tests |>
  dplyr::filter(.data$significant_BH_0_05 %in% TRUE) |>
  dplyr::mutate(
    marker_class = factor(.data$marker_class, levels = fidelity_class_order),
    x = match(.data$dataset_1, compartment_datasets),
    xend = match(.data$dataset_2, compartment_datasets),
    xmid = (.data$x + .data$xend) / 2
  ) |>
  dplyr::filter(!is.na(.data$x), !is.na(.data$xend))

if (nrow(stat_annotations)) {
  annotation_ranges <- plot_scores |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::summarise(
      data_y_min = min(.data$marker_score, na.rm = TRUE),
      data_y_max = max(.data$marker_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      data_y_span = pmax(.data$data_y_max - .data$data_y_min, 0.5)
    )

  stat_annotations <- stat_annotations |>
    dplyr::left_join(annotation_ranges, by = "marker_class") |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::arrange(.data$p_adj_BH, .by_group = TRUE) |>
    dplyr::mutate(
      annotation_rank = dplyr::row_number(),
      bracket_y = .data$data_y_max + .data$data_y_span * (0.14 + 0.17 * (.data$annotation_rank - 1)),
      tick_y = .data$bracket_y - .data$data_y_span * 0.045,
      label_y = .data$bracket_y + .data$data_y_span * 0.024
    ) |>
    dplyr::ungroup()

  stat_annotation_limits <- stat_annotations |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::summarise(y = max(.data$label_y + .data$data_y_span * 0.08, na.rm = TRUE), .groups = "drop")
} else {
  stat_annotation_limits <- data.frame(marker_class = character(), y = numeric())
}

plot_mean_sem <- plot_scores |>
  dplyr::group_by(.data$marker_class, .data$dataset_key, .data$dataset_x) |>
  dplyr::summarise(
    mean_marker_score = mean(.data$marker_score, na.rm = TRUE),
    sem_marker_score = stats::sd(.data$marker_score, na.rm = TRUE) / sqrt(sum(is.finite(.data$marker_score))),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    sem_marker_score = ifelse(is.nan(.data$sem_marker_score), NA_real_, .data$sem_marker_score),
    ymin = .data$mean_marker_score - .data$sem_marker_score,
    ymax = .data$mean_marker_score + .data$sem_marker_score,
    x = .data$dataset_x - 0.15,
    xend = .data$dataset_x + 0.15
  )

summary_plot_scores <- summary_scores |>
  dplyr::mutate(
    dataset_label = factor(format_dataset_labels(dataset), levels = format_dataset_labels(compartment_datasets)),
    marker_class = factor(as.character(marker_class), levels = fidelity_class_order)
  )

scatter_panel_width_pt <- 92
scatter_panel_height_pt <- 126
scatter_n_panels <- max(1L, dplyr::n_distinct(plot_scores$marker_class))
scatter_width_mm <- max(170, pt_to_mm(scatter_panel_width_pt * scatter_n_panels + 54))
scatter_height_mm <- max(88, pt_to_mm(scatter_panel_height_pt + 34))

# Distribution plot: light raw data + mean/SEM summary in the same visual
# language as the compact grouped dot plots used elsewhere in the project.
p_box <- ggplot2::ggplot(
  plot_scores,
  ggplot2::aes(x = dataset_x, y = marker_score, color = dataset_key)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "longdash",
    linewidth = 0.36,
    color = "grey30"
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.10, height = 0, seed = 42),
    alpha = 0.62,
    size = 1.75,
    stroke = 0
  ) +
  ggplot2::geom_errorbar(
    data = plot_mean_sem,
    ggplot2::aes(x = dataset_x, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    width = 0.20,
    linewidth = 0.46,
    color = "grey25"
  ) +
  ggplot2::geom_segment(
    data = plot_mean_sem,
    ggplot2::aes(x = x, xend = xend, y = mean_marker_score, yend = mean_marker_score),
    inherit.aes = FALSE,
    linewidth = 0.46,
    color = "grey25"
  ) +
  ggplot2::facet_wrap(~ marker_class, nrow = 1, scales = "free_y") +
  ggplot2::geom_blank(
    data = stat_annotation_limits,
    ggplot2::aes(x = 1, y = y),
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(values = dataset_palette, guide = "none", drop = FALSE) +
  ggplot2::scale_x_continuous(
    breaks = dataset_x_breaks,
    labels = dataset_x_labels,
    expand = ggplot2::expansion(mult = c(0.06, 0.08))
  ) +
  ggplot2::labs(title = "Compartment Marker Fidelity", x = NULL, y = "Mean marker\nabundance") +
  ggplot2::coord_cartesian(clip = "off") +
  compact_plot_theme(base_size = 7) +
  ggplot2::theme(
    axis.title.y = ggplot2::element_text(size = 7.2, color = "grey25", margin = ggplot2::margin(r = 5)),
    axis.text.x = ggplot2::element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5, color = "grey25"),
    axis.text.y = ggplot2::element_text(size = 7, color = "grey25"),
    axis.line = ggplot2::element_line(linewidth = 0.42, color = "grey25"),
    axis.ticks = ggplot2::element_line(linewidth = 0.42, color = "grey25"),
    axis.ticks.length = grid::unit(0.9, "mm"),
    strip.text = ggplot2::element_text(size = 7.2, color = "grey25", margin = ggplot2::margin(b = 3.5)),
    plot.title = ggplot2::element_text(size = 8, hjust = 0.5, color = "grey25", margin = ggplot2::margin(b = 8)),
    plot.margin = ggplot2::margin(5, 8, 5, 6),
    panel.spacing.x = grid::unit(6, "mm")
  )

if (nrow(stat_annotations)) {
  p_box <- p_box +
    ggplot2::geom_segment(
      data = stat_annotations,
      ggplot2::aes(x = x, xend = xend, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.38,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = stat_annotations,
      ggplot2::aes(x = x, xend = x, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.38,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = stat_annotations,
      ggplot2::aes(x = xend, xend = xend, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.38,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = stat_annotations,
      ggplot2::aes(x = xmid, y = label_y, label = p_label),
      inherit.aes = FALSE,
      size = 2.2,
      color = "grey25",
      vjust = 0
    )
}
save_compact_figure(
  p_box,
  "compartment_marker_fidelity_boxplots",
  width_mm = scatter_width_mm,
  height_mm = scatter_height_mm
)

heat_nx <- max(1L, dplyr::n_distinct(summary_plot_scores$dataset_label))
heat_ny <- max(1L, dplyr::n_distinct(summary_plot_scores$marker_class))
heat_tile_mm <- 11.0
heat_panel_width_mm <- heat_tile_mm * heat_nx
heat_panel_height_mm <- heat_tile_mm * heat_ny
heat_total_width_mm <- heat_panel_width_mm + 28
heat_total_height_mm <- heat_panel_height_mm + 14

# Summary heatmap: keep the figure visual and avoid tile labels; exact values are
# already written to compartment_marker_fidelity_summary.csv.
p_heat <- ggplot2::ggplot(
  summary_plot_scores,
  ggplot2::aes(x = dataset_label, y = marker_class, fill = median_marker_score)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::scale_fill_gradientn(
    colors = c("#F4F4F4", "#DBE6EF", "#86A9C5", "#3B6F9E"),
    na.value = "grey92"
  ) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Median\nscore") +
  compact_plot_theme(base_size = 7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(14, "mm"),
    legend.key.width = grid::unit(2.8, "mm"),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(3, 3, 3, 3)
  )
save_compact_figure(p_heat, "compartment_marker_fidelity_heatmap", width_mm = heat_total_width_mm, height_mm = heat_total_height_mm)

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

z_plot_scores <- z_scores |>
  dplyr::mutate(
    dataset_label = factor(format_dataset_labels(dataset), levels = format_dataset_labels(compartment_datasets)),
    marker_class = factor(as.character(marker_class), levels = fidelity_class_order)
  )

z_limit <- finite_abs_limit(z_plot_scores$median_marker_score_z)

p_z <- ggplot2::ggplot(
  z_plot_scores,
  ggplot2::aes(x = dataset_label, y = marker_class, fill = median_marker_score_z)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::scale_fill_gradient2(
    low = "#4E79A7",
    mid = "#F4F4F4",
    high = "#C44E52",
    midpoint = 0,
    limits = c(-z_limit, z_limit),
    na.value = "grey92"
  ) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Median\nz-score") +
  compact_plot_theme(base_size = 7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(14, "mm"),
    legend.key.width = grid::unit(2.8, "mm"),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(3, 3, 3, 3)
  )
save_compact_figure(p_z, "compartment_marker_fidelity_zscore_heatmap", width_mm = heat_total_width_mm, height_mm = heat_total_height_mm)

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
  "Protein/source audit tables:",
  "- compartment_marker_fidelity_protein_provenance.csv: cross-dataset 04c deduplicated marker proteins used/missing for each fidelity class, including source provenance when available.",
  "- compartment_marker_fidelity_source_summary.csv: compact counts by dataset, fidelity class, source_name/source_type/evidence_level/confidence.", "",
  "Statistics:",
  "- compartment_marker_fidelity_pairwise_tests.csv: pairwise Wilcoxon rank-sum tests between datasets within each marker class, with BH adjustment within marker class.",
  "- The distribution plot annotates only BH-adjusted significant comparisons (adjusted P < 0.05) to avoid visual clutter.", "",
  paste0("Datasets requested: ", paste(compartment_datasets, collapse = ", ")),
  paste0("Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", ")),
  paste0("Marker classes found: ", paste(unique(as.character(compartment_scores$marker_class)), collapse = ", ")),
  paste0("Source table(s): ", paste(unique(compartment_scores$source_table), collapse = ", ")),
  paste0("Source substep: ", SOURCE_SUBSTEP_ID)
)
writeLines(notes, file.path(PATHS$reports, "compartment_marker_fidelity_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(
    preferred_fidelity_files = as.list(fidelity_files),
    fallback_legacy_score_files = as.list(legacy_score_files),
    protein_provenance_files = as.list(fidelity_protein_files)
  ),
  outputs = list(
    figures = PATHS$figures,
    tables = PATHS$tables,
    reports = PATHS$reports,
    protein_provenance = file.path(PATHS$tables, "compartment_marker_fidelity_protein_provenance.csv"),
    source_summary = file.path(PATHS$tables, "compartment_marker_fidelity_source_summary.csv"),
    pairwise_tests = file.path(PATHS$tables, "compartment_marker_fidelity_pairwise_tests.csv")
  ),
  parameters = list(dataset = DATASET, compartment_datasets = compartment_datasets, marker_classes = fidelity_class_order),
  notes = "Cross-dataset soma/neuropil/microglia-PVM marker fidelity QC from 04c marker_fidelity_by_sample.csv outputs, with legacy fallback."
)

message("Compartment marker fidelity QC complete. Datasets found: ", paste(unique(as.character(compartment_scores$dataset)), collapse = ", "))
