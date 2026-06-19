#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/08_wgcna_score_publication_summary.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: results/tables/06_modules_WGCNA/module_score/<dataset>/wgcna/supermodule_*.csv; cleaned WGCNA labels from 06/07.
# Produces: results/figures/06_modules_WGCNA/score_publication_summary/<dataset>/publication_supermodule_*.svg/pdf.
# Notes: Final publication rendering for score-derived supermodule robustness/correlation plots; statistics remain produced by 03_score_module_activity.R.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr", "stringr", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) {
  suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
}

run <- wgcna_cli(allow_all = TRUE)
DATASET_ARG <- run$dataset
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}
module_source_arg <- tolower(arg_value("--module-source", Sys.getenv("PROTEOMICS_MODULE_DEFINITION_SOURCE", unset = "wgcna")))
if (!module_source_arg %in% c("wgcna")) {
  stop("08_wgcna_score_publication_summary.R currently renders WGCNA-source supermodule score outputs only.", call. = FALSE)
}

analysis_primary <- "primary_all_replicates"
analysis_qc_sensitivity <- "sensitivity_flagged_replicates_removed"
analysis_expected <- c(analysis_primary, analysis_qc_sensitivity)

analysis_display_label <- function(x) {
  dplyr::case_when(
    x == analysis_primary ~ "Primary: all technical replicates",
    x == analysis_qc_sensitivity ~ "Sensitivity: flagged technical replicates removed",
    TRUE ~ as.character(x)
  )
}

analysis_file_label <- function(x) {
  dplyr::case_when(
    x == analysis_primary ~ analysis_primary,
    x == analysis_qc_sensitivity ~ "sensitivity",
    TRUE ~ as.character(x)
  )
}

dataset_display_label <- function(dataset) {
  dplyr::case_when(
    dataset == "neuron_neuropil" ~ "Neuron neuropil",
    dataset == "neuron_soma" ~ "Neuron soma",
    dataset == "microglia" ~ "Microglia",
    TRUE ~ stringr::str_to_sentence(gsub("_", " ", as.character(dataset)))
  )
}

col_or_na <- function(df, nm) {
  if (nm %in% names(df)) return(df[[nm]])
  rep(NA_character_, nrow(df))
}

clean_label_value <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "NAN", "NULL")] <- NA_character_
  x
}

is_nonfinal_fallback_label <- function(x) {
  z <- stringr::str_squish(as.character(x))
  is.na(z) | !nzchar(z) |
    grepl("Hub-supported cluster|Hub-supported module cluster|Unresolved module cluster|Mixed / unresolved|mixed / low-specificity", z, ignore.case = TRUE)
}

supermodule_label_sep <- function() paste0(" ", intToUtf8(183), " ")

strip_supermodule_prefix <- function(x) {
  stringr::str_squish(sub(paste0("^SM[0-9]+\\s*(", intToUtf8(183), "|:|-)\\s*"), "", as.character(x), ignore.case = TRUE))
}

add_supermodule_id_prefix <- function(id, label) {
  id <- clean_label_value(id)
  label <- clean_label_value(label)
  out <- label
  needs_id <- !is.na(id) & !is.na(out) & !grepl(paste0("^SM[0-9]+\\s*(", intToUtf8(183), "|:|-)"), out)
  out[needs_id] <- paste0(id[needs_id], supermodule_label_sep(), out[needs_id])
  out
}

label_wrap <- function(x, width = 34) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

ordered_levels <- function(x, preferred = character()) {
  observed <- unique(as.character(x[!is.na(x)]))
  c(preferred[preferred %in% observed], sort(setdiff(observed, preferred)))
}

contrast_order <- c("RES vs CON", "SUS vs CON", "SUS vs RES")
heat_cols <- c(low = "#3B6FB6", mid = "#F8FAFC", high = "#C84C5A", ink = "#3F3F3F", grid = "#FFFFFF")

theme_publication_heat <- function(base_size = 7) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 6.5, color = heat_cols["ink"], angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.8, color = heat_cols["ink"], lineheight = 0.9),
      strip.text = ggplot2::element_text(size = 7.2, color = heat_cols["ink"], margin = ggplot2::margin(b = 3)),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 9.2, face = "plain", color = heat_cols["ink"], hjust = 0, margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = 6.8, color = "#6A6A6A", hjust = 0, margin = ggplot2::margin(b = 4)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6.5, color = heat_cols["ink"]),
      legend.text = ggplot2::element_text(size = 6, color = heat_cols["ink"]),
      legend.key.height = grid::unit(13, "mm"),
      legend.key.width = grid::unit(2.4, "mm"),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

scale_effect_fill <- function(name, limits) {
  ggplot2::scale_fill_gradient2(
    low = heat_cols["low"],
    mid = heat_cols["mid"],
    high = heat_cols["high"],
    midpoint = 0,
    limits = limits,
    oob = scales::squish,
    name = name,
    guide = ggplot2::guide_colorbar(
      frame.colour = heat_cols["ink"],
      frame.linewidth = 0.25,
      ticks.colour = heat_cols["ink"],
      ticks.linewidth = 0.25,
      barheight = grid::unit(20, "mm"),
      barwidth = grid::unit(2.5, "mm")
    )
  )
}

read_csv_required <- function(path, label) {
  if (!file.exists(path)) stop("Missing required ", label, ": ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

first_existing_csv <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else NA_character_
}

legacy_multifile_supermodule_labels_frozen <- function(dataset) {
  stop("Legacy multi-file label fallback is disabled; use WGCNA_final_label_lookup.csv.", call. = FALSE)
  candidates <- c(
    path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_supermodule_label_audit.csv"),
    path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_supermodule_plot_label_qc.csv"),
    path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_supermodule_group_effects_interpretable.csv"),
    path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_supermodule_biological_annotation.csv")
  )
  tabs <- lapply(candidates[file.exists(candidates)], function(path) {
    x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    x$.__label_source_file <- basename(path)
    x
  })
  if (!length(tabs)) stop("No cleaned supermodule label inputs found for dataset ", dataset, call. = FALSE)
  df <- dplyr::bind_rows(tabs)
  id <- dplyr::coalesce(
    clean_label_value(col_or_na(df, "SupermoduleID")),
    clean_label_value(col_or_na(df, "supermodule_id")),
    clean_label_value(col_or_na(df, "Supermodule_DataDriven"))
  )
  raw <- dplyr::coalesce(
    clean_label_value(col_or_na(df, "Supermodule_FinalLabel")),
    clean_label_value(col_or_na(df, "Supermodule_DisplayLabel")),
    clean_label_value(col_or_na(df, "Supermodule_PlotLabel")),
    id
  )
  candidates <- list(
    Supermodule_CleanPlotLabel = clean_label_value(col_or_na(df, "Supermodule_CleanPlotLabel")),
    Supermodule_CompositionDisplayLabel = add_supermodule_id_prefix(id, clean_label_value(col_or_na(df, "Supermodule_CompositionDisplayLabel"))),
    Supermodule_CompositionLabel = add_supermodule_id_prefix(id, clean_label_value(col_or_na(df, "Supermodule_CompositionLabel"))),
    Supermodule_PlotLabel = clean_label_value(col_or_na(df, "Supermodule_PlotLabel")),
    Supermodule_DisplayLabel = clean_label_value(col_or_na(df, "Supermodule_DisplayLabel")),
    Supermodule_FinalLabel = add_supermodule_id_prefix(id, clean_label_value(col_or_na(df, "Supermodule_FinalLabel"))),
    supermodule_id = id
  )
  best_label <- rep(NA_character_, length(id))
  best_source <- rep(NA_character_, length(id))
  for (nm in names(candidates)) {
    cand <- candidates[[nm]]
    good <- is.na(best_label) & !is.na(cand) & !is_nonfinal_fallback_label(cand)
    best_label[good] <- cand[good]
    best_source[good] <- nm
  }
  for (nm in names(candidates)) {
    cand <- candidates[[nm]]
    fill <- is.na(best_label) & !is.na(cand)
    best_label[fill] <- cand[fill]
    best_source[fill] <- nm
  }
  dplyr::tibble(
    supermodule_id = id,
    final_plot_label = best_label,
    final_plot_label_short = label_wrap(best_label, width = 34),
    final_label_source = best_source,
    raw_label = raw,
    raw_label_source = col_or_na(df, ".__label_source_file"),
    raw_label_rationale = dplyr::coalesce(
      clean_label_value(col_or_na(df, "Supermodule_CompositionRationale")),
      clean_label_value(col_or_na(df, "Supermodule_LabelRationale")),
      clean_label_value(col_or_na(df, "supermodule_theme_label_qc_warning"))
    ),
    cleaned_label_available = !is.na(dplyr::coalesce(
      clean_label_value(col_or_na(df, "Supermodule_CleanPlotLabel")),
      clean_label_value(col_or_na(df, "Supermodule_CompositionDisplayLabel")),
      clean_label_value(col_or_na(df, "Supermodule_CompositionLabel"))
    )),
    source_file = col_or_na(df, ".__label_source_file")
  ) |>
    dplyr::filter(!is.na(.data$supermodule_id), nzchar(.data$supermodule_id)) |>
    dplyr::arrange(.data$supermodule_id, dplyr::desc(.data$final_label_source == "Supermodule_CleanPlotLabel")) |>
    dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
}

legacy_multifile_module_labels_frozen <- function(dataset) {
  stop("Legacy multi-file label fallback is disabled; use WGCNA_final_label_lookup.csv.", call. = FALSE)
  candidates <- c(
    path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_module_plot_label_qc.csv"),
    path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_module_group_effects_interpretable.csv"),
    path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_module_biological_annotation.csv")
  )
  tabs <- lapply(candidates[file.exists(candidates)], function(path) {
    x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    x$.__label_source_file <- basename(path)
    x
  })
  if (!length(tabs)) return(tibble::tibble())
  df <- dplyr::bind_rows(tabs)
  id <- dplyr::coalesce(clean_label_value(col_or_na(df, "module_id")), clean_label_value(col_or_na(df, "ModuleID")))
  candidates <- list(
    Module_CleanPlotLabel = clean_label_value(col_or_na(df, "Module_CleanPlotLabel")),
    module_label_display = clean_label_value(col_or_na(df, "module_label_display")),
    cleaned_biological_label = clean_label_value(col_or_na(df, "cleaned_biological_label")),
    module_biological_label = clean_label_value(col_or_na(df, "module_biological_label")),
    ModuleLabel_Final = clean_label_value(col_or_na(df, "ModuleLabel_Final")),
    module_label = clean_label_value(col_or_na(df, "module_label")),
    module_id = id
  )
  best_label <- rep(NA_character_, length(id))
  best_source <- rep(NA_character_, length(id))
  for (nm in names(candidates)) {
    cand <- candidates[[nm]]
    good <- is.na(best_label) & !is.na(cand) & !is_nonfinal_fallback_label(cand)
    best_label[good] <- cand[good]
    best_source[good] <- nm
  }
  dplyr::tibble(
    module_id = id,
    final_plot_label = best_label,
    final_label_source = best_source,
    raw_label = dplyr::coalesce(clean_label_value(col_or_na(df, "ModuleLabel_Final")), clean_label_value(col_or_na(df, "module_label")), id),
    raw_label_source = col_or_na(df, ".__label_source_file"),
    raw_label_rationale = dplyr::coalesce(clean_label_value(col_or_na(df, "cleaned_biological_label_rationale")), clean_label_value(col_or_na(df, "GO_label_relevance_rationale")))
  ) |>
    dplyr::filter(!is.na(.data$module_id), nzchar(.data$module_id)) |>
    dplyr::distinct(.data$module_id, .keep_all = TRUE)
}

canonical_final_label_lookup <- function(dataset) {
  path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_final_label_lookup.csv")
  lookup <- read_csv_required(path, "canonical WGCNA final label lookup produced by 07_wgcna_interpretable_summary.r")
  wgcna_validate_label_lookup(lookup)
  if (!all(lookup$dataset == dataset)) {
    stop("Canonical WGCNA final label lookup contains rows for the wrong dataset: ", path, call. = FALSE)
  }
  lookup
}

canonical_supermodule_labels <- function(lookup) {
  lookup |>
    dplyr::filter(.data$level == "supermodule") |>
    dplyr::transmute(
      supermodule_id = .data$entity_id,
      final_plot_label = .data$final_plot_label,
      final_plot_label_short = label_wrap(.data$final_plot_label, width = 34),
      final_label_source = "WGCNA_final_label_lookup.csv",
      raw_label = dplyr::coalesce(.data$raw_top_GO_label, .data$entity_id),
      raw_label_source = "WGCNA_final_label_lookup.csv",
      raw_label_rationale = .data$label_rationale,
      cleaned_label_available = TRUE,
      source_file = "WGCNA_final_label_lookup.csv"
    )
}

canonical_module_labels <- function(lookup) {
  lookup |>
    dplyr::filter(.data$level == "module") |>
    dplyr::transmute(
      module_id = .data$entity_id,
      final_plot_label = .data$final_plot_label,
      final_label_source = "WGCNA_final_label_lookup.csv",
      raw_label = dplyr::coalesce(.data$raw_top_GO_label, .data$entity_id),
      raw_label_source = "WGCNA_final_label_lookup.csv",
      raw_label_rationale = .data$label_rationale
    )
}

validate_render_ids_in_lookup <- function(ids, labels, id_col, context) {
  ids <- unique(as.character(ids[!is.na(ids) & nzchar(as.character(ids))]))
  available <- unique(as.character(labels[[id_col]]))
  missing <- setdiff(ids, available)
  if (length(missing)) {
    stop("Rendered ", context, " IDs are absent from WGCNA_final_label_lookup.csv: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

attach_supermodule_labels <- function(df, labels, module_col = "Module") {
  validate_render_ids_in_lookup(df[[module_col]], labels, "supermodule_id", module_col)
  out <- df |>
    dplyr::left_join(labels, by = stats::setNames("supermodule_id", module_col)) |>
    dplyr::mutate(
      raw_score_module_id = .data[[module_col]],
      final_plot_label_short = label_wrap(.data$final_plot_label, width = 34)
    )
  out
}

score_effect_title <- function(dataset) {
  paste(dataset_display_label(dataset), "WGCNA supermodule eigengene effect sizes")
}

plot_supermodule_effect_heatmap <- function(df, analysis_label, dataset) {
  plot_df <- df |>
    dplyr::filter(.data$Analysis == analysis_label, is.finite(.data$Cohen_d))
  if (!nrow(plot_df)) return(NULL)

  region_levels <- ordered_levels(plot_df$RegionLayer)
  module_levels <- labels_by_id(plot_df, "Module")
  contrast_levels <- ordered_levels(plot_df$contrast_label, contrast_order)

  plot_df <- plot_df |>
    dplyr::mutate(
      RegionLayer = factor(.data$RegionLayer, levels = region_levels),
      Module_label = factor(.data$final_plot_label_short, levels = rev(module_levels$label)),
      contrast_label = factor(.data$contrast_label, levels = contrast_levels),
      show_label = abs(.data$Cohen_d) >= 0.8,
      d_label = sprintf("%.1f", .data$Cohen_d)
    ) |>
    dplyr::filter(!is.na(.data$RegionLayer), !is.na(.data$Module_label), !is.na(.data$contrast_label))

  max_abs_d <- stats::quantile(abs(plot_df$Cohen_d), 0.95, na.rm = TRUE)
  max_abs_d <- min(max(1, max_abs_d), 2.5)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$RegionLayer, y = .data$Module_label, fill = .data$Cohen_d)) +
    ggplot2::geom_tile(color = heat_cols["grid"], linewidth = 0.55, width = 0.94, height = 0.94) +
    ggplot2::geom_text(data = plot_df |> dplyr::filter(.data$show_label), ggplot2::aes(label = .data$d_label), size = 1.85, color = heat_cols["ink"]) +
    ggplot2::geom_point(data = plot_df |> dplyr::filter(.data$within_BH_significant), shape = 21, size = 1.55, color = heat_cols["ink"], fill = "white", stroke = 0.28) +
    ggplot2::facet_wrap(~ contrast_label, nrow = 1) +
    scale_effect_fill("Cohen's d", limits = c(-max_abs_d, max_abs_d)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = score_effect_title(dataset),
      subtitle = paste0(analysis_display_label(analysis_label), " | dots indicate BH-adjusted p <= 0.05; numbers shown for |d| >= 0.8")
    ) +
    theme_publication_heat(base_size = 7)
}

labels_by_id <- function(df, id_col) {
  df |>
    dplyr::transmute(id = as.character(.data[[id_col]]), label = .data$final_plot_label_short) |>
    dplyr::distinct(.data$id, .keep_all = TRUE) |>
    dplyr::arrange(.data$id)
}

plot_supermodule_correlation_heatmap <- function(df, analysis_label) {
  plot_df <- df |>
    dplyr::filter(.data$Analysis == analysis_label, is.finite(.data$rho))
  if (!nrow(plot_df)) return(NULL)

  label_map <- dplyr::bind_rows(
    plot_df |> dplyr::transmute(id = .data$module1, label = .data$module1_label_short),
    plot_df |> dplyr::transmute(id = .data$module2, label = .data$module2_label_short)
  ) |>
    dplyr::distinct(.data$id, .keep_all = TRUE) |>
    dplyr::arrange(.data$id)
  module_label_levels <- label_map$label
  region_levels <- ordered_levels(plot_df$RegionLayer)

  plot_df <- plot_df |>
    dplyr::transmute(
      RegionLayer,
      module1 = .data$module1_label_short,
      module2 = .data$module2_label_short,
      rho,
      p_adj_BH
    )

  mirror_df <- plot_df |>
    dplyr::transmute(RegionLayer, module1 = module2, module2 = module1, rho, p_adj_BH)

  diag_df <- expand.grid(
    RegionLayer = unique(plot_df$RegionLayer),
    module1 = module_label_levels,
    module2 = module_label_levels,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(.data$module1 == .data$module2) |>
    dplyr::mutate(rho = 1, p_adj_BH = 0)

  full_df <- dplyr::bind_rows(plot_df, mirror_df, diag_df) |>
    dplyr::mutate(
      RegionLayer = factor(.data$RegionLayer, levels = region_levels),
      module1 = factor(.data$module1, levels = module_label_levels),
      module2 = factor(.data$module2, levels = rev(module_label_levels))
    ) |>
    dplyr::filter(!is.na(.data$RegionLayer), !is.na(.data$module1), !is.na(.data$module2))

  ggplot2::ggplot(full_df, ggplot2::aes(x = .data$module1, y = .data$module2, fill = .data$rho)) +
    ggplot2::geom_tile(color = heat_cols["grid"], linewidth = 0.55, width = 0.94, height = 0.94) +
    ggplot2::geom_point(data = full_df |> dplyr::filter(.data$module1 != .data$module2, .data$p_adj_BH <= 0.05), shape = 21, size = 1.45, color = heat_cols["ink"], fill = "white", stroke = 0.28) +
    ggplot2::facet_wrap(~ RegionLayer) +
    scale_effect_fill("Spearman\nrho", limits = c(-1, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Supermodule correlation structure",
      subtitle = paste0(analysis_display_label(analysis_label), " | dots indicate BH-adjusted p <= 0.05")
    ) +
    theme_publication_heat(base_size = 6.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 5.5, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 5.5, lineheight = 0.9)
    )
}

plot_supermodule_consistency <- function(df, analysis_label) {
  plot_df <- df |>
    dplyr::filter(.data$Analysis == analysis_label, is.finite(.data$direction_consistency))
  if (!nrow(plot_df)) return(NULL)

  module_levels <- labels_by_id(plot_df, "Module")
  contrast_levels <- ordered_levels(plot_df$contrast_label, contrast_order)
  plot_df <- plot_df |>
    dplyr::mutate(
      Module_label = factor(.data$final_plot_label_short, levels = rev(module_levels$label)),
      contrast_label = factor(.data$contrast_label, levels = contrast_levels),
      evidence_class = factor(
        .data$evidence_class,
        levels = c("convergent_with_variable_power", "directionally_stable_trend", "limited_evidence", "mixed_or_unstable_direction")
      )
    ) |>
    dplyr::filter(!is.na(.data$Module_label), !is.na(.data$contrast_label))

  evidence_cols <- c(
    convergent_with_variable_power = heat_cols["high"],
    directionally_stable_trend = heat_cols["low"],
    limited_evidence = "#8A8A84",
    mixed_or_unstable_direction = "#7A5A8A"
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$direction_consistency, y = .data$Module_label, color = .data$evidence_class, size = .data$mean_abs_Cohen_d)) +
    ggplot2::geom_vline(xintercept = 0.80, linetype = "dashed", linewidth = 0.30, color = heat_cols["ink"]) +
    ggplot2::geom_point(alpha = 0.90, stroke = 0) +
    ggplot2::facet_wrap(~ contrast_label, nrow = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::scale_color_manual(values = evidence_cols, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(1.4, 3.4), name = "|mean d|") +
    ggplot2::labs(
      title = "Supermodule directional consistency across regions",
      subtitle = paste0(analysis_display_label(analysis_label), " | dashed line marks 80% directional agreement"),
      x = "Fraction in dominant direction",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linewidth = 0.22, color = "#EAE6DD"),
      panel.spacing = grid::unit(3.5, "mm"),
      plot.title = ggplot2::element_text(size = 9.2, color = heat_cols["ink"], hjust = 0, margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = 6.8, color = "#6A6A6A", hjust = 0, margin = ggplot2::margin(b = 4)),
      axis.title.x = ggplot2::element_text(size = 7, color = heat_cols["ink"], margin = ggplot2::margin(t = 4)),
      axis.text = ggplot2::element_text(size = 6.5, color = heat_cols["ink"]),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7.2, color = heat_cols["ink"], margin = ggplot2::margin(b = 3)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6.5, color = heat_cols["ink"]),
      legend.text = ggplot2::element_text(size = 6, color = heat_cols["ink"])
    )
}

save_publication_plot <- function(plot, out_dir, stem, dataset, module_source, analysis_label, width, height) {
  if (is.null(plot)) return(character())
  dir_create(out_dir)
  out <- character()
  for (ext in c("svg", "pdf")) {
    target <- file.path(out_dir, paste0(stem, "_", dataset, "_", module_source, "_", analysis_file_label(analysis_label), ".", ext))
    ggplot2::ggsave(target, plot = plot, width = width, height = height, units = "mm")
    out <- c(out, target)
  }
  out
}

validate_unchanged_numeric <- function(original, labelled, key_cols, numeric_cols) {
  key_cols <- intersect(key_cols, intersect(names(original), names(labelled)))
  numeric_cols <- intersect(numeric_cols, intersect(names(original), names(labelled)))
  if (!length(key_cols) || !length(numeric_cols)) {
    return(tibble::tibble(metric = character(), max_abs_diff = numeric(), status = character()))
  }
  original_keyed <- original |> dplyr::arrange(dplyr::across(dplyr::all_of(key_cols)))
  labelled_keyed <- labelled |> dplyr::arrange(dplyr::across(dplyr::all_of(key_cols)))
  key_ok <- identical(original_keyed[key_cols], labelled_keyed[key_cols])
  dplyr::bind_rows(lapply(numeric_cols, function(nm) {
    a <- suppressWarnings(as.numeric(original_keyed[[nm]]))
    b <- suppressWarnings(as.numeric(labelled_keyed[[nm]]))
    tibble::tibble(
      metric = nm,
      max_abs_diff = if (length(a) && length(b)) max(abs(a - b), na.rm = TRUE) else NA_real_,
      status = if (isTRUE(key_ok) && isTRUE(all.equal(a, b, tolerance = 0, check.attributes = FALSE))) "unchanged" else if (isTRUE(key_ok) && max(abs(a - b), na.rm = TRUE) <= 1e-12) "unchanged_roundoff" else "changed_or_key_mismatch"
    )
  }))
}

run_one_dataset <- function(dataset, module_source = "wgcna") {
  score_dir <- path_results("tables", "06_modules_WGCNA", "module_score", dataset, module_source)
  paths <- list(
    directional_effects = file.path(score_dir, "supermodule_directional_effects.csv"),
    directional_consistency = file.path(score_dir, "supermodule_directional_consistency.csv"),
    correlation_structure = file.path(score_dir, "supermodule_score_correlation_structure.csv"),
    eigengene_scores = file.path(score_dir, "supermodule_scores_eigengene_per_sample.csv"),
    composition = file.path(score_dir, "supermodule_score_composition.csv"),
    full_workbook = file.path(score_dir, "WGCNA_supermodule_score_statistics_full.xlsx")
  )

  out_paths <- create_module_dirs("06_modules_WGCNA", file.path("score_publication_summary", dataset))
  out_fig <- out_paths$figures
  out_tables <- out_paths$tables
  out_source <- out_paths$source_data
  out_logs <- out_paths$logs

  if (run$dry_run) {
    dry_run_line("Script", "06_modules_WGCNA/08_wgcna_score_publication_summary.R")
    dry_run_line("Dataset", dataset)
    dry_run_line("Module source", module_source)
    dry_run_line("Score directional effects", paths$directional_effects, if (file.exists(paths$directional_effects)) "PASS" else "FAIL")
    label_lookup_path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_final_label_lookup.csv")
    dry_run_line("Canonical final label lookup", label_lookup_path, if (file.exists(label_lookup_path)) "PASS" else "FAIL")
    dry_run_line("Output figures", out_fig)
    return(invisible(NULL))
  }

  effects_raw <- read_csv_required(paths$directional_effects, "supermodule directional effects")
  consistency_raw <- read_csv_required(paths$directional_consistency, "supermodule directional consistency")
  corr_raw <- read_csv_required(paths$correlation_structure, "supermodule score correlation structure")
  final_label_lookup <- canonical_final_label_lookup(dataset)
  labels <- canonical_supermodule_labels(final_label_lookup)
  module_labels <- canonical_module_labels(final_label_lookup)
  validate_render_ids_in_lookup(effects_raw$Module, labels, "supermodule_id", "effect")
  validate_render_ids_in_lookup(consistency_raw$Module, labels, "supermodule_id", "consistency")
  validate_render_ids_in_lookup(c(corr_raw$module1, corr_raw$module2), labels, "supermodule_id", "correlation")

  effects <- attach_supermodule_labels(effects_raw, labels, "Module")
  consistency <- attach_supermodule_labels(consistency_raw, labels, "Module")
  corr <- corr_raw |>
    dplyr::left_join(labels |> dplyr::select(module1 = "supermodule_id", module1_label = "final_plot_label", module1_label_short = "final_plot_label_short"), by = "module1") |>
    dplyr::left_join(labels |> dplyr::select(module2 = "supermodule_id", module2_label = "final_plot_label", module2_label_short = "final_plot_label_short"), by = "module2") |>
    dplyr::mutate(
      module1_label_short = label_wrap(.data$module1_label, width = 34),
      module2_label_short = label_wrap(.data$module2_label, width = 34)
    )

  analysis_seen <- sort(unique(c(effects$Analysis, consistency$Analysis, corr$Analysis)))
  missing_analysis <- setdiff(analysis_expected, analysis_seen)
  if (length(missing_analysis)) {
    stop("Missing expected primary/sensitivity analyses in score-publication inputs for ", dataset, ": ", paste(missing_analysis, collapse = ", "), call. = FALSE)
  }

  readr::write_csv(labels, file.path(out_tables, "WGCNA_score_publication_supermodule_label_lookup.csv"), na = "")
  readr::write_csv(labels, file.path(out_source, "WGCNA_score_publication_supermodule_label_lookup.csv"), na = "")
  if (nrow(module_labels)) {
    readr::write_csv(module_labels, file.path(out_tables, "WGCNA_score_publication_module_label_lookup.csv"), na = "")
    readr::write_csv(module_labels, file.path(out_source, "WGCNA_score_publication_module_label_lookup.csv"), na = "")
  }
  readr::write_csv(effects, file.path(out_source, "publication_supermodule_effect_heatmap_source.csv"), na = "")
  readr::write_csv(consistency, file.path(out_source, "publication_supermodule_consistency_source.csv"), na = "")
  readr::write_csv(corr, file.path(out_source, "publication_supermodule_correlation_source.csv"), na = "")

  rendered <- character()
  for (analysis_label in analysis_expected) {
    rendered <- c(rendered, save_publication_plot(plot_supermodule_effect_heatmap(effects, analysis_label, dataset), out_fig, "publication_supermodule_effect_heatmap", dataset, module_source, analysis_label, 185, 85))
    rendered <- c(rendered, save_publication_plot(plot_supermodule_consistency(consistency, analysis_label), out_fig, "publication_supermodule_consistency", dataset, module_source, analysis_label, 175, 75))
    rendered <- c(rendered, save_publication_plot(plot_supermodule_correlation_heatmap(corr, analysis_label), out_fig, "publication_supermodule_correlation", dataset, module_source, analysis_label, 185, 105))
  }

  label_qc <- dplyr::bind_rows(
    effects |> dplyr::distinct(plot_type = "effect_heatmap", supermodule_id = .data$Module, final_plot_label, raw_label, cleaned_label_available),
    consistency |> dplyr::distinct(plot_type = "consistency", supermodule_id = .data$Module, final_plot_label, raw_label, cleaned_label_available),
    corr |> dplyr::distinct(plot_type = "correlation_module1", supermodule_id = .data$module1, final_plot_label = .data$module1_label, raw_label = .data$module1, cleaned_label_available = !is_nonfinal_fallback_label(.data$module1_label)),
    corr |> dplyr::distinct(plot_type = "correlation_module2", supermodule_id = .data$module2, final_plot_label = .data$module2_label, raw_label = .data$module2, cleaned_label_available = !is_nonfinal_fallback_label(.data$module2_label))
  ) |>
    dplyr::mutate(
      fallback_visible_with_clean_label = .data$cleaned_label_available & is_nonfinal_fallback_label(.data$final_plot_label),
      final_label_matches_lookup = .data$final_plot_label %in% labels$final_plot_label
    )

  numeric_qc <- dplyr::bind_rows(
    validate_unchanged_numeric(
      effects_raw,
      effects,
      c("Analysis", "RegionLayer", "Module", "group1", "group2", "contrast_label"),
      c("estimate", "SE", "df", "t.ratio", "p.value", "p_adj_within_model_BH", "p_adj_global_BH", "Cohen_d", "rank_biserial", "estimate_group2_minus_group1", "p_nominal", "p_within_BH", "p_global_BH")
    ) |> dplyr::mutate(table = "supermodule_directional_effects"),
    validate_unchanged_numeric(
      consistency_raw,
      consistency,
      c("Analysis", "Module", "group1", "group2", "contrast_label"),
      c("n_tests", "n_positive", "n_negative", "direction_consistency", "mean_Cohen_d", "median_Cohen_d", "mean_abs_Cohen_d", "fisher_p", "stouffer_signed_z", "stouffer_p")
    ) |> dplyr::mutate(table = "supermodule_directional_consistency"),
    validate_unchanged_numeric(
      corr_raw,
      corr,
      c("Analysis", "RegionLayer", "module1", "module2"),
      c("n", "rho", "p", "p_adj_BH")
    ) |> dplyr::mutate(table = "supermodule_score_correlation_structure")
  ) |>
    dplyr::select("table", dplyr::everything())

  rendered_qc <- tibble::tibble(
    dataset = dataset,
    module_source = module_source,
    rendered_file = rendered,
    analysis = dplyr::case_when(
      grepl("_primary_all_replicates\\.", basename(rendered)) ~ analysis_primary,
      grepl("_sensitivity\\.", basename(rendered)) ~ analysis_qc_sensitivity,
      TRUE ~ NA_character_
    ),
    plot_type = dplyr::case_when(
      grepl("effect_heatmap", basename(rendered)) ~ "effect_heatmap",
      grepl("correlation", basename(rendered)) ~ "correlation",
      grepl("consistency", basename(rendered)) ~ "consistency",
      TRUE ~ "unknown"
    )
  )

  validation <- tibble::tibble(
    validation_check = c(
      "no_nonfinal_fallback_labels_with_clean_labels",
      "final_labels_match_cleaned_lookup",
      "primary_and_sensitivity_rendered",
      "numeric_statistics_unchanged"
    ),
    validation_status = c(
      if (any(label_qc$fallback_visible_with_clean_label, na.rm = TRUE)) "fail" else "ok",
      if (all(label_qc$final_label_matches_lookup, na.rm = TRUE)) "ok" else "fail",
      if (all(analysis_expected %in% rendered_qc$analysis)) "ok" else "fail",
      if (all(numeric_qc$status %in% c("unchanged", "unchanged_roundoff"))) "ok" else "fail"
    ),
    validation_message = c(
      paste(label_qc$supermodule_id[label_qc$fallback_visible_with_clean_label], collapse = "; "),
      paste(unique(label_qc$supermodule_id[!label_qc$final_label_matches_lookup]), collapse = "; "),
      paste(setdiff(analysis_expected, rendered_qc$analysis), collapse = "; "),
      paste(unique(numeric_qc$table[!numeric_qc$status %in% c("unchanged", "unchanged_roundoff")]), collapse = "; ")
    )
  )

  readr::write_csv(label_qc, file.path(out_tables, "WGCNA_score_publication_label_qc.csv"), na = "")
  readr::write_csv(label_qc, file.path(out_source, "WGCNA_score_publication_label_qc.csv"), na = "")
  readr::write_csv(numeric_qc, file.path(out_tables, "WGCNA_score_publication_numeric_qc.csv"), na = "")
  readr::write_csv(rendered_qc, file.path(out_tables, "WGCNA_score_publication_rendered_files.csv"), na = "")
  readr::write_csv(validation, file.path(out_tables, "WGCNA_score_publication_validation.csv"), na = "")

  if (any(validation$validation_status != "ok")) {
    stop("Score-publication validation failed for ", dataset, ": ", paste(validation$validation_check[validation$validation_status != "ok"], collapse = ", "), call. = FALSE)
  }

  write_run_manifest(
    file.path(out_logs, "run_manifest.yml"),
    inputs = c(paths, list(
      final_label_lookup = path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_final_label_lookup.csv")
    )),
    outputs = list(figures = out_fig, tables = out_tables, source_data = out_source),
    parameters = list(dataset = dataset, module_source = module_source, analyses = paste(analysis_expected, collapse = ";")),
    notes = "Final score-derived WGCNA supermodule publication plots. Score/statistics computation stays in 03_score_module_activity.R; every rendered label comes only from the canonical lookup produced by 07."
  )

  message("WGCNA score publication summary complete for dataset: ", dataset)
  invisible(list(labels = labels, validation = validation, rendered = rendered_qc))
}

datasets_to_run <- if (identical(DATASET_ARG, "all")) valid_datasets() else DATASET_ARG
invisible(lapply(datasets_to_run, run_one_dataset, module_source = module_source_arg))
