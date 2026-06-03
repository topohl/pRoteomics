#!/usr/bin/env Rscript
#
# Combine WGCNA group effects and biological annotation into interpretable summaries.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- wgcna_cli(allow_all = TRUE)
DATASET_ARG <- run$dataset

if (run$dry_run) {
  datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
  for (ds in datasets) {
    paths <- wgcna_downstream_paths("interpretable_summary", ds)
    invisible(lapply(unlist(paths), dir_create))
    dry_run_line("Dataset", ds)
    dry_run_line("Module effects", path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"), if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))) "PASS" else "WARN")
    dry_run_line("Supermodule effects", path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"), if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))) "PASS" else "WARN")
    dry_run_line("Module annotation", path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"), if (file.exists(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))) "PASS" else "WARN")
  }
  if (DATASET_ARG == "all") dry_run_line("Cross-dataset output", path_results("tables", "06_modules_WGCNA", "interpretable_summary", "all"))
  quit(status = 0, save = "no")
}

dataset_label <- function(ds) {
  switch(ds, neuron_neuropil = "Neuron neuropil", neuron_soma = "Neuron soma", microglia = "Microglia ROI", ds)
}

effect_sentence <- function(row, level = "supermodule") {
  row_get <- function(name, default = NA) {
    if (name %in% names(row)) row[[name]][[1]] else default
  }
  ds_label <- dataset_label(row$dataset)
  id <- if (level == "supermodule") row_get("supermodule_id") else row_get("module_id")
  label <- if (level == "supermodule") row_get("Supermodule_FinalLabel", row_get("supermodule_label")) else row_get("module_label", row_get("ModuleID"))
  cls <- row_get("dominant_microenvironment_class", row_get("microenvironment_class"))
  estimate <- suppressWarnings(as.numeric(row_get("estimate")))
  direction <- ifelse(is.na(estimate), "is altered", ifelse(estimate > 0, "is higher", "is lower"))
  spatial_unit <- row_get("spatial_unit")
  spatial <- ifelse(is.na(spatial_unit) || spatial_unit == "global_spatial_adjusted", "after spatial adjustment", paste0("in ", spatial_unit))
  effect_scope <- row_get("effect_scope")
  scope <- ifelse(!is.na(effect_scope), paste0(" (", effect_scope, ")"), "")
  base <- paste(ds_label, level, id, "annotated as", label, direction, "for", row_get("contrast"), spatial, scope)
  base <- gsub("[[:space:]]+\\(", " (", base)
  if (row$dataset != "microglia" || is.na(cls)) return(paste0(base, "."))
  if (cls == "microglia_supported") return(paste0(base, "; this ROI signal has microglia-supported ROI signal evidence."))
  if (cls == "microglia_state_or_activation_supported") return(paste0(base, "; this ROI signal has microglia state/phagolysosomal support."))
  if (cls == "shared_microenvironment") return(paste0(base, "; this ROI signal is shared local microenvironment signal, not purified microglial regulation."))
  if (cls == "neuropil_sensitive") return(paste0(base, "; this ROI signal is neuropil-sensitive; do not interpret as purified microglial regulation."))
  if (cls == "other_cellular_or_vascular_sensitive") return(paste0(base, "; this ROI signal shows other cellular/vascular marker support; interpret cautiously."))
  paste0(base, "; microenvironment support is ambiguous.")
}

add_interpretation_sentences <- function(df, level) {
  df$interpretation_sentence <- vapply(seq_len(nrow(df)), function(i) effect_sentence(df[i, , drop = FALSE], level), character(1))
  df
}

make_dataset_summary <- function(ds) {
  paths <- wgcna_downstream_paths("interpretable_summary", ds)
  module_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))
  super_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))
  module_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))
  super_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"))
  overlap <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", ds, "WGCNA_vs_DE_GSEA_overlap.csv"))

  if (is.null(module_effects)) module_effects <- empty_group_effects(ds, "module", "missing module_group_effects.csv")
  if (is.null(super_effects)) super_effects <- empty_group_effects(ds, "supermodule", "missing supermodule_group_effects.csv")
  if (is.null(module_annot)) module_annot <- data.frame(dataset = ds, ModuleID = character(), microenvironment_class = character())
  if (is.null(super_annot)) super_annot <- data.frame(dataset = ds, SupermoduleID = character(), Supermodule_FinalLabel = character(), dominant_microenvironment_class = character())
  for (nm in c(
    "dominant_microenvironment_class", "fraction_modules_microglia_supported",
    "fraction_modules_microglia_state_or_activation_supported", "fraction_modules_shared_microenvironment",
    "fraction_modules_neuropil_sensitive", "fraction_modules_other_cellular_or_vascular_sensitive",
    "fraction_modules_ambiguous", "marker_registry_version", "empirical_marker_set_version",
    "classification_rationale", "Supermodule_LabelSource", "Supermodule_LabelConfidence", "ManualReviewRequired"
  )) {
    if (!nm %in% names(super_annot)) super_annot[[nm]] <- NA
  }
  super_annot <- super_annot |>
    dplyr::arrange(.data$SupermoduleID) |>
    dplyr::distinct(.data$dataset, .data$SupermoduleID, .keep_all = TRUE)

  module_join <- module_effects |>
    dplyr::left_join(module_annot, by = c("dataset" = "dataset", "module_id" = "ModuleID")) |>
    add_interpretation_sentences("module")

  super_join <- super_effects |>
    dplyr::left_join(super_annot, by = c("dataset" = "dataset", "supermodule_id" = "SupermoduleID")) |>
    add_interpretation_sentences("supermodule")

  if (!is.null(overlap) && nrow(overlap) && "ModuleID" %in% names(overlap)) {
    overlap_summary <- overlap |>
      dplyr::group_by(.data$ModuleID) |>
      dplyr::summarise(best_wgcna_de_gsea_overlap = paste(utils::head(unique(.data$contrast), 5), collapse = ";"), .groups = "drop")
    module_join <- module_join |> dplyr::left_join(overlap_summary, by = c("module_id" = "ModuleID"))
  }

  top_super <- super_join |>
    dplyr::filter(!is.na(.data$p_value)) |>
    dplyr::arrange(.data$FDR_global, .data$p_value, dplyr::desc(abs(.data$estimate))) |>
    dplyr::slice_head(n = 50)

  write_table_and_source(super_join, paths$tables, paths$source_data, "WGCNA_supermodule_group_effects_interpretable.csv")
  write_table_and_source(module_join, paths$tables, paths$source_data, "WGCNA_module_group_effects_interpretable.csv")
  write_table_and_source(top_super, paths$tables, paths$source_data, "WGCNA_top_changed_supermodules.csv")
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(list(supermodules = super_join, modules = module_join, top_supermodules = top_super), file.path(paths$tables, "WGCNA_interpretable_summary.xlsx"))
  }

  plot_super <- super_join |> dplyr::filter(!is.na(.data$p_value), .data$spatial_unit != "global_spatial_adjusted")
  if (nrow(plot_super)) {
    plot_super <- plot_super |>
      dplyr::mutate(sig = dplyr::case_when(.data$FDR_global < 0.05 ~ "**", .data$FDR_global < 0.10 ~ "*", TRUE ~ ""))
    p <- ggplot2::ggplot(plot_super, ggplot2::aes(x = .data$contrast, y = dplyr::coalesce(.data$Supermodule_FinalLabel, .data$supermodule_id), fill = .data$estimate)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = .data$sig), size = 2.2) +
      ggplot2::facet_grid(effect_scope ~ spatial_unit, scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "grey96", high = "#C84C5A") +
      ggplot2::labs(x = NULL, y = NULL, fill = "Estimate") +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    ggplot2::ggsave(file.path(paths$figures, "interpretable_supermodule_effect_heatmap.svg"), p, width = 170, height = 120, units = "mm", device = svglite::svglite)
    p2 <- ggplot2::ggplot(top_super, ggplot2::aes(x = .data$estimate, y = reorder(.data$supermodule_id, abs(.data$estimate)), color = .data$contrast, size = -log10(.data$p_value))) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Estimate", y = NULL, color = "Contrast", size = "-log10 P") +
      ggplot2::theme_classic(base_size = 8)
    ggplot2::ggsave(file.path(paths$figures, "top_changed_supermodules_dotplot.svg"), p2, width = 130, height = 100, units = "mm", device = svglite::svglite)
  }
  plot_module <- module_join |> dplyr::filter(!is.na(.data$p_value), !is.na(.data$supermodule_label))
  if (nrow(plot_module)) {
    p3 <- ggplot2::ggplot(plot_module, ggplot2::aes(x = .data$contrast, y = .data$module_id, color = .data$estimate, size = -log10(.data$p_value))) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::facet_wrap(~ supermodule_label, scales = "free_y") +
      ggplot2::scale_color_gradient2(low = "#3B6FB6", mid = "grey92", high = "#C84C5A") +
      ggplot2::labs(x = NULL, y = NULL, color = "Estimate", size = "-log10 P") +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    ggplot2::ggsave(file.path(paths$figures, "module_effects_by_supermodule.svg"), p3, width = 170, height = 130, units = "mm", device = svglite::svglite)
  }

  if (ds == "microglia" && nrow(super_join)) {
    micro_modules <- module_annot
    effect_panel <- super_join |>
      dplyr::filter(!is.na(.data$p_value)) |>
      dplyr::mutate(panel = "A group effect", x = .data$contrast, y = .data$supermodule_id, value = .data$estimate)
    comp_panel <- super_annot |>
      dplyr::select("SupermoduleID", dplyr::starts_with("fraction_modules_")) |>
      tidyr::pivot_longer(cols = -dplyr::all_of("SupermoduleID"), names_to = "class", values_to = "value") |>
      dplyr::mutate(panel = "B microenvironment composition", x = .data$class, y = .data$SupermoduleID)
    x_col <- if ("empirical_neuropil_enriched_fraction" %in% names(micro_modules)) "empirical_neuropil_enriched_fraction" else if ("canonical_neuronal_synaptic_neuropil_fraction" %in% names(micro_modules)) "canonical_neuronal_synaptic_neuropil_fraction" else "neuropil_synaptic_neuronal_marker_fraction"
    y_col <- if ("empirical_microglia_roi_enriched_fraction" %in% names(micro_modules)) "empirical_microglia_roi_enriched_fraction" else if ("canonical_microglia_homeostatic_fraction" %in% names(micro_modules)) "canonical_microglia_homeostatic_fraction" else "microglia_marker_fraction"
    module_panel <- micro_modules |>
      dplyr::mutate(panel = "C module marker evidence", x = as.character(round(.data[[x_col]], 2)), y = as.character(round(.data[[y_col]], 2)), value = as.numeric(.data[[y_col]]))
    panel_df <- dplyr::bind_rows(
      effect_panel[, c("panel", "x", "y", "value")],
      comp_panel[, c("panel", "x", "y", "value")],
      module_panel[, c("panel", "x", "y", "value")]
    )
    p_micro <- ggplot2::ggplot(panel_df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.85) +
      ggplot2::facet_wrap(~ panel, scales = "free") +
      ggplot2::scale_color_gradient2(low = "#3B6FB6", mid = "grey85", high = "#C84C5A", na.value = "grey60") +
      ggplot2::labs(x = NULL, y = NULL, color = "Value") +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "bottom")
    ggplot2::ggsave(file.path(paths$figures, "microglia_supermodule_interpretation_panel.svg"), p_micro, width = 190, height = 120, units = "mm", device = svglite::svglite)
  }

  write_run_manifest(
    file.path(paths$logs, "run_manifest.yml"),
    inputs = list(module_effects = "group_effects/module_group_effects.csv", supermodule_effects = "group_effects/supermodule_group_effects.csv"),
    outputs = list(tables = paths$tables, source_data = paths$source_data, figures = paths$figures),
    parameters = list(dataset = ds),
    notes = "Interpretable layer; cross-dataset comparisons should remain at broad supermodule/program label level."
  )
  list(super = super_join, module = module_join, top = top_super)
}

datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
summaries <- lapply(datasets, make_dataset_summary)
names(summaries) <- datasets

if (DATASET_ARG == "all") {
  paths_all <- wgcna_downstream_paths("interpretable_summary", "all")
  cross <- dplyr::bind_rows(lapply(summaries, `[[`, "super")) |>
    dplyr::mutate(program_label = dplyr::coalesce(.data$Supermodule_FinalLabel, .data$supermodule_label, .data$supermodule_id)) |>
    dplyr::group_by(.data$dataset, .data$program_label, .data$contrast, .data$spatial_unit) |>
    dplyr::summarise(
      min_FDR_global = min(.data$FDR_global, na.rm = TRUE),
      max_abs_estimate = max(abs(.data$estimate), na.rm = TRUE),
      representative_direction = .data$direction[which.max(abs(.data$estimate))],
      interpretation_sentence = .data$interpretation_sentence[which.max(abs(.data$estimate))],
      .groups = "drop"
    )
  write_table_and_source(cross, paths_all$tables, paths_all$source_data, "WGCNA_cross_dataset_supermodule_program_summary.csv")
  if (nrow(cross)) {
    p <- ggplot2::ggplot(cross, ggplot2::aes(x = .data$dataset, y = .data$program_label, fill = log10(.data$min_FDR_global))) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ contrast) +
      ggplot2::labs(x = NULL, y = NULL, fill = "log10 FDR") +
      ggplot2::theme_classic(base_size = 8)
    ggplot2::ggsave(file.path(paths_all$figures, "cross_dataset_supermodule_program_effect_heatmap.svg"), p, width = 170, height = 120, units = "mm", device = svglite::svglite)
  }
}

message("WGCNA interpretable summary complete for dataset argument: ", DATASET_ARG)
