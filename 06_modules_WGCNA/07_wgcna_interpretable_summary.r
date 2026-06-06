#!/usr/bin/env Rscript
#
# Combine WGCNA group effects and biological annotation into interpretable summaries.
#
# Design intent:
#   - keep the primary inference from 05_module_supermodule_group_effects.r
#   - join biological annotation from 06_annotate_module_microenvironment.r
#   - produce readable, publication-ready summary figures
#   - avoid overloading one figure with spatial-unit, contrast, dataset, and annotation dimensions at once

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

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

dataset_label <- function(ds) {
  switch(ds,
    neuron_neuropil = "Neuron neuropil",
    neuron_soma = "Neuron soma",
    microglia = "Microglia ROI",
    ds
  )
}

spatial_unit_label <- function(ds) {
  if (identical(ds, "neuron_neuropil")) "region-layer" else "region"
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

first_existing_col <- function(df, candidates, fallback = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else fallback
}

coalesce_chr <- function(...) {
  x <- dplyr::coalesce(...)
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Unlabelled"
  x
}

label_wrap <- function(x, width = 36) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

clip_p <- function(p) pmax(safe_num(p), 1e-300)

col_or_na <- function(df, nm) {
  if (nm %in% names(df)) return(df[[nm]])
  rep(NA_character_, nrow(df))
}

ensure_columns <- function(df, cols) {
  for (nm in cols) {
    if (!nm %in% names(df)) df[[nm]] <- NA
  }
  df
}

effect_sentence <- function(row, level = "supermodule") {
  row_get <- function(name, default = NA) {
    if (name %in% names(row)) row[[name]][[1]] else default
  }

  ds_label <- dataset_label(row$dataset)
  id <- if (level == "supermodule") row_get("supermodule_id") else row_get("module_id")
  label <- if (level == "supermodule") {
    row_get("Supermodule_DisplayLabel", row_get("Supermodule_FinalLabel", row_get("Macroprogram_Display", row_get("supermodule_label", row_get("supermodule_id")))))
  } else {
    row_get("module_label", row_get("ModuleID", row_get("module_id")))
  }
  cls <- row_get("dominant_microenvironment_class", row_get("microenvironment_class"))
  estimate <- safe_num(row_get("estimate"))
  direction <- ifelse(is.na(estimate), "is altered", ifelse(estimate > 0, "is higher", "is lower"))
  spatial_unit <- row_get("spatial_unit")
  spatial <- ifelse(
    is.na(spatial_unit) || spatial_unit == "global_spatial_adjusted",
    "after spatial adjustment",
    paste0("in ", spatial_unit)
  )
  effect_scope <- row_get("effect_scope")
  scope <- ifelse(!is.na(effect_scope), paste0(" (", effect_scope, ")"), "")
  base <- paste(ds_label, level, id, "annotated as", label, direction, "for", row_get("contrast"), spatial, scope)
  base <- gsub("[[:space:]]+\\(", " (", base)

  if (row$dataset != "microglia" || is.na(cls)) return(paste0(base, "."))
  if (cls == "microglia_supported") return(paste0(base, "; this ROI signal has microglia-supported ROI evidence."))
  if (cls == "microglia_state_or_activation_supported") return(paste0(base, "; this ROI signal has microglia state/phagolysosomal support."))
  if (cls == "shared_microenvironment") return(paste0(base, "; this ROI signal is shared local microenvironment signal, not purified microglial regulation."))
  if (cls == "neuropil_sensitive") return(paste0(base, "; this ROI signal is neuropil-sensitive; do not interpret as purified microglial regulation."))
  if (cls == "other_cellular_or_vascular_sensitive") return(paste0(base, "; this ROI signal shows other cellular/vascular marker support; interpret cautiously."))
  paste0(base, "; microenvironment support is ambiguous.")
}

add_interpretation_sentences <- function(df, level) {
  if (!nrow(df)) {
    df$interpretation_sentence <- character()
    return(df)
  }
  df$interpretation_sentence <- vapply(seq_len(nrow(df)), function(i) {
    effect_sentence(df[i, , drop = FALSE], level)
  }, character(1))
  df
}

theme_clean <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = base_size)
    )
}

add_plot_metrics <- function(df, fdr_col = "FDR_global", p_col = "p_value") {
  if (!fdr_col %in% names(df)) df[[fdr_col]] <- NA_real_
  if (!p_col %in% names(df)) df[[p_col]] <- NA_real_
  if (!"estimate" %in% names(df)) df$estimate <- NA_real_

  df |>
    dplyr::mutate(
      estimate = safe_num(.data$estimate),
      p_value = safe_num(.data[[p_col]]),
      FDR_global = safe_num(.data[[fdr_col]]),
      neg_log10_FDR = -log10(clip_p(.data$FDR_global)),
      neg_log10_P = -log10(clip_p(.data$p_value)),
      signed_FDR_score = sign(.data$estimate) * .data$neg_log10_FDR,
      signed_P_score = sign(.data$estimate) * .data$neg_log10_P,
      sig_label = dplyr::case_when(
        !is.na(.data$FDR_global) & .data$FDR_global < 0.05 ~ "**",
        !is.na(.data$FDR_global) & .data$FDR_global < 0.10 ~ "*",
        TRUE ~ ""
      ),
      evidence_rank = dplyr::case_when(
        .data$evidence_status == "robust_FDR" ~ 1L,
        .data$evidence_status == "suggestive_FDR10" ~ 2L,
        .data$evidence_status == "nominal_only" ~ 3L,
        .data$evidence_status == "model_unstable" ~ 4L,
        TRUE ~ 5L
      )
    )
}

main_effect_rows <- function(df) {
  if (!nrow(df)) return(df)
  out <- df |> dplyr::filter(!is.na(.data$p_value))
  if ("spatial_unit" %in% names(out) && any(out$spatial_unit == "global_spatial_adjusted", na.rm = TRUE)) {
    out <- out |> dplyr::filter(.data$spatial_unit == "global_spatial_adjusted")
  }
  out
}

spatial_effect_rows <- function(df) {
  if (!nrow(df) || !"spatial_unit" %in% names(df)) return(df[0, , drop = FALSE])
  df |>
    dplyr::filter(!is.na(.data$p_value), !is.na(.data$spatial_unit), .data$spatial_unit != "global_spatial_adjusted")
}

program_label_col <- function(df, level = "supermodule") {
  if (level == "supermodule") {
    coalesce_chr(
      col_or_na(df, "Supermodule_DisplayLabel"),
      col_or_na(df, "Supermodule_FinalLabel"),
      col_or_na(df, "Supermodule_FinalLabel.y"),
      col_or_na(df, "Supermodule_FinalLabel.x"),
      col_or_na(df, "Macroprogram_Display"),
      col_or_na(df, "supermodule_label_for_module"),
      col_or_na(df, "module_supermodule_label"),
      col_or_na(df, "supermodule_label"),
      col_or_na(df, "supermodule_label.y"),
      col_or_na(df, "supermodule_label.x"),
      col_or_na(df, "SupermoduleID"),
      col_or_na(df, "supermodule_id_for_module"),
      col_or_na(df, "module_supermodule_id"),
      col_or_na(df, "supermodule_id")
    )
  } else {
    # Display labels should prefer the annotation table from 06_annotate_module_microenvironment.r.
    # The group-effect table from 05 usually only knows the raw WGCNA module label/color.
    coalesce_chr(
      col_or_na(df, "module_display_label"),
      col_or_na(df, "ModuleLabel_Final"),
      col_or_na(df, "module_label.y"),
      col_or_na(df, "module_label"),
      col_or_na(df, "endpoint_label"),
      col_or_na(df, "ModuleID"),
      col_or_na(df, "module_id"),
      col_or_na(df, "ModuleColor")
    )
  }
}

module_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\s*[—|-].*$", "", x)  # strip display suffix if present
  x <- sub("^ME", "", x, ignore.case = FALSE)
  x <- sub("^WGCNA[_-]?", "", x, ignore.case = TRUE)
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

module_join_key <- function(df) {
  # Raw IDs/colors only. Do not use module_display_label here because it may contain biology text.
  candidates <- c(
    "module_id", "ModuleID", "ModuleColor", "module_eigengene",
    "endpoint_id", "module_label.x", "module_label", "module_label.y"
  )
  vals <- rep(NA_character_, nrow(df))
  for (nm in candidates) {
    if (nm %in% names(df)) vals <- dplyr::coalesce(vals, as.character(df[[nm]]))
  }
  module_key(vals)
}

pretty_module_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^WGCNA[_-]?", "", x, ignore.case = TRUE)
  x <- sub("^ME", "", x, ignore.case = FALSE)
  x[is.na(x) | !nzchar(x)] <- "unlabelled"
  x
}

pretty_program_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- "Unassigned"
  x
}

candidate_cols <- function(df, candidates) intersect(candidates, names(df))

build_module_supermodule_map <- function(module_effects, module_to_supermodule_map = NULL, supermodule_composition = NULL, super_annot = NULL) {
  out <- data.frame(
    module_key = character(),
    supermodule_id = character(),
    supermodule_label = character(),
    map_source = character(),
    stringsAsFactors = FALSE
  )

  add_map_rows <- function(keys, sid, slabel, source) {
    keys <- unique(module_key(keys))
    keys <- keys[nzchar(keys) & !is.na(keys)]
    if (!length(keys)) return(NULL)
    data.frame(
      module_key = keys,
      supermodule_id = as.character(sid %||% NA_character_),
      supermodule_label = as.character(slabel %||% sid %||% NA_character_),
      map_source = source,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(module_to_supermodule_map) && nrow(module_to_supermodule_map)) {
    map_df <- as.data.frame(module_to_supermodule_map, stringsAsFactors = FALSE)

    module_cols <- candidate_cols(map_df, c(
      "ModuleID", "module_id", "module", "Module", "ModuleColor", "module_color",
      "module_eigengene", "moduleEigengene", "module_eigengene_col"
    ))
    sid_cols <- candidate_cols(map_df, c(
      "Supermodule_DataDrivenID", "Supermodule_DataDriven", "SupermoduleID",
      "supermodule_id", "Supermodule", "Supermodule_DataDrivenLabel"
    ))
    slabel_cols <- candidate_cols(map_df, c(
      "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleLabel", "Supermodule_DataDrivenLabel",
      "supermodule_label", "Supermodule", "SupermoduleID"
    ))

    if (length(module_cols) && length(sid_cols)) {
      rows <- lapply(seq_len(nrow(map_df)), function(i) {
        keys <- unlist(map_df[i, module_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- unlist(map_df[i, sid_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- as.character(sid_vals[!is.na(sid_vals) & nzchar(as.character(sid_vals))])
        sid <- if (length(sid_vals)) sid_vals[[1]] else NA_character_
        if (length(slabel_cols)) {
          lab_vals <- unlist(map_df[i, slabel_cols, drop = TRUE], use.names = FALSE)
          lab_vals <- as.character(lab_vals[!is.na(lab_vals) & nzchar(as.character(lab_vals))])
          slabel <- if (length(lab_vals)) lab_vals[[1]] else sid
        } else {
          slabel <- sid
        }
        add_map_rows(keys, sid, slabel, "module_to_supermodule_map")
      })
      out <- dplyr::bind_rows(out, dplyr::bind_rows(rows))
    }
  }

  if (!is.null(supermodule_composition) && nrow(supermodule_composition) && "member_modules" %in% names(supermodule_composition)) {
    comp <- as.data.frame(supermodule_composition, stringsAsFactors = FALSE)
    for (nm in c("supermodule_id", "supermodule_label", "SupermoduleID", "SupermoduleLabel")) if (!nm %in% names(comp)) comp[[nm]] <- NA_character_
    rows <- lapply(seq_len(nrow(comp)), function(i) {
      members <- unlist(strsplit(as.character(comp$member_modules[[i]]), "[;|, ]+"), use.names = FALSE)
      sid <- dplyr::coalesce(as.character(comp$supermodule_id[[i]]), as.character(comp$SupermoduleID[[i]]))
      slabel <- dplyr::coalesce(as.character(comp$supermodule_label[[i]]), as.character(comp$SupermoduleLabel[[i]]), sid)
      add_map_rows(members, sid, slabel, "supermodule_composition")
    })
    out <- dplyr::bind_rows(out, dplyr::bind_rows(rows))
  }

  if (!nrow(out)) {
    return(data.frame(module_key = character(), module_supermodule_id = character(), module_supermodule_label = character(), module_supermodule_map_source = character()))
  }

  if (!is.null(super_annot) && nrow(super_annot)) {
    ann <- as.data.frame(super_annot, stringsAsFactors = FALSE)
    sid_cols <- candidate_cols(ann, c("SupermoduleID", "supermodule_id", "Supermodule_DataDrivenID", "Supermodule_DataDriven", "Supermodule"))
    label_cols <- candidate_cols(ann, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleLabel", "Supermodule_DataDrivenLabel", "supermodule_label"))
    if (length(sid_cols) && length(label_cols)) {
      ann2 <- dplyr::bind_rows(lapply(seq_len(nrow(ann)), function(i) {
        sid_vals <- unlist(ann[i, sid_cols, drop = TRUE], use.names = FALSE)
        lab_vals <- unlist(ann[i, label_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- as.character(sid_vals[!is.na(sid_vals) & nzchar(as.character(sid_vals))])
        lab_vals <- as.character(lab_vals[!is.na(lab_vals) & nzchar(as.character(lab_vals))])
        data.frame(
          supermodule_id = if (length(sid_vals)) sid_vals[[1]] else NA_character_,
          curated_supermodule_label = if (length(lab_vals)) lab_vals[[1]] else NA_character_,
          stringsAsFactors = FALSE
        )
      })) |>
        dplyr::filter(!is.na(.data$supermodule_id), nzchar(.data$supermodule_id)) |>
        dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
      out <- out |>
        dplyr::left_join(ann2, by = "supermodule_id") |>
        dplyr::mutate(supermodule_label = dplyr::coalesce(.data$curated_supermodule_label, .data$supermodule_label)) |>
        dplyr::select(-dplyr::any_of("curated_supermodule_label"))
    }
  }

  out |>
    dplyr::filter(!is.na(.data$module_key), nzchar(.data$module_key)) |>
    dplyr::filter(!is.na(.data$supermodule_id), nzchar(.data$supermodule_id)) |>
    dplyr::distinct(.data$module_key, .keep_all = TRUE) |>
    dplyr::rename(
      module_supermodule_id = "supermodule_id",
      module_supermodule_label = "supermodule_label",
      module_supermodule_map_source = "map_source"
    )
}

attach_module_supermodules <- function(module_join, module_super_map) {
  if (!nrow(module_join)) return(module_join)
  module_join$module_key <- module_join_key(module_join)
  if (!is.null(module_super_map) && nrow(module_super_map)) {
    module_join <- module_join |>
      dplyr::left_join(module_super_map, by = "module_key")
  } else {
    module_join$module_supermodule_id <- NA_character_
    module_join$module_supermodule_label <- NA_character_
    module_join$module_supermodule_map_source <- NA_character_
  }
  module_join |>
    dplyr::mutate(
      supermodule_id_for_module = dplyr::coalesce(.data$module_supermodule_id, .data$supermodule_id),
      supermodule_label_for_module = dplyr::coalesce(.data$module_supermodule_label, .data$supermodule_label, "Unassigned"),
      supermodule_map_source_for_module = dplyr::coalesce(.data$module_supermodule_map_source, "not_mapped")
    )
}

save_svg <- function(plot, filename, width = 170, height = 110) {
  ggplot2::ggsave(filename, plot, width = width, height = height, units = "mm", device = svglite::svglite)
}

plot_supermodule_main_heatmap <- function(super_join, paths, ds) {
  plot_df <- main_effect_rows(super_join) |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(program_label_col(plot_df, "supermodule"), 34)
  plot_df <- plot_df |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = unique(.data$contrast)),
      program_label = stats::reorder(.data$program_label, .data$estimate, FUN = function(z) max(abs(z), na.rm = TRUE))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$program_label, fill = .data$signed_FDR_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.4, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "grey96", high = "#C84C5A",
      midpoint = 0, name = "signed\n-log10 FDR",
      na.value = "grey92"
    ) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": supermodule group effects"),
      subtitle = "Primary spatial-adjusted model; asterisks mark FDR < 0.10 (*) and FDR < 0.05 (**)",
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "supermodule_group_effects_main_heatmap.svg"), width = 170, height = 120)
}

plot_supermodule_spatial_heatmap <- function(super_join, paths, ds) {
  plot_df <- spatial_effect_rows(super_join) |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(program_label_col(plot_df, "supermodule"), 30)
  plot_df <- plot_df |>
    dplyr::mutate(
      spatial_unit = factor(.data$spatial_unit, levels = unique(.data$spatial_unit)),
      contrast = factor(.data$contrast, levels = unique(.data$contrast)),
      panel = paste(.data$contrast, .data$spatial_unit, sep = " | ")
    )

  # Keep the spatial plot readable by showing strongest rows when the matrix is too large.
  max_rows <- 45L
  keep_labels <- plot_df |>
    dplyr::group_by(.data$program_label) |>
    dplyr::summarise(rank_metric = min(.data$FDR_global, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$rank_metric) |>
    dplyr::slice_head(n = max_rows) |>
    dplyr::pull(.data$program_label)

  plot_df <- plot_df |> dplyr::filter(.data$program_label %in% keep_labels)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$program_label, fill = .data$signed_FDR_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.20) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.0, color = "black") +
    ggplot2::facet_wrap(~ contrast, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "grey96", high = "#C84C5A",
      midpoint = 0, name = "signed\n-log10 FDR",
      na.value = "grey92"
    ) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": spatially resolved supermodule effects"),
      subtitle = paste0("Spatial rows use ", spatial_unit_label(ds), "-specific tests; strongest programs shown if crowded"),
      x = spatial_unit_label(ds), y = NULL
    ) +
    theme_clean(base_size = 7.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "supermodule_group_effects_spatial_heatmap.svg"), width = 190, height = 130)
}


plot_module_main_heatmap <- function(module_join, paths, ds) {
  plot_df <- main_effect_rows(module_join) |>
    add_plot_metrics() |>
    dplyr::filter(!is.na(.data$p_value))

  if (!nrow(plot_df)) return(invisible(NULL))

  raw_id <- pretty_module_label(coalesce_chr(col_or_na(plot_df, "module_id"), col_or_na(plot_df, "ModuleID"), col_or_na(plot_df, "ModuleColor")))
  ann_lab <- shorten_supermodule_label(coalesce_chr(col_or_na(plot_df, "microenvironment_label"), program_label_col(plot_df, "module")), max_chars = 28)
  plot_df$module_label_plot <- paste0(raw_id, " | ", ann_lab)
  plot_df$module_label_plot <- label_wrap(plot_df$module_label_plot, 26)
  plot_df$supermodule_plot <- label_wrap(shorten_supermodule_label(pretty_program_label(plot_df$supermodule_label_for_module), max_chars = 45), 30)

  max_modules <- 45L
  keep_modules <- plot_df |>
    dplyr::group_by(.data$module_label_plot) |>
    dplyr::summarise(best_p = min(.data$p_value, na.rm = TRUE), best_abs_estimate = max(abs(.data$estimate), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_p, dplyr::desc(.data$best_abs_estimate)) |>
    dplyr::slice_head(n = max_modules) |>
    dplyr::pull(.data$module_label_plot)

  plot_df <- plot_df |>
    dplyr::filter(.data$module_label_plot %in% keep_modules) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "RES - SUS")),
      module_label_plot = stats::reorder(.data$module_label_plot, .data$estimate, FUN = function(z) max(abs(z), na.rm = TRUE)),
      evidence_label = dplyr::case_when(
        .data$FDR_global < 0.05 ~ "FDR < 0.05",
        .data$FDR_global < 0.10 ~ "FDR < 0.10",
        .data$p_value < 0.05 ~ "P < 0.05",
        TRUE ~ "n.s."
      )
    )
  write_csv_safe2(plot_df, file.path(paths$source_data, "module_group_effects_main_dotplot_source.csv"))
  fig_h <- max(120, 45 + 6.0 * length(unique(plot_df$module_label_plot)) + 5.0 * length(unique(plot_df$supermodule_plot)))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$module_label_plot)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$estimate, size = .data$neg_log10_P, alpha = .data$evidence_label),
      shape = 21, color = "grey20", stroke = 0.20
    ) +
    ggplot2::facet_grid(supermodule_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "grey96", high = "#C84C5A",
      midpoint = 0, name = "Estimate"
    ) +
    ggplot2::scale_size_continuous(range = c(1.4, 4.2), name = "-log10 P") +
    ggplot2::scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR < 0.10" = 0.9, "P < 0.05" = 0.70, "n.s." = 0.30), name = "Evidence") +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module-level group effects"),
      subtitle = paste0("One row per WGCNA module; showing strongest ", min(max_modules, length(unique(plot_df$module_label_plot))), " modules by P value/effect size"),
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 7.6) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "module_group_effects_main_dotplot.svg"), width = 185, height = fig_h)

  # Also write a heatmap with effect-size fill for readers who prefer a matrix view.
  p_heat <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$module_label_plot, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.0, color = "black") +
    ggplot2::facet_grid(supermodule_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "grey96", high = "#C84C5A", midpoint = 0, name = "Estimate") +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module-level group effects"),
      subtitle = paste0("Effect-size heatmap; showing strongest ", min(max_modules, length(unique(plot_df$module_label_plot))), " modules; asterisks mark FDR < 0.10 (*) and FDR < 0.05 (**)"),
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 7.4) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      legend.position = "right"
    )

  save_svg(p_heat, file.path(paths$figures, "module_group_effects_main_heatmap.svg"), width = 185, height = fig_h)
}

plot_module_spatial_heatmap <- function(module_join, paths, ds) {
  plot_df <- spatial_effect_rows(module_join) |>
    add_plot_metrics() |>
    dplyr::filter(!is.na(.data$p_value))

  if (!nrow(plot_df)) return(invisible(NULL))

  raw_id <- pretty_module_label(coalesce_chr(col_or_na(plot_df, "module_id"), col_or_na(plot_df, "ModuleID"), col_or_na(plot_df, "ModuleColor")))
  ann_lab <- shorten_supermodule_label(coalesce_chr(col_or_na(plot_df, "microenvironment_label"), program_label_col(plot_df, "module")), max_chars = 26)
  plot_df$module_label_plot <- label_wrap(paste0(raw_id, " | ", ann_lab), 25)
  plot_df$supermodule_plot <- label_wrap(shorten_supermodule_label(pretty_program_label(plot_df$supermodule_label_for_module), max_chars = 45), 30)

  max_modules <- 35L
  keep_modules <- plot_df |>
    dplyr::group_by(.data$module_label_plot) |>
    dplyr::summarise(best_p = min(.data$p_value, na.rm = TRUE), best_abs_estimate = max(abs(.data$estimate), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_p, dplyr::desc(.data$best_abs_estimate)) |>
    dplyr::slice_head(n = max_modules) |>
    dplyr::pull(.data$module_label_plot)

  plot_df <- plot_df |>
    dplyr::filter(.data$module_label_plot %in% keep_modules) |>
    dplyr::mutate(
      spatial_unit = factor(.data$spatial_unit, levels = unique(.data$spatial_unit)),
      contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "RES - SUS"))
    )
  write_csv_safe2(plot_df, file.path(paths$source_data, "module_group_effects_spatial_heatmap_source.csv"))
  fig_h <- max(135, 50 + 6.5 * length(unique(plot_df$module_label_plot)) + 5.0 * length(unique(plot_df$supermodule_plot)))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$module_label_plot, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.18) +
    ggplot2::facet_grid(supermodule_plot ~ contrast, scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "grey96", high = "#C84C5A",
      midpoint = 0, name = "Estimate",
      na.value = "grey92"
    ) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": spatially resolved module effects"),
      subtitle = paste0("Effect estimates across ", spatial_unit_label(ds), " tests; showing strongest ", min(max_modules, length(unique(plot_df$module_label_plot))), " modules"),
      x = spatial_unit_label(ds), y = NULL
    ) +
    theme_clean(base_size = 6.7) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "module_group_effects_spatial_heatmap.svg"), width = 210, height = fig_h)
}

plot_supermodule_membership_overview <- function(module_join, super_join, paths, ds) {
  if (!nrow(module_join)) return(invisible(NULL))

  df <- module_join |>
    dplyr::distinct(.data$module_id, .data$supermodule_id_for_module, .data$supermodule_label_for_module, .keep_all = TRUE)

  df$supermodule_plot <- label_wrap(pretty_program_label(df$supermodule_label_for_module), 30)

  overview <- df |>
    dplyr::count(.data$supermodule_plot, name = "n_modules") |>
    dplyr::arrange(dplyr::desc(.data$n_modules))

  if (!nrow(overview)) return(invisible(NULL))

  p <- ggplot2::ggplot(overview, ggplot2::aes(x = .data$n_modules, y = stats::reorder(.data$supermodule_plot, .data$n_modules))) +
    ggplot2::geom_col(width = 0.64, fill = "grey45") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_modules), hjust = -0.18, size = 2.5) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module-to-supermodule compression"),
      subtitle = "Number of WGCNA modules assigned to each supermodule; not an effect-size plot",
      x = "Modules per supermodule", y = NULL
    ) +
    theme_clean(base_size = 8)

  save_svg(p, file.path(paths$figures, "supermodule_membership_overview.svg"), width = 150, height = 95)
}

plot_top_supermodules <- function(top_super, paths, ds) {
  plot_df <- top_super |>
    add_plot_metrics() |>
    dplyr::filter(!is.na(.data$p_value))

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(program_label_col(plot_df, "supermodule"), 36)
  plot_df <- plot_df |>
    dplyr::mutate(
      program_label = stats::reorder(.data$program_label, abs(.data$estimate))
    ) |>
    dplyr::slice_head(n = 30)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$estimate, y = .data$program_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
    ggplot2::geom_point(ggplot2::aes(size = .data$neg_log10_P, fill = .data$contrast), shape = 21, color = "black", stroke = 0.20, alpha = 0.90) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": strongest supermodule contrasts"),
      subtitle = "Point size reflects nominal evidence; use tables for exact adjusted evidence",
      x = "Model estimate", y = NULL,
      size = "-log10 P", fill = "Contrast"
    ) +
    theme_clean(base_size = 8) +
    ggplot2::theme(legend.position = "right")

  save_svg(p, file.path(paths$figures, "top_supermodule_effects_dotplot.svg"), width = 150, height = 115)
}

plot_module_effects_by_supermodule <- function(module_join, paths, ds) {
  plot_df <- module_join |>
    dplyr::filter(!is.na(.data$p_value), !is.na(.data$supermodule_label)) |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$module_label_plot <- program_label_col(plot_df, "module")
  plot_df$supermodule_label_plot <- label_wrap(as.character(plot_df$supermodule_label), 24)

  # This is a supplementary diagnostic figure. Limit each facet to the strongest module rows.
  plot_df <- plot_df |>
    dplyr::group_by(.data$supermodule_label_plot, .data$module_id) |>
    dplyr::mutate(module_rank = min(.data$p_value, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$supermodule_label_plot) |>
    dplyr::arrange(.data$module_rank, .by_group = TRUE) |>
    dplyr::filter(dplyr::row_number() <= 30L) |>
    dplyr::ungroup()

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$module_id, fill = .data$estimate, size = .data$neg_log10_P)) +
    ggplot2::geom_point(shape = 21, color = "black", stroke = 0.15, alpha = 0.85) +
    ggplot2::facet_wrap(~ supermodule_label_plot, scales = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "grey96", high = "#C84C5A", midpoint = 0, name = "Estimate") +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module effects within supermodules"),
      subtitle = "Supplementary view; strongest modules per supermodule shown",
      x = NULL, y = NULL, size = "-log10 P"
    ) +
    theme_clean(base_size = 7.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  save_svg(p, file.path(paths$figures, "module_effects_by_supermodule_dotplot.svg"), width = 180, height = 140)
}

plot_microglia_effect_heatmap <- function(super_join, paths) {
  plot_df <- main_effect_rows(super_join) |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(program_label_col(plot_df, "supermodule"), 34)
  if (!"dominant_microenvironment_class" %in% names(plot_df)) plot_df$dominant_microenvironment_class <- NA_character_
  plot_df <- plot_df |>
    dplyr::mutate(
      microenvironment_class = dplyr::coalesce(.data$dominant_microenvironment_class, "unclassified")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$program_label, fill = .data$signed_FDR_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.4, color = "black") +
    ggplot2::facet_grid(microenvironment_class ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "grey96", high = "#C84C5A", midpoint = 0, name = "signed\n-log10 FDR") +
    ggplot2::labs(
      title = "Microglia ROI: supermodule effects by microenvironment class",
      subtitle = "Facet labels distinguish microglia-supported from shared/vascular/neuropil-sensitive ROI signals",
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 7.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  save_svg(p, file.path(paths$figures, "microglia_supermodule_effects_by_microenvironment.svg"), width = 175, height = 130)
}

plot_microglia_composition_heatmap <- function(super_annot, paths) {
  frac_cols <- grep("^fraction_modules_", names(super_annot), value = TRUE)
  if (!length(frac_cols) || !nrow(super_annot)) return(invisible(NULL))

  plot_df <- super_annot |>
    dplyr::mutate(
      program_label = coalesce_chr(
        col_or_na(super_annot, "Supermodule_DisplayLabel"),
        col_or_na(super_annot, "Supermodule_FinalLabel"),
        col_or_na(super_annot, "Macroprogram_Display"),
        col_or_na(super_annot, "SupermoduleID")
      ),
      program_label = label_wrap(shorten_supermodule_label(.data$program_label, max_chars = 45), 34)
    ) |>
    dplyr::select(dplyr::any_of(c("SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display")), "program_label", dplyr::all_of(frac_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(frac_cols), names_to = "class", values_to = "fraction") |>
    dplyr::mutate(
      fraction = safe_num(.data$fraction),
      class = gsub("^fraction_modules_", "", .data$class),
      class = gsub("_", " ", .data$class)
    ) |>
    dplyr::filter(is.finite(.data$fraction), .data$fraction > 0)

  if (!nrow(plot_df)) return(invisible(NULL))
  write_csv_safe2(plot_df, file.path(paths$source_data, "microglia_supermodule_annotation_composition_source.csv"))
  fig_h <- max(105, 16 * length(unique(plot_df$program_label)) + 40)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$fraction, y = stats::reorder(.data$program_label, .data$SupermoduleID), fill = .data$class)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(
      title = "Microglia ROI: supermodule annotation composition",
      subtitle = "Fractions summarize the module classes contributing to each supermodule",
      x = "Fraction of member modules", y = NULL, fill = NULL
    ) +
    theme_clean(base_size = 7.5) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 7.2)
    )

  save_svg(p, file.path(paths$figures, "microglia_supermodule_annotation_composition.svg"), width = 180, height = fig_h)
}

plot_microglia_marker_evidence <- function(module_annot, paths) {
  if (!nrow(module_annot)) return(invisible(NULL))

  x_col <- first_existing_col(
    module_annot,
    c("empirical_neuropil_enriched_fraction", "canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction")
  )
  y_col <- first_existing_col(
    module_annot,
    c("empirical_microglia_roi_enriched_fraction", "canonical_microglia_homeostatic_fraction", "microglia_marker_fraction")
  )
  if (is.null(x_col) || is.null(y_col)) return(invisible(NULL))

  class_col <- first_existing_col(module_annot, c("microenvironment_class", "dominant_microenvironment_class"), fallback = NULL)

  plot_df <- module_annot |>
    dplyr::mutate(
      x_value = safe_num(.data[[x_col]]),
      y_value = safe_num(.data[[y_col]])
    )
  plot_df$module_label_plot <- program_label_col(plot_df, "module")

  if (!is.null(class_col)) {
    plot_df$marker_class <- as.character(plot_df[[class_col]])
  } else {
    plot_df$marker_class <- "unclassified"
  }

  plot_df <- plot_df |> dplyr::filter(!is.na(.data$x_value), !is.na(.data$y_value))
  if (!nrow(plot_df)) return(invisible(NULL))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x_value, y = .data$y_value)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.25, color = "grey60", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$marker_class), shape = 21, color = "black", stroke = 0.20, size = 2.0, alpha = 0.85) +
    ggplot2::labs(
      title = "Microglia ROI: module marker evidence",
      subtitle = "Numeric scatter; diagonal separates neuropil-skewed from microglia-skewed marker evidence",
      x = x_col,
      y = y_col,
      fill = "Class"
    ) +
    theme_clean(base_size = 8) +
    ggplot2::theme(legend.position = "right")

  save_svg(p, file.path(paths$figures, "microglia_module_marker_evidence_scatter.svg"), width = 145, height = 115)
}

make_dataset_summary <- function(ds) {
  paths <- wgcna_downstream_paths("interpretable_summary", ds)

  # Correct upstream inputs:
  #   - group_effects/<dataset>/module_group_effects.csv
  #   - group_effects/<dataset>/supermodule_group_effects.csv
  #   - module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
  #   - module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv
  #   - 04_wgcna_de_gsea_overlap/<dataset>/WGCNA_vs_DE_GSEA_overlap.csv
  module_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))
  super_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))
  module_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))
  super_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"))
  overlap <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", ds, "WGCNA_vs_DE_GSEA_overlap.csv"))
  super_comp <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_composition.csv"))
  module_super_map_file <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_to_supermodule_map_with_annotations.csv"))

  if (is.null(module_effects)) module_effects <- empty_group_effects(ds, "module", "missing module_group_effects.csv")
  if (is.null(super_effects)) super_effects <- empty_group_effects(ds, "supermodule", "missing supermodule_group_effects.csv")
  if (is.null(module_annot)) module_annot <- data.frame(dataset = ds, ModuleID = character(), microenvironment_class = character())
  if (is.null(super_annot)) super_annot <- data.frame(dataset = ds, SupermoduleID = character(), Supermodule_FinalLabel = character(), dominant_microenvironment_class = character())

  needed_super_cols <- c(
    "Supermodule_DisplayLabel",
    "Supermodule_LongLabel",
    "Macroprogram_Display",
    "dominant_microenvironment_class",
    "fraction_modules_microglia_supported",
    "fraction_modules_microglia_state_or_activation_supported",
    "fraction_modules_shared_microenvironment",
    "fraction_modules_neuropil_sensitive",
    "fraction_modules_other_cellular_or_vascular_sensitive",
    "fraction_modules_ambiguous",
    "marker_registry_version",
    "empirical_marker_set_version",
    "classification_rationale",
    "Supermodule_LabelSource",
    "Supermodule_LabelConfidence",
    "GO_label_confidence_class",
    "annotation_scope",
    "manual_label_status",
    "ManualReviewRequired"
  )
  for (nm in needed_super_cols) {
    if (!nm %in% names(super_annot)) super_annot[[nm]] <- NA
  }

  super_annot <- super_annot |>
    dplyr::arrange(.data$SupermoduleID) |>
    dplyr::distinct(.data$dataset, .data$SupermoduleID, .keep_all = TRUE)

  module_join <- module_effects |>
    dplyr::left_join(module_annot, by = c("dataset" = "dataset", "module_id" = "ModuleID")) |>
    add_interpretation_sentences("module") |>
    ensure_columns(c("module_label", "module_label.x", "module_label.y", "ModuleID", "ModuleColor", "module_id", "module_display_label", "microenvironment_label", "microenvironment_class", "microenvironment_confidence", "supermodule_label", "p_value", "estimate", "FDR_global", "contrast"))

  super_join <- super_effects |>
    dplyr::left_join(super_annot, by = c("dataset" = "dataset", "supermodule_id" = "SupermoduleID")) |>
    add_interpretation_sentences("supermodule") |>
    ensure_columns(c("Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display", "Supermodule_FinalLabel", "supermodule_label", "supermodule_id", "dominant_microenvironment_class", "p_value", "estimate", "FDR_global", "contrast", "spatial_unit", "effect_scope", "evidence_status", "direction"))

  module_super_map <- build_module_supermodule_map(
    module_effects = module_join,
    module_to_supermodule_map = module_super_map_file,
    supermodule_composition = super_comp,
    super_annot = super_annot
  )
  module_join <- attach_module_supermodules(module_join, module_super_map)

  module_supermodule_join_qc <- module_join |>
    dplyr::distinct(.data$module_id, .data$module_key, .data$module_display_label, .data$microenvironment_label, .data$supermodule_label_for_module, .data$supermodule_map_source_for_module)

  write_table_and_source(
    module_supermodule_join_qc,
    paths$tables,
    paths$source_data,
    "WGCNA_module_supermodule_join_qc.csv"
  )

  if (!is.null(overlap) && nrow(overlap) && "ModuleID" %in% names(overlap)) {
    overlap_summary <- overlap |>
      dplyr::group_by(.data$ModuleID) |>
      dplyr::summarise(
        best_wgcna_de_gsea_overlap = paste(utils::head(unique(.data$contrast), 5), collapse = ";"),
        .groups = "drop"
      )
    module_join <- module_join |> dplyr::left_join(overlap_summary, by = c("module_id" = "ModuleID"))
  }

  top_super <- super_join |>
    dplyr::filter(!is.na(.data$p_value)) |>
    add_plot_metrics() |>
    dplyr::arrange(.data$evidence_rank, .data$FDR_global, .data$p_value, dplyr::desc(abs(.data$estimate))) |>
    dplyr::select(-dplyr::any_of(c("neg_log10_FDR", "neg_log10_P", "signed_FDR_score", "signed_P_score", "sig_label", "evidence_rank"))) |>
    dplyr::slice_head(n = 50)

  write_table_and_source(super_join, paths$tables, paths$source_data, "WGCNA_supermodule_group_effects_interpretable.csv")
  write_table_and_source(module_join, paths$tables, paths$source_data, "WGCNA_module_group_effects_interpretable.csv")
  write_table_and_source(top_super, paths$tables, paths$source_data, "WGCNA_top_changed_supermodules.csv")

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(supermodules = super_join, modules = module_join, top_supermodules = top_super),
      file.path(paths$tables, "WGCNA_interpretable_summary.xlsx")
    )
  }

  plot_supermodule_main_heatmap(super_join, paths, ds)
  plot_supermodule_spatial_heatmap(super_join, paths, ds)
  plot_supermodule_membership_overview(module_join, super_join, paths, ds)
  plot_module_main_heatmap(module_join, paths, ds)
  plot_module_spatial_heatmap(module_join, paths, ds)
  plot_top_supermodules(top_super, paths, ds)
  plot_module_effects_by_supermodule(module_join, paths, ds)

  if (identical(ds, "microglia")) {
    plot_microglia_effect_heatmap(super_join, paths)
    plot_microglia_composition_heatmap(super_annot, paths)
    plot_microglia_marker_evidence(module_annot, paths)
  }

  write_run_manifest(
    file.path(paths$logs, "run_manifest.yml"),
    inputs = list(
      module_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
      supermodule_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
      module_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"),
      supermodule_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
      de_gsea_overlap = path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", ds, "WGCNA_vs_DE_GSEA_overlap.csv"),
      supermodule_composition = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_composition.csv"),
      module_to_supermodule_map = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_to_supermodule_map_with_annotations.csv")
    ),
    outputs = list(tables = paths$tables, source_data = paths$source_data, figures = paths$figures),
    parameters = list(dataset = ds),
    notes = "Interpretable layer. Supermodules are treated as compressed overview units; module-level dotplots/heatmaps are exported as the biological-resolution view after explicitly joining the module-to-supermodule map; WGCNA_module_supermodule_join_qc.csv records whether labels came from 06 annotations and whether modules mapped to supermodules. Main heatmaps use global spatial-adjusted rows where available; spatial rows are plotted separately."
  )

  list(super = super_join, module = module_join, top = top_super)
}

make_cross_dataset_summary <- function(summaries) {
  paths_all <- wgcna_downstream_paths("interpretable_summary", "all")

  all_super <- dplyr::bind_rows(lapply(summaries, `[[`, "super"))
  if (!nrow(all_super)) {
    cross <- data.frame()
    write_table_and_source(cross, paths_all$tables, paths_all$source_data, "WGCNA_cross_dataset_supermodule_program_summary.csv")
    return(invisible(cross))
  }

  cross_base <- all_super |>
    add_plot_metrics() |>
    dplyr::mutate(
      program_label = coalesce_chr(
        col_or_na(all_super, "Supermodule_DisplayLabel"),
        col_or_na(all_super, "Supermodule_FinalLabel"),
        col_or_na(all_super, "Macroprogram_Display"),
        col_or_na(all_super, "supermodule_label"),
        col_or_na(all_super, "supermodule_id")
      ),
      is_global = !is.na(.data$spatial_unit) & .data$spatial_unit == "global_spatial_adjusted",
      row_priority = dplyr::case_when(
        .data$is_global ~ 1L,
        TRUE ~ 2L
      )
    ) |>
    dplyr::filter(!is.na(.data$p_value))

  # One row per dataset x broad program x contrast.
  # Prefer global spatial-adjusted rows; if absent, use the strongest available spatial row.
  cross <- cross_base |>
    dplyr::group_by(.data$dataset, .data$program_label, .data$contrast) |>
    dplyr::arrange(.data$row_priority, .data$FDR_global, .data$p_value, dplyr::desc(abs(.data$estimate)), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dataset,
      program_label,
      contrast,
      selected_spatial_unit = .data$spatial_unit,
      selected_effect_scope = .data$effect_scope,
      estimate,
      p_value,
      FDR_global,
      neg_log10_FDR,
      signed_FDR_score,
      evidence_status,
      dominant_microenvironment_class,
      interpretation_sentence
    )

  write_table_and_source(cross, paths_all$tables, paths_all$source_data, "WGCNA_cross_dataset_supermodule_program_summary.csv")

  if (nrow(cross)) {
    plot_df <- cross |>
      dplyr::mutate(
        dataset = factor(.data$dataset, levels = c("microglia", "neuron_soma", "neuron_neuropil")),
        dataset_label = factor(dataset_label(as.character(.data$dataset)), levels = dataset_label(c("microglia", "neuron_soma", "neuron_neuropil"))),
        program_label_plot = label_wrap(.data$program_label, 34),
        sig_label = dplyr::case_when(
          !is.na(.data$FDR_global) & .data$FDR_global < 0.05 ~ "**",
          !is.na(.data$FDR_global) & .data$FDR_global < 0.10 ~ "*",
          TRUE ~ ""
        )
      )

    # Limit crowded cross-dataset plot to programs with at least nominal evidence somewhere.
    keep_programs <- plot_df |>
      dplyr::group_by(.data$program_label_plot) |>
      dplyr::summarise(best_p = min(.data$p_value, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data$best_p) |>
      dplyr::slice_head(n = 40) |>
      dplyr::pull(.data$program_label_plot)

    plot_df <- plot_df |> dplyr::filter(.data$program_label_plot %in% keep_programs)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$dataset_label, y = .data$program_label_plot, fill = .data$signed_FDR_score)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.2, color = "black") +
      ggplot2::facet_wrap(~ contrast, nrow = 1) +
      ggplot2::scale_fill_gradient2(
        low = "#3B6FB6", mid = "grey96", high = "#C84C5A",
        midpoint = 0, name = "signed\n-log10 FDR",
        na.value = "grey92"
      ) +
      ggplot2::labs(
        title = "Cross-dataset supermodule program effects",
        subtitle = "One selected row per dataset, program, and contrast; global spatial-adjusted rows preferred",
        x = NULL, y = NULL
      ) +
      theme_clean(base_size = 7.5) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

    save_svg(p, file.path(paths_all$figures, "cross_dataset_supermodule_program_effect_heatmap.svg"), width = 190, height = 130)
  }

  invisible(cross)
}

if (run$dry_run) {
  datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
  for (ds in datasets) {
    paths <- wgcna_downstream_paths("interpretable_summary", ds)
    invisible(lapply(unlist(paths), dir_create))
    dry_run_line("Dataset", ds)
    dry_run_line(
      "Module effects",
      path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Supermodule effects",
      path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Module annotation",
      path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Supermodule annotation",
      path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"))) "PASS" else "WARN"
    )
  }
  if (DATASET_ARG == "all") {
    dry_run_line("Cross-dataset output", path_results("tables", "06_modules_WGCNA", "interpretable_summary", "all"))
  }
  quit(status = 0, save = "no")
}

datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
summaries <- lapply(datasets, make_dataset_summary)
names(summaries) <- datasets

if (DATASET_ARG == "all") {
  make_cross_dataset_summary(summaries)
}

message("WGCNA interpretable summary complete for dataset argument: ", DATASET_ARG)
