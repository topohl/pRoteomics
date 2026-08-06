#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/09b_microglia_neuropil_independence_figures.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: Existing microglia matched-neuropil independence diagnostic tables.
# Produces: Publication-style figures and source data only.
# Notes: Rendering layer only; no WGCNA recomputation and no model refitting.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "plotting_nature.R"))

SCRIPT_ID <- "06_modules_WGCNA/09b_microglia_neuropil_independence_figures.R"
runtime <- init_script_runtime(SCRIPT_ID, stage = "modules_downstream", default_dataset = "microglia")
if (!identical(runtime$dataset, "microglia") && !isTRUE(runtime$dry_run)) {
  stop("This figure layer is microglia-only. Use --dataset microglia.", call. = FALSE)
}
parse_bool_arg <- function(flag, default = FALSE, args = runtime$args) {
  value <- script_arg_value(flag, if (isTRUE(default)) "TRUE" else "FALSE", args = args)
  value <- tolower(as.character(value))
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Expected TRUE/FALSE for ", flag, ", got: ", value, call. = FALSE)
}
COMPACT <- parse_bool_arg("--compact", TRUE)
INCLUDE_EXPLORATORY <- parse_bool_arg("--include-exploratory", FALSE)
LABEL_ALL <- parse_bool_arg("--label-all", FALSE)

required_pkgs <- c("dplyr", "tidyr", "tibble", "readr", "ggplot2", "scales", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !isTRUE(runtime$dry_run)) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

dataset <- "microglia"
in_dir <- path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", dataset)
paths <- create_module_dirs("06_modules_WGCNA", file.path("microglia_neuropil_independence", dataset))
fig_dir <- paths$figures
source_dir <- paths$source_data
log_dir <- paths$logs

inputs <- list(
  classification = file.path(in_dir, "microglia_module_neuropil_independence_classification.csv"),
  effects = file.path(in_dir, "microglia_neuropil_independence_effects.csv"),
  handoff_with_independence = file.path(in_dir, "WGCNA_module_inferential_handoff_with_neuropil_independence.csv"),
  module_annotation_optional = path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_module_biological_annotation.csv"),
  final_label_lookup_optional = path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset, "WGCNA_final_label_lookup.csv")
)
required_inputs <- inputs[c("classification", "effects", "handoff_with_independence")]

figure_outputs <- c(
  diagnostic_svg = file.path(fig_dir, "microglia_neuropil_diagnostic_screen.svg"),
  diagnostic_pdf = file.path(fig_dir, "microglia_neuropil_diagnostic_screen.pdf"),
  heatmap_svg = file.path(fig_dir, "microglia_neuropil_retained_effect_heatmap.svg"),
  heatmap_pdf = file.path(fig_dir, "microglia_neuropil_retained_effect_heatmap.pdf"),
  dumbbell_svg = file.path(fig_dir, "WGCNA_4D4D4D_neuropil_adjustment_dumbbell.svg"),
  dumbbell_pdf = file.path(fig_dir, "WGCNA_4D4D4D_neuropil_adjustment_dumbbell.pdf")
)
source_outputs <- c(
  diagnostic = file.path(source_dir, "microglia_neuropil_diagnostic_screen_source.csv"),
  heatmap = file.path(source_dir, "microglia_neuropil_retained_effect_heatmap_source.csv"),
  dumbbell = file.path(source_dir, "WGCNA_4D4D4D_neuropil_adjustment_dumbbell_source.csv")
)
manifest_path <- file.path(log_dir, "run_manifest.yml")

if (isTRUE(runtime$dry_run)) {
  dry_run_line("Script", SCRIPT_ID)
  dry_run_line("Dataset", dataset)
  dry_run_line("Compact", COMPACT)
  dry_run_line("Include exploratory", INCLUDE_EXPLORATORY)
  dry_run_line("Label all diagnostic modules", LABEL_ALL)
  for (nm in names(required_inputs)) dry_run_line(nm, required_inputs[[nm]], if (file.exists(required_inputs[[nm]])) "PASS" else "FAIL")
  for (nm in setdiff(names(inputs), names(required_inputs))) dry_run_line(nm, inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "WARN")
  dry_run_line("Figures", fig_dir)
  dry_run_line("Source data", source_dir)
  dry_run_line("Manifest", manifest_path)
  quit(status = if (all(file.exists(unlist(required_inputs, use.names = FALSE)))) 0 else 1, save = "no")
}

missing_required <- required_inputs[!file.exists(unlist(required_inputs, use.names = FALSE))]
if (length(missing_required)) {
  stop(
    "Missing required input table(s) for the matched-neuropil figure layer:\n",
    paste0(" - ", names(missing_required), ": ", unname(missing_required), collapse = "\n"),
    "\nRun 06_modules_WGCNA/09_microglia_neuropil_independence.R first; this script only visualizes its outputs.",
    call. = FALSE
  )
}

read_csv_safe <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
classification <- read_csv_safe(inputs$classification)
effects <- read_csv_safe(inputs$effects)
invisible(read_csv_safe(inputs$handoff_with_independence))

require_columns <- function(dat, cols, label) {
  missing <- setdiff(cols, names(dat))
  if (length(missing)) stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

require_columns(
  classification,
  c("module_id", "endpoint_label", "claim_gate_eligible", "best_contrast", "min_p_before", "max_percent_attenuation"),
  "classification table"
)
require_columns(
  effects,
  c("module_id", "endpoint_type", "contrast", "covariate_family", "effect_before", "effect_after", "direction_preserved"),
  "effects table"
)

first_present <- function(dat, cols) intersect(cols, names(dat))[1]
clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}

label_lookup <- tibble::tibble(module_id = character(), biological_label = character())
if (file.exists(inputs$final_label_lookup_optional)) {
  lookup <- read_csv_safe(inputs$final_label_lookup_optional)
  if (all(c("entity_id", "level") %in% names(lookup))) {
    label_col <- first_present(lookup, c("best_data_driven_label", "final_plot_label", "parent_program", "raw_top_GO_label"))
    if (!is.na(label_col)) {
      label_lookup <- lookup |>
        dplyr::filter(.data$level == "module") |>
        dplyr::transmute(module_id = .data$entity_id, biological_label = .data[[label_col]]) |>
        dplyr::distinct()
    }
  }
}
if (!nrow(label_lookup) && file.exists(inputs$module_annotation_optional)) {
  annot <- read_csv_safe(inputs$module_annotation_optional)
  id_col <- first_present(annot, c("ModuleID", "module_id"))
  label_col <- first_present(annot, c("cleaned_biological_label_short", "module_biological_label_short", "cleaned_biological_label", "module_biological_label", "ModuleLabel_Final", "module_label"))
  if (!is.na(id_col) && !is.na(label_col)) {
    label_lookup <- annot |>
      dplyr::transmute(module_id = .data[[id_col]], biological_label = .data[[label_col]]) |>
      dplyr::distinct()
  }
}

compact_module_id <- function(x) sub("^WGCNA_", "", as.character(x))
strip_module_from_label_one <- function(label, module_id) {
  label <- clean_chr(label)
  id <- as.character(module_id)
  hex <- compact_module_id(id)
  label <- gsub("^WGCNA_#[0-9A-Fa-f]{6}\\s*[|.\u00b7-]+\\s*", "", label)
  label <- ifelse(!is.na(label), trimws(gsub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", id), "\\s*[|.\u00b7-]+\\s*"), "", label)), label)
  label <- ifelse(!is.na(label), trimws(gsub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", hex), "\\s*[|.\u00b7-]+\\s*"), "", label)), label)
  label
}
make_compact_label <- function(module_id, biological_label) {
  mapply(
    function(mid, lab) {
      bio <- strip_module_from_label_one(lab, mid)
      if (is.na(bio) || !nzchar(bio)) bio <- clean_chr(mid)
      paste(compact_module_id(mid), "\u00b7", bio)
    },
    module_id,
    biological_label,
    USE.NAMES = FALSE
  )
}
compact_module_label_lookup <- c(
  `WGCNA_#7F2704` = "#7F2704 \u00b7 Mito/energy",
  `WGCNA_#3F007D` = "#3F007D \u00b7 Mito/energy",
  `WGCNA_#4D4D4D` = "#4D4D4D \u00b7 ECM/adhesion",
  `WGCNA_#A6611A` = "#A6611A \u00b7 Endomembrane",
  `WGCNA_#41AB5D` = "#41AB5D \u00b7 Ion homeostasis",
  `WGCNA_#08519C` = "#08519C \u00b7 Syn/cytosk.",
  `WGCNA_#8C510A` = "#8C510A \u00b7 Proteostasis",
  `WGCNA_#2B8CBE` = "#2B8CBE \u00b7 Vesicle",
  `WGCNA_#969696` = "#969696 \u00b7 Actin/cytosk.",
  `WGCNA_#9E9AC8` = "#9E9AC8 \u00b7 Translation",
  `WGCNA_#737373` = "#737373 \u00b7 Barrier/junction",
  `WGCNA_#1F4E79` = "#1F4E79 \u00b7 Mito gene expr.",
  `WGCNA_#006D2C` = "#006D2C \u00b7 RNA/RNP"
)
compact_module_label <- function(module_id, fallback_label) {
  out <- unname(compact_module_label_lookup[as.character(module_id)])
  missing <- is.na(out) | !nzchar(out)
  out[missing] <- fallback_label[missing]
  out
}

class_for_plot <- classification |>
  dplyr::mutate(
    min_p_before = suppressWarnings(as.numeric(.data$min_p_before)),
    max_percent_attenuation = suppressWarnings(as.numeric(.data$max_percent_attenuation)),
    claim_gate_eligible = as.logical(.data$claim_gate_eligible)
  ) |>
  dplyr::left_join(label_lookup, by = "module_id") |>
  dplyr::mutate(
    biological_label = dplyr::coalesce(clean_chr(.data$biological_label), clean_chr(.data$endpoint_label), clean_chr(.data$module_id)),
    compact_label = make_compact_label(.data$module_id, .data$biological_label),
    compact_plot_label = compact_module_label(.data$module_id, .data$compact_label),
    highlight_4D4D4D = .data$module_id == "WGCNA_#4D4D4D"
  )

diagnostic_source <- class_for_plot |>
  dplyr::filter(is.finite(.data$min_p_before), .data$min_p_before > 0, is.finite(.data$max_percent_attenuation)) |>
  dplyr::mutate(
    neg_log10_min_p_before = -log10(.data$min_p_before),
    diagnostic_stability_score = 100 - .data$max_percent_attenuation,
    figure_text_scope = dplyr::if_else(.data$claim_gate_eligible, "claim_gate_eligible", "matched-neuropil adjustment diagnostic")
  ) |>
  dplyr::select(
    "module_id", "compact_label", "compact_plot_label", "biological_label", "claim_gate_eligible", "best_contrast",
    "min_p_before", "neg_log10_min_p_before", "max_percent_attenuation", "diagnostic_stability_score",
    "neuropil_independence_classification", "primary_effect_status_summary", "figure_text_scope", "highlight_4D4D4D"
  ) |>
  dplyr::arrange(.data$min_p_before, .data$max_percent_attenuation, .data$module_id)

covariate_order <- c(
  "global_neuropil_score",
  "mitochondrial_neuropil_score",
  "ECM_adhesion_neuropil_score",
  "RNA_translation_neuropil_score",
  "synaptic_neuropil_score",
  "cytoskeletal_neuropil_score",
  "exploratory_best_spearman"
)
covariate_labels <- c(
  global_neuropil_score = "Global",
  mitochondrial_neuropil_score = "Mito",
  ECM_adhesion_neuropil_score = "ECM/vascular",
  RNA_translation_neuropil_score = "RNA/transl.",
  synaptic_neuropil_score = "Synaptic",
  cytoskeletal_neuropil_score = "Cytoskeletal",
  exploratory_best_spearman = "Best Spearman"
)
compact_covariate_labels <- c(
  global_neuropil_score = "Global",
  mitochondrial_neuropil_score = "Mito",
  ECM_adhesion_neuropil_score = "ECM",
  RNA_translation_neuropil_score = "RNA",
  synaptic_neuropil_score = "Syn.",
  cytoskeletal_neuropil_score = "Cytosk.",
  exploratory_best_spearman = "Best Spearman"
)
diagnostic_key_labels <- c(
  `WGCNA_#4D4D4D` = "#4D4D4D \u00b7 ECM/adhesion",
  `WGCNA_#3F007D` = "#3F007D \u00b7 mito/energy",
  `WGCNA_#7F2704` = "#7F2704 \u00b7 Acetyl-CoA"
)

row_order <- class_for_plot |>
  dplyr::arrange(.data$min_p_before, .data$max_percent_attenuation, .data$module_id) |>
  dplyr::pull(.data$module_id) |>
  unique()

heatmap_source <- effects |>
  dplyr::filter(.data$endpoint_type == "module_eigengene") |>
  dplyr::inner_join(
    class_for_plot |> dplyr::select("module_id", "best_contrast", "min_p_before", "max_percent_attenuation", "biological_label", "compact_label", "compact_plot_label", "module_claim_gate_eligible" = "claim_gate_eligible"),
    by = "module_id"
  ) |>
  dplyr::filter(.data$contrast == .data$best_contrast, .data$covariate_family %in% covariate_order) |>
  dplyr::mutate(
    effect_before = suppressWarnings(as.numeric(.data$effect_before)),
    effect_after = suppressWarnings(as.numeric(.data$effect_after)),
    direction_preserved = as.logical(.data$direction_preserved)
  ) |>
  dplyr::filter(is.finite(.data$effect_before), is.finite(.data$effect_after)) |>
  dplyr::mutate(
    retained_abs_effect_pct = 100 * abs(.data$effect_after) / abs(.data$effect_before),
    retained_abs_effect_pct_display_capped = pmin(.data$retained_abs_effect_pct, 250),
    covariate_display = unname(covariate_labels[.data$covariate_family]),
    compact_covariate_display = unname(compact_covariate_labels[.data$covariate_family]),
    covariate_display = factor(.data$covariate_display, levels = unname(covariate_labels[covariate_order])),
    compact_covariate_display = factor(.data$compact_covariate_display, levels = unname(compact_covariate_labels[covariate_order])),
    module_display = factor(.data$compact_label, levels = rev(class_for_plot$compact_label[match(row_order, class_for_plot$module_id)])),
    compact_module_display = factor(.data$compact_plot_label, levels = rev(class_for_plot$compact_plot_label[match(row_order, class_for_plot$module_id)])),
    highlight_4D4D4D = .data$module_id == "WGCNA_#4D4D4D",
    claim_gate_eligible = .data$module_claim_gate_eligible,
    cell_label = paste0(round(.data$retained_abs_effect_pct), dplyr::if_else(!dplyr::coalesce(.data$direction_preserved, TRUE), "\u00d7", ""))
  ) |>
  dplyr::filter(is.finite(.data$retained_abs_effect_pct)) |>
  dplyr::select(
    "module_id", "module_display", "compact_module_display", "endpoint_label", "biological_label", "compact_label", "compact_plot_label", "best_contrast", "covariate_family", "covariate_display", "compact_covariate_display",
    "effect_before", "effect_after", "retained_abs_effect_pct", "retained_abs_effect_pct_display_capped",
    "direction_preserved", "cell_label", "claim_gate_eligible", "min_p_before", "max_percent_attenuation", "highlight_4D4D4D"
  ) |>
  dplyr::arrange(factor(.data$module_id, levels = row_order), factor(.data$covariate_family, levels = covariate_order))

dumbbell_source <- effects |>
  dplyr::filter(.data$module_id == "WGCNA_#4D4D4D", .data$endpoint_type == "module_eigengene") |>
  dplyr::mutate(
    effect_before = suppressWarnings(as.numeric(.data$effect_before)),
    effect_after = suppressWarnings(as.numeric(.data$effect_after)),
    direction_preserved = as.logical(.data$direction_preserved),
    covariate_display = unname(covariate_labels[.data$covariate_family]),
    covariate_display = dplyr::coalesce(.data$covariate_display, .data$covariate_family),
    contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "SUS - RES")),
    row_label = paste(as.character(.data$contrast), .data$covariate_display, sep = " | ")
  ) |>
  dplyr::filter(is.finite(.data$effect_before), is.finite(.data$effect_after), !is.na(.data$contrast)) |>
  dplyr::arrange(.data$contrast, factor(.data$covariate_family, levels = covariate_order)) |>
  dplyr::mutate(row_label = factor(.data$row_label, levels = rev(unique(.data$row_label)))) |>
  dplyr::select(
    "module_id", "contrast", "covariate_family", "covariate_display", "row_label",
    "effect_before", "effect_after", "direction_preserved", "p_before", "p_after",
    "claim_gate_eligible", "independence_classification", "downgrade_reason"
  )

dir_create(fig_dir)
dir_create(source_dir)
dir_create(log_dir)
readr::write_csv(diagnostic_source, source_outputs["diagnostic"], na = "")
readr::write_csv(heatmap_source, source_outputs["heatmap"], na = "")
if (nrow(dumbbell_source)) readr::write_csv(dumbbell_source, source_outputs["dumbbell"], na = "")

save_svg_pdf <- function(plot, stem, width_mm, height_mm) {
  svg_path <- file.path(fig_dir, paste0(stem, ".svg"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(svg_path, plot = plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = svg_device, limitsize = FALSE)
  ggplot2::ggsave(pdf_path, plot = plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = pdf_device, limitsize = FALSE)
  c(svg_path, pdf_path)
}

poster_teal <- "#176B57"
poster_mint <- "#A8DCCF"
poster_grey <- "#5F6568"
light_grey <- "#E8ECEC"
dark_grey <- "#2D3436"
poster_red <- "#C94A5A"
poster_blue <- "#6D8FB3"

plot_base_size <- if (isTRUE(COMPACT)) 6.8 else 7.5
diagnostic_dims <- if (isTRUE(COMPACT)) c(45, 50) else c(120, 92)
heatmap_dims <- if (isTRUE(COMPACT)) c(105, 80) else c(128, 104)
dumbbell_dims <- if (isTRUE(COMPACT)) c(105, 65) else c(118, 92)

base_theme <- theme_nature_base(base_size = plot_base_size, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", colour = dark_grey, size = 8.5),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(hjust = 0, size = 5.5, colour = poster_grey, lineheight = 0.95),
    axis.text = ggplot2::element_text(colour = dark_grey, size = 6),
    axis.title = ggplot2::element_text(colour = dark_grey, size = 7),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(2, 3, 2, 2, "pt")
  )

diagnostic_label_source <- diagnostic_source |>
  dplyr::filter(isTRUE(LABEL_ALL) | .data$module_id %in% names(diagnostic_key_labels)) |>
  dplyr::mutate(
    diagnostic_label = dplyr::coalesce(unname(diagnostic_key_labels[.data$module_id]), .data$compact_plot_label),
    label_colour = dplyr::if_else(.data$module_id == "WGCNA_#4D4D4D", poster_teal, poster_grey),
    point_fill = dplyr::case_when(
      .data$module_id == "WGCNA_#4D4D4D" ~ poster_teal,
      .data$module_id == "WGCNA_#3F007D" ~ "#7EA99D",
      .data$module_id == "WGCNA_#7F2704" ~ poster_grey,
      TRUE ~ poster_grey
    )
  )

diagnostic_plot <- ggplot2::ggplot(diagnostic_source, ggplot2::aes(x = .data$neg_log10_min_p_before, y = .data$diagnostic_stability_score)) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.18, colour = light_grey) +
  ggplot2::geom_point(shape = 21, size = 1.0, stroke = 0.12, colour = "white", fill = light_grey) +
  ggplot2::geom_point(
    data = diagnostic_label_source,
    ggplot2::aes(fill = .data$point_fill),
    shape = 21,
    size = 2.3,
    stroke = 0.25,
    colour = "white"
  ) +
    ggplot2::scale_fill_identity() +
  ggplot2::labs(
    title = "Module-level neuropil adjustment",
    x = "-log10(best nominal p before adjustment)",
    y = "Stability score",
    caption = "Matched-neuropil diagnostics from AnimalID + Region adjusted models;\nall modules are diagnostic unless claim_gate_eligible = TRUE."
  ) +
  base_theme

if (requireNamespace("ggrepel", quietly = TRUE)) {
  diagnostic_plot <- diagnostic_plot +
    ggrepel::geom_text_repel(
      data = diagnostic_label_source,
      ggplot2::aes(label = .data$diagnostic_label),
      size = 2.0,
      colour = dark_grey,
      min.segment.length = 0,
      segment.size = 0.18,
      segment.colour = poster_grey,
      box.padding = 0.16,
      point.padding = 0.08,
      max.overlaps = Inf,
      seed = 9
    )
} else {
  diagnostic_plot <- diagnostic_plot +
    ggplot2::geom_text(
      data = diagnostic_label_source,
      ggplot2::aes(label = .data$diagnostic_label),
      size = 2.0,
      colour = dark_grey,
      hjust = -0.05,
      vjust = -0.25,
      check_overlap = TRUE
    )
}
diagnostic_files <- save_svg_pdf(diagnostic_plot, "microglia_neuropil_diagnostic_screen", diagnostic_dims[[1]], diagnostic_dims[[2]])

heatmap_visual_covariates <- if (isTRUE(INCLUDE_EXPLORATORY)) {
  c(
    "global_neuropil_score",
    "mitochondrial_neuropil_score",
    "ECM_adhesion_neuropil_score",
    "RNA_translation_neuropil_score",
    "exploratory_best_spearman"
  )
} else {
  c(
    "global_neuropil_score",
    "mitochondrial_neuropil_score",
    "ECM_adhesion_neuropil_score",
    "RNA_translation_neuropil_score"
  )
}
heatmap_plot_source <- heatmap_source |>
  dplyr::filter(.data$covariate_family %in% heatmap_visual_covariates) |>
  dplyr::mutate(
    compact_covariate_display = factor(.data$compact_covariate_display, levels = unname(compact_covariate_labels[heatmap_visual_covariates])),
    compact_module_display = factor(.data$compact_module_display, levels = levels(heatmap_source$compact_module_display))
  )

outline_dat <- heatmap_plot_source |>
  dplyr::filter(.data$highlight_4D4D4D) |>
  dplyr::distinct(.data$compact_module_display)

heatmap_plot <- ggplot2::ggplot(heatmap_plot_source, ggplot2::aes(x = .data$compact_covariate_display, y = .data$compact_module_display)) +
  ggplot2::geom_tile(ggplot2::aes(fill = .data$retained_abs_effect_pct_display_capped), colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = .data$cell_label), size = 2.0, colour = dark_grey) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7F8F8", light_grey, poster_mint, poster_teal),
    values = scales::rescale(c(0, 75, 150, 250)),
    limits = c(0, 250),
    oob = scales::squish,
    name = "Retained (%)"
  ) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(30, "mm"), barwidth = grid::unit(3, "mm"))) +
  ggplot2::coord_fixed(ratio = 0.78, clip = "off") +
  ggplot2::labs(
    title = "Retained effect after adjustment",
    x = "Covariate",
    y = NULL,
    caption = "Values are retained absolute effect (%); \u00d7 marks direction flip."
  ) +
  theme_nature_heatmap(base_size = plot_base_size, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", colour = dark_grey, size = 8.5),
    plot.caption = ggplot2::element_text(hjust = 0, size = 5.5, colour = poster_grey, lineheight = 0.95),
    axis.title.x = ggplot2::element_text(colour = dark_grey, size = 7),
    axis.text = ggplot2::element_text(colour = dark_grey, size = 6),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
    legend.title = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 5.5),
    legend.key.height = grid::unit(8, "mm"),
    legend.box.spacing = grid::unit(0.5, "mm"),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = ggplot2::margin(2, 3, 2, 2, "pt")
  )

if (nrow(outline_dat)) {
  x_min <- 0.5
  x_max <- length(levels(heatmap_plot_source$compact_covariate_display)) + 0.5
  outline_rect <- outline_dat |>
    dplyr::mutate(y_num = as.numeric(.data$compact_module_display), xmin = x_min, xmax = x_max, ymin = .data$y_num - 0.5, ymax = .data$y_num + 0.5)
  heatmap_plot <- heatmap_plot +
    ggplot2::geom_rect(
      data = outline_rect,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = poster_teal,
      linewidth = 0.35
    )
}
heatmap_files <- save_svg_pdf(heatmap_plot, "microglia_neuropil_retained_effect_heatmap", heatmap_dims[[1]], heatmap_dims[[2]])

dumbbell_files <- character()
if (nrow(dumbbell_source)) {
  dumbbell_visual_covariates <- if (isTRUE(INCLUDE_EXPLORATORY)) {
    c("global_neuropil_score", "mitochondrial_neuropil_score", "ECM_adhesion_neuropil_score", "RNA_translation_neuropil_score", "exploratory_best_spearman")
  } else {
    c("global_neuropil_score", "mitochondrial_neuropil_score", "ECM_adhesion_neuropil_score", "RNA_translation_neuropil_score")
  }
  dumbbell_plot_source <- dumbbell_source |>
    dplyr::filter(.data$covariate_family %in% dumbbell_visual_covariates) |>
    dplyr::mutate(
      covariate_display = factor(.data$covariate_display, levels = rev(unname(covariate_labels[dumbbell_visual_covariates]))),
      contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "SUS - RES"))
    )
  dumbbell_long <- dumbbell_source |>
    dplyr::filter(.data$covariate_family %in% dumbbell_visual_covariates) |>
    tidyr::pivot_longer(
      cols = c("effect_before", "effect_after"),
      names_to = "timepoint",
      values_to = "effect"
    ) |>
    dplyr::mutate(
      timepoint = factor(.data$timepoint, levels = c("effect_before", "effect_after"), labels = c("Before", "After")),
      covariate_display = factor(.data$covariate_display, levels = levels(dumbbell_plot_source$covariate_display)),
      contrast = factor(.data$contrast, levels = levels(dumbbell_plot_source$contrast))
    )
  direction_all_preserved <- all(dumbbell_source$direction_preserved, na.rm = TRUE)
  dumbbell_caption <- if (isTRUE(direction_all_preserved)) {
    "Direction is preserved across adjustments; effects remain nominal/non-claim-gated."
  } else {
    "Direction is not preserved for all adjustments; interpret as diagnostic."
  }
  dumbbell_plot <- ggplot2::ggplot(dumbbell_plot_source, ggplot2::aes(y = .data$covariate_display)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.18, colour = light_grey) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$effect_before, xend = .data$effect_after, yend = .data$covariate_display), colour = "#C9CECE", linewidth = 0.28) +
    ggplot2::geom_point(data = dumbbell_long, ggplot2::aes(x = .data$effect, colour = .data$timepoint), size = 1.9) +
    ggplot2::facet_wrap(~ contrast, ncol = 3) +
    ggplot2::scale_colour_manual(values = c(Before = poster_blue, After = poster_red), name = NULL) +
    ggplot2::labs(
      title = "ECM/adhesion module sensitivity",
      x = "Contrast estimate",
      y = NULL,
      caption = dumbbell_caption
    ) +
    base_theme +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "center",
      legend.key.size = grid::unit(2.4, "mm"),
      legend.spacing.x = grid::unit(1.2, "mm"),
      strip.text = ggplot2::element_text(size = 6, face = "bold", colour = dark_grey),
      panel.spacing.x = grid::unit(1.4, "mm"),
      axis.text.y = ggplot2::element_text(size = 5.8),
      plot.margin = ggplot2::margin(2, 3, 2, 2, "pt")
    )
  dumbbell_files <- save_svg_pdf(dumbbell_plot, "WGCNA_4D4D4D_neuropil_adjustment_dumbbell", dumbbell_dims[[1]], dumbbell_dims[[2]])
}

write_run_manifest(
  manifest_path,
  inputs = inputs,
  outputs = list(
    figures = c(diagnostic_files, heatmap_files, dumbbell_files),
    source_data = c(source_outputs[c("diagnostic", "heatmap")], if (nrow(dumbbell_source)) source_outputs["dumbbell"] else character())
  ),
  parameters = list(
    dataset = dataset,
    source_tables = "Existing outputs from 06_modules_WGCNA/09_microglia_neuropil_independence.R",
    model_refit = FALSE,
    wgcna_recomputed = FALSE,
    neuropil_matching = "neuron_neuropil covariates aggregated by AnimalID + Region",
    claim_wording_guard = "Exploratory_best_spearman rows and non-claim-gated modules are diagnostic only; supermodules are not directly tested here."
  ),
  notes = "Publication-style rendering of existing microglia matched-neuropil adjustment diagnostics; no model refitting and no WGCNA recomputation."
)

message("Microglia matched-neuropil diagnostic figures complete.")
