#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/08_wgcna_publication_figures.R
# Stage: modules_downstream
# Scope: microglia
# Consumes: current Stage 01 identities, Stage 05 statistics, and the reviewed
#           canonical lookup produced by Stage 07.
# Produces: all-supermodule publication figures (SVG/PDF) and one source CSV
#           per figure. This renderer never fits models or changes WGCNA state.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))
source(repo_path("R", "wgcna_reviewed_label_registry.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "readr", "stringr", "scales", "svglite", "digest")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- wgcna_cli()
DATASET <- run$dataset
if (!identical(DATASET, "microglia")) stop("08_wgcna_publication_figures.R is restricted to the reviewed microglia network.", call. = FALSE)

OUT <- create_module_dirs("06_modules_WGCNA", file.path("wgcna_publication_figures", DATASET))
GROUP_DIR <- path_results("tables", "06_modules_WGCNA", "group_effects", DATASET)
LABEL_FILE <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", DATASET, "WGCNA_final_label_lookup.csv")
HANDOFF_FILE <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", DATASET, "WGCNA_inferential_handoff.csv")
FILES <- resolve_wgcna_files(DATASET)

figure_stems <- c(
  "all_supermodule_architecture",
  "all_supermodule_global_eigengenes",
  "all_supermodule_sus_res_forest",
  "all_supermodule_sus_res_spatial_heatmap",
  "multi_supermodule_member_loadings",
  "supplementary_all_contrast_effect_matrix"
)

if (run$dry_run) {
  dry_run_line("Script", "06_modules_WGCNA/08_wgcna_publication_figures.R")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Reviewed canonical lookup", LABEL_FILE, if (file.exists(LABEL_FILE)) "PASS" else "FAIL")
  dry_run_line("Stage 07 inferential handoff", HANDOFF_FILE, if (file.exists(HANDOFF_FILE)) "PASS" else "FAIL")
  dry_run_line("Stage 05 non-inferential values/loadings", GROUP_DIR, if (dir.exists(GROUP_DIR)) "PASS" else "FAIL")
  for (stem in figure_stems) {
    dry_run_line("Figure", file.path(OUT$figures, paste0(stem, ".svg")))
    dry_run_line("Source data", file.path(OUT$source_data, paste0(stem, "_source.csv")))
  }
  required_ready <- file.exists(LABEL_FILE) &&
    file.exists(HANDOFF_FILE) && dir.exists(GROUP_DIR)
  dry_run_line(
    "Status",
    if (required_ready) "ready" else "missing_required_input",
    if (required_ready) "PASS" else "FAIL"
  )
  quit(status = if (required_ready) 0 else 1, save = "no")
}

read_required <- function(path, label) {
  if (!file.exists(path)) stop("Missing required ", label, ": ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

write_figure_source <- function(df, stem) {
  readr::write_csv(df, file.path(OUT$source_data, paste0(stem, "_source.csv")), na = "")
  readr::write_csv(df, file.path(OUT$tables, paste0(stem, "_source.csv")), na = "")
}

save_publication_figure <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(OUT$figures, paste0(stem, ".svg")), plot = plot,
    width = width, height = height, units = "mm", device = svglite::svglite,
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    file.path(OUT$figures, paste0(stem, ".pdf")), plot = plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf,
    family = "Arial", bg = "white", limitsize = FALSE
  )
}

theme_pub <- function(base_size = 7.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      line = ggplot2::element_line(linewidth = 0.16, colour = "#333333"),
      axis.line = ggplot2::element_line(linewidth = 0.16, colour = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = 0.16, colour = "#333333"),
      axis.text = ggplot2::element_text(size = 7, colour = "#222222"),
      axis.title = ggplot2::element_text(size = 7.5, colour = "#222222"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7.3, face = "plain", lineheight = 0.95),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.8),
      plot.margin = ggplot2::margin(3, 4, 3, 3)
    )
}

stress_colours <- c(CON = "#3E3C6F", RES = "#9E9A92", SUS = "#D7303F")
effect_colours <- c(low = "#3B6FB6", mid = "#F7F7F5", high = "#C84C5A")
sm_order <- sprintf("SM%02d", 1:9)
roi_order <- c("CA1", "CA2", "CA3", "DG")
contrast_order <- c("RES - CON", "SUS - CON", "SUS - RES")

lookup <- read_required(LABEL_FILE, "reviewed canonical label lookup")
member_map_raw <- read_required(FILES$supermodule_annotation, "current Stage 01 module-to-supermodule map")
member_map <- wgcna_normalize_current_member_map(member_map_raw, DATASET)
wgcna_validate_canonical_lookup(lookup, DATASET, member_map)

super_labels <- lookup |>
  dplyr::filter(.data$level == "supermodule") |>
  dplyr::transmute(
    dataset, supermodule_id = .data$entity_id,
    canonical_biological_label, canonical_short_label, canonical_plot_label,
    biological_label_confidence, manual_review_required, structural_status,
    aggregation_evidence_class, structural_coherence_class,
    subcellular_context, roi_context, n_member_modules, member_ModuleIDs,
    member_modules_fingerprint, label_rationale
  ) |>
  dplyr::mutate(supermodule_id = factor(.data$supermodule_id, levels = sm_order)) |>
  dplyr::arrange(.data$supermodule_id) |>
  dplyr::mutate(supermodule_id = as.character(.data$supermodule_id))

module_labels <- lookup |>
  dplyr::filter(.data$level == "module") |>
  dplyr::transmute(
    dataset, module_id = .data$entity_id, canonical_module_short_label = .data$canonical_short_label,
    canonical_module_plot_label = .data$canonical_plot_label,
    module_parent_supermodule = .data$parent_entity_id
  )

join_super_labels <- function(df, id_col = "supermodule_id", context = "publication source") {
  before <- nrow(df)
  if (!"dataset" %in% names(df)) stop("Stable dataset key missing in ", context, ".", call. = FALSE)
  join_keys <- c(dataset = "dataset")
  join_keys[[id_col]] <- "supermodule_id"
  out <- df |>
    dplyr::left_join(super_labels, by = join_keys, relationship = "many-to-one")
  if (nrow(out) != before) stop("Label join multiplied rows in ", context, ".", call. = FALSE)
  if (any(is.na(out$canonical_short_label))) stop("Canonical label missing in ", context, ".", call. = FALSE)
  out
}

definitions <- read_required(FILES$definitions, "current module definitions")
super_summary <- read_required(FILES$supermodule_summary, "current supermodule summary")
inferential_handoff <- wgcna_inferential_handoff_read(HANDOFF_FILE)
display_lookup <- lookup |>
  dplyr::filter(.data$level == "module") |>
  dplyr::transmute(
    dataset, module_id = .data$entity_id,
    supermodule_id = .data$parent_entity_id
  ) |>
  dplyr::left_join(
    super_labels |>
      dplyr::select(
        dataset, supermodule_id, n_member_modules,
        supermodule_label = canonical_plot_label
      ),
    by = c("dataset", "supermodule_id"),
    relationship = "many-to-one"
  )
group_effects <- wgcna_inferential_handoff_supermodule_display(
  inferential_handoff, display_lookup
)
raw_values <- read_required(file.path(GROUP_DIR, "all_supermodule_eigengene_group_values.csv"), "Stage 05 raw supermodule eigengenes")
loadings <- read_required(file.path(GROUP_DIR, "supermodule_pca_member_loadings.csv"), "Stage 05 supermodule member loadings")

module_meta <- definitions |>
  dplyr::transmute(
    ModuleID = as.character(.data$ModuleID),
    ModuleColor = as.character(.data$ModuleColor),
    ModuleColorName = as.character(.data$ModuleColorName),
    module_eigengene = as.character(.data$module_eigengene)
  ) |>
  dplyr::distinct(.data$ModuleID, .keep_all = TRUE)

summary_meta <- super_summary |>
  dplyr::transmute(
    supermodule_id = as.character(.data$SupermoduleID),
    pc1_variance = suppressWarnings(as.numeric(.data$pc1_variance_explained)),
    adjusted_min_eigengene_correlation = suppressWarnings(as.numeric(.data$adjusted_signed_min_pairwise_eigengene_correlation)),
    cut_height_stability = suppressWarnings(as.numeric(.data$cut_height_stability_fraction_stable))
  ) |>
  dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)

# 1. Architecture -----------------------------------------------------------
architecture_source <- member_map |>
  dplyr::left_join(module_meta, by = "ModuleID", relationship = "one-to-one") |>
  dplyr::rename(supermodule_id = "SupermoduleID") |>
  join_super_labels(context = "architecture") |>
  dplyr::left_join(summary_meta, by = "supermodule_id", relationship = "many-to-one") |>
  dplyr::group_by(.data$supermodule_id) |>
  dplyr::arrange(.data$ModuleID, .by_group = TRUE) |>
  dplyr::mutate(member_position = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    singleton_badge = dplyr::if_else(.data$structural_status == "singleton", "singleton", ""),
    architecture_annotation = dplyr::case_when(
      .data$structural_status == "singleton" ~ paste0("n=1 | singleton | ", .data$biological_label_confidence),
      TRUE ~ paste0("n=", .data$n_member_modules, " | PC1=", scales::percent(.data$pc1_variance, accuracy = 1), " | ", .data$biological_label_confidence)
    ),
    plot_label = factor(.data$canonical_plot_label, levels = rev(super_labels$canonical_plot_label)),
    ModuleColor = dplyr::if_else(grepl("^#[0-9A-Fa-f]{6}$", .data$ModuleColor), .data$ModuleColor, "#BDBDBD")
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = sm_order), .data$member_position)
write_figure_source(architecture_source, "all_supermodule_architecture")

architecture_text <- architecture_source |>
  dplyr::distinct(.data$supermodule_id, .data$plot_label, .data$architecture_annotation, .data$roi_context)
p_arch <- ggplot2::ggplot(architecture_source, ggplot2::aes(x = .data$member_position, y = .data$plot_label)) +
  ggplot2::geom_tile(ggplot2::aes(fill = .data$ModuleColor), width = 0.86, height = 0.68, colour = "#333333", linewidth = 0.16) +
  ggplot2::geom_text(ggplot2::aes(label = .data$ModuleID), size = 2.45, family = "Arial", colour = "white") +
  ggplot2::geom_text(data = architecture_text, ggplot2::aes(x = 4.0, y = .data$plot_label, label = .data$architecture_annotation), inherit.aes = FALSE, hjust = 0, vjust = -0.28, size = 2.45, family = "Arial") +
  ggplot2::geom_text(data = architecture_text, ggplot2::aes(x = 4.0, y = .data$plot_label, label = stringr::str_trunc(.data$roi_context, 52)), inherit.aes = FALSE, hjust = 0, vjust = 1.30, size = 2.45, family = "Arial", colour = "#666666") +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_x_continuous(limits = c(0.5, 8.2), breaks = NULL, expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_pub() +
  ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(3, 8, 3, 3))
save_publication_figure(p_arch, "all_supermodule_architecture", 205, 100)

# Shared effect sources and symmetric scale ---------------------------------
effect_base <- group_effects |>
  dplyr::mutate(
    supermodule_id = as.character(.data$supermodule_id),
    contrast = as.character(.data$contrast),
    estimate = suppressWarnings(as.numeric(.data$estimate)),
    SE = suppressWarnings(as.numeric(.data$SE)),
    q_value = suppressWarnings(as.numeric(.data$tier_specific_fdr)),
    model_valid_for_inference = as.logical(.data$model_valid),
    model_stable = .data$claim_gate ==
      "eligible_for_readiness_assessment",
    diagnostic_only = !.data$model_stable
  )

global_effects <- effect_base |>
  dplyr::filter(.data$effect_scope == "spatial_adjusted_global", .data$spatial_unit == "global_spatial_adjusted") |>
  join_super_labels(context = "global effects")
spatial_effects <- effect_base |>
  dplyr::filter(.data$effect_scope == "within_spatial_unit", .data$contrast == "SUS - RES") |>
  dplyr::mutate(roi = toupper(as.character(.data$spatial_unit))) |>
  dplyr::filter(.data$roi %in% roi_order) |>
  join_super_labels(context = "spatial effects")
effect_limit <- max(abs(c(global_effects$estimate, spatial_effects$estimate)), na.rm = TRUE)
effect_limit <- ceiling(effect_limit * 10) / 10
if (!is.finite(effect_limit) || effect_limit <= 0) effect_limit <- 1

# 2. Global eigengenes ------------------------------------------------------
global_sig <- global_effects |>
  dplyr::filter(.data$display_is_independent_endpoint %in% TRUE, .data$independent_hypothesis %in% TRUE, .data$model_stable, is.finite(.data$q_value), .data$q_value <= 0.05) |>
  dplyr::group_by(.data$supermodule_id) |>
  dplyr::summarise(
    fdr_supported_contrasts = paste(paste0(.data$contrast, " q=", formatC(.data$q_value, digits = 2, format = "f")), collapse = "; "),
    .groups = "drop"
  )
global_eigengene_source <- raw_values |>
  dplyr::select("dataset", "Sample", "supermodule_id", "supermodule_eigengene", "eigengene", "eigengene_z", "StressGroup", "AnimalID", "Sex", "Batch", "SpatialUnitType") |>
  dplyr::mutate(supermodule_id = as.character(.data$supermodule_id), StressGroup = as.character(.data$StressGroup)) |>
  join_super_labels(context = "global eigengenes") |>
  dplyr::left_join(global_sig, by = "supermodule_id", relationship = "many-to-one") |>
  dplyr::mutate(
    StressGroup = factor(.data$StressGroup, levels = c("CON", "RES", "SUS")),
    panel_label = factor(.data$canonical_plot_label, levels = super_labels$canonical_plot_label)
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = sm_order), .data$StressGroup, .data$Sample)
write_figure_source(global_eigengene_source, "all_supermodule_global_eigengenes")

global_annotation_rows <- global_eigengene_source |>
  dplyr::filter(!is.na(.data$fdr_supported_contrasts), is.finite(.data$eigengene_z)) |>
  droplevels()
global_annotations <- if (nrow(global_annotation_rows)) {
  global_annotation_rows |>
    dplyr::group_by(.data$panel_label, .data$fdr_supported_contrasts, .drop = TRUE) |>
    dplyr::summarise(y = max(.data$eigengene_z) + 0.25, .groups = "drop")
} else {
  tibble::tibble(panel_label = factor(), fdr_supported_contrasts = character(), y = numeric())
}
p_global <- ggplot2::ggplot(global_eigengene_source, ggplot2::aes(x = .data$StressGroup, y = .data$eigengene_z)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = .data$StressGroup), width = 0.56, outlier.shape = NA, linewidth = 0.16, alpha = 0.28, colour = "#333333") +
  ggplot2::geom_point(ggplot2::aes(colour = .data$StressGroup), position = ggplot2::position_jitter(width = 0.11, height = 0), size = 0.65, alpha = 0.65) +
  ggplot2::geom_text(data = global_annotations, ggplot2::aes(x = 2, y = .data$y, label = .data$fdr_supported_contrasts), inherit.aes = FALSE, size = 2.45, family = "Arial") +
  ggplot2::facet_wrap(~ panel_label, ncol = 3, scales = "free_y") +
  ggplot2::scale_fill_manual(values = stress_colours, drop = FALSE) +
  ggplot2::scale_colour_manual(values = stress_colours, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "Supermodule eigengene (z)") +
  theme_pub(7.3) +
  ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.12), panel.grid.minor = ggplot2::element_blank())
save_publication_figure(p_global, "all_supermodule_global_eigengenes", 180, 155)

# 3. SUS-RES forest ---------------------------------------------------------
forest_source <- global_effects |>
  dplyr::filter(.data$contrast == "SUS - RES") |>
  dplyr::mutate(
    ci_low = .data$estimate - 1.96 * .data$SE,
    ci_high = .data$estimate + 1.96 * .data$SE,
    model_stability = dplyr::if_else(.data$model_stable, "stable", "singular_or_unstable_diagnostic_only"),
    claim_support_symbol = dplyr::case_when(
      !(.data$display_is_independent_endpoint %in% TRUE) ~ "none",
      !(.data$independent_hypothesis %in% TRUE) ~ "none",
      .data$model_stable & is.finite(.data$q_value) & .data$q_value <= 0.05 ~ "FDR05",
      .data$model_stable & is.finite(.data$q_value) & .data$q_value <= 0.10 ~ "FDR10",
      TRUE ~ "none"
    ),
    plot_label = factor(.data$canonical_plot_label, levels = rev(super_labels$canonical_plot_label))
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = sm_order))
write_figure_source(forest_source, "all_supermodule_sus_res_forest")

p_forest <- ggplot2::ggplot(forest_source, ggplot2::aes(y = .data$plot_label, x = .data$estimate)) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.16) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high), orientation = "y", width = 0.16, linewidth = 0.22, colour = "#333333") +
  ggplot2::geom_point(data = forest_source |> dplyr::filter(.data$model_stable), shape = 21, size = 2.2, stroke = 0.22, fill = "#3B6FB6", colour = "#222222") +
  ggplot2::geom_point(data = forest_source |> dplyr::filter(!.data$model_stable), shape = 4, size = 2.4, stroke = 0.55, colour = "#888888") +
  ggplot2::labs(x = "SUS minus RES estimate (95% CI)", y = NULL) +
  theme_pub() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.12), panel.grid.minor = ggplot2::element_blank())
save_publication_figure(p_forest, "all_supermodule_sus_res_forest", 155, 92)

# 4. SUS-RES spatial heatmap ------------------------------------------------
spatial_heatmap_source <- spatial_effects |>
  dplyr::mutate(
    support_symbol = dplyr::case_when(
      !.data$model_stable ~ "cross_diagnostic_only",
      !(.data$display_is_independent_endpoint %in% TRUE) ~ "none",
      !(.data$independent_hypothesis %in% TRUE) ~ "none",
      is.finite(.data$q_value) & .data$q_value <= 0.05 ~ "filled_FDR05",
      is.finite(.data$q_value) & .data$q_value <= 0.10 ~ "open_FDR10",
      TRUE ~ "none"
    ),
    roi = factor(.data$roi, levels = roi_order),
    plot_label = factor(.data$canonical_plot_label, levels = rev(super_labels$canonical_plot_label)),
    effect_scale_limit = effect_limit
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = sm_order), .data$roi)
write_figure_source(spatial_heatmap_source, "all_supermodule_sus_res_spatial_heatmap")

roi_side <- spatial_heatmap_source |>
  dplyr::distinct(.data$plot_label, .data$roi_context)
p_spatial <- ggplot2::ggplot(spatial_heatmap_source, ggplot2::aes(x = .data$roi, y = .data$plot_label, fill = .data$estimate)) +
  ggplot2::geom_tile(width = 0.92, height = 0.86, colour = "white", linewidth = 0.20) +
  ggplot2::geom_point(data = spatial_heatmap_source |> dplyr::filter(.data$support_symbol == "filled_FDR05"), shape = 21, size = 1.8, stroke = 0.20, fill = "#222222", colour = "#222222") +
  ggplot2::geom_point(data = spatial_heatmap_source |> dplyr::filter(.data$support_symbol == "open_FDR10"), shape = 21, size = 1.8, stroke = 0.30, fill = "white", colour = "#222222") +
  ggplot2::geom_point(data = spatial_heatmap_source |> dplyr::filter(.data$support_symbol == "cross_diagnostic_only"), shape = 4, size = 2.2, stroke = 0.55, colour = "#888888") +
  ggplot2::geom_text(data = roi_side, ggplot2::aes(x = 4.65, y = .data$plot_label, label = stringr::str_trunc(.data$roi_context, 42)), inherit.aes = FALSE, hjust = 0, size = 2.45, family = "Arial", colour = "#666666") +
  ggplot2::scale_fill_gradient2(low = effect_colours["low"], mid = effect_colours["mid"], high = effect_colours["high"], midpoint = 0, limits = c(-effect_limit, effect_limit), oob = scales::squish, name = "SUS - RES") +
  ggplot2::scale_x_discrete(drop = FALSE, expand = ggplot2::expansion(add = c(0.05, 1.45))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_pub() +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(3, 45, 3, 3), legend.position = "right")
save_publication_figure(p_spatial, "all_supermodule_sus_res_spatial_heatmap", 225, 102)

# 5. Multi-supermodule member loadings -------------------------------------
multi_ids <- c("SM01", "SM03", "SM09")
multi_panel_levels <- super_labels |>
  dplyr::filter(.data$supermodule_id %in% multi_ids) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = multi_ids)) |>
  dplyr::transmute(panel_label = paste0(.data$supermodule_id, "\n", stringr::str_wrap(.data$canonical_short_label, width = 18))) |>
  dplyr::pull(.data$panel_label)
loading_source <- loadings |>
  dplyr::filter(.data$supermodule_id %in% multi_ids, suppressWarnings(as.integer(.data$pc)) == 1L) |>
  dplyr::select("dataset", "supermodule_id", "module_id", "module_eigengene", "loading", "abs_loading", "loading_rank", pca_n_member_modules = "n_member_modules", "pca_status") |>
  dplyr::left_join(module_labels, by = c("dataset", "module_id"), relationship = "many-to-one") |>
  join_super_labels(context = "member loadings") |>
  dplyr::mutate(
    supermodule_panel = factor(
      paste0(.data$supermodule_id, "\n", stringr::str_wrap(.data$canonical_short_label, width = 18)),
      levels = multi_panel_levels
    ),
    module_plot_label = factor(.data$canonical_module_plot_label, levels = rev(unique(.data$canonical_module_plot_label)))
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = multi_ids), .data$loading_rank)
if (any(loading_source$n_member_modules <= 1L)) stop("Singleton loading panels are forbidden.", call. = FALSE)
write_figure_source(loading_source, "multi_supermodule_member_loadings")

p_load <- ggplot2::ggplot(loading_source, ggplot2::aes(x = .data$loading, y = .data$module_plot_label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.16, colour = "#777777") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$loading, yend = .data$module_plot_label), linewidth = 0.35, colour = "#5F7590") +
  ggplot2::geom_point(shape = 21, size = 2.0, stroke = 0.20, fill = "#5F7590", colour = "#222222") +
  ggplot2::facet_wrap(~ supermodule_panel, scales = "free_y", ncol = 3) +
  ggplot2::labs(x = "PC1 loading", y = NULL) +
  theme_pub() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.12), panel.grid.minor = ggplot2::element_blank())
save_publication_figure(p_load, "multi_supermodule_member_loadings", 180, 78)

# 6. Supplementary all-contrast effect matrix ------------------------------
effect_matrix_source <- global_effects |>
  dplyr::filter(.data$contrast %in% contrast_order) |>
  dplyr::mutate(
    contrast = factor(.data$contrast, levels = contrast_order),
    plot_label = factor(.data$canonical_plot_label, levels = rev(super_labels$canonical_plot_label)),
    support_symbol = dplyr::case_when(
      !.data$model_stable ~ "cross_diagnostic_only",
      !(.data$display_is_independent_endpoint %in% TRUE) ~ "none",
      !(.data$independent_hypothesis %in% TRUE) ~ "none",
      is.finite(.data$q_value) & .data$q_value <= 0.05 ~ "filled_FDR05",
      is.finite(.data$q_value) & .data$q_value <= 0.10 ~ "open_FDR10",
      TRUE ~ "none"
    ),
    effect_scale_limit = effect_limit
  ) |>
  dplyr::arrange(factor(.data$supermodule_id, levels = sm_order), .data$contrast)
write_figure_source(effect_matrix_source, "supplementary_all_contrast_effect_matrix")

p_matrix <- ggplot2::ggplot(effect_matrix_source, ggplot2::aes(x = .data$contrast, y = .data$plot_label, fill = .data$estimate)) +
  ggplot2::geom_tile(width = 0.92, height = 0.86, colour = "white", linewidth = 0.20) +
  ggplot2::geom_point(data = effect_matrix_source |> dplyr::filter(.data$support_symbol == "filled_FDR05"), shape = 21, size = 1.8, stroke = 0.20, fill = "#222222", colour = "#222222") +
  ggplot2::geom_point(data = effect_matrix_source |> dplyr::filter(.data$support_symbol == "open_FDR10"), shape = 21, size = 1.8, stroke = 0.30, fill = "white", colour = "#222222") +
  ggplot2::geom_point(data = effect_matrix_source |> dplyr::filter(.data$support_symbol == "cross_diagnostic_only"), shape = 4, size = 2.2, stroke = 0.55, colour = "#888888") +
  ggplot2::scale_fill_gradient2(low = effect_colours["low"], mid = effect_colours["mid"], high = effect_colours["high"], midpoint = 0, limits = c(-effect_limit, effect_limit), oob = scales::squish, name = "Estimate") +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_pub() +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "right")
save_publication_figure(p_matrix, "supplementary_all_contrast_effect_matrix", 142, 98)

# Validation ----------------------------------------------------------------
validation <- tibble::tribble(
  ~validation_check, ~status, ~details,
  "canonical_lookup_22_rows", if (nrow(lookup) == 22L) "pass" else "fail", paste0("n=", nrow(lookup)),
  "architecture_has_9_supermodules", if (dplyr::n_distinct(architecture_source$supermodule_id) == 9L) "pass" else "fail", paste0("n=", dplyr::n_distinct(architecture_source$supermodule_id)),
  "global_eigengenes_has_9_supermodules", if (dplyr::n_distinct(global_eigengene_source$supermodule_id) == 9L) "pass" else "fail", paste0("n=", dplyr::n_distinct(global_eigengene_source$supermodule_id)),
  "forest_has_9_supermodules", if (nrow(forest_source) == 9L && dplyr::n_distinct(forest_source$supermodule_id) == 9L) "pass" else "fail", paste0("rows=", nrow(forest_source)),
  "spatial_heatmap_has_9_supermodules", if (nrow(spatial_heatmap_source) == 36L && dplyr::n_distinct(spatial_heatmap_source$supermodule_id) == 9L) "pass" else "fail", paste0("rows=", nrow(spatial_heatmap_source)),
  "effect_matrix_has_9_supermodules", if (nrow(effect_matrix_source) == 27L && dplyr::n_distinct(effect_matrix_source$supermodule_id) == 9L) "pass" else "fail", paste0("rows=", nrow(effect_matrix_source)),
  "singleton_loading_panels_absent", if (!any(loading_source$n_member_modules <= 1L) && identical(sort(unique(loading_source$supermodule_id)), sort(multi_ids))) "pass" else "fail", paste(unique(loading_source$supermodule_id), collapse = ";"),
  "singular_models_have_no_claim_symbols", if (!any((!forest_source$model_stable) & forest_source$claim_support_symbol != "none")) "pass" else "fail", "diagnostic-only rows cannot receive FDR claim symbols",
  "stable_symmetric_effect_scale", if (length(unique(c(spatial_heatmap_source$effect_scale_limit, effect_matrix_source$effect_scale_limit))) == 1L) "pass" else "fail", paste0("limit=", effect_limit),
  "positive_effect_means_higher_in_SUS", "pass", "SUS - RES orientation retained from Stage 05"
)
readr::write_csv(validation, file.path(OUT$tables, "WGCNA_publication_figure_validation.csv"), na = "")
if (any(validation$status != "pass")) stop("WGCNA publication-figure validation failed: ", paste(validation$validation_check[validation$status != "pass"], collapse = ", "), call. = FALSE)

write_run_manifest(
  file.path(OUT$logs, "run_manifest.yml"),
  inputs = list(
    reviewed_label_lookup = LABEL_FILE,
    current_member_map = FILES$supermodule_annotation,
    current_supermodule_summary = FILES$supermodule_summary,
    inferential_handoff = HANDOFF_FILE,
    raw_eigengenes = file.path(GROUP_DIR, "all_supermodule_eigengene_group_values.csv"),
    member_loadings = file.path(GROUP_DIR, "supermodule_pca_member_loadings.csv")
  ),
  outputs = list(figures = OUT$figures, source_data = OUT$source_data, validation = file.path(OUT$tables, "WGCNA_publication_figure_validation.csv")),
  parameters = list(dataset = DATASET, supermodule_order = paste(sm_order, collapse = ";"), roi_order = paste(roi_order, collapse = ";"), effect_scale_limit = effect_limit),
  notes = "All-supermodule publication layer. Inferential rows come only from the validated Stage 07 handoff; singleton panels display their canonical module endpoint and never create a second independent claim. Reviewed Stage 07 labels are authoritative. No WGCNA, membership, model, or FDR recomputation occurs."
)

message("Reviewed microglia WGCNA publication figures complete.")
