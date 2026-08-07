#!/usr/bin/env Rscript
# Corrected additive microglia WGCNA publication layer.  It uses immutable
# Stage 05 statistics and Stage 12 audit outputs; historic Stage 08 figures
# are intentionally not overwritten.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))
source(repo_path("R", "wgcna_reviewed_label_registry.R"))
pkgs <- c("dplyr", "readr", "tidyr", "ggplot2", "svglite", "scales", "stringr")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) && !is_dry_run()) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
if (!length(missing)) suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
run <- wgcna_cli(); if (!identical(run$dataset, "microglia")) stop("08b is microglia-only", call. = FALSE)
DATASET <- "microglia"; sm_order <- sprintf("SM%02d", 1:9); roi_order <- c("CA1", "CA2", "CA3", "DG")
OUT <- create_module_dirs("06_modules_WGCNA", file.path("wgcna_publication_figures_corrected", DATASET))
FILES <- resolve_wgcna_files(DATASET); TBL <- path_results("tables", "06_modules_WGCNA")
AUDIT <- path_results("reviewer_audit", "microglia_wgcna_nature_readiness")
paths <- list(lookup = file.path(TBL, "interpretable_summary", DATASET, "WGCNA_final_label_lookup.csv"),
  values = file.path(TBL, "group_effects", DATASET, "all_supermodule_eigengene_group_values.csv"),
  effects = file.path(TBL, "interpretable_summary", DATASET, "WGCNA_inferential_handoff.csv"),
  member = FILES$supermodule_annotation, summary = FILES$supermodule_summary,
  definitions = FILES$definitions,
  loadings = file.path(TBL, "group_effects", DATASET, "supermodule_pca_member_loadings.csv"),
  modules = file.path(AUDIT, "module_robustness_consensus.csv"), blocks = file.path(AUDIT, "higher_order_block_readiness_summary.csv"),
  variance = file.path(AUDIT, "module_eigengene_variance_partition_fixed_group.csv"))
if (run$dry_run) { for (p in unlist(paths)) dry_run_line("Input", p, if (file.exists(p)) "PASS" else "FAIL"); dry_run_line("Output", file.path(OUT$figures, "wgcna_readiness_summary.svg")); quit(status = 0) }
if (any(!file.exists(unlist(paths)))) stop("Missing corrected-figure input(s): ", paste(unlist(paths)[!file.exists(unlist(paths))], collapse = ", "), call. = FALSE)
read_csv <- function(p) readr::read_csv(p, show_col_types = FALSE, progress = FALSE)
save <- function(p, stem, w, h) { for (ext in c("svg", "pdf")) ggplot2::ggsave(file.path(OUT$figures, paste0(stem, ".", ext)), p, width = w, height = h, units = "in", device = if (ext == "svg") svglite::svglite else grDevices::cairo_pdf) }
write_source <- function(x, stem) readr::write_csv(x, file.path(OUT$source_data, paste0(stem, "_source.csv")), na = "")
theme_pub <- function(base = 7.5) ggplot2::theme_minimal(base_family = "Arial", base_size = base) + ggplot2::theme(axis.line = ggplot2::element_line(linewidth = .45, colour = "#333333"), panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = base))
lookup <- read_csv(paths$lookup); member <- wgcna_normalize_current_member_map(read_csv(paths$member), DATASET); wgcna_validate_canonical_lookup(lookup, DATASET, member)
definitions <- read_csv(paths$definitions) |> dplyr::transmute(ModuleID = as.character(.data$ModuleID), ModuleColor = as.character(.data$ModuleColor)) |> dplyr::distinct(.data$ModuleID, .keep_all = TRUE)
super_summary <- read_csv(paths$summary) |> dplyr::transmute(supermodule_id = as.character(.data$SupermoduleID), pc1_variance = as.numeric(.data$pc1_variance_explained)) |> dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
loadings <- read_csv(paths$loadings)
labels <- lookup |> dplyr::filter(.data$level == "supermodule") |> dplyr::transmute(supermodule_id = .data$entity_id, canonical_short_label, canonical_plot_label, structural_status, biological_label_confidence, roi_context)
join_labels <- function(x) { n <- nrow(x); label_columns <- setdiff(names(labels), "supermodule_id"); z <- x |> dplyr::select(-dplyr::any_of(label_columns)) |> dplyr::left_join(labels, by = "supermodule_id", relationship = "many-to-one"); if (nrow(z) != n) stop("Label join multiplied rows"); z }
values <- read_csv(paths$values) |> dplyr::mutate(supermodule_id = as.character(.data$supermodule_id), roi = factor(toupper(.data$canonical_spatial_unit), levels = roi_order), StressGroup = factor(.data$StressGroup, levels = c("CON", "RES", "SUS"))) |> join_labels()
inferential_handoff <- wgcna_inferential_handoff_read(paths$effects)
display_lookup <- lookup |>
  dplyr::filter(.data$level == "module") |>
  dplyr::transmute(
    dataset, module_id = .data$entity_id,
    supermodule_id = .data$parent_entity_id
  ) |>
  dplyr::left_join(
    lookup |>
      dplyr::filter(.data$level == "supermodule") |>
      dplyr::transmute(
        dataset, supermodule_id = .data$entity_id,
        n_member_modules,
        supermodule_label = canonical_plot_label
      ),
    by = c("dataset", "supermodule_id"),
    relationship = "many-to-one"
  )
effects <- wgcna_inferential_handoff_supermodule_display(
  inferential_handoff, display_lookup
) |>
  dplyr::mutate(
    supermodule_id = as.character(.data$supermodule_id),
    estimate = as.numeric(.data$estimate),
    raw_SE = as.numeric(.data$SE),
    q_value = as.numeric(.data$tier_specific_fdr),
    stable = .data$claim_gate == "eligible_for_readiness_assessment"
  ) |>
  join_labels()
response_sd_global <- values |> dplyr::group_by(.data$supermodule_id) |> dplyr::summarise(response_SD_global = stats::sd(.data$eigengene), model_analysis_N_global = dplyr::n(), animal_N_global = dplyr::n_distinct(.data$AnimalID), .groups = "drop")
response_sd_spatial <- values |> dplyr::mutate(spatial_unit = tolower(as.character(.data$roi))) |> dplyr::group_by(.data$supermodule_id, .data$spatial_unit) |> dplyr::summarise(response_SD_spatial = stats::sd(.data$eigengene), model_analysis_N_spatial = dplyr::n(), animal_N_spatial = dplyr::n_distinct(.data$AnimalID), .groups = "drop")
effects_std <- effects |> dplyr::mutate(spatial_unit = tolower(as.character(.data$spatial_unit))) |> dplyr::left_join(response_sd_global, by = "supermodule_id", relationship = "many-to-one") |> dplyr::left_join(response_sd_spatial, by = c("supermodule_id", "spatial_unit"), relationship = "many-to-one") |> dplyr::mutate(response_SD = dplyr::if_else(.data$effect_scope == "within_spatial_unit", .data$response_SD_spatial, .data$response_SD_global), model_analysis_N = dplyr::if_else(.data$effect_scope == "within_spatial_unit", .data$model_analysis_N_spatial, .data$model_analysis_N_global), animal_N = dplyr::if_else(.data$effect_scope == "within_spatial_unit", .data$animal_N_spatial, .data$animal_N_global), raw_CI_low = .data$estimate - 1.96 * .data$raw_SE, raw_CI_high = .data$estimate + 1.96 * .data$raw_SE, standardized_estimate = .data$estimate / .data$response_SD, standardized_SE = .data$raw_SE / .data$response_SD, standardized_CI_low = .data$raw_CI_low / .data$response_SD, standardized_CI_high = .data$raw_CI_high / .data$response_SD, standardization_scope = "SD of exact Stage 05 eigengene response analysis subset")
readr::write_csv(effects_std, file.path(OUT$tables, "supermodule_group_effects_standardized.csv"), na = "")

# Architecture: blocks and singleton compatibility identities are visually separated.
architecture <- member |> dplyr::left_join(definitions, by = "ModuleID", relationship = "one-to-one") |> dplyr::rename(supermodule_id = "SupermoduleID") |> join_labels() |> dplyr::left_join(super_summary, by = "supermodule_id", relationship = "many-to-one") |> dplyr::group_by(.data$supermodule_id) |> dplyr::arrange(.data$ModuleID, .by_group = TRUE) |> dplyr::mutate(member_count = dplyr::n(), member_position = dplyr::row_number()) |> dplyr::ungroup() |> dplyr::mutate(structure_class = ifelse(.data$member_count > 1L, "multi-module higher-order block", "standalone module / singleton compatibility ID"), singleton_badge = ifelse(.data$structural_status == "singleton", "singleton", ""), ModuleColor = ifelse(grepl("^#[0-9A-Fa-f]{6}$", .data$ModuleColor), .data$ModuleColor, "#BDBDBD"), y = factor(paste(.data$supermodule_id, .data$canonical_short_label, sep = " | "), levels = rev(paste(sm_order, labels$canonical_short_label[match(sm_order, labels$supermodule_id)], sep = " | "))))
write_source(architecture, "corrected_all_supermodule_architecture")
p_arch_text <- architecture |> dplyr::distinct(.data$supermodule_id, .data$y, .data$member_count, .data$pc1_variance, .data$structural_status, .data$biological_label_confidence, .data$roi_context) |> dplyr::mutate(annotation = dplyr::if_else(.data$structural_status == "singleton", paste0("standalone; ", .data$biological_label_confidence), paste0("block n=", .data$member_count, "; PC1=", scales::percent(.data$pc1_variance, accuracy = 1), "; ", .data$biological_label_confidence)))
p_arch <- ggplot2::ggplot(architecture, ggplot2::aes(x = .data$member_position, y = .data$y)) + ggplot2::geom_tile(ggplot2::aes(fill = .data$ModuleColor), width = .85, height = .66, colour = "#333333", linewidth = .16) + ggplot2::geom_text(ggplot2::aes(label = .data$ModuleID), size = 2.25, family = "Arial", colour = "white") + ggplot2::geom_text(data = p_arch_text, ggplot2::aes(x = 4.0, y = .data$y, label = .data$annotation), inherit.aes = FALSE, hjust = 0, vjust = -.3, size = 2.2, family = "Arial") + ggplot2::geom_text(data = p_arch_text, ggplot2::aes(x = 4.0, y = .data$y, label = stringr::str_trunc(.data$roi_context, 52)), inherit.aes = FALSE, hjust = 0, vjust = 1.25, size = 2.2, family = "Arial", colour = "#666666") + ggplot2::scale_fill_identity() + ggplot2::scale_x_continuous(limits = c(.5, 8.4), breaks = NULL) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(x = NULL, y = NULL) + theme_pub() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(3, 60, 3, 3))
save(p_arch, "corrected_all_supermodule_architecture", 8.2, 4.3)

# Global panels: ROI points are faint; one prominent mean per AnimalID is the biological-unit summary.
animal_means <- values |> dplyr::group_by(.data$supermodule_id, .data$canonical_plot_label, .data$StressGroup, .data$AnimalID) |> dplyr::summarise(animal_mean_eigengene = mean(.data$eigengene), n_roi = dplyr::n(), .groups = "drop")
global_source <- values |> dplyr::select(supermodule_id, canonical_plot_label, StressGroup, AnimalID, roi, roi_eigengene = eigengene) |> dplyr::left_join(animal_means, by = c("supermodule_id", "canonical_plot_label", "StressGroup", "AnimalID"), relationship = "many-to-one")
write_source(global_source, "corrected_all_supermodule_global_eigengenes")
p_global <- ggplot2::ggplot(global_source, ggplot2::aes(.data$StressGroup, .data$roi_eigengene, colour = .data$StressGroup)) + ggplot2::geom_point(position = ggplot2::position_jitter(width = .12, height = 0, seed = 1), alpha = .22, size = .55) + ggplot2::geom_point(data = animal_means, ggplot2::aes(y = .data$animal_mean_eigengene), inherit.aes = TRUE, shape = 21, fill = "white", size = 1.75, stroke = .35) + ggplot2::facet_wrap(~canonical_plot_label, ncol = 3, scales = "free_y") + ggplot2::scale_colour_manual(values = c(CON = "#3E3C6F", RES = "#9E9A92", SUS = "#D7303F")) + ggplot2::labs(x = NULL, y = "Eigengene; faint ROI observations, n=3 animal means/group") + theme_pub(7.2) + ggplot2::theme(legend.position = "none")
save(p_global, "corrected_all_supermodule_global_eigengenes", 7.2, 6.1)

global <- effects_std |> dplyr::filter(.data$effect_scope == "spatial_adjusted_global", .data$contrast == "SUS - RES") |> dplyr::mutate(plot_label = factor(.data$canonical_plot_label, levels = rev(labels$canonical_plot_label)), claim_symbol = dplyr::case_when(!.data$stable ~ "none", !(.data$display_is_independent_endpoint %in% TRUE) ~ "none", !(.data$independent_hypothesis %in% TRUE) ~ "none", .data$q_value <= .05 ~ "FDR05", .data$q_value <= .10 ~ "FDR10", TRUE ~ "none"))
write_source(global, "corrected_all_supermodule_sus_res_forest")
p_forest <- ggplot2::ggplot(global, ggplot2::aes(.data$standardized_estimate, .data$plot_label)) + ggplot2::geom_vline(xintercept = 0, linewidth = .4) + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$standardized_CI_low, xmax = .data$standardized_CI_high), orientation = "y", linewidth = .35, width = .16) + ggplot2::geom_point(ggplot2::aes(shape = .data$stable), size = 1.8) + ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), labels = c(`TRUE` = "stable", `FALSE` = "singular/unstable diagnostic")) + ggplot2::labs(x = "SUS minus RES (response-SD units, 95% CI)", y = NULL, shape = NULL) + theme_pub()
save(p_forest, "corrected_all_supermodule_sus_res_forest", 6.6, 4.2)

spatial <- effects_std |> dplyr::filter(.data$effect_scope == "within_spatial_unit", .data$contrast == "SUS - RES") |> dplyr::mutate(roi = factor(toupper(.data$spatial_unit), levels = roi_order), plot_label = factor(.data$canonical_plot_label, levels = rev(labels$canonical_plot_label)), support = dplyr::case_when(!.data$stable ~ "cross", !(.data$display_is_independent_endpoint %in% TRUE) ~ "none", !(.data$independent_hypothesis %in% TRUE) ~ "none", .data$q_value <= .05 ~ "filled", .data$q_value <= .10 ~ "open", TRUE ~ "none"))
scale_limit <- max(abs(c(spatial$standardized_estimate, global$standardized_estimate)), na.rm = TRUE); if (!is.finite(scale_limit)) scale_limit <- 1
write_source(spatial, "corrected_all_supermodule_sus_res_spatial_heatmap")
roi_side <- spatial |> dplyr::distinct(.data$plot_label, .data$roi_context)
p_spatial <- ggplot2::ggplot(spatial, ggplot2::aes(.data$roi, .data$plot_label, fill = .data$standardized_estimate)) + ggplot2::geom_tile(colour = "white", linewidth = .25) + ggplot2::geom_point(data = spatial |> dplyr::filter(.data$support == "filled"), shape = 16, size = 1.5) + ggplot2::geom_point(data = spatial |> dplyr::filter(.data$support == "open"), shape = 1, size = 1.9) + ggplot2::geom_point(data = spatial |> dplyr::filter(.data$support == "cross"), shape = 4, size = 2, colour = "grey45") + ggplot2::geom_text(data = roi_side, ggplot2::aes(x = 4.6, y = .data$plot_label, label = stringr::str_trunc(.data$roi_context, 42)), inherit.aes = FALSE, hjust = 0, size = 2.1, family = "Arial", colour = "#666666") + ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84C5A", midpoint = 0, limits = c(-scale_limit, scale_limit), name = "SUS - RES\nresponse SD") + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(.05, 1.6))) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(x = NULL, y = NULL) + theme_pub() + ggplot2::theme(plot.margin = ggplot2::margin(3, 50, 3, 3))
save(p_spatial, "corrected_all_supermodule_sus_res_spatial_heatmap", 6.2, 4.3)

# Only genuine multi-module blocks receive member-loading panels.
module_labels <- lookup |> dplyr::filter(.data$level == "module") |> dplyr::transmute(module_id = .data$entity_id, canonical_module_short_label = .data$canonical_short_label)
loading_source <- loadings |> dplyr::filter(.data$supermodule_id %in% c("SM01", "SM03", "SM09"), as.integer(.data$pc) == 1L) |> dplyr::transmute(supermodule_id = as.character(.data$supermodule_id), module_id = as.character(.data$module_id), loading = as.numeric(.data$loading), loading_rank = as.integer(.data$loading_rank), n_member_modules = as.integer(.data$n_member_modules)) |> dplyr::left_join(module_labels, by = "module_id", relationship = "many-to-one") |> join_labels() |> dplyr::mutate(panel = factor(paste(.data$supermodule_id, .data$canonical_short_label, sep = " | "), levels = paste(c("SM01", "SM03", "SM09"), labels$canonical_short_label[match(c("SM01", "SM03", "SM09"), labels$supermodule_id)], sep = " | ")), module_label = factor(.data$canonical_module_short_label, levels = rev(unique(.data$canonical_module_short_label))))
if (any(loading_source$n_member_modules <= 1L)) stop("Singleton loading panels are forbidden", call. = FALSE)
write_source(loading_source, "corrected_multi_supermodule_member_loadings")
p_load <- ggplot2::ggplot(loading_source, ggplot2::aes(.data$loading, .data$module_label)) + ggplot2::geom_vline(xintercept = 0, linewidth = .35) + ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$loading, yend = .data$module_label), linewidth = .35, colour = "#5F7590") + ggplot2::geom_point(size = 1.8, colour = "#5F7590") + ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 3) + ggplot2::labs(x = "PC1 loading", y = NULL) + theme_pub()
save(p_load, "corrected_multi_supermodule_member_loadings", 7.0, 2.8)

matrix <- effects_std |> dplyr::filter(.data$effect_scope == "spatial_adjusted_global", .data$contrast %in% c("RES - CON", "SUS - CON", "SUS - RES")) |> dplyr::mutate(contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "SUS - RES")), plot_label = factor(.data$canonical_plot_label, levels = rev(labels$canonical_plot_label)))
write_source(matrix, "corrected_supplementary_all_contrast_effect_matrix")
p_matrix <- ggplot2::ggplot(matrix, ggplot2::aes(.data$contrast, .data$plot_label, fill = .data$standardized_estimate)) + ggplot2::geom_tile(colour = "white", linewidth = .25) + ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84C5A", midpoint = 0, limits = c(-scale_limit, scale_limit), name = "Response SD") + ggplot2::labs(x = NULL, y = NULL) + theme_pub()
save(p_matrix, "corrected_supplementary_all_contrast_effect_matrix", 5.8, 4.3)

# Audit-derived readiness figure is intentionally separate from Stage 08 effects.
mod <- read_csv(paths$modules); block <- read_csv(paths$blocks); vp <- read_csv(paths$variance) |> dplyr::filter(.data$component %in% c("StressGroup_fixed", "Region_fixed", "AnimalID_random", "Residuals"))
readiness_source <- mod |> dplyr::transmute(ModuleID, bilateral = .data$bilateral_reproducibility_status, spatial_adjusted = .data$spatial_adjusted_robustness_status, animal = .data$animal_stability_status, leave_one_region_out = .data$leave_one_region_out_status, strict_nonspatial = .data$strict_nonspatial_sensitivity_status, primary_architecture_status, spatial_dependence_class) |> dplyr::left_join(lookup |> dplyr::filter(.data$level == "module") |> dplyr::transmute(ModuleID = .data$entity_id, roi_context), by = "ModuleID", relationship = "one-to-one")
write_source(readiness_source, "wgcna_readiness_summary")
readiness_long <- readiness_source |> tidyr::pivot_longer(c("bilateral", "spatial_adjusted", "animal", "leave_one_region_out", "strict_nonspatial"), names_to = "check", values_to = "status") |> dplyr::mutate(score = dplyr::case_when(.data$status %in% c("pass", "hippocampus_wide_stable") ~ 2, .data$status %in% c("suggestive", "partially_region_dependent") ~ 1, TRUE ~ 0), ModuleID = factor(.data$ModuleID, levels = rev(sprintf("WGCNA_m%02d", 1:13))))
p_ready <- ggplot2::ggplot(readiness_long, ggplot2::aes(.data$check, .data$ModuleID, fill = .data$score)) + ggplot2::geom_tile(colour = "white", linewidth = .25) + ggplot2::scale_fill_gradient(low = "#D9D9D9", high = "#2166AC", breaks = 0:2, labels = c("fail/diagnostic", "suggestive", "pass"), name = "Robustness") + ggplot2::facet_grid(. ~ "Corrected fixed-membership readiness") + ggplot2::labs(x = NULL, y = NULL) + theme_pub()
save(p_ready, "wgcna_readiness_summary", 7.5, 5.0)
validation <- tibble::tribble(~check, ~status, ~detail,
 "nine_supermodules_in_architecture", if (dplyr::n_distinct(architecture$supermodule_id) == 9) "pass" else "fail", as.character(dplyr::n_distinct(architecture$supermodule_id)),
 "nine_animal_means_per_supermodule", if (all(table(animal_means$supermodule_id) == 9)) "pass" else "fail", "3 groups x 3 animals",
 "raw_and_standardized_effects", if (all(c("estimate", "raw_SE", "response_SD", "standardized_estimate", "standardized_SE", "standardized_CI_low", "standardized_CI_high") %in% names(effects_std))) "pass" else "fail", "denominators retained",
 "no_singular_claim_symbols", if (!any(global$claim_symbol != "none" & !global$stable)) "pass" else "fail", "diagnostic rows unpromoted",
 "strict_nonspatial_separate", if ("strict_nonspatial" %in% unique(readiness_long$check)) "pass" else "fail", "shown separately")
readr::write_csv(validation, file.path(OUT$tables, "WGCNA_corrected_publication_figure_validation.csv")); if (any(validation$status != "pass")) stop("Corrected figure validation failed")
write_run_manifest(file.path(OUT$logs, "run_manifest.yml"), inputs = paths, outputs = list(figures = OUT$figures, source_data = OUT$source_data), parameters = list(dataset = DATASET, effect_standardization = "estimate / SD(exact Stage 05 eigengene response analysis subset)"), notes = "Additive corrected publication layer; inference comes only from the Stage 07 handoff. Singleton panels display their canonical module endpoint without a duplicate independent claim. Existing primary publication figures were not modified.")
message("Corrected additive microglia WGCNA figures complete.")
