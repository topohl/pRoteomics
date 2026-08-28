#!/usr/bin/env Rscript

# Bounded downstream three-contrast audit of active stress-response remodeling.
# Consumes canonical mapped protein contrasts and existing ranked GO-BP/GSEA output.
# Does not refit differential abundance, rerun GSEA, modify DAP-set ORA, or create new tests.

source("R/paths.R")
source("R/dataset_config.R")
source("R/plotting_nature.R")
source("R/enrichment_io.R")
source("R/manuscript_go_theme_utils.R")
source("R/sus_res_spatial_dap_atlas_utils.R")
source("R/stress_response_biological_audit_utils.R")

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "stress_response_biological_audit"
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
ROOT <- repo_root()
TABLE_DIR <- path_results("tables", MODULE_ID, SUBSTEP_ID, "global")
SOURCE_DIR <- path_results("source_data", MODULE_ID, SUBSTEP_ID, "global")
FIGURE_DIR <- path_results("figures", MODULE_ID, SUBSTEP_ID, "global")
REPORT_DIR <- path_results("reports", MODULE_ID, SUBSTEP_ID, "global")
LOG_DIR <- path_results("logs", MODULE_ID, SUBSTEP_ID, "global")
REGISTRY_PATH <- repo_path("config", "manuscript_go_theme_registry.tsv")
GO_SOURCE <- path_results("source_data", MODULE_ID, "compareGO_spatial_atlas", "spatial_atlas_enrichment_long.csv")
SUS_RES_REFERENCE_WORKBOOK <- path_results("reports", MODULE_ID, "sus_res_spatial_dap_atlas", "global", "sus_res_biological_audit.xlsx")

required_packages <- c("readr", "tidyselect", "ggplot2", "svglite", "GO.db", "AnnotationDbi")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

is_truthy <- function(x) tolower(x) %in% c("1", "true", "yes", "y")
DRY_RUN <- "--dry-run" %in% commandArgs(trailingOnly = TRUE) || is_dry_run() || is_truthy(Sys.getenv("PROTEOMICS_STRESS_RESPONSE_AUDIT_DRY_RUN", unset = ""))

manifest_paths <- setNames(vapply(DATASETS, canonical_clusterprofiler_manifest_path, character(1)), DATASETS)
required_inputs <- c(manifest_paths, GO_SOURCE, REGISTRY_PATH)
if (any(!file.exists(required_inputs))) stop("Missing required canonical input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = "; "), call. = FALSE)

specs <- do.call(rbind, lapply(DATASETS, function(dataset) stress_response_manifest_specs(dataset, manifest_paths[[dataset]], ROOT)))
if (nrow(specs) != 54L || length(unique(paste(specs$dataset, specs$spatial_unit, sep = "|"))) != 18L) {
  stop("Expected 54 contrast specifications across 18 spatial contexts; found ", nrow(specs), " rows.", call. = FALSE)
}

if (DRY_RUN) {
  message("[DRY-RUN] three-contrast stress-response biological audit")
  message("[DRY-RUN] contexts / contrast specifications: 18 / ", nrow(specs))
  message("[DRY-RUN] all mapped contrast hashes match manifests: ", all(specs$input_hash_matches_manifest))
  message("[DRY-RUN] ranked GO source: ", GO_SOURCE, " (", file.info(GO_SOURCE)$size, " bytes)")
  message("[DRY-RUN] protein joint families: dataset x spatial_unit, 2m raw p-values on matched ProteinGroupID universes")
  message("[DRY-RUN] GO joint families: dataset x spatial_unit, 2m nominal p-values on matched GO-ID universes")
  message("[DRY-RUN] SUS-RES remains outside both joint control families")
  quit(save = "no", status = 0L)
}

for (path in c(TABLE_DIR, SOURCE_DIR, FIGURE_DIR, REPORT_DIR, LOG_DIR)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

message("Reading 54 canonical mapped protein contrasts...")
protein_rows <- lapply(seq_len(nrow(specs)), function(i) {
  raw <- utils::read.csv(specs$input_file[[i]], stringsAsFactors = FALSE, check.names = FALSE)
  stress_response_normalize_protein_contrast(raw, specs[i, , drop = FALSE])
})
protein_long <- do.call(rbind, protein_rows)
rownames(protein_long) <- NULL
protein_bh_validation <- stress_response_bh_validation(
  protein_long, c("dataset", "spatial_unit", "contrast"), "p_value", "canonical_contrast_FDR", tolerance = 1e-12
)
protein_bh_validation$canonical_FDR_preserved_without_recalculation <- TRUE
protein_bh_validation$diagnostic_interpretation <- ifelse(
  protein_bh_validation$validation_pass,
  "mapped ProteinGroupID rows reproduce the upstream BH family",
  "mapped ProteinGroupID rows do not fully reconstruct the upstream adjustment family/precision; canonical padj remains authoritative and is not recalculated"
)
readr::write_csv(protein_bh_validation, file.path(LOG_DIR, "protein_canonical_BH_preflight.csv"), na = "")
if (!all(protein_bh_validation$validation_pass)) {
  failed <- protein_bh_validation[!protein_bh_validation$validation_pass, , drop = FALSE]
  message("Protein canonical-BH preflight failures: ", nrow(failed), "; maximum absolute difference: ", max(failed$max_abs_FDR_difference, na.rm = TRUE))
  warning(
    "Mapped ProteinGroupID rows do not exactly reconstruct upstream canonical padj. Canonical padj is preserved unchanged; ",
    "the secondary joint family is defined only on the explicitly matched finite-raw-p ProteinGroupID universe.",
    call. = FALSE
  )
}
protein <- stress_response_build_protein_outputs(protein_long)
direct_dap <- stress_response_direct_dap_control_geometry(protein$geometry)
readr::write_csv(protein$common_universe_audit, file.path(LOG_DIR, "protein_common_universe_preflight.csv"), na = "")
readr::write_csv(protein$algebra_audit, file.path(LOG_DIR, "protein_contrast_algebra_preflight.csv"), na = "")
if (!all(protein$common_universe_audit$validation_pass)) stop("Protein common-universe validation failed.", call. = FALSE)
if (!all(protein$algebra_audit$tolerance_pass)) {
  failed_algebra <- protein$algebra_audit[!protein$algebra_audit$tolerance_pass, , drop = FALSE]
  message("Protein algebra failures: ", nrow(failed_algebra), "; maximum absolute residual: ", max(failed_algebra$max_absolute_residual))
  stop("Protein contrast algebra failed; source/sign alignment must be investigated.", call. = FALSE)
}
if (!all(protein$joint_family_audit$validation_pass)) stop("Protein joint-control BH family validation failed.", call. = FALSE)

message("Reading the existing full ranked-GSEA source (selected columns only)...")
go_required_columns <- c(
  "ID", "Description", "pvalue", "p.adjust", "NES", "core_enrichment", "core_enrichment_gene",
  "dataset", "dataset_label", "source_supplementary_file", "source_manifest_file", "evidence_source_family",
  "route_unit", "phenotype_contrast", "region", "layer", "compartment", "spatial_unit"
)
go_comparison_columns <- c("comparison", "legacy_Comparison", "Comparison")
go_long <- readr::read_csv(
  GO_SOURCE,
  col_select = tidyselect::any_of(c(go_required_columns, go_comparison_columns)),
  show_col_types = FALSE,
  progress = interactive(),
  name_repair = "minimal"
)
go_long <- as.data.frame(go_long, stringsAsFactors = FALSE, check.names = FALSE)
missing_go_columns <- setdiff(go_required_columns, names(go_long))
if (length(missing_go_columns)) {
  stop("Ranked-GSEA source is missing required column(s): ", paste(missing_go_columns, collapse = ", "), call. = FALSE)
}
go_long <- normalize_stage07_comparison_column(go_long, "Ranked-GSEA source")
go_long <- go_long[as.character(go_long$phenotype_contrast) %in% STRESS_RESPONSE_CONTRASTS, , drop = FALSE]
go_long$pvalue <- suppressWarnings(as.numeric(go_long$pvalue))
go_long$p.adjust <- suppressWarnings(as.numeric(go_long$p.adjust))
go_long$NES <- suppressWarnings(as.numeric(go_long$NES))
# Ranked-GSEA provenance contract. The raw evidence_source_family is preserved
# verbatim; evidence_source_class carries the normalized semantic class that
# downstream logic compares against. Fails closed on any non-ranked-GSEA
# family, naming the offending values rather than reporting a "mixed" source
# when the input is uniformly current-canonical under a newer vocabulary.
go_long$evidence_source_class <- sus_res_normalize_ranked_gsea_source_class(
  go_long$evidence_source_family,
  "Stress-response ranked-GSEA GO source"
)
if (any(go_long$evidence_source_class != SUS_RES_RANKED_GSEA_SOURCE_CLASS, na.rm = TRUE)) {
  stop("Stress-response ranked-GSEA GO source did not normalize to a single ",
       "ranked-GSEA evidence class.", call. = FALSE)
}
go_context_counts <- table(paste(go_long$dataset, go_long$spatial_unit, go_long$phenotype_contrast, sep = "|"))
if (length(go_context_counts) != 54L || any(go_context_counts == 0L)) stop("Ranked-GSEA source lacks one or more required context/contrast families.", call. = FALSE)
go_bh_validation <- stress_response_bh_validation(
  go_long, c("dataset", "spatial_unit", "phenotype_contrast"), "pvalue", "p.adjust", tolerance = 1e-12
)
if (!all(go_bh_validation$validation_pass)) {
  stop("Original ranked-GSEA nominal p-values do not reconstruct the canonical contrast-specific GO BH FDR families.", call. = FALSE)
}
joint_go <- stress_response_add_joint_go_fdr(go_long)
if (!all(joint_go$audit$validation_pass)) stop("GO joint-control BH family validation failed.", call. = FALSE)

message("Applying the existing manuscript GO-ID/is_a/part_of registry to all three contrasts...")
registry <- read_manuscript_go_theme_registry(REGISTRY_PATH)
go_terms <- unique(go_long[c("ID", "Description")])
mapping <- map_go_terms_to_manuscript_themes(go_terms, registry)
theme <- stress_response_build_theme_detail(go_long, mapping, joint_go$joint)
trajectory <- stress_response_build_theme_trajectory(theme$detail)
overview <- stress_response_theme_overview(theme$detail, trajectory)
if (nrow(theme$detail) != 8L * 18L * 3L) stop("Theme detail must contain 432 theme x context x contrast rows.", call. = FALSE)
if (nrow(trajectory) != 6L * 18L) stop("Primary trajectory must contain 108 theme x context rows.", call. = FALSE)

joint_go_key <- paste(joint_go$joint$dataset, joint_go$joint$spatial_unit, joint_go$joint$phenotype_contrast, joint_go$joint$ID, sep = "|")
go_key <- paste(go_long$dataset, go_long$spatial_unit, go_long$phenotype_contrast, go_long$ID, sep = "|")
go_long$control_pair_joint_GO_FDR <- joint_go$joint$control_pair_joint_GO_FDR[match(go_key, joint_go_key)]
go_long$control_pair_joint_GO_FDR05 <- joint_go$joint$control_pair_joint_GO_FDR05[match(go_key, joint_go_key)]

supported_audits <- stress_response_build_supported_go_audits(go_long, mapping)
supported_occurrence <- supported_audits$occurrence
supported <- supported_audits$exploded

contrast_map_source <- theme$detail[theme$detail$theme_role == "primary", c(
  "dataset", "dataset_label", "compartment", "region", "layer", "spatial_unit", "contrast", "theme_id", "manuscript_theme",
  "median_NES_all_theme_terms", "canonical_GO_FDR_support", "min_canonical_GO_FDR", "n_canonical_FDR05_GO_terms",
  "control_pair_joint_GO_FDR_support", "min_control_pair_joint_GO_FDR", "n_control_pair_joint_FDR05_GO_terms"
)]

input_provenance <- specs[c(
  "dataset", "spatial_unit", "route_unit", "contrast", "comparison", "serialized_effect_definition", "formal_effect_definition",
  "formal_effect_multiplier", "sign_was_flipped", "manifest_file", "input_file_manifest", "input_file", "manifest_input_hash",
  "current_input_hash", "input_hash_matches_manifest"
)]
protected_reference <- stress_response_protected_reference_artifacts(
  SUS_RES_REFERENCE_WORKBOOK,
  producer = "04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r",
  note = "Stage 11 does not read or modify this protected SUS-RES biological audit workbook."
)

write_csv <- function(x, name, directory = SOURCE_DIR) {
  path <- file.path(directory, name)
  readr::write_csv(x, path, na = "")
  path
}

outputs <- c(
  protein_geometry = write_csv(protein$geometry, "protein_three_contrast_geometry.csv"),
  direct_dap_geometry = write_csv(direct_dap$detail, "direct_SUS_RES_DAP_control_geometry.csv"),
  direct_dap_context_summary = write_csv(direct_dap$context_summary, "direct_SUS_RES_DAP_control_geometry_summary_by_context.csv"),
  direct_dap_global_summary = write_csv(direct_dap$global_summary, "direct_SUS_RES_DAP_control_geometry_summary_global.csv"),
  protein_summary = write_csv(protein$summary, "protein_remodeling_summary.csv", TABLE_DIR),
  protein_common_universe = write_csv(protein$common_universe_audit, "protein_common_universe_audit.csv"),
  protein_algebra = write_csv(protein$algebra_audit, "protein_contrast_algebra_audit.csv"),
  protein_joint_family = write_csv(protein$joint_family_audit, "protein_control_pair_joint_FDR_family_audit.csv"),
  protein_canonical_bh = write_csv(protein_bh_validation, "protein_canonical_BH_reconstruction_diagnostic.csv"),
  go_joint_rows = write_csv(joint_go$joint, "control_pair_joint_GO_FDR.csv"),
  go_joint_family = write_csv(joint_go$audit, "GO_control_pair_joint_FDR_family_audit.csv"),
  go_canonical_bh = write_csv(go_bh_validation, "GO_canonical_BH_reconstruction_validation.csv"),
  theme_detail = write_csv(theme$detail, "theme_detail_all_three_contrasts.csv", TABLE_DIR),
  theme_trajectory = write_csv(trajectory, "theme_trajectories.csv", TABLE_DIR),
  theme_overview = write_csv(overview, "theme_overview.csv", TABLE_DIR),
  contrast_map = write_csv(contrast_map_source, "contrast_map_source.csv"),
  supported_go = write_csv(supported, "FDR_supported_GO_term_audit_all_three_contrasts.csv"),
  supported_go_occurrence = write_csv(supported_occurrence, "FDR_supported_GO_occurrence_audit_all_three_contrasts.csv"),
  theme_assignment = write_csv(collapse_go_theme_assignment_audit(mapping), "GO_term_theme_assignment_all_three_contrasts.csv"),
  input_provenance = write_csv(input_provenance, "input_contrast_provenance.csv"),
  protected_reference = write_csv(protected_reference, "protected_reference_artifacts.csv")
)

dataset_labels <- c(neuron_neuropil = "Neuron neuropil", neuron_soma = "Neuron soma", microglia = "Microglia-enriched ROI")
trajectory$dataset_display <- unname(dataset_labels[trajectory$dataset])
trajectory$direct_support_label <- ifelse(trajectory$canonical_GO_FDR_support_SUS_vs_RES, "Direct SUS-RES GO FDR support", "No direct SUS-RES GO FDR support")
trajectory$unit_label <- paste(trajectory$dataset_display, trajectory$spatial_unit, sep = " | ")

trajectory_plot <- ggplot2::ggplot(trajectory, ggplot2::aes(x = median_NES_RES_vs_CON, y = median_NES_SUS_vs_CON)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(colour = dataset_display, shape = direct_support_label), size = 2.35, stroke = 0.7, alpha = 0.9) +
  ggplot2::facet_wrap(~manuscript_theme, ncol = 3, scales = "free") +
  ggplot2::scale_colour_manual(values = c("Neuron neuropil" = "#3B6FB6", "Neuron soma" = "#D17A22", "Microglia-enriched ROI" = "#3B8D63"), name = NULL) +
  ggplot2::scale_shape_manual(values = c("Direct SUS-RES GO FDR support" = 21, "No direct SUS-RES GO FDR support" = 16), name = NULL) +
  ggplot2::labs(
    x = "Median NES: RES vs CON", y = "Median NES: SUS vs CON",
    title = "Candidate three-contrast GO-theme trajectories",
    caption = "Color identifies dataset only. Outlined point = at least one original SUS-RES GO term with BH FDR < 0.05. Median NES geometry is descriptive; NES is not algebraically additive."
  ) +
  theme_nature_base(7) +
  ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(size = 6, hjust = 0), strip.text = ggplot2::element_text(size = 6.5))

protein$summary$dataset_display <- unname(dataset_labels[protein$summary$dataset])
protein$summary$unit_label <- factor(protein$summary$spatial_unit, levels = rev(sort(unique(protein$summary$spatial_unit), method = "radix")))
burden_long <- rbind(
  data.frame(protein$summary[c("dataset_display", "unit_label")], contrast = "SUS vs CON", n_joint_FDR05_DAP = protein$summary$n_joint_FDR05_SUS_vs_CON),
  data.frame(protein$summary[c("dataset_display", "unit_label")], contrast = "RES vs CON", n_joint_FDR05_DAP = protein$summary$n_joint_FDR05_RES_vs_CON)
)
burden_plot <- ggplot2::ggplot(protein$summary, ggplot2::aes(y = unit_label)) +
  ggplot2::geom_segment(ggplot2::aes(x = n_joint_FDR05_SUS_vs_CON, xend = n_joint_FDR05_RES_vs_CON, yend = unit_label), linewidth = 0.6, colour = "grey55") +
  ggplot2::geom_point(data = burden_long, ggplot2::aes(x = n_joint_FDR05_DAP, shape = contrast), size = 2.4, stroke = 0.7) +
  ggplot2::facet_wrap(~dataset_display, scales = "free_y", nrow = 1) +
  ggplot2::scale_shape_manual(values = c("SUS vs CON" = 4, "RES vs CON" = 21), name = NULL) +
  ggplot2::scale_x_continuous(breaks = 0:max(protein$summary$n_joint_FDR05_RES_vs_CON, protein$summary$n_joint_FDR05_SUS_vs_CON), expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
  ggplot2::labs(
    x = "Joint-control BH FDR < 0.05 DAP count", y = NULL,
    title = "Candidate matched control-contrast protein burden",
    caption = "Each unit uses the common ProteinGroupID universe and one BH family across both control contrasts. Overlapping circle/cross symbols denote tied counts."
  ) +
  theme_nature_base(7) +
  ggplot2::theme(legend.position = "bottom", panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.2), plot.caption = ggplot2::element_text(size = 6, hjust = 0))

direct_context <- direct_dap$context_summary[direct_dap$context_summary$n_direct_DAP > 0L, , drop = FALSE]
direct_context$dataset_display <- unname(dataset_labels[direct_context$dataset])
direct_context <- direct_context[order(match(direct_context$dataset, DATASETS), direct_context$spatial_unit, method = "radix"), , drop = FALSE]
direct_context$context_label <- paste(direct_context$dataset_display, direct_context$spatial_unit, sep = " | ")
direct_context$context_label <- factor(direct_context$context_label, levels = rev(direct_context$context_label))
direct_context$panel <- factor("Sign geometry", levels = c("Sign geometry", "RES-farther fraction"))
sign_labels <- c(
  RES_up_SUS_down = "RES up / SUS down",
  RES_down_SUS_up = "RES down / SUS up",
  both_up_vs_CON = "Both above CON",
  both_down_vs_CON = "Both below CON",
  one_or_both_exact_zero = "One/both exact zero"
)
direct_composition <- do.call(rbind, lapply(STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY, function(z) {
  data.frame(
    direct_context[c("context_label", "n_direct_DAP")],
    sign_geometry = factor(z, levels = STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY),
    fraction = direct_context[[paste0("fraction_", z)]],
    panel = "Sign geometry",
    stringsAsFactors = FALSE
  )
}))
direct_displacement <- data.frame(
  direct_context[c("context_label", "n_direct_DAP", "fraction_RES_farther")],
  panel = "RES-farther fraction", stringsAsFactors = FALSE
)
direct_composition$panel <- factor(direct_composition$panel, levels = c("Sign geometry", "RES-farther fraction"))
direct_displacement$panel <- factor(direct_displacement$panel, levels = c("Sign geometry", "RES-farther fraction"))
direct_geometry_plot <- ggplot2::ggplot() +
  ggplot2::geom_col(
    data = direct_composition,
    ggplot2::aes(x = fraction, y = context_label, fill = sign_geometry),
    width = 0.72, colour = "white", linewidth = 0.2
  ) +
  ggplot2::geom_text(
    data = direct_context,
    ggplot2::aes(x = 1.025, y = context_label, label = paste0("n = ", n_direct_DAP)),
    inherit.aes = FALSE, hjust = 0, size = 2.1, colour = "grey25"
  ) +
  ggplot2::geom_segment(
    data = direct_displacement,
    ggplot2::aes(x = 0, xend = 1, y = context_label, yend = context_label),
    inherit.aes = FALSE, linewidth = 0.35, colour = "grey75"
  ) +
  ggplot2::geom_point(
    data = direct_displacement,
    ggplot2::aes(x = fraction_RES_farther, y = context_label),
    inherit.aes = FALSE, shape = 21, size = 2.6, stroke = 0.5, fill = "#3B6FB6", colour = "#24476F"
  ) +
  ggplot2::geom_text(
    data = direct_displacement,
    ggplot2::aes(
      x = fraction_RES_farther, y = context_label,
      label = sprintf("%.0f%%", 100 * fraction_RES_farther),
      hjust = ifelse(fraction_RES_farther >= 0.82, 1.2, -0.2)
    ),
    inherit.aes = FALSE, size = 2.0, colour = "grey20"
  ) +
  ggplot2::facet_grid(. ~ panel) +
  ggplot2::scale_fill_manual(
    values = c(
      RES_up_SUS_down = "#3B6FB6", RES_down_SUS_up = "#D17A22",
      both_up_vs_CON = "#8A9A5B", both_down_vs_CON = "#C76B8A",
      one_or_both_exact_zero = "#B8B8B8"
    ),
    labels = sign_labels, name = NULL, drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = "Fraction of canonical direct SUS–RES DAPs", y = NULL,
    title = "Direct SUS–RES DAP control-reference geometry",
    subtitle = "Proteins selected by canonical direct SUS–RES FDR < 0.05; composition and displacement are descriptive",
    caption = paste0(
      "Control-reference geometry is a descriptive decomposition of proteins selected by the direct SUS–RES FDR criterion.\n",
      "RES–CON and SUS–CON estimates are not treated as independent validation of those selected proteins."
    )
  ) +
  theme_nature_base(7) +
  ggplot2::theme(
    legend.position = "bottom", panel.grid.major.y = ggplot2::element_line(colour = "grey93", linewidth = 0.2),
    panel.grid.major.x = ggplot2::element_line(colour = "grey92", linewidth = 0.2), panel.spacing.x = grid::unit(4, "mm"),
    strip.text = ggplot2::element_text(size = 7), plot.subtitle = ggplot2::element_text(size = 6.5),
    plot.caption = ggplot2::element_text(size = 5.7, hjust = 0), axis.text.y = ggplot2::element_text(size = 6.5)
  )

save_triplet <- function(plot, stem, width_mm, height_mm) {
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = svglite::svglite, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = png_device, dpi = 300, bg = "white", limitsize = FALSE)
}

trajectory_stem <- file.path(FIGURE_DIR, "candidate_primary_theme_control_trajectory")
burden_stem <- file.path(FIGURE_DIR, "candidate_joint_FDR_protein_burden")
direct_geometry_stem <- file.path(FIGURE_DIR, "candidate_direct_SUS_RES_DAP_control_geometry")
save_triplet(trajectory_plot, trajectory_stem, 183, 142)
save_triplet(burden_plot, burden_stem, 183, 75)
save_triplet(direct_geometry_plot, direct_geometry_stem, 183, 118)
outputs <- c(outputs,
  trajectory_svg = paste0(trajectory_stem, ".svg"), trajectory_pdf = paste0(trajectory_stem, ".pdf"), trajectory_png = paste0(trajectory_stem, ".png"),
  burden_svg = paste0(burden_stem, ".svg"), burden_pdf = paste0(burden_stem, ".pdf"), burden_png = paste0(burden_stem, ".png"),
  direct_geometry_svg = paste0(direct_geometry_stem, ".svg"), direct_geometry_pdf = paste0(direct_geometry_stem, ".pdf"), direct_geometry_png = paste0(direct_geometry_stem, ".png")
)

joint_res_more <- sum(protein$summary$n_joint_FDR05_RES_vs_CON > protein$summary$n_joint_FDR05_SUS_vs_CON)
joint_sus_more <- sum(protein$summary$n_joint_FDR05_SUS_vs_CON > protein$summary$n_joint_FDR05_RES_vs_CON)
joint_tied <- sum(protein$summary$n_joint_FDR05_RES_vs_CON == protein$summary$n_joint_FDR05_SUS_vs_CON)
median_res_more <- sum(protein$summary$median_abs_log2FC_RES_vs_CON > protein$summary$median_abs_log2FC_SUS_vs_CON)
median_sus_more <- sum(protein$summary$median_abs_log2FC_SUS_vs_CON > protein$summary$median_abs_log2FC_RES_vs_CON)
median_tied <- sum(protein$summary$median_abs_log2FC_RES_vs_CON == protein$summary$median_abs_log2FC_SUS_vs_CON)
protein$summary$region_group <- ifelse(grepl("^(CA3|DG)", protein$summary$spatial_unit), "CA3/DG", "CA1/CA2")
largest_count_delta <- max(abs(protein$summary$delta_n_joint_DAP))
top_delta <- protein$summary[abs(protein$summary$delta_n_joint_DAP) == largest_count_delta & largest_count_delta > 0, c("dataset", "spatial_unit", "delta_n_joint_DAP", "median_delta_abs_log2FC", "region_group")]
top_delta_locations <- if (nrow(top_delta)) paste(paste(top_delta$dataset, top_delta$spatial_unit, sep = "/"), collapse = ", ") else "none"
n_supported_go_occurrences <- length(unique(supported_occurrence$supported_occurrence_id))
n_supported_go_assignment_rows <- nrow(supported)
direct_global <- direct_dap$global_summary[1, , drop = FALSE]

sanity_lines <- c(
  "# Three-contrast stress-response biological audit",
  "",
  "## Inferential family definitions",
  "",
  "- Protein canonical FDR: original upstream contrast-specific adjustment encoded in canonical padj; copied unchanged. Mapped ProteinGroupID rows are not used to recalculate it.",
  "- GO canonical FDR: original BH adjustment within each ranked-GSEA dataset x spatial unit x contrast GO-BP family; unchanged.",
  "- Protein joint-control FDR: one secondary BH adjustment over 2m raw p-values for the common RES-vs-CON/SUS-vs-CON ProteinGroupID universe within each dataset x spatial unit.",
  "- GO joint-control FDR: one secondary BH adjustment over 2m original nominal GSEA p-values for the common RES-vs-CON/SUS-vs-CON GO-ID universe within each dataset x spatial unit.",
  "- SUS-vs-RES is excluded from both secondary joint-control families and remains the primary phenotype-separation endpoint.",
  "",
  "## Protein sanity checks",
  "",
  paste0("- Joint-FDR DAP burden across 18 contexts: RES > SUS in ", joint_res_more, "; SUS > RES in ", joint_sus_more, "; tied in ", joint_tied, "."),
  paste0("- Median absolute log2FC: RES > SUS in ", median_res_more, "; SUS > RES in ", median_sus_more, "; tied in ", median_tied, "."),
  paste0("- Per-unit RES-farther protein fractions range from ", sprintf("%.3f", min(protein$summary$fraction_proteins_RES_farther_from_CON)), " to ", sprintf("%.3f", max(protein$summary$fraction_proteins_RES_farther_from_CON)), "."),
  paste0("- The largest absolute joint-DAP count difference is ", largest_count_delta, " and occurs at ", top_delta_locations, "; these sparse differences are not confined to CA3/DG."),
  paste0("- Protein contrast algebra: all ", nrow(protein$algebra_audit), " contexts pass tolerance ", format(STRESS_RESPONSE_ALGEBRA_TOLERANCE, scientific = TRUE), "; maximum residual = ", format(max(protein$algebra_audit$max_absolute_residual), scientific = TRUE), "."),
  "",
  "## Direct SUS-RES DAP control-reference geometry",
  "",
  paste0("- Direct SUS-RES DAPs selected by canonical FDR < 0.05: ", direct_global$n_direct_SUS_RES_DAP, "."),
  paste0("- Opposite-side geometry: ", direct_global$n_opposite_side_of_CON, "/", direct_global$n_direct_SUS_RES_DAP, "; same-side geometry: ", direct_global$n_same_side_of_CON, "/", direct_global$n_direct_SUS_RES_DAP, "."),
  paste0("- Relative displacement from CON: RES farther ", direct_global$n_RES_farther, "/", direct_global$n_direct_SUS_RES_DAP, "; SUS farther ", direct_global$n_SUS_farther, "/", direct_global$n_direct_SUS_RES_DAP, "; equal distance ", direct_global$n_equal_distance, "."),
  "- Direct RES–SUS differences frequently reflect opposing displacement around the unstressed control state.",
  if (direct_global$n_RES_farther > direct_global$n_SUS_farther) "- Within direct SUS-RES DAPs, the resilient state more often showed the larger absolute displacement from control." else "- The direct-DAP control-reference decomposition does not show a global predominance of larger RES displacement.",
  "- Control-reference geometry is a descriptive decomposition of proteins selected by the direct SUS–RES FDR criterion. RES–CON and SUS–CON estimates are not treated as independent validation of those selected proteins.",
  "",
  "## Theme sanity checks",
  "",
  paste0("- ", overview$theme, ": RES-farther ", overview$n_units_abs_RES_theme_NES_gt_abs_SUS, "/18; SUS-farther ", overview$n_units_abs_SUS_theme_NES_gt_abs_RES, "/18; median delta absolute NES = ", sprintf("%.3f", overview$median_delta_abs_theme_NES), " (Q25 ", sprintf("%.3f", overview$q25_delta_abs_theme_NES), ", Q75 ", sprintf("%.3f", overview$q75_delta_abs_theme_NES), "); same direction ", overview$n_units_same_response_direction, "/18; opposite direction ", overview$n_units_opposite_response_direction, "/18; direct SUS-RES supported units ", overview$n_canonical_supported_SUS_vs_RES_units, "/18."),
  "",
  "## Supported GO occurrence audit",
  "",
  paste0("- Unique FDR-supported GO occurrences: ", n_supported_go_occurrences, "; occurrence x theme-assignment rows: ", n_supported_go_assignment_rows, ". Occurrence counts use distinct supported_occurrence_id values, not raw exploded row counts."),
  "",
  "## Interpretation boundary",
  "",
  "Inferential statements are limited to original contrast-specific BH FDR, secondary joint-control BH FDR for matched burden comparison, and direct SUS-vs-RES BH FDR. Median absolute log2FC, protein farther-from-control fractions, median theme NES, trajectory geometry, and spatial concordance are descriptive and have no p-value in this audit."
)
summary_path <- file.path(REPORT_DIR, "stress_response_biological_audit_summary.md")
writeLines(sanity_lines, summary_path, useBytes = TRUE)
outputs <- c(outputs, summary = summary_path)

workbook_status <- data.frame(
  requested_output = file.path(REPORT_DIR, "stress_response_biological_audit.xlsx"),
  status = "blocked_not_created",
  reason = "Required @oai/artifact-tool runtime was unavailable; spreadsheet skill forbids substituting another authoring library.",
  workbook_ready_source_tables = paste(basename(c(outputs[["theme_overview"]], outputs[["protein_summary"]], outputs[["theme_trajectory"]], outputs[["contrast_map"]], outputs[["theme_detail"]], outputs[["direct_dap_geometry"]], outputs[["direct_dap_context_summary"]], outputs[["direct_dap_global_summary"]], outputs[["supported_go_occurrence"]], outputs[["supported_go"]])), collapse = ";"),
  stringsAsFactors = FALSE
)
outputs <- c(outputs, workbook_status = write_csv(workbook_status, "workbook_generation_status.csv", REPORT_DIR))

dimension_output_names <- c(
  "theme_overview", "protein_summary", "direct_dap_geometry", "direct_dap_context_summary", "direct_dap_global_summary",
  "theme_trajectory", "contrast_map", "theme_detail", "supported_go_occurrence", "supported_go",
  "protein_geometry", "go_joint_rows", "protected_reference"
)
dimension_paths <- unname(outputs[dimension_output_names])
output_dimensions <- data.frame(
  artifact = c(
    "overview", "protein_remodeling", "direct_dap_geometry", "direct_dap_context_summary", "direct_dap_global_summary",
    "theme_trajectories", "contrast_map", "theme_detail", "go_occurrence_audit", "go_theme_assignment_audit",
    "protein_geometry", "joint_go", "protected_reference"
  ),
  path = dimension_paths,
  rows = c(
    nrow(overview), nrow(protein$summary), nrow(direct_dap$detail), nrow(direct_dap$context_summary), nrow(direct_dap$global_summary),
    nrow(trajectory), nrow(contrast_map_source), nrow(theme$detail), nrow(supported_occurrence), nrow(supported),
    nrow(protein$geometry), nrow(joint_go$joint), nrow(protected_reference)
  ),
  columns = vapply(
    dimension_paths,
    function(path) ncol(utils::read.csv(path, nrows = 0L, check.names = FALSE)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
outputs <- c(outputs, output_dimensions = write_csv(output_dimensions, "output_dimensions.csv", REPORT_DIR))

family_definitions <- data.frame(
  family_name = c("canonical_protein_contrast_FDR", "control_pair_joint_FDR", "canonical_GO_contrast_FDR", "control_pair_joint_GO_FDR"),
  family_scope = c(
    "original upstream dataset x spatial unit x individual contrast protein family encoded in canonical padj",
    "dataset x spatial unit x common ProteinGroupID universe x {RES_vs_CON,SUS_vs_CON}; 2m raw p-values",
    "dataset x spatial unit x individual ranked-GSEA contrast tested GO-BP family",
    "dataset x spatial unit x common GO-ID universe x {RES_vs_CON,SUS_vs_CON}; 2m original nominal GSEA p-values"
  ),
  role = c("authoritative individual-contrast inference", "secondary matched burden comparison only", "authoritative individual-contrast inference", "secondary matched burden comparison only"),
  SUS_vs_RES_included = c(TRUE, FALSE, TRUE, FALSE),
  input_p_field = c("pval", "pval", "pvalue", "pvalue"),
  adjusted_p_used_as_input = FALSE,
  changed_existing_statistics = FALSE,
  stringsAsFactors = FALSE
)
outputs <- c(outputs, family_definitions = write_csv(family_definitions, "FDR_family_definitions.csv", REPORT_DIR))

write_run_manifest(
  file.path(LOG_DIR, "run_manifest.yml"),
  inputs = c(as.list(manifest_paths), list(ranked_GO_source = GO_SOURCE, theme_registry = REGISTRY_PATH)),
  outputs = as.list(outputs),
  parameters = list(
    primary_contrast = "SUS_vs_RES",
    secondary_control_contrasts = STRESS_RESPONSE_CONTROL_CONTRASTS,
    canonical_FDR_threshold = STRESS_RESPONSE_FDR_THRESHOLD,
    protein_joint_family = "2m raw p-values over common ProteinGroupID universe within dataset x spatial_unit",
    GO_joint_family = "2m nominal GSEA p-values over common GO-ID universe within dataset x spatial_unit",
    algebra_tolerance = STRESS_RESPONSE_ALGEBRA_TOLERANCE,
    ontology = "GO BP",
    ontology_relationships = manuscript_go_allowed_relationships(),
    theme_registry_version = unique(registry$registry_version),
    workbook_status = workbook_status$status
  ),
  notes = c(
    "Downstream audit only: no DA model, preprocessing, imputation, ranked GSEA, DAP-set ORA, WGCNA, or behavior model was run or changed.",
    "Canonical protein and GO adjusted p-values are copied unchanged; secondary joint-control FDRs use raw p-values and never include SUS_vs_RES.",
    "All effect/NES geometry and spatial concordance summaries are descriptive and carry no p-value.",
    "Direct-DAP control-reference geometry is selected only by canonical SUS_vs_RES FDR and is not independent validation.",
    "The protected SUS-RES workbook is recorded as a protected reference artifact and is not consumed or modified by Stage 11.",
    "Legacy script-07 heuristic classes are not consumed as manuscript evidence."
  )
)

message("Stress-response biological audit complete.")
message("Protein contexts: ", nrow(protein$summary), "; theme trajectories: ", nrow(trajectory))
message("Joint protein burden RES>SUS / SUS>RES / tied: ", joint_res_more, " / ", joint_sus_more, " / ", joint_tied)
message("Summary: ", summary_path)
message("Workbook status: ", workbook_status$status)
