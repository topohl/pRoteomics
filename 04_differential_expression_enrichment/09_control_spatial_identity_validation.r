# Figure 2 control-anatomy validation only; does not consume or modify stress/WGCNA outputs.
source("R/paths.R"); source("R/dataset_config.R"); source("R/dataset_inputs.R")
source("R/validation_utils.R"); source("R/qc_exploration_utils.R")
source("R/protein_group_enrichment_utils.R"); source("R/control_spatial_identity_utils.R")
source("R/clusterprofiler_reproducibility.R")
source("R/plotting_nature.R")
suppressPackageStartupMessages({ library(limma); library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2) })

control_spatial_publication_contrast_label <- function(x) {
  labels <- control_spatial_contrast_label(x)
  labels[x == "DG_neuropil_vs_mean_non_DG_regions"] <- "DG neuropil vs non-DG"
  labels
}

control_spatial_prepare_figure2f <- function(
    go, contrasts, panel_contract, interpretation_note,
    analytical_source_path = NULL, require_dg_layers = FALSE) {
  universe <- sort(unique(as.character(go$contrast)))
  missing <- setdiff(contrasts, universe)
  if (length(missing)) {
    stop("Completed GO-BP table is missing display contrast(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  dg_layers <- c("DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers")
  if (require_dg_layers && !all(dg_layers %in% universe)) {
    stop("Completed GO-BP table no longer retains both DG-layer contrasts.", call. = FALSE)
  }
  out <- control_spatial_select_go_display(go, max_terms = 2L, contrasts = contrasts)
  out <- control_spatial_complete_go_display_grid(go, out, contrasts)
  out$internal_contrast_label <- control_spatial_publication_contrast_label(out$contrast)
  out$semantic_simplification_applied <- FALSE
  out$display_selection_rule <-
    paste(
      "completed + positive NES + adjusted p < 0.05;",
      "adjusted p-value, descending NES, GO ID; maximum two per contrast"
    )
  out$panel_contract <- panel_contract
  out$interpretation_note <- interpretation_note(out$contrast)
  if (!is.null(analytical_source_path)) {
    out$analytical_source_file <- normalizePath(
      analytical_source_path, winslash = "/", mustWork = TRUE
    )
    out$analytical_contrast_universe <- paste(universe, collapse = ";")
    out$DG_layer_results_retained_in_analytical_source <-
      all(dg_layers %in% universe)
  }
  out$.contrast_order <- match(out$contrast, contrasts)
  out <- out[order(
    out$.contrast_order, !out$display_term_available,
    out$p_adjust, -out$NES, out$ID, na.last = TRUE
  ), , drop = FALSE]
  out$.contrast_order <- NULL
  out
}

control_spatial_build_figure2f <- function(
    go_display, contrasts, grouped_layout = FALSE, label_wrap = 18L) {
  manuscript_text <- nature_manuscript_text_sizes_pt()
  required <- c(
    "contrast", "internal_contrast_label", "ID", "Description", "NES",
    "p_adjust", "display_evidence_status", "display_term_available"
  )
  missing_columns <- setdiff(required, names(go_display))
  if (length(missing_columns)) {
    stop("Figure 2f display source is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (grouped_layout) {
    layout <- control_spatial_figure2f_grouped_layout()
    idx <- match(go_display$contrast, layout$contrast)
    if (anyNA(idx) || !identical(layout$contrast, contrasts)) {
      stop("Grouped Figure 2f layout does not match the seven-contrast contract.", call. = FALSE)
    }
    go_display$contrast_plot_label <- factor(
      layout$short_label[idx], levels = layout$short_label
    )
    go_display$contrast_group <- factor(
      layout$group[idx], levels = unique(layout$group)
    )
  } else {
    wrapped <- function(z) paste(strwrap(z, width = label_wrap), collapse = "\n")
    go_display$contrast_plot_label <- factor(
      vapply(go_display$internal_contrast_label, wrapped, character(1)),
      levels = vapply(
        control_spatial_publication_contrast_label(contrasts), wrapped, character(1)
      )
    )
  }
  term_rows <- go_display[
    order(
      !go_display$display_term_available, go_display$p_adjust,
      -go_display$NES, go_display$ID, na.last = TRUE
    ),
    c("ID", "Description"), drop = FALSE
  ]
  term_rows <- term_rows[!duplicated(term_rows$ID), , drop = FALSE]
  go_display$term_plot_key <- factor(as.character(go_display$ID), levels = rev(term_rows$ID))
  term_labels <- stats::setNames(
    vapply(term_rows$Description,
           function(z) paste(strwrap(z, width = 38), collapse = "\n"), character(1)),
    term_rows$ID
  )
  go_display$minus_log10_adjusted_p <- -log10(
    pmax(go_display$p_adjust, .Machine$double.xmin)
  )
  evidence <- go_display[go_display$display_term_available %in% TRUE, , drop = FALSE]
  missing_evidence <- go_display[!(go_display$display_term_available %in% TRUE), , drop = FALSE]
  plot <- ggplot(
    go_display,
    aes(contrast_plot_label, term_plot_key)
  ) +
    geom_point(
      data = evidence, aes(fill = NES, size = minus_log10_adjusted_p),
      shape = 21, colour = "#FFFFFF", stroke = 0.24, alpha = 0.96
    ) +
    geom_text(
      data = missing_evidence, label = "x", colour = "#777777",
      size = 2.5, family = "Arial"
    ) +
    scale_y_discrete(labels = term_labels, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_gradientn(colours = nature_palette("support"), name = "NES") +
    scale_size_continuous(
      name = expression(-log[10]("adjusted p-value")), range = c(2.6, 5.4)
    ) +
    guides(
      fill = guide_colourbar(
        order = 1, title.position = "top",
        barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.5, "mm")
      ),
      size = guide_legend(order = 2, title.position = "top", nrow = 1)
    ) +
    labs(x = NULL, y = NULL) +
    theme_nature_manuscript_panel(
      base_size = manuscript_text[["normal"]], base_family = "Arial",
      axes = FALSE, publication_legible = TRUE
    ) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_blank(),
      panel.spacing = grid::unit(0.6, "mm"),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(size = manuscript_text[["dense"]], lineheight = 0.94),
      axis.text.x = element_text(
        size = manuscript_text[["dense"]], lineheight = 0.88, margin = margin(t = 1)
      ),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = element_text(size = manuscript_text[["normal"]]),
      legend.text = element_text(size = manuscript_text[["normal"]]),
      legend.spacing.x = grid::unit(1, "mm"),
      plot.margin = margin(1.2, 1.5, 1.0, 1.2)
    )
  if (grouped_layout) {
    plot <- plot +
      facet_grid(
        cols = vars(contrast_group), scales = "free_x", space = "free_x"
      ) +
      theme(
        panel.spacing.x = grid::unit(1.6, "mm"),
        strip.background = element_rect(fill = "#F7F7F7", colour = NA),
        strip.text.x = element_text(
          face = "plain", size = manuscript_text[["normal"]],
          margin = margin(1, 2, 1, 2)
        )
      )
  }
  list(plot = plot, source = go_display, term_rows = term_rows)
}

control_spatial_save_figure2f <- function(rendered, path, n_contrasts) {
  ggsave(
    path, rendered$plot,
    width = max(154, 54 + 18 * n_contrasts),
    height = max(78, 34 + 4 * nrow(rendered$term_rows)),
    units = "mm", device = svglite::svglite, bg = "white", limitsize = FALSE
  )
}

control_spatial_require_candidate_writable <- function(paths, allow_env, label) {
  existing <- paths[file.exists(paths)]
  allowed <- tolower(trimws(Sys.getenv(allow_env, ""))) %in% c("1", "true", "yes")
  if (length(existing) && !allowed) {
    stop("Refusing to overwrite existing ", label, " output(s): ",
         paste(existing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

control_spatial_render_figures <- function(ks, go, source_data, figures) {
  manuscript_text <- nature_manuscript_text_sizes_pt()
  ks <- ks[ks$internal_contrast %in% control_spatial_figure2e_external_contrasts(), , drop = FALSE]
  ks <- control_spatial_apply_signature_fdr(ks)
  ks$panel_contract <- "Figure 2e: External anatomical validation against matched or explicitly approximate Kaulich regional / CA1-stratum reference signatures."
  ks$matched_exactly <- ks$match_type == "exact"
  ks$external_reference_context <- ks$external_signature == "SP" & !ks$matched_exactly
  utils::write.csv(ks, source_data[["e"]], row.names = FALSE)
  figure2f_contrasts <- control_spatial_figure2f_display_contrasts()
  go_display <- control_spatial_prepare_figure2f(
    go, figure2f_contrasts,
    "Figure 2f: Internal CON-only positive target-enriched GO-BP terms.",
    function(x) ifelse(
      x %in% c("DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers"),
      "Internal DG layer contrast; no matched DG-layer Kaulich reference signature available.",
      "Internal CON-only positive target-enriched GO-BP terms."
    )
  )
  utils::write.csv(go_display, source_data[["f"]], row.names = FALSE)

  kaulich_plot_data <- control_spatial_figure2e_plot_data(ks)
  if (nrow(kaulich_plot_data)) {
    kaulich_plot_data$internal_contrast_label <- control_spatial_publication_contrast_label(
      kaulich_plot_data$internal_contrast
    )
    kaulich_plot_data$external_signature_plot_label <- ifelse(
      kaulich_plot_data$external_reference_context,
      paste0(kaulich_plot_data$external_signature, " (reference only)"),
      as.character(kaulich_plot_data$external_signature)
    )
    kaulich_plot_data$internal_contrast_label <- vapply(
      kaulich_plot_data$internal_contrast_label,
      function(z) paste(strwrap(z, width = 24), collapse = "\n"),
      character(1)
    )
    domain_plot_labels <- c(
      soma_tissue = "Soma region \u2014 tissue reference",
      neuropil_subregion = "Neuropil region \u2014 synaptosome reference",
      ca1_strata = "CA1 strata \u2014 synaptosome reference"
    )
    kaulich_plot_data$validation_domain_plot_label <- factor(
      unname(domain_plot_labels[kaulich_plot_data$signature_fdr_family]),
      levels = unname(domain_plot_labels)
    )
    colour_limit <- max(abs(kaulich_plot_data$NES), na.rm = TRUE)
    significant <- kaulich_plot_data$outline_status == "FDR < 0.05"
    figure2e <- ggplot(kaulich_plot_data, aes(internal_contrast_label, external_signature_plot_label, fill = NES)) +
      geom_point(
        data = kaulich_plot_data[!significant, , drop = FALSE],
        shape = 21, size = 3.2, colour = "#D9D9D9", stroke = 0.35
      ) +
      geom_point(
        data = kaulich_plot_data[significant, , drop = FALSE],
        aes(colour = outline_status),
        shape = 21, size = 3.2, stroke = 1
      ) +
      facet_wrap(~validation_domain_plot_label, ncol = 1, scales = "free") +
      scale_fill_gradient2(
        low = nature_palette("signed")[["low"]],
        mid = nature_palette("signed")[["mid"]],
        high = nature_palette("signed")[["high"]],
        midpoint = 0, limits = c(-colour_limit, colour_limit), name = "NES"
      ) +
      scale_colour_manual(
        values = c("FDR < 0.05" = "black"),
        breaks = "FDR < 0.05",
        labels = c("FDR < 0.05" = "Black outline"),
        name = "Signature FDR < 0.05"
      ) +
      guides(
        fill = guide_colourbar(
          order = 1, title.position = "top",
          barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.5, "mm"),
          theme = theme(legend.text = element_text(size = manuscript_text[["normal"]] / 0.7))
        ),
        colour = guide_legend(
          order = 2, title.position = "top",
          override.aes = list(shape = 21, fill = "white", size = 3, stroke = 1)
        )
      ) +
      labs(x = NULL, y = NULL) +
      theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
      theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.spacing.y = grid::unit(0.8, "mm"),
        strip.background = element_rect(fill = "#F7F7F7", colour = NA),
        strip.text = element_text(face = "plain", hjust = 0, size = manuscript_text[["normal"]], margin = margin(1, 2, 1, 2)),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = manuscript_text[["dense"]], margin = margin(t = 1)),
        axis.text.y = element_text(size = manuscript_text[["normal"]]),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = manuscript_text[["normal"]]),
        legend.text = element_text(size = manuscript_text[["normal"]]),
        legend.spacing.x = grid::unit(1, "mm"),
        plot.margin = margin(1.2, 1.2, 1.0, 1.2)
      )
    ggsave(
      figures[["e"]], figure2e, width = 110, height = 92, units = "mm",
      device = svglite::svglite, bg = "white", limitsize = FALSE
    )
  } else {
    control_spatial_write_empty_state_svg(
      figures[["e"]], "No Kaulich signatures met the prespecified mapping threshold"
    )
  }

  if (nrow(go_display)) {
    control_spatial_save_figure2f(
      control_spatial_build_figure2f(go_display, figure2f_contrasts),
      figures[["f"]], length(figure2f_contrasts)
    )
  } else {
    control_spatial_write_empty_state_svg(
      figures[["f"]], "No GO-BP terms met the prespecified display criteria"
    )
  }
  invisible(list(figure2e_source = ks, figure2f_source = go_display))
}

control_spatial_render_figure2f_regions_ca1layers_candidate <- function(
    go, source_path, figure_path, analytical_source_path) {
  contrasts <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  go_display <- control_spatial_prepare_figure2f(
    go, contrasts,
    paste(
      "Figure 2f candidate: Internal CON-only soma-region and CA1-laminar",
      "positive target-enriched GO-BP terms; presentation-only selection from completed results."
    ),
    function(x) rep(
      "Internal CON-only positive target-enriched GO-BP terms.", length(x)
    ),
    analytical_source_path = analytical_source_path,
    require_dg_layers = TRUE
  )
  utils::write.csv(go_display, source_path, row.names = FALSE)
  control_spatial_save_figure2f(
    control_spatial_build_figure2f(go_display, contrasts, label_wrap = 13L),
    figure_path, length(contrasts)
  )
  invisible(go_display)
}

control_spatial_render_figure2f_grouped_layout_candidate <- function(
    source_path, figure_path) {
  go_display <- utils::read.csv(source_path, check.names = FALSE)
  contrasts <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  counts <- table(factor(go_display$contrast, levels = contrasts))
  if (!setequal(unique(as.character(go_display$contrast)), contrasts) ||
      any(counts < 1L | counts > 2L) ||
      !all(c("display_evidence_status", "display_term_available") %in% names(go_display))) {
    stop(
      "Grouped Figure 2f input must retain every contrast with zero to two supported terms.",
      call. = FALSE
    )
  }
  control_spatial_save_figure2f(
    control_spatial_build_figure2f(
      go_display, contrasts, grouped_layout = TRUE
    ),
    figure_path, length(contrasts)
  )
  invisible(go_display)
}
control_spatial_identity_main <- function() {
message("Starting control spatial identity validation")
root <- repo_path(); wb <- repo_path("data", "external", "kaulich_2025", "kaulich_supplementary_data_2.xlsx")
message("Resolved project root: ", root)
if (!file.exists(wb)) stop("Required Kaulich workbook not found: ", wb, call. = FALSE)
expected_sheets <- c("subregion tissue FCs", "subregion tissue GO", "strata tissue FCs", "strata tissue GO", "subregion syn FCs", "strata syn FCs")
missing_sheets <- setdiff(expected_sheets, readxl::excel_sheets(wb))
if (length(missing_sheets)) stop("Kaulich workbook is missing expected sheets: ", paste(missing_sheets, collapse = ", "), call. = FALSE)
wb_hash <- file_hash_sha256(wb)
message("Validated Kaulich workbook: ", wb, " (SHA-256 ", wb_hash, ")")
if (is_dry_run()) {
  message("[dry-run] Inputs validated; no scientific outputs will be written.")
  return(invisible(list(status = "dry_run", project_root = root, workbook = wb)))
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required for GSEA reproducibility settings.", call. = FALSE)
}
gsea_cfg <- validate_clusterprofiler_gsea_config(yaml::read_yaml(
  repo_path("config", "clusterProfiler_config.yml")
))
gsea_seed_base <- gsea_cfg$analysis$gsea_seed_base
n_perm_simple <- gsea_cfg$analysis$n_perm_simple
validation_output_root <- Sys.getenv(
  "PROTEOMICS_CONTROL_SPATIAL_OUTPUT_ROOT", unset = ""
)
if (!nzchar(validation_output_root)) validation_output_root <- root
validation_output_root <- normalizePath(
  validation_output_root, winslash = "/", mustWork = FALSE
)
out <- function(kind, name) file.path(validation_output_root, "results", kind, "04_differential_expression_enrichment", "control_spatial_identity_validation", "global", name)
tables <- c(protein=out("tables","anatomical_protein_contrasts.csv"), mapping=out("tables","kaulich_signature_mapping.csv"), kaulich=out("tables","kaulich_signature_gsea.csv"), go=out("tables","control_anatomical_go_bp_gsea.csv"), matching_audit=out("tables","figure2e_matching_audit.csv"), status=out("tables","analysis_status.csv"))
source_data <- c(e=out("source_data","figure2e_source_data.csv"), f=out("source_data","figure2f_source_data.csv")); figures <- c(e=out("figures","figure2e_kaulich_validation.svg"), f=out("figures","figure2f_control_anatomical_GO.svg")); manifest <- out("logs","run_manifest.yml")
candidate_source <- out("source_data", "figure2f_regions_CA1layers_source_data.csv")
candidate_figure <- out("figures", "figure2f_control_anatomical_GO_regions_CA1layers.svg")
grouped_candidate_figure <- out(
  "figures", "figure2f_control_anatomical_GO_regions_CA1layers_grouped.svg"
)
dir.create(dirname(tables[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(source_data[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(figures[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(manifest), recursive=TRUE, showWarnings=FALSE)

render_only <- tolower(trimws(Sys.getenv("PROTEOMICS_CONTROL_SPATIAL_RENDER_ONLY", ""))) %in% c("1", "true", "yes")
candidate_render_only <- tolower(trimws(Sys.getenv(
  "PROTEOMICS_CONTROL_SPATIAL_FIGURE2F_CA1LAYERS_RENDER_ONLY", ""
))) %in% c("1", "true", "yes")
grouped_candidate_render_only <- tolower(trimws(Sys.getenv(
  "PROTEOMICS_CONTROL_SPATIAL_FIGURE2F_GROUPED_LAYOUT_RENDER_ONLY", ""
))) %in% c("1", "true", "yes")
if (grouped_candidate_render_only) {
  message("Starting presentation-only grouped Figure 2f candidate rendering")
  control_spatial_validate_output_bundle(candidate_source)
  control_spatial_require_candidate_writable(
    grouped_candidate_figure,
    "PROTEOMICS_CONTROL_SPATIAL_FIGURE2F_GROUPED_LAYOUT_ALLOW_OVERWRITE",
    "grouped Figure 2f candidate"
  )
  control_spatial_render_figure2f_grouped_layout_candidate(
    candidate_source, grouped_candidate_figure
  )
  control_spatial_validate_output_bundle(grouped_candidate_figure)
  message("Grouped Figure 2f candidate rendering complete")
  return(invisible(list(
    source_data = candidate_source, figure = grouped_candidate_figure,
    status = "presentation_only_grouped_layout"
  )))
}
if (candidate_render_only) {
  message("Starting presentation-only Figure 2f regions/CA1-layers candidate rendering")
  control_spatial_validate_output_bundle(tables[["go"]])
  candidate_paths <- c(candidate_source, candidate_figure)
  control_spatial_require_candidate_writable(
    candidate_paths,
    "PROTEOMICS_CONTROL_SPATIAL_FIGURE2F_CA1LAYERS_ALLOW_OVERWRITE",
    "Figure 2f candidate"
  )
  go <- utils::read.csv(tables[["go"]], check.names = FALSE)
  control_spatial_render_figure2f_regions_ca1layers_candidate(
    go, candidate_source, candidate_figure, tables[["go"]]
  )
  control_spatial_validate_output_bundle(candidate_paths)
  message("Figure 2f regions/CA1-layers candidate rendering complete")
  return(invisible(list(
    analytical_source = tables[["go"]], source_data = candidate_source,
    figure = candidate_figure, status = "presentation_only_candidate"
  )))
}
if (render_only) {
  message("Starting presentation-only control spatial identity rendering")
  control_spatial_validate_output_bundle(c(tables, manifest))
  ks <- utils::read.csv(tables[["kaulich"]], check.names = FALSE)
  go <- utils::read.csv(tables[["go"]], check.names = FALSE)
  ks <- ks[ks$internal_contrast %in% control_spatial_figure2e_external_contrasts(), , drop = FALSE]
  ks <- control_spatial_apply_signature_fdr(ks)
  utils::write.csv(ks, tables[["kaulich"]], row.names = FALSE)
  utils::write.csv(control_spatial_figure2e_matching_audit(ks), tables[["matching_audit"]], row.names = FALSE)
  control_spatial_render_figures(ks, go, source_data, figures)
  control_spatial_validate_output_bundle(c(tables, source_data, figures, manifest))
  message("Control spatial identity presentation rendering complete")
  return(invisible(list(
    tables = tables, source_data = source_data, figures = figures,
    manifest = manifest, status = "presentation_only"
  )))
}

regional_signatures <- list(
  tissue = control_spatial_kaulich_region_signatures(wb, "subregion tissue FCs", wb_hash),
  synaptosome = control_spatial_kaulich_region_signatures(wb, "subregion syn FCs", wb_hash)
)
stratum_signatures <- control_spatial_kaulich_signatures(wb, "strata syn FCs", wb_hash)
all_protein <- list(); all_map <- list(); all_ks <- list(); all_go <- list()
statuses <- list(); details <- list(); kaulich_jobs <- list()
go_announced <- FALSE
for (dataset in c("neuron_soma", "neuron_neuropil")) {
  message("Loading ", dataset)
  inputs <- resolve_dataset_inputs(dataset, purpose="wgcna", script="09_control_spatial_identity_validation.r", stage="differential_expression_enrichment")
  canonical_metadata <- path_processed("01_preprocessing", "06_merged_metadata_module_score", dataset, "sample_metadata_merged_clean_for_module_scores.xlsx")
  if (!file.exists(canonical_metadata)) stop("Canonical merged metadata not found: ", canonical_metadata, call. = FALSE)
  # The canonical loader keeps quantitative processing untouched while reconstructing the existing ProteinGroupID/mapping contract.
  canonical <- qc_load_canonical_expression(inputs$expression_file, canonical_metadata, dataset=dataset, strict=TRUE)
  meta <- control_spatial_prepare_metadata(canonical$meta); keep <- colnames(canonical$mat) %in% meta$Sample
  meta <- meta[match(colnames(canonical$mat)[keep], meta$Sample),,drop=FALSE]; mat <- canonical$mat[,keep,drop=FALSE]
  if (anyNA(match(colnames(mat), meta$Sample))) stop("Canonical matrix and metadata failed to align.", call.=FALSE)
  d <- control_spatial_design(meta); design <- d$design
  message("Fitting ", dataset, " model")
  corfit <- limma::duplicateCorrelation(mat, design, block=meta$AnimalID)
  fit <- limma::lmFit(mat, design, block=meta$AnimalID, correlation=corfit$consensus)
  unit_levels <- sub("^anatomical_unit_", "", colnames(design)[grepl("^anatomical_unit_", colnames(design))])
  make_contrast <- function(name, weights) { v <- stats::setNames(rep(0,ncol(design)),colnames(design)); v[paste0("anatomical_unit_", names(weights))] <- weights; list(name=name, weights=v) }
  contrasts <- list()
  if (dataset == "neuron_soma") for (target in sort(unique(toupper(meta$Region)))) contrasts[[length(contrasts)+1L]] <- make_contrast(paste0(target,"_vs_mean_other_soma_regions"), control_spatial_target_rest_weights(unit_levels,target))
  if (dataset == "neuron_neuropil") {
    regions <- sub("_.*$", "", unit_levels); contrasts[[1]] <- make_contrast("DG_neuropil_vs_mean_non_DG_regions", control_spatial_region_mean_weights(unit_levels,regions,"DG"))
    if (all(c("CA1_SO","CA3_SO") %in% unit_levels)) contrasts[[length(contrasts)+1L]] <- make_contrast("CA1_SO_vs_CA3_SO", stats::setNames(c(1,-1),c("CA1_SO","CA3_SO"))) else statuses[[length(statuses)+1L]] <- control_spatial_empty_status(dataset,"CA1_SO_vs_CA3_SO","skipped","required spatial units absent")
    ca1 <- unit_levels[grepl("^CA1_",unit_levels)]; if(length(ca1)>=3L) for (target in ca1) contrasts[[length(contrasts)+1L]] <- make_contrast(paste0(target,"_vs_mean_other_CA1_strata"),control_spatial_target_rest_weights(ca1,target))
    dg <- unit_levels[grepl("^DG_", unit_levels)]
    if (length(dg) >= 2L) for (target in dg) contrasts[[length(contrasts) + 1L]] <- make_contrast(paste0(target, "_vs_mean_other_DG_layers"), control_spatial_target_rest_weights(dg, target))
  }
  for (ct in contrasts) {
    cf <- limma::contrasts.fit(fit, matrix(ct$weights,ncol=1,dimnames=list(names(ct$weights),ct$name))) |> limma::eBayes(robust=TRUE,trend=TRUE)
    tt <- limma::topTable(cf, number=Inf, sort.by="none"); tt$ProteinGroupID <- rownames(tt); tt <- cbind(tt, canonical$feature_table[match(tt$ProteinGroupID,canonical$feature_table$ProteinGroupID), setdiff(names(canonical$feature_table),"ProteinGroupID"),drop=FALSE]); tt$dataset <- dataset; tt$contrast <- ct$name
    all_protein[[length(all_protein)+1L]] <- tt
    gi <- build_enrichment_gene_inputs(tt, strict=TRUE); ranked <- gi$ranked
    if (ct$name %in% control_spatial_figure2e_external_contrasts()) {
      signatures <- if(dataset=="neuron_soma") regional_signatures$tissue else if(grepl("CA1_.*strata",ct$name)) stratum_signatures else regional_signatures$synaptosome
      for (sg in split(signatures, signatures$external_signature)) {
      mapped_signature <- control_spatial_map_signature(sg, names(ranked), mapped_gene_threshold = 5L)
      interpretation <- control_spatial_signature_interpretation(dataset, ct$name, sg$external_signature[1])
      fdr_family <- control_spatial_signature_family(interpretation$validation_domain)
      identity <- data.frame(
        dataset = dataset,
        internal_contrast = ct$name,
        external_signature = sg$external_signature[1],
        validation_domain = interpretation$validation_domain,
        external_source_compartment = interpretation$external_source_compartment,
        expected_match = interpretation$expected_match,
        match_type = interpretation$match_type,
        interpretation_note = interpretation$interpretation_note,
        signature_fdr_family = fdr_family,
        source_compartment = interpretation$external_source_compartment,
        anatomical_match_type = interpretation$match_type,
        stringsAsFactors = FALSE
      )
      mapping_summary <- cbind(
        identity,
        mapped_signature$summary
      )
      candidate_data <- mapped_signature$candidates[
        , setdiff(names(mapped_signature$candidates), "external_signature"), drop = FALSE
      ]
      mapping_candidates <- cbind(
        identity,
        candidate_data,
        mapped_signature$summary[rep(1L, nrow(mapped_signature$candidates)), , drop = FALSE]
      )
      all_map[[length(all_map) + 1L]] <- mapping_candidates
      gsea_prefix <- cbind(
        mapping_summary,
        data.frame(
          source_sheet = sg$source_sheet[1],
          workbook_sha256 = wb_hash,
          stringsAsFactors = FALSE
        )
      )
      kaulich_jobs[[length(kaulich_jobs) + 1L]] <- list(
        prefix = gsea_prefix,
        ranked = ranked,
        mapped = mapped_signature$mapped_official_symbols,
        external_signature = sg$external_signature[1]
      )
      }
    }
    if (!go_announced) { message("Running GO-BP GSEA"); go_announced <- TRUE }
    go_repro <- control_spatial_gsea_reproducibility(
      dataset, ct$name, "gseGO_BP", gsea_seed_base, n_perm_simple
    )
    go <- as.data.frame(run_seeded_clusterprofiler_gsea(
      clusterProfiler::gseGO,
      gsea_seed = go_repro$gsea_seed,
      n_perm_simple = go_repro$n_perm_simple,
      geneList = ranked, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
      ont = "BP", pvalueCutoff = 1, verbose = FALSE
    ))
    if(nrow(go)) {
      go$dataset<-dataset;go$contrast<-ct$name;go$status<-"completed"
      go$p_value <- go$pvalue; go$p_adjust <- go$p.adjust
      if (!"setSize" %in% names(go)) stop("Successful GO-BP GSEA rows are missing setSize.", call. = FALSE)
      go$mapped_unique_genes <- go$setSize
      go <- cbind(go, go_repro[rep(1L, nrow(go)), , drop = FALSE])
      all_go[[length(all_go)+1L]] <- go
    } else all_go[[length(all_go)+1L]]<-cbind(data.frame(dataset=dataset,contrast=ct$name,status="completed_zero_terms"), go_repro)
    statuses[[length(statuses)+1L]] <- control_spatial_empty_status(dataset,ct$name,"completed",paste0("duplicateCorrelation=",round(corfit$consensus,4),"; hemisphere=",if(d$hemisphere_included)"included" else d$hemisphere_omission_reason)); details[[length(details)+1L]]<-list(dataset=dataset,formula=if(d$hemisphere_included)"abundance ~ 0 + anatomical_unit + hemisphere" else "abundance ~ 0 + anatomical_unit",n_samples=ncol(mat),n_animals=length(unique(meta$AnimalID)),contrast=ct$name,weights=ct$weights,duplicate_correlation=corfit$consensus)
  }
}
largest_mapped_signature_size <- max(vapply(kaulich_jobs, function(z) length(z$mapped), integer(1)))
min_gs_size <- 5L
max_gs_size <- control_spatial_signature_max_gs_size(largest_mapped_signature_size)
message("Running Kaulich signature GSEA")
for (job in kaulich_jobs) {
  job_repro <- control_spatial_gsea_reproducibility(
    job$prefix$dataset[[1]], job$prefix$internal_contrast[[1]],
    "GSEA_Kaulich_signature", gsea_seed_base, n_perm_simple,
    external_signature = job$external_signature
  )
  job_prefix <- cbind(job$prefix, job_repro)
  if (length(job$mapped) < min_gs_size) {
    all_ks[[length(all_ks) + 1L]] <- cbind(
      job_prefix,
      data.frame(
        status = "skipped_lt5_mapped_genes",
        NES = NA_real_, p_value = NA_real_, p_adjust = NA_real_,
        leading_edge_genes = NA_character_, stringsAsFactors = FALSE
      )
    )
    next
  }
  z <- run_seeded_clusterprofiler_gsea(
    clusterProfiler::GSEA,
    gsea_seed = job_repro$gsea_seed,
    n_perm_simple = job_repro$n_perm_simple,
    geneList = job$ranked,
    TERM2GENE = data.frame(term = job$external_signature, gene = job$mapped),
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    pvalueCutoff = 1,
    verbose = FALSE
  )
  q <- as.data.frame(z)
  all_ks[[length(all_ks) + 1L]] <- if (nrow(q)) {
    cbind(
      job_prefix,
      data.frame(
        status = "completed",
        NES = q$NES,
        p_value = q$pvalue,
        p_adjust = q$p.adjust,
        leading_edge_genes = q$core_enrichment,
        stringsAsFactors = FALSE
      )
    )
  } else {
    cbind(
      job_prefix,
      data.frame(
        status = "completed_zero_terms",
        NES = NA_real_, p_value = NA_real_, p_adjust = NA_real_,
        leading_edge_genes = NA_character_, stringsAsFactors = FALSE
      )
    )
  }
}
protein<-control_spatial_bind_rows_fill(all_protein)
mapping<-control_spatial_bind_rows_fill(all_map)
ks<-control_spatial_apply_signature_fdr(control_spatial_bind_rows_fill(all_ks))
go<-control_spatial_bind_rows_fill(all_go)
status<-control_spatial_bind_rows_fill(statuses)
message("Writing output bundle")
utils::write.csv(protein,tables[["protein"]],row.names=FALSE); utils::write.csv(mapping,tables[["mapping"]],row.names=FALSE); utils::write.csv(ks,tables[["kaulich"]],row.names=FALSE); utils::write.csv(go,tables[["go"]],row.names=FALSE); utils::write.csv(control_spatial_figure2e_matching_audit(ks),tables[["matching_audit"]],row.names=FALSE); utils::write.csv(status,tables[["status"]],row.names=FALSE)
control_spatial_render_figures(ks, go, source_data, figures)
write_run_manifest(
  manifest,
  inputs = c(workbook = wb),
  outputs = c(tables, source_data, figures),
  parameters = list(
    git_commit = system("git rev-parse HEAD", intern = TRUE),
    external_workbook_sha256 = wb_hash,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    largest_mapped_prespecified_signature_size = largest_mapped_signature_size,
    gsea_seed_base = gsea_seed_base,
    n_perm_simple = n_perm_simple,
    gsea_rng_kind = "L'Ecuyer-CMRG/Inversion/Rejection",
    gsea_seed_identity = "script + dataset + contrast + enrichment type + optional external signature",
    signature_fdr_family_sizes = as.list(control_spatial_signature_family_sizes()),
    figure2f_semantic_simplification = "No standalone repository simplification helper was readily reusable; deterministic top two by adjusted p-value, descending NES, and GO ID.",
    models = details,
    session = utils::sessionInfo()
  ),
  notes = "Supportive external concordance only; three CON animals with repeated spatial observations. Canonical input hashes and fingerprints are retained by qc_load_canonical_expression input manifests."
)
control_spatial_validate_output_bundle(c(tables, source_data, figures, manifest))
message("Control spatial identity validation complete")
invisible(list(tables = tables, source_data = source_data, figures = figures, manifest = manifest, models = details))
}

if (control_spatial_direct_execution(sys.nframe())) {
  options(error = function() {
    traceback(2)
    quit(save = "no", status = 1, runLast = FALSE)
  })
  control_spatial_identity_main()
}
