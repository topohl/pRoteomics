# =====================================================================
# Candidate figure panel renderers.
#
# EVERY function here is a DOWNSTREAM RENDERER. It reads validated canonical
# outputs and joins, filters, reshapes, annotates and plots them. None of them
# fits a model, runs an enrichment, or computes an adjusted p-value. Summary
# statistics that already exist upstream are READ, not recomputed; where a
# renderer needs a count it counts rows of an already-adjudicated status
# column. cf_assert_no_model_fitting() scans this file and a test pins it.
#
# Each renderer has the signature (panel, svg_path, source_csv_path) and must
# write both the SVG and its source-data CSV.
# =====================================================================

# ==================================================== FIGURE 2 CANDIDATES

# F2-C1. Paired hemisphere reproducibility of anatomical contrast effects.
# The correlation/sign-agreement statistics are READ from the validated
# summary table, never recomputed here.
cfp_render_bilateral_concordance <- function(panel, svg_path, source_csv) {
  prot <- cf_read_csv(repo_path(panel$primary_source))
  summ <- cf_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  cf_require_columns(prot, as.character(unlist(panel$required_columns)), panel$id)

  # One representative anatomical contrast per dataset, chosen by a fixed rule:
  # the contrast with the MEDIAN reproducibility for that dataset, so the panel
  # is neither a best-case nor a worst-case showcase.
  pick <- do.call(rbind, lapply(split(summ, summ$dataset), function(d) {
    d <- d[order(d$pearson_r), , drop = FALSE]
    d[ceiling(nrow(d) / 2), , drop = FALSE]
  }))
  pick$facet <- sprintf("%s\n%s", pick$dataset, pick$contrast)

  dat <- merge(prot, pick[, c("dataset", "contrast", "facet", "pearson_r",
                              "spearman_rho", "sign_agreement_fraction",
                              "n_evaluable_proteins")],
               by = c("dataset", "contrast"))
  dat <- dat[is.finite(dat$estimate_L) & is.finite(dat$estimate_R), , drop = FALSE]

  lab <- unique(dat[, c("facet", "pearson_r", "spearman_rho",
                        "sign_agreement_fraction", "n_evaluable_proteins")])
  lab$text <- sprintf("r = %.2f\nrho = %.2f\nsign agree = %.0f%%\nn = %s",
                      lab$pearson_r, lab$spearman_rho,
                      100 * lab$sign_agreement_fraction,
                      format(lab$n_evaluable_proteins, big.mark = ","))
  rng <- range(c(dat$estimate_L, dat$estimate_R), na.rm = TRUE)
  lab$x <- rng[1]; lab$y <- rng[2]

  p <- ggplot2::ggplot(dat, ggplot2::aes(estimate_L, estimate_R)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                         linewidth = 0.3, colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(colour = sign_agreement), size = 0.25,
                        alpha = 0.35) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#23384D", "FALSE" = "#D7303F"),
                                 name = "sign agrees") +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = x, y = y, label = text),
                       hjust = 0, vjust = 1, size = 1.9, lineheight = 1.05,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ facet, nrow = 1) +
    ggplot2::labs(
      title = "Paired hemisphere reproducibility of anatomical contrast effects",
      subtitle = paste0("left versus right estimate of the SAME anatomical contrast; ",
                        "paired measurements within an animal, not independent replicates"),
      x = "left-hemisphere effect", y = "right-hemisphere effect") +
    cf_theme()

  write_csv_safe(dat[, c("dataset", "contrast", "ProteinGroupID", "estimate_L",
                         "estimate_R", "sign_agreement", "pearson_r",
                         "spearman_rho", "sign_agreement_fraction",
                         "n_evaluable_proteins")], source_csv)
  cf_save_panel(p, svg_path, 183, 62)
  invisible(list(status = "ok"))
}

# F2-C2. Precision benefit of bilateral sampling: single-side ICC -> bilateral ICC.
cfp_render_precision_gain <- function(panel, svg_path, source_csv) {
  d <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(d, as.character(unlist(panel$required_columns)), panel$id)
  d <- d[is.finite(d$ICC_single_side) & is.finite(d$ICC_bilateral_mean), , drop = FALSE]
  d$endpoint_class <- factor(d$endpoint_class,
    levels = c("wgcna_module_eigengene", "reference_marker_score",
               "empirical_compartment_score"),
    labels = c("WGCNA modules", "reference marker scores", "empirical compartment scores"))
  long <- rbind(
    data.frame(d[, c("dataset", "endpoint_class", "endpoint_id")],
               sampling = "single side", ICC = d$ICC_single_side),
    data.frame(d[, c("dataset", "endpoint_class", "endpoint_id")],
               sampling = "bilateral mean", ICC = d$ICC_bilateral_mean))
  long$sampling <- factor(long$sampling, levels = c("single side", "bilateral mean"))

  p <- ggplot2::ggplot(long, ggplot2::aes(sampling, ICC)) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(dataset, endpoint_id)),
                       colour = "grey72", linewidth = 0.22) +
    ggplot2::geom_point(ggplot2::aes(colour = sampling), size = 0.5, alpha = 0.75) +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.42,
                          linewidth = 0.28, colour = "#23384D") +
    ggplot2::scale_colour_manual(values = c("single side" = "#9E9A92",
                                            "bilateral mean" = "#23384D"),
                                 guide = "none") +
    ggplot2::facet_grid(dataset ~ endpoint_class) +
    ggplot2::labs(
      title = "Measurement reliability gained by bilateral sampling",
      subtitle = paste0("each line is one endpoint; crossbar is the median. ICC is a ",
                        "reliability descriptor of the measurement, not a stress effect"),
      x = NULL, y = "intraclass correlation") +
    cf_theme()

  out <- d[, c("dataset", "endpoint_class", "endpoint_id", "n_animals",
               "ICC_single_side", "ICC_bilateral_mean",
               "relative_measurement_variance_reduction", "expected_reliability_gain",
               "assumption_status", "is_singular")]
  write_csv_safe(out, source_csv)
  cf_save_panel(p, svg_path, 89, 78)
  invisible(list(status = "ok"))
}

# F2-C3. WGCNA module bilateral reproducibility, two axes only.
cfp_render_wgcna_bilateral <- function(panel, svg_path, source_csv) {
  d <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(d, as.character(unlist(panel$required_columns)), panel$id)
  d <- d[d$level == "module", , drop = FALSE]
  d$reproducible <- d$bilateral_reproducibility_class == "reproducible_level_and_pattern"

  p <- ggplot2::ggplot(d, ggplot2::aes(pearson_r, median_profile_Pearson)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 3, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = reproducible), size = 1.1) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#23384D", "FALSE" = "#D7303F"),
      labels = c("TRUE" = "reproducible level and pattern",
                 "FALSE" = "poor bilateral reproducibility"), name = NULL) +
    ggplot2::facet_wrap(~ dataset, nrow = 1) +
    ggplot2::labs(
      title = "Molecular programs keep their level and their spatial pattern across hemispheres",
      subtitle = "x: absolute left/right module reproducibility.  y: reproducibility of the spatial profile",
      x = "absolute L/R reproducibility (Pearson r)",
      y = "spatial-profile reproducibility") +
    ggplot2::coord_cartesian(xlim = c(-0.2, 1), ylim = c(-0.2, 1)) +
    cf_theme()

  write_csv_safe(d[, c("dataset", "level", "endpoint_id", "n_pairs", "pearson_r",
                       "spearman_rho", "median_profile_Pearson",
                       "sign_agreement_fraction", "bilateral_reproducibility_class")],
                 source_csv)
  cf_save_panel(p, svg_path, 89, 62)
  invisible(list(status = "ok"))
}

# F2-C4. Compartment identity, with cross-hemisphere concordance alongside.
cfp_render_compartment_plus_bilateral <- function(panel, svg_path, source_csv) {
  d <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(d, as.character(unlist(panel$required_columns)), panel$id)
  bil_p <- repo_path(as.character(unlist(panel$input_dependencies))[1])
  bil <- cf_read_csv(bil_p, required = FALSE)

  d$marker_class <- as.character(d$marker_class)
  d$intended_dataset <- as.character(d$intended_dataset)
  agg <- stats::aggregate(intended_minus_comparator_log2 ~ marker_class + intended_dataset,
                          data = d, FUN = stats::median)
  names(agg)[3] <- "median_intended_minus_comparator_log2"

  p1 <- ggplot2::ggplot(d, ggplot2::aes(stats::reorder(marker_class,
                                                       intended_minus_comparator_log2),
                                        intended_minus_comparator_log2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = intended_dataset), outlier.size = 0.2,
                          linewidth = 0.2, alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
    ggplot2::labs(
      title = "Compartment identity: intended minus comparator marker abundance",
      subtitle = "positive means the marker is higher in the compartment it should mark",
      x = NULL, y = "intended - comparator (log2)") +
    cf_theme()

  if (!is.null(bil) && all(c("dataset", "contrast", "pearson_r") %in% names(bil))) {
    b <- bil[, c("dataset", "contrast", "pearson_r", "sign_agreement_fraction")]
    p2 <- ggplot2::ggplot(b, ggplot2::aes(stats::reorder(contrast, pearson_r), pearson_r)) +
      ggplot2::geom_col(ggplot2::aes(fill = dataset), width = 0.68) +
      ggplot2::coord_flip(ylim = c(0, 1)) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
      ggplot2::labs(subtitle = "cross-hemisphere concordance of the same compartment contrast",
                    x = NULL, y = "paired hemisphere Pearson r") +
      cf_theme()
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1.35, 1))
    } else {
      p <- p1
    }
    write_csv_safe(merge(agg, b, by.x = "intended_dataset", by.y = "dataset",
                         all.x = TRUE), source_csv)
  } else {
    p <- p1
    write_csv_safe(agg, source_csv)
  }
  cf_save_panel(p, svg_path, 89, 92)
  invisible(list(status = "ok"))
}

# ==================================================== FIGURE 3 CANDIDATES

# Shared loader for the canonical ranked-GSEA theme assignments. Reads only.
cfp_gsea_theme_table <- function(panel) {
  d <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(d, as.character(unlist(panel$required_columns)), panel$id)
  d <- d[d$assignment_status %in% c("single_theme", "multi_theme", "qc_review"), , drop = FALSE]
  d <- d[nzchar(as.character(d$theme_id)), , drop = FALSE]
  d
}

# Collapse constituent GO terms to a theme cell. This is a DESCRIPTIVE summary
# of canonical values: the median NES of the theme's constituent terms, plus a
# count of how many of those terms are FDR-supported IN THE CANONICAL GSEA.
# No new statistic is produced and no p-value is adjusted.
cfp_theme_cells <- function(d, fdr_cut = 0.05) {
  key <- paste(d$dataset, d$spatial_unit, d$contrast, d$theme_id, sep = "\r")
  sp <- split(seq_len(nrow(d)), key)
  out <- do.call(rbind, lapply(sp, function(ix) {
    z <- d[ix, , drop = FALSE]
    data.frame(
      dataset = z$dataset[1], spatial_unit = z$spatial_unit[1],
      contrast = z$contrast[1], theme_id = z$theme_id[1],
      manuscript_theme = z$manuscript_theme[1], theme_role = z$theme_role[1],
      n_terms = nrow(z),
      n_terms_FDR_supported = sum(is.finite(z$GSEA_FDR) & z$GSEA_FDR < fdr_cut),
      median_NES = stats::median(z$NES, na.rm = TRUE),
      min_constituent_FDR = suppressWarnings(min(z$GSEA_FDR, na.rm = TRUE)),
      stringsAsFactors = FALSE)
  }))
  out$min_constituent_FDR[!is.finite(out$min_constituent_FDR)] <- NA_real_
  out$has_FDR_support <- out$n_terms_FDR_supported > 0L
  out$summary_basis <- paste0(
    "median NES of the theme's constituent canonical GO-BP GSEA terms; ",
    "FDR support counted from those same terms at FDR < ", fdr_cut,
    "; the theme summary is NOT itself a test")
  rownames(out) <- NULL
  out
}

cfp_unit_order <- function(u) {
  lev <- c("CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
           "CA3_so", "CA3_sr", "DG_mo", "DG_po",
           "CA1_sp", "CA2_sp", "CA3_sp", "DG_sg", "CA1", "CA2", "CA3", "DG")
  factor(u, levels = intersect(lev, unique(u)))
}

cfp_gsea_atlas_plot <- function(cells, contrasts, title, subtitle) {
  z <- cells[cells$contrast %in% contrasts, , drop = FALSE]
  z$spatial_unit <- cfp_unit_order(z$spatial_unit)
  # qc_review themes are QC context, never claims. Mark them in the axis label
  # so a reader cannot mistake one for a primary program.
  z$theme_label <- ifelse(z$theme_role == "qc_review",
                          paste0("[QC context] ", z$manuscript_theme),
                          as.character(z$manuscript_theme))
  all_labels <- ifelse(cells$theme_role == "qc_review",
                       paste0("[QC context] ", cells$manuscript_theme),
                       as.character(cells$manuscript_theme))
  z$theme <- factor(z$theme_label,
                    levels = rev(unique(all_labels[
                      order(cells$theme_role, cells$manuscript_theme)])))
  z$contrast <- factor(z$contrast, levels = contrasts)
  lim <- max(abs(z$median_NES), na.rm = TRUE)
  ggplot2::ggplot(z, ggplot2::aes(spatial_unit, theme, fill = median_NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_point(data = z[z$has_FDR_support, , drop = FALSE],
                        ggplot2::aes(size = n_terms_FDR_supported),
                        shape = 21, fill = "black", colour = "white",
                        stroke = 0.15, show.legend = TRUE) +
    ggplot2::scale_size_continuous(range = c(0.35, 1.7),
                                   name = "constituent GO terms\nwith canonical FDR < 0.05") +
    ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                                  midpoint = 0, limits = c(-lim, lim),
                                  name = "median NES") +
    ggplot2::facet_grid(. ~ contrast + dataset, scales = "free_x", space = "free_x") +
    ggplot2::labs(title = if (is.null(title)) NULL else cf_wrap(title, 110),
                  subtitle = if (is.null(subtitle)) NULL else cf_wrap(subtitle, 125),
                  x = NULL, y = NULL) +
    cf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 4.6),
                   axis.text.y = ggplot2::element_text(size = 5))
}

# F3-C1 variant A: SUS-RES primary matrix with the other contrasts as strips.
cfp_render_gsea_atlas_susres <- function(panel, svg_path, source_csv) {
  d <- cfp_gsea_theme_table(panel)
  cells <- cfp_theme_cells(d)
  main <- cfp_gsea_atlas_plot(cells, "SUS - RES",
    "Ranked GSEA program atlas across canonical spatial contexts",
    paste0("cell = median NES of the theme's constituent canonical GO-BP terms; ",
           "dot = number of those terms with canonical FDR < 0.05. ",
           "qc_review themes are QC context, not claims"))
  strip <- cfp_gsea_atlas_plot(cells, c("RES - CON", "SUS - CON"), NULL,
    "reference contrasts")
  strip <- strip + ggplot2::theme(legend.position = "none")
  p <- if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(main, strip, ncol = 1, heights = c(1, 0.92)) else main

  write_csv_safe(cells, source_csv)
  cf_save_panel(p, svg_path, 183, 108)
  invisible(list(status = "ok"))
}

# F3-C1 variant B: three aligned contrast blocks.
cfp_render_gsea_atlas_blocks <- function(panel, svg_path, source_csv) {
  d <- cfp_gsea_theme_table(panel)
  cells <- cfp_theme_cells(d)
  p <- cfp_gsea_atlas_plot(cells, c("RES - CON", "SUS - CON", "SUS - RES"),
    "Ranked GSEA program atlas, three aligned contrast blocks",
    paste0("same evidence as the SUS-RES atlas, laid out so the contrast ",
           "geometry is directly comparable"))
  write_csv_safe(cells, source_csv)
  cf_save_panel(p, svg_path, 183, 92)
  invisible(list(status = "ok"))
}

# F3-C2. Global WGCNA module effects, REUSED verbatim, with side annotations.
cfp_render_wgcna_annotated <- function(panel, svg_path, source_csv) {
  eff <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(eff, as.character(unlist(panel$required_columns)), panel$id)
  deps <- as.character(unlist(panel$input_dependencies))
  aff <- cf_read_csv(repo_path(deps[1]), required = FALSE)
  bil <- cf_read_csv(repo_path(deps[2]), required = FALSE)

  lab_col <- intersect(c("module_display_label", "ModuleLabel_Final", "module_label"),
                       names(eff))[1]
  sm_col <- intersect(c("supermodule_label_for_module", "Supermodule_CleanPlotLabel",
                        "SupermoduleID", "supermodule_id"), names(eff))[1]
  e <- data.frame(
    dataset = eff$dataset, module_id = eff$module_id, contrast = eff$contrast,
    estimate = eff$estimate, tier_specific_fdr = eff$tier_specific_fdr,
    model_valid_for_inference = eff$model_valid_for_inference,
    module_label = if (is.na(lab_col)) eff$module_id else eff[[lab_col]],
    supermodule = if (is.na(sm_col)) NA_character_ else as.character(eff[[sm_col]]),
    stringsAsFactors = FALSE)

  if (!is.null(aff)) {
    a <- aff[aff$dataset == "neuron_neuropil", , drop = FALSE]
    ext <- intersect(c("external_celltype_all", "external_celltype_core_kME06"), names(a))[1]
    e$external_celltype <- if (is.na(ext)) NA_character_ else
      a[[ext]][match(e$module_id, a$ModuleID)]
    e$bilateral_support <- if ("bilateral_reproducibility_class" %in% names(a))
      a$bilateral_reproducibility_class[match(e$module_id, a$ModuleID)] else NA_character_
  }
  if (!is.null(bil) && all(is.na(e$bilateral_support))) {
    b <- bil[bil$dataset == "neuron_neuropil" & bil$level == "module", , drop = FALSE]
    e$bilateral_support <- b$bilateral_reproducibility_class[match(e$module_id, b$endpoint_id)]
  }

  e$fdr_symbol <- ifelse(is.finite(e$tier_specific_fdr) & e$tier_specific_fdr < 0.05,
                         "*", "")
  ord <- unique(e[order(e$supermodule, e$module_id), c("module_id", "module_label")])
  e$module_label <- factor(e$module_label, levels = rev(ord$module_label))
  lim <- max(abs(e$estimate), na.rm = TRUE)

  main <- ggplot2::ggplot(e, ggplot2::aes(contrast, module_label, fill = estimate)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = fdr_symbol), size = 2.6, vjust = 0.78) +
    ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                                  midpoint = 0, limits = c(-lim, lim),
                                  name = "model effect") +
    ggplot2::labs(
      title = cf_wrap("Global WGCNA module effects with program and context annotation", 52),
      subtitle = cf_wrap(paste0("colour = model effect, reused verbatim from the canonical ",
                                "panel source. * marks tier-specific FDR < 0.05; colour alone ",
                                "is NOT significance"), 60),
      x = NULL, y = NULL) +
    cf_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4.6),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 4.6))

  ann <- unique(e[, c("module_label", "supermodule", "external_celltype",
                      "bilateral_support")])
  annl <- rbind(
    data.frame(module_label = ann$module_label, track = "supermodule",
               value = as.character(ann$supermodule)),
    data.frame(module_label = ann$module_label, track = "cell-type affinity",
               value = as.character(ann$external_celltype)),
    data.frame(module_label = ann$module_label, track = "bilateral",
               value = as.character(ann$bilateral_support)))
  annl$value[is.na(annl$value) | annl$value == ""] <- "not available"
  # A discrete colour legend over this many distinct annotation values swamps
  # the panel at assembly scale, so the annotation is labelled IN the tile and
  # carries no legend. The full strings stay in the source data.
  annl$short <- vapply(annl$value, function(v) {
    v <- gsub("_", " ", v)
    if (nchar(v) <= 10) v else paste0(substr(v, 1, 9), "…")
  }, character(1))
  side <- ggplot2::ggplot(annl, ggplot2::aes(track, module_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white",
                       linewidth = 0.3, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = short), size = 1.35, colour = "black") +
    ggplot2::scale_fill_viridis_d(option = "D", alpha = 0.35, na.value = "grey92") +
    ggplot2::labs(x = NULL, y = NULL) +
    cf_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 4.6),
                   panel.grid = ggplot2::element_blank())

  p <- if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(main, side, nrow = 1, widths = c(1, 0.9)) else main

  e$descriptive_geometry_note <- paste0(
    "colour encodes the model effect only; the descriptive RES > CON > SUS ",
    "ordering is never drawn as statistical support")
  write_csv_safe(e, source_csv)
  cf_save_panel(p, svg_path, 110, 82)
  invisible(list(status = "ok"))
}

# F3-C3. Cross-method program convergence. An evidence inventory, not a test.
cfp_render_program_convergence <- function(panel, svg_path, source_csv) {
  d <- cfp_gsea_theme_table(panel)
  cells <- cfp_theme_cells(d)
  deps <- as.character(unlist(panel$input_dependencies))
  map <- cf_read_csv(repo_path(deps[1]), required = FALSE)
  eff <- cf_read_csv(repo_path(deps[2]), required = FALSE)
  bil <- cf_read_csv(repo_path(deps[3]), required = FALSE)
  aff <- cf_read_csv(repo_path(deps[4]), required = FALSE)

  themes <- unique(cells[, c("theme_id", "manuscript_theme", "theme_role")])
  rows <- list()
  for (i in seq_len(nrow(themes))) {
    tid <- themes$theme_id[i]
    z <- cells[cells$theme_id == tid, , drop = FALSE]

    # 1. ranked GSEA support
    n_units_fdr <- sum(z$has_FDR_support & z$contrast == "SUS - RES")
    rows[[length(rows) + 1L]] <- data.frame(
      theme_id = tid, program = themes$manuscript_theme[i],
      theme_role = themes$theme_role[i], evidence = "ranked GSEA",
      value = n_units_fdr,
      label = paste0(n_units_fdr, " units"),
      source_artifact = "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
      stringsAsFactors = FALSE)

    # 2. WGCNA module support, via the curated theme -> entity mapping
    n_mod <- if (is.null(map)) NA_integer_ else
      sum(map$theme_id == tid & map$approved_for_manuscript_interpretation %in% TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      theme_id = tid, program = themes$manuscript_theme[i],
      theme_role = themes$theme_role[i], evidence = "WGCNA module support",
      value = n_mod, label = paste0(n_mod, " approved"),
      source_artifact = "config/gsea_wgcna_theme_module_mapping.csv",
      stringsAsFactors = FALSE)

    # 3. bilateral reproducibility of the mapped modules.
    # Module ids are DATASET-SCOPED: WGCNA_m05 exists in all three datasets and
    # they are unrelated modules, so the join must carry dataset or a
    # neuropil-only mapping silently counts soma and microglia modules too.
    n_bil <- NA_integer_
    if (!is.null(map) && !is.null(bil)) {
      mm <- map[map$theme_id == tid & map$entity_level == "module", , drop = FALSE]
      key <- paste(mm$dataset, mm$entity_id)
      b <- bil[bil$level == "module" &
                 paste(bil$dataset, bil$endpoint_id) %in% key, , drop = FALSE]
      n_bil <- sum(b$bilateral_reproducibility_class == "reproducible_level_and_pattern")
    }
    rows[[length(rows) + 1L]] <- data.frame(
      theme_id = tid, program = themes$manuscript_theme[i],
      theme_role = themes$theme_role[i], evidence = "bilateral reproducibility",
      value = n_bil, label = paste0(n_bil, " reproducible"),
      source_artifact = "WGCNA_module_bilateral_reproducibility.csv",
      stringsAsFactors = FALSE)

    # 4. external cell-type context for the mapped modules
    n_ext <- NA_integer_
    if (!is.null(map) && !is.null(aff)) {
      mm <- map[map$theme_id == tid & map$entity_level == "module", , drop = FALSE]
      key <- paste(mm$dataset, mm$entity_id)
      a <- aff[paste(aff$dataset, aff$ModuleID) %in% key, , drop = FALSE]
      n_ext <- if ("external_any_significant" %in% names(a))
        sum(a$external_any_significant %in% TRUE) else NA_integer_
    }
    rows[[length(rows) + 1L]] <- data.frame(
      theme_id = tid, program = themes$manuscript_theme[i],
      theme_role = themes$theme_role[i], evidence = "external cell-type context",
      value = n_ext, label = paste0(n_ext, " modules"),
      source_artifact = "WGCNA_module_spatial_cell_affinity.csv",
      stringsAsFactors = FALSE)

    # 5. spatial breadth of the GSEA signal
    n_ctx <- length(unique(z$spatial_unit[z$has_FDR_support]))
    rows[[length(rows) + 1L]] <- data.frame(
      theme_id = tid, program = themes$manuscript_theme[i],
      theme_role = themes$theme_role[i], evidence = "spatial breadth",
      value = n_ctx, label = paste0(n_ctx, " contexts"),
      source_artifact = "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
      stringsAsFactors = FALSE)
  }
  conv <- dplyr::bind_rows(rows)
  conv$evidence <- factor(conv$evidence,
    levels = c("ranked GSEA", "WGCNA module support", "bilateral reproducibility",
               "external cell-type context", "spatial breadth"))
  conv$program <- factor(conv$program,
                         levels = rev(sort(unique(conv$program))))
  conv$present <- is.finite(conv$value) & conv$value > 0

  p <- ggplot2::ggplot(conv, ggplot2::aes(evidence, program)) +
    ggplot2::geom_tile(ggplot2::aes(fill = present), colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 1.8) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#CFE3D4", "FALSE" = "grey94"),
                               guide = "none") +
    ggplot2::labs(
      title = "Cross-method program convergence",
      subtitle = paste0("an evidence inventory, NOT a statistical test. Each cell names ",
                        "its source artifact in the accompanying source data; evidence ",
                        "is never collapsed into one score"),
      x = NULL, y = NULL) +
    cf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 4.8),
                   axis.text.y = ggplot2::element_text(size = 5))

  conv$interpretation <- "evidence inventory; no cell is a significance call"
  write_csv_safe(conv, source_csv)
  cf_save_panel(p, svg_path, 89, 72)
  invisible(list(status = "ok"))
}

# F3-C4. Selected QC-clean program examples, chosen by a recorded rule.
cfp_render_program_examples <- function(panel, svg_path, source_csv) {
  d <- cfp_gsea_theme_table(panel)
  cells <- cfp_theme_cells(d)

  # PRESPECIFIED SELECTION RULE, recorded in the source data:
  #   (1) theme_role must be primary (qc_review themes are excluded);
  #   (2) the theme must have genuine ranked-GSEA FDR support in SUS - RES;
  #   (3) rank by the number of spatial contexts with support, so the example
  #       is not chosen for its biological appeal;
  #   (4) keep the top 3 distinct programs.
  rule <- paste0("primary themes only; require >=1 constituent GO term with ",
                 "canonical FDR < 0.05 in SUS - RES; rank by number of spatial ",
                 "contexts with such support; keep the top 3. CA2-SLM is not ",
                 "privileged and no selection uses the DAP burden")
  sr <- cells[cells$contrast == "SUS - RES" & cells$theme_role == "primary", , drop = FALSE]
  brk <- stats::aggregate(has_FDR_support ~ theme_id + manuscript_theme, data = sr, FUN = sum)
  names(brk)[3] <- "n_contexts_with_support"
  brk <- brk[brk$n_contexts_with_support > 0, , drop = FALSE]
  brk <- brk[order(-brk$n_contexts_with_support, brk$manuscript_theme), , drop = FALSE]
  keep <- utils::head(brk$theme_id, 3)

  z <- cells[cells$theme_id %in% keep & cells$contrast == "SUS - RES", , drop = FALSE]
  z$spatial_unit <- cfp_unit_order(z$spatial_unit)
  z$program <- factor(z$manuscript_theme,
                      levels = brk$manuscript_theme[brk$theme_id %in% keep])

  p <- ggplot2::ggplot(z, ggplot2::aes(spatial_unit, median_NES)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
    ggplot2::geom_col(ggplot2::aes(fill = has_FDR_support), width = 0.68) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#23384D", "FALSE" = "grey82"),
      labels = c("TRUE" = "constituent GO term at FDR < 0.05", "FALSE" = "no FDR support"),
      name = NULL) +
    ggplot2::facet_grid(program ~ dataset, scales = "free_x", space = "free_x") +
    ggplot2::labs(
      title = "Selected QC-clean program examples (SUS - RES)",
      subtitle = paste0("selection rule is prespecified and recorded in the source ",
                        "data; bar height is the median NES of constituent terms"),
      x = NULL, y = "median NES") +
    cf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 4.4),
                   strip.text.y = ggplot2::element_text(size = 4.2))

  z$selection_rule <- rule
  z$no_running_enrichment_note <- paste0(
    "no running-enrichment curve is drawn: no gseaResult object is stored in ",
    "the repository and re-running GSEA is not permitted in a renderer")
  write_csv_safe(z, source_csv)
  cf_save_panel(p, svg_path, 89, 84)
  invisible(list(status = "ok"))
}

# F3-C5. Stress effects relative to baseline spatial affinity.
cfp_render_stress_identity <- function(panel, svg_path, source_csv) {
  cmp <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(cmp, as.character(unlist(panel$required_columns)), panel$id)
  atlas <- cf_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]),
                       required = FALSE)

  cmp$subset_label <- gsub("_", " ", cmp$subset)
  cmp$subset_label <- factor(cmp$subset_label, levels = rev(cmp$subset_label))
  cmp$pct <- 100 * cmp$fraction_outside_baseline_affinity
  cmp$primary <- cmp$subset == "CA2_SLM_robustness_qualified"

  bars <- ggplot2::ggplot(cmp, ggplot2::aes(pct, subset_label)) +
    ggplot2::geom_col(ggplot2::aes(fill = primary), width = 0.66) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/%d", effect_outside_baseline_affinity,
                                                    n_hits)),
                       hjust = -0.12, size = 1.9) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#23384D", "FALSE" = "#9E9A92"),
                               guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, 118)) +
    ggplot2::labs(
      title = cf_wrap("Stress effects sit outside the protein's dominant baseline spatial niche", 58),
      subtitle = cf_wrap(paste0("share of FDR-supported proteins whose strongest effect falls ",
                                "outside their own baseline affinity; robustness-qualified ",
                                "subset highlighted"), 66),
      x = "% outside baseline affinity", y = NULL) +
    cf_theme()

  p <- bars
  if (!is.null(atlas)) {
    a <- atlas[atlas$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
    a$claim <- ifelse(is.na(a$QC_claimability) | a$QC_claimability == "",
                      "not_audited", a$QC_claimability)
    ins <- ggplot2::ggplot(a, ggplot2::aes(baseline_rank_of_strongest_effect_unit,
                                           fill = claim)) +
      ggplot2::geom_bar(width = 0.82) +
      ggplot2::scale_fill_manual(values = c(
        claimable = "#23384D", claimable_with_caveat = "#6E7FA0",
        not_claimable = "#D7303F", not_evaluable = "#C9C6BF",
        not_audited = "grey80"), name = NULL) +
      ggplot2::scale_x_continuous(breaks = 1:10) +
      ggplot2::labs(subtitle = "baseline rank of the strongest effect, by claim status",
                    x = "baseline rank (1 = the protein's own peak unit)",
                    y = "proteins") +
      cf_theme() + ggplot2::theme(legend.text = ggplot2::element_text(size = 4))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(bars, ins, ncol = 1, heights = c(1, 0.95))
    }
  }

  cmp$interpretation_note <- paste0(
    "this is where an effect sits relative to the protein's own baseline ",
    "ranking. It is NOT protein redistribution. The extreme rank-10 tail did ",
    "not survive the QC audit and is deliberately not emphasised")
  write_csv_safe(cmp, source_csv)
  cf_save_panel(p, svg_path, 89, 92)
  invisible(list(status = "ok"))
}

# Optional candidate replacement for the canonical DAP-count panel.
cfp_render_dap_status <- function(panel, svg_path, source_csv) {
  counts <- cf_read_csv(repo_path(panel$primary_source))
  cf_require_columns(counts, as.character(unlist(panel$required_columns)), panel$id)
  deps <- as.character(unlist(panel$input_dependencies))
  member <- cf_read_csv(repo_path(deps[1]), required = FALSE)
  atlas <- cf_read_csv(repo_path(deps[2]), required = FALSE)

  # canonical FDR counts are preserved exactly; status is overlaid, never
  # substituted, and only CA2-SLM carries an audited status
  base <- unique(counts[, c("dataset", "spatial_unit", "n_DAP_FDR05")])
  stat <- NULL
  if (!is.null(member) && !is.null(atlas)) {
    a <- atlas[, c("ProteinGroupID", "QC_claimability")]
    m <- merge(member[, c("dataset", "spatial_unit", "ProteinGroupID")], a,
               by = "ProteinGroupID", all.x = TRUE)
    m$QC_claimability[is.na(m$QC_claimability) | m$QC_claimability == ""] <- "not_audited"
    stat <- as.data.frame(table(dataset = m$dataset, spatial_unit = m$spatial_unit,
                                status = m$QC_claimability), stringsAsFactors = FALSE)
    stat <- stat[stat$Freq > 0, , drop = FALSE]
    names(stat)[names(stat) == "Freq"] <- "n_proteins"
  }
  if (is.null(stat) || !nrow(stat)) {
    stat <- data.frame(dataset = base$dataset, spatial_unit = base$spatial_unit,
                       status = "not_audited", n_proteins = base$n_DAP_FDR05,
                       stringsAsFactors = FALSE)
  }
  stat <- stat[stat$n_proteins > 0, , drop = FALSE]
  stat$plot_unit <- paste0(sub("neuron_", "", stat$dataset), " ", stat$spatial_unit)
  stat$status <- factor(stat$status,
    levels = c("claimable", "claimable_with_caveat", "not_claimable",
               "not_evaluable", "not_audited"))

  p <- ggplot2::ggplot(stat, ggplot2::aes(stats::reorder(plot_unit, n_proteins,
                                                         FUN = sum),
                                          n_proteins, fill = status)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      claimable = "#23384D", claimable_with_caveat = "#6E7FA0",
      not_claimable = "#D7303F", not_evaluable = "#C9C6BF",
      not_audited = "grey80"), name = NULL, drop = FALSE) +
    ggplot2::labs(
      title = "Spatial DAP counts with interpretation status",
      subtitle = paste0("canonical BH FDR < 0.05 counts are preserved exactly; ",
                        "status is overlaid, never substituted. Only CA2-SLM has ",
                        "been audited - every other unit reads not_audited"),
      x = NULL, y = "FDR-supported proteins") +
    cf_theme()

  out <- merge(stat, base, by = c("dataset", "spatial_unit"), all.x = TRUE)
  out$canonical_counts_preserved <- TRUE
  out$status_source <- "protein_spatial_cell_affinity.csv QC_claimability (CA2-SLM audit only)"
  write_csv_safe(out, source_csv)
  cf_save_panel(p, svg_path, 89, 62)
  invisible(list(status = "ok"))
}
