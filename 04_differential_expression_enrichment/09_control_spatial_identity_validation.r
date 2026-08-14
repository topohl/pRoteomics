# Figure 2 control-anatomy validation only; does not consume or modify stress/WGCNA outputs.
source("R/paths.R"); source("R/dataset_config.R"); source("R/dataset_inputs.R")
source("R/validation_utils.R"); source("R/qc_exploration_utils.R")
source("R/protein_group_enrichment_utils.R"); source("R/control_spatial_identity_utils.R")
suppressPackageStartupMessages({ library(limma); library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2) })

control_spatial_render_figures <- function(ks, go, source_data, figures) {
  publication_contrast_label <- function(x) {
    labels <- control_spatial_contrast_label(x)
    labels[x == "DG_neuropil_vs_mean_non_DG_regions"] <- "DG neuropil vs non-DG"
    labels
  }
  ks <- ks[ks$internal_contrast %in% control_spatial_figure2e_external_contrasts(), , drop = FALSE]
  ks <- control_spatial_apply_signature_fdr(ks)
  ks$panel_contract <- "Figure 2e: External anatomical validation against matched or explicitly approximate Kaulich regional / CA1-stratum reference signatures."
  ks$matched_exactly <- ks$match_type == "exact"
  ks$external_reference_context <- ks$external_signature == "SP" & !ks$matched_exactly
  utils::write.csv(ks, source_data[["e"]], row.names = FALSE)
  go_display <- control_spatial_select_go_display(
    go, max_terms = 2L, contrasts = control_spatial_figure2f_display_contrasts()
  )
  if (nrow(go_display)) {
    go_display$internal_contrast_label <- publication_contrast_label(go_display$contrast)
    go_display$semantic_simplification_applied <- FALSE
    go_display$display_selection_rule <- "adjusted p-value, descending NES, GO ID; maximum two per contrast"
    go_display$panel_contract <- "Figure 2f: Internal CON-only anatomical GO-BP differentiation."
    go_display$interpretation_note <- ifelse(
      go_display$contrast %in% c("DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers"),
      "Internal DG layer contrast; no matched DG-layer Kaulich reference signature available.",
      "Internal CON-only anatomical GO-BP differentiation."
    )
  }
  utils::write.csv(go_display, source_data[["f"]], row.names = FALSE)

  kaulich_plot_data <- control_spatial_figure2e_plot_data(ks)
  if (nrow(kaulich_plot_data)) {
    kaulich_plot_data$internal_contrast_label <- publication_contrast_label(
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
        low = "#3568A8", mid = "white", high = "#C43C39",
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
          barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.5, "mm")
        ),
        colour = guide_legend(
          order = 2, title.position = "top",
          override.aes = list(shape = 21, fill = "white", size = 3, stroke = 1)
        )
      ) +
      labs(
        x = NULL, y = NULL, tag = "a",
        subtitle = "External Kaulich validation; CA2/3 correspondence is approximate; CA1 SP is reference-only."
      ) +
      theme_minimal(base_size = 8.5) +
      theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#BDBDBD", linewidth = 0.35),
        panel.spacing.y = grid::unit(4.5, "mm"),
        strip.background = element_rect(fill = "#F3F3F3", colour = NA),
        strip.text = element_text(face = "bold", hjust = 0, size = 8.3, margin = margin(3, 4, 3, 4)),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 7.2, margin = margin(t = 3)),
        axis.text.y = element_text(size = 7.5),
        plot.subtitle = element_text(size = 7.7, colour = "#444444", margin = margin(b = 6)),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.spacing.x = grid::unit(3, "mm"),
        plot.margin = margin(8, 8, 6, 8)
      )
    ggsave(figures[["e"]], figure2e, width = 7.2, height = 7.4)
  } else {
    control_spatial_write_empty_state_svg(
      figures[["e"]], "No Kaulich signatures met the prespecified mapping threshold"
    )
  }

  if (nrow(go_display)) {
    contrast_order <- control_spatial_figure2f_display_contrasts()
    go_display$internal_contrast_plot_label <- vapply(
      go_display$internal_contrast_label,
      function(z) paste(strwrap(z, width = 21), collapse = "\n"),
      character(1)
    )
    go_display$internal_contrast_plot_label <- factor(
      go_display$internal_contrast_plot_label,
      levels = vapply(
        publication_contrast_label(contrast_order),
        function(z) paste(strwrap(z, width = 21), collapse = "\n"),
        character(1)
      )
    )
    go_display$term_plot_key <- paste(go_display$contrast, go_display$ID, sep = "\r")
    go_display$term_plot_key <- factor(
      go_display$term_plot_key, levels = rev(unique(go_display$term_plot_key))
    )
    term_labels <- stats::setNames(
      vapply(
        go_display$Description,
        function(z) paste(strwrap(z, width = 32), collapse = "\n"),
        character(1)
      ),
      as.character(go_display$term_plot_key)
    )
    go_display$minus_log10_adjusted_p <- -log10(
      pmax(go_display$p_adjust, .Machine$double.xmin)
    )
    figure2f <- ggplot(
      go_display,
      aes(NES, term_plot_key, colour = NES, size = minus_log10_adjusted_p)
    ) +
      geom_point(alpha = 0.92) +
      facet_wrap(~internal_contrast_plot_label, ncol = 3, scales = "free_y") +
      scale_y_discrete(labels = term_labels) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.1))) +
      scale_colour_gradient(low = "#4C78A8", high = "#C43C39", name = "NES") +
      scale_size_continuous(
        name = expression(-log[10]("adjusted p-value")), range = c(2.3, 5.2)
      ) +
      guides(
        colour = guide_colourbar(
          order = 1, title.position = "top",
          barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.5, "mm")
        ),
        size = guide_legend(order = 2, title.position = "top", nrow = 1)
      ) +
      labs(
        x = "Normalized enrichment score (NES)", y = NULL, tag = "b",
        subtitle = "Internal CON-only anatomical GO-BP differentiation."
      ) +
      theme_minimal(base_size = 8) +
      theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major.x = element_line(colour = "#E7E7E7", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = grid::unit(4, "mm"),
        strip.background = element_rect(fill = "#F3F3F3", colour = NA),
        strip.text = element_text(face = "bold", size = 7.3, margin = margin(3, 3, 3, 3)),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 6.4, lineheight = 0.9),
        axis.text.x = element_text(size = 6.8),
        axis.title.x = element_text(size = 7.5, margin = margin(t = 5)),
        plot.subtitle = element_text(size = 7.4, colour = "#444444", margin = margin(b = 5)),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 7.2),
        legend.text = element_text(size = 6.8),
        legend.spacing.x = grid::unit(3, "mm"),
        plot.margin = margin(8, 8, 6, 8)
      )
    ggsave(figures[["f"]], figure2f, width = 9.2, height = 8.2)
  } else {
    control_spatial_write_empty_state_svg(
      figures[["f"]], "No GO-BP terms met the prespecified display criteria"
    )
  }
  invisible(list(figure2e_source = ks, figure2f_source = go_display))
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
out <- function(kind, name) file.path(root, "results", kind, "04_differential_expression_enrichment", "control_spatial_identity_validation", "global", name)
tables <- c(protein=out("tables","anatomical_protein_contrasts.csv"), mapping=out("tables","kaulich_signature_mapping.csv"), kaulich=out("tables","kaulich_signature_gsea.csv"), go=out("tables","control_anatomical_go_bp_gsea.csv"), matching_audit=out("tables","figure2e_matching_audit.csv"), status=out("tables","analysis_status.csv"))
source_data <- c(e=out("source_data","figure2e_source_data.csv"), f=out("source_data","figure2f_source_data.csv")); figures <- c(e=out("figures","figure2e_kaulich_validation.svg"), f=out("figures","figure2f_control_anatomical_GO.svg")); manifest <- out("logs","run_manifest.yml")
dir.create(dirname(tables[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(source_data[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(figures[[1]]), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(manifest), recursive=TRUE, showWarnings=FALSE)

render_only <- tolower(trimws(Sys.getenv("PROTEOMICS_CONTROL_SPATIAL_RENDER_ONLY", ""))) %in% c("1", "true", "yes")
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
    go <- as.data.frame(clusterProfiler::gseGO(ranked,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP",pvalueCutoff=1,verbose=FALSE))
    if(nrow(go)) {
      go$dataset<-dataset;go$contrast<-ct$name;go$status<-"completed"
      go$p_value <- go$pvalue; go$p_adjust <- go$p.adjust
      if (!"setSize" %in% names(go)) stop("Successful GO-BP GSEA rows are missing setSize.", call. = FALSE)
      go$mapped_unique_genes <- go$setSize
      all_go[[length(all_go)+1L]] <- go
    } else all_go[[length(all_go)+1L]]<-data.frame(dataset=dataset,contrast=ct$name,status="completed_zero_terms")
    statuses[[length(statuses)+1L]] <- control_spatial_empty_status(dataset,ct$name,"completed",paste0("duplicateCorrelation=",round(corfit$consensus,4),"; hemisphere=",if(d$hemisphere_included)"included" else d$hemisphere_omission_reason)); details[[length(details)+1L]]<-list(dataset=dataset,formula=if(d$hemisphere_included)"abundance ~ 0 + anatomical_unit + hemisphere" else "abundance ~ 0 + anatomical_unit",n_samples=ncol(mat),n_animals=length(unique(meta$AnimalID)),contrast=ct$name,weights=ct$weights,duplicate_correlation=corfit$consensus)
  }
}
largest_mapped_signature_size <- max(vapply(kaulich_jobs, function(z) length(z$mapped), integer(1)))
min_gs_size <- 5L
max_gs_size <- control_spatial_signature_max_gs_size(largest_mapped_signature_size)
message("Running Kaulich signature GSEA")
for (job in kaulich_jobs) {
  if (length(job$mapped) < min_gs_size) {
    all_ks[[length(all_ks) + 1L]] <- cbind(
      job$prefix,
      data.frame(
        status = "skipped_lt5_mapped_genes",
        NES = NA_real_, p_value = NA_real_, p_adjust = NA_real_,
        leading_edge_genes = NA_character_, stringsAsFactors = FALSE
      )
    )
    next
  }
  z <- clusterProfiler::GSEA(
    job$ranked,
    TERM2GENE = data.frame(term = job$external_signature, gene = job$mapped),
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    pvalueCutoff = 1,
    verbose = FALSE
  )
  q <- as.data.frame(z)
  all_ks[[length(all_ks) + 1L]] <- if (nrow(q)) {
    cbind(
      job$prefix,
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
      job$prefix,
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
