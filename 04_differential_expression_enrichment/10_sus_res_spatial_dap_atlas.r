#!/usr/bin/env Rscript

# Manuscript-facing spatial atlas of SUS-RES differential protein changes and
# GO-ID ontology-mapped ranked-GSEA themes.
#
# Direction-aware Jaccard for units A and B is
#   |{ProteinGroupID|direction}_A intersect {ProteinGroupID|direction}_B| /
#   |{ProteinGroupID|direction}_A union     {ProteinGroupID|direction}_B|.
# Unsigned Jaccard uses ProteinGroupID alone. Both sets are first restricted to
# ProteinGroupIDs with finite FDR and formal SUS-RES effect in both units.

source("R/paths.R")
source("R/dataset_config.R")
source("R/plotting_nature.R")
source("R/protein_group_enrichment_utils.R")
source("R/enrichment_io.R")
source("R/manuscript_go_theme_utils.R")
source("R/sus_res_spatial_dap_atlas_utils.R")
source("R/sus_res_biological_audit_workbook.R")

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- file.path("sus_res_spatial_dap_atlas", "global")
DATASETS <- valid_datasets()
AUDIT_ONLY <- "--audit-only" %in% commandArgs(trailingOnly = TRUE)
RENDER_ONLY <- "--render-only" %in% commandArgs(trailingOnly = TRUE)
DRY_RUN <- is_dry_run()

if (AUDIT_ONLY && RENDER_ONLY) {
  stop("Choose either --audit-only or --render-only, not both.", call. = FALSE)
}

root_parts <- c(MODULE_ID, "sus_res_spatial_dap_atlas", "global")
TABLE_DIR <- do.call(path_results, as.list(c("tables", root_parts)))
SOURCE_DIR <- do.call(path_results, as.list(c("source_data", root_parts)))
FIGURE_DIR <- do.call(path_results, as.list(c("figures", root_parts)))
LOG_DIR <- do.call(path_results, as.list(c("logs", root_parts)))
REPORT_DIR <- do.call(path_results, as.list(c("reports", root_parts)))

files <- list(
  counts = file.path(TABLE_DIR, "sus_res_dap_counts.csv"),
  membership = file.path(SOURCE_DIR, "sus_res_dap_membership.csv"),
  similarity = file.path(TABLE_DIR, "sus_res_pairwise_similarity.csv"),
  go = file.path(TABLE_DIR, "sus_res_go_bp_ora.csv"),
  go_display = file.path(SOURCE_DIR, "sus_res_go_display_terms.csv"),
  input_audit = file.path(SOURCE_DIR, "sus_res_input_audit.csv"),
  status = file.path(TABLE_DIR, "analysis_status.csv"),
  panel_a_source = file.path(SOURCE_DIR, "panel_a_sus_res_dap_counts.csv"),
  panel_b_source = file.path(SOURCE_DIR, "panel_b_sus_res_dap_jaccard.csv"),
  panel_c_ora_source = file.path(SOURCE_DIR, "panel_c_sus_res_dap_GO_BP.csv"),
  panel_c_source = file.path(SOURCE_DIR, "panel_c_sus_res_ranked_GSEA_themes.csv"),
  panel_c_source_compatibility_alias = file.path(SOURCE_DIR, "panel_c_sus_res_ranked_GSEA_programs.csv"),
  ranked_theme_inventory = file.path(SOURCE_DIR, "sus_res_ranked_GSEA_theme_inventory.csv"),
  theme_overview_source = file.path(SOURCE_DIR, "sus_res_manuscript_theme_overview.csv"),
  theme_overview = file.path(TABLE_DIR, "sus_res_manuscript_theme_overview.csv"),
  biological_summary = file.path(TABLE_DIR, "sus_res_manuscript_theme_summary.csv"),
  biological_summary_source = file.path(SOURCE_DIR, "sus_res_manuscript_theme_summary.csv"),
  biological_summary_compatibility_alias = file.path(TABLE_DIR, "sus_res_biological_program_summary.csv"),
  biological_summary_source_compatibility_alias = file.path(SOURCE_DIR, "sus_res_biological_program_summary.csv"),
  biological_audit_workbook = file.path(REPORT_DIR, "sus_res_biological_audit.xlsx"),
  biological_audit_workbook_validation = file.path(REPORT_DIR, "sus_res_biological_audit_validation.csv"),
  manifest = file.path(LOG_DIR, "run_manifest.yml")
)

RANKED_PROGRAM_SOURCE <- path_results(
  "source_data", MODULE_ID, "compareGO_spatial_atlas",
  "source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv"
)
SUPPORTED_GO_AUDIT_SOURCE <- path_results(
  "source_data", MODULE_ID, "compareGO_spatial_atlas",
  "sus_res_supported_go_term_theme_audit.csv"
)
MANUSCRIPT_GO_THEME_REGISTRY <- repo_path("config", "manuscript_go_theme_registry.tsv")

manifest_paths <- setNames(vapply(DATASETS, canonical_clusterprofiler_manifest_path, character(1)), DATASETS)
specs <- do.call(rbind, lapply(DATASETS, function(dataset) {
  sus_res_manifest_contrasts(dataset, manifest_paths[[dataset]], repository_root = repo_root())
}))
rownames(specs) <- NULL

if (isTRUE(DRY_RUN)) {
  message("[DRY-RUN] SUS-RES spatial DAP atlas")
  message("[DRY-RUN] canonical manifests: ", paste(manifest_paths, collapse = "; "))
  message("[DRY-RUN] selected phenotype-within-unit SUS-RES contrasts: ", nrow(specs))
  message("[DRY-RUN] datasets: ", paste(sprintf("%s=%d", DATASETS, table(factor(specs$dataset, levels = DATASETS))), collapse = ", "))
  message("[DRY-RUN] all manifest hashes match current canonical mapped files: ", all(specs$input_hash_matches_manifest))
  message("[DRY-RUN] ranked GO/GSEA Panel C source: ", RANKED_PROGRAM_SOURCE)
  message("[DRY-RUN] ranked GO/GSEA Panel C source exists: ", file.exists(RANKED_PROGRAM_SOURCE))
  quit(status = if (file.exists(RANKED_PROGRAM_SOURCE)) 0 else 1, save = "no")
}

contrast_rows <- vector("list", nrow(specs))
for (i in seq_len(nrow(specs))) {
  raw <- utils::read.csv(specs$input_file[[i]], stringsAsFactors = FALSE, check.names = FALSE)
  normalized <- sus_res_normalize_contrast(raw, specs$comparison[[i]], SUS_RES_DAP_FDR_THRESHOLD)
  mapped <- sus_res_attach_gene_mapping(normalized)$data
  mapped$dataset <- specs$dataset[[i]]
  mapped$dataset_label <- specs$dataset_label[[i]]
  mapped$route_category <- specs$route_category[[i]]
  mapped$route_unit <- specs$route_unit[[i]]
  mapped$spatial_unit <- specs$spatial_unit[[i]]
  mapped$input_file <- specs$input_file[[i]]
  mapped$input_hash <- specs$current_input_hash[[i]]
  contrast_rows[[i]] <- mapped
}
membership_all <- do.call(rbind, contrast_rows)
rownames(membership_all) <- NULL

counts <- sus_res_count_summary(membership_all, MIN_ORA_DAP_GENES)
dataset_order <- match(counts$dataset, DATASETS)
unit_order <- ave(seq_len(nrow(counts)), counts$dataset, FUN = function(idx) {
  match(counts$spatial_unit[idx], anatomical_spatial_unit_levels(counts$spatial_unit[idx]))
})
counts <- counts[order(dataset_order, unit_order, method = "radix"), , drop = FALSE]
rownames(counts) <- NULL

membership_columns <- c(
  "dataset", "dataset_label", "comparison", "route_category", "route_unit", "spatial_unit",
  "ProteinGroupID", "original_identifier", "member_accessions", "member_gene_symbols",
  "official_gene_symbol", "official_entrez_id", "protein_group_gene_annotation_status",
  "protein_group_ambiguity_class", "gene_level_claim_allowed", "protein_level_claim_allowed",
  "mapping_status", "protein_FDR", "protein_p_value", "serialized_effect_column",
  "serialized_effect", "formal_SUS_minus_RES_effect", "is_DAP_FDR05", "DAP_direction",
  "eligible_official_gene_symbol", "eligible_official_entrez_id",
  "gene_mapping_eligibility_status", "gene_mapping_exclusion_reason", "input_file", "input_hash"
)
missing_membership_columns <- setdiff(membership_columns, names(membership_all))
if (length(missing_membership_columns)) {
  stop("DAP membership source data is missing columns: ", paste(missing_membership_columns, collapse = ", "), call. = FALSE)
}
dap_membership <- membership_all[membership_all$is_DAP_FDR05, membership_columns, drop = FALSE]

input_audit <- specs[, c(
  "dataset", "dataset_label", "comparison", "route_category", "route_unit", "spatial_unit",
  "manifest_file", "input_file_manifest", "input_file", "input_hash", "current_input_hash",
  "input_hash_matches_manifest", "serialized_left_phenotype", "serialized_right_phenotype",
  "serialized_effect_definition", "formal_effect_definition", "formal_effect_multiplier", "sign_was_flipped"
), drop = FALSE]
input_audit$protein_fdr_threshold <- SUS_RES_DAP_FDR_THRESHOLD
input_audit$absolute_log2fc_threshold <- 0

ora_plan <- sus_res_build_ora_plan(membership_all, MIN_ORA_DAP_GENES)
analysis_status <- ora_plan$summary
analysis_status$analysis_mode <- if (AUDIT_ONLY) "audit_only" else "full"

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!RENDER_ONLY) {
  utils::write.csv(counts, files$counts, row.names = FALSE)
  utils::write.csv(counts, files$panel_a_source, row.names = FALSE)
  utils::write.csv(dap_membership, files$membership, row.names = FALSE)
  utils::write.csv(input_audit, files$input_audit, row.names = FALSE)
  utils::write.csv(analysis_status, files$status, row.names = FALSE)
} else {
  preserved_ora_outputs <- unlist(files[c("go", "go_display", "panel_c_ora_source")], use.names = FALSE)
  missing_preserved_ora <- preserved_ora_outputs[!file.exists(preserved_ora_outputs)]
  if (length(missing_preserved_ora)) {
    stop(
      "--render-only requires the existing DAP-set ORA outputs and will not recompute them. Missing: ",
      paste(missing_preserved_ora, collapse = ", "),
      call. = FALSE
    )
  }
}

audit_print <- counts[, c(
  "dataset", "spatial_unit", "n_tested_ProteinGroupID", "n_DAP_FDR05",
  "n_higher_in_SUS", "n_higher_in_RES",
  "n_gene_level_eligible_DAP_genes_higher_in_SUS",
  "n_gene_level_eligible_DAP_genes_higher_in_RES",
  "higher_in_SUS_meets_MIN_ORA_DAP_GENES",
  "higher_in_RES_meets_MIN_ORA_DAP_GENES"
), drop = FALSE]
message("SUS-RES spatial DAP audit (BH protein FDR <= 0.05; MIN_ORA_DAP_GENES = 10):")
print(audit_print, row.names = FALSE)

if (isTRUE(AUDIT_ONLY)) {
  write_run_manifest(
    files$manifest,
    inputs = c(as.list(manifest_paths), as.list(setNames(specs$input_file, paste(specs$dataset, specs$comparison, sep = "::")))),
    outputs = files[c("counts", "membership", "input_audit", "status", "panel_a_source")],
    parameters = list(
      audit_only = TRUE,
      protein_FDR_threshold = SUS_RES_DAP_FDR_THRESHOLD,
      absolute_log2fc_threshold = 0,
      MIN_ORA_DAP_GENES = MIN_ORA_DAP_GENES,
      formal_effect = "SUS - RES"
    ),
    notes = "Audit-only mode writes counts, DAP membership, provenance, and ORA eligibility status; it does not run GO ORA or write scientific figures."
  )
  message("Audit-only complete. No GO ORA results or scientific figures were written.")
  quit(status = 0, save = "no")
}

required_packages <- c("ggplot2", "scales", "svglite", "cowplot", "patchwork")
if (!RENDER_ONLY) required_packages <- c("AnnotationDbi", "clusterProfiler", "org.Mm.eg.db", required_packages)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required package(s) for figure generation: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

similarity <- sus_res_pairwise_by_dataset(membership_all)
if (!nrow(similarity)) stop("No within-dataset spatial-unit pairs were available for similarity analysis.", call. = FALSE)
similarity <- similarity[order(match(similarity$dataset, DATASETS), similarity$spatial_unit_a, similarity$spatial_unit_b, method = "radix"), , drop = FALSE]
if (!RENDER_ONLY) utils::write.csv(similarity, files$similarity, row.names = FALSE)

empty_go_results <- function() {
  data.frame(
    dataset = character(), dataset_label = character(), comparison = character(),
    route_unit = character(), spatial_unit = character(), direction = character(),
    n_DAP_ProteinGroupID = integer(), n_gene_level_eligible_ProteinGroupID = integer(),
    n_unique_query_genes = integer(), n_universe_genes = integer(),
    ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), geneID = character(),
    Count = integer(), stringsAsFactors = FALSE
  )
}

if (RENDER_ONLY) {
  go_display <- utils::read.csv(files$go_display, stringsAsFactors = FALSE, check.names = FALSE)
  message("[RENDER-ONLY] Reusing existing DAP-set ORA outputs without running clusterProfiler.")
} else {
  org_db <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  org_symbols <- AnnotationDbi::keys(org_db, keytype = "SYMBOL")
  go_rows <- list()
  analysis_status$n_GO_BP_terms_returned <- NA_integer_
  for (i in seq_len(nrow(analysis_status))) {
    if (!identical(analysis_status$ORA_status[[i]], "eligible_for_ORA")) next
    key <- paste(analysis_status$dataset[[i]], analysis_status$spatial_unit[[i]], analysis_status$direction[[i]], sep = "|")
    query <- ora_plan$queries[[key]]
    universe <- ora_plan$universes[[key]]
    invalid <- setdiff(universe, org_symbols)
    if (length(invalid)) {
      stop("Canonical eligible official symbols are absent from the current org.Mm.eg.db for ", key, ": ", paste(head(invalid, 10L), collapse = ", "), call. = FALSE)
    }
    ora <- clusterProfiler::enrichGO(
      gene = query,
      universe = universe,
      OrgDb = org_db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10L,
      maxGSSize = 500L,
      readable = FALSE
    )
    result <- enrichment_result_table(ora, paste("SUS-RES GO-BP ORA", key), required = TRUE)
    analysis_status$n_GO_BP_terms_returned[[i]] <- nrow(result)
    analysis_status$ORA_status[[i]] <- if (nrow(result)) "completed_with_terms" else "completed_zero_terms"
    if (nrow(result)) {
      meta <- analysis_status[i, c(
        "dataset", "dataset_label", "comparison", "route_unit", "spatial_unit", "direction",
        "n_DAP_ProteinGroupID", "n_gene_level_eligible_ProteinGroupID", "n_unique_query_genes", "n_universe_genes"
      ), drop = FALSE]
      meta <- meta[rep(1L, nrow(result)), , drop = FALSE]
      go_rows[[length(go_rows) + 1L]] <- cbind(meta, result, stringsAsFactors = FALSE)
    }
  }
  go_results <- if (length(go_rows)) do.call(rbind, go_rows) else empty_go_results()
  rownames(go_results) <- NULL
  go_display <- sus_res_select_go_display(go_results, max_terms = 2L, fdr_threshold = 0.05)
  utils::write.csv(go_results, files$go, row.names = FALSE)
  utils::write.csv(go_display, files$go_display, row.names = FALSE)
  utils::write.csv(analysis_status, files$status, row.names = FALSE)
}

dataset_levels <- vapply(DATASETS, sus_res_publication_dataset_label, character(1))
plot_unit_levels <- unlist(lapply(DATASETS, function(dataset) {
  units <- counts$spatial_unit[counts$dataset == dataset]
  paste(dataset, anatomical_spatial_unit_levels(units), sep = "|")
}), use.names = FALSE)
plot_unit_labels <- setNames(
  clean_spatial_unit_label(sub("^[^|]+\\|", "", plot_unit_levels)),
  plot_unit_levels
)
dataset_panel_labels <- c(
  Neuropil = "Neuropil",
  Soma = "Soma",
  `Microglia-enriched ROI` = "Microglia-\nenriched ROI"
)
dataset_panel_labeller <- ggplot2::as_labeller(dataset_panel_labels)

panel_a_source <- rbind(
  transform(counts, direction = "Higher in SUS", signed_count = n_higher_in_SUS),
  transform(counts, direction = "Higher in RES", signed_count = -n_higher_in_RES)
)
panel_a_source$dataset_label <- factor(panel_a_source$dataset_label, levels = dataset_levels)
panel_a_source$plot_unit <- factor(paste(panel_a_source$dataset, panel_a_source$spatial_unit, sep = "|"), levels = rev(plot_unit_levels))
panel_a_source$direction <- factor(panel_a_source$direction, levels = c("Higher in RES", "Higher in SUS"))
panel_a_source$panel <- "a"
panel_a_source$statistical_view <- "individual proteins passing BH FDR <= 0.05"
panel_a_source$panel_data_basis <- "FDR_supported_DAPs"
utils::write.csv(panel_a_source, files$panel_a_source, row.names = FALSE)
manuscript_text <- nature_manuscript_text_sizes_pt()

panel_a <- ggplot2::ggplot(panel_a_source, ggplot2::aes(x = signed_count, y = plot_unit, fill = direction)) +
  ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.18) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "black") +
  ggplot2::facet_grid(
    dataset_label ~ ., scales = "free_y", space = "free_y",
    labeller = dataset_panel_labeller
  ) +
  ggplot2::scale_y_discrete(labels = plot_unit_labels) +
  ggplot2::scale_x_continuous(labels = function(x) abs(x), expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
  ggplot2::scale_fill_manual(
    values = c(
      "Higher in RES" = nature_palette("group")[["RES"]],
      "Higher in SUS" = nature_palette("group")[["SUS"]]
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Differential ProteinGroupIDs\n(BH FDR <= 0.05)", y = NULL, fill = NULL
  ) +
  theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", publication_legible = TRUE) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = manuscript_text[["dense"]]),
    axis.title.x = ggplot2::element_text(size = manuscript_text[["axis_title"]]),
    strip.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
    strip.text.y.right = ggplot2::element_text(angle = 0, size = manuscript_text[["normal"]], lineheight = 0.9),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
    plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
  )

panel_b_source <- rbind(
  transform(similarity, plot_unit_x = paste(dataset, spatial_unit_a, sep = "|"), plot_unit_y = paste(dataset, spatial_unit_b, sep = "|")),
  transform(similarity, plot_unit_x = paste(dataset, spatial_unit_b, sep = "|"), plot_unit_y = paste(dataset, spatial_unit_a, sep = "|"))
)
diagonal <- transform(
  counts,
  comparison_a = comparison,
  comparison_b = comparison,
  spatial_unit_a = spatial_unit,
  spatial_unit_b = spatial_unit,
  n_tested_common = n_tested_ProteinGroupID,
  n_DAP_A_in_common_universe = n_DAP_FDR05,
  n_DAP_B_in_common_universe = n_DAP_FDR05,
  n_shared_same_direction = n_DAP_FDR05,
  n_shared_opposite_direction = 0L,
  direction_aware_jaccard = ifelse(n_DAP_FDR05 > 0, 1, NA_real_),
  unsigned_protein_jaccard = ifelse(n_DAP_FDR05 > 0, 1, NA_real_),
  direction_concordance_among_shared = ifelse(n_DAP_FDR05 > 0, 1, NA_real_),
  plot_unit_x = paste(dataset, spatial_unit, sep = "|"),
  plot_unit_y = paste(dataset, spatial_unit, sep = "|")
)
keep_b <- names(panel_b_source)
for (column in setdiff(keep_b, names(diagonal))) diagonal[[column]] <- NA
panel_b_source <- rbind(panel_b_source, diagonal[, keep_b, drop = FALSE])
panel_b_source$dataset_label <- factor(panel_b_source$dataset_label, levels = dataset_levels)
panel_b_source$plot_unit_x <- factor(panel_b_source$plot_unit_x, levels = plot_unit_levels)
panel_b_source$plot_unit_y <- factor(panel_b_source$plot_unit_y, levels = rev(plot_unit_levels))
panel_b_source$panel <- "b"
panel_b_source$statistical_view <- "overlap of FDR-supported directional protein sets"
panel_b_source$panel_data_basis <- "FDR_supported_DAPs"
utils::write.csv(panel_b_source, files$panel_b_source, row.names = FALSE)

build_jaccard_block <- function(dataset) {
  unit_keys <- paste(
    dataset,
    anatomical_spatial_unit_levels(unique(counts$spatial_unit[counts$dataset == dataset])),
    sep = "|"
  )
  block <- panel_b_source[panel_b_source$dataset == dataset, , drop = FALSE]
  block$plot_unit_x <- factor(as.character(block$plot_unit_x), levels = unit_keys)
  block$plot_unit_y <- factor(as.character(block$plot_unit_y), levels = rev(unit_keys))
  block$dataset_label <- factor(block$dataset_label, levels = dataset_levels)
  block_labels <- setNames(clean_spatial_unit_label(sub("^[^|]+\\|", "", unit_keys)), unit_keys)
  ggplot2::ggplot(block, ggplot2::aes(x = plot_unit_x, y = plot_unit_y, fill = direction_aware_jaccard)) +
    ggplot2::geom_tile(width = 0.98, height = 0.98, colour = NA) +
    ggplot2::facet_wrap(~dataset_label, nrow = 1, labeller = dataset_panel_labeller) +
    ggplot2::scale_x_discrete(labels = block_labels, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = block_labels, drop = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = nature_palette("jaccard"), limits = c(0, 1), na.value = "#FFFFFF",
      name = "Jaccard",
      guide = ggplot2::guide_colourbar(
        title.position = "top", barwidth = grid::unit(27, "mm"), barheight = grid::unit(2.1, "mm")
      )
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = manuscript_text[["dense"]]),
      axis.text.y = ggplot2::element_text(size = manuscript_text[["dense"]]),
      strip.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = manuscript_text[["normal"]]),
      legend.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
      plot.margin = ggplot2::margin(0.5, 0.8, 0.5, 0.5)
    )
}
jaccard_units_per_dataset <- vapply(DATASETS, function(dataset) {
  length(unique(counts$spatial_unit[counts$dataset == dataset]))
}, integer(1))
panel_b <- patchwork::wrap_plots(
  lapply(DATASETS, build_jaccard_block), nrow = 1,
  widths = jaccard_units_per_dataset, guides = "collect"
) & ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

if (nrow(go_display)) {
  go_display$dataset_label <- factor(go_display$dataset_label, levels = dataset_levels)
  go_display$plot_unit <- factor(paste(go_display$dataset, go_display$spatial_unit, sep = "|"), levels = plot_unit_levels)
  go_display$minus_log10_GO_FDR <- -log10(pmax(go_display$p.adjust, .Machine$double.xmin))
  go_display$facet_label <- paste(as.character(go_display$dataset_label), go_display$direction, sep = " | ")
  go_display$term_label <- vapply(go_display$Description, function(x) paste(strwrap(x, width = 36), collapse = "\n"), character(1))
  go_display$term_key <- paste(go_display$dataset, go_display$direction, go_display$ID, sep = "|")
  term_order <- unique(go_display$term_key[order(go_display$dataset, go_display$direction, go_display$p.adjust, go_display$Description, method = "radix")])
  term_labels <- setNames(go_display$term_label[match(term_order, go_display$term_key)], term_order)
  go_display$term_key <- factor(go_display$term_key, levels = rev(term_order))
  if (!RENDER_ONLY) utils::write.csv(go_display, files$panel_c_ora_source, row.names = FALSE)
  panel_c_ora <- ggplot2::ggplot(go_display, ggplot2::aes(x = plot_unit, y = term_key)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, fill = minus_log10_GO_FDR, shape = FDR_supported), colour = "black", stroke = 0.3) +
    ggplot2::facet_wrap(~facet_label, scales = "free", ncol = 2) +
    ggplot2::scale_x_discrete(labels = plot_unit_labels) +
    ggplot2::scale_y_discrete(labels = term_labels) +
    ggplot2::scale_fill_gradient(low = "#F4E3A1", high = "#A23B72", name = expression(-log[10]("GO BH FDR"))) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 21),
      labels = c(`FALSE` = "Not FDR-supported", `TRUE` = "FDR-supported"),
      name = "GO-BP support"
    ) +
    ggplot2::scale_size_continuous(range = c(1.4, 4.2), name = "Overlap genes") +
    ggplot2::labs(x = NULL, y = NULL, title = "GO-BP ORA of significant DAP gene sets") +
    theme_nature_dotplot(base_size = 6.5, base_family = "sans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.5), axis.text.y = ggplot2::element_text(size = 5.5))
} else {
  if (!RENDER_ONLY) utils::write.csv(sus_res_empty_go_display(), files$panel_c_ora_source, row.names = FALSE)
  panel_c_ora <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No GO-BP terms met BH FDR <= 0.05\nfor eligible DAP gene sets", size = 3) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = "GO-BP ORA of significant DAP gene sets") +
    theme_nature_base(base_size = 7, base_family = "sans") +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
}

if (!file.exists(RANKED_PROGRAM_SOURCE)) {
  stop(
    "Missing ranked GO/GSEA Panel C source: ", RANKED_PROGRAM_SOURCE,
    ". Run 04_differential_expression_enrichment/07_compareGO_spatial_program_atlas.r first; script 10 does not recompute ranked GO/GSEA statistics.",
    call. = FALSE
  )
}
ranked_program_raw <- utils::read.csv(RANKED_PROGRAM_SOURCE, stringsAsFactors = FALSE, check.names = FALSE)
ranked_theme_inventory <- sus_res_prepare_ranked_program_panel(ranked_program_raw, expected_datasets = DATASETS)
ranked_theme_inventory <- sus_res_attach_dap_counts_to_ranked_programs(ranked_theme_inventory, counts)
biological_summary <- sus_res_biological_inspection_summary(ranked_theme_inventory)
theme_overview <- sus_res_theme_overview(ranked_theme_inventory)
panel_c_source <- ranked_theme_inventory[
  as.character(ranked_theme_inventory$theme_role) == "primary" &
    suppressWarnings(as.integer(ranked_theme_inventory$n_theme_terms_tested)) > 0L,
  , drop = FALSE
]
utils::write.csv(ranked_theme_inventory, files$ranked_theme_inventory, row.names = FALSE)
utils::write.csv(theme_overview, files$theme_overview_source, row.names = FALSE)
utils::write.csv(theme_overview, files$theme_overview, row.names = FALSE)
utils::write.csv(panel_c_source, files$panel_c_source, row.names = FALSE)
utils::write.csv(panel_c_source, files$panel_c_source_compatibility_alias, row.names = FALSE)
utils::write.csv(biological_summary, files$biological_summary, row.names = FALSE)
utils::write.csv(biological_summary, files$biological_summary_source, row.names = FALSE)
utils::write.csv(biological_summary, files$biological_summary_compatibility_alias, row.names = FALSE)
utils::write.csv(biological_summary, files$biological_summary_source_compatibility_alias, row.names = FALSE)
if (!file.exists(SUPPORTED_GO_AUDIT_SOURCE)) {
  stop("Missing supported GO-term audit for the human-facing workbook: ", SUPPORTED_GO_AUDIT_SOURCE, call. = FALSE)
}
supported_go_audit <- utils::read.csv(SUPPORTED_GO_AUDIT_SOURCE, stringsAsFactors = FALSE, check.names = FALSE)
manuscript_registry <- read_manuscript_go_theme_registry(MANUSCRIPT_GO_THEME_REGISTRY)
write_sus_res_biological_audit_workbook(
  overview = theme_overview,
  detail = biological_summary,
  supported_go_audit = supported_go_audit,
  dap_counts = counts,
  registry = manuscript_registry,
  output_file = files$biological_audit_workbook
)
workbook_validation <- validate_sus_res_biological_audit_workbook(files$biological_audit_workbook, biological_summary)
utils::write.csv(workbook_validation, files$biological_audit_workbook_validation, row.names = FALSE)

primary_theme <- as.character(panel_c_source$theme_role) == "primary"
panel_c_axis_source <- panel_c_source[primary_theme, , drop = FALSE]
panel_c_plot_source <- panel_c_source[primary_theme, , drop = FALSE]
if (!nrow(panel_c_plot_source)) {
  stop("Ranked GO/GSEA Panel C source has no ontology-mapped primary manuscript theme rows.", call. = FALSE)
}
panel_c_axis_source <- panel_c_axis_source[order(panel_c_axis_source$display_order, panel_c_axis_source$program, method = "radix"), , drop = FALSE]
program_levels <- unique(as.character(panel_c_axis_source$program))
panel_c_axis_source$program <- factor(panel_c_axis_source$program, levels = rev(program_levels))
panel_c_axis_source$spatial_unit_label <- factor(
  panel_c_axis_source$spatial_unit_label,
  levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(panel_c_axis_source$spatial_unit)))
)
panel_c_axis_source$dataset_compartment_label <- factor(
  panel_c_axis_source$dataset_compartment_label,
  levels = dataset_levels
)
panel_c_plot_source$program <- factor(panel_c_plot_source$program, levels = levels(panel_c_axis_source$program))
panel_c_plot_source$spatial_unit_label <- factor(panel_c_plot_source$spatial_unit_label, levels = levels(panel_c_axis_source$spatial_unit_label))
panel_c_plot_source$dataset_compartment_label <- factor(panel_c_plot_source$dataset_compartment_label, levels = dataset_levels)
panel_c_plot_source$unit_direction_status <- ifelse(
  as.character(panel_c_plot_source$direction_consistency) == "mixed_direction",
  "mixed supported directions", "supported direction not mixed"
)
nes_limit <- suppressWarnings(stats::quantile(abs(panel_c_plot_source$median_NES_all_theme_terms), probs = 0.98, na.rm = TRUE, names = FALSE))
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- suppressWarnings(max(abs(panel_c_plot_source$median_NES_all_theme_terms), na.rm = TRUE))
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- 1
nes_limit <- min(nes_limit, 2.5)
panel_c_descriptive_only <- panel_c_plot_source[!as.logical(panel_c_plot_source$FDR_support_present), , drop = FALSE]
panel_c_supported <- panel_c_plot_source[as.logical(panel_c_plot_source$FDR_support_present), , drop = FALSE]
fdr_size_values <- panel_c_plot_source$representative_minus_log10_FDR[
  is.finite(panel_c_plot_source$representative_minus_log10_FDR)
]
fdr_size_breaks <- if (length(fdr_size_values)) {
  unique(signif(stats::quantile(fdr_size_values, probs = c(0, 0.5, 1), names = FALSE), 2))
} else {
  c(-log10(0.05), 2, 3)
}

panel_c <- ggplot2::ggplot(
  panel_c_axis_source,
  ggplot2::aes(x = spatial_unit_label, y = program)
) +
  ggplot2::geom_blank() +
  ggplot2::geom_point(
    data = panel_c_descriptive_only,
    ggplot2::aes(fill = median_NES_all_theme_terms),
    shape = 21, size = 1.4, colour = "transparent", alpha = 0.45
  ) +
  ggplot2::geom_point(
    data = panel_c_supported,
    ggplot2::aes(fill = median_NES_all_theme_terms, size = representative_minus_log10_FDR, shape = unit_direction_status),
    colour = "black", stroke = 0.35, alpha = 0.96
  ) +
  ggplot2::facet_wrap(
    ~dataset_compartment_label, ncol = 1, scales = "free_x",
    strip.position = "right", labeller = dataset_panel_labeller
  ) +
  ggplot2::scale_fill_gradient2(
    low = nature_palette("signed")[["low"]],
    mid = nature_palette("signed")[["mid"]],
    high = nature_palette("signed")[["high"]], midpoint = 0,
    limits = c(-nes_limit, nes_limit), oob = scales::squish, name = "Median NES"
  ) +
  ggplot2::scale_size_continuous(range = c(2.2, 4.5), breaks = fdr_size_breaks, name = expression(-log[10]("GO FDR"))) +
  ggplot2::scale_shape_manual(values = c("supported direction not mixed" = 21, "mixed supported directions" = 23), name = "FDR-supported\nterm directions") +
  ggplot2::guides(
    fill = ggplot2::guide_colourbar(
      order = 1, title.position = "top",
      barwidth = grid::unit(18, "mm"), barheight = grid::unit(2.1, "mm")
    ),
    size = ggplot2::guide_legend(order = 2, title.position = "top", nrow = 1),
    shape = if (length(unique(panel_c_supported$unit_direction_status)) > 1L) {
      ggplot2::guide_legend(order = 3, title.position = "top", nrow = 1)
    } else {
      "none"
    }
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = manuscript_text[["dense"]]),
    axis.text.y = ggplot2::element_text(size = manuscript_text[["dense"]]),
    strip.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
    strip.text.y.right = ggplot2::element_text(angle = 0, size = manuscript_text[["normal"]], lineheight = 0.9),
    panel.spacing.y = grid::unit(0.45, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(size = manuscript_text[["normal"]]),
    legend.text = ggplot2::element_text(size = manuscript_text[["normal"]]),
    legend.spacing.x = grid::unit(0.6, "mm"),
    plot.margin = ggplot2::margin(0.5, 1.2, 0.5, 0.8)
  )

save_plot_triplet <- function(plot, stem, width_mm, height_mm) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = svglite::svglite, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = png_device, dpi = 300, bg = "white", limitsize = FALSE)
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
panel_c_terms <- if (nrow(go_display)) length(unique(paste(go_display$dataset, go_display$direction, go_display$ID, sep = "|"))) else 0L
panel_c_height <- max(80, min(210, 50 + 5 * panel_c_terms))
jaccard_width_mm <- 5 + sum(14 + 3.8 * jaccard_units_per_dataset)
jaccard_height_mm <- max(60, 22 + 3.8 * max(jaccard_units_per_dataset))
save_plot_triplet(panel_a, file.path(FIGURE_DIR, "panel_a_sus_res_dap_counts"), 50, 48)
save_plot_triplet(panel_b, file.path(FIGURE_DIR, "panel_b_sus_res_dap_jaccard"), jaccard_width_mm, jaccard_height_mm)
if (!RENDER_ONLY) save_plot_triplet(panel_c_ora, file.path(FIGURE_DIR, "panel_c_sus_res_dap_GO_BP"), 183, panel_c_height)
save_plot_triplet(panel_c, file.path(FIGURE_DIR, "panel_c_sus_res_ranked_GSEA_themes"), 71, 80)
# Preserve the established figure stem as a presentation-only compatibility alias.
save_plot_triplet(panel_c, file.path(FIGURE_DIR, "panel_c_sus_res_ranked_GSEA_programs"), 71, 80)

combined_readable <- TRUE
panel_a_combined <- panel_a + ggplot2::labs(title = "DAP counts") + ggplot2::theme(plot.title = ggplot2::element_text(size = manuscript_text[["title"]], face = "bold"))
panel_b_combined <- panel_b + ggplot2::labs(title = "Direction-aware similarity") + ggplot2::theme(plot.title = ggplot2::element_text(size = manuscript_text[["title"]], face = "bold"))
panel_c_combined <- panel_c + ggplot2::theme(plot.title = ggplot2::element_text(size = manuscript_text[["title"]], face = "bold"))
top <- patchwork::wrap_plots(
  panel_a_combined, panel_b_combined,
  nrow = 1, widths = c(0.78, 1.42)
)
combined <- patchwork::wrap_plots(
  top, panel_c_combined,
  ncol = 1, heights = c(1.0, 1.55)
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", size = manuscript_text[["panel_letter"]], face = "bold"))
  )
combined_stem <- file.path(FIGURE_DIR, "figure_sus_res_spatial_dap_atlas_183mm")
save_plot_triplet(combined, combined_stem, 183, 150)
combined_files <- paste0(combined_stem, c(".svg", ".pdf", ".png"))

analysis_status$combined_183mm_figure_readable <- combined_readable
analysis_status$n_panel_c_union_terms <- length(unique(panel_c_plot_source$theme_id))
analysis_status$panel_c_evidence_basis <- "ranked_proteome_wide_GO_GSEA"
analysis_status$panel_c_effect_display <- "descriptive_median_NES_all_mapped_tested_GO_terms_with_original_GO_FDR_support_emphasis"
if (!RENDER_ONLY) utils::write.csv(analysis_status, files$status, row.names = FALSE)

figure_outputs <- list.files(FIGURE_DIR, full.names = TRUE)
write_run_manifest(
  files$manifest,
    inputs = c(
      as.list(manifest_paths),
      as.list(setNames(specs$input_file, paste(specs$dataset, specs$comparison, sep = "::"))),
      ranked_program_source = RANKED_PROGRAM_SOURCE
    ),
  outputs = c(files, as.list(setNames(figure_outputs, basename(figure_outputs)))),
  parameters = list(
      audit_only = FALSE,
      render_only = RENDER_ONLY,
    protein_FDR_threshold = SUS_RES_DAP_FDR_THRESHOLD,
    absolute_log2fc_threshold = 0,
    MIN_ORA_DAP_GENES = MIN_ORA_DAP_GENES,
    formal_effect = "SUS - RES",
    similarity_scope = "within dataset only; pairwise-common tested ProteinGroupID universe",
    GO_ontology = "BP",
    GO_adjustment = "BH",
    GO_query = "eligible unique official genes from directional ProteinGroupIDs with protein BH FDR <= 0.05",
      GO_universe = "eligible unique official genes among all tested ProteinGroupIDs in the same contrast",
      panel_ab_basis = "individual proteins passing BH FDR <= 0.05 and their direction-aware set overlap",
      panel_c_basis = "existing canonical ranked proteome-wide GO/GSEA terms filtered to SUS_vs_RES; GO-ID ontology themes use only is_a/part_of ancestry; color is the descriptive median NES across all mapped tested terms; original GO-term FDR support is emphasized separately",
      panel_c_representative_selection = "p.adjust < 0.05; lowest p.adjust; largest abs(NES); GO ID; Description",
      panel_c_recurrence_role = "descriptive_only_not_a_significance_or_display_gate",
      panel_c_NES_interpretation = "positive = higher in SUS; negative = higher in RES",
      combined_183mm_figure_readable = combined_readable
    ),
    notes = paste(
      "Downstream synthesis only. Differential abundance and ranked GO/GSEA statistics were not recomputed.",
      if (RENDER_ONLY) "Existing DAP-set ORA outputs were reused unchanged; clusterProfiler was not run." else "DAP-set GO-BP ORA was generated by this script.",
      "Panel C is ranked proteome-wide GO/GSEA evidence and is not calculated from the Panels A/B DAP sets. Manuscript themes are ontology-mapped summaries in the same ranked_GSEA evidence family. Median NES across mapped tested terms is descriptive only; representative term statistics remain the original individual GSEA term statistics; no theme-level p-value or FDR is constructed."
    )
  )

message("SUS-RES spatial DAP atlas complete: ", TABLE_DIR)
