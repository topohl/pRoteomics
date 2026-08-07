#!/usr/bin/env Rscript

# Manuscript-facing spatial atlas of SUS-RES differential protein changes.
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
source("R/sus_res_spatial_dap_atlas_utils.R")

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
  panel_c_source = file.path(SOURCE_DIR, "panel_c_sus_res_ranked_GSEA_programs.csv"),
  manifest = file.path(LOG_DIR, "run_manifest.yml")
)

RANKED_PROGRAM_SOURCE <- path_results(
  "source_data", MODULE_ID, "compareGO_spatial_atlas",
  "source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv"
)

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

required_packages <- c("ggplot2", "scales", "svglite", "cowplot")
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

panel_a <- ggplot2::ggplot(panel_a_source, ggplot2::aes(x = signed_count, y = plot_unit, fill = direction)) +
  ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.18) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "black") +
  ggplot2::facet_grid(dataset_label ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_y_discrete(labels = plot_unit_labels) +
  ggplot2::scale_x_continuous(labels = function(x) abs(x), expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
  ggplot2::scale_fill_manual(values = c("Higher in RES" = "#4E79A7", "Higher in SUS" = "#B04A7A"), drop = FALSE) +
  ggplot2::labs(
    x = "Differential ProteinGroupIDs (BH FDR <= 0.05)", y = NULL,
    title = "SUS-RES differential proteins by spatial unit", fill = NULL
  ) +
  theme_nature_base(base_size = 7, base_family = "sans") +
  ggplot2::theme(legend.position = "top", panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.2))

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

panel_b <- ggplot2::ggplot(panel_b_source, ggplot2::aes(x = plot_unit_x, y = plot_unit_y, fill = direction_aware_jaccard)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::facet_wrap(~dataset_label, scales = "free") +
  ggplot2::scale_x_discrete(labels = plot_unit_labels) +
  ggplot2::scale_y_discrete(labels = plot_unit_labels) +
  ggplot2::scale_fill_gradient(low = "#D9E2E8", high = "#1F5A83", limits = c(0, 1), na.value = "white") +
  ggplot2::labs(
    x = NULL, y = NULL, title = "Direction-aware DAP similarity", fill = "Jaccard",
    caption = "White = non-estimable: both DAP sets are empty in the pairwise-common tested universe."
  ) +
  theme_nature_heatmap(base_size = 6.5, base_family = "sans") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.5), axis.text.y = ggplot2::element_text(size = 5.5))

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
panel_c_source <- sus_res_prepare_ranked_program_panel(ranked_program_raw, expected_datasets = DATASETS)
utils::write.csv(panel_c_source, files$panel_c_source, row.names = FALSE)

include_ranked <- if (is.logical(panel_c_source$publication_include)) {
  panel_c_source$publication_include
} else {
  tolower(as.character(panel_c_source$publication_include)) %in% c("true", "1", "yes")
}
panel_c_plot_source <- panel_c_source[include_ranked, , drop = FALSE]
if (!nrow(panel_c_plot_source)) {
  stop("Ranked GO/GSEA Panel C source has no rows passing the existing publication inclusion rule.", call. = FALSE)
}
program_levels <- unique(as.character(panel_c_plot_source$program))
panel_c_plot_source$program <- factor(panel_c_plot_source$program, levels = rev(program_levels))
panel_c_plot_source$spatial_unit_label <- factor(
  panel_c_plot_source$spatial_unit_label,
  levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(panel_c_plot_source$spatial_unit)))
)
panel_c_plot_source$dataset_compartment_label <- factor(
  panel_c_plot_source$dataset_compartment_label,
  levels = dataset_levels
)
nes_limit <- suppressWarnings(stats::quantile(abs(panel_c_plot_source$mean_NES), probs = 0.98, na.rm = TRUE, names = FALSE))
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- suppressWarnings(max(abs(panel_c_plot_source$mean_NES), na.rm = TRUE))
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- 1
nes_limit <- min(nes_limit, 2.5)
support_values <- sort(unique(panel_c_plot_source$significant_term_count[panel_c_plot_source$significant_term_count > 0]))
support_breaks <- if (length(support_values)) unique(round(stats::quantile(support_values, probs = c(0, 0.5, 1), names = FALSE))) else c(1, 2, 3)

panel_c <- ggplot2::ggplot(
  panel_c_plot_source,
  ggplot2::aes(x = spatial_unit_label, y = program, colour = mean_NES, size = significant_term_count)
) +
  ggplot2::geom_point(alpha = 0.92) +
  ggplot2::facet_grid(dataset_compartment_label ~ ., scales = "free_x", space = "free_x") +
  ggplot2::scale_colour_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
    limits = c(-nes_limit, nes_limit), oob = scales::squish, name = "Mean NES"
  ) +
  ggplot2::scale_size_continuous(range = c(1.0, 4.2), breaks = support_breaks, name = "Sig. terms") +
  ggplot2::labs(
    x = NULL, y = NULL, title = "Ranked GO/GSEA program shifts",
    subtitle = "Positive NES = higher in SUS; negative NES = higher in RES",
    caption = "Proteome-wide ranked GO/GSEA evidence; Panel c is not calculated from Panels a/b FDR-supported DAP sets."
  ) +
  theme_nature_dotplot(base_size = 6.5, base_family = "sans") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.5),
    axis.text.y = ggplot2::element_text(size = 5.5),
    plot.caption = ggplot2::element_text(size = 5.2, hjust = 0)
  )

save_plot_triplet <- function(plot, stem, width_mm, height_mm) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = svglite::svglite, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = grDevices::pdf, useDingbats = FALSE, bg = "white", limitsize = FALSE)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = png_device, dpi = 300, bg = "white", limitsize = FALSE)
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
panel_c_terms <- if (nrow(go_display)) length(unique(paste(go_display$dataset, go_display$direction, go_display$ID, sep = "|"))) else 0L
panel_c_height <- max(80, min(210, 50 + 5 * panel_c_terms))
save_plot_triplet(panel_a, file.path(FIGURE_DIR, "panel_a_sus_res_dap_counts"), 105, 135)
save_plot_triplet(panel_b, file.path(FIGURE_DIR, "panel_b_sus_res_dap_jaccard"), 183, 120)
if (!RENDER_ONLY) save_plot_triplet(panel_c_ora, file.path(FIGURE_DIR, "panel_c_sus_res_dap_GO_BP"), 183, panel_c_height)
save_plot_triplet(panel_c, file.path(FIGURE_DIR, "panel_c_sus_res_ranked_GSEA_programs"), 183, 110)

combined_readable <- TRUE
panel_a_combined <- panel_a + ggplot2::labs(title = "DAP counts") + ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold"))
panel_b_combined <- panel_b + ggplot2::labs(title = "Direction-aware similarity") + ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold"))
panel_c_combined <- panel_c + ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold"))
top <- cowplot::plot_grid(panel_a_combined, panel_b_combined, nrow = 1, rel_widths = c(0.78, 1.42), labels = c("a", "b"), label_size = 10)
combined <- cowplot::plot_grid(top, panel_c_combined, ncol = 1, rel_heights = c(1.0, 1.15), labels = c("", "c"), label_size = 10)
combined_stem <- file.path(FIGURE_DIR, "figure_sus_res_spatial_dap_atlas_183mm")
save_plot_triplet(combined, combined_stem, 183, 225)
combined_files <- paste0(combined_stem, c(".svg", ".pdf", ".png"))

analysis_status$combined_183mm_figure_readable <- combined_readable
analysis_status$n_panel_c_union_terms <- length(unique(panel_c_plot_source$program))
analysis_status$panel_c_evidence_basis <- "ranked_proteome_wide_GO_GSEA"
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
      panel_c_basis = "existing canonical ranked proteome-wide GO/GSEA program summary filtered to SUS_vs_RES",
      panel_c_NES_interpretation = "positive = higher in SUS; negative = higher in RES",
      combined_183mm_figure_readable = combined_readable
    ),
    notes = paste(
      "Downstream synthesis only. Differential abundance and ranked GO/GSEA statistics were not recomputed.",
      if (RENDER_ONLY) "Existing DAP-set ORA outputs were reused unchanged; clusterProfiler was not run." else "DAP-set GO-BP ORA was generated by this script.",
      "Panel C is ranked proteome-wide GO/GSEA evidence and is not calculated from the Panels A/B DAP sets."
    )
  )

message("SUS-RES spatial DAP atlas complete: ", TABLE_DIR)
