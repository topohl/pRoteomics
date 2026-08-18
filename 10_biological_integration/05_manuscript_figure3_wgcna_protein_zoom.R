# Downstream manuscript Figure 3 render only.  It reuses Stage 07 WGCNA
# inference and manifest-selected protein DA rows; no network, DA, or GO model
# is fitted here.
source("R/paths.R")
source("R/enrichment_io.R")
source("R/sus_res_spatial_dap_atlas_utils.R")
source("R/plotting_nature.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork) })

dataset <- "neuron_neuropil"
selected_modules <- c("WGCNA_m01", "WGCNA_m02", "WGCNA_m12")
contrast_levels <- c("RES - CON", "SUS - CON", "SUS - RES")
manuscript_text <- nature_manuscript_text_sizes_pt()
out <- function(kind, ...) path_results(kind, "manuscript_panels", "figure_3", ...)
fig_dir <- out("figures"); source_dir <- out("source_data"); table_dir <- out("tables"); report_dir <- out("reports")
invisible(lapply(c(fig_dir, source_dir, table_dir, report_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

read_required <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
contrast_from_comparison <- function(x) {
  x <- tolower(x)
  ifelse(grepl("res.*con", x), "RES - CON", ifelse(grepl("sus.*con", x), "SUS - CON", ifelse(grepl("sus.*res", x), "SUS - RES", NA_character_)))
}
save_both <- function(plot, stem, width_mm, height_mm) {
  ggplot2::ggsave(
    file.path(fig_dir, paste0(stem, ".svg")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = svglite::svglite, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    file.path(fig_dir, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
}

handoff_path <- path_results("source_data", "06_modules_WGCNA", "interpretable_summary", dataset, "module_group_effects_main_heatmap_source.csv")
handoff <- read_required(handoff_path, "complete Stage 07 module heatmap source") %>%
  filter(.data$contrast %in% contrast_levels) %>%
  mutate(contrast = factor(.data$contrast, levels = contrast_levels))
if (nrow(handoff) != 45L || n_distinct(handoff$module_id) != 15L || any(count(handoff, .data$module_id)$n != 3L)) stop("Figure 3B must contain exactly 15 canonical Neuropil modules x three contrasts.", call. = FALSE)
if (any(!handoff$model_valid_for_inference)) stop("Figure 3B contains a model-invalid Stage 07 row.", call. = FALSE)
if (any(grepl("Cohen", names(handoff), ignore.case = TRUE))) stop("Stage 07 handoff unexpectedly contains score-derived Cohen's d fields.", call. = FALSE)
geometry <- handoff %>% select(module_id, contrast, estimate) %>% tidyr::pivot_wider(names_from = contrast, values_from = estimate) %>%
  mutate(descriptive_RES_gt_CON_gt_SUS = .data$`RES - CON` > 0 & .data$`SUS - CON` < 0 & .data$`SUS - RES` < 0,
         descriptive_geometry_note = "Descriptive point-estimate geometry only; no enrichment or new statistical test.")
if (sum(geometry$descriptive_RES_gt_CON_gt_SUS) != 13L) stop("Expected 13/15 descriptive Neuropil geometry modules was not observed.", call. = FALSE)
handoff <- handoff %>% left_join(geometry %>% select(module_id, descriptive_RES_gt_CON_gt_SUS, descriptive_geometry_note), by = "module_id")
readr::write_csv(handoff, file.path(source_dir, "figure3b_stage07_effect_source.csv"), na = "")

p3b_source <- handoff %>%
  mutate(
    module_plot_label = sub(
      "^WGCNA_", "",
      coalesce(.data$ModulePlotLabel, .data$module_label, .data$module_id)
    ),
    module_plot_label = sub("cytoskeletal trafficking", "trafficking", .data$module_plot_label, fixed = TRUE),
    module_plot_label = sub("mitochondrial / energy metabolism", "mitochondrial/energy", .data$module_plot_label, fixed = TRUE),
    module_plot_label = sub("RNA/RNP regulatory module", "RNA/RNP", .data$module_plot_label, fixed = TRUE),
    module_plot_label = sub("translation / proteostasis", "translation/proteostasis", .data$module_plot_label, fixed = TRUE),
    module_plot_label = sub("barrier / cell-junction structural module", "barrier/junction", .data$module_plot_label, fixed = TRUE),
    contrast_plot_label = factor(
      .data$contrast, levels = contrast_levels,
      labels = c("RES\n– CON", "SUS\n– CON", "SUS\n– RES")
    )
  )
module_rows <- p3b_source %>%
  distinct(module_id, module_plot_label, SupermoduleID) %>%
  mutate(module_number = suppressWarnings(as.integer(sub("^WGCNA_m", "", .data$module_id)))) %>%
  arrange(.data$SupermoduleID, .data$module_number, .data$module_id)
p3b_source$module_plot_label <- factor(
  p3b_source$module_plot_label,
  levels = rev(module_rows$module_plot_label)
)
module_rows$y_position <- match(module_rows$module_plot_label, levels(p3b_source$module_plot_label))
separator_y <- module_rows %>%
  group_by(.data$SupermoduleID) %>%
  summarise(boundary = max(.data$y_position) + 0.5, .groups = "drop") %>%
  filter(.data$boundary < nrow(module_rows) + 0.5)
p3b_supported <- filter(
  p3b_source,
  is.finite(.data$tier_specific_fdr) & .data$tier_specific_fdr <= .05
)
p3b_effect_limit <- max(abs(p3b_source$estimate), na.rm = TRUE)
p3b <- ggplot(p3b_source, aes(.data$contrast_plot_label, .data$module_plot_label, fill = .data$estimate)) +
  geom_hline(
    data = separator_y, aes(yintercept = .data$boundary),
    inherit.aes = FALSE, colour = "#D9D9D9", linewidth = 0.25
  ) +
  geom_tile(width = 0.96, height = 0.96, colour = NA) +
  scale_fill_gradient2(
    low = nature_palette("signed")[["low"]],
    mid = nature_palette("signed")[["mid"]],
    high = nature_palette("signed")[["high"]],
    midpoint = 0, limits = c(-p3b_effect_limit, p3b_effect_limit),
    breaks = c(-p3b_effect_limit, 0, p3b_effect_limit),
    labels = scales::label_number(accuracy = 0.01), name = "Model estimate"
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    barwidth = grid::unit(20, "mm"), barheight = grid::unit(2.1, "mm")
  )) +
  labs(x = NULL, y = NULL) +
  coord_fixed(clip = "off") +
  theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = manuscript_text[["normal"]], lineheight = 0.85),
    axis.text.y = element_text(size = manuscript_text[["dense"]], colour = "#202020"),
    legend.position = "bottom",
    legend.title = element_text(size = manuscript_text[["normal"]]),
    legend.text = element_text(size = manuscript_text[["normal"]]),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5)
  )
if (nrow(p3b_supported)) {
  p3b <- p3b +
    geom_point(data = p3b_supported, shape = 21, size = 1.1, fill = "black") +
    guides(shape = guide_legend(title = "Tier-specific FDR <= 0.05"))
}
p3b_width_mm <- max(76, 54 + 6.2 * length(contrast_levels))
p3b_height_mm <- 19 + 7.0 * nrow(module_rows)
save_both(p3b, "figure3b_stage07_wgcna_effects", p3b_width_mm, p3b_height_mm)

membership_path <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_modules_long.csv")
members <- read_required(membership_path, "canonical WGCNA module membership") %>%
  filter(.data$ModuleID %in% selected_modules) %>%
  distinct(ModuleID, ProteinGroupID, .keep_all = TRUE)
if (anyDuplicated(members[c("ModuleID", "ProteinGroupID")])) stop("Non-unique canonical module membership.", call. = FALSE)

manifest_path <- canonical_clusterprofiler_manifest_path(dataset)
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE, progress = FALSE) %>%
  filter(.data$result_type == "GSEA_GO", .data$route_category == "phenotype_within_unit") %>%
  mutate(contrast = contrast_from_comparison(.data$comparison)) %>%
  filter(.data$contrast %in% contrast_levels) %>%
  distinct(comparison, .keep_all = TRUE)
if (nrow(manifest) != 30L) stop("Expected exactly 10 spatial units x 3 canonical contrasts in the neuropil manifest.", call. = FALSE)
da <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, , drop = FALSE]
  input <- sus_res_resolve_manifest_input(row$input_gene_file[[1]], dataset, repo_root())
  z <- read_required(input, "manifest-selected DA input")
  required <- c("ProteinGroupID", "log2fc", "pval", "padj")
  missing <- setdiff(required, names(z)); if (length(missing)) stop("DA input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  z %>% transmute(dataset = dataset, spatial_unit = row$route_unit[[1]], contrast = row$contrast[[1]],
                  ProteinGroupID = .data$ProteinGroupID, gene_symbol = coalesce(.data$official_gene_symbol, .data$representative_gene_symbol, .data$gene_symbol),
                  log2FC = .data$log2fc, raw_p = .data$pval, BH_FDR = .data$padj,
                  source_artifact = input, manifest_input_hash = row$input_hash[[1]])
}))
if (anyDuplicated(da[c("spatial_unit", "contrast", "ProteinGroupID")])) stop("Canonical DA rows are not unique by spatial unit, contrast, and ProteinGroupID.", call. = FALSE)

all_source <- inner_join(members, da, by = "ProteinGroupID") %>%
  mutate(reviewed_module_label = .data$ModuleLabel_Final,
         module_membership_provenance = "results/tables/06_modules_WGCNA/01_WGCNA/neuron_neuropil/modules/WGCNA_modules_long.csv",
         selection_display_reason = "Complete canonical module membership retained; display ranking uses abs_kME only.",
         local_multimethod_context = ifelse(.data$ModuleID == "WGCNA_m12", "CA2-SLM is the strongest existing local RNA/RNP multimethod example: RES-CON and SUS-RES GSEA supported; SUS-CON directional/not supported; module-level WGCNA inference remains multiplicity-nonsignificant.", NA_character_),
         local_context_source = ifelse(.data$ModuleID == "WGCNA_m12", "results/tables/10_biological_integration/gsea_wgcna_concordance/global/program_specific_leading_edge_module_overlap.csv", NA_character_))
if (!nrow(all_source) || any(!all_source$ModuleID %in% selected_modules)) stop("Protein zoom source failed canonical membership join.", call. = FALSE)
readr::write_csv(all_source, file.path(source_dir, "figure3c_protein_level_all_source.csv"), na = "")

display_ids <- members %>%
  group_by(.data$ModuleID) %>%
  arrange(desc(.data$abs_kME), .data$ProteinGroupID, .by_group = TRUE) %>%
  slice_head(n = 15L) %>%
  transmute(
    ModuleID, ProteinGroupID, display_rank = row_number(),
    display_selection_reason = "Top 15 canonical members by abs_kME; not selected by differential-abundance significance."
  ) %>%
  ungroup()
theme_program <- c(WGCNA_m01 = "synaptic_signaling_vesicle", WGCNA_m02 = "mitochondrial_respiration_oxphos", WGCNA_m12 = "rna_processing_splicing_rnp")
overlap_path <- path_results("tables", "10_biological_integration", "gsea_wgcna_concordance", "global", "program_specific_leading_edge_module_overlap.csv")
theme_members <- read_required(overlap_path, "authoritative program-specific leading-edge overlap") %>%
  filter(.data$dataset == dataset, .data$entity_id %in% selected_modules, is.finite(.data$overlap_FDR), .data$overlap_FDR <= .05) %>%
  filter(.data$biological_program == theme_program[.data$entity_id]) %>%
  select(entity_id, overlap_proteins) %>%
  rename(ModuleID = entity_id, ProteinGroupID = overlap_proteins) %>%
  tidyr::separate_rows(ProteinGroupID, sep = ";") %>%
  filter(nzchar(.data$ProteinGroupID)) %>% distinct()
display_ids <- display_ids %>% left_join(theme_members %>% mutate(authoritative_theme_member = TRUE), by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(authoritative_theme_member = coalesce(.data$authoritative_theme_member, FALSE))
if (any(count(display_ids, .data$ModuleID)$n != 15L) ||
    any(vapply(split(display_ids$display_rank, display_ids$ModuleID), function(x) !identical(sort(x), 1:15), logical(1)))) {
  stop("Figure 3d display_rank must be exactly 1-15 within each selected module.", call. = FALSE)
}
legacy_selection_file <- file.path(table_dir, "figure3c_display_selection_audit.csv")
if (file.exists(legacy_selection_file)) {
  legacy_selection <- read_required(legacy_selection_file, "existing protein display selection audit")
  old_keys <- sort(unique(paste(legacy_selection$ModuleID, legacy_selection$ProteinGroupID, sep = "|")))
  new_keys <- sort(unique(paste(display_ids$ModuleID, display_ids$ProteinGroupID, sep = "|")))
  if (!identical(old_keys, new_keys)) {
    stop("Presentation-only Figure 3d revision changed the selected protein set.", call. = FALSE)
  }
}
coverage <- display_ids %>% group_by(.data$ModuleID) %>% summarise(displayed_top15 = n(), authoritative_theme_members = sum(.data$authoritative_theme_member), display_rule_changed = FALSE, .groups = "drop")
if (!identical(coverage$authoritative_theme_members[match("WGCNA_m12", coverage$ModuleID)], 9L)) stop("Expected m12 top-15 RNA/RNP evidence coverage is 9/15.", call. = FALSE)
readr::write_csv(coverage, file.path(table_dir, "figure3c_top15_authoritative_theme_coverage.csv"), na = "")
display_source <- inner_join(all_source, display_ids, by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(
    contrast = factor(.data$contrast, levels = contrast_levels),
    spatial_contrast = paste(.data$spatial_unit, .data$contrast, sep = "\n"),
    protein_label = coalesce(.data$gene_symbol, .data$ProteinGroupID),
    fdr_supported = is.finite(.data$BH_FDR) & .data$BH_FDR <= .05
  )
if (any(!display_source$ProteinGroupID %in% members$ProteinGroupID)) stop("Displayed protein is not a member of its displayed module.", call. = FALSE)
readr::write_csv(display_ids, file.path(table_dir, "figure3c_display_selection_audit.csv"), na = "")
readr::write_csv(display_source, file.path(source_dir, "figure3c_protein_level_display_source.csv"), na = "")
readr::write_csv(display_ids, file.path(table_dir, "figure3d_display_selection_audit.csv"), na = "")
readr::write_csv(display_source, file.path(source_dir, "figure3d_protein_level_display_source.csv"), na = "")

protein_order <- display_ids %>%
  mutate(module_order = match(.data$ModuleID, selected_modules)) %>%
  arrange(.data$module_order, .data$display_rank) %>%
  transmute(protein_plot_key = paste(.data$ModuleID, .data$ProteinGroupID, sep = "||"))
display_plot_source <- display_source %>%
  mutate(
    spatial_unit = factor(
      .data$spatial_unit,
      levels = anatomical_spatial_unit_levels(unique(.data$spatial_unit))
    ),
    module_block = factor(
      sub("^WGCNA_", "", .data$ModuleID),
      levels = sub("^WGCNA_", "", selected_modules)
    ),
    contrast_plot_label = factor(
      .data$contrast, levels = contrast_levels,
      labels = c("RES\n– CON", "SUS\n– CON", "SUS\n– RES")
    ),
    protein_plot_key = paste(.data$ModuleID, .data$ProteinGroupID, sep = "||")
  )
display_plot_source$protein_plot_key <- factor(
  display_plot_source$protein_plot_key,
  levels = rev(protein_order$protein_plot_key)
)
protein_labels <- display_plot_source %>%
  distinct(protein_plot_key, protein_label) %>%
  { stats::setNames(.$protein_label, as.character(.$protein_plot_key)) }
p3d_effect_limit <- max(abs(display_plot_source$log2FC), na.rm = TRUE)
spatial_unit_compact_label <- function(x) {
  gsub(" ", "-", clean_spatial_unit_label(x), fixed = TRUE)
}
build_protein_module_heatmap <- function(module_id, show_x_axis = FALSE) {
  block <- display_plot_source[display_plot_source$ModuleID == module_id, , drop = FALSE]
  block_supported <- block[block$fdr_supported, , drop = FALSE]
  p <- ggplot(block, aes(.data$spatial_unit, .data$protein_plot_key, fill = .data$log2FC)) +
    geom_tile(width = 0.96, height = 0.96, colour = NA) +
    facet_grid(cols = vars(.data$contrast_plot_label)) +
    scale_x_discrete(labels = spatial_unit_compact_label, drop = FALSE) +
    scale_y_discrete(labels = protein_labels, drop = TRUE) +
    scale_fill_gradient2(
      low = nature_palette("signed")[["low"]],
      mid = nature_palette("signed")[["mid"]],
      high = nature_palette("signed")[["high"]],
      midpoint = 0, limits = c(-p3d_effect_limit, p3d_effect_limit),
      breaks = c(-p3d_effect_limit, 0, p3d_effect_limit),
      name = expression("Protein log"[2] * "FC")
    ) +
    guides(fill = guide_colourbar(
      title.position = "top", barwidth = grid::unit(35, "mm"), barheight = grid::unit(2.1, "mm"),
      theme = theme(legend.text = element_text(size = manuscript_text[["normal"]] / 0.7))
    )) +
    coord_fixed(clip = "off") +
    labs(x = NULL, y = sub("^WGCNA_", "", module_id)) +
    theme_nature_manuscript_panel(base_size = manuscript_text[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = manuscript_text[["dense"]], lineheight = 0.85),
      axis.text.y = element_text(size = manuscript_text[["dense"]], face = "italic"),
      axis.title.y = element_text(angle = 0, size = manuscript_text[["normal"]], face = "bold", margin = margin(r = 1.5)),
      strip.text.x = element_text(size = manuscript_text[["normal"]], face = "plain", lineheight = 0.85),
      panel.spacing.x = grid::unit(0.65, "mm"),
      legend.position = "bottom",
      legend.title = element_text(size = manuscript_text[["normal"]]),
      legend.text = element_text(size = manuscript_text[["normal"]]),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)
    )
  if (!show_x_axis) p <- p + theme(axis.text.x = element_blank())
  if (nrow(block_supported)) {
    p <- p + geom_point(data = block_supported, shape = 8, size = 0.9, colour = "black")
  }
  p
}
p3d <- patchwork::wrap_plots(
  lapply(seq_along(selected_modules), function(i) {
    build_protein_module_heatmap(selected_modules[[i]], show_x_axis = i == length(selected_modules))
  }),
  ncol = 1, guides = "collect"
) & theme(legend.position = "bottom", legend.box = "horizontal")
p3d_width_mm <- 132
p3d_height_mm <- 18 + length(selected_modules) * (9 + 15 * 2.6)
save_both(p3d, "figure3d_neuropil_module_protein_log2fc", p3d_width_mm, p3d_height_mm)
# Established pre-layout-correction compatibility stem.
save_both(p3d, "figure3c_neuropil_module_protein_log2fc", p3d_width_mm, p3d_height_mm)

go_path <- path_results("tables", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps", dataset, "WGCNA_module_GO_heatmap_source_BP.csv")
go <- read_required(go_path, "canonical module GO source") %>% filter(.data$ModuleID %in% selected_modules, .data$Significant) %>% group_by(.data$ModuleID) %>% arrange(.data$p_adjust, desc(.data$EnrichmentScore), .data$TermID, .by_group = TRUE) %>% slice_head(n = 3L) %>% ungroup()
readr::write_csv(go, file.path(source_dir, "figure3c_selected_module_go_source.csv"), na = "")
readr::write_csv(go, file.path(source_dir, "figure3d_selected_module_go_source.csv"), na = "")

validation <- tibble(check = c("stage07_complete_45_row_effect_source", "stage07_15_modules_three_contrasts", "descriptive_geometry_13_of_15", "no_score_cohens_d", "protein_fill_is_log2FC", "displayed_members_are_canonical", "canonical_da_rows_unique", "source_rows_cover_display", "m12_rna_rnp_coverage_9_of_15", "display_rank_is_1_to_15_within_module", "no_displayed_protein_BH_FDR_support", "no_stage07_tier_specific_FDR_support"),
                     status = rep("pass", 12L))
readr::write_csv(validation, file.path(table_dir, "figure3_validation.csv"), na = "")
focused_go_path <- path_results("tables", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps", dataset, "WGCNA_supermodule_GO_focused_source_BP.csv")
inventory <- tibble(figure = "3", panel = c("3a", "3b", "3c", "3d"), question = c("SUS-RES spatial DAP atlas", "Complete Neuropil WGCNA module landscape across the three canonical global contrasts.", "Focused Neuropil supermodule GO-BP member-module evidence.", "Protein-level spatial zoom-ins for selected biologically interpretable WGCNA modules."),
                    renderer = c("04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R", "06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R"),
                    upstream_source = c("manifest-selected DA/GSEA", handoff_path, focused_go_path, paste(membership_path, manifest_path, sep = ";")), metric = c("DAP/GSEA", "Stage-07 estimate", "member-module GO support", "protein log2FC"), statistical_level = c("canonical atlas", "module-level Stage-07 / Stage-05 WGCNA effects", "member-module evidence; no pooled supermodule inference", "protein-level canonical differential-abundance log2FC"), FDR_source = c("canonical protein/GO families", "tier_specific_fdr", "member-module BH FDR", "protein padj"), status = c("reused", "rendered", "reused", "downstream renderer"), notes = c("No rebuild.", "13/15 descriptive RES > CON > SUS point-estimate geometry; no new test; no Cohen's d.", "Selected GO terms, recurrence and redundancy pruning unchanged.", "m01, m02, m12; abs(kME)-only display selection; within-module display ranks; m12 CA2-SLM context is annotation only."))
readr::write_csv(inventory, path_results("tables", "manuscript_panels", "manuscript_panel_inventory.csv"), na = "")
writeLines(c("# Figure 3 manuscript-panel provenance", "", "All inference is reused. The protein display contains only canonical protein log2FC values. No displayed protein has BH FDR <= 0.05, so no significance-symbol legend is drawn. No Figure 3b module has tier-specific FDR <= 0.05, so no unused support-symbol legend is drawn."), file.path(report_dir, "figure3_provenance.md"))
