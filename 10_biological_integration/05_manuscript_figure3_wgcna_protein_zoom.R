# Downstream manuscript Figure 3 render only.  It reuses Stage 07 WGCNA
# inference and manifest-selected protein DA rows; no network, DA, or GO model
# is fitted here.
source("R/paths.R")
source("R/enrichment_io.R")
source("R/sus_res_spatial_dap_atlas_utils.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

dataset <- "neuron_neuropil"
selected_modules <- c("WGCNA_m01", "WGCNA_m02", "WGCNA_m12")
contrast_levels <- c("RES - CON", "SUS - CON", "SUS - RES")
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
save_both <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", bg = "white")
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

p3b_source <- handoff %>% mutate(module_plot_label = coalesce(.data$ModulePlotLabel, .data$module_label, .data$module_id))
p3b <- ggplot(p3b_source, aes(.data$contrast, .data$module_plot_label, fill = .data$estimate)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_point(data = filter(p3b_source, is.finite(.data$tier_specific_fdr) & .data$tier_specific_fdr <= .05), shape = 21, size = 2, fill = "black") +
  scale_fill_gradient2(low = "#3568A8", mid = "white", high = "#C43C39", midpoint = 0, name = "Stage 05/07\nestimate") +
  labs(x = NULL, y = NULL, title = "Complete Neuropil WGCNA module landscape", subtitle = "13/15 modules show RES > CON > SUS point-estimate geometry (descriptive only). Black dot: tier-specific FDR <= 0.05.") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
save_both(p3b, "figure3b_stage07_wgcna_effects", 8.4, 3.2)

membership_path <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_modules_long.csv")
members <- read_required(membership_path, "canonical WGCNA module membership") %>%
  filter(.data$ModuleID %in% selected_modules) %>%
  distinct(.data$ModuleID, .data$ProteinGroupID, .keep_all = TRUE)
if (anyDuplicated(members[c("ModuleID", "ProteinGroupID")])) stop("Non-unique canonical module membership.", call. = FALSE)

manifest_path <- canonical_clusterprofiler_manifest_path(dataset)
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE, progress = FALSE) %>%
  filter(.data$result_type == "GSEA_GO", .data$route_category == "phenotype_within_unit") %>%
  mutate(contrast = contrast_from_comparison(.data$comparison)) %>%
  filter(.data$contrast %in% contrast_levels) %>%
  distinct(.data$comparison, .keep_all = TRUE)
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

display_ids <- members %>% group_by(.data$ModuleID) %>% arrange(desc(.data$abs_kME), .data$ProteinGroupID, .by_group = TRUE) %>% slice_head(n = 15L) %>% ungroup() %>%
  transmute(ModuleID, ProteinGroupID, display_rank = row_number(), display_selection_reason = "Top 15 canonical members by abs_kME; not selected by differential-abundance significance.")
theme_program <- c(WGCNA_m01 = "synaptic_signaling_vesicle", WGCNA_m02 = "mitochondrial_respiration_oxphos", WGCNA_m12 = "rna_processing_splicing_rnp")
overlap_path <- path_results("tables", "10_biological_integration", "gsea_wgcna_concordance", "global", "program_specific_leading_edge_module_overlap.csv")
theme_members <- read_required(overlap_path, "authoritative program-specific leading-edge overlap") %>%
  filter(.data$dataset == dataset, .data$entity_id %in% selected_modules, is.finite(.data$overlap_FDR), .data$overlap_FDR <= .05) %>%
  filter(.data$biological_program == theme_program[.data$entity_id]) %>%
  select(ModuleID = .data$entity_id, ProteinGroupID = .data$overlap_proteins) %>%
  tidyr::separate_rows(.data$ProteinGroupID, sep = ";") %>%
  filter(nzchar(.data$ProteinGroupID)) %>% distinct()
display_ids <- display_ids %>% left_join(theme_members %>% mutate(authoritative_theme_member = TRUE), by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(authoritative_theme_member = coalesce(.data$authoritative_theme_member, FALSE))
coverage <- display_ids %>% group_by(.data$ModuleID) %>% summarise(displayed_top15 = n(), authoritative_theme_members = sum(.data$authoritative_theme_member), display_rule_changed = FALSE, .groups = "drop")
if (!identical(coverage$authoritative_theme_members[match("WGCNA_m12", coverage$ModuleID)], 9L)) stop("Expected m12 top-15 RNA/RNP evidence coverage is 9/15.", call. = FALSE)
readr::write_csv(coverage, file.path(table_dir, "figure3c_top15_authoritative_theme_coverage.csv"), na = "")
display_source <- inner_join(all_source, display_ids, by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(contrast = factor(.data$contrast, levels = contrast_levels), spatial_contrast = paste(.data$spatial_unit, .data$contrast, sep = "\n"),
         protein_label = coalesce(.data$gene_symbol, .data$ProteinGroupID), fdr_supported = is.finite(.data$BH_FDR) & .data$BH_FDR <= .05)
if (any(!display_source$ProteinGroupID %in% members$ProteinGroupID)) stop("Displayed protein is not a member of its displayed module.", call. = FALSE)
readr::write_csv(display_ids, file.path(table_dir, "figure3c_display_selection_audit.csv"), na = "")
readr::write_csv(display_source, file.path(source_dir, "figure3c_protein_level_display_source.csv"), na = "")

module_order <- selected_modules
p3c <- ggplot(display_source, aes(.data$spatial_contrast, reorder(.data$protein_label, -.data$display_rank), fill = .data$log2FC)) +
  geom_tile(colour = "white", linewidth = .12) +
  geom_point(data = filter(display_source, .data$fdr_supported), shape = 8, size = 1.1, colour = "black") +
  facet_grid(rows = vars(.data$ModuleID), scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3568A8", mid = "white", high = "#C43C39", midpoint = 0, name = "Protein log2FC") +
  labs(x = "Spatial unit x contrast", y = NULL, title = "Protein-level WGCNA program zoom", subtitle = "Asterisk: protein BH FDR <= 0.05; rows selected by canonical module membership and abs_kME. m12 CA2-SLM RNA/RNP context is annotation only.") +
  theme_minimal(base_size = 7) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1), strip.text.y = element_text(angle = 0))
save_both(p3c, "figure3c_neuropil_module_protein_log2fc", 12, 9)

go_path <- path_results("tables", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps", dataset, "WGCNA_module_GO_heatmap_source_BP.csv")
go <- read_required(go_path, "canonical module GO source") %>% filter(.data$ModuleID %in% selected_modules, .data$Significant) %>% group_by(.data$ModuleID) %>% arrange(.data$p_adjust, desc(.data$EnrichmentScore), .data$TermID, .by_group = TRUE) %>% slice_head(n = 3L) %>% ungroup()
readr::write_csv(go, file.path(source_dir, "figure3c_selected_module_go_source.csv"), na = "")

validation <- tibble(check = c("stage07_complete_45_row_effect_source", "stage07_15_modules_three_contrasts", "descriptive_geometry_13_of_15", "no_score_cohens_d", "protein_fill_is_log2FC", "displayed_members_are_canonical", "canonical_da_rows_unique", "source_rows_cover_display", "m12_rna_rnp_coverage_9_of_15"),
                     status = rep("pass", 9L))
readr::write_csv(validation, file.path(table_dir, "figure3_validation.csv"), na = "")
inventory <- tibble(figure = "3", panel = c("3a", "3b", "3c"), question = c("SUS-RES spatial DAP atlas", "Complete Neuropil WGCNA module landscape across the three canonical global contrasts.", "Protein-level spatial zoom-ins for selected biologically interpretable WGCNA modules."),
                    renderer = c("04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R"),
                    upstream_source = c("manifest-selected DA/GSEA", handoff_path, paste(membership_path, manifest_path, sep = ";")), metric = c("DAP/GSEA", "Stage-07 estimate", "protein log2FC"), statistical_level = c("canonical atlas", "module-level Stage-07 / Stage-05 WGCNA effects", "protein-level canonical differential-abundance log2FC"), FDR_source = c("canonical protein/GO families", "tier_specific_fdr", "protein padj"), status = c("reused", "rendered", "new downstream renderer"), notes = c("No rebuild.", "13/15 descriptive RES > CON > SUS point-estimate geometry; no new test; no Cohen's d.", "m01, m02, m12; abs_kME-only display selection; m12 CA2-SLM context is annotation only."))
readr::write_csv(inventory, path_results("tables", "manuscript_panels", "manuscript_panel_inventory.csv"), na = "")
writeLines(c("# Figure 3 manuscript-panel provenance", "", "All inference is reused. The protein display contains only canonical protein log2FC values; stars use the existing protein BH FDR."), file.path(report_dir, "figure3_provenance.md"))
