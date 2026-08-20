# Downstream manuscript Figure 3 render only.  It reuses Stage 07 WGCNA
# inference and manifest-selected protein DA rows; no network, DA, or GO model
# is fitted here.
source("R/paths.R")
source("R/enrichment_io.R")
source("R/sus_res_spatial_dap_atlas_utils.R")
source("R/plotting_nature.R")
source("R/manuscript_figure3_utils.R")
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
handoff_hash <- file_hash_sha256(handoff_path)
handoff <- read_required(handoff_path, "complete Stage 07 module heatmap source") %>%
  filter(.data$contrast %in% contrast_levels) %>%
  mutate(contrast = factor(.data$contrast, levels = contrast_levels))
validation_records <- list()
stage07_row_count <- nrow(handoff)
stage07_module_count <- n_distinct(handoff$module_id)
stage07_counts <- count(handoff, .data$module_id)
stage07_contrasts <- split(as.character(handoff$contrast), handoff$module_id)
stage07_three_contrasts <-
  length(stage07_contrasts) == 15L &&
  all(stage07_counts$n == 3L) &&
  all(vapply(
    stage07_contrasts,
    function(x) identical(sort(unique(x)), sort(contrast_levels)),
    logical(1)
  ))
stage07_model_valid_rows <- sum(handoff$model_valid_for_inference %in% TRUE)
stage07_unique_rows <- !anyDuplicated(handoff[c("module_id", "contrast")])
stage07_no_cohens_d <- !any(grepl("Cohen", names(handoff), ignore.case = TRUE))
validation_records <- c(validation_records, list(
  figure3_validation_record(
    "stage07_row_count", "structural_contract", stage07_row_count, 45L,
    stage07_row_count == 45L, "stop", handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "stage07_module_count", "structural_contract", stage07_module_count, 15L,
    stage07_module_count == 15L, "stop", handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "stage07_three_contrasts_per_module", "structural_contract",
    stage07_three_contrasts, TRUE, stage07_three_contrasts, "stop",
    handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "stage07_model_valid_rows", "structural_contract", stage07_model_valid_rows,
    stage07_row_count, stage07_model_valid_rows == stage07_row_count, "stop",
    handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "stage07_module_contrast_rows_unique", "structural_contract",
    stage07_unique_rows, TRUE, stage07_unique_rows, "stop", handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "stage07_no_score_cohens_d_fields", "structural_contract",
    stage07_no_cohens_d, TRUE, stage07_no_cohens_d, "stop", handoff_path, handoff_hash
  )
))
figure3_assert_structural_checks(do.call(rbind, validation_records))
geometry <- handoff %>% select(module_id, contrast, estimate) %>% tidyr::pivot_wider(names_from = contrast, values_from = estimate) %>%
  mutate(descriptive_RES_gt_CON_gt_SUS = .data$`RES - CON` > 0 & .data$`SUS - CON` < 0 & .data$`SUS - RES` < 0,
         descriptive_geometry_note = "Descriptive point-estimate geometry only; no enrichment or new statistical test.")
geometry_count <- sum(geometry$descriptive_RES_gt_CON_gt_SUS %in% TRUE)
validation_records <- c(validation_records, list(figure3_validation_record(
  "descriptive_geometry_modules", "frozen_manuscript_result_regression",
  geometry_count, 13L, geometry_count == 13L, "report_only", handoff_path,
  handoff_hash,
  "Snapshot-linked descriptive point-estimate geometry; not an analytical validity condition."
)))
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
membership_hash <- file_hash_sha256(membership_path)
members_raw <- read_required(membership_path, "canonical WGCNA module membership") %>%
  filter(.data$ModuleID %in% selected_modules)
canonical_membership_unique <-
  !anyDuplicated(members_raw[c("ModuleID", "ProteinGroupID")])
validation_records <- c(validation_records, list(figure3_validation_record(
  "canonical_module_membership_rows_unique", "structural_contract",
  canonical_membership_unique, TRUE, canonical_membership_unique, "stop",
  membership_path, membership_hash
)))
figure3_assert_structural_checks(do.call(rbind, validation_records))
members <- distinct(members_raw, ModuleID, ProteinGroupID, .keep_all = TRUE)

manifest_path <- canonical_clusterprofiler_manifest_path(dataset)
manifest_hash <- file_hash_sha256(manifest_path)
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
                  ProteinGroupID = .data$ProteinGroupID,
                  gene_symbol = figure3_display_gene_symbol(z),
                  log2FC = .data$log2fc, raw_p = .data$pval, BH_FDR = .data$padj,
                  source_artifact = input, manifest_input_hash = row$input_hash[[1]])
}))
canonical_da_rows_unique <- !anyDuplicated(da[c("spatial_unit", "contrast", "ProteinGroupID")])
validation_records <- c(validation_records, list(figure3_validation_record(
  "canonical_da_rows_unique", "structural_contract", canonical_da_rows_unique,
  TRUE, canonical_da_rows_unique, "stop", manifest_path, manifest_hash
)))
figure3_assert_structural_checks(do.call(rbind, validation_records))

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
overlap_hash <- file_hash_sha256(overlap_path)
theme_members <- read_required(overlap_path, "authoritative program-specific leading-edge overlap") %>%
  filter(.data$dataset == dataset, .data$entity_id %in% selected_modules, is.finite(.data$overlap_FDR), .data$overlap_FDR <= .05) %>%
  filter(.data$biological_program == theme_program[.data$entity_id]) %>%
  select(entity_id, overlap_proteins) %>%
  rename(ModuleID = entity_id, ProteinGroupID = overlap_proteins) %>%
  tidyr::separate_rows(ProteinGroupID, sep = ";") %>%
  filter(nzchar(.data$ProteinGroupID)) %>% distinct()
display_ids <- display_ids %>% left_join(theme_members %>% mutate(authoritative_theme_member = TRUE), by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(authoritative_theme_member = coalesce(.data$authoritative_theme_member, FALSE))
display_counts <- count(display_ids, .data$ModuleID)
exactly_15_selected <-
  setequal(display_counts$ModuleID, selected_modules) && all(display_counts$n == 15L)
display_rank_exact <-
  length(split(display_ids$display_rank, display_ids$ModuleID)) == length(selected_modules) &&
  all(vapply(
    split(display_ids$display_rank, display_ids$ModuleID),
    function(x) identical(sort(as.integer(x)), 1:15), logical(1)
  ))
validation_records <- c(validation_records, list(
  figure3_validation_record(
    "exactly_15_selected_proteins_per_module", "structural_contract",
    paste(display_counts$ModuleID, display_counts$n, sep = "="),
    paste(selected_modules, 15L, sep = "="), exactly_15_selected, "stop",
    membership_path, membership_hash
  ),
  figure3_validation_record(
    "display_rank_is_1_to_15_within_module", "structural_contract",
    display_rank_exact, TRUE, display_rank_exact, "stop",
    membership_path, membership_hash
  )
))
figure3_assert_structural_checks(do.call(rbind, validation_records))
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
m12_coverage <- coverage$authoritative_theme_members[match("WGCNA_m12", coverage$ModuleID)]
validation_records <- c(validation_records, list(figure3_validation_record(
  "m12_authoritative_rna_rnp_overlap_coverage",
  "frozen_manuscript_result_regression", m12_coverage, 9L,
  identical(as.integer(m12_coverage), 9L), "report_only",
  paste(membership_path, overlap_path, sep = ";"),
  paste(membership_hash, overlap_hash, sep = ";"),
  "Snapshot-linked top-15 overlap coverage; not an analytical validity condition."
)))
readr::write_csv(coverage, file.path(table_dir, "figure3c_top15_authoritative_theme_coverage.csv"), na = "")
display_source <- inner_join(all_source, display_ids, by = c("ModuleID", "ProteinGroupID")) %>%
  mutate(
    contrast = factor(.data$contrast, levels = contrast_levels),
    spatial_contrast = paste(.data$spatial_unit, .data$contrast, sep = "\n"),
    protein_label = coalesce(.data$gene_symbol, .data$ProteinGroupID),
    fdr_supported = is.finite(.data$BH_FDR) & .data$BH_FDR <= .05
  )
member_keys <- paste(members$ModuleID, members$ProteinGroupID, sep = "\r")
display_keys <- paste(display_ids$ModuleID, display_ids$ProteinGroupID, sep = "\r")
displayed_members_are_canonical <- all(display_keys %in% member_keys)
source_display_keys <- unique(paste(
  display_source$ModuleID, display_source$ProteinGroupID, sep = "\r"
))
source_rows_cover_display <- setequal(unique(display_keys), source_display_keys)
validation_records <- c(validation_records, list(
  figure3_validation_record(
    "displayed_proteins_are_canonical_module_members", "structural_contract",
    displayed_members_are_canonical, TRUE, displayed_members_are_canonical, "stop",
    membership_path, membership_hash
  ),
  figure3_validation_record(
    "source_rows_cover_all_displayed_proteins", "structural_contract",
    length(source_display_keys), length(unique(display_keys)),
    source_rows_cover_display, "stop", manifest_path, manifest_hash
  )
))
figure3_assert_structural_checks(do.call(rbind, validation_records))
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

displayed_protein_keys <- unique(paste(
  display_source$ModuleID, display_source$ProteinGroupID, sep = "\r"
))
supported_display_protein_keys <- unique(paste(
  display_source$ModuleID[display_source$fdr_supported],
  display_source$ProteinGroupID[display_source$fdr_supported], sep = "\r"
))
displayed_proteins_with_bh_support <- length(supported_display_protein_keys)
displayed_da_rows_with_bh_support <- sum(display_source$fdr_supported %in% TRUE)
stage07_effects_with_tier_fdr_support <- nrow(p3b_supported)
validation_records <- c(validation_records, list(
  figure3_validation_record(
    "displayed_proteins_with_BH_FDR_le_0_05", "observed_manuscript_result",
    displayed_proteins_with_bh_support, "observed_value_reported",
    NA, "report_only", manifest_path, manifest_hash
  ),
  figure3_validation_record(
    "displayed_DA_rows_with_BH_FDR_le_0_05", "observed_manuscript_result",
    displayed_da_rows_with_bh_support, "observed_value_reported",
    NA, "report_only", manifest_path, manifest_hash
  ),
  figure3_validation_record(
    "figure3b_module_effects_with_tier_specific_FDR_le_0_05",
    "observed_manuscript_result", stage07_effects_with_tier_fdr_support,
    "observed_value_reported", NA, "report_only", handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "descriptive_RES_gt_CON_gt_SUS_geometry_count", "observed_manuscript_result",
    geometry_count, "observed_value_reported", NA, "report_only",
    handoff_path, handoff_hash
  ),
  figure3_validation_record(
    "current_m12_authoritative_RNA_RNP_overlap_coverage",
    "observed_manuscript_result", paste0(m12_coverage, "/15"),
    "observed_value_reported", NA, "report_only",
    paste(membership_path, overlap_path, sep = ";"),
    paste(membership_hash, overlap_hash, sep = ";")
  )
))
validation <- do.call(rbind, validation_records)
figure3_assert_structural_checks(validation)
figure3_warn_snapshot_mismatches(validation)
readr::write_csv(validation, file.path(table_dir, "figure3_validation.csv"), na = "")
focused_go_path <- path_results("tables", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps", dataset, "WGCNA_supermodule_GO_focused_source_BP.csv")
inventory <- tibble(figure = "3", panel = c("3a", "3b", "3c", "3d"), question = c("SUS-RES spatial DAP atlas", "Complete Neuropil WGCNA module landscape across the three canonical global contrasts.", "Focused Neuropil supermodule GO-BP member-module evidence.", "Protein-level spatial zoom-ins for selected biologically interpretable WGCNA modules."),
                    renderer = c("04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R", "06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R", "10_biological_integration/05_manuscript_figure3_wgcna_protein_zoom.R"),
                    upstream_source = c("manifest-selected DA/GSEA", handoff_path, focused_go_path, paste(membership_path, manifest_path, sep = ";")), metric = c("DAP/GSEA", "Stage-07 estimate", "member-module GO support", "protein log2FC"), statistical_level = c("canonical atlas", "module-level Stage-07 / Stage-05 WGCNA effects", "member-module evidence; no pooled supermodule inference", "protein-level canonical differential-abundance log2FC"), FDR_source = c("canonical protein/GO families", "tier_specific_fdr", "member-module BH FDR", "protein padj"), status = c("reused", "rendered", "reused", "downstream renderer"), notes = c("No rebuild.", sprintf("%d/%d descriptive RES > CON > SUS point-estimate geometry; no new test; no Cohen's d.", geometry_count, stage07_module_count), "Selected GO terms, recurrence and redundancy pruning unchanged.", sprintf("m01, m02, m12; abs(kME)-only display selection; within-module display ranks; current m12 authoritative RNA/RNP overlap coverage %d/15; CA2-SLM context is annotation only.", m12_coverage)))
readr::write_csv(inventory, path_results("tables", "manuscript_panels", "manuscript_panel_inventory.csv"), na = "")
protein_support_text <- if (displayed_proteins_with_bh_support > 0L) {
  sprintf(
    "%d of %d displayed canonical module proteins have at least one displayed DA row with BH FDR <= 0.05 (%d supported rows); those rows carry the plotted support symbol.",
    displayed_proteins_with_bh_support, length(displayed_protein_keys),
    displayed_da_rows_with_bh_support
  )
} else {
  sprintf(
    "0 of %d displayed canonical module proteins has a displayed DA row with BH FDR <= 0.05; no protein-level support symbols are present.",
    length(displayed_protein_keys)
  )
}
module_support_text <- if (stage07_effects_with_tier_fdr_support > 0L) {
  sprintf(
    "%d of %d Figure 3b module-effect rows have tier-specific FDR <= 0.05 and carry the plotted support marker.",
    stage07_effects_with_tier_fdr_support, stage07_row_count
  )
} else {
  sprintf(
    "0 of %d Figure 3b module-effect rows has tier-specific FDR <= 0.05; no module-level support markers are present.",
    stage07_row_count
  )
}
writeLines(c(
  "# Figure 3 manuscript-panel provenance", "",
  "All inference is reused. The protein display contains only canonical protein log2FC values.",
  protein_support_text,
  module_support_text,
  sprintf(
    "%d/%d modules show the descriptive RES > CON > SUS point-estimate geometry; this is descriptive only and is not a new statistical test.",
    geometry_count, stage07_module_count
  ),
  sprintf(
    "Current WGCNA_m12 authoritative RNA/RNP overlap coverage is %d/15 displayed proteins.",
    m12_coverage
  )
), file.path(report_dir, "figure3_provenance.md"))
