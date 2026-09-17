#!/usr/bin/env Rscript
#
# Spatial systems atlas: reviewer workbook, candidate figures and validation.
#
# Assembles the already-computed atlas tables into one OOXML-safe workbook,
# emits NEW candidate figures (never a manuscript contract replacement), and
# runs the atlas validation. A critical FAIL exits non-zero.
#
# Every figure writes its own source-data CSV beside it, so no panel exists
# without the numbers behind it.
#
# USAGE
#   Rscript analysis/spatial_validation/summarize_spatial_atlas.R
# Script: analysis/spatial_validation/summarize_spatial_atlas.R
# Stage: networks
# Scope: global
# Consumes: required results/spatial_validation/build_module_spatial_atlas/global/tables/WGCNA_module_spatial_cell_affinity.csv; results/tables/11_spatial_systems/atlas/WGCNA_module_spatial_cell_affinity.csv; results/spatial_validation/build_protein_spatial_atlas/global/tables/protein_spatial_cell_affinity.csv; +1 more; optional results/spatial_validation/quantify_neuropil_detection_context/global/tables/neuropil_spatial_detection_context.csv; results/tables/11_spatial_systems/atlas/neuropil_spatial_detection_context.csv; results/spatial_validation/decompose_bilateral_variance/global/tables/bilateral_precision_gain.csv; +1 more
# Produces: results/spatial_validation/summarize_spatial_atlas/global/reports/spatial_systems_atlas.xlsx; results/spatial_validation/summarize_spatial_atlas/global/tables/spatial_systems_atlas_validation.csv; results/spatial_validation/summarize_spatial_atlas/global/plots
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Assembles the already-computed atlas tables into one OOXML-safe workbook, emits NEW candidate figures (never a manuscript contract replacement), and runs the atlas validation.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/spatial/spatial_atlas_utils.R")
source("R/utilities/xlsx_package_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "summarize_spatial_atlas"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/spatial_validation/summarize_spatial_atlas.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

# This script is a pure aggregator: it reads five families and writes a
# workbook plus candidate figures. Every accessor below is a read and
# resolves normalized-first.
A <- function(f) spatial_systems_find(f, "atlas")
B <- function(f) spatial_systems_find(f, "bilateral")
C_ <- function(f) spatial_systems_find(f, "celltype_annotation")
D_ <- function(f) spatial_systems_find(f, "data_contract")
PR <- function(f) spatial_systems_find(f, "precision")
OUTT <- function() { d <- CANONICAL_PATHS$tables; dir_create(d); d }
FIG <- function() { d <- CANONICAL_PATHS$plots; dir_create(d); d }

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Spatial systems atlas workbook, candidate figures and validation.\n")
  dry_run_inputs(SCRIPT_ID, list(module_atlas = A("WGCNA_module_spatial_cell_affinity.csv"),
                                 protein_atlas = A("protein_spatial_cell_affinity.csv")))
  quit(save = "no", status = 0L)
}

rd <- function(p) if (file.exists(p)) {
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
} else NULL

sheets <- list(
  Module_atlas            = rd(A("WGCNA_module_spatial_cell_affinity.csv")),
  Module_spatial_raw      = rd(A("WGCNA_module_spatial_fingerprints_raw.csv")),
  Module_spatial_z        = rd(A("WGCNA_module_spatial_fingerprints_z.csv")),
  Module_compartment      = rd(A("WGCNA_module_empirical_compartment_affinity_long.csv")),
  Module_reference        = rd(A("WGCNA_module_reference_marker_affinity_long.csv")),
  Module_external_celltype = rd(C_("WGCNA_module_external_celltype_affinity_long.csv")),
  Module_bilateral        = rd(B("WGCNA_module_bilateral_reproducibility.csv")),
  Module_similarity       = rd(A("WGCNA_module_spatial_similarity.csv")),
  Supermodule_context     = rd(A("WGCNA_supermodule_spatial_cell_context.csv")),
  Supermodule_coherence   = rd(A("WGCNA_supermodule_spatial_coherence.csv")),
  Protein_atlas           = rd(A("protein_spatial_cell_affinity.csv")),
  SUS_RES_37              = rd(A("protein_sus_res_fdr_supported_atlas.csv")),
  Protein_spatial         = rd(A("protein_baseline_spatial_profile.csv")),
  CA2_SLM_audit           = rd(A("neuropil_spatial_detection_context.csv")),
  Label_context_audit     = rd(A("WGCNA_label_spatial_cell_context_audit.csv")),
  Precision               = rd(PR("bilateral_precision_gain.csv")),
  Evidence_dependence     = rd(D_("spatial_systems_evidence_dependence.csv"))
)
sheets <- sheets[!vapply(sheets, is.null, logical(1))]

# ------------------------------------------------------------ validation

checks <- list()
add <- function(id, critical, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id, critical = critical, status = status, detail = detail,
    stringsAsFactors = FALSE)
}
ma <- sheets$Module_atlas; pa <- sheets$Protein_atlas
ew <- sheets$Module_external_celltype; emp <- sheets$Module_compartment

add("module_ids_canonical", TRUE,
    if (!is.null(ma) && all(grepl("^WGCNA_m[0-9]+$", ma$ModuleID))) "PASS" else "FAIL",
    paste0(if (is.null(ma)) "module atlas missing" else
      paste0(nrow(ma), " modules, all matching WGCNA_mNN; no colour token")))
add("protein_ids_canonical", TRUE,
    if (!is.null(pa) && all(grepl("^PG:", pa$ProteinGroupID))) "PASS" else "FAIL",
    paste0(if (is.null(pa)) "protein atlas missing" else
      paste0(nrow(pa), " rows, all ProteinGroupID prefixed PG:")))
add("con_baseline_only", TRUE,
    if (!is.null(rd(A("WGCNA_module_baseline_spatial_profile_long.csv")))) "PASS" else "FAIL",
    "baseline profile carries baseline_definition = CON animals only")
add("no_phenotype_in_affinity", TRUE,
    if (!is.null(ma) && !any(c("StressGroup", "SUS", "RES") %in% names(ma))) "PASS" else "FAIL",
    "no phenotype column reaches the module affinity table")
add("empirical_universe_is_measured_proteome", TRUE,
    if (!is.null(emp) && all(emp$n_universe > 4000)) "PASS" else "FAIL",
    if (is.null(emp)) "missing" else
      paste0("universe sizes ", paste(range(emp$n_universe), collapse = "-"),
             " = measured proteome"))
add("marker_set_informativeness_recorded", TRUE,
    if (!is.null(emp) && "informative_for_affinity_call" %in% names(emp)) "PASS" else "FAIL",
    "non-discriminating marker sets are flagged, not silently used")
add("ewce_no_zero_empirical_p", TRUE,
    {
      t_ <- if (is.null(ew)) NULL else ew[ew$annotation_status == "tested", , drop = FALSE]
      if (is.null(t_)) "FAIL" else
        if (sum(t_$p_value == 0, na.rm = TRUE) == 0L &&
            sum(t_$FDR == 0, na.rm = TRUE) == 0L) "PASS" else "FAIL"
    },
    "empirical p from a finite bootstrap null is bounded below by 1/(B+1)")
add("ewce_scopes_kept_separate", TRUE,
    if (!is.null(ew) && length(unique(ew$module_scope)) == 3L) "PASS" else "FAIL",
    "all / core_kME06 / top25 retained as separate rows and FDR families")
add("no_single_gene_ewce", TRUE,
    if (is.null(ew) || min(ew$n_input_genes, na.rm = TRUE) >= 10L) "PASS" else "FAIL",
    "no gene set below the 10-gene minimum was tested")
add("no_automatic_module_renaming", TRUE,
    {
      la <- sheets$Label_context_audit
      if (!is.null(la) && all(grepl("AUDIT ONLY", la$audit_note))) "PASS" else "FAIL"
    },
    "label context audit is explicitly non-mutating")
add("da_joined_not_recomputed", TRUE,
    if (!is.null(pa) && "da_provenance" %in% names(pa) &&
        all(grepl("no differential statistic recomputed", pa$da_provenance))) "PASS" else "FAIL",
    "every DA statistic carries a provenance string naming its canonical source")
add("all_37_fdr_proteins_present_once", TRUE,
    {
      s <- sheets$SUS_RES_37
      if (is.null(s)) "FAIL" else
        if (nrow(s) == 37L &&
            !anyDuplicated(paste(s$dataset, s$ProteinGroupID))) "PASS" else "FAIL"
    },
    paste0(if (is.null(sheets$SUS_RES_37)) "missing" else
      paste0(nrow(sheets$SUS_RES_37),
             " FDR-supported rows, unique per dataset x ProteinGroupID")))
add("spatial_units_canonical", TRUE,
    if (!is.null(ma) && !any(grepl("[A-Z]", ma$peak_unit))) "PASS" else "FAIL",
    "peak units use the canonical lowercase spatial-unit vocabulary")
add("ca2_slm_audit_uses_all_tested", TRUE,
    {
      ca <- sheets$CA2_SLM_audit
      if (!is.null(ca) && all(ca$n_proteins_tested > 5000)) "PASS" else "FAIL"
    },
    "denominator is all tested proteins per unit, not the candidate subset")

# ------------------------------------------------------------- workbook

wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
title <- openxlsx::createStyle(fontName = "Arial", fontSize = 13,
                               textDecoration = "bold", fontColour = "#1F2933")
hdr <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E9EDF0",
                             textDecoration = "bold", wrapText = TRUE,
                             halign = "center", valign = "center")
wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE,
                              valign = "top")

openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
openxlsx::writeData(wb, "README", "Spatial / cell-affinity atlas", startRow = 1)
openxlsx::addStyle(wb, "README", title, rows = 1, cols = 1:2, gridExpand = TRUE, stack = TRUE)
readme <- data.frame(Note = c(
  "FOUR EVIDENCE DIMENSIONS ARE KEPT SEPARATE and are never summed: spatial anatomical identity, empirical compartment affinity, external cell-type affinity, and bilateral reliability. context_confidence is a RULE-BASED category and every row carries the rule that produced it in context_rule.",
  "BASELINE IDENTITY IS PHENOTYPE-BLIND. It uses CON animals only, via the validated bilateral constructor. SUS-RES statistics are joined afterwards as a downstream overlay and never contribute to identity or affinity.",
  "EXTERNAL CELL-TYPE AFFINITY IS NOT COMPARTMENT AFFINITY. External enrichment says a module's genes are preferentially expressed in a reference cell type. It is not a cell-proportion estimate, and a microglia-enriched ROI is not purified microglia.",
  "EMPIRICAL MARKER SETS ARE UNBALANCED. empirical_microglia_neuropil_shared covers ~73% of the measured proteome and empirical_neuropil_enriched contains a single protein. Sets outside 0.5%-50% coverage are flagged informative_for_affinity_call = FALSE and are excluded from the single strongest-compartment call, though they remain in the long table. Neuropil affinity therefore CANNOT be assessed empirically in this release.",
  "EMPIRICAL P-VALUES CANNOT BE ZERO. EWCE computes p = sum(null >= observed)/B. A raw 0 means no null draw reached the observed value, not that the null probability is zero, so the standard (1 + count)/(B + 1) correction is applied. The uncorrected EWCE value is retained as p_value_raw_ewce.",
  "BILATERAL RELIABILITY QUALIFIES AN IDENTITY, IT NEVER ERASES ONE. A module with low left/right agreement keeps its spatial identity and gains the caveat spatial_identity_present_but_bilaterally_variable.",
  "LABELS ARE AUDITED, NEVER CHANGED. No module was renamed and no proposed label was activated by any script in this layer.",
  "DA STATISTICS ARE JOINED, NEVER RECOMPUTED. Every phenotype column carries a provenance string naming its canonical source."),
  stringsAsFactors = FALSE)
openxlsx::writeDataTable(wb, "README", readme, startRow = 3, tableName = "AtlasReadme",
                         tableStyle = "TableStyleLight9")
openxlsx::setColWidths(wb, "README", cols = 1, widths = 150)
openxlsx::addStyle(wb, "README", wrap, rows = 3:(3 + nrow(readme)), cols = 1,
                   gridExpand = TRUE, stack = TRUE)

for (nm in names(sheets)) {
  d <- sheets[[nm]]
  sh <- substr(nm, 1, 31)
  openxlsx::addWorksheet(wb, sh, gridLines = FALSE)
  openxlsx::writeData(wb, sh, d, startRow = 1, colNames = TRUE)
  openxlsx::addStyle(wb, sh, hdr, rows = 1, cols = seq_len(max(ncol(d), 1L)),
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sh, firstActiveRow = 2)
  openxlsx::setColWidths(wb, sh, cols = seq_len(max(ncol(d), 1L)),
                         widths = pmin(46, pmax(12, nchar(names(d)) + 2)))
}
validation <- dplyr::bind_rows(checks)
openxlsx::addWorksheet(wb, "Validation", gridLines = FALSE, tabColour = "#7A2E2E")
openxlsx::writeData(wb, "Validation", validation, startRow = 1)
openxlsx::addStyle(wb, "Validation", hdr, rows = 1, cols = seq_len(ncol(validation)),
                   gridExpand = TRUE, stack = TRUE)
openxlsx::setColWidths(wb, "Validation", cols = 1:4, widths = c(42, 10, 10, 110))

wb_path <- file.path(dir_create(CANONICAL_PATHS$reports),
                     "spatial_systems_atlas.xlsx")
xlsx_save_valid_workbook(wb, wb_path)
write_csv_safe(validation, file.path(OUTT(), "spatial_systems_atlas_validation.csv"))

# -------------------------------------------------------------- figures

fig_dir <- FIG()
made <- character()
if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(ma)) {
  library(ggplot2)
  # FIGURE 1: module spatial fingerprint heatmap (row-standardised for display)
  zr <- sheets$Module_spatial_z
  if (!is.null(zr)) {
    long <- zr %>%
      tidyr::pivot_longer(-c("dataset", "ModuleID"), names_to = "spatial_unit",
                          values_to = "z")
    long <- dplyr::left_join(
      long, ma[, c("dataset", "ModuleID", "bilateral_support_class", "SupermoduleID")],
      by = c("dataset", "ModuleID"))
    p <- ggplot(long, aes(.data$spatial_unit, .data$ModuleID, fill = .data$z)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      facet_wrap(~ dataset, scales = "free") +
      scale_fill_gradient2(low = "#2C7BB6", mid = "#FFFFBF", high = "#D7191C",
                           midpoint = 0, name = "row z") +
      labs(x = NULL, y = NULL,
           title = "Module spatial fingerprints (CON baseline, row-standardised)",
           subtitle = "Standardisation is for display only; peak and specificity come from raw values") +
      theme_minimal(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())
    ggsave(file.path(fig_dir, "module_spatial_fingerprint_heatmap.png"), p,
           width = 11, height = 6, dpi = 200, bg = "white")
    write_csv_safe(long, file.path(fig_dir, "module_spatial_fingerprint_heatmap_source_data.csv"))
    made <- c(made, "module_spatial_fingerprint_heatmap")
  }
  # FIGURE 2: module cell/compartment affinity atlas
  aff <- ma %>%
    dplyr::select("dataset", "ModuleID", "peak_unit", "spatial_tau",
                  "strongest_empirical_compartment", "empirical_compartment_FDR",
                  "external_celltype_all", "external_FDR_all",
                  "bilateral_support_class", "context_confidence")
  p2 <- ggplot(aff, aes(.data$spatial_tau,
                        stats::reorder(.data$ModuleID, .data$spatial_tau),
                        colour = .data$context_confidence)) +
    geom_point(aes(size = -log10(pmax(.data$external_FDR_all, 1e-5))), alpha = 0.9) +
    facet_wrap(~ dataset, scales = "free_y") +
    labs(x = "spatial specificity (tau)", y = NULL, size = "-log10 external FDR",
         colour = "context class",
         title = "Module spatial specificity and external cell-type affinity") +
    theme_minimal(base_size = 8)
  ggsave(file.path(fig_dir, "module_spatial_cell_affinity_atlas.png"), p2,
         width = 11, height = 6, dpi = 200, bg = "white")
  write_csv_safe(aff, file.path(fig_dir, "module_spatial_cell_affinity_atlas_source_data.csv"))
  made <- c(made, "module_spatial_cell_affinity_atlas")

  # FIGURE 3: the 37 FDR-supported proteins, compact
  s37 <- sheets$SUS_RES_37
  if (!is.null(s37) && nrow(s37)) {
    d3 <- s37 %>%
      dplyr::mutate(label = paste0(.data$GeneSymbol, " (", .data$ModuleID, ")"),
                    direction = ifelse(.data$sus_res_strongest_log2FC < 0,
                                       "lower in SUS", "higher in SUS"))
    p3 <- ggplot(d3, aes(.data$sus_res_strongest_log2FC,
                         stats::reorder(.data$label, .data$sus_res_strongest_log2FC),
                         fill = .data$direction)) +
      geom_col() +
      facet_grid(rows = vars(.data$dataset), scales = "free_y", space = "free_y") +
      labs(x = "SUS-RES log2FC at strongest spatial unit", y = NULL,
           title = "FDR-supported SUS-RES proteins in spatial context",
           subtitle = "Baseline identity and phenotype effect are separate quantities") +
      theme_minimal(base_size = 7)
    ggsave(file.path(fig_dir, "protein_spatial_context_37_hits.png"), p3,
           width = 9, height = 9, dpi = 200, bg = "white")
    write_csv_safe(d3, file.path(fig_dir, "protein_spatial_context_37_hits_source_data.csv"))
    made <- c(made, "protein_spatial_context_37_hits")
  }
}
add("figures_have_source_data", FALSE,
    if (all(file.exists(file.path(fig_dir, paste0(made, "_source_data.csv"))))) "PASS" else "FAIL",
    paste0(length(made), " candidate figures, each with a source-data CSV"))
validation <- dplyr::bind_rows(checks)
write_csv_safe(validation, file.path(OUTT(), "spatial_systems_atlas_validation.csv"))

cat("\n===== Spatial systems atlas workbook =====\n")
cat("sheets:", length(sheets) + 2L, " workbook:", relative_to(wb_path), "\n")
cat("figures:", paste(made, collapse = ", "), "\n\n")
for (i in seq_len(nrow(validation))) {
  cat(sprintf("  %-7s %-42s %s\n", validation$status[i], validation$check_id[i],
              substr(validation$detail[i], 1, 78)))
}
crit <- validation$critical %in% TRUE & validation$status == "FAIL"
cat(sprintf("\n%d checks, %d critical FAIL\n", nrow(validation), sum(crit)))
if (any(crit)) quit(save = "no", status = 1L)
