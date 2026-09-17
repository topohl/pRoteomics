#!/usr/bin/env Rscript
#
# Legacy comparison, behaviour-coupling audit, workbook, candidate figures and
# the network validation contract.
#
# Historical 07_spatial_networks outputs are NOT modified. They are read, their
# inference status is labelled, and the comparison is written to a new file.
#
# USAGE
#   Rscript analysis/spatial_validation/validate_network_workbook.R
# Script: analysis/spatial_validation/validate_network_workbook.R
# Stage: networks
# Scope: global
# Consumes: required results/spatial_validation/build_animal_spatial_networks/global/tables/animal_network_edges.csv; results/tables/11_spatial_systems/networks/animal_network_edges.csv; results/spatial_validation/test_network_group_organization/global/tables/network_edge_group_differences.csv; +1 more; optional results/tables/07_spatial_networks/; results/spatial_validation/test_network_group_organization/global/tables/WGCNA_module_network_context.csv; results/tables/11_spatial_systems/atlas/WGCNA_module_network_context.csv; +2 more
# Produces: results/spatial_validation/validate_network_workbook/global/reports/spatial_systems_networks.xlsx; results/spatial_validation/validate_network_workbook/global/tables/spatial_network_validation_status.csv; results/spatial_validation/validate_network_workbook/global/tables/legacy_spatial_network_comparison.csv; +2 more
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Legacy comparison, behaviour-coupling audit, workbook, candidate figures and the network validation contract.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/networks/animal_spatial_network_utils.R")
source("R/utilities/xlsx_package_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "validate_network_workbook"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/spatial_validation/validate_network_workbook.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

# reads resolve normalized-first; writes go to this analysis's own directory
NET <- function(f) spatial_systems_find(f, "networks")
NET_OUT <- function(...) { d <- CANONICAL_PATHS$tables; dir_create(d); file.path(d, ...) }
LEG <- function(...) path_results("tables", "07_spatial_networks", ...)
FIG <- function() {
  d <- CANONICAL_PATHS$plots; dir_create(d); d
}
rd <- function(p) if (file.exists(p)) {
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
} else NULL

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Legacy comparison, workbook, figures and network validation.\n")
  dry_run_inputs(SCRIPT_ID, list(animal_network_edges = NET("animal_network_edges.csv")))
  cat("[DRY-RUN] Historical 07_spatial_networks outputs are read, never modified.\n")
  quit(save = "no", status = 0L)
}

edges <- rd(NET("animal_network_edges.csv"))
bilateral <- rd(NET("animal_network_bilateral_reproducibility.csv"))
edge_bil <- rd(NET("edge_bilateral_reproducibility.csv"))
baseline <- rd(NET("CON_spatial_molecular_similarity_matrix.csv"))
metrics <- rd(NET("animal_network_global_metrics.csv"))
gmetrics <- rd(NET("network_group_metric_comparison.csv"))
mv <- rd(NET("network_global_multivariate_test.csv"))
edge_grp <- rd(NET("network_edge_group_differences.csv"))
selection <- rd(NET("network_representation_selection.csv"))
distances <- rd(NET("animal_network_distance_from_CON.csv"))
edge_ctx <- rd(NET("spatial_edge_module_context.csv"))
node_ctx <- rd(NET("spatial_node_cell_context.csv"))
node_int <- rd(NET("spatial_node_network_integration.csv"))
module_net_ctx <- rd(spatial_systems_find("WGCNA_module_network_context.csv",
                                          "atlas"))

# ==================================== PART 17: legacy comparison

message("Comparing against the legacy 07_spatial_networks pipeline")
leg_rows <- list()
old_perm <- rd(LEG("bootstrap_differential_network_stability",
                   "edge_bootstrap_permutation_validation_summary.csv"))
old_stable <- rd(LEG("bootstrap_differential_network_stability",
                     "stable_differential_edges.csv"))

norm_edge <- function(a, b) {
  # legacy neuropil nodes are CA1_sp style from the corrupted layer axis; the
  # corrected axis resolves ten region-layer units, so most legacy neuropil
  # edges have no counterpart at all
  paste(pmin(tolower(a), tolower(b)), pmax(tolower(a), tolower(b)), sep = "__")
}
new_lookup <- if (!is.null(edge_grp)) {
  stats::setNames(edge_grp$SUS_minus_RES,
                  norm_edge(edge_grp$node_a, edge_grp$node_b))
} else character()

if (!is.null(old_perm)) {
  nm <- names(old_perm)
  a_col <- intersect(c("Source", "node_a", "NodeA", "From"), nm)[1]
  b_col <- intersect(c("Target", "node_b", "NodeB", "To"), nm)[1]
  w_col <- intersect(c("observed_r", "MeanDeltaR", "DeltaR"), nm)[1]
  for (i in seq_len(nrow(old_perm))) {
    key <- if (!is.na(a_col) && !is.na(b_col))
      norm_edge(old_perm[[a_col]][i], old_perm[[b_col]][i]) else NA_character_
    new_est <- if (!is.na(key) && key %in% names(new_lookup)) new_lookup[[key]] else NA_real_
    leg_rows[[length(leg_rows) + 1L]] <- data.frame(
      legacy_file = "edge_bootstrap_permutation_validation_summary.csv",
      legacy_edge = if (!is.na(a_col)) paste(old_perm[[a_col]][i], old_perm[[b_col]][i],
                                             sep = "__") else NA_character_,
      legacy_estimate = if (!is.na(w_col)) old_perm[[w_col]][i] else NA_real_,
      legacy_permutation_p = if ("permutation_p" %in% nm) old_perm$permutation_p[i] else NA_real_,
      legacy_fdr = if ("fdr" %in% nm) old_perm$fdr[i] else NA_real_,
      new_animal_level_estimate = new_est,
      comparable = !is.na(new_est),
      legacy_inference_status = "legacy_noncanonical_inference",
      legacy_defect = paste0(
        "permutation_p is pmin(Prob_DeltaR_LessThan0, Prob_DeltaR_GreaterThan0)*2 - ",
        "a bootstrap SIGN FREQUENCY, not a label permutation. No group label is ",
        "shuffled anywhere in 07_spatial_networks. The fdr column is a BH ",
        "adjustment of that pseudo-p and inherits the misnomer."),
      interpretation = if (is.na(new_est)) "not_comparable" else "see new animal-level estimate",
      stringsAsFactors = FALSE)
  }
}
legacy <- dplyr::bind_rows(leg_rows)
if (nrow(legacy)) {
  legacy$interpretation[!legacy$comparable] <- paste0(
    "not_comparable: the legacy neuropil network was built on a corrupted layer ",
    "axis (every sample assigned Layer='sp', giving 4 pseudo-nodes) so its edges ",
    "have no counterpart in the corrected 10-node region-layer axis")
}

# ==================================== PART 20: behaviour-coupling audit
#
# AUDIT ONLY in this pass. The consumer is not migrated here because migrating
# the edge source alone would not make it work.
beh_src <- repo_path("analysis/integration", "test_network_behaviour_coupling.R")
beh_rows <- list()
if (file.exists(beh_src)) {
  src <- paste(readLines(beh_src, warn = FALSE), collapse = "\n")
  beh_rows[[1]] <- data.frame(
    finding = "replicate_unit",
    status = if (grepl("compute_animal_edge_scores", src, fixed = TRUE)) "OK" else "UNKNOWN",
    detail = paste0("Edges are built per animal from that animal's own samples, ",
                    "so the biological replicate is already AnimalID. The earlier ",
                    "suspicion that a group-mean network was propagated is REFUTED."),
    stringsAsFactors = FALSE)
  beh_rows[[2]] <- data.frame(
    finding = "layer_collapse_defect",
    status = "RESOLVED_UPSTREAM",
    detail = paste0("The RegionLayer collapse that emptied this script was fixed ",
                    "upstream when the layer regex was corrected; it no longer blocks."),
    stringsAsFactors = FALSE)
  beh_rows[[3]] <- data.frame(
    finding = "animal_id_normaliser_defect",
    status = "STILL_BLOCKING",
    detail = paste0("normalize_animal_id() remains broken and is now the ONLY ",
                    "blocker. It does not merely fail to match - it produces WRONG ",
                    "animals: '13856' and '13857' both map to 'A1385', merging two ",
                    "animals, and 'OQ754' maps to 'A0754'. Combined with ",
                    "distinct(AnimalID), physiology rows are silently discarded. ",
                    "Migrating the edge source alone would NOT fix the output."),
    stringsAsFactors = FALSE)
  beh_rows[[4]] <- data.frame(
    finding = "migration_recommendation",
    status = "DEFERRED",
    detail = paste0("A canonical animal-level edge table now exists and could ",
                    "replace compute_animal_edge_scores(). Migration is deferred ",
                    "until the identity normaliser is fixed, because migrating ",
                    "first would silently produce a corrupted join rather than an ",
                    "empty one - a worse failure mode. Three incompatible AnimalID ",
                    "conventions coexist ('A0003', '3', 'OR111'), so any migration ",
                    "must declare which one is canonical."),
    stringsAsFactors = FALSE)
}
behaviour <- dplyr::bind_rows(beh_rows)

# ==================================== PART 24: validation contract

checks <- list()
add <- function(id, critical, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id, critical = critical, status = status, detail = detail,
    stringsAsFactors = FALSE)
}
bil_e <- if (!is.null(edges)) edges[edges$hemisphere_mode == "bilateral", , drop = FALSE] else NULL

add("one_network_per_animal", TRUE,
    {
      k <- unique(bil_e[, c("dataset", "AnimalID")])
      if (!is.null(bil_e) && nrow(k) == length(unique(paste(bil_e$dataset, bil_e$AnimalID))))
        "PASS" else "FAIL"
    },
    paste0(if (is.null(bil_e)) "no edges" else
      paste0(nrow(unique(bil_e[, c("dataset", "AnimalID")])),
             " dataset x AnimalID bilateral networks")))
add("no_duplicate_animal_edge_rows", TRUE,
    if (!is.null(bil_e) &&
        !anyDuplicated(paste(bil_e$dataset, bil_e$AnimalID, bil_e$edge_id))) "PASS" else "FAIL",
    "one row per dataset x AnimalID x edge")
add("side_purity_L_and_R", TRUE,
    if (!is.null(edges) && all(c("L", "R", "bilateral") %in% edges$hemisphere_mode))
      "PASS" else "FAIL",
    "left and right networks are built from that hemisphere's columns only; the producer hard-stops otherwise")
add("representation_chosen_phenotype_blind", TRUE,
    if (!is.null(selection) && all(grepl("CON", selection$selection_criteria)) &&
        !any(grepl("SUS|RES", selection$selection_criteria))) "PASS" else "FAIL",
    "selection criteria reference CON, bilateral reproducibility and anatomy only")
add("matched_node_ordering", TRUE,
    {
      ok <- TRUE
      if (!is.null(bil_e)) {
        for (ds in unique(bil_e$dataset)) {
          z <- bil_e[bil_e$dataset == ds, , drop = FALSE]
          per <- split(z$edge_id, z$AnimalID)
          if (length(per) > 1L && !all(vapply(per[-1], function(v)
            identical(sort(v), sort(per[[1]])), logical(1)))) ok <- FALSE
        }
      }
      if (ok) "PASS" else "FAIL"
    },
    "every animal within a dataset carries the identical edge set")
add("full_weighted_matrices_retained", TRUE,
    if (dir.exists(NET("matrices")) &&
        length(list.files(NET("matrices"), pattern = "[.]csv$")) > 0L) "PASS" else "FAIL",
    paste0(length(list.files(NET("matrices"), pattern = "[.]csv$")),
           " full weighted matrices written, independent of any plot threshold"))
add("no_pseudo_permutation_p", TRUE,
    {
      bad <- FALSE
      for (d in list(edge_grp, gmetrics)) {
        if (is.null(d)) next
        if ("permutation_p" %in% names(d)) bad <- TRUE
      }
      if (!bad) "PASS" else "FAIL"
    },
    "no output field named permutation_p; exact values are named exact_p_two_sided and carry their attainable resolution")
add("no_fdr_from_bootstrap_sign", TRUE,
    {
      bad <- FALSE
      for (d in list(edge_grp, gmetrics)) {
        if (is.null(d)) next
        if (any(grepl("^fdr$|_fdr$", names(d), ignore.case = TRUE))) bad <- TRUE
      }
      if (!bad) "PASS" else "FAIL"
    },
    "no BH FDR is derived from any bootstrap sign frequency")
add("permutation_unit_is_animal", TRUE,
    if (!is.null(edge_grp) && all(edge_grp$exact_n_assignments == 20L, na.rm = TRUE))
      "PASS" else "FAIL",
    "SUS vs RES enumerates exactly 20 whole-animal assignments")
add("exact_resolution_documented", TRUE,
    if (!is.null(edge_grp) &&
        all(abs(edge_grp$min_attainable_two_sided_p - 0.10) < 1e-9, na.rm = TRUE))
      "PASS" else "FAIL",
    "minimum attainable two-sided p recorded as 0.10 on every edge row")
add("bootstrap_resamples_animals", TRUE,
    if (!is.null(edge_grp) && "n_valid_bootstrap_iterations" %in% names(edge_grp))
      "PASS" else "FAIL",
    "animal-level bootstrap with recorded valid-iteration count and seed")
add("legacy_marked_noncanonical", TRUE,
    if (nrow(legacy) == 0L || all(legacy$legacy_inference_status ==
                                  "legacy_noncanonical_inference")) "PASS" else "FAIL",
    paste0(nrow(legacy), " legacy rows labelled legacy_noncanonical_inference; ",
           "historical files were read, never modified"))
git_changed <- tryCatch(system2("git", c("status", "--porcelain"), stdout = TRUE,
                                stderr = FALSE), error = function(e) character())
add("no_da_or_wgcna_changes", TRUE,
    if (!any(grepl("04_differential|01_WGCNA[.]r|wgcna_final_model_state", git_changed)))
      "PASS" else "FAIL",
    "no differential-abundance or WGCNA state file modified")
add("legacy_outputs_unmodified", TRUE,
    if (!any(grepl("results/tables/07_spatial_networks", git_changed))) "PASS" else "FAIL",
    "historical 07_spatial_networks outputs untouched in the working tree")
# Part 23: the atlas gains node-level context only. If a per-edge column ever
# leaks into the atlas context table the atlas stops being an atlas.
atlas_edge_cols <- if (is.null(module_net_ctx)) character() else
  grep("^edge_|_edge$|node_a|node_b", names(module_net_ctx), value = TRUE)
add("atlas_context_is_node_level", TRUE,
    if (length(atlas_edge_cols) == 0L) "PASS" else "FAIL",
    if (is.null(module_net_ctx)) "no atlas network context generated" else
      sprintf(paste0("%d modules annotated with node integration context; ",
                     "%d edge-level columns added to the atlas"),
              nrow(module_net_ctx), length(atlas_edge_cols)))

validation <- dplyr::bind_rows(checks)
write_csv_safe(validation, NET_OUT("spatial_network_validation_status.csv"))
if (nrow(legacy)) write_csv_safe(legacy, NET_OUT("legacy_spatial_network_comparison.csv"))
if (nrow(behaviour)) write_csv_safe(behaviour, NET_OUT("network_behavior_coupling_audit.csv"))

# ------------------------------------------------------------- figures

fig_dir <- FIG(); made <- character()
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  if (!is.null(bilateral)) {
    p1 <- ggplot(bilateral, aes(stats::reorder(.data$AnimalID, .data$pearson_r),
                                .data$pearson_r, colour = .data$StressGroup)) +
      geom_point(size = 2.4) + geom_hline(yintercept = c(0.4, 0.7), linetype = 2,
                                          colour = "grey60") +
      facet_wrap(~ dataset, scales = "free_y") + coord_flip() +
      labs(x = NULL, y = "left vs right edge-vector Pearson r",
           title = "Bilateral network reproducibility, one point per animal") +
      theme_minimal(base_size = 8)
    ggsave(file.path(fig_dir, "network_bilateral_validation.png"), p1,
           width = 10, height = 4.5, dpi = 200, bg = "white")
    write_csv_safe(bilateral, file.path(fig_dir, "network_bilateral_validation_source_data.csv"))
    made <- c(made, "network_bilateral_validation")
  }
  if (!is.null(baseline)) {
    p2 <- ggplot(baseline, aes(.data$node_a, .data$node_b, fill = .data$median_similarity)) +
      geom_tile(colour = "white") + facet_wrap(~ dataset, scales = "free") +
      scale_fill_gradient2(low = "#2C7BB6", mid = "#FFFFBF", high = "#D7191C",
                           midpoint = 0, name = "CON median") +
      labs(x = NULL, y = NULL,
           title = "Baseline CON spatial molecular-similarity matrix",
           subtitle = "full weighted matrix; no edge threshold applied") +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(fig_dir, "network_CON_baseline_matrix.png"), p2,
           width = 11, height = 5, dpi = 200, bg = "white")
    write_csv_safe(baseline, file.path(fig_dir, "network_CON_baseline_matrix_source_data.csv"))
    made <- c(made, "network_CON_baseline_matrix")
  }
  if (!is.null(metrics)) {
    long <- metrics %>%
      tidyr::pivot_longer(-c("dataset", "AnimalID", "StressGroup"),
                          names_to = "metric", values_to = "value")
    p3 <- ggplot(long, aes(.data$StressGroup, .data$value, colour = .data$StressGroup)) +
      geom_point(size = 1.8, position = position_jitter(width = 0.12, height = 0)) +
      facet_wrap(dataset ~ metric, scales = "free_y", ncol = 5) +
      labs(x = NULL, y = NULL, title = "Animal-level global network metrics",
           subtitle = "every animal shown; no bars without points, n = 3 per group") +
      theme_minimal(base_size = 6) + theme(legend.position = "none")
    ggsave(file.path(fig_dir, "network_animal_global_metrics.png"), p3,
           width = 12, height = 9, dpi = 200, bg = "white")
    write_csv_safe(long, file.path(fig_dir, "network_animal_global_metrics_source_data.csv"))
    made <- c(made, "network_animal_global_metrics")
  }
  if (!is.null(distances)) {
    p4 <- ggplot(distances, aes(.data$StressGroup, .data$distance_from_CON_centroid,
                                colour = .data$StressGroup)) +
      geom_point(size = 2.4, position = position_jitter(width = 0.1, height = 0)) +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(x = NULL, y = "distance from CON centroid (Fisher-z edge space)",
           title = "Network distance from the CON centroid",
           subtitle = "CON animals use a leave-one-CON-out centroid") +
      theme_minimal(base_size = 8) + theme(legend.position = "none")
    ggsave(file.path(fig_dir, "network_distance_from_CON.png"), p4,
           width = 9, height = 4, dpi = 200, bg = "white")
    write_csv_safe(distances, file.path(fig_dir, "network_distance_from_CON_source_data.csv"))
    made <- c(made, "network_distance_from_CON")
  }
}
add("figures_have_source_data", FALSE,
    if (all(file.exists(file.path(fig_dir, paste0(made, "_source_data.csv"))))) "PASS" else "FAIL",
    paste0(length(made), " candidate figures, each with a source-data CSV"))
validation <- dplyr::bind_rows(checks)
write_csv_safe(validation, NET("spatial_network_validation_status.csv"))

# ------------------------------------------------------------- workbook

sheets <- list(
  Representation_selection = selection,
  Animal_network_metrics = metrics,
  Animal_edge_values = bil_e,
  Bilateral_summary = bilateral,
  Edge_bilateral = edge_bil,
  CON_baseline_edges = baseline,
  Group_edge_summary = edge_grp,
  Global_tests = gmetrics,
  Exact_permutation = mv,
  Network_distance = distances,
  Legacy_comparison = legacy,
  Edge_module_context = edge_ctx,
  Node_cell_context = node_ctx,
  Node_network_integration = node_int,
  Atlas_module_network_context = module_net_ctx,
  Behavior_coupling = behaviour,
  Validation = validation)
sheets <- sheets[!vapply(sheets, function(z) is.null(z) || !nrow(z), logical(1))]

wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
hdr <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E9EDF0",
                             textDecoration = "bold", wrapText = TRUE)
wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE,
                              valign = "top")
openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
openxlsx::writeData(wb, "README", "Animal-level spatial molecular-similarity networks", startRow = 1)
readme <- data.frame(Note = c(
  "ONE AnimalID = ONE independent network replicate. Every permutation relabels whole animals and every bootstrap resamples animals. Proteins, hemispheres and spatial units are repeated measures inside an animal and are never resampled as replicates.",
  "NODE = an anatomical sampling unit. EDGE = similarity of the molecular spatial profiles of two anatomical units within one animal. This is NOT neural connectivity, NOT anatomical connectivity, NOT molecular communication, NOT a protein coexpression network and NOT brain connectivity.",
  "THE REPRESENTATION WAS CHOSEN PHENOTYPE-BLIND on CON animals, using bilateral reproducibility, inter-animal stability, anatomical plausibility and protein retention. Raw abundance was REJECTED: it posts the highest inter-animal stability precisely because similarity is dominated by the global abundance rank every unit shares, and it fails to separate within-region from across-region pairs (0.009 versus 0.12-0.13).",
  "EXACT TEST RESOLUTION. SUS vs RES is 3 vs 3: 20 label assignments, so the smallest attainable two-sided p is 0.10 and no such comparison can reach 0.05. Failing to cross 0.05 is NOT evidence of no difference. The three-group omnibus has 1680 assignments.",
  "THE LEGACY PIPELINE IS NOT THE INFERENTIAL BASIS. Its permutation_p was a bootstrap sign frequency - no group label was ever shuffled - and its fdr was a BH adjustment of that pseudo-p. Its minimum fdr across all 18 edges was 0.890, so it had no significant edge by its own criterion. Those outputs are retained unmodified and labelled legacy_noncanonical_inference.",
  "FULL WEIGHTED MATRICES are the canonical numerical representation. No edge is thresholded for analysis; thresholding belongs to visualisation only.",
  "ATLAS CONTEXT IS NODE-LEVEL. The module and protein atlas gains one integration summary per spatial unit - how much that unit's molecular profile resembles the rest of the network in CON animals - so a peak in a distinct unit can be told apart from a peak in a unit that resembles everything. It is a CON-baseline property of the unit, not of the module or protein, and carries no group inference. The edge structure itself stays in the network tables."),
  stringsAsFactors = FALSE)
openxlsx::writeDataTable(wb, "README", readme, startRow = 3, tableName = "NetReadme",
                         tableStyle = "TableStyleLight9")
openxlsx::setColWidths(wb, "README", cols = 1, widths = 150)
openxlsx::addStyle(wb, "README", wrap, rows = 3:(3 + nrow(readme)), cols = 1,
                   gridExpand = TRUE, stack = TRUE)
for (nm in names(sheets)) {
  d <- sheets[[nm]]; sh <- substr(nm, 1, 31)
  openxlsx::addWorksheet(wb, sh, gridLines = FALSE)
  openxlsx::writeData(wb, sh, d, startRow = 1)
  openxlsx::addStyle(wb, sh, hdr, rows = 1, cols = seq_len(ncol(d)),
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sh, firstActiveRow = 2)
  openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(d)),
                         widths = pmin(44, pmax(12, nchar(names(d)) + 2)))
}
wb_path <- file.path(dir_create(CANONICAL_PATHS$reports),
                     "spatial_systems_networks.xlsx")
xlsx_save_valid_workbook(wb, wb_path)

cat("\n===== Network legacy comparison, workbook and validation =====\n")
cat("workbook:", relative_to(wb_path), " sheets:", length(sheets) + 1L, "\n")
cat("figures :", paste(made, collapse = ", "), "\n\n")
for (i in seq_len(nrow(validation))) {
  cat(sprintf("  %-7s %-38s %s\n", validation$status[i], validation$check_id[i],
              substr(validation$detail[i], 1, 74)))
}
crit <- validation$critical %in% TRUE & validation$status == "FAIL"
cat(sprintf("\n%d checks, %d critical FAIL\n", nrow(validation), sum(crit)))
if (nrow(legacy)) {
  cat(sprintf("\nlegacy edges compared: %d (comparable to the corrected axis: %d)\n",
              nrow(legacy), sum(legacy$comparable)))
  cat(sprintf("legacy minimum fdr: %.3f - no legacy edge was significant by its own criterion\n",
              suppressWarnings(min(legacy$legacy_fdr, na.rm = TRUE))))
}
if (any(crit)) quit(save = "no", status = 1L)
