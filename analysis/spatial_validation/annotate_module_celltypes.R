#!/usr/bin/env Rscript
#
#
# Uses the shared engine in R/enrichment/ewce_gene_set_engine.R - the SAME EWCE method and
# the SAME specificity reference as the canonical analysis. Nothing about the
# test is reimplemented.
#
# PHENOTYPE-BLIND BY CONSTRUCTION
#   A module's gene set is its membership. No StressGroup, contrast, CON/RES/SUS
#   label or DAP status is read anywhere in this script, and the engine rejects
#   gene-set frames that carry phenotype columns.
#
# FDR FAMILIES ARE ISOLATED
#   Each row carries ewce_module_annotation_<dataset>_<scope>_<level> and BH is
#   applied strictly within it, across ModuleID x CellType. Adding or removing
#   rows from the phenotype arm can never move these numbers.
#
# THIS IS EXTERNAL EXPRESSION ENRICHMENT, NOT COMPARTMENT AFFINITY
#   It says a module's genes are preferentially expressed in a reference cell
#   type. It is not a cell-proportion estimate, and a microglia-enriched ROI is
#   not purified microglia.
#
# WGCNA LABELS ARE NOT ALTERED BY THIS SCRIPT.
#
# USAGE
#   Rscript analysis/spatial_validation/annotate_module_celltypes.R
#   Rscript analysis/spatial_validation/annotate_module_celltypes.R --dry-run
#   Rscript analysis/spatial_validation/annotate_module_celltypes.R --reps 1000
# Script: analysis/spatial_validation/annotate_module_celltypes.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; optional none declared in pipeline.yml
# Produces: results/tables/11_spatial_systems/celltype_annotation/WGCNA_module_external_celltype_affinity_long.csv; results/tables/11_spatial_systems/celltype_annotation/WGCNA_module_external_celltype_affinity_summary.csv
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Phenotype-blind external cell-type annotation of WGCNA modules.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/enrichment/ewce_gene_set_engine.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/spatial_validation/annotate_module_celltypes.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

args <- commandArgs(trailingOnly = TRUE)
REPS <- {
  i <- match("--reps", args)
  if (!is.na(i) && length(args) > i) as.integer(args[[i + 1L]]) else 10000L
}
LEVELS <- c(1L, 2L)
SEED <- 20260101L

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "celltype_annotation"); dir_create(d); d
}
membership_path <- function(ds) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules",
               "WGCNA_modules_long.csv")
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- stats::setNames(lapply(DATASETS, membership_path),
                            paste0("wgcna_membership__", DATASETS))
  cat("[DRY-RUN] Phenotype-blind WGCNA module external cell-type annotation.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Scopes: ", paste(ewce_module_scopes(), collapse = ", "), "\n")
  cat("[DRY-RUN] reps=", REPS, " levels=", paste(LEVELS, collapse = ","), "\n", sep = "")
  cat("[DRY-RUN] No StressGroup, contrast or DAP status is read.\n")
  quit(save = "no", status = 0L)
}

if (!requireNamespace("EWCE", quietly = TRUE) ||
    !requireNamespace("ewceData", quietly = TRUE)) {
  stop("missing_required_input: EWCE and ewceData are required.", call. = FALSE)
}

message("Loading external specificity reference")
ctd <- ewceData::ctd()

# ----------------------------------------------------- module gene sets

# The WGCNA membership tables store gene symbols UPPERCASED. The specificity
# reference is mouse Title-case, and an uppercased symbol is not a valid
# org.Mm.eg.db key, so passing them straight through yields zero reference
# matches. Resolve the canonical mouse symbol once, here.
sym_of <- function(x) {
  s <- if ("official_gene_symbol" %in% names(x)) x$official_gene_symbol else x$GeneSymbol
  ewce_to_mouse_symbols(s)
}

build_scope_sets <- function(m, scope) {
  keep <- switch(scope,
                 all = rep(TRUE, nrow(m)),
                 core_kME06 = m[["is_core_kME_0.6"]] %in% TRUE,
                 top25 = m[["is_top_hub_25"]] %in% TRUE)
  z <- m[keep, , drop = FALSE]
  s <- sym_of(z)
  ok <- !is.na(s) & nzchar(s)
  split(s[ok], as.character(z$ModuleID[ok]))
}

all_rows <- list()
for (ds in DATASETS) {
  p <- membership_path(ds)
  if (!file.exists(p)) stop("missing_required_input: WGCNA membership: ", p, call. = FALSE)
  m <- as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                     guess_max = Inf))
  # PHENOTYPE-BLIND GUARD: drop any group/contrast column before use so a
  # phenotype value cannot reach the gene-set construction even accidentally.
  m <- m[, !names(m) %in% c("StressGroup", "ExpGroup", "Group", "contrast",
                            "Contrast", "Direction", "Condition"), drop = FALSE]

  # Measured-proteome background: every gene measured in THIS dataset.
  background <- unique(sym_of(m))
  background <- background[!is.na(background) & nzchar(background)]
  message(sprintf("%-16s %d modules, measured background %d genes",
                  ds, length(unique(m$ModuleID)), length(background)))

  for (scope in ewce_module_scopes()) {
    sets <- build_scope_sets(m, scope)
    sets <- sets[vapply(sets, function(g) length(unique(g)) > 0L, logical(1))]
    if (!length(sets)) next
    for (lvl in LEVELS) {
      message(sprintf("  scope=%-10s level=%d  (%d module sets)",
                      scope, lvl, length(sets)))
      res <- run_ewce_gene_set_annotation(
        gene_sets = sets,
        background = background,
        reference = ctd,
        dataset = ds,
        celltype_level = lvl,
        n_boot = REPS,
        seed = SEED,
        fdr_family = ewce_fdr_family_module_annotation(ds, scope, lvl),
        provenance = list(module_scope = scope, reps = REPS,
                          membership_source = relative_to(p)))
      res$module_scope <- scope
      all_rows[[length(all_rows) + 1L]] <- res
    }
  }
}

long <- dplyr::bind_rows(all_rows)
long <- ewce_apply_family_fdr(long)
long <- long %>%
  dplyr::rename(ModuleID = "gene_set_id") %>%
  dplyr::select("dataset", "ModuleID", "module_scope", "level", "cell_type",
                "n_input_genes", "n_mapped_genes", "n_background",
                "observed_statistic", "null_mean", "null_sd", "z_score",
                "fold_change", "p_value_raw_ewce", "p_value", "p_correction",
                "min_attainable_p", "FDR", "fdr_family",
                "annotation_status", "provenance")

# --------------------------------------------------------------- summary

tested <- long[long$annotation_status == "tested" & is.finite(long$FDR), , drop = FALSE]
best <- tested %>%
  dplyr::filter(.data$z_score > 0) %>%
  dplyr::group_by(.data$dataset, .data$ModuleID, .data$module_scope, .data$level) %>%
  dplyr::slice_min(.data$FDR, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup()

summary_tbl <- best %>%
  dplyr::group_by(.data$dataset, .data$ModuleID, .data$level) %>%
  dplyr::summarise(
    strongest_cell_type = .data$cell_type[which.min(.data$FDR)],
    best_FDR = min(.data$FDR, na.rm = TRUE),
    best_z = .data$z_score[which.min(.data$FDR)],
    scope_of_best = .data$module_scope[which.min(.data$FDR)],
    n_scopes_tested = dplyr::n(),
    scopes_agree = dplyr::n_distinct(.data$cell_type) == 1L,
    cell_types_by_scope = paste(paste0(.data$module_scope, ":", .data$cell_type),
                                collapse = "; "),
    any_scope_significant = any(.data$FDR < 0.05, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(
    ambiguity_flag = dplyr::case_when(
      !.data$any_scope_significant ~ "no_significant_affinity",
      !.data$scopes_agree ~ "mixed_scopes_disagree",
      .data$n_scopes_tested < length(ewce_module_scopes()) ~ "partial_scope_coverage",
      TRUE ~ "consistent"),
    annotation_note = paste0(
      "External cell-type expression enrichment against a reference dataset. ",
      "NOT empirical compartment affinity and NOT a cell-proportion estimate."))

root <- OUT()
write_csv_safe(long, file.path(root, "WGCNA_module_external_celltype_affinity_long.csv"))
write_csv_safe(summary_tbl, file.path(root, "WGCNA_module_external_celltype_affinity_summary.csv"))

cat("\n===== WGCNA module external cell-type annotation =====\n")
cat(sprintf("rows: %d   modules x scope x level tested: %d   FDR families: %d\n",
            nrow(long), nrow(best), length(unique(long$fdr_family))))
cat("reps:", REPS, " levels:", paste(LEVELS, collapse = ","), " seed:", SEED, "\n")
for (ds in DATASETS) {
  s <- summary_tbl[summary_tbl$dataset == ds & summary_tbl$level == 1L, , drop = FALSE]
  if (!nrow(s)) next
  cat(sprintf("\n--- %s (level 1) ---\n", ds))
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-12s %-28s FDR=%-9.2e %-22s %s\n", s$ModuleID[i],
                substr(s$strongest_cell_type[i], 1, 28), s$best_FDR[i],
                s$ambiguity_flag[i], s$scope_of_best[i]))
  }
}
cat("\nambiguity flags:\n"); print(table(summary_tbl$ambiguity_flag))
cat("\nOutputs:", relative_to(root), "\n")
cat("Phenotype-blind: no StressGroup, contrast or DAP status was read.\n")
cat("WGCNA labels are unchanged by this script.\n")
