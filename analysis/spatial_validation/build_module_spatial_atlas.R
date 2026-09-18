#!/usr/bin/env Rscript
#
#
# CONSUMES validated canonical outputs; recomputes none of them. In particular
# the hemisphere contract, bilateral validations, variance decomposition and
# module EWCE annotation are read, not rebuilt.
#
# FOUR EVIDENCE DIMENSIONS, KEPT SEPARATE
#   A spatial anatomical identity    (CON baseline only)
#   B empirical compartment affinity (this experiment's ROIs)
#   C external cell-type affinity    (outside reference, via EWCE)
#   D bilateral reliability
# There is no composite score. A categorical context class is derived from
# explicit rules and every row carries the rule that produced it.
#
# PHENOTYPE-BLIND. No SUS/RES value enters any affinity or identity calculation.
# WGCNA labels are audited against the atlas but never changed.
#
# USAGE
#   Rscript analysis/spatial_validation/build_module_spatial_atlas.R
#   Rscript analysis/spatial_validation/build_module_spatial_atlas.R --dry-run
# Script: analysis/spatial_validation/build_module_spatial_atlas.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_hemisphere_values.csv; results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; results/spatial_validation/annotate_module_celltypes/global/tables/WGCNA_module_external_celltype_affinity_long.csv; +1 more; optional results/spatial_validation/quantify_module_bilateral_identity/global/tables/WGCNA_module_bilateral_reproducibility.csv; results/tables/11_spatial_systems/bilateral/WGCNA_module_bilateral_reproducibility.csv; results/spatial_validation/decompose_bilateral_variance/global/tables/bilateral_precision_gain.csv; +2 more
# Produces: results/spatial_validation/build_module_spatial_atlas/global/tables/WGCNA_module_spatial_cell_affinity.csv; results/spatial_validation/build_module_spatial_atlas/global/tables/WGCNA_module_baseline_spatial_profile_long.csv; results/spatial_validation/build_module_spatial_atlas/global/tables/WGCNA_module_spatial_fingerprints_raw.csv; +7 more
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Integrated spatial / cell-affinity atlas for WGCNA modules.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/spatial/spatial_atlas_utils.R")
source("R/enrichment/ewce_gene_set_engine.R")
source(repo_path("R", "spatial_systems_paths.R"))
source(repo_path("R", "qc_result_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "build_module_spatial_atlas"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/spatial_validation/build_module_spatial_atlas.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- CANONICAL_PATHS$tables; dir_create(d); d
}
P <- list(
  membership = function(ds) path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds,
                                         "modules", "WGCNA_modules_long.csv"),
  universe   = function(ds) path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds,
                                         "modules", "WGCNA_feature_universe.csv"),
  contract   = function(ds) path_results("tables", "06_modules_WGCNA", "identity_contract",
                                         ds, "WGCNA_module_supermodule_membership_contract.csv"),
  lookup     = function(ds) path_results("tables", "06_modules_WGCNA", "interpretable_summary",
                                         ds, "WGCNA_final_label_lookup.csv"),
  bilateral  = spatial_systems_find("WGCNA_module_bilateral_reproducibility.csv",
                                    "bilateral"),
  bil_prof   = spatial_systems_find("WGCNA_module_bilateral_spatial_profile_summary.csv",
                                    "bilateral"),
  precision  = spatial_systems_find("bilateral_precision_gain.csv", "precision"),
  ewce_long  = spatial_systems_find("WGCNA_module_external_celltype_affinity_long.csv",
                                    "celltype_annotation"),
  emp_sets   = qc_find("empirical_roi_marker_sets.csv",
                       owner = "discover_empirical_roi_markers",
                       legacy_substep = "05_empirical_roi_marker_discovery"),
  emp_bil    = spatial_systems_find("bilateral_empirical_compartment_summary.csv",
                                    "bilateral"),
  emp_trans  = spatial_systems_find("bilateral_empirical_marker_transfer.csv",
                                    "bilateral"),
  ref_panels = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
)
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- list(ewce_module_annotation = P$ewce_long,
                 module_bilateral = P$bilateral,
                 empirical_marker_sets = P$emp_sets,
                 reference_marker_panels = P$ref_panels)
  for (ds in DATASETS) inputs[[paste0("membership__", ds)]] <- P$membership(ds)
  cat("[DRY-RUN] Module spatial/cell-affinity atlas.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Consumes validated outputs; recomputes no WGCNA, DA or EWCE.\n")
  quit(save = "no", status = 0L)
}

read_req <- function(p, lbl) {
  if (!file.exists(p)) stop("missing_required_input: ", lbl, ": ", p, call. = FALSE)
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}
read_opt <- function(p) if (file.exists(p)) {
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
} else NULL

# ===================================================== PART 2 + 11: baseline

message("Building CON baseline module spatial profiles")
baseline_long <- dplyr::bind_rows(lapply(DATASETS, sat_con_bilateral_module_profiles))
prof_metrics <- list(); raw_rows <- list(); z_rows <- list(); sim_rows <- list()

for (ds in DATASETS) {
  b <- baseline_long[baseline_long$dataset == ds, , drop = FALSE]
  m <- sat_profile_matrix(b)
  # peak/second/ranks come from the RAW matrix, never the z-scored one
  prof_metrics[[ds]] <- dplyr::bind_rows(lapply(rownames(m), function(r)
    cbind(data.frame(dataset = ds, ModuleID = r, stringsAsFactors = FALSE),
          sat_profile_metrics(m[r, ], colnames(m)))))
  raw_rows[[ds]] <- cbind(data.frame(dataset = ds, ModuleID = rownames(m),
                                     stringsAsFactors = FALSE),
                          as.data.frame(m, check.names = FALSE))
  z <- sat_row_z(m)
  z_rows[[ds]] <- cbind(data.frame(dataset = ds, ModuleID = rownames(z),
                                   stringsAsFactors = FALSE),
                        as.data.frame(z, check.names = FALSE))
  # PART 11: module-module spatial-profile correlation, within dataset
  if (nrow(m) >= 2L) {
    cm <- suppressWarnings(stats::cor(t(m), use = "pairwise.complete.obs"))
    idx <- which(upper.tri(cm), arr.ind = TRUE)
    sim_rows[[ds]] <- data.frame(
      dataset = ds, module_a = rownames(cm)[idx[, 1]],
      module_b = colnames(cm)[idx[, 2]],
      spatial_profile_pearson = cm[idx],
      spatial_distance = 1 - cm[idx], stringsAsFactors = FALSE)
  }
}
prof_metrics <- dplyr::bind_rows(prof_metrics)
similarity <- dplyr::bind_rows(sim_rows)

# ===================================== PART 4 + 6: enrichment against markers

message("Testing empirical compartment and reference marker affinity")
emp_sets <- read_req(P$emp_sets, "empirical marker sets")
ref_panels <- read_req(P$ref_panels, "reference marker panels")

SCOPES <- c("all", "core_kME06", "top25")
scope_filter <- function(m, scope) switch(
  scope,
  all = rep(TRUE, nrow(m)),
  core_kME06 = m[["is_core_kME_0.6"]] %in% TRUE,
  top25 = m[["is_top_hub_25"]] %in% TRUE)

emp_rows <- list(); ref_rows <- list()
for (ds in DATASETS) {
  memb <- read_req(P$membership(ds), "WGCNA membership")
  univ <- read_req(P$universe(ds), "measured feature universe")
  universe_ids <- unique(as.character(univ$ProteinGroupID))

  # empirical compartment marker sets, mapped to THIS dataset's protein groups
  id_col <- intersect(c(paste0("ProteinGroupID_", ds), "ProteinGroupID"),
                      names(emp_sets))
  # reference panels are gene-symbol based; map through the membership table
  sym <- toupper(trimws(as.character(memb$GeneSymbol)))
  sym_lut <- data.frame(sym = sym, pg = as.character(memb$ProteinGroupID),
                        stringsAsFactors = FALSE)
  sym_lut <- sym_lut[nzchar(sym_lut$sym) & !is.na(sym_lut$sym), , drop = FALSE]

  for (scope in SCOPES) {
    keep <- scope_filter(memb, scope)
    mod_by_id <- split(as.character(memb$ProteinGroupID[keep]),
                       as.character(memb$ModuleID[keep]))

    for (mid in names(mod_by_id)) {
      for (set in sort(unique(emp_sets$marker_set))) {
        mk <- unique(unlist(lapply(id_col, function(cc)
          as.character(emp_sets[[cc]][emp_sets$marker_set == set]))))
        mk <- mk[!is.na(mk) & nzchar(mk)]
        if (!length(mk)) next
        e <- sat_fisher_enrichment(mod_by_id[[mid]], mk, universe_ids)
        emp_rows[[length(emp_rows) + 1L]] <- cbind(
          data.frame(dataset = ds, ModuleID = mid, module_scope = scope,
                     empirical_compartment = set,
                     fdr_family = sat_fdr_family_empirical_compartment(ds, scope),
                     universe_definition = "measured proteome (WGCNA feature universe)",
                     stringsAsFactors = FALSE), e)
      }
      for (panel in sort(unique(ref_panels$marker_set))) {
        syms <- toupper(trimws(ref_panels$gene_symbol[ref_panels$marker_set == panel]))
        pg <- unique(sym_lut$pg[sym_lut$sym %in% syms])
        if (length(pg) < 3L) next
        e <- sat_fisher_enrichment(mod_by_id[[mid]], pg, universe_ids)
        ref_rows[[length(ref_rows) + 1L]] <- cbind(
          data.frame(dataset = ds, ModuleID = mid, module_scope = scope,
                     reference_panel = panel,
                     fdr_family = sat_fdr_family_reference_marker(ds, scope),
                     universe_definition = "measured proteome (WGCNA feature universe)",
                     stringsAsFactors = FALSE), e)
      }
    }
  }
}
empirical <- sat_apply_family_fdr(dplyr::bind_rows(emp_rows))
reference <- sat_apply_family_fdr(dplyr::bind_rows(ref_rows))

# Flag which marker sets can discriminate at all. Nothing is dropped; the flag
# governs only the single "strongest compartment" call.
for (nm in c("empirical", "reference")) {
  d <- get(nm)
  d$marker_set_universe_fraction <- d$n_markers_in_universe / d$n_universe
  d$informative_for_affinity_call <- sat_marker_set_is_informative(
    d$n_markers_in_universe, d$n_universe)
  d$informativeness_rule <- sprintf(
    "set must cover %.1f%%-%.0f%% of the measured universe and have >= %d members",
    100 * sat_marker_set_bounds()$min_fraction,
    100 * sat_marker_set_bounds()$max_fraction,
    sat_marker_set_bounds()$min_members)
  assign(nm, d)
}

# ===================================================== PART 7: external EWCE

message("Joining validated external cell-type annotation")
ewce <- read_req(P$ewce_long, "module EWCE annotation")
ewce_tested <- ewce[ewce$annotation_status == "tested" & is.finite(ewce$FDR) &
                      ewce$z_score > 0, , drop = FALSE]
best_by_scope <- ewce_tested %>%
  dplyr::group_by(.data$dataset, .data$ModuleID, .data$module_scope, .data$level) %>%
  dplyr::slice_min(.data$FDR, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup()

ext_wide <- best_by_scope %>%
  dplyr::filter(.data$level == 1L) %>%
  dplyr::select("dataset", "ModuleID", "module_scope", "cell_type", "FDR") %>%
  tidyr::pivot_wider(names_from = "module_scope",
                     values_from = c("cell_type", "FDR"))
names(ext_wide) <- sub("^cell_type_", "external_celltype_", names(ext_wide))
names(ext_wide) <- sub("^FDR_", "external_FDR_", names(ext_wide))

ext_l2 <- best_by_scope %>%
  dplyr::filter(.data$level == 2L, .data$module_scope == "all") %>%
  dplyr::select("dataset", "ModuleID",
                external_celltype_level2 = "cell_type",
                external_FDR_level2 = "FDR")

scope_concord <- best_by_scope %>%
  dplyr::filter(.data$level == 1L) %>%
  dplyr::group_by(.data$dataset, .data$ModuleID) %>%
  dplyr::summarise(
    n_scopes_with_result = dplyr::n(),
    external_celltype_scope_concordance = dplyr::case_when(
      dplyr::n_distinct(.data$cell_type) == 1L & dplyr::n() == length(SCOPES) ~ "consistent",
      dplyr::n_distinct(.data$cell_type) > 1L ~ "mixed_scopes_disagree",
      TRUE ~ "partial_scope_coverage"),
    external_any_significant = any(.data$FDR < 0.05, na.rm = TRUE),
    .groups = "drop")

# ===================================================== PART 3: bilateral

message("Joining bilateral reliability")
bil <- read_req(P$bilateral, "module bilateral reproducibility")
bil <- bil[bil$level == "module", , drop = FALSE]
bil_prof <- read_opt(P$bil_prof)
prec <- read_opt(P$precision)

bil_join <- bil %>%
  dplyr::select("dataset", ModuleID = "endpoint_id",
                bilateral_absolute_pearson = "pearson_r",
                bilateral_absolute_spearman = "spearman_rho",
                bilateral_MAE = "MAE",
                bilateral_median_abs_difference = "median_abs_difference",
                bilateral_n_pairs = "n_pairs",
                bilateral_n_animals = "n_animals",
                bilateral_reproducibility_class = "bilateral_reproducibility_class",
                bilateral_spatial_profile_pearson = "median_profile_Pearson")
if (!is.null(prec)) {
  pj <- prec %>%
    dplyr::filter(.data$endpoint_class == "wgcna_module_eigengene") %>%
    dplyr::select("dataset", ModuleID = "endpoint_id",
                  ICC_single_side = "ICC_single_side",
                  ICC_bilateral_mean = "ICC_bilateral_mean",
                  precision_assumption_status = "assumption_status")
  bil_join <- dplyr::left_join(bil_join, pj, by = c("dataset", "ModuleID"))
}
bil_join$bilateral_support_class <- sat_bilateral_support_class(
  bil_join$bilateral_absolute_pearson, bil_join$bilateral_spatial_profile_pearson)

# ===================================================== PART 8: integrate

message("Assembling the module atlas")
identity_rows <- dplyr::bind_rows(lapply(DATASETS, function(ds) {
  memb <- read_req(P$membership(ds), "WGCNA membership")
  ct <- read_opt(P$contract(ds))
  lk <- read_opt(P$lookup(ds))
  size <- memb %>% dplyr::count(.data$ModuleID, name = "module_size")
  hubs <- memb %>%
    dplyr::filter(.data[["is_top_hub_25"]] %in% TRUE) %>%
    dplyr::arrange(.data$ModuleID, dplyr::desc(.data$abs_kME)) %>%
    dplyr::group_by(.data$ModuleID) %>%
    dplyr::summarise(top_hubs = paste(utils::head(.data$GeneSymbol, 10), collapse = ", "),
                     .groups = "drop")
  out <- data.frame(dataset = ds, ModuleID = size$ModuleID,
                    module_size = size$module_size, stringsAsFactors = FALSE) %>%
    dplyr::left_join(hubs, by = "ModuleID")
  if (!is.null(ct)) {
    out <- dplyr::left_join(out, ct %>%
      dplyr::select(ModuleID = "module_id", SupermoduleID = "supermodule_id"),
      by = "ModuleID")
  }
  if (!is.null(lk)) {
    lab <- lk %>% dplyr::filter(.data$level == "module") %>%
      dplyr::select(ModuleID = "entity_id", canonical_display_label = "final_plot_label")
    out <- dplyr::left_join(out, lab, by = "ModuleID")
  }
  out
}))

# strongest empirical compartment per module, taken at the `all` scope so the
# choice is prospective rather than the most attractive across scopes
emp_best <- empirical %>%
  dplyr::filter(.data$module_scope == "all",
                .data$informative_for_affinity_call %in% TRUE) %>%
  dplyr::group_by(.data$dataset, .data$ModuleID) %>%
  dplyr::slice_min(.data$FDR, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select("dataset", "ModuleID",
                strongest_empirical_compartment = "empirical_compartment",
                empirical_compartment_FDR = "FDR",
                empirical_compartment_odds_ratio = "odds_ratio",
                empirical_compartment_overlap = "n_overlap")
emp_scope_support <- empirical %>%
  dplyr::filter(.data$informative_for_affinity_call %in% TRUE) %>%
  dplyr::group_by(.data$dataset, .data$ModuleID, .data$module_scope) %>%
  dplyr::summarise(sig = any(.data$FDR < 0.05, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = "module_scope", values_from = "sig",
                     names_prefix = "empirical_significant_")

ref_best <- reference %>%
  dplyr::filter(.data$module_scope == "all",
                .data$informative_for_affinity_call %in% TRUE) %>%
  dplyr::group_by(.data$dataset, .data$ModuleID) %>%
  dplyr::slice_min(.data$FDR, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select("dataset", "ModuleID",
                strongest_reference_panel = "reference_panel",
                reference_panel_FDR = "FDR",
                reference_panel_overlap = "n_overlap")

atlas <- identity_rows %>%
  dplyr::left_join(prof_metrics, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(bil_join, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(emp_best, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(emp_scope_support, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(ref_best, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(ext_wide, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(ext_l2, by = c("dataset", "ModuleID")) %>%
  dplyr::left_join(scope_concord, by = c("dataset", "ModuleID"))

# PRESPECIFIED evidence thresholds for the rule-based class
SPATIAL_TAU_MIN <- 0.60
atlas$has_spatial_identity <- atlas$spatial_tau >= SPATIAL_TAU_MIN
atlas$has_empirical_affinity <- atlas$empirical_compartment_FDR < 0.05
atlas$has_external_affinity <- atlas$external_any_significant %in% TRUE
ctx <- sat_context_class(
  atlas$has_spatial_identity, atlas$has_empirical_affinity,
  atlas$has_external_affinity, atlas$bilateral_support_class,
  external_mixed = atlas$external_celltype_scope_concordance %in% "mixed_scopes_disagree")
atlas <- cbind(atlas, ctx)

atlas$spatial_context_summary <- ifelse(
  atlas$has_spatial_identity %in% TRUE,
  sprintf("peak %s (tau=%.2f, peak-vs-rest=%.3f)", atlas$peak_unit,
          atlas$spatial_tau, atlas$peak_minus_rest),
  sprintf("no dominant unit (tau=%.2f)", atlas$spatial_tau))
atlas$cell_context_summary <- sprintf(
  "empirical: %s (FDR=%.2g); external L1: %s",
  ifelse(is.na(atlas$strongest_empirical_compartment), "none",
         atlas$strongest_empirical_compartment),
  atlas$empirical_compartment_FDR,
  ifelse(is.na(atlas$external_celltype_all), "none", atlas$external_celltype_all))
atlas$caveats <- ifelse(
  atlas$bilateral_support_class %in% c("low_bilateral_support",
                                       "insufficient_bilateral_data"),
  "spatial_identity_present_but_bilaterally_variable", NA_character_)
atlas$contract_version <- sat_contract_version()

# ===================================================== PART 9: label audit

label_audit <- atlas %>%
  dplyr::select("dataset", "ModuleID", "canonical_display_label",
                "spatial_context_summary", "strongest_empirical_compartment",
                "empirical_compartment_FDR", "external_celltype_all",
                "external_celltype_scope_concordance", "context_confidence")
rel <- sat_label_context_relationship(
  label_audit$canonical_display_label, label_audit$spatial_context_summary,
  label_audit$strongest_empirical_compartment, label_audit$external_celltype_all)
label_audit <- cbind(label_audit, rel)

# The ACTIVE canonical label and the PROPOSED label from the approval table can
# differ, and the external evidence may side with one against the other. That is
# precisely the case a reviewer needs to see, so both are audited.
approval <- read_opt(path_results("reviewer_audit", "wgcna_label_approval",
                                  "WGCNA_final_label_approval_table.csv"))
if (!is.null(approval)) {
  prop <- approval %>%
    dplyr::filter(.data$level == "module") %>%
    dplyr::select("dataset", ModuleID = "entity_id",
                  proposed_final_label = "proposed_final_label",
                  proposed_label_confidence = "confidence",
                  proposed_recommended_for_activation = "recommended_for_activation")
  label_audit <- dplyr::left_join(label_audit, prop, by = c("dataset", "ModuleID"))
  rel_prop <- sat_label_context_relationship(
    label_audit$proposed_final_label, label_audit$spatial_context_summary,
    label_audit$strongest_empirical_compartment, label_audit$external_celltype_all)
  label_audit$proposed_label_context_relationship <- rel_prop$label_context_relationship
  label_audit$proposed_label_context_rule <- rel_prop$label_context_rule
  # where the external evidence backs the proposal but not the active label,
  # that is independent support for a relabel a human has not yet approved
  label_audit$external_supports_proposed_over_active <-
    label_audit$proposed_label_context_relationship == "context_corroborates_label" &
    label_audit$label_context_relationship != "context_corroborates_label"
}
label_audit$audit_note <- "AUDIT ONLY: no module was renamed and no proposed label was activated."

# ===================================================== PART 10 + 12: supermodule

message("Summarising supermodule context")
sm_rows <- list(); coh_rows <- list()
for (ds in DATASETS) {
  a <- atlas[atlas$dataset == ds & !is.na(atlas$SupermoduleID), , drop = FALSE]
  if (!nrow(a)) next
  sim <- similarity[similarity$dataset == ds, , drop = FALSE]
  for (sm in sort(unique(a$SupermoduleID))) {
    mem <- a[a$SupermoduleID == sm, , drop = FALSE]
    within <- sim[sim$module_a %in% mem$ModuleID & sim$module_b %in% mem$ModuleID, ]
    between <- sim[xor(sim$module_a %in% mem$ModuleID, sim$module_b %in% mem$ModuleID), ]
    # MEMBER MODULES ARE WEIGHTED EQUALLY: a large module must not dominate.
    sm_rows[[length(sm_rows) + 1L]] <- data.frame(
      dataset = ds, SupermoduleID = sm, n_member_modules = nrow(mem),
      member_modules = paste(mem$ModuleID, collapse = ";"),
      member_spatial_peaks = paste(mem$peak_unit, collapse = ";"),
      n_distinct_peaks = length(unique(mem$peak_unit)),
      mean_member_tau = mean(mem$spatial_tau, na.rm = TRUE),
      member_empirical_affinities = paste(mem$strongest_empirical_compartment, collapse = ";"),
      member_external_celltypes = paste(mem$external_celltype_all, collapse = ";"),
      n_distinct_external_celltypes = length(unique(stats::na.omit(mem$external_celltype_all))),
      mean_bilateral_absolute_pearson = mean(mem$bilateral_absolute_pearson, na.rm = TRUE),
      bilateral_classes = paste(unique(mem$bilateral_support_class), collapse = ";"),
      weighting = "member modules weighted equally; module size does not influence any summary",
      stringsAsFactors = FALSE)
    coh_rows[[length(coh_rows) + 1L]] <- data.frame(
      dataset = ds, SupermoduleID = sm, n_member_modules = nrow(mem),
      median_within_supermodule_similarity = if (nrow(within)) stats::median(within$spatial_profile_pearson, na.rm = TRUE) else NA_real_,
      median_between_supermodule_similarity = if (nrow(between)) stats::median(between$spatial_profile_pearson, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
supermodule <- dplyr::bind_rows(sm_rows)
coherence <- dplyr::bind_rows(coh_rows)
coherence$spatial_coherence_delta <- coherence$median_within_supermodule_similarity -
  coherence$median_between_supermodule_similarity
coherence$interpretation <- paste0(
  "Descriptive only. Spatial divergence does NOT invalidate a WGCNA supermodule: ",
  "modules are defined by eigengene covariation, not by anatomy.")

# discordant members: module whose profile is farther from its own supermodule
# peers than the dataset median
disc_rows <- list()
for (ds in DATASETS) {
  a <- atlas[atlas$dataset == ds & !is.na(atlas$SupermoduleID), , drop = FALSE]
  sim <- similarity[similarity$dataset == ds, , drop = FALSE]
  if (!nrow(a) || !nrow(sim)) next
  for (i in seq_len(nrow(a))) {
    mid <- a$ModuleID[i]; sm <- a$SupermoduleID[i]
    peers <- setdiff(a$ModuleID[a$SupermoduleID == sm], mid)
    if (!length(peers)) next
    s <- sim[(sim$module_a == mid & sim$module_b %in% peers) |
               (sim$module_b == mid & sim$module_a %in% peers), ]
    disc_rows[[length(disc_rows) + 1L]] <- data.frame(
      dataset = ds, ModuleID = mid, SupermoduleID = sm,
      median_similarity_to_supermodule_peers = stats::median(s$spatial_profile_pearson, na.rm = TRUE),
      n_peers = length(peers), stringsAsFactors = FALSE)
  }
}
discordant <- dplyr::bind_rows(disc_rows)
if (nrow(discordant)) {
  discordant <- discordant %>%
    dplyr::group_by(.data$dataset) %>%
    dplyr::mutate(dataset_median = stats::median(.data$median_similarity_to_supermodule_peers, na.rm = TRUE),
                  spatially_discordant_from_peers =
                    .data$median_similarity_to_supermodule_peers < .data$dataset_median) %>%
    dplyr::ungroup()
}

# ------------------------------------------------------------------ write

root <- OUT()
write_csv_safe(baseline_long, file.path(root, "WGCNA_module_baseline_spatial_profile_long.csv"))
write_csv_safe(dplyr::bind_rows(raw_rows), file.path(root, "WGCNA_module_spatial_fingerprints_raw.csv"))
write_csv_safe(dplyr::bind_rows(z_rows), file.path(root, "WGCNA_module_spatial_fingerprints_z.csv"))
write_csv_safe(similarity, file.path(root, "WGCNA_module_spatial_similarity.csv"))
write_csv_safe(empirical, file.path(root, "WGCNA_module_empirical_compartment_affinity_long.csv"))
write_csv_safe(reference, file.path(root, "WGCNA_module_reference_marker_affinity_long.csv"))
write_csv_safe(atlas, file.path(root, "WGCNA_module_spatial_cell_affinity.csv"))
write_csv_safe(label_audit, file.path(root, "WGCNA_label_spatial_cell_context_audit.csv"))
write_csv_safe(supermodule, file.path(root, "WGCNA_supermodule_spatial_cell_context.csv"))
write_csv_safe(coherence, file.path(root, "WGCNA_supermodule_spatial_coherence.csv"))
if (nrow(discordant)) {
  write_csv_safe(discordant, file.path(root, "WGCNA_module_supermodule_spatial_discordance.csv"))
}

cat("\n===== Module spatial / cell-affinity atlas =====\n")
for (ds in DATASETS) {
  a <- atlas[atlas$dataset == ds, , drop = FALSE]
  if (!nrow(a)) next
  cat(sprintf("\n--- %s (%d modules) ---\n", ds, nrow(a)))
  cat(sprintf("  %-11s %-9s %5s %-22s %-20s %-11s %s\n",
              "module", "peak", "tau", "empirical", "external L1", "bilateral", "context"))
  for (i in seq_len(nrow(a))) {
    cat(sprintf("  %-11s %-9s %5.2f %-22s %-20s %-11s %s\n",
                a$ModuleID[i], substr(a$peak_unit[i], 1, 9), a$spatial_tau[i],
                substr(ifelse(is.na(a$strongest_empirical_compartment[i]), "-",
                              a$strongest_empirical_compartment[i]), 1, 22),
                substr(ifelse(is.na(a$external_celltype_all[i]), "-",
                              a$external_celltype_all[i]), 1, 20),
                sub("_bilateral_support", "", a$bilateral_support_class[i]),
                a$context_confidence[i]))
  }
}
cat("\ncontext classes:\n"); print(table(atlas$context_confidence))
cat("\nlabel-context relationship:\n"); print(table(label_audit$label_context_relationship))
cat("\nOutputs:", relative_to(root), "\n")
cat("Phenotype-blind: no SUS/RES value defined any identity or affinity.\n")
cat("WGCNA labels audited, never changed.\n")
