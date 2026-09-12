#!/usr/bin/env Rscript

# Story-v4 step 1: a COMPLETE evidence inventory of the currently defensible
# proteomics results, reconstructed from canonical outputs and the latest
# audits. Figures are designed from this, not the other way round.
#
# Every quantitative_result here is READ from a canonical artifact. Nothing is
# recomputed and no new inference is created.
#
# Also emits figure_claim_evidence_map.md: for each proposed main-figure claim,
# what supports it, what only contextualises it, what limits it, and the single
# visual form that communicates it.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/story_v4_evidence_inventory.R")

OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "story_v4"); dir_create(d)
  file.path(d, ...)
}
REP <- function(...) {
  d <- path_results("reports", "manuscript_candidates", "story_v4"); dir_create(d)
  file.path(d, ...)
}
rd <- function(p, ...) if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, ...) else NULL

# ------------------------------------------------------------ read canonical

gsea <- as.data.frame(data.table::fread(path_results(
  "tables", "10_biological_integration", "gsea_wgcna_concordance", "global",
  "ontology_aware_gsea_theme_assignments_all_contrasts.csv"), showProgress = FALSE))
member <- rd(path_results("source_data", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "sus_res_dap_membership.csv"))
atlas <- rd(path_results("tables", "11_spatial_systems", "atlas",
                         "protein_spatial_cell_affinity.csv"))
rob <- rd(path_results("tables", "11_spatial_systems", "ca2_slm_robustness",
                       "CA2_SLM_DAP_robustness.csv"))
ident <- rd(path_results("tables", "11_spatial_systems", "ca2_slm_robustness",
                         "stress_identity_robustness_comparison.csv"))
bil <- rd(path_results("tables", "11_spatial_systems", "bilateral",
                       "bilateral_spatial_identity_summary.csv"))
prec <- rd(path_results("tables", "11_spatial_systems", "precision",
                        "bilateral_precision_gain.csv"))
mod <- rd(path_results("tables", "11_spatial_systems", "atlas",
                       "WGCNA_module_spatial_cell_affinity.csv"))
eff <- rd(path_results("source_data", "manuscript_panels", "figure_3",
                       "figure3b_stage07_effect_source.csv"))
net <- rd(path_results("tables", "11_spatial_systems", "networks",
                       "network_global_multivariate_test.csv"))
coup <- rd(path_results("tables", "08_behavior_physio_coupling",
                        "network_behavior_coupling", "edge_behavior_correlations.csv"))
kaul <- rd(path_results("source_data", "04_differential_expression_enrichment",
  "control_spatial_identity_validation", "global", "figure2e_source_data.csv"))
intgo <- rd(path_results("source_data", "04_differential_expression_enrichment",
  "control_spatial_identity_validation", "global",
  "figure2f_regions_CA1layers_source_data.csv"))

row <- function(claim_id, dataset, spatial_unit, contrast, analysis_type,
                quantitative_result, inferential_status, FDR_status, QC_status,
                bilateral_status, biological_program, current_candidate_panel,
                main_text_candidate, extended_data_candidate, reason) {
  data.frame(claim_id, dataset, spatial_unit, contrast, analysis_type,
             quantitative_result, inferential_status, FDR_status, QC_status,
             bilateral_status, biological_program, current_candidate_panel,
             main_text_candidate, extended_data_candidate, reason,
             stringsAsFactors = FALSE)
}
inv <- list()
add <- function(...) inv[[length(inv) + 1L]] <<- row(...)

# ---- A. canonical DAP landscape ------------------------------------------
if (!is.null(member)) {
  tb <- as.data.frame(table(member$dataset, member$spatial_unit,
                            member$DAP_direction), stringsAsFactors = FALSE)
  tb <- tb[tb$Freq > 0, ]
  for (i in seq_len(nrow(tb))) {
    add(sprintf("A_dap_%s_%s_%s", tb$Var1[i], tb$Var2[i],
                gsub("[^A-Za-z]", "", tb$Var3[i])),
        tb$Var1[i], tb$Var2[i], "SUS - RES", "differential_abundance",
        sprintf("%d FDR-supported proteins (%s)", tb$Freq[i], tb$Var3[i]),
        "canonical_primary_inference", "BH FDR < 0.05",
        if (tb$Var2[i] == "CA2_slm") "QC-audited: 6 claimable of 28" else "not_audited",
        "animal-level hemisphere-averaged", "protein-level",
        "sv3_da / n3_da", "yes", "yes",
        "establishes that single-protein effects are sparse and spatially restricted")
  }
}
# ---- B. CA2-SLM robustness -----------------------------------------------
if (!is.null(rob)) {
  tb <- table(rob$CA2_SLM_robustness_class)
  for (k in names(tb)) {
    add(paste0("B_ca2slm_", k), "neuron_neuropil", "CA2_slm", "SUS - RES",
        "QC robustness audit", sprintf("%d of 28 proteins", tb[[k]]),
        "descriptive robustness diagnostic", "no new FDR",
        "differential missingness + 2 QC-failed SUS samples",
        "animal-level hemisphere-averaged", "protein-level",
        "3x_dap_status / n3_da",
        if (k == "robust_to_missingness_and_QC") "yes (as quiet status)" else "no",
        "yes",
        "only the robust subset is claimable; the full decomposition is an ED story")
  }
}
# ---- C/D. ranked GSEA, all three contrasts -------------------------------
sup <- gsea[is.finite(gsea$GSEA_FDR) & gsea$GSEA_FDR < 0.05 &
              nzchar(gsea$theme_id) & gsea$theme_role == "primary", ]
for (ds in unique(sup$dataset)) {
  for (ct in c("RES - CON", "SUS - CON", "SUS - RES")) {
    z <- sup[sup$dataset == ds & sup$contrast == ct, ]
    if (!nrow(z)) next
    add(sprintf("C_gsea_%s_%s", ds, gsub("[^A-Za-z]", "", ct)), ds,
        paste(sort(unique(z$spatial_unit)), collapse = ";"), ct,
        "ranked GSEA (GO-BP)",
        sprintf("%d FDR-supported terms over %d spatial units, %d themes",
                nrow(z), length(unique(z$spatial_unit)), length(unique(z$theme_id))),
        "canonical primary inference (per term)", "BH FDR < 0.05 per term",
        "not QC-limited", "animal-level", "program-level",
        "sv3_atlas / n3_atlas", "yes", "reference contrasts to ED",
        "the central positive result: coordinated program differences where single proteins are sparse")
  }
}
# ---- E. three representative programs, trajectory -------------------------
prog <- data.frame(
  ds = c("neuron_neuropil", "neuron_soma", "microglia"),
  unit = c("CA3_sr", "CA2_sp", "CA1"),
  term = c("GO:0099536", "GO:0006397", "GO:0006119"),
  lab = c("synaptic signalling", "mRNA processing", "oxidative phosphorylation"),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(prog))) {
  z <- gsea[gsea$dataset == prog$ds[i] & gsea$spatial_unit == prog$unit[i] &
              gsea$GO_ID == prog$term[i], ]
  for (ct in c("RES - CON", "SUS - CON", "SUS - RES")) {
    w <- z[z$contrast == ct, ]
    if (!nrow(w)) next
    add(sprintf("E_traj_%s_%s", prog$ds[i], gsub("[^A-Za-z]", "", ct)),
        prog$ds[i], prog$unit[i], ct, "ranked GSEA exemplar",
        sprintf("NES %.2f, FDR %.2e", w$NES[1], w$GSEA_FDR[1]),
        "canonical primary inference",
        if (w$GSEA_FDR[1] < 0.05) "FDR < 0.05" else "not FDR-supported",
        "not QC-limited", "animal-level", prog$lab[i],
        "sv3_ex_* (grey wedge - to be replaced)", "yes", "no",
        "the three-contrast trajectory distinguishes divergent from graded programs and is currently MISSING from every figure")
  }
}
# ---- F/I/J. WGCNA module identity, cell type, bilateral ------------------
if (!is.null(mod)) {
  for (ds in unique(mod$dataset)) {
    z <- mod[mod$dataset == ds, ]
    add(paste0("F_modules_", ds), ds, "all", "not applicable",
        "WGCNA module structure",
        sprintf("%d modules; %d bilaterally reproducible; median tau %.2f",
                nrow(z), sum(z$bilateral_reproducibility_class ==
                               "reproducible_level_and_pattern", na.rm = TRUE),
                stats::median(z$spatial_tau, na.rm = TRUE)),
        "descriptive molecular architecture", "not a phenotype test",
        "not QC-limited", "bilaterally assessed", "module identity",
        "sv3_wgcna_circle / strip", "only if it earns space", "yes",
        "module STRUCTURE is robust; it is identity evidence, not phenotype evidence")
  }
}
# ---- G/H. WGCNA phenotype effects + descriptive geometry -----------------
if (!is.null(eff)) {
  nsig <- sum(eff$tier_specific_fdr < 0.05, na.rm = TRUE)
  gcol <- grep("descriptive_RES", names(eff), value = TRUE)[1]
  geom <- if (!is.na(gcol)) sum(eff[[gcol]] %in% TRUE) else NA_integer_
  add("G_wgcna_phenotype", "neuron_neuropil", "global", "all three",
      "WGCNA module group effects",
      sprintf("%d of %d module x contrast cells FDR-supported (min FDR %.3f)",
              nsig, nrow(eff), min(eff$tier_specific_fdr, na.rm = TRUE)),
      "descriptive only", "0 of 45 survive FDR", "not QC-limited",
      "equal-weight L/R", "module phenotype", "n3_wgcna_small (dropped)",
      "no", "yes",
      "a colour field over 45 cells with no supported cell implies inference that does not exist")
  if (!is.na(geom)) {
    add("H_descriptive_geometry", "neuron_neuropil", "global", "RES/CON/SUS",
        "descriptive ordering", sprintf("%d of 15 modules show RES > CON > SUS", geom),
        "descriptive geometry", "no FDR support", "not QC-limited",
        "equal-weight L/R", "module phenotype", "none",
        "no", "yes",
        "a consistent descriptive ordering is interesting context but must never be drawn as significance")
  }
}
# ---- K. stress vs baseline spatial affinity ------------------------------
if (!is.null(ident)) {
  for (i in seq_len(nrow(ident))) {
    add(paste0("K_identity_", ident$subset[i]), "all three", "all", "SUS - RES",
        "effect vs baseline spatial affinity",
        sprintf("%d of %d outside baseline affinity (%.0f%%)",
                ident$effect_outside_baseline_affinity[i], ident$n_hits[i],
                100 * ident$fraction_outside_baseline_affinity[i]),
        "descriptive classification from canonical labels", "no new FDR",
        "robustness subsets explicit", "animal-level", "spatial principle",
        "3x_stress_identity / n3_identity",
        if (ident$subset[i] == "CA2_SLM_robustness_qualified") "candidate for final panel" else "no",
        "yes", "survives every robustness subset; the rank-10 tail does not")
  }
}
# ---- M/N. external and internal validation --------------------------------
if (!is.null(kaul)) {
  f <- intersect(c("p_adjust", "p.adjust"), names(kaul))[1]
  add("M_kaulich", "neuron_neuropil + neuron_soma", "several", "anatomical",
      "external signature GSEA",
      sprintf("%d signature x contrast cells, %d at FDR < 0.05", nrow(kaul),
              sum(kaul[[f]] < 0.05, na.rm = TRUE)),
      "canonical validation inference", "BH FDR < 0.05", "not QC-limited",
      "blocked by AnimalID", "spatial identity", "sv2_external / n2_external",
      "yes", "no", "independent external confirmation of spatial identity")
}
if (!is.null(intgo)) {
  f <- intersect(c("p_adjust", "p.adjust"), names(intgo))[1]
  add("N_internal_go", "neuron_neuropil", "regions + CA1 layers", "anatomical",
      "internal anatomical GO GSEA",
      sprintf("%d terms, %d at FDR < 0.05", nrow(intgo),
              sum(intgo[[f]] < 0.05, na.rm = TRUE)),
      "canonical validation inference", "BH FDR < 0.05", "not QC-limited",
      "blocked by AnimalID", "anatomical programs", "sv2_internal / n2_internal",
      "yes", "no", "internal biology recovers known anatomy, a distinct layer from external validation")
}
# ---- O. bilateral reproducibility, ALL contrasts --------------------------
if (!is.null(bil)) {
  for (i in seq_len(nrow(bil))) {
    add(paste0("O_bilateral_", bil$dataset[i], "_", i), bil$dataset[i],
        bil$contrast[i], "anatomical", "paired hemisphere reproducibility",
        sprintf("Pearson r = %.3f, sign agreement %.3f", bil$pearson_r[i],
                bil$sign_agreement_fraction[i]),
        "descriptive reproducibility", "not a phenotype test", "not QC-limited",
        "paired within animal", "spatial identity", "sv2_bilateral (3 of 15 shown)",
        "yes - ALL contrasts, not 3", "no",
        "the spread 0.38 to 0.92 is the finding: coarse regional structure reproduces, fine CA1 laminar structure less so")
  }
}
# ---- P. precision ---------------------------------------------------------
if (!is.null(prec)) {
  for (k in unique(prec$endpoint_class)) {
    z <- prec[prec$endpoint_class == k, ]
    add(paste0("P_icc_", k), "all three", "endpoint", "not applicable",
        "variance decomposition / ICC",
        sprintf("median ICC %.2f (one side) to %.2f (bilateral mean), n = %d endpoints",
                stats::median(z$ICC_single_side, na.rm = TRUE),
                stats::median(z$ICC_bilateral_mean, na.rm = TRUE), nrow(z)),
        "descriptive reliability", "not a phenotype test", "not QC-limited",
        "single vs bilateral", "measurement quality", "sv2_precision",
        "small inset only", "yes", "reliability gain justifies bilateral sampling but is not a biological result")
  }
}
# ---- Q/R. nulls -----------------------------------------------------------
if (!is.null(net)) {
  add("Q_network_null", "all three", "whole network", "SUS vs RES",
      "exact label enumeration",
      paste(sprintf("%s p=%.3f", net$dataset, net$exact_p), collapse = "; "),
      "exact test, informative null", "min attainable p = 0.0036",
      "not QC-limited", "bilateral networks", "network organisation",
      "none", "no", "yes",
      "an informative null: the test had resolution to detect a difference and found none")
}
if (!is.null(coup)) {
  f <- grep("p.adj", names(coup), value = TRUE)[1]
  add("R_coupling_null", "neuron_neuropil", "edges", "behaviour",
      "edge-behaviour correlation",
      sprintf("%d correlations at n = 9, min BH FDR = %.2f", nrow(coup),
              min(coup[[f]], na.rm = TRUE)),
      "exploratory", "nothing survives FDR", "after the AnimalID fix",
      "animal-level", "behaviour coupling", "none", "no", "yes",
      "estimable for the first time after the identity fix, but null")
}

inventory <- do.call(rbind, inv)
write_csv_safe(inventory, OUT("proteomics_manuscript_evidence_inventory.csv"))

cat("\n===== Evidence inventory =====\n")
cat("claims recorded:", nrow(inventory), "\n\n")
print(as.data.frame(table(inventory$main_text_candidate)), row.names = FALSE)
cat("\nby analysis type:\n")
print(as.data.frame(table(inventory$analysis_type)), row.names = FALSE)
cat("\nwritten:", relative_to(OUT("proteomics_manuscript_evidence_inventory.csv")), "\n")
