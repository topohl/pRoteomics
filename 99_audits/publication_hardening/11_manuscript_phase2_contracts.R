#!/usr/bin/env Rscript

# Manuscript phase 2: the contract tables.
#
# Figure 1 could not be drafted. This script records WHY, in a form that can be
# checked, and resolves the three Methods questions that ARE answerable from this
# repository: the WGCNA "0 of 45" scope, the external-validation p_adjust field
# (MT-04), and the PH-009 terminology.
#
# WRITES ONLY manuscript/*.csv. Touches no scientific output.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
MS <- "manuscript"
dir.create(MS, recursive = TRUE, showWarnings = FALSE)

# ================================================ Figure 1 source inventory
#
# The question this table answers is not "where is Figure 1" but "is Figure 1
# constructible from this repository at all". Each row is a thing phase 2 was
# asked to reconstruct, with what was actually found.
F1 <- function(item, status, location, evidence, consequence)
  data.frame(required_item = item, status = status,
             location_in_repo = location, evidence = evidence,
             consequence_for_manuscript = consequence,
             stringsAsFactors = FALSE)
f1 <- rbind(
  F1("Figure 1 panel contract", "ABSENT",
     "figures/figure_final_truth_v9_contract.yml",
     "The frozen contract defines F2_NATURE_FINAL_V9 and F3_NATURE_FINAL_V9 only. No F1 figure key, no F1 panels, no F1 renderers.",
     "There is no frozen Figure 1 to write Results against."),
  F1("Figure 1 output directory", "EMPTY", "results/manuscript/figure_1/",
     "Contains only .gitkeep; no source data, no SVG, no PDF.",
     "No released Figure 1 artefact exists."),
  F1("Behavioural / physiological raw data", "PRESENT (external input)",
     "data/external/behavior/E9_Behavior_Data.xlsx",
     "34 sheets; the zScore sheet carries 117 animals with NOR, sucrose_pref, weight_dev, delta_cort, adrenal_weight, spleen_weight and a precomputed CombZ column.",
     "The data exist but arrive as a finished product."),
  F1("Home-cage / RFID movement AUC", "PRESENT (external input)",
     "data/external/behavior/auc_individual_animals_firstChangeActive.csv (and _all.csv)",
     "322 animals; columns AnimalNum, Batch, Group, metric, Change, Sex, Phase, window, AUC, AUC_norm, AUC_per_hour, prediction_type = subject_specific_gamm_observed_grid.",
     "The GAMM-derived AUC is an INPUT here; the GAMM itself was fitted elsewhere."),
  F1("CombZ / outcome-score construction", "NOT IN REPOSITORY",
     "read at 01_preprocessing/06_merged_metadata_module_score.r:399-427",
     "The script READS sheet 'zScore' and renames comb_z -> CombZ. It performs no z-scoring, no sign orientation, no aggregation, no sex or batch handling. The formula is not present in any tracked file.",
     "Sign orientation, z-scoring population, aggregation formula and the treatment of sex and batch CANNOT be stated from this repository."),
  F1("RES / SUS classification rule", "NOT IN REPOSITORY",
     "01_preprocessing/06_merged_metadata_module_score.r StressGroup = case_when(!is.na(group_auc_all) ~ group_auc_all, ...)",
     "The group label is taken from the Group column already present in the external AUC csv. No threshold, no cut-point and no classifier is computed here.",
     "The RES/SUS boundary definition cannot be stated from this repository."),
  F1("GAMM model specification", "NOT IN REPOSITORY",
     "implied by prediction_type = subject_specific_gamm_observed_grid",
     "No mgcv/gamm4 call exists in any tracked script. Only the fitted per-animal AUC survives.",
     "The GAMM formula, smoothers, random-effect structure and grid cannot be stated."),
  F1("Out-of-sample prediction / cross-validation", "ABSENT",
     "repository-wide search",
     "No glmnet, caret, pROC, randomForest, cv.glmnet, trainControl or createFolds anywhere. Every 'leave-one-animal-out' hit is a PROTEOMIC stability analysis (CA2-SLM robustness, WGCNA readiness, GSEA-WGCNA concordance).",
     "The word 'predicts' cannot be licensed. No AUC, no CV fold structure, no permutation null exists to report."),
  F1("HMM / behavioural state architecture", "ABSENT",
     "repository-wide search",
     "No state model, no state definitions, no state-fraction outputs.",
     "Cannot be described."),
  F1("Experimental timeline (ages, windows)", "NOT RECOVERABLE",
     "data/external/behavior/*",
     "The AUC csv carries Phase (Active) and window (all) and a halfhour grid 0-23, but no age, no date and no cage-change timestamp. The workbook sheets were not searched exhaustively, but no tracked script states ages.",
     "P22/P25/P36/P56-style ages cannot be verified and must not be asserted."),
  F1("Sex-stratified / interaction statistics", "NOT IN REPOSITORY",
     "repository-wide search",
     "Sex is carried as a column in both external inputs. No sex-stratified model, no sex x predictor interaction test is fitted in any tracked script.",
     "No sex claim of any strength can be made."),
  F1("Proteomics-behaviour coupling (the part that IS here)", "PRESENT",
     "08_behavior_physio_coupling/01-04",
     "Correlates the proteome with behavioural summaries. Its outputs already appear in the drafted Results 3 as the edge-behaviour null (0 of 48 tests, 8 neuropil pairs).",
     "Already reported, correctly scoped, in Results 3."))
utils::write.csv(f1, file.path(MS, "figure1_authoritative_source_inventory.csv"),
                 row.names = FALSE)

# ================================================== WGCNA "0 of 45" contract
#
# Recomputed from results/tables/06_modules_WGCNA/group_effects/*/
# module_group_effects.csv. Each dataset carries its own BH families.
W <- function(family, scope, n, n_sup, minfdr, note)
  data.frame(fdr_family = family, scope = scope, n_tests = n,
             n_fdr_supported = n_sup, smallest_fdr = minfdr, note = note,
             stringsAsFactors = FALSE)
wg <- rbind(
  W("FDR_primary_global", "neuron_neuropil; 15 modules x SUS-RES", 15, 0,
    0.245023, "the primary module-phenotype endpoint for the neuropil"),
  W("FDR_secondary_global", "neuron_neuropil; 15 modules x (RES-CON, SUS-CON)",
    30, 0, 0.611653, "the two secondary contrasts"),
  W("primary + secondary", "neuron_neuropil; 15 modules x 3 group contrasts",
    45, 0, 0.245023,
    "THIS IS THE '0 of 45'. It is NEUROPIL-ONLY and covers the three group contrasts. It excludes the interaction omnibus."),
  W("FDR_primary_global", "microglia; 13 modules x SUS-RES", 13, 0, 0.442632, ""),
  W("FDR_primary_global", "neuron_soma; 7 modules x SUS-RES", 7, 0, 0.155629,
    "smallest primary FDR anywhere in the module layer"),
  W("primary + secondary, all datasets",
    "35 modules (neuropil 15, microglia 13, soma 7) x 3 group contrasts", 105, 0,
    0.155629,
    "the all-compartment equivalent of the 45; smallest FDR is the soma primary value"),
  W("FDR_interaction_omnibus",
    "35 modules across three datasets; StressGroup x SpatialUnit", 35, 0,
    0.274144,
    "THIS IS THE '0 of 35'. One omnibus test per module, all three compartments; smallest value is the neuropil one."),
  W("FDR_conservative_all_tests",
    "140 module x contrast cells across three datasets", 140, 0, 0.521571,
    "the most conservative pooling; still zero supported"))
utils::write.csv(wg, file.path(MS, "wgcna_45_contract.csv"), row.names = FALSE)

# ============================================== methods statement provenance
M <- function(id, sec, fact, src, loc, pname, pval, verified, notes)
  data.frame(methods_id = id, section = sec, sentence_or_fact = fact,
             authoritative_source = src, source_location = loc,
             parameter_name = pname, parameter_value = pval,
             verified = verified, notes = notes, stringsAsFactors = FALSE)
mp <- rbind(
  M("M-01", "Animal-level bilateral aggregation",
    "Hemispheres are aggregated to one animal-level value before any group comparison",
    "manuscript_methods_contract.csv", "row preprocessing_bilateral",
    "analysis unit", "protein group x animal x spatial unit", "VERIFIED",
    "the biological replicate is the animal; hemispheres are repeated tissue"),
  M("M-02", "Differential protein abundance",
    "Moderated linear model on animal-level values within each spatial unit",
    "manuscript_methods_contract.csv", "row differential_abundance",
    "contract id", "animal_level_protigy_da_v1", "VERIFIED",
    "BH within each comparison; SUS-RES is the primary contrast"),
  M("M-03", "Ranked GO enrichment",
    "Genes ranked by the median moderated t per official gene symbol",
    "manuscript_methods_contract.csv", "row gsea", "ranking statistic",
    "moderated t, median per SYMBOL", "VERIFIED",
    "54 comparisons; no rank-statistic fallback"),
  M("M-04", "Ranked GO enrichment",
    "clusterProfiler::gseGO over GO biological process",
    "04_differential_expression_enrichment/01_clusterProfiler.r", "gseGO call",
    "minGSSize / maxGSSize / pAdjustMethod", "10 / 800 / BH", "VERIFIED", ""),
  M("M-05", "Ranked GO enrichment", "Enrichment p-values are floored at eps",
    "clusterProfiler default", "gseGO", "eps", "1e-10", "VERIFIED",
    "never overridden; a term at the floor has an unresolved true p"),
  M("M-06", "Curated GO-program atlas",
    "Seven ontology-defined program families, anchor-based",
    "config/manuscript_go_theme_registry.tsv", "theme_role == primary",
    "registry version", "manuscript_go_themes_v3", "VERIFIED",
    "27 rows; match_scope in {anchor_and_descendants, exact_go_id, exclude_anchor_and_descendants}"),
  M("M-07", "Curated GO-program atlas",
    "The mitochondrial family excludes the glycolysis sub-DAG",
    "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
    "theme_id == mitochondrial_respiration_oxphos",
    "constituent GO terms / glycolytic members", "16 / 0", "VERIFIED",
    "PH-012: two reader-facing artefacts still say 20 terms including glycolysis. 324 glycolytic rows carry theme_id empty and assignment_status unclassified. PDH GO:0006086 and TCA GO:0006099 retained."),
  M("M-08", "CAMERA sensitivity analysis",
    "limma::cameraPR preranked on the exact canonical ranked statistic",
    "manuscript_methods_contract.csv", "row camera_sensitivity",
    "inter.gene.cor", "0.01 (prespecified; band 0, 0.01, 0.05)", "VERIFIED",
    "sensitivity only; never validation"),
  M("M-09", "WGCNA",
    "Module eigengenes related to outcome with a linear mixed model",
    "manuscript_methods_contract.csv", "row wgcna", "model",
    "lmerTest eigengene ~ StressGroup + SpatialUnit + (1|AnimalID)", "VERIFIED",
    "random intercept for animal"),
  M("M-10", "WGCNA", "The displayed module-phenotype family is neuropil-only",
    "results/tables/06_modules_WGCNA/group_effects/*/module_group_effects.csv",
    "FDR_primary_global + FDR_secondary_global", "n tests / supported",
    "45 / 0 (neuropil); 105 / 0 across all three compartments", "VERIFIED",
    "resolves the previously ambiguous 0 of 45; see wgcna_45_contract.csv"),
  M("M-11", "External-reference validation",
    "Each internal contrast is tested against one external signature at a time",
    "04_differential_expression_enrichment/09_control_spatial_identity_validation.r",
    "line 591-601, TERM2GENE = data.frame(term = job$external_signature, gene = job$mapped)",
    "gene-set collection size per call", "1", "VERIFIED",
    "MT-04: because the collection has one term, clusterProfiler's BH adjustment is a no-op, so the field named p_adjust is numerically identical to the raw p. This is a scope artefact of the per-pairing design, NOT a coding error."),
  M("M-12", "External-reference validation",
    "The genuine multi-test correction is applied within signature families",
    "04_differential_expression_enrichment/09_control_spatial_identity_validation.r",
    "control_spatial_signature_family() at line 503", "family sizes",
    "soma_tissue 12, neuropil_subregion 6, ca1_laminar 12", "VERIFIED",
    "stored as signature_FDR; this is the value Methods must cite, not p_adjust"),
  M("M-13", "External-reference validation",
    "Nine control-only internal contrasts against seven published signatures",
    "results/source_data/.../v9_ed_external_full_source_data.csv", "all 30 rows",
    "contrasts / signatures / pairings", "9 / 7 / 30", "VERIFIED",
    "5 neuropil + 4 neuronal soma; zero microglia contrasts"),
  M("M-14", "Spatial molecular-identity analyses",
    "Bilateral concordance is computed on prespecified contrasts in control animals",
    "results/tables/11_spatial_systems/bilateral/bilateral_spatial_identity_summary.csv",
    "manuscript_locked column", "prespecified contrasts",
    "11 locked (7 neuropil + 4 soma); 4 microglia regional-context, not locked",
    "VERIFIED",
    "MT-13: DG_MO and DG_PO are algebraic mirror images (DG neuropil has two layers) and are identical to 7 dp"),
  M("M-15", "Statistics", "The biological replicate is the animal",
    "docs/MANUSCRIPT_STATISTICAL_CONTRACT.md", "section 5", "n per group", "3",
    "VERIFIED",
    "9 animals; acquisitions, hemispheres and animal x dataset network instances are never biological replicates"),
  M("M-16", "Behavioural outcome definition", "CombZ composite outcome score",
    "data/external/behavior/E9_Behavior_Data.xlsx", "sheet zScore",
    "construction formula", "[METHOD DETAIL UNRESOLVED]", "NOT_VERIFIABLE",
    "the column is precomputed in the workbook; no tracked code constructs it. Components visible: NOR, sucrose_pref, weight_dev, delta_cort, adrenal_weight, spleen_weight"),
  M("M-17", "Phenotype classification", "RES / SUS assignment",
    "data/external/behavior/auc_individual_animals_*.csv", "Group column",
    "threshold / cut-point", "[METHOD DETAIL UNRESOLVED]", "NOT_VERIFIABLE",
    "the label arrives already assigned; no classifier exists in this repository"),
  M("M-18", "Home-cage monitoring", "Per-animal movement AUC",
    "data/external/behavior/auc_individual_animals_firstChangeActive.csv",
    "prediction_type column", "model", "subject_specific_gamm_observed_grid",
    "PARTIAL",
    "[METHOD DETAIL UNRESOLVED] the GAMM formula, smoothers and random-effect structure are not in this repository; only the fitted AUC is"))
utils::write.csv(mp, file.path(MS, "methods_statement_provenance.csv"),
                 row.names = FALSE)

cat("figure1_authoritative_source_inventory.csv rows:", nrow(f1), "\n")
print(table(f1$status))
cat("\nwgcna_45_contract.csv rows:", nrow(wg), "\n")
cat("  '0 of 45' =", wg$n_tests[wg$scope == "neuron_neuropil; 15 modules x 3 group contrasts"],
    "neuropil tests, smallest FDR",
    wg$smallest_fdr[wg$scope == "neuron_neuropil; 15 modules x 3 group contrasts"], "\n")
cat("\nmethods_statement_provenance.csv rows:", nrow(mp), "\n")
print(table(mp$verified))
