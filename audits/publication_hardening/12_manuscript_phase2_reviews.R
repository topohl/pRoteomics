#!/usr/bin/env Rscript

# Manuscript phase 2: the review and contract tables that record what could NOT
# be established, alongside the Methods red-team.
#
# WRITES ONLY manuscript/*.csv. Touches no scientific output.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
MS <- "manuscript"
dir.create(MS, recursive = TRUE, showWarnings = FALSE)

# ================================================== Methods red-team review
RT <- function(id, defect, check, verdict, evidence)
  data.frame(review_id = id, reproducibility_defect = defect,
             what_was_checked = check, verdict = verdict, evidence = evidence,
             stringsAsFactors = FALSE)
mrt <- rbind(
  RT("MRT-01", "wrong biological n", "every n stated in Methods", "CLEAN",
     "n = 3 per group, 9 animals, stated once at the head and never contradicted; acquisitions and hemispheres explicitly excluded as replicates"),
  RT("MRT-02", "wrong analysis unit", "DA, GSEA, WGCNA, network units", "CLEAN",
     "DA and GSEA on animal-level values within spatial unit; WGCNA on module eigengenes with a per-animal random intercept; network on animal x dataset instances, named as such"),
  RT("MRT-03", "omitted repeated-measures structure",
     "WGCNA model and the bilateral design", "CLEAN",
     "(1|AnimalID) stated; bilateral aggregation stated before any group comparison; the two DG contrasts flagged as algebraic mirror images"),
  RT("MRT-04", "wrong contrast orientation",
     "primary contrast and the algebraic relation", "CLEAN",
     "SUS-RES named primary; SUS-RES = SUS-CON minus RES-CON stated in both Methods and Results"),
  RT("MRT-05", "wrong BH family", "every declared family",
     "DEFECT FOUND AND FIXED",
     "The external-validation section would have cited p_adjust as the correction. Methods now states that p_adjust is a single-set no-op and that signature_FDR, applied within families of 12, 6 and 12, is the operative correction. MT-04."),
  RT("MRT-06", "software function named incorrectly",
     "gseGO, cameraPR, lmerTest, GSEA", "CLEAN",
     "each checked against the calling script; cameraPR is the preranked form, not CAMERA on an expression matrix"),
  RT("MRT-07", "parameter copied from superseded code",
     "atlas membership", "DEFECT AVOIDED",
     "PH-012: two reader-facing artefacts still describe a 20-term mitochondrial theme containing glycolysis. Methods states the verified frozen membership: 16 terms, glycolysis excluded, PDH and TCA retained."),
  RT("MRT-08", "method described more strongly than implemented",
     "CAMERA, external validation, atlas", "CLEAN",
     "CAMERA is sensitivity only; the atlas carries no theme-level FDR; the off-target pairings are explicitly not a specificity test"),
  RT("MRT-09", "figure post-processing presented as inference",
     "atlas median NES; CA2-SLM qualification", "CLEAN",
     "median NES named descriptive; the robustness qualification explicitly not a second FDR family and not an additional test"),
  RT("MRT-10", "ambiguous denominator", "the WGCNA 0 of 45",
     "DEFECT FOUND AND FIXED",
     "Recomputed: 45 = 15 neuropil modules x 3 group contrasts, NEUROPIL ONLY, smallest FDR 0.245. The all-compartment equivalent is 105 tests, smallest FDR 0.156. Both now stated in Methods and in Results 3."),
  RT("MRT-11", "unrecoverable parameter filled from convention",
     "CombZ, RES/SUS cut-point, GAMM specification",
     "CLEAN - MARKED UNRESOLVED",
     "three [METHOD DETAIL UNRESOLVED] markers rather than invented values; rows M-16, M-17, M-18"),
  RT("MRT-12", "version mismatch", "package versions, seeds, ontology release",
     "OUTSTANDING",
     "METHODS TODO MT-01: still to be inserted from the run manifests"))
utils::write.csv(mrt, file.path(MS, "methods_red_team_review.csv"),
                 row.names = FALSE)

# ================================================= Figure 1 red-team review
#
# The checklist was applied in full. Every question is unanswerable from this
# repository, which is itself the finding.
F <- function(id, q, verdict, evidence)
  data.frame(review_id = id, question = q, verdict = verdict,
             evidence = evidence, stringsAsFactors = FALSE)
f1rt <- rbind(
  F("F1RT-01", "Is this expected because of how RES and SUS were defined?",
    "UNANSWERABLE",
    "The classification rule is not in this repository, so circularity cannot be assessed."),
  F("F1RT-02", "Is this independent evidence?", "UNANSWERABLE",
    "Requires the outcome-score construction, which is precomputed upstream."),
  F("F1RT-03", "Does 'predict' mean held-out prediction?",
    "NO - AND NO PREDICTION EXISTS",
    "No cross-validation or prediction machinery anywhere in the repository. The verb is not usable."),
  F("F1RT-04", "Was the predictor collected before the outcome?",
    "UNVERIFIABLE",
    "No ages, dates or cage-change timestamps appear in any tracked file."),
  F("F1RT-05", "Is the cross-validation unit the animal?", "NOT APPLICABLE",
    "No cross-validation exists."),
  F("F1RT-06", "Is sex-specific wording licensed?", "NO",
    "Sex is carried as a column in both external inputs, but no stratified or interaction model is fitted here."),
  F("F1RT-07", "Could cage or batch explain the result?", "UNANSWERABLE",
    "Batch is present in both external inputs; no model in this repository adjusts for it."),
  F("F1RT-08", "Is the stated n biological animals?", "NOT APPLICABLE",
    "No Figure 1 claim is made, so no n is asserted."))
utils::write.csv(f1rt, file.path(MS, "figure1_red_team_review.csv"),
                 row.names = FALSE)

# ============================== timeline / outcome score / prediction contracts
tl <- data.frame(
  stage = c("early home-cage monitoring", "later outcome assessment",
            "tissue collection"),
  age_or_window = "[UNRESOLVED - not recorded in this repository]",
  measure = c("RFID movement, GAMM-derived AUC per animal",
              "NOR, sucrose preference, weight deviation, delta corticosterone, adrenal and spleen weight",
              "laser-capture microdissection across 18 spatial units"),
  analysis_role = c("predictor in the upstream analysis", "outcome components",
                    "proteomic measurement"),
  used_in_outcome_score = c("no", "yes", "no"),
  used_as_predictor = c("yes, upstream only", "no", "no"),
  notes = c("prediction_type = subject_specific_gamm_observed_grid; the model is not in this repository",
            "CombZ is precomputed in the source workbook",
            "9 animals, 3 per group"),
  stringsAsFactors = FALSE)
utils::write.csv(tl, file.path(MS, "figure1_timeline_contract.csv"),
                 row.names = FALSE)

os <- data.frame(
  element = c("canonical manuscript term", "component measures",
              "sign orientation", "z-scoring population", "aggregation formula",
              "treatment of sex", "treatment of batch", "RES/SUS boundary",
              "do controls enter score construction",
              "leave-one-out score construction"),
  value = c("CombZ (composite outcome score)",
            "NOR, sucrose_pref, weight_dev, delta_cort, adrenal_weight, spleen_weight",
            rep("[UNRESOLVED]", 8)),
  source = c("data/external/behavior/E9_Behavior_Data.xlsx sheet zScore",
             "same workbook, column names",
             rep("NOT IN REPOSITORY - precomputed upstream", 8)),
  stringsAsFactors = FALSE)
utils::write.csv(os, file.path(MS, "outcome_score_contract.csv"),
                 row.names = FALSE)

bp <- data.frame(
  element = c("CV type", "grouping unit", "held-out animals",
              "sex stratification", "batch or cage grouping",
              "feature-selection location", "scaling location",
              "hyperparameter tuning", "outcome thresholding",
              "performance metric", "uncertainty interval", "permutation null"),
  value = "NO PREDICTION ANALYSIS EXISTS IN THIS REPOSITORY",
  evidence = "Repository-wide search found no glmnet, caret, pROC, randomForest, cv.glmnet, trainControl or createFolds. Every leave-one-animal-out match is a proteomic stability analysis.",
  consequence = "The verb 'predicts' is not licensed. 'Was associated with' becomes available only if an association analysis is recovered from the upstream repository.",
  stringsAsFactors = FALSE)
utils::write.csv(bp, file.path(MS, "behavior_prediction_contract.csv"),
                 row.names = FALSE)

cat("methods_red_team_review.csv rows:", nrow(mrt), "\n")
print(table(mrt$verdict))
cat("\nfigure1_red_team_review.csv rows:", nrow(f1rt), "\n")
print(table(f1rt$verdict))
cat("\nfigure1_timeline_contract.csv rows:", nrow(tl),
    "| outcome_score_contract.csv rows:", nrow(os),
    "| behavior_prediction_contract.csv rows:", nrow(bp), "\n")
