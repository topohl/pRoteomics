#!/usr/bin/env Rscript
# Review workbooks for the normalisation decision (2026-10-21). Read-only for the pipeline; reads the tables written by
# tools/normalisation_review/*.R and writes three house-style workbooks:
#
#   Proteomics_normalisation_comparison.xlsx    the two centrings and independent estimators, side by side
#   Proteomics_WGCNA_annotation_decisions.xlsx  editable decision table for every module and supermodule
#   Proteomics_normalisation_provenance.xlsx    every analysis against the two versions; checks; decisions
#
# Numbers quoted in the Read me sheets are computed from the result tables below. Typed exceptions: the original archive's
# file count and size (its ARCHIVE_SUMMARY.txt) and the verification bounds 6e-15 / 2e-14 (validation/verification_additive_offsets.csv).
# USAGE  Rscript tools/normalisation_review/build_normalisation_review_workbooks.R      (working directory: repo root)

source(if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R"))
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(readr); library(tidyr) })
source(repo_path("R", "review_workbook_style.R"))
REV <- path_results("reviewer_audit", "normalisation_choice_2026-10-21")
VAL <- file.path(REV, "validation"); ANN <- file.path(REV, "annotation_review"); INV <- file.path(REV, "inventory")
OUT <- file.path(REV, "workbooks"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(dir, f) if (file.exists(file.path(dir, f))) read_csv(file.path(dir, f), show_col_types = FALSE, guess_max = 1e5) else NULL
CMP <- c(neuron_neuropil = "Neuropil", neuron_soma = "Soma", microglia = "Microglia")
EST_LABEL <- c(original_all_observed = "Original: median of all observed values",
               new_complete_case = "New: median of proteins quantified in every sample",
               median_ratio_complete_case = "Median ratio to the protein mean, same anchor proteins",
               median_ratio_shared50 = "Median ratio, proteins observed in at least half of the samples",
               joint_shared_core = "Joint shared-core median across compartments (QC stage; fixed, so in-sample)",
               none_diann_only = "No further centring (DIA-NN output as delivered)")
lab_est <- function(x) unname(ifelse(x %in% names(EST_LABEL), EST_LABEL[x], x))
lab_cmp <- function(x) unname(ifelse(x %in% names(CMP), CMP[x], x))
num_cols <- function(df) which(vapply(df, is.numeric, logical(1)))
sig3 <- function(df) {   # 3 significant figures for measured values; counts (whole numbers) are left as they are
  nc <- num_cols(df); nc <- nc[!vapply(df[nc], function(v) all(is.na(v) | v == round(v)), logical(1))]
  df[nc] <- lapply(df[nc], function(v) signif(v, 3)); df }
rng <- function(v, d = 2) { v <- round(v[is.finite(v)], d) + 0; if (!length(v)) return("NA"); f <- function(x) sub("^-(0[.]0+)$", "\\1", formatC(x, format = "f", digits = d)); a <- f(min(v)); b <- f(max(v)); if (a == b) a else paste0(a, " to ", b) }
fmt <- function(x) format(x, big.mark = ",")
put <- function(wb, sheet, df, ...) { if (is.null(df) || !nrow(df)) return(invisible()); put_table(wb, sheet, df, ...) }

# ----------------------------------------------------------------------------------- numbers for the text
hf <- rd(VAL, "addendum_heldout_by_functional.csv"); wu <- rd(VAL, "addendum_offset_vs_missing_within_unit.csv")
iu <- rd(VAL, "addendum_intensity_dependence_within_unit.csv"); gp <- rd(VAL, "addendum_group_offsets_plate_adjusted.csv")
mc <- rd(VAL, "addendum_microglia_anchor_carryover.csv"); rv <- rd(VAL, "residual_variance_by_estimator.csv")
cv <- rd(VAL, "ca2slm_variants_and_unit_only_model.csv"); cs <- rd(VAL, "ca2slm_sample_detail.csv")
hits <- rd(VAL, "ca2slm_hits_sample_dependence.csv"); conc <- rd(VAL, "offset_concordance_between_estimators.csv")
ws <- rd(VAL, "wgcna_stability_by_estimator.csv"); wc <- rd(VAL, "wgcna_stability_common_power.csv")
wp <- rd(VAL, "wgcna_leave_one_animal_out_common_power_paired.csv"); ew <- rd(VAL, "ewce_normalisation_vs_background.csv")
foc <- rd(VAL, "heldout_focus_samples.csv"); dz <- rd(VAL, "diann_level_summary.csv")
hemi <- rd(VAL, "hemisphere_concordance.csv"); rc <- rd(VAL, "ca2slm_removal_control.csv"); lp <- rd(VAL, "limma_vs_protigy.csv"); pqc <- rd(VAL, "precentring_qc_variance_explained.csv")
anc <- rd(VAL, "anchor_set_characteristics.csv"); design <- rd(VAL, "design_plate_by_animal.csv"); comp <- rd(VAL, "anchor_functional_composition.csv")
dia_all <- readxl::read_excel(file.path("data", "raw", "pg_matrix", "quicksearch.stats.annotated.xlsx"))
dia_used <- rd(VAL, "diann_level_samples.csv"); dia_excl <- dia_all[!dia_all$sample_id %in% dia_used$sample_id, ]
ANCHOR_EST <- c("new_complete_case", "median_ratio_complete_case")
held <- hf |> filter(set == "anchor", startsWith(mode, "re"))                     # re-estimated without the held-out proteins
o_hold <- held |> filter(estimator == "original_all_observed"); n_hold <- held |> filter(estimator == "new_complete_case")
m_hold <- held |> filter(estimator == "median_ratio_complete_case")
rowc <- held |> filter(functional != "raw_median")
a755 <- cs[cs$label == "A755_L_CA2_slm", ]
f_part <- foc |> filter(set == "held_out_partially_observed", label == "A755_L_CA2_slm", estimator == "new_complete_case")
cvn <- function(est, var) cv$n_fdr05_pipeline_model[cv$estimator == est & cv$variant == var]
pct <- function(a, b) rng(100 * (b - a) / b, 0)                                    # relative reduction of a against b, in %
rv_w <- rv |> select(dataset, estimator, median_residual_sd) |> pivot_wider(names_from = estimator, values_from = median_residual_sd)
gp_w <- function(col) gp |> select(dataset, estimator, v = all_of(col)) |> pivot_wider(names_from = estimator, values_from = v)
plate_tab <- design |> mutate(plate = ifelse(B > 0, "B", "C")) |> count(plate, group) |> group_by(plate) |>
  summarise(txt = paste(sprintf("%s %d", group, n), collapse = "/"), .groups = "drop")
sft_r2 <- function(ds, est, pw) { f <- file.path(VAL, "wgcna_builds", sprintf("%s_%s_softpower.rds", ds, est)); if (!file.exists(f)) return(NA_real_)
  x <- readRDS(f); x <- if (is.data.frame(x)) x else x$fitIndices; v <- x$SFT.R.sq * -sign(x$slope); v[x$Power == pw][1] }
new_hits <- hits |> filter(version == "new"); new_up <- new_hits |> filter(log2fc > 0)
rc_rem <- rc |> filter(removed_group != "none")
rc_with <- rc_rem |> filter(grepl("A755_L", removed)); rc_without <- rc_rem |> filter(!grepl("A755_L", removed))
rc_n <- function(x, est) rng(x$SUS_vs_RES[x$estimator == est], 0)
fam_win <- function(fun, md) {                     # which anchor-family estimator has the lower residual SD, per compartment
  x <- hf |> filter(set == "anchor", functional == fun, mode == md, estimator %in% ANCHOR_EST) |> select(dataset, estimator, residual_sample_sd) |>
    pivot_wider(names_from = estimator, values_from = residual_sample_sd)
  paste(sprintf("%s %s", tolower(lab_cmp(x$dataset)), ifelse(x$median_ratio_complete_case < x$new_complete_case, "median ratio", "anchor median")), collapse = ", ") }
RE <- "re-estimated without held-out proteins"; DEP <- "as deployed (in-sample)"
fold_red <- held |> filter(estimator %in% c("original_all_observed", "new_complete_case")) |> select(dataset, functional, estimator, residual_sample_sd) |>
  pivot_wider(names_from = estimator, values_from = residual_sample_sd) |> mutate(f = original_all_observed / new_complete_case)
mg_comp <- comp |> filter(dataset == "microglia", set == "anchor", grepl("myelin|synap|spine", Description))
# held-out residual of the original as deployed (its own all-observed median, not re-estimated on the depleted set) against the new re-estimated
fold_dep <- hf |> filter(set == "anchor", mode == DEP, estimator == "original_all_observed") |> select(dataset, functional, o = residual_sample_sd) |>
  inner_join(n_hold |> select(dataset, functional, n = residual_sample_sd), by = c("dataset", "functional")) |> mutate(f = o / n)
rawc <- held |> filter(functional == "raw_median")
# SUS minus RES offset per region-layer unit (hemispheres averaged per animal, then animals per group): what each unit's logFC is shifted by
offs <- rd(VAL, "sample_offsets_all_estimators.csv")
pu <- offs |> mutate(unit = sub("^[^_]+_[LR]_", "", label)) |> select(dataset, animal, unit, group, starts_with("rel_")) |>
  pivot_longer(starts_with("rel_"), names_to = "estimator", values_to = "o") |> mutate(estimator = sub("^rel_", "", estimator)) |>
  group_by(dataset, unit, estimator, group, animal) |> summarise(o = mean(o), .groups = "drop") |>
  group_by(dataset, unit, estimator, group) |> summarise(o = mean(o), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = o) |> mutate(sus_res = SUS - RES)
pu_w <- pu |> select(dataset, unit, estimator, sus_res) |> pivot_wider(names_from = estimator, values_from = sus_res)
pu_neutral <- pu_w |> filter(abs(original_all_observed) < abs(new_complete_case))
pu_neutral_txt <- paste(sprintf("%s %s %.3f vs %.3f", tolower(lab_cmp(pu_neutral$dataset)), sub("_microglia$", "", pu_neutral$unit),
                                pu_neutral$original_all_observed, pu_neutral$new_complete_case), collapse = "; ")
dshift <- rd(VAL, "da_shift_vs_new_by_estimator.csv") |> filter(estimator == "original_all_observed", comparison == "SUS vs RES")
dshift_max <- dshift[which.max(abs(dshift$logFC_shift_vs_new)), ]
# partially observed proteins: is the original lower than both anchor estimators, per statistic (re-estimated)?
po <- hf |> filter(set == "partially_observed", mode == RE, estimator %in% c("original_all_observed", ANCHOR_EST)) |>
  select(dataset, functional, estimator, residual_sample_sd) |> pivot_wider(names_from = estimator, values_from = residual_sample_sd) |>
  mutate(orig_lower = original_all_observed < pmin(new_complete_case, median_ratio_complete_case))
po_not <- po |> filter(!orig_lower)
ew_c <- if (!is.null(ew)) ew |> filter(grepl("^centring", effect), analysis_type == "Differential") |>
  mutate(jac = significant_in_both / (significant_1 + significant_2 - significant_in_both)) else NULL
new_up <- new_up |> mutate(any_imputed = imputed_A755L | imputed_A764R | n_imputed_SUS > 0 | n_imputed_RES > 0)
rc_grp <- function(g, est, a755 = FALSE) { x <- rc_rem |> filter(removed_group == g, estimator == est, grepl("A755_L", removed) == a755); rng(x$SUS_vs_RES, 0) }
ag <- rd(VAL, "anchor_definition_by_group.csv")
pqc_raw <- pqc |> filter(response == "raw_original_all_observed", predictor == "missing_fraction")
FUN_LAB <- c(raw_median = "raw median", row_centred_median = "row-centred median", row_centred_mean = "row-centred mean")
oc <- bind_rows(lapply(c("SUS_minus_CON_plate_adjusted", "RES_minus_CON_plate_adjusted"), function(col)
  gp_w(col) |> mutate(contrast = sub("_plate_adjusted$", "", sub("_minus_", " minus ", col))))) |> filter(abs(original_all_observed) < abs(new_complete_case))
oc_txt <- if (nrow(oc)) paste(sprintf("%s %s (%.3f vs %.3f)", tolower(lab_cmp(oc$dataset)), oc$contrast, oc$original_all_observed, oc$new_complete_case), collapse = "; ") else "none"
imp_parts <- function(r) paste(c(if (r$imputed_A755L) "A755 L itself", if (r$imputed_A764R) "A764 R",
                                 if (r$n_imputed_SUS - r$imputed_A755L - r$imputed_A764R > 0) sprintf("%d other SUS", r$n_imputed_SUS - r$imputed_A755L - r$imputed_A764R),
                                 if (r$n_imputed_RES > 0) sprintf("%d RES", r$n_imputed_RES)), collapse = " and ")
imp_txt <- paste(vapply(which(new_up$any_imputed), function(i) sprintf("%s: %s", new_up$gene_symbol[i], imp_parts(new_up[i, ])), ""), collapse = "; ")
ari_family <- c(wc$ari_vs_current_new_network[wc$estimator == "median_ratio_complete_case"], ws$ari_vs_current_new_network[ws$estimator == "median_ratio_complete_case"])
txt_wgcna <- if (!is.null(wc)) {
  pr <- wp |> group_by(dataset) |> summarise(better = sum(new_minus_original > 0), n = n(), .groups = "drop")
  s <- wc |> select(dataset, estimator, common_power, loao_ari_mean) |> pivot_wider(names_from = estimator, values_from = loao_ari_mean)
  s$r2_orig <- mapply(sft_r2, s$dataset, "original_all_observed", s$common_power)
  me <- ws |> select(dataset, estimator, median_abs_cor_ME_vs_missing) |> pivot_wider(names_from = estimator, values_from = median_abs_cor_ME_vs_missing)
  mg_own <- ws |> filter(dataset == "microglia")
  paste0("WGCNA: at a common soft power, the new network's own power (", paste(sprintf("%s %d", tolower(lab_cmp(s$dataset)), s$common_power), collapse = ", "),
         "; there the original's scale-free fit R2 is ", paste(sprintf("%.2f", s$r2_orig), collapse = ", "), ", against the 0.80 rule), leave-one-animal-out stability (mean ARI) original / new / median ratio is ",
         paste(sprintf("%s %.2f / %.2f / %.2f", lab_cmp(s$dataset), s$original_all_observed, s$new_complete_case, s$median_ratio_complete_case), collapse = "; "),
         "; the new centring is more stable than the original for ", paste(sprintf("%d of %d animals in %s", pr$better, pr$n, tolower(lab_cmp(pr$dataset))), collapse = ", "),
         ". Microglia stability favours the original at both the common and the own powers (",
         sprintf("%.2f vs %.2f", mg_own$loao_ari_mean[mg_own$estimator == "original_all_observed"], mg_own$loao_ari_mean[mg_own$estimator == "new_complete_case"]),
         " at own powers). The median ratio and the anchor median, whose offsets differ by at most ", rng(conc$max_abs_difference[conc$estimator_a == "median_ratio_complete_case" & conc$estimator_b == "new_complete_case"]),
         " log2, give partitions that agree by ARI ", rng(ari_family), " (common and own powers); with one animal left out the mean ARI is ", rng(wc$loao_ari_mean),
         " at the common power and ", rng(ws$loao_ari_mean), " at the own powers (lowest single animal ", rng(min(c(wc$loao_ari_min, ws$loao_ari_min))), "). Module identities are therefore fragile under every centring. At the selected powers module eigengenes track depth less under the new centring in microglia and neuropil (median |r| ",
         sprintf("%.2f vs %.2f and %.2f vs %.2f", me$new_complete_case[me$dataset == "microglia"], me$original_all_observed[me$dataset == "microglia"],
                 me$new_complete_case[me$dataset == "neuron_neuropil"], me$original_all_observed[me$dataset == "neuron_neuropil"]),
         "); in soma they are equal (", sprintf("%.2f vs %.2f", me$new_complete_case[me$dataset == "neuron_soma"], me$original_all_observed[me$dataset == "neuron_soma"]), ").")
} else "Common-power WGCNA stability is pending."
txt_ewce <- if (!is.null(ew)) {
  b <- ew |> filter(grepl("^background", effect)); c <- ew |> filter(grepl("^centring", effect)); cb <- c |> filter(analysis_type == "Baseline")
  cd <- c |> filter(analysis_type == "Differential") |> mutate(jac = significant_in_both / (significant_1 + significant_2 - significant_in_both))
  paste0("With the same measured-proteome background, the centring changes differential-target scores little (score r ",
         rng(cd$cor_sd_from_mean), ", call agreement ", rng(100 * cd$call_agreement, 0),
         "% over all tests), but only ", rng(100 * cd$jac, 0), "% of the significant differential calls are shared (",
         paste(sprintf("%s %d of %d/%d", tolower(lab_cmp(cd$dataset)), cd$significant_in_both, cd$significant_1, cd$significant_2), collapse = "; "),
         ", original/new); baseline scores are almost unchanged (r ", rng(cb$cor_sd_from_mean), ", call agreement ", rng(100 * cb$call_agreement, 1), "%; ",
         rng(cb$significant_1 + cb$significant_2 - 2 * cb$significant_in_both, 0),
         " calls differ through the global FDR). The change from the archived run comes from the step from run A to run C (differential r ",
         rng(b$cor_sd_from_mean[b$analysis_type == "Differential"]), "), which changes the background, the species space and the code since b478a13 together.")
} else "EWCE comparison pending."

# =================================================================================== 1. normalisation comparison
wb <- createWorkbook()
add_readme(wb, "Proteomics: choice of per-sample centring", list(
  list(heading = "What is compared", lines = c(
    "Every sample is median-centred on the log2 scale before the 70% missingness filter and the down-shift imputation. Two centrings were used:",
    "  Original: median of all observed values of the sample (analysis before 2026-10-19).",
    "  New: median of the proteins quantified in every sample of the compartment (the anchors; re-normalisation of 2026-10-19).",
    "They are compared with a median ratio to each protein's mean on the same anchors, a median ratio on proteins observed in at least half of the samples, the joint shared-core median the QC stage uses, and no further centring (DIA-NN output as delivered).",
    "The filter uses missingness only and the imputation is relative to each sample's own mean and SD, so every downstream matrix under any estimator equals the current matrix plus one constant per sample; this is verified (sheet Checks), so no estimator needed a pipeline rerun.",
    "None of the criteria uses the number of significant proteins; DA counts are shown as consequences only (sheet DA consequences).")),
  list(heading = "Upstream quantification", lines = c(
    "DIA-NN (protein-group matrix and its stats.tsv columns, e.g. Normalisation.Instability, Median.RT.Prediction.Acc). The sequence sheets name a timsTOF Ultra with a dia-PASEF low-sample-amount method and the LC method WhisperZoom_SPD40_32.5 min, while the run names read Evo2_80SPDzoom: the throughput is to be confirmed with the facility. Two acquisition sequences: plate B 2025-07-03 and plate C 2025-07-07.",
    sprintf("%d runs in one matrix named 'quicksearch'. DIA-NN reported a normalisation instability for every run (%d analysed: median %.3f, max %.2f; the %d excluded runs %s), which indicates that its cross-run normalisation was on (to confirm from report.log.txt).",
            nrow(dia_all), dz$n[dz$dataset == "all compartments"], dz$median_norm_instability[dz$dataset == "all compartments"], dz$max_norm_instability[dz$dataset == "all compartments"],
            nrow(dia_excl), paste(formatC(sort(dia_excl$Normalisation.Instability), format = "f", digits = 3), collapse = ", ")),
    sprintf("Not in the delivered files: the DIA-NN version, the normalisation mode (global or RT-dependent), whether it ran across all %d runs jointly, and whether 'quicksearch' is the final search. The search log (report.log.txt) from the facility would answer all four.", nrow(dia_all)),
    sprintf("Each animal was acquired on one plate (animal nested in plate): %s on plate B, %s on plate C. Groups are unbalanced across plates (B: %s; C: %s), which is why group offsets are also given adjusted for plate (design: sheet Design; adjusted offsets: sheet Offsets by group).",
            paste(design$animal[design$B > 0], collapse = ", "), paste(design$animal[design$C > 0], collapse = ", "),
            plate_tab$txt[plate_tab$plate == "B"], plate_tab$txt[plate_tab$plate == "C"]))),
  list(heading = "Findings", lines = c(
    sprintf("1. The original offset tracks detection depth: given the region-layer unit, the missing fraction explains R2 %s of the raw original offset (within-unit r %s; pooled r %s); the anchor-based offsets (new, median ratio) correlate %s within units (%s including the shared-50%% median ratio). The missing fraction explains R2 %s of the original-minus-new difference, by construction (dropping low-abundance values raises an all-observed median). In microglia, depth itself co-varies with neuronal/myelin carry-over (finding 8), so part of it is composition.",
            rng(pqc_raw$partial_r2_given_unit), rng(wu$r_within_unit[wu$estimator == "original_all_observed"]), rng(wu$r_pooled[wu$estimator == "original_all_observed"]),
            rng(wu$r_within_unit[wu$estimator %in% ANCHOR_EST]), rng(wu$r_within_unit[wu$estimator %in% c(ANCHOR_EST, "median_ratio_shared50")]),
            rng(pqc$r2_alone[pqc$response == "original_minus_new" & pqc$predictor == "missing_fraction"])),
    sprintf("2. On proteins quantified in every sample, the original leaves a depth-correlated sample effect under every evaluation statistic. Re-estimated without the held-out proteins, the residual sample SD is %s log2 for the original (r with missingness %s), %s for the anchor median and %s for the median ratio: %s-fold lower under the new centring. Against the original as deployed (in-sample for the original, which avoids the extra error of re-estimating an all-observed median on a depleted protein set) the ratio is %s-fold. Under row-centred statistics the new centring leaves a smaller depth-correlated residual of opposite sign (r %s) and the median ratio does not (r %s); under the raw median the reverse holds (median ratio r %s, new r %s). Caveat: held-out anchors are the protein class the anchor estimators are fitted on.",
            rng(o_hold$residual_sample_sd, 3), rng(o_hold$cor_residual_vs_missing), rng(n_hold$residual_sample_sd, 3), rng(m_hold$residual_sample_sd, 3),
            rng(fold_red$f, 1), rng(fold_dep$f, 1),
            rng(rowc$cor_residual_vs_missing[rowc$estimator == "new_complete_case"]), rng(rowc$cor_residual_vs_missing[rowc$estimator == "median_ratio_complete_case"]),
            rng(rawc$cor_residual_vs_missing[rawc$estimator == "median_ratio_complete_case"]), rng(rawc$cor_residual_vs_missing[rawc$estimator == "new_complete_case"])),
    sprintf("3. Within-group residual SD (median limma sigma over all filtered proteins) is lower under every anchor-based estimator than under the original: new vs original %s (%s%% lower). The criterion is partly in-sample for the anchor estimators, because the anchors are %s%% of the proteins scored. With no centring at all it is also lower than the original in microglia and neuropil (%s), not in soma (%s).",
            paste(sprintf("%s %.3f vs %.3f", lab_cmp(rv_w$dataset), rv_w$new_complete_case, rv_w$original_all_observed), collapse = "; "), pct(rv_w$new_complete_case, rv_w$original_all_observed),
            rng(100 * anc$anchor_fraction_of_filtered, 0),
            paste(sprintf("%.3f vs %.3f", rv_w$none_diann_only[rv_w$dataset != "neuron_soma"], rv_w$original_all_observed[rv_w$dataset != "neuron_soma"]), collapse = " and "),
            sprintf("%.3f vs %.3f", rv_w$none_diann_only[rv_w$dataset == "neuron_soma"], rv_w$original_all_observed[rv_w$dataset == "neuron_soma"])),
    sprintf("4. SUS minus RES offsets. Averaged over region-layer units at animal level, the SUS-minus-RES offset difference is %s log2 under the original and %s under the new centring (plate-adjusted, as a measure of group neutrality: %s and %s). The deployed model tests each unit separately, so each unit's SUS-vs-RES logFC is shifted by that unit's own difference: %s log2 under the original and %s under the new centring (largest original-minus-new logFC shift %+.2f, %s %s; sheet DA shift vs new). In %d of %d units the original is the more neutral (%s). Among the SUS-vs-CON and RES-vs-CON contrasts (plate-adjusted), the original is more neutral only for %s; soma RES minus CON is far from zero under both (%.3f original, %.3f new). No offset difference is significant (p %s; 1-2 animals per group and plate), which is weak evidence of absence.",
            rng(gp$SUS_minus_RES[gp$estimator == "original_all_observed"], 3), rng(gp$SUS_minus_RES[gp$estimator == "new_complete_case"], 3),
            rng(gp$SUS_minus_RES_plate_adjusted[gp$estimator == "original_all_observed"], 3), rng(gp$SUS_minus_RES_plate_adjusted[gp$estimator == "new_complete_case"], 3),
            rng(pu_w$original_all_observed, 3), rng(pu_w$new_complete_case, 3),
            dshift_max$logFC_shift_vs_new, tolower(lab_cmp(dshift_max$dataset)), dshift_max$unit,
            nrow(pu_neutral), nrow(pu_w), pu_neutral_txt, oc_txt,
            gp_w("RES_minus_CON_plate_adjusted")$original_all_observed[gp_w("RES_minus_CON_plate_adjusted")$dataset == "neuron_soma"],
            gp_w("RES_minus_CON_plate_adjusted")$new_complete_case[gp_w("RES_minus_CON_plate_adjusted")$dataset == "neuron_soma"],
            rng(gp$p_group_given_plate[gp$estimator %in% c("original_all_observed", ANCHOR_EST)])),
    sprintf("5. Partially observed proteins (detected in 50-95%% of samples) are left-censored in low-depth samples, so no residual criterion can show which centring is right for them. The original gives a lower residual there than both anchor estimators in %d of %d compartment-by-statistic combinations (re-estimated; not in %s), because its over-correction offsets the censoring; the anchor estimators leave the censoring inflation in place (A755 L %+.2f log2 on held-out partially observed proteins under the new centring).",
            sum(po$orig_lower), nrow(po), if (nrow(po_not)) paste(sprintf("%s under the %s", tolower(lab_cmp(po_not$dataset)), FUN_LAB[po_not$functional]), collapse = ", ") else "none",
            f_part$mean_residual),
    sprintf("6. Within the anchor family neither estimator is better. Lower held-out residual (re-estimated): raw median %s; row-centred median, the median ratio's own statistic: %s; row-centred mean %s. As deployed (in-sample), each is lower under its own statistic. Results on partially observed proteins are mixed and not a correctness test (finding 5). Their offsets correlate r %s (largest difference %s log2).",
            fam_win("raw_median", RE), fam_win("row_centred_median", RE), fam_win("row_centred_mean", RE),
            rng(conc$pearson_r[conc$estimator_a == "median_ratio_complete_case" & conc$estimator_b == "new_complete_case"]),
            rng(conc$max_abs_difference[conc$estimator_a == "median_ratio_complete_case" & conc$estimator_b == "new_complete_case"])),
    sprintf("7. Sample effects are partly intensity-dependent: within units the high-minus-low intensity-tertile residual has SD %s log2 (%s%% of the within-unit original-minus-new SD, %s), of which %s (%s%%) is sample-specific (from hemisphere pairs). A single constant per sample cannot remove that part under either centring.",
            rng(iu$sd_high_minus_low_within_unit, 3), rng(100 * iu$sd_high_minus_low_within_unit / iu$sd_original_minus_new_within_unit, 0),
            rng(iu$sd_original_minus_new_within_unit, 3), rng(iu$sd_from_hemisphere_pairs, 3), rng(100 * iu$sd_from_hemisphere_pairs / iu$sd_original_minus_new_within_unit, 0)),
    sprintf("8. Microglia: all %d neuronal/myelin marker proteins tested are anchors, and neuronal/myelin GO terms are over-represented among the anchors (%s-fold). The markers' level relative to the anchor median (a carry-over score) correlates r %.2f with the missing fraction and r %.2f with the original-minus-new difference. Anchoring on the %d microglia-enriched anchors instead gives different offsets (r %.2f with the new offsets, up to %.1f log2 apart). In microglia-enriched ROIs, microglial loading and neighbourhood carry-over cannot be separated by any of the estimators.",
            mc$n_markers_in_anchors, rng(mg_comp$fold, 1), mc$cor_marker_score_vs_missing, mc$cor_marker_score_vs_original_minus_new,
            mc$n_microglia_enriched_anchors, mc$cor_microglia_enriched_anchor_offset_vs_new, mc$max_abs_diff_microglia_enriched_vs_new),
    sprintf("9. CA2-SLM SUS vs RES (%d proteins original, %d new) depends on single SUS hemispheres under every centring (A755 L and A764 R, the two SUS CA2-SLM samples with the highest missing fractions): without A755 L %d (original) and %d (new) remain, without A764 R %d and %d, without both %d and %d, with A755 represented by A755 L alone %s and %d. Under the new centring the result therefore depends on A764 R as much as on A755 L. A CA2-SLM-only model (fitted for the original, the new and the median ratio) finds none. A755 L is the lowest-depth neuropil acquisition before any centring (%s proteins identified, robust z %.1f; missing fraction z %+.1f; normalisation instability z %+.1f). Under the original, %d of %d hits are decreases, consistent with the offset over-correction of a low-depth sample (finding 1). Under the new centring, %d of %d hits are increases; all of them are detected in at least %d%% of neuropil samples (%d anchors), yet %d of the %d rest partly on imputed CA2-SLM values (%s). On them A755 L lies %s log2 above the mean of the four SUS hemispheres other than A755 L and A764 R (current imputed matrix). The mechanism of these increases is not established. %d of the %d new hits have padj 0.035-0.05.",
            cvn("original_all_observed", "all samples"), cvn("new_complete_case", "all samples"),
            cvn("original_all_observed", "without A755 L"), cvn("new_complete_case", "without A755 L"),
            cvn("original_all_observed", "without A764 R"), cvn("new_complete_case", "without A764 R"),
            cvn("original_all_observed", "without A755 L and A764 R"), cvn("new_complete_case", "without A755 L and A764 R"),
            format(cvn("original_all_observed", "without A755 R (A755 L alone)"), big.mark = ","), cvn("new_complete_case", "without A755 R (A755 L alone)"),
            fmt(a755$Proteins.Identified), a755$z_proteins_identified_neuropil, a755$z_missing_neuropil, a755$z_norm_instability_neuropil,
            sum(hits$version == "original" & hits$log2fc < 0), sum(hits$version == "original"), nrow(new_up), nrow(new_hits),
            as.integer(floor(100 * min(new_up$detection_fraction_neuropil))), sum(grepl("anchor", new_up$detection_class)),
            sum(new_up$any_imputed), nrow(new_up), imp_txt, rng(new_up$A755L_minus_other_SUS),
            sum(new_hits$padj >= 0.035), nrow(new_hits)),
    sprintf("9a. The CA2-SLM count is unstable in general. Removing one hemisphere from each of two animals of one group (%d removals per centring, %d of them without A755 L) changes the SUS-vs-RES count under the new centring (14 with all samples) to %s when two CON hemispheres are removed (they are not in the SUS-vs-RES contrast but enter the variance estimate), %s for two RES, %s for two SUS without A755 L and %s with A755 L; under the original (28) to %s, %s, %s and %s.",
            sum(rc_rem$estimator == "new_complete_case"), sum(rc_without$estimator == "new_complete_case"),
            rc_grp("CON", "new_complete_case"), rc_grp("RES", "new_complete_case"), rc_grp("SUS", "new_complete_case"), rc_grp("SUS", "new_complete_case", TRUE),
            rc_grp("CON", "original_all_observed"), rc_grp("RES", "original_all_observed"), rc_grp("SUS", "original_all_observed"), rc_grp("SUS", "original_all_observed", TRUE)),
    paste("10.", txt_wgcna), paste("11. EWCE:", txt_ewce))),
  list(heading = "Recommendation", lines = c(
    "Recommendation only: the choice of centring is pending your decision (config/normalisation.yml, decision: pending).",
    sprintf("The new centring is recommended over the original for proteins quantified in every sample (%s%% of the filtered proteins) and for neutrality of the SUS-vs-RES offsets: it lowers the residual sample SD on held-out fully quantified proteins %s-fold (%s-fold against the original as deployed), and the unit-averaged SUS-vs-RES offset difference is %s log2 under the original and %s under the new centring. The evidence is strongest in neuropil and microglia on residual criteria and weakest in soma.",
            rng(100 * anc$anchor_fraction_of_filtered, 0), rng(fold_red$f, 1), rng(fold_dep$f, 1),
            rng(gp$SUS_minus_RES[gp$estimator == "original_all_observed"]), rng(abs(gp$SUS_minus_RES[gp$estimator == "new_complete_case"]), 3)),
    sprintf("Criteria where the original does as well or better: microglia WGCNA leave-one-animal-out stability (finding 10); microglia hemisphere agreement in the tail of non-anchor proteins (90th percentile gap %.3f vs %.3f); soma median hemisphere gap on non-anchor proteins (%.3f vs %.3f); SUS-vs-RES offset neutrality in %d of %d units and soma SUS-vs-CON offset neutrality (finding 4); module eigengene depth correlation in soma (finding 10); and the residual on partially observed proteins (finding 5, not a correctness test). Neither centring is shown correct for partially observed proteins.",
            hemi$p90_abs_LR_gap_nonanchor[hemi$dataset == "microglia" & hemi$estimator == "original_all_observed"], hemi$p90_abs_LR_gap_nonanchor[hemi$dataset == "microglia" & hemi$estimator == "new_complete_case"],
            hemi$median_abs_LR_gap_nonanchor_observed_in_both[hemi$dataset == "neuron_soma" & hemi$estimator == "original_all_observed"],
            hemi$median_abs_LR_gap_nonanchor_observed_in_both[hemi$dataset == "neuron_soma" & hemi$estimator == "new_complete_case"],
            nrow(pu_neutral), nrow(pu_w)),
    sprintf("Not settled by the data: the choice within the anchor family (anchor median vs median ratio; offsets differ by at most %s log2, yet their WGCNA partitions agree only by ARI %s), and the treatment of partially observed proteins in low-depth samples, for which no constant centring is correct.",
            rng(conc$max_abs_difference[conc$estimator_a == "median_ratio_complete_case" & conc$estimator_b == "new_complete_case"]), rng(ari_family)),
    sprintf("Assumptions of the anchor-based centring: (1) the anchors carry no net group effect. This is not testable with 3 animals per group; group-specific anchor sets only show that the offsets depend little on which group defines completeness (r %s with the new offsets; largest difference %s log2, neuropil largest). (2) Per-sample loading differences are global rather than intensity- or protein-class-specific; this is partly violated (finding 7). (3) In microglia ROIs the anchor median reflects microglial loading; because neuronal/myelin proteins are over-represented among the anchors, loading and neighbourhood carry-over are confounded there (finding 8).",
            rng(ag$cor_with_new_offset), rng(ag$max_abs_difference_relative_offset)),
    sprintf("Alternatives for a sensitivity analysis: a missingness-aware differential model (detection-probability weighting or a hurdle model) and down-weighting or excluding the lowest-depth acquisitions (A755 L, A764 R). Independent of the choice, CA2-SLM protein-level results need such a sensitivity analysis before any claim; WGCNA module identities should be described as one partition among several plausible ones (ARI %s under the perturbations of finding 10).",
            rng(c(ari_family, wc$loao_ari_mean, wc$loao_ari_min, ws$loao_ari_mean, ws$loao_ari_min)))))))
# Every criterion evaluated, including those where the original does as well or better, in the order listed here.
# In-sample rows are labelled: the anchor estimators are fitted on the anchors, so they favour those estimators.
hf_rows <- bind_rows(lapply(names(FUN_LAB), function(fu) bind_rows(
  hf |> filter(set == "anchor", functional == fu, mode == RE) |>
    transmute(dataset, estimator, criterion = sprintf("Fully quantified proteins: residual sample SD, %s, held out (log2)", FUN_LAB[fu]), value = residual_sample_sd),
  hf |> filter(set == "anchor", functional == fu, mode == DEP) |>
    transmute(dataset, estimator, criterion = sprintf("Fully quantified proteins: residual sample SD, %s, as deployed (in-sample; log2)", FUN_LAB[fu]), value = residual_sample_sd),
  hf |> filter(set == "anchor", functional == fu, mode == RE) |>
    transmute(dataset, estimator, criterion = sprintf("Fully quantified proteins: residual vs missing fraction (r), %s, held out (0 is best)", FUN_LAB[fu]), value = cor_residual_vs_missing))))
crit <- bind_rows(
  hf_rows,
  wu |> transmute(dataset, estimator, criterion = "Offset vs missing fraction within units (partial r; 0 is best)", value = r_within_unit),
  rv |> transmute(dataset, estimator, criterion = "Within-group residual SD, animal level (median over proteins; partly in-sample)", value = median_residual_sd),
  gp |> transmute(dataset, estimator, criterion = "Offset SUS minus RES, animal level, unadjusted, averaged over units (log2; 0 is best)", value = SUS_minus_RES),
  pu |> group_by(dataset, estimator) |> summarise(value = max(abs(sus_res)), .groups = "drop") |>
    transmute(dataset, estimator, criterion = "Offset SUS minus RES per unit: largest absolute value, the largest shift of a unit's logFC (log2; 0 is best)", value),
  gp |> transmute(dataset, estimator, criterion = "Offset SUS minus RES, animal level, plate-adjusted (log2; 0 is best)", value = SUS_minus_RES_plate_adjusted),
  gp |> transmute(dataset, estimator, criterion = "Offset SUS minus CON, animal level, plate-adjusted (log2; 0 is best)", value = SUS_minus_CON_plate_adjusted),
  gp |> transmute(dataset, estimator, criterion = "Offset RES minus CON, animal level, plate-adjusted (log2; 0 is best)", value = RES_minus_CON_plate_adjusted),
  hemi |> transmute(dataset, estimator, criterion = "Hemisphere agreement: median |left-right gap|, anchors (in-sample; log2)", value = median_abs_LR_gap_anchor),
  hemi |> transmute(dataset, estimator, criterion = "Hemisphere agreement: median |left-right gap|, non-anchor proteins observed in both (log2)", value = median_abs_LR_gap_nonanchor_observed_in_both),
  hemi |> transmute(dataset, estimator, criterion = "Hemisphere agreement: 90th percentile |left-right gap|, non-anchor proteins (log2)", value = p90_abs_LR_gap_nonanchor),
  wc |> transmute(dataset, estimator, criterion = "WGCNA: leave-one-animal-out stability, mean ARI, common power (higher is better)", value = loao_ari_mean),
  ws |> transmute(dataset, estimator, criterion = "WGCNA: leave-one-animal-out stability, mean ARI, own power (higher is better)", value = loao_ari_mean),
  ws |> transmute(dataset, estimator, criterion = "WGCNA: module eigengene vs missing fraction, median |r|, own power (0 is best)", value = median_abs_cor_ME_vs_missing),
  hf |> filter(set == "partially_observed", functional == "raw_median", mode == RE) |>
    transmute(dataset, estimator, criterion = "Partially observed proteins: residual sample SD (left-censored; not a correctness test)", value = residual_sample_sd))
crit_order <- unique(crit$criterion)
score <- crit |> mutate(estimator = factor(estimator, levels = names(EST_LABEL))) |> arrange(estimator) |>
  mutate(estimator = as.character(estimator), dataset = lab_cmp(dataset)) |>
  pivot_wider(names_from = estimator, values_from = value) |> arrange(match(criterion, crit_order), dataset) |> relocate(criterion, dataset) |> sig3()
names(score) <- c("Criterion", "Compartment", lab_est(names(score)[-(1:2)]))
put(wb, "Scorecard", score, n2_cols = 3:ncol(score), left_cols = 1, widths = c(62, 13, rep(22, ncol(score) - 2)))
put(wb, "Design", rd(VAL, "design_plate_by_animal.csv"), widths = c(16, 10, 10, 10, 10))
put(wb, "DIA-NN level", sig3(rd(VAL, "diann_level_summary.csv") |> mutate(dataset = lab_cmp(dataset))), widths = c(18, rep(15, 11)))
off <- rd(VAL, "sample_offsets_all_estimators.csv") |>
  transmute(Compartment = lab_cmp(dataset), Sample = label, Group = group, Plate = plate, `Run order` = run_order,
            `Proteins identified (DIA-NN)` = Proteins.Identified, `Missing fraction` = missing_fraction,
            Original = rel_original_all_observed, New = rel_new_complete_case, `Median ratio, anchors` = rel_median_ratio_complete_case,
            `Median ratio, shared 50%` = rel_median_ratio_shared50, `Joint shared core` = rel_joint_shared_core, `Original minus new (raw)` = original_minus_new) |> sig3()
put(wb, "Sample offsets", off, d_cols = 8:13, int_cols = 5:6, n2_cols = 7, lim = 0.3, widths = c(12, 20, 8, 7, 10, rep(12, 8)))
put(wb, "Estimator agreement", sig3(conc |> mutate(dataset = lab_cmp(dataset), estimator_a = lab_est(estimator_a), estimator_b = lab_est(estimator_b))), widths = c(12, 44, 44, 12, 14, 14))
put(wb, "Held-out by statistic", sig3(hf |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 18, 20, 34, 50, 9, 14, 14, 14))
put(wb, "Held-out A755 A764", sig3(foc |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 26, 50, 18, 14, 14))
put(wb, "Depth within units", sig3(wu |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, 12, 14))
put(wb, "Intensity within units", sig3(iu |> mutate(dataset = lab_cmp(dataset))), widths = c(12, rep(22, 4)))
put(wb, "Hemispheres", sig3(rd(VAL, "hemisphere_concordance.csv") |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, 8, 16, 16, 16, 14, 50))
put(wb, "Residual variance", sig3(rv |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, rep(15, 4)))
put(wb, "Offsets by group", sig3(gp |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, rep(14, 6)))
put(wb, "Pre-centring QC", sig3(rd(VAL, "precentring_qc_variance_explained.csv") |>
  filter(response %in% c("raw_new_complete_case", "original_minus_new", "missing_fraction")) |>
  transmute(Compartment = lab_cmp(dataset), Response = response_description, Predictor = predictor, `R2 alone` = r2_alone,
            `Partial R2 given region/layer` = partial_r2_given_unit, Slope = slope_if_numeric)), widths = c(12, 48, 28, 12, 16, 12))
put(wb, "Anchor proteins", sig3(rd(VAL, "anchor_set_characteristics.csv") |> mutate(dataset = lab_cmp(dataset))), widths = c(12, rep(13, 15)))
put(wb, "Anchors by group", sig3(rd(VAL, "anchor_definition_by_group.csv") |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 26, 12, 16, 14, 18))
put(wb, "Anchor composition", sig3(rd(VAL, "anchor_functional_composition.csv") |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 12, 9, 9, 8, 12, 46, 10, 11, 10, 9))
put(wb, "Microglia carry-over", sig3(mc), left_cols = 1:2, widths = c(12, 70, 14, 14, 14, 14, 14, 14))
chk <- rd(VAL, "verification_additive_offsets.csv") |> mutate(dataset = lab_cmp(dataset), file = basename(file))
put(wb, "Checks", chk, widths = c(12, 80, 60, 18))
writeData(wb, "Checks", data.frame(x = c("Reproduction through config/normalisation.yml: the imputation rerun with each centring equals the archived (original, 20260601 and 20260526) and current (new, 20261019) workbooks value for value in all three compartments (reproduction/ in this review folder). Equality is of values, not files (dates are part of the file names).",
                                         "Animal-level inputs: archived = current + mean hemisphere offset change (2e-14).",
                                         sprintf("Scripted limma vs the archived ProTigy statistics on the identical archived input (validation/limma_vs_protigy.csv): %d contrasts, largest difference %.0e in logFC, %.0e in t, %.0e in adjusted P; FDR 5%% calls identical in %d of %d.",
                                                 nrow(lp), max(lp$max_abs_diff_logFC), max(lp$max_abs_diff_t), max(lp$max_abs_diff_adjP), sum(lp$same_FDR05_calls), nrow(lp)))),
          startRow = nrow(chk) + 3, colNames = FALSE)
dac <- rd(VAL, "da_counts_by_estimator.csv") |> select(dataset, unit, comparison, estimator, n_fdr05) |>
  mutate(estimator = factor(estimator, levels = names(EST_LABEL))) |> arrange(estimator) |> mutate(estimator = as.character(estimator)) |>
  pivot_wider(names_from = estimator, values_from = n_fdr05) |> mutate(dataset = lab_cmp(dataset)) |> arrange(dataset, comparison, unit)
names(dac)[-(1:3)] <- lab_est(names(dac)[-(1:3)])
put(wb, "DA consequences", dac, int_cols = 4:ncol(dac), widths = c(12, 14, 14, rep(20, ncol(dac) - 3)))
put(wb, "DA shift vs new", sig3(rd(VAL, "da_shift_vs_new_by_estimator.csv") |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), d_cols = 6, lim = 0.1, widths = c(12, 50, 30, 12, 14, 14, 14, 12, 12, 14))
put(wb, "CA2-SLM samples", sig3(cs), widths = c(18, 12, rep(12, 30)))
put(wb, "CA2-SLM variants", sig3(cv |> mutate(estimator = lab_est(estimator))), widths = c(50, 30, 14, 14, 14, 14, 14, 14))
put(wb, "CA2-SLM removal control", rd(VAL, "ca2slm_removal_control.csv") |> mutate(estimator = lab_est(estimator)), widths = c(50, 12, 40, 12, 18, 12, 12, 12))
put(wb, "CA2-SLM hits", sig3(hits), widths = c(10, 18, 12, 10, 10, 12, 12, 10, 10, 14, 14))
put(wb, "WGCNA stability", sig3(ws |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, rep(13, 14)))
put(wb, "WGCNA common power", sig3(wc |> mutate(dataset = lab_cmp(dataset), estimator = lab_est(estimator))), widths = c(12, 50, rep(13, 6)))
put(wb, "WGCNA paired animals", sig3(wp |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 14, rep(14, 5)))
att <- bind_rows(lapply(names(CMP), function(ds) rd(file.path("results", "wgcna", "attribute_wgcna_change", ds, "tables"), "attribution_steps.csv")))
put(wb, "WGCNA attribution", sig3(att |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 6, 60, 8, 10, 10, 10, 14, 14, 10))
# Run A (archived) used the reference-transcriptome background; its stamped background size is the length of the protein list, not the background used.
ew_show <- if (!is.null(ew)) sig3(ew |> mutate(dataset = lab_cmp(dataset),
  background_size_1 = ifelse(grepl("A vs", runs), NA, background_size_1),
  background_note = ifelse(grepl("A vs", runs), "run A: reference-transcriptome background (size not recorded)", "")))
put(wb, "EWCE", ew_show, widths = c(12, 12, 36, 8, 9, 12, 12, 12, 11, 11, 11, 11, 50))
saveWorkbook(wb, file.path(OUT, "Proteomics_normalisation_comparison.xlsx"), overwrite = TRUE)

# =================================================================================== 2. annotation decisions
# Both networks side by side, classified by the same scripted rules (tools/normalisation_review/label_modules_by_rule.R);
# only the wording of each label is written by hand (annotation_review/rule_based/labels.tsv, checked by check_rule_labels.R).
RB <- file.path(ANN, "rule_based")
rmod <- rd(RB, "module_rule_classes.csv"); rsm <- rd(RB, "supermodule_rule_classes.csv"); rsig <- rd(RB, "module_signals.csv")
# This workbook joins tables of the 2026-10-22 review by supermodule ID. The cut-height rule of 2026-10-26 merged and
# renumbered supermodules, so the workbook is not rewritten while the IDs differ (current labels and the ID map:
# Proteomics_WGCNA_module_labels.xlsx and annotation_review/label_changes_2026-10-23.csv).
sm22 <- read_csv(file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_review_2026-10-22",
                           "results/reviewer_audit/normalisation_choice_2026-10-21/annotation_review/rule_based/supermodule_rule_classes.csv"), show_col_types = FALSE)
sm_key <- function(d) sort(paste(d$version, d$dataset, d$SupermoduleID, d$member_modules))
ids_moved <- !identical(sm_key(rsm), sm_key(sm22))
rlab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = ""); rchk <- rd(RB, "label_form_check.csv")
rules <- yaml::read_yaml(file.path(RB, "rules.yml"))
disp <- rd(ANN, "displayed_labels.csv") |> select(dataset, entity_id, displayed_label)
VLAB <- c(original = "original centring (median of all observed values)", new = "complete-case centring (median of proteins quantified in every sample)")
old_labels <- bind_rows(lapply(names(CMP), function(ds) {
  s <- readRDS(file.path(if (dir.exists("Q:/data")) "Q:" else "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19",
                         "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"))
  a <- as.data.frame(s$WGCNA_modules_long) |> distinct(entity_id = ModuleID, auto = ModuleLabel_Final) |> mutate(dataset = ds, level = "module")
  f <- file.path("config", "wgcna_labels", "reference_pre_renormalisation_2026-10-19", paste0(ds, ".csv"))
  rv <- if (file.exists(f)) read_csv(f, show_col_types = FALSE) |> transmute(dataset, level, entity_id, reviewed = reviewed_biological_label) else NULL
  out <- a; if (!is.null(rv)) out <- full_join(out, rv, by = c("dataset", "level", "entity_id")) else out$reviewed <- NA_character_
  out })) |> mutate(shown = ifelse(!is.na(reviewed), paste0(reviewed, " (reviewed)"), auto))
old_sm_auto <- rsm |> filter(version == "original") |> transmute(dataset, level = "supermodule", entity_id = SupermoduleID, auto_sm = pipeline_label)
pctx <- function(x, d = 0) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f%%"), 100 * x))
# Label review of 2026-10-22, corrections applied to labels.tsv on 2026-10-23 (annotation_review/label_changes_2026-10-23.csv): one
# row per label; where the wording was corrected, the previous label and the reason; group and phenotype statistics are context
OA <- file.path(REV, "original_assessment")
sug <- rd(ANN, "label_review_suggestions.csv"); dli <- rd(ANN, "display_label_issues.csv"); f3o <- rd(ANN, "figure3_candidates_original.csv")
mb <- rd(OA, "modules_both_networks.csv"); sb <- rd(OA, "supermodules_both_networks.csv"); f3x <- rd(OA, "figure3_exemplar_check.csv")
SEV <- c(high = 1, medium = 2, low = 3)
sug_cols <- function(x, id) {   # x carries version, dataset and the entity id in column `id`
  if (is.null(sug)) return(x |> mutate(sug_label = "", sug_why = "", sug_sev = "", sug_with = "", sug_stab = ""))
  x |> left_join(sug |> transmute(version, dataset, !!id := entity_id, sug_label = ifelse(grepl("^applied", change), coalesce(previous_label, ""), ""),
                                  sug_why = ifelse(grepl("^applied", change), coalesce(reason, ""), ""), sug_sev = ifelse(grepl("^applied", change), coalesce(severity, ""), ""),
                                  sug_with = coalesce(depends_on_or_changes_with, ""), sug_stab = coalesce(stability_note, "")),
                 by = c("version", "dataset", id))
}
stab_flag <- function(conf, k) ifelse(!is.na(k) & k <= 3 & conf %in% c("high", "moderate"), sprintf("%s confidence but recovered in %d of 9", conf, as.integer(k)), "")
mod_sheet <- function(v) {
  x <- rmod |> filter(version == v) |> left_join(rlab |> filter(level == "module") |> select(version, dataset, ModuleID = entity_id, label), by = c("version", "dataset", "ModuleID"))
  shown <- if (v == "new") disp |> transmute(dataset, ModuleID = entity_id, shown = displayed_label) else
    old_labels |> filter(level == "module") |> transmute(dataset, ModuleID = entity_id, shown)
  x <- x |> left_join(shown, by = c("dataset", "ModuleID")) |> sug_cols("ModuleID") |>
    left_join(mb |> select(version, dataset, ModuleID, group_min_FDR_global, group_best_contrast, phenotype_min_FDR_family, phenotype_min_FDR_all),
              by = c("version", "dataset", "ModuleID"))
  x |> transmute(Compartment = lab_cmp(dataset), ID = ModuleID, Proteins = n_proteins, `Gene-mapped` = n_gene_mapped,
    `Signal class` = signal_class, Confidence = rule_confidence, `Current label` = label,
    `Previous label (before 2026-10-23)` = sug_label, `Why changed` = sug_why, Priority = sug_sev, `Changed with` = sug_with,
    `Dominant signal (top terms)` = ifelse(nzchar(coalesce(dominant_compartment, "")), paste0(dominant_signal, " [", dominant_compartment, "]"), coalesce(dominant_signal, "")),
    `Second signal` = coalesce(second_signal, ""), `Third signal` = coalesce(third_signal, ""),
    `Weak association` = ifelse(nzchar(coalesce(weak_association, "")), sprintf("%s (%s; %d proteins)", weak_association, weak_ontologies, as.integer(weak_size)), ""),
    `Signals from` = coalesce(signal_basis, ""), Coverage = pctx(coverage), `Hub support` = pctx(hub_support), `Hubs (top 10 kME)` = hubs,
    `Stability: recovered (of rebuilds)` = ifelse(is.na(stability_n_recovered), "", sprintf("%d of %d", as.integer(stability_n_recovered), as.integer(stability_n_rebuilds))),
    `Stability: median Jaccard` = round(stability_median_jaccard, 2), `Stability caveat` = stab_flag(rule_confidence, stability_n_recovered),
    `Group effect: lowest FDR` = group_min_FDR_global, `Group effect: best contrast` = coalesce(group_best_contrast, ""),
    `Phenotype: lowest FDR within family` = phenotype_min_FDR_family, `Phenotype: lowest FDR over all tests` = phenotype_min_FDR_all,
    Lineage = if (v == "new") coalesce(lineage_class, "") else "",
    `Corresponding module (share)` = ifelse(is.na(lineage_module), "", sprintf("%s %s (%s)", if (v == "new") "from original" else "to complete-case",
                                                                                 lineage_module, pctx(lineage_share))),
    `Label shown so far` = coalesce(shown, ""), `Reviewer decision` = "", `Reviewer label` = "", `Reviewer notes` = "") |>
    (\(d) if (v != "new") select(d, -Lineage) else d)()
}
sm_sheet <- function(v) {
  x <- rsm |> filter(version == v) |> left_join(rlab |> filter(level == "supermodule") |> select(version, dataset, SupermoduleID = entity_id, label), by = c("version", "dataset", "SupermoduleID"))
  shown <- if (v == "new") disp |> transmute(dataset, SupermoduleID = entity_id, shown = displayed_label) else
    old_labels |> filter(level == "supermodule") |> transmute(dataset, SupermoduleID = entity_id, shown) |>
      full_join(old_sm_auto |> transmute(dataset, SupermoduleID = entity_id, auto_sm), by = c("dataset", "SupermoduleID")) |>
      transmute(dataset, SupermoduleID, shown = coalesce(shown, auto_sm))
  x |> left_join(shown, by = c("dataset", "SupermoduleID")) |> sug_cols("SupermoduleID") |>
    left_join(sb |> select(version, dataset, SupermoduleID, group_min_FDR_global, group_best_contrast, phenotype_min_FDR_family, phenotype_min_FDR_all),
              by = c("version", "dataset", "SupermoduleID")) |>
    transmute(Compartment = lab_cmp(dataset), ID = SupermoduleID, `Members (proteins): class, confidence` = members,
              `Members' dominant signals` = member_dominant_signals, Class = rule_class, Confidence = rule_confidence, `Current label` = label,
              `Previous label (before 2026-10-23)` = sug_label, `Why changed` = sug_why, Priority = sug_sev, `Changed with` = sug_with,
              `Members' stability` = sug_stab,
              `Group effect: lowest FDR` = group_min_FDR_global, `Group effect: best contrast` = coalesce(group_best_contrast, ""),
              `Phenotype: lowest FDR within family` = phenotype_min_FDR_family, `Phenotype: lowest FDR over all tests` = phenotype_min_FDR_all,
              `Label shown so far` = coalesce(shown, ""), `Reviewer decision` = "", `Reviewer label` = "", `Reviewer notes` = "")
}
# review queue: every label corrected on 2026-10-23 and every stability caveat, primary (original) network first
queue <- if (is.null(sug)) NULL else sug |>
  left_join(rmod |> select(version, dataset, entity_id = ModuleID, rule_confidence, stability_n_recovered), by = c("version", "dataset", "entity_id")) |>
  mutate(caveat = ifelse(level == "module", stab_flag(rule_confidence, stability_n_recovered), ""),
         kind = case_when(grepl("^applied", change) & nzchar(caveat) ~ "label corrected + stability caveat", grepl("^applied", change) ~ "label corrected",
                          nzchar(caveat) ~ "stability caveat (label unchanged)", TRUE ~ NA_character_)) |>
  filter(!is.na(kind)) |>
  arrange(factor(version, c("original", "new")), coalesce(SEV[severity], 4), factor(dataset, names(CMP)), level, entity_id) |>
  transmute(`#` = row_number(), Network = ifelse(version == "original", "original (primary)", "complete-case (sensitivity)"), Compartment = lab_cmp(dataset),
            Level = level, ID = entity_id, What = kind, Priority = coalesce(severity, ""), `Label (corrected)` = current_label,
            `Previous label` = coalesce(previous_label, ""), Why = reason,
            Stability = ifelse(nzchar(caveat), caveat, coalesce(stability_note, "")), `Changes with` = coalesce(depends_on_or_changes_with, ""),
            `Decide in sheet` = paste(ifelse(level == "module", "Modules", "Supermodules"), ifelse(version == "original", "(original)", "(complete-case)")))
corr <- bind_rows(lapply(names(CMP), function(ds) read_csv(file.path("results", "wgcna", "map_wgcna_module_lineage", ds, "tables", "module_lineage_pairs.csv"), show_col_types = FALSE))) |>
  filter(share_of_old >= 0.2 | share_of_new >= 0.2) |>
  left_join(rlab |> filter(version == "original", level == "module") |> select(dataset, old_module = entity_id, original_label = label), by = c("dataset", "old_module")) |>
  left_join(rlab |> filter(version == "new", level == "module") |> select(dataset, new_module = entity_id, new_label = label), by = c("dataset", "new_module")) |>
  transmute(Compartment = lab_cmp(dataset), `Original module` = old_module, `Original label` = original_label, `Complete-case module` = new_module,
            `Complete-case label` = new_label, `Shared proteins` = n_shared, `Share of original` = pctx(share_of_old), `Share of complete-case` = pctx(share_of_new),
            Jaccard = round(jaccard, 2)) |> arrange(Compartment, `Original module`, desc(`Shared proteins`))
n_mod <- rmod |> count(version); n_sm <- rsm |> count(version)
cls_txt <- function(v) { x <- rmod |> filter(version == v) |> count(signal_class); paste(sprintf("%s %d", x$signal_class, x$n), collapse = ", ") }
stab_txt <- function(v) { x <- rmod |> filter(version == v, !is.na(stability_median_jaccard))
  sprintf("median Jaccard %s (median %.2f); recovered in at least 5 of 9 rebuilds: %d of %d modules", rng(x$stability_median_jaccard), median(x$stability_median_jaccard),
          sum(x$stability_n_recovered >= 5), nrow(x)) }
wb <- createWorkbook()
n_sug <- if (is.null(sug)) NULL else sug |> group_by(version) |> summarise(n = n(), wording = sum(grepl("^applied", change)), .groups = "drop")
n_cav <- if (is.null(queue)) 0 else sum(grepl("caveat", queue$What))
sug_txt <- if (is.null(n_sug)) "Label review suggestions: not yet available." else
  sprintf("Review of 2026-10-22: every one of the %d labels was re-checked against its own network's GO evidence and the wording rules, and an independent check confirmed or refined each correction. On 2026-10-23 the corrections were applied to the labels (annotation_review/rule_based/labels.tsv; log: annotation_review/label_changes_2026-10-23.csv): %d in the original network, %d in the complete-case network; the previous wording is in column Previous label. %d high- or moderate-confidence modules are recovered in 3 or fewer of 9 rebuilds (stability caveats; their labels stay).",
          sum(n_sug$n), n_sug$wording[n_sug$version == "original"], n_sug$wording[n_sug$version == "new"], n_cav)
add_readme(wb, "WGCNA modules and supermodules: label review (original centring primary)", list(
  list(heading = "Purpose", lines = c(
    sprintf("Primary analysis: the original-centring network (%d modules, %d supermodules); sensitivity analysis: the complete-case network (%d modules, %d supermodules) (decision of 2026-10-22). The same scripted rules classify every module and supermodule of both networks; only the wording of each label is written by hand. Each network's labels rest on its own GO evidence; module identities and labels are never carried between the networks (sheet Correspondence is for orientation only).",
            n_mod$n[n_mod$version == "original"], n_sm$n[n_sm$version == "original"], n_mod$n[n_mod$version == "new"], n_sm$n[n_sm$version == "new"]),
    sug_txt,
    "Labels name what the proteins have in common; they do not claim a demonstrated function, a cellular origin, stress responsiveness or causality. In microglia-enriched ROIs, synaptic and neuronal proteins are most likely neighbourhood proteins; their cellular origin is not resolved by these data, and labels must not imply microglial expression.")),
  list(heading = "How to review (suggested order)", lines = c(
    "1. Sheet Review queue: every corrected label and every stability caveat, original network first, then by priority. It only lists; decisions go into the module and supermodule sheets.",
    "2. Sheets Modules (original) and Supermodules (original): one row per label. Fill Reviewer decision for every row (accept, or 'accept with edited label' with your wording in Reviewer label; to return to a previous wording, copy it from Previous label). A supermodule changes with its members (column Changed with).",
    "3. Sheet Figure 3 candidates: the zoom module for each exemplar in the original network, with the complete-case choice for comparison. Selection rests on enrichment of the exemplar term and stability only.",
    "4. Then the complete-case sheets (sensitivity analysis) the same way. Sheet Display labels lists automatic pipeline display labels (both networks) that misstate their module or were carried between networks; figures and tables should use the reviewed labels instead.",
    "Nothing is promoted until the decisions are filled in; promotion to the label registry is a separate, recorded step.")),
  list(heading = "Column groups (kept separate on purpose)", lines = c(
    "Functional annotation: signal class, confidence, current label, previous label and why it changed, the signals with their top terms, coverage, hub support and hubs. Only this group determines a label.",
    "Network stability: recovered in k of 9 leave-one-animal-out rebuilds and median Jaccard; 'Stability caveat' marks high- or moderate-confidence modules recovered in 3 or fewer rebuilds. Report stability next to a label, never inside it.",
    "Group effect and phenotype: lowest FDR of the eigengene group test (all units) with its best contrast, and the lowest phenotype-association FDR within its family and over all tests (original network only; the complete-case network has the module-behaviour coupling test instead, 0 of 528 at FDR < 0.05). Context only: never used to choose or word a label.")),
  list(heading = "How a module is classified (scripted)", lines = c(
    sprintf("1. Evidence: the module's GO over-representation from the network build (BP, CC, MF; BH q < %.2f; background = all proteins tested in the compartment, identical for both networks).", rules$q),
    sprintf("2. Signals: the %d most significant BP terms are grouped by shared proteins (overlap coefficient >= %.1f, the EnrichmentMap default). Groups whose proteins are mostly the same, or that sit mostly (at least half) in one significant CC term of the module (e.g. mitochondrial inner membrane), are one signal. A module without a passing BP signal is classified on its CC and MF terms (a location- or complex-defined module; column 'Signals from').",
            rules$top_terms, rules$overlap),
    sprintf("3. Label threshold: the module needs significant terms in at least %d ontologies, and a signal passes when its leading term is carried by more than %d proteins. Without a passing signal the module is 'unannotated'; with exactly one (small) signal, 'unannotated (weak ... association)'.", rules$min_ontologies, rules$min_proteins),
    sprintf("4. Signal class: the largest passing signal is the dominant one; every passing signal at least %.0f%% its size is comparable. single signal / two signals (both named) / several signals.", 100 * rules$comparable_ratio),
    sprintf("5. Confidence: coverage = proteins in the dominant signal's largest term / gene-mapped module proteins; hub support = share of the 10 highest-kME hubs in the dominant signal. Single signal: high with coverage >= %.0f%% and hub support >= %.0f%%; moderate with coverage >= %.0f%% (or the module holds at least half of a term with at least %d proteins); low otherwise. Two signals: moderate or low by the same coverage rule. Several signals and weak associations: low. No signal: none.",
            100 * rules$high_cov, 100 * rules$hub_share, 100 * rules$mod_cov, rules$complex_n),
    sprintf("6. Lineage (complete-case network): carried when one original module supplies at least %.0f%% of the proteins and its dominant signal is found again; revised when it supplies that share but the signal changed; new below that share.", 100 * rules$carry_share),
    "7. Supermodules: singleton (one member; its class and confidence), unannotated (no member labelled), mixed (one member unannotated, another labelled), concordant (members share the same dominant signal), compound (members with different signals). Confidence: the lowest member's; mixed: low.",
    sprintf("8. Stability (descriptive, not a filter): each module is matched in the 9 leave-one-animal-out rebuilds of its network (best Jaccard per rebuild); only rebuilds that reproduce the network exactly are used (ARI >= %.2f; for the original neuropil network the rebuilds start from its archived input). Original: %s. Complete-case: %s.",
            rules$rebuild_ari, stab_txt("original"), stab_txt("new")))),
  list(heading = "How a label is worded (by hand, checked)", lines = c(
    "single signal: 'X-enriched' at high confidence, 'X-associated' otherwise; two signals: 'X- and Y-associated'; several signals: 'Several signals: X, Y, Z'; 'unannotated (weak X association)'; 'unannotated'.",
    "Supermodules: a singleton carries its member's label with the ID; 'Compound: X (mA) + Y (mB)'; 'Mixed: X (mA) + unannotated member (mB)'; concordant: the shared label with both IDs. A member's part uses the member's own wording: its dominant signal, plus its second signal when it has two or several.",
    "Terms keep their meaning when shortened (a 'regulation of' term stays regulation; a process is not renamed after a structure), signals are named after their leading or largest term, and a module carried between the networks with the same signals carries the same words.",
    if (!is.null(rchk)) sprintf("Check: the form of all %d labels follows from their class and confidence, and every supermodule label names all members (%d of %d pass; sheet Label form check).", nrow(rchk), sum(rchk$problem == "" | is.na(rchk$problem)), nrow(rchk)) else "Label form check pending.")),
  list(heading = "Results", lines = c(
    sprintf("Original network: %s. Complete-case network: %s.", cls_txt("original"), cls_txt("new")),
    "Sheet Correspondence pairs the modules of the two networks (shared proteins, share of each, Jaccard) with both labels.")),
  list(heading = "Columns", lines = c(
    "Dominant, second and third signal: their most significant terms, with the shared CC term in brackets where the signal was grouped by compartment. Label shown so far: original network = the reviewed label where one exists (microglia; neuropil m11), otherwise the automatic label; complete-case network = the label registry's canonical display label.")),
  list(heading = "Reviewer decision values", lines = c("accept | accept with edited label | reject (keep unannotated) | needs discussion. Put an edited label in Reviewer label.")),
  list(heading = "Zoom figure", lines = c(
    "Sheet Figure 3 candidates gives the zoom module for each exemplar in both networks, original (primary) first. The sheets named Zoom ... below are the detailed checks for the complete-case network (sensitivity analysis).",
    "Sheet Zoom choices re-checks the exemplar-aligned modules of the complete-case network (neuropil m04, soma m01, microglia m09) against two rules (most significant BP terms in the exemplar's theme; strongest enrichment of the exemplar term) and against how the exemplar term's GSEA leading edge distributes over modules (contrasts where the term is significant; a module without leading-edge proteins counts as 0%).",
    "The previous zoom used three neuropil modules of the original network (m01, m02, m12); the RNA-processing and respiration exemplars now sit in soma and microglia, so the two choices are not like-for-like. The original choices are kept; nothing is swapped."))))
dv <- '"accept,accept with edited label,reject (keep unannotated),needs discussion"'
W_MOD <- c(Compartment = 11, ID = 11, Proteins = 9, `Gene-mapped` = 9, `Signal class` = 13, Confidence = 11, `Current label` = 44,
           `Previous label (before 2026-10-23)` = 44, `Why changed` = 60, Priority = 9, `Changed with` = 22, `Dominant signal (top terms)` = 60,
           `Second signal` = 44, `Third signal` = 36, `Weak association` = 40, `Signals from` = 16, Coverage = 9, `Hub support` = 9, `Hubs (top 10 kME)` = 40,
           `Stability: recovered (of rebuilds)` = 12, `Stability: median Jaccard` = 11, `Stability caveat` = 24, `Members' stability` = 30,
           `Group effect: lowest FDR` = 11, `Group effect: best contrast` = 18, `Phenotype: lowest FDR within family` = 12, `Phenotype: lowest FDR over all tests` = 12,
           Lineage = 12, `Corresponding module (share)` = 30, `Members (proteins): class, confidence` = 50, `Members' dominant signals` = 50, Class = 13,
           `Label shown so far` = 36, `Reviewer decision` = 18, `Reviewer label` = 30, `Reviewer notes` = 30)
NUMERIC <- c("Proteins", "Gene-mapped", "Stability: median Jaccard", "Group effect: lowest FDR", "Phenotype: lowest FDR within family", "Phenotype: lowest FDR over all tests")
put_entity <- function(df, sn, freeze) {
  pc <- col_idx(df, intersect(c("Group effect: lowest FDR", "Phenotype: lowest FDR within family", "Phenotype: lowest FDR over all tests"), names(df)))
  put(wb, sn, df, int_cols = col_idx(df, intersect(c("Proteins", "Gene-mapped"), names(df))), p_cols = pc, sig_cols = pc,
      n2_cols = col_idx(df, intersect("Stability: median Jaccard", names(df))), left_cols = which(!names(df) %in% c(NUMERIC, "Compartment", "ID", "Priority")),
      widths = unname(W_MOD[names(df)]), freeze_col = freeze)
  dataValidation(wb, sn, cols = which(names(df) == "Reviewer decision"), rows = 2:(nrow(df) + 1), type = "list", value = dv)
  addStyle(wb, sn, WRAP, rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
}
if (!is.null(queue)) {
  put(wb, "Review queue", queue, int_cols = 1, left_cols = c(2:6, 8:13), widths = c(5, 22, 11, 12, 11, 24, 9, 44, 44, 70, 30, 24, 24), freeze_col = 6)
  addStyle(wb, "Review queue", WRAP, rows = 2:(nrow(queue) + 1), cols = seq_len(ncol(queue)), gridExpand = TRUE, stack = TRUE)
}
for (v in c("original", "new")) {
  put_entity(mod_sheet(v), if (v == "new") "Modules (complete-case)" else "Modules (original)", 3)
  put_entity(sm_sheet(v), if (v == "new") "Supermodules (complete-case)" else "Supermodules (original)", 3)
}
# Figure 3 zoom candidates: original network (primary) and complete-case network (sensitivity), by exemplar-term enrichment and stability
if (!is.null(f3x)) {
  pick_new <- c(neuron_neuropil = "WGCNA_m04", neuron_soma = "WGCNA_m01", microglia = "WGCNA_m09")
  pick_orig <- if (!is.null(f3o)) setNames(f3o$module, f3o$dataset) else c(neuron_neuropil = "WGCNA_m01", neuron_soma = "WGCNA_m02", microglia = "WGCNA_m01")
  f3 <- f3x |> filter((version == "original" & ModuleID == pick_orig[dataset]) | (version == "new" & ModuleID == pick_new[dataset])) |>
    left_join(mb |> select(version, dataset, ModuleID, group_min_FDR_global), by = c("version", "dataset", "ModuleID")) |>
    left_join(if (!is.null(f3o)) f3o |> transmute(dataset, version = "original", alternatives, caveat) else tibble(dataset = character(), version = character(), alternatives = character(), caveat = character()),
              by = c("dataset", "version")) |>
    arrange(factor(dataset, names(CMP)), factor(version, c("original", "new"))) |>
    transmute(Compartment = lab_cmp(dataset), Network = ifelse(version == "original", "original (primary)", "complete-case (sensitivity)"), Exemplar = exemplar,
              Module = ModuleID, Label = label, `Signal class` = signal_class, Confidence = rule_confidence,
              `Exemplar term: q` = exemplar_padj, `Exemplar term: overlap` = exemplar_overlap, `Exemplar term: fold` = round(exemplar_fold, 2),
              `Theme terms (significant BP)` = n_sig_BP_in_theme, `Stability: recovered (of 9)` = stability_n_recovered, `Stability: median Jaccard` = round(stability_median_jaccard, 2),
              `Group effect: lowest FDR` = group_min_FDR_global, Alternatives = coalesce(alternatives, ""), Caveat = coalesce(caveat, ""))
  put(wb, "Figure 3 candidates", f3, p_cols = c(8, 14), int_cols = c(11, 12), n2_cols = c(10, 13), left_cols = c(1:7, 9, 15, 16),
      widths = c(11, 22, 34, 11, 44, 13, 11, 10, 12, 9, 10, 10, 10, 11, 60, 60), freeze_col = 3)
  addStyle(wb, "Figure 3 candidates", WRAP, rows = 2:(nrow(f3) + 1), cols = seq_len(ncol(f3)), gridExpand = TRUE, stack = TRUE)
}
if (!is.null(dli)) {
  put(wb, "Display labels", dli |> transmute(Network = ifelse(network == "original", "original (primary)", "complete-case (sensitivity)"), Compartment = lab_cmp(dataset),
                                            ID = entity_id, `Display label (pipeline)` = displayed_label, Issue = issue, `Use instead` = use_instead),
      left_cols = 1:6, widths = c(22, 11, 11, 50, 60, 50))
  addStyle(wb, "Display labels", WRAP, rows = 2:(nrow(dli) + 1), cols = 1:6, gridExpand = TRUE, stack = TRUE)
}
put(wb, "Correspondence", corr, int_cols = 6, left_cols = c(2:5), widths = c(11, 11, 50, 13, 50, 10, 10, 12, 9))
addStyle(wb, "Correspondence", WRAP, rows = 2:(nrow(corr) + 1), cols = c(3, 5), gridExpand = TRUE, stack = TRUE)
put(wb, "Signals", sig3(rsig |> mutate(dataset = lab_cmp(dataset))), left_cols = c(1:5, 9, 11, 17, 18), widths = c(10, 11, 11, 24, 8, 9, 11, 34, 9, 10, 9, 34, 9, 10, 8, 12, 80, 10))
if (!is.null(rchk)) put(wb, "Label form check", rchk |> mutate(dataset = lab_cmp(dataset)), left_cols = 1:8, widths = c(10, 11, 12, 11, 16, 11, 70, 30))
zt <- bind_rows(lapply(names(CMP), function(ds) rd(ANN, sprintf("zoom_exemplar_term_by_module_%s.csv", ds)) |> mutate(dataset = ds)))
zl <- rd(ANN, "zoom_exemplar_leading_edge_by_module.csv")
pick <- function(ds, m) zt |> filter(dataset == ds, ModuleID == m)
share <- function(ds, m) {   # over the contrasts where the exemplar term is significant; a module absent from a leading edge counts as 0%
  sig <- unique(zl$contrast[zl$dataset == ds & zl$p.adjust < 0.05]); if (!length(sig)) return("0")
  rng(100 * vapply(sig, function(ct) { v <- zl$share_of_leading_edge[zl$dataset == ds & zl$contrast == ct & zl$module == m]; if (length(v)) v[1] else 0 }, 0), 0) }
fmt_term <- function(ds, m, extra = "") { x <- pick(ds, m); sprintf("%sq %.1e; %s; %.2f-fold%s", if (nzchar(extra)) "" else "", x$exemplar_padj, x$exemplar_overlap, x$exemplar_fold, extra) }
zoom <- tibble(
  Compartment = c("Neuropil", "Soma", "Microglia"),
  `Exemplar term (unit; theme)` = c("synaptic signaling GO:0099536 (CA3 sr; synaptic_signaling_vesicle)", "mRNA processing GO:0006397 (CA2 sp; rna_processing_splicing_rnp)", "oxidative phosphorylation GO:0006119 (CA1; mitochondrial_respiration_oxphos)"),
  `Chosen module` = c("WGCNA_m04", "WGCNA_m01", "WGCNA_m09"),
  `In-theme significant BP terms (chosen; next)` = c(sprintf("%d; m10 %d", pick("neuron_neuropil", "WGCNA_m04")$n_sig_BP_terms_in_theme, pick("neuron_neuropil", "WGCNA_m10")$n_sig_BP_terms_in_theme),
                                                    sprintf("%d; no other module", pick("neuron_soma", "WGCNA_m01")$n_sig_BP_terms_in_theme),
                                                    sprintf("%d; m16 %d", pick("microglia", "WGCNA_m09")$n_sig_BP_terms_in_theme, pick("microglia", "WGCNA_m16")$n_sig_BP_terms_in_theme)),
  `Exemplar term in chosen module (q; overlap; fold)` = c(fmt_term("neuron_neuropil", "WGCNA_m04"), fmt_term("neuron_soma", "WGCNA_m01"), fmt_term("microglia", "WGCNA_m09")),
  `Next module for the term` = c(paste("m10", fmt_term("neuron_neuropil", "WGCNA_m10")), "none (q > 0.98)", paste("m16", fmt_term("microglia", "WGCNA_m16", " (same SM06)"))),
  `GSEA leading edge (significant contrasts)` = c(sprintf("m04 %s%%; m10 %s%%", share("neuron_neuropil", "WGCNA_m04"), share("neuron_neuropil", "WGCNA_m10")),
                                                  sprintf("m01 %s%%", share("neuron_soma", "WGCNA_m01")),
                                                  sprintf("m09 %s%%; m16 %s%%; m02 %s%%", share("microglia", "WGCNA_m09"), share("microglia", "WGCNA_m16"), share("microglia", "WGCNA_m02"))),
  Assessment = c("Both rules pick m04 (more in-theme terms; smallest q), but m04 and m10 have equivalent fold enrichment and the leading edge splits between them about in proportion to where the term's proteins sit in the network: it does not single out either module. Show m04 with m10 as context, or name both. Not swapped.",
                 "Unambiguous by both rules and the leading edge.",
                 "Both rules pick m09; m09 holds more of the leading edge than of the term's proteins, so the support is real. SM06 (m09 + m16) is the fuller unit for this term."))
put(wb, "Zoom choices", zoom, left_cols = 1:8, widths = c(12, 44, 12, 22, 34, 34, 34, 70))
addStyle(wb, "Zoom choices", WRAP, rows = 2:4, cols = 1:8, gridExpand = TRUE, stack = TRUE)
put(wb, "Zoom previous choices", rd(ANN, "zoom_previous_choices.csv"), left_cols = 1:3, widths = c(14, 10, 80, 24, 10, 12, 10))
put(wb, "Zoom leading edge", sig3(zl |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 14, 22, 12, 12, 30, 14, 9, 11, 10, 12, 12, 12, 12, 12))
put(wb, "Zoom term by module", sig3(zt |> mutate(dataset = lab_cmp(dataset))), widths = c(12, 14, 12, 12, 12, 12, 12, 12, 30))
if (ids_moved) message("Proteomics_WGCNA_annotation_decisions.xlsx not rewritten: supermodule IDs differ from those of 2026-10-22") else
  saveWorkbook(wb, file.path(OUT, "Proteomics_WGCNA_annotation_decisions.xlsx"), overwrite = TRUE)

# =================================================================================== 3. provenance inventory
inv <- rd(INV, "normalisation_dependence_inventory.csv"); fam <- rd(INV, "analysis_families.csv"); cov <- rd(INV, "ledger_coverage.csv")
# Stale rows of the inventory, with the reason recorded in their status; reads outside the ledger
stale_rows <- if (!is.null(inv)) grepl("^stale", inv$status) else logical()
stale_reasons <- if (!is.null(inv)) paste(unique(sprintf("%s (%s)", sub("[.][Rr]$", "", basename(inv$script[stale_rows])),
  sub("^stale( [(]audit[)])?: ", "", inv$status[stale_rows]))), collapse = "; ") else ""
stale_outside_ledger <- if (!is.null(inv)) paste(unique(sub("[.][Rr]$", "", basename(inv$script[stale_rows & coalesce(inv$ledger_inputs, 0L) == 0]))), collapse = ", ") else ""
# Archive and test evidence, read from the files the Checks sheet cites
# absolute: the repository is usually opened through subst P:, where ".." of the root is the root itself
ANALYSIS_DIR <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
ARCH <- setNames(file.path(ANALYSIS_DIR, c("_archive_proteomics_pre_renormalisation_2026-10-19", "_archive_proteomics_complete_case_2026-10-21",
                                           "_archive_proteomics_pride_submission_2026-10-21")), c("original", "new", "pride"))
stopifnot(all(file.exists(file.path(ARCH[c("new", "pride")], "SHA256_MANIFEST.csv"))))
man_status <- function(a) { f <- file.path(a, "SHA256_MANIFEST.csv"); if (!file.exists(f)) return(NULL)
  x <- read_csv(f, show_col_types = FALSE, col_select = any_of("status")); if ("status" %in% names(x)) count(x, status) else tibble(status = "listed", n = nrow(x)) }
st_new <- man_status(ARCH["new"]); st_pride <- man_status(ARCH["pride"])
code_line <- grep("^Code:", readLines(file.path(ARCH["new"], "ARCHIVE_SUMMARY.txt"), warn = FALSE), value = TRUE)
code_state <- if (length(code_line) && grepl("[(]refreshed ", code_line)) sub("^.*[(](refreshed [^;)]*).*$", "\\1", code_line) else "see ARCHIVE_SUMMARY.txt"
cnt_of <- function(s, k) if (is.null(s)) NA_integer_ else sum(s$n[s$status %in% k])
spot <- rd(VAL, "archive_spot_check_2026-10-21.csv"); repro <- rd(file.path(REV, "reproduction"), "reproduction_check.csv")
tst <- if (file.exists(file.path(VAL, "test_run_2026-10-21.txt"))) readLines(file.path(VAL, "test_run_2026-10-21.txt")) else character()
tst <- tst[grepl("tests=", tst)]; tst_n <- sum(as.integer(sub(".*tests= *([0-9]+).*", "\\1", tst))); tst_f <- sum(as.integer(sub(".*failed=([0-9]+).*", "\\1", tst)))
n_guarded <- 18L   # files calling assert_normalisation_decided(): 18 export, freeze and deposit entry points (listed in R/data_contracts/normalisation_contract.R)
wb <- createWorkbook()
add_readme(wb, "Proteomics: analyses against the two centring versions", list(
  list(heading = "Versions", lines = c(
    sprintf("Original: median of all observed values per sample. Outputs preserved and hash-verified in Analysis/_archive_proteomics_pre_renormalisation_2026-10-19: 455,556 files, 38 GB (35.6 GiB) of exports/, data/processed and results/ (spot-checked again 2026-10-21: %s random files re-hashed, %s mismatches). Its pride_submission/ tree (not in that archive) is preserved separately in Analysis/_archive_proteomics_pride_submission_2026-10-21 (%s files, %s match). Code pRoteomics b478a13; differential model fitted in ProTigy.",
            fmt(nrow(spot)), fmt(sum(spot$status != "match")), fmt(cnt_of(st_pride, c("match", "mismatch", "source_error", "listed"))), fmt(cnt_of(st_pride, "match"))),
    sprintf("New: median of the proteins quantified in every sample of the compartment. Analysis outputs in the working tree (branch renormalisation-2026-10-19, uncommitted), preserved 2026-10-21 as a hash-verified copy in Analysis/_archive_proteomics_complete_case_2026-10-21: %s files, %s match; %s Excel lock files were copied but could not be verified. The code state (patch against b478a13 and the untracked code files) is in its code/ folder (%s). This review folder is not part of that copy.",
            fmt(sum(st_new$n)), fmt(cnt_of(st_new, "match")), fmt(sum(st_new$n) - cnt_of(st_new, "match")), code_state),
    "The frozen export bundle (exports/, results/manuscript, pride_submission) in the working tree still holds the original version; Exp9_manuscript currently imports from it. The import tools in Exp9_manuscript are not guarded (other repository).",
    "Legacy original-centring trees are still in the working tree (results/tables, data/processed/02-10, gct_extractR, protigy_output). The guard checks the imputed workbooks, the metadata-joined workbook, the animal-level GCT and the limma fit; it does not check these trees or the export bundle.",
    "One setting selects the version: config/normalisation.yml (sample_centring, decision, approved_sample_centring). The guards are listed in sheet Checks. A version is regenerated in a separate working tree at the same code (config/normalisation_versions.yml); regeneration in place is guarded but not the recommended route.")),
  list(heading = "Inventory columns", lines = c(
    "One row per step and dataset: every pipeline.yml step, the legacy preprocessing chain (impute, metadata join, animal-level input), and the steps outside the registry (tools/ exporters and review tools, export trees, Exp9_manuscript, discovery audits). writes_into_legacy_tree marks new-version files written into an original-centring tree.",
    "Outputs: n_outputs_archive_same_path counts the archived files at the same path; archive_has_same_named_outputs counts the step's current file names found anywhere in the archive manifest (most archived originals sit at legacy paths).",
    "ran_before_upstream_regenerated: the step's oldest output predates the regeneration of what it reads (network, labels, differential tables, enrichment), so it ran on previous inputs. legacy_paths_in_code: hard-coded legacy result paths in the script (may be fallbacks).",
    if (!is.null(cov)) sprintf("Ledger columns use the rerun window only and the last read of each input by dataset. The ledger covers %d of %d centring-dependent steps; %d of %d rerun rows carry no script ID.",
            cov$dependent_steps_with_ledger_reads, cov$dependent_steps, cov$unattributed_rows, cov$ledger_rows_in_rerun) else "Ledger coverage pending.",
    "original_reproducible / new_reproducible: whether each version can be regenerated through the setting, and whether the code changed since b478a13 (then the archived outputs are not byte-reproducible).")),
  list(heading = "Special cases", lines = c(
    paste("EWCE.", txt_ewce),
    "Neuropil overlap scores: the curated overlap programmes were defined from proteins recurring among the top GO-driving proteins of the original differential results. Scoring them on the new data, or regenerating them from the new results, is outcome-selected and cannot serve as confirmation (sheet Decisions)."))))
if (!is.null(inv)) put(wb, "Inventory", inv |> select(-any_of("notes")), left_cols = c(1, 2, 6, 8, 13:15, 20:23, 25, 26),
    widths = c(16, 52, 14, 9, 11, 48, 11, 16, 11, 18, 18, 11, 22, 34, 24, 18, 9, 9, 11, 40, 40, 24, 50, 11, 40, 50), freeze_col = 3)
if (!is.null(inv)) { summ <- inv |> count(stage, status) |> pivot_wider(names_from = status, values_from = n, values_fill = 0)
  put(wb, "Inventory summary", summ, int_cols = 2:ncol(summ), widths = c(30, rep(18, ncol(summ) - 1))) }
put(wb, "Analysis families", fam, left_cols = c(1, 3, 4), widths = c(34, 9, 90, 80))
checks <- tibble(
  Check = c("Imputed matrices, original setting vs archived workbooks", "Imputed matrices, new setting vs working tree",
            "Archived minus current matrices = one constant per sample", "Animal-level inputs: archived = current + offset change",
            "Scripted limma vs archived ProTigy statistics", "Back-filled stamp vs hashes recorded on 2026-10-19",
            "Archive unchanged", "Launcher: setting vs stamp mismatch", "Launcher: inputs not made from the stamped workbooks",
            "Freeze/export/import while decision pending", "Approval tied to one centring", "Regeneration over an unpreserved version",
            "Imputed workbook from another run in the same directory", "Last read of every ledgered input in the rerun",
            "Steps run before their inputs were regenerated (ordering check)", "Static tests touching the changed code",
            "Snapshot of the new version"),
  Result = c(if (!is.null(repro)) sprintf("%s: %d of %d compartments identical in every value and annotation (run %s with the final guarded code)", "original setting", sum(repro$identical_values & repro$ann_identical & grepl("^original", repro$check)), sum(grepl("^original", repro$check)), max(repro$run_at)) else "pending",
             if (!is.null(repro)) sprintf("%s: %d of %d compartments identical", "new setting", sum(repro$identical_values & repro$ann_identical & grepl("^new", repro$check)), sum(grepl("^new", repro$check))) else "pending",
             "max spread 6e-15", "max difference 2e-14",
             sprintf("%d contrasts: largest difference %.0e (logFC), %.0e (t); FDR 5%% calls identical in %d of %d", nrow(lp), max(lp$max_abs_diff_logFC), max(lp$max_abs_diff_t), sum(lp$same_FDR05_calls), nrow(lp)),
             "equal (animal-level builder source hashes)",
             sprintf("%s random files re-hashed: %s mismatches (2026-10-21)", fmt(nrow(spot)), fmt(sum(spot$status != "match"))),
             "refused; nothing run",
             "refused when the newest metadata-joined workbook's name does not match a stamped workbook (a name check only: a joined workbook built from an earlier same-day imputation is not detected), when the source hash recorded in aggregation_summary.csv differs from the stamp (it fails when that summary is missing but a GCT exists), or when the limma fit's recorded input md5 differs from the GCT on disk. The GCT itself is not hashed against the summary",
             sprintf("%d entry points stop; the dry runs of 08, 09, RUN_EXPORT and the migrate preview (DRY=1) only report. The Exp9_manuscript import tools are not guarded.", n_guarded),
             "decision approved AND approved_sample_centring equal to the setting in the file; the environment override never satisfies it",
             "refused unless the stamped (or, for an unstamped directory, the on-disk) workbook hashes are a preserved set in config/normalisation_versions.yml with a manifest; the check does not test that the preserved copy is complete, and a directory holding imputed workbooks but neither a stamp nor imputation_qc.csv is not checked",
             "workbooks with another date are refused before any workbook is written; same-day workbooks are not refused",
             if (!is.null(cov)) sprintf("within the ledger (%d of %d centring-dependent steps write ledger reads), only the curated overlap programmes are an original-centring input (score_module_activity; decision needed). Original-centring reads outside the ledger: %s (see Inventory status)", cov$dependent_steps_with_ledger_reads, cov$dependent_steps, stale_outside_ledger) else "pending",
             if (!is.null(inv)) sprintf("%d step x dataset rows are stale (Inventory status 'stale'): %d pipeline rows and %d discovery audits. All were identified by review; the ordering check (ran_before_upstream_regenerated) flagged %d of them, because it misses steps without code markers for their inputs.", sum(stale_rows), sum(stale_rows & !grepl("discovery audits", inv$stage)), sum(stale_rows & grepl("discovery audits", inv$stage)), sum(stale_rows & nzchar(coalesce(as.character(inv$ran_before_upstream_regenerated), "")) & !inv$ran_before_upstream_regenerated %in% c("FALSE", "NA"))) else "pending",
             sprintf("%d suites, %d tests, %d failures (2026-10-21)", length(tst), tst_n, tst_f),
             sprintf("%s files: %s match; %s Excel lock files copied but not verifiable", fmt(sum(st_new$n)), fmt(cnt_of(st_new, "match")), fmt(sum(st_new$n) - cnt_of(st_new, "match")))),
  Evidence = c("reproduction/reproduction_check.csv", "reproduction/reproduction_check.csv", "validation/verification_additive_offsets.csv",
               "validation/verification_additive_offsets.csv", "validation/limma_vs_protigy.csv", "data/processed/01_preprocessing/impute/normalisation_stamp.yml",
               "validation/archive_spot_check_2026-10-21.csv", "R/data_contracts/normalisation_contract.R", "R/data_contracts/normalisation_contract.R",
               "assert_normalisation_decided()", "config/normalisation.yml", "assert_working_version_preserved()", "archive/01_preprocessing/01_impute.r",
               "inventory/last_reads_of_pre_rerun_derived_inputs.csv", "inventory/normalisation_dependence_inventory.csv", "validation/test_run_2026-10-21.txt",
               "Analysis/_archive_proteomics_complete_case_2026-10-21/SHA256_MANIFEST.csv"))
put(wb, "Checks", checks, left_cols = 1:3, widths = c(56, 80, 70))
addStyle(wb, "Checks", WRAP, rows = 2:(nrow(checks) + 1), cols = 1:3, gridExpand = TRUE, stack = TRUE)
dec <- tibble(
  Item = c("Normalisation", "Within the anchor family", "CA2-SLM protein results", "Module and supermodule labels", "Zoom figure",
           "Neuropil overlap scores", "EWCE background", "Stale steps found by review", "Legacy original-centring trees in the working tree",
           "Preserve both versions", "Exports, Exp9_manuscript import and the import tools", "DIA-NN search details"),
  `Decision needed` = c("original | new (complete-case anchor median) | neither conclusively favoured",
                        "keep the anchor median, or adopt the median ratio on the same anchors",
                        "report with a sensitivity analysis without A755 L (and A764 R), restrict claims, or use a missingness-aware model",
                        if (ids_moved) "review the current labels in Proteomics_WGCNA_module_labels.xlsx (supermodule IDs of 2026-10-22 are mapped in annotation_review/label_changes_2026-10-23.csv)" else
                          "accept / edit / reject each row of Proteomics_WGCNA_annotation_decisions.xlsx (both networks)",
                        "keep the exemplar-aligned modules (neuropil m04 with m10 context, soma m01, microglia m09/SM06), or choose other exemplars",
                        "retire, or keep as a predefined descriptive analysis labelled as defined on the original results",
                        "measured-proteome background as primary (the archived reference-transcriptome tables would then be historical), or keep the reference-transcriptome background",
                        if (!is.null(inv)) paste0("rerun after the decision, or retire: ", stale_reasons) else "pending",
                        "move the pre-rerun derived trees (results/tables, data/processed/02-10, gct_extractR, protigy_output) out of the working tree (they are in the archive), or make the resolvers strict, before any further rerun or export. First relocate the new-version files written into them (summarize_ca2_slm_robustness; the animal-level aggregation tables in results/tables/01_preprocessing, which the guard reads)",
                        "no decision needed (done 2026-10-21): the complete-case copy (Analysis/_archive_proteomics_complete_case_2026-10-21) and the original pride_submission/ tree (Analysis/_archive_proteomics_pride_submission_2026-10-21) are hash-verified and registered in config/normalisation_versions.yml",
                        "set decision: approved and approved_sample_centring; then rework export, freeze and deposit (the export config still packages the legacy original-centring DA trees), and add a centring check to the Exp9_manuscript import tools",
                        "ask the facility for report.log.txt (version, normalisation mode, joint search, final vs quick search)"),
  Recommendation = c("recommendation only, decision pending: new over original for proteins quantified in every sample, with the qualifications in Proteomics_normalisation_comparison.xlsx (Read me, Recommendation)",
                     "either; the data do not separate them", "sensitivity analysis required before any claim",
                     "see the proposals: classes, confidence and stability come from the same scripted rules for both networks; only the wording is hand-written (form-checked)",
                     "keep, with the m10 note", "retire as confirmation; at most a historical descriptive record",
                     sprintf("proteome background (the change from the archived EWCE comes with the step from run A to run C, which changes the background, the species space and the code together; the centring changes the scores little, but %s%% of the significant differential calls are not shared)", rng(100 * (1 - ew_c$jac), 0)),
                     "repoint the hard-coded reads and remove the EWCE fallbacks of test_microglia_targeted_signatures, then rerun with the chosen version; rerun the discovery audits only if they are to be cited", "move them (reversible, with a hash record)", "none",
                     "after all of the above", "request it"))
put(wb, "Decisions", dec, left_cols = 1:3, widths = c(34, 80, 70))
addStyle(wb, "Decisions", WRAP, rows = 2:(nrow(dec) + 1), cols = 1:3, gridExpand = TRUE, stack = TRUE)
st <- rd(INV, "statements_review.csv")
put(wb, "Statements", st, left_cols = seq_len(ncol(st)), widths = rep(30, ncol(st)))
if (!is.null(st)) addStyle(wb, "Statements", WRAP, rows = 2:(nrow(st) + 1), cols = seq_len(ncol(st)), gridExpand = TRUE, stack = TRUE)
put(wb, "EWCE", ew_show, widths = c(12, 12, 36, 8, 9, 12, 12, 12, 11, 11, 11, 11, 50))
put(wb, "Pre-rerun inputs read", rd(INV, "last_reads_of_pre_rerun_derived_inputs.csv"), left_cols = 1:5, widths = c(30, 16, 30, 24, 100))
put(wb, "Ledger coverage", cov, widths = rep(22, 4))
saveWorkbook(wb, file.path(OUT, "Proteomics_normalisation_provenance.xlsx"), overwrite = TRUE)
cat("written:", OUT, "\n")
