#!/usr/bin/env Rscript
#
# CA2-SLM DAP robustness audit.
#
# THE QUESTION
#   Which of the 28 CA2-SLM SUS-vs-RES DAPs remain supported independently of
#   differential missingness, imputation, and the two QC-failed SUS samples?
#
#   This is NOT an attempt to rescue the CA2-SLM signal, and NOT an attempt to
#   eliminate it. Every rule below is fixed before any protein is inspected.
#
# WHAT IS AND IS NOT RECOMPUTED
#   The canonical differential abundance is READ, never refitted. log2fc, pval
#   and padj are carried through verbatim from the canonical contrast file. No
#   new FDR is computed anywhere in this script. Every quantity produced here
#   is a DESCRIPTIVE robustness diagnostic on the same animal-level values the
#   canonical model saw, and is labelled as such.
#
# THE DATA THE CANONICAL MODEL ACTUALLY SAW
#   Each animal-level value is the unweighted mean of that animal's Left and
#   Right POST-imputation log2 values (aggregation_audit.csv records
#   equal_weight_mean_LR_on_existing_imputed_log2_values for all 90 units).
#   So the contrast is n = 3 vs 3 on animal-level columns, and an animal whose
#   hemisphere was never detected still carries a drawn value.
#
#   Imputation status is therefore derived from the PRE-imputation raw matrix
#   (data/raw/pg_matrix/quicksearch.pg_matrix.tsv, where missing is a literal
#   NA), never from the post-imputation matrix, which has no NAs left at all.
#
# PRESPECIFIED RULES (fixed before looking at any protein)
#
#   OBSERVED-ONLY ESTIMABILITY. An animal contributes to an observed-only
#   estimate only if BOTH of its hemispheres were detected, because its value
#   is the mean of the two and a half-imputed mean is partly synthetic. The
#   observed-only effect is estimable only with at least MIN_OBS_ANIMALS = 2
#   fully observed animals in EACH of SUS and RES. A 1-versus-3 comparison is
#   not reported as an effect.
#
#   IMPUTATION DEPENDENCE. Relative to the canonical effect magnitude:
#     <= 0.25 change  -> minimal
#     <= 0.50 change  -> moderate
#     >  0.50 change or a sign flip -> strong
#
#   MAGNITUDE COLLAPSE. Losing more than HALF the canonical effect magnitude
#   in a sensitivity analysis is a collapse.
#
# USAGE
#   Rscript analysis/03_spatial_validation/16_ca2_slm_robustness_audit.R
#   Rscript analysis/03_spatial_validation/16_ca2_slm_robustness_audit.R --dry-run

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/spatial/spatial_atlas_utils.R")
source("R/spatial/ca2_slm_robustness_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/03_spatial_validation/16_ca2_slm_robustness_audit.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "neuron_neuropil")

DS <- "neuron_neuropil"
UNIT <- "CA2_slm"
DA_DIR <- path_processed("02_id_mapping", "mapped", DS, "forward", "per_file")
DA_FILE <- file.path(DA_DIR, "CA2slmsus_CA2slmres.csv")
RAW <- repo_path("data", "raw", "pg_matrix", "quicksearch.pg_matrix.tsv")
META <- repo_path("data", "metadata", "TPE9_sample_metadata_males.xlsx")
GCT <- path_processed("01_preprocessing", "protigy_input_animal_level", DS,
                      "neuron_neuropil_animal_level.gct")
QC <- path_results("tables", "03_qc_exploration", "00_dataset_qc_report", DS,
                   "dataset_qc_outlier_flags.csv")
IMPUTED <- path_processed("01_preprocessing", "impute",
  "20260601_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct.xlsx")
OUT <- function(...) {
  d <- path_results("tables", "11_spatial_systems", "ca2_slm_robustness")
  dir_create(d); file.path(d, ...)
}

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] CA2-SLM DAP robustness audit.\n")
  dry_run_inputs(SCRIPT_ID, list(canonical_da = DA_FILE, raw_pg_matrix = RAW,
                                 sample_metadata = META, animal_level_gct = GCT,
                                 imputed_sample_matrix = IMPUTED,
                                 qc_outlier_flags = QC))
  cat("[DRY-RUN] Canonical DA is READ, never refitted; no new FDR is computed.\n")
  cat("[DRY-RUN] Imputation status comes from the PRE-imputation raw matrix.\n")
  quit(save = "no", status = 0L)
}

for (p in c(DA_FILE, RAW, META, GCT, IMPUTED)) {
  if (!file.exists(p)) stop("missing_required_input: ", p, call. = FALSE)
}

# ===================================================== canonical DA (read only)

message("Reading the canonical CA2-SLM SUS-vs-RES contrast")
da <- as.data.frame(readr::read_csv(DA_FILE, show_col_types = FALSE,
                                    progress = FALSE, guess_max = Inf))
daps <- da[!is.na(da$padj) & da$padj < csr_fdr_threshold(), , drop = FALSE]
daps <- daps[order(daps$padj), , drop = FALSE]
message("  FDR-supported CA2-SLM DAPs: ", nrow(daps))

# ================================================ sample map and QC flags

meta <- as.data.frame(readxl::read_excel(META))
meta$AnimalID <- as.character(meta$AnimalID)
keep <- meta$celltype_layer == DS & !(meta$exclude %in% TRUE)
neuro <- meta[keep, , drop = FALSE]
ca2 <- neuro[neuro$region == "CA2" & neuro$layer == "slm", , drop = FALSE]
ca2$StressGroup <- csr_expgroup_to_stress(ca2$ExpGroup)
ca2$Hemi <- ifelse(tolower(ca2$ReplicateGroup) == "left", "L", "R")
if (nrow(ca2) != 18L) stop("expected 18 CA2_slm samples, found ", nrow(ca2), call. = FALSE)

qc_flag <- stats::setNames(rep("PASS", nrow(ca2)), ca2$sample_id)
if (file.exists(QC)) {
  q <- utils::read.csv(QC, stringsAsFactors = FALSE)
  m <- match(ca2$sample_id, q$Sample)
  qc_flag[!is.na(m)] <- q$outlier_flag[m[!is.na(m)]]
}
ca2$qc_flag <- unname(qc_flag[ca2$sample_id])
qc_fail_samples <- ca2$sample_id[ca2$qc_flag == "FAIL"]
qc_fail_animals <- sort(unique(ca2$AnimalID[ca2$qc_flag == "FAIL"]))
message("  QC in CA2_slm: ", paste(sprintf("%s_%s=%s", ca2$AnimalID, ca2$Hemi,
        ca2$qc_flag)[ca2$qc_flag != "PASS"], collapse = " "))

# ==================================== PRE-imputation observation mask

message("Deriving the pre-imputation observation mask from the raw matrix")
mask <- csr_preimputation_mask(RAW, neuro$sample_id, ca2$sample_id)
message("  analysed rows after the 70% filter: ", mask$n_analysed_rows,
        " (pre-imputation NAs over the 180 neuropil samples: ",
        mask$n_missing_total, ")")

# ================== sample-level values, and proof they rebuild the DA input

message("Reading the sample-level imputed matrix")
imputed <- csr_read_imputed_samples(IMPUTED, ca2$sample_id)

sample_meta <- data.frame(
  sample_id = ca2$sample_id,
  AnimalID = csr_bare_animal(ca2$AnimalID),
  hemisphere = ca2$Hemi,
  StressGroup = ca2$StressGroup,
  plate = ca2$plate,
  qc_flag = ca2$qc_flag,
  stringsAsFactors = FALSE)
animals <- sort(unique(sample_meta$AnimalID))
stress <- stats::setNames(
  sample_meta$StressGroup[match(animals, sample_meta$AnimalID)], animals)

ctx <- list(
  sample_meta = sample_meta, animals = animals, stress = stress,
  qc_fail_samples = qc_fail_samples,
  sample_of = function(a, h) sample_meta$sample_id[
    sample_meta$AnimalID == a & sample_meta$hemisphere == h][1])

# The animal-level value the differential model consumed is the unweighted
# mean of Left and Right. Prove the reconstruction is exact rather than
# assuming it: if this drifts, every sensitivity number below is meaningless.
gct <- csr_read_gct(GCT)
# the GCT keeps the raw A-prefixed animal labels; take them from the metadata
# rather than reconstructing the padding, which differs between A0003 and A111
raw_label <- ca2$AnimalID[match(animals, csr_bare_animal(ca2$AnimalID))]
unit_cols <- paste0(raw_label, "_", UNIT)
missing_cols <- setdiff(unit_cols, colnames(gct$mat))
if (length(missing_cols)) {
  stop("animal-level column(s) absent from the GCT: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}
shared <- intersect(rownames(imputed), rownames(gct$mat))
recon <- vapply(animals, function(a) {
  cols <- sample_meta$sample_id[sample_meta$AnimalID == a]
  rowMeans(imputed[shared, cols, drop = FALSE])
}, numeric(length(shared)))
gct_sub <- gct$mat[shared, unit_cols, drop = FALSE]
max_recon_err <- max(abs(recon - gct_sub), na.rm = TRUE)
message(sprintf("  L/R mean reproduces the animal-level matrix to %.3g over %d proteins",
                max_recon_err, length(shared)))
if (!is.finite(max_recon_err) || max_recon_err > 1e-8) {
  stop("the sample-level reconstruction does not reproduce the animal-level ",
       "matrix the differential model consumed (max abs error ",
       signif(max_recon_err, 3), "); sensitivity analysis would not be ",
       "describing the canonical values", call. = FALSE)
}

# =========================================== PARTS 7-13: per protein

message("Computing robustness diagnostics for ", nrow(daps), " DAPs")
rows <- list(); loo_rows <- list()
for (i in seq_len(nrow(daps))) {
  pid <- daps$original_identifier[i]
  if (!pid %in% rownames(imputed)) {
    stop("canonical DAP absent from the imputed matrix: ", pid, call. = FALSE)
  }
  mrow <- mask$unit_missing[pid, sample_meta$sample_id]
  names(mrow) <- sample_meta$sample_id
  r <- csr_protein_robustness(pid, daps[i, , drop = FALSE], imputed, mrow, ctx)
  rows[[length(rows) + 1L]] <- r$summary
  loo_rows[[length(loo_rows) + 1L]] <- r$loo
}
rob <- dplyr::bind_rows(rows)
loo_long <- dplyr::bind_rows(loo_rows)

# ------------------------------------------------- module / tier annotation
atlas_p <- path_results("tables", "11_spatial_systems", "atlas",
                        "protein_spatial_cell_affinity.csv")
if (file.exists(atlas_p)) {
  at <- as.data.frame(readr::read_csv(atlas_p, show_col_types = FALSE,
                                      progress = FALSE, guess_max = Inf))
  at <- at[at$dataset == DS, c("ProteinGroupID", "GeneSymbol", "ModuleID",
                               "module_display_label", "candidate_tier",
                               "abs_kME"), drop = FALSE]
  rob <- dplyr::left_join(rob, at, by = "ProteinGroupID")
}

rob <- rob |> dplyr::relocate(dplyr::any_of(c(
  "ProteinGroupID", "original_identifier", "GeneSymbol", "gene_symbol",
  "ModuleID", "module_display_label", "candidate_tier", "abs_kME",
  "canonical_log2FC_SUS_minus_RES", "canonical_p_value", "canonical_BH_FDR")))

# ================== normalisation-bias context (protein-independent shift)
#
# Per-sample median centring uses only OBSERVED values, so a heavily missing
# sample gets an inflated median and every one of its proteins is displaced
# downwards - including proteins with no imputed value at all. Measured here
# on proteins observed in all 18 samples, so imputation cannot contribute.

message("Measuring the missingness-driven normalisation displacement")
bias <- csr_normalization_bias_context(RAW, neuro$sample_id, sample_meta)
message(sprintf("  cor(missingness, observed median) = %+.3f", bias$cor_missing_vs_median))
message(sprintf("  cor(missingness, centred value of always-observed proteins) = %+.3f",
                bias$cor_missing_vs_centred))
message(sprintf("  implied SUS-RES displacement on EVERY protein = %+.4f log2",
                bias$sus_minus_res_displacement))

# How much of each canonical effect is already accounted for by that global,
# protein-independent shift? This does not refit anything; it compares the
# published effect with a displacement measured from always-observed proteins.
rob$global_normalization_displacement_SUS_minus_RES <- bias$sus_minus_res_displacement
rob$effect_beyond_global_displacement <-
  rob$canonical_log2FC_SUS_minus_RES - bias$sus_minus_res_displacement
rob$fraction_of_effect_explained_by_displacement <- ifelse(
  rob$canonical_log2FC_SUS_minus_RES == 0, NA_real_,
  bias$sus_minus_res_displacement / rob$canonical_log2FC_SUS_minus_RES)
rob$exceeds_global_displacement <-
  abs(rob$canonical_log2FC_SUS_minus_RES) > abs(bias$sus_minus_res_displacement)
rob$displacement_note <- paste0(
  "the displacement is a protein-independent consequence of per-sample median ",
  "centring under unequal missingness, measured on the ", bias$n_always_observed,
  " proteins observed in all 18 CA2_slm samples; it is NOT a refitted statistic")

# ======================================= PART 13: the fully observed subset

fully <- rob[rob$fully_observed %in% TRUE, , drop = FALSE]
fo_audit <- csr_fully_observed_audit(fully, loo_long, sample_meta)

# ------------------------------------------------------- sample context

sample_ctx <- sample_meta
sample_ctx$n_proteins_analysed <- mask$n_analysed_rows
sample_ctx$n_missing_preimputation <- as.integer(
  colSums(mask$unit_missing[, sample_meta$sample_id, drop = FALSE]))
sample_ctx$fraction_missing_preimputation <- as.numeric(
  colMeans(mask$unit_missing[, sample_meta$sample_id, drop = FALSE]))
sample_ctx$missingness_note <- paste0(
  "computed over the ", mask$n_analysed_rows, " rows that survive the 70% ",
  "filter, i.e. the rows actually analysed, not the full raw matrix")

# ----------------------------------------------------------------- write

write_csv_safe(rob, OUT("CA2_SLM_DAP_robustness.csv"))
write_csv_safe(loo_long, OUT("CA2_SLM_leave_one_animal_out_long.csv"))
write_csv_safe(fo_audit, OUT("CA2_SLM_fully_observed_DAP_audit.csv"))
write_csv_safe(sample_ctx, OUT("CA2_SLM_sample_context.csv"))
write_csv_safe(csr_thresholds(), OUT("CA2_SLM_prespecified_thresholds.csv"))
write_csv_safe(bias$per_sample, OUT("CA2_SLM_normalization_bias_context.csv"))

cat("\n===== CA2-SLM robustness audit =====\n")
cat(sprintf("canonical DAPs audited: %d (read, never refitted)\n", nrow(rob)))
cat(sprintf("fully observed in all 12 SUS+RES samples: %d\n", sum(rob$fully_observed)))
cat(sprintf("observed-only effect estimable (>=%d observed animals per group): %d\n",
            csr_min_observed_animals(), sum(rob$observed_only_estimable)))

cat("\n--- imputation dependence ---\n")
print(table(rob$imputation_dependence_class))

cat("\n--- final robustness classification ---\n")
print(table(rob$CA2_SLM_robustness_class))

cat("\n--- QC-sample dependence ---\n")
cat(sprintf("  A755-sensitive: %d   A764-sensitive: %d   single-animal-sensitive: %d\n",
            sum(rob$A755_sensitive), sum(rob$A764_sensitive),
            sum(rob$single_animal_sensitive)))
cat(sprintf("  LOO sign stable: %d of %d\n", sum(rob$loo_sign_stable), nrow(rob)))

cat("\n--- missingness-driven normalisation displacement ---\n")
cat(sprintf("  cor(missingness, observed median)                 = %+.3f\n",
            bias$cor_missing_vs_median))
cat(sprintf("  cor(missingness, centred always-observed value)   = %+.3f\n",
            bias$cor_missing_vs_centred))
cat(sprintf("  SUS-RES displacement applied to EVERY protein     = %+.4f log2\n",
            bias$sus_minus_res_displacement))
cat(sprintf("  DAPs whose effect exceeds that displacement       = %d of %d\n",
            sum(rob$exceeds_global_displacement), nrow(rob)))
cat(sprintf("  median share of effect explained by displacement  = %.0f%%\n",
            100 * stats::median(rob$fraction_of_effect_explained_by_displacement[
              rob$canonical_log2FC_SUS_minus_RES < 0], na.rm = TRUE)))

cat("\n--- the 8 fully observed DAPs ---\n")
if (nrow(fo_audit)) {
  for (i in seq_len(nrow(fo_audit))) {
    cat(sprintf("  %-10s log2FC=%6.3f  LOO[%6.3f,%6.3f]  signStable=%-5s  QClev=%.2f  %s\n",
                substr(fo_audit$GeneSymbol[i], 1, 10),
                fo_audit$canonical_log2FC_SUS_minus_RES[i],
                fo_audit$loo_min_effect[i], fo_audit$loo_max_effect[i],
                fo_audit$loo_sign_stable[i], fo_audit$max_qc_animal_leverage[i],
                fo_audit$CA2_SLM_robustness_class[i]))
  }
}
cat("\nOutputs:", relative_to(dirname(OUT("x"))), "\n")
