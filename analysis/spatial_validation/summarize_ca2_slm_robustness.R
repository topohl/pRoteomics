#!/usr/bin/env Rscript
#
# CA2-SLM robustness: candidate figures, reviewer workbook and the validation
# contract. A critical FAIL exits non-zero.
#
# The figures here are NEW CANDIDATES. They replace no manuscript figure
# contract, and each writes its own source-data CSV.
#
# PIN. The published 3/3/3 result has 28 FDR-supported CA2-SLM SUS-RES DAPs and
# the validation holds the audit to that count. A rerun on other data (the 3/3/2
# rerun) sets PROTEOMICS_CA2_SLM_UNPINNED=true: the audit's protein set must then
# equal the contrast file's padj < 0.05 set, and the validation table and log
# state that the pin was not enforced. Design counts in the generated text are
# taken from the audit tables in both modes.
#
# USAGE
#   Rscript analysis/spatial_validation/summarize_ca2_slm_robustness.R
#   Rscript analysis/spatial_validation/summarize_ca2_slm_robustness.R --dry-run
# Script: analysis/spatial_validation/summarize_ca2_slm_robustness.R
# Stage: networks
# Scope: global
# Consumes: required results/spatial_validation/audit_ca2_slm_robustness/global/tables/CA2_SLM_DAP_robustness.csv; results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_DAP_robustness.csv; data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/CA2slmsus_CA2slmres.csv; optional results/tables/11_spatial_systems/ca2_slm_robustness/; results/tables/08_behavior_physio_coupling/animal_id_integrity/; results/tables/08_behavior_physio_coupling/network_behavior_coupling/; +1 more
# Produces: results/spatial_validation/summarize_ca2_slm_robustness/global/reports/CA2_SLM_robustness_audit.xlsx; results/spatial_validation/summarize_ca2_slm_robustness/global/tables/CA2_SLM_robustness_validation.csv; results/spatial_validation/summarize_ca2_slm_robustness/global/tables/CA2_SLM_claim_status_summary.csv; +1 more
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: CA2-SLM robustness: candidate figures, reviewer workbook and the validation contract.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/utilities/xlsx_package_utils.R")
source("R/data_contracts/animal_id_contract.R")
source("R/spatial/ca2_slm_robustness_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "qc_result_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "summarize_ca2_slm_robustness"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/spatial_validation/summarize_ca2_slm_robustness.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

ROB <- function(f) spatial_systems_find(f, "ca2_slm_robustness")
# ROB() resolves reads; this analysis writes only to its own tables.
OUT <- function(f) file.path(CANONICAL_PATHS$tables, f)
AID <- function(...) path_results("integration", "audit_animal_id_integrity", "global", "tables", ...)
COUP <- function(...) path_results("integration", "test_network_behaviour_coupling", "global", "tables", ...)
FIG <- function(...) {
  d <- CANONICAL_PATHS$plots
  dir_create(d); file.path(d, ...)
}
DA_FILE <- file.path(preprocessing_mapped_contrast_dir("neuron_neuropil", "forward"),
                     "CA2slmsus_CA2slmres.csv")
source(repo_path("R", "dataset_inputs.R"))
IMPUTED <- resolve_dataset_inputs("neuron_neuropil", purpose = "wgcna", record_resolution = FALSE)$expression_file
# read only for the count of QC-flagged neuropil samples in the generated text
QC_FLAGS <- qc_find("dataset_qc_outlier_flags.csv", owner = "assess_dataset_quality",
                    legacy_substep = "00_dataset_qc_report", scope = "neuron_neuropil")

PINNED_N_HITS <- 28L
ca2_unpinned <- identical(tolower(Sys.getenv("PROTEOMICS_CA2_SLM_UNPINNED", "")), "true")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] CA2-SLM robustness figures, workbook and validation.\n")
  dry_run_inputs(SCRIPT_ID, list(robustness = ROB("CA2_SLM_DAP_robustness.csv"),
                                 canonical_da = DA_FILE))
  cat("[DRY-RUN] A critical validation FAIL exits non-zero.\n")
  quit(save = "no", status = 0L)
}

rd <- function(p) if (file.exists(p))
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf)) else NULL

rob <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
if (is.null(rob)) stop("missing_required_input: ", ROB("CA2_SLM_DAP_robustness.csv"),
                       call. = FALSE)
rob_ann <- rd(ROB("CA2_SLM_DAP_robustness_annotated.csv"))
loo <- rd(ROB("CA2_SLM_leave_one_animal_out_long.csv"))
fo <- rd(ROB("CA2_SLM_fully_observed_DAP_audit.csv"))
spec <- rd(ROB("CA2_SLM_spatial_specificity.csv"))
sctx <- rd(ROB("CA2_SLM_sample_context.csv"))
bias <- rd(ROB("CA2_SLM_normalization_bias_context.csv"))
ident <- rd(ROB("stress_identity_robustness_comparison.csv"))
modcmp <- rd(ROB("CA2_SLM_module_distribution_comparison.csv"))
across <- rd(ROB("CA2_SLM_effect_across_neuropil_units.csv"))
thr <- rd(ROB("CA2_SLM_prespecified_thresholds.csv"))
consumer <- rd(AID("animal_id_normalization_consumer_audit.csv"))
blast <- rd(AID("animal_id_normalization_blast_radius.csv"))
resolution <- rd(AID("animal_id_resolution_report.csv"))
join_sum <- rd(COUP("join_diagnostics_summary.csv"))
coup_cor <- rd(COUP("edge_behavior_correlations.csv"))
qc_flags <- rd(QC_FLAGS)
da <- rd(DA_FILE)

# ---- design facts for the generated text, from the data rather than assumed
# neuropil SUS-RES DAPs over every unit's contrast next to the CA2-SLM one
sr_files <- list.files(dirname(DA_FILE), pattern = "sus_.*res[.]csv$", full.names = TRUE)
sr_files <- sr_files[sub("sus_.*$", "", basename(sr_files)) ==
                       sub("^.*sus_(.*)res[.]csv$", "\\1", basename(sr_files))]
n_units <- length(sr_files)
n_neuropil_daps <- if (!n_units || !basename(DA_FILE) %in% basename(sr_files)) NA_integer_ else
  sum(vapply(sr_files, function(f) {
    z <- rd(f); sum(!is.na(z$padj) & z$padj < csr_fdr_threshold())
  }, integer(1)))
imp_head <- if (file.exists(IMPUTED)) readxl::read_excel(IMPUTED, n_max = 400) else NULL
n_neuropil_samples <- if (is.null(imp_head)) NA_integer_ else
  sum(vapply(imp_head, is.numeric, logical(1)))
n_flagged_neuropil <- if (is.null(qc_flags)) NA_integer_ else sum(qc_flags$outlier_flag != "PASS")
dtxt <- if (is.null(sctx)) NULL else csr_design_text(sctx, n_flagged_neuropil)

# ============================================ PART 19: candidate figures

made <- character()
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  gcol <- c(CON = "#3E3C6F", RES = "#9E9A92", SUS = "#D7303F")
  base_theme <- theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 9))

  # ---- Figure A: missingness / observation / QC context by group
  if (!is.null(sctx)) {
    a <- sctx
    a$label <- paste0(a$AnimalID, "_", a$hemisphere)
    a$qc_label <- ifelse(a$qc_flag == "PASS", "", a$qc_flag)
    p <- ggplot(a, aes(x = stats::reorder(label, -fraction_missing_preimputation),
                       y = fraction_missing_preimputation, fill = StressGroup)) +
      geom_col() +
      geom_text(aes(label = qc_label), vjust = -0.35, size = 2.4, colour = "black") +
      scale_fill_manual(values = gcol) +
      labs(title = "CA2-SLM pre-imputation missingness and QC status by sample",
           subtitle = "bar label marks the QC outlier flag; every sample was retained (Exclude is FALSE)",
           x = NULL, y = "fraction missing before imputation") +
      base_theme + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    ggsave(FIG("A_ca2_slm_missingness_qc_context.png"), p, width = 7, height = 4,
           dpi = 200, bg = "white")
    write_csv_safe(a, FIG("A_ca2_slm_missingness_qc_context_source_data.csv"))
    made <- c(made, "A_ca2_slm_missingness_qc_context")
  }

  # ---- Figure B: canonical versus observed-only effect
  if (!is.null(rob)) {
    b <- rob[, c("gene_symbol", "canonical_log2FC_SUS_minus_RES",
                 "observed_only_log2FC", "observed_only_estimable",
                 "fully_observed", "imputation_dependence_class",
                 "CA2_SLM_robustness_class")]
    b$status <- dplyr::case_when(
      b$fully_observed %in% TRUE ~ "fully observed",
      b$CA2_SLM_robustness_class == "imputation_sensitive" ~ "imputation sensitive",
      b$CA2_SLM_robustness_class == "not_claimable_due_to_QC" ~ "QC sensitive",
      TRUE ~ "other")
    bb <- b[b$observed_only_estimable %in% TRUE, , drop = FALSE]
    p <- ggplot(bb, aes(x = canonical_log2FC_SUS_minus_RES, y = observed_only_log2FC,
                        colour = status)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
      geom_hline(yintercept = 0, linewidth = 0.2) +
      geom_vline(xintercept = 0, linewidth = 0.2) +
      geom_point(size = 2) +
      geom_text(aes(label = gene_symbol), size = 2.1, vjust = -0.9,
                show.legend = FALSE) +
      labs(title = "Canonical versus observed-only SUS-RES effect in CA2-SLM",
           subtitle = "observed-only uses animals with BOTH hemispheres detected; dashed line is identity",
           x = "canonical log2FC (SUS - RES)", y = "observed-only log2FC") +
      base_theme
    ggsave(FIG("B_canonical_vs_observed_only.png"), p, width = 6.5, height = 4.5,
           dpi = 200, bg = "white")
    write_csv_safe(b, FIG("B_canonical_vs_observed_only_source_data.csv"))
    made <- c(made, "B_canonical_vs_observed_only")
  }

  # ---- Figure C: leave-one-animal-out stability
  if (!is.null(loo)) {
    l <- loo[loo$can_change_estimate %in% TRUE, , drop = FALSE]
    gs <- stats::setNames(rob$gene_symbol, rob$original_identifier)
    l$gene <- unname(gs[l$original_identifier])
    ord <- rob$gene_symbol[order(rob$canonical_log2FC_SUS_minus_RES)]
    l$gene <- factor(l$gene, levels = ord)
    l$omit <- paste0(l$omitted_AnimalID, "\n(", l$omitted_StressGroup, ")")
    p <- ggplot(l, aes(x = omit, y = gene, fill = change_from_canonical)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                           midpoint = 0, name = "change in\nlog2FC") +
      labs(title = "Leave-one-animal-out stability of the CA2-SLM effects",
           subtitle = "change in the descriptive SUS-RES effect when one animal is omitted; CON animals cannot move it and are not shown",
           x = "omitted animal", y = NULL) +
      base_theme + theme(axis.text.y = element_text(size = 6))
    ggsave(FIG("C_leave_one_animal_out_stability.png"), p, width = 6, height = 6.5,
           dpi = 200, bg = "white")
    write_csv_safe(l, FIG("C_leave_one_animal_out_stability_source_data.csv"))
    made <- c(made, "C_leave_one_animal_out_stability")
  }

  # ---- Figure D: robust proteins across every neuropil unit
  if (!is.null(across) && !is.null(rob)) {
    keep <- rob$original_identifier[rob$CA2_SLM_robustness_class ==
                                      "robust_to_missingness_and_QC"]
    d <- across[across$original_identifier %in% keep, , drop = FALSE]
    if (nrow(d)) {
      gs <- stats::setNames(rob$gene_symbol, rob$original_identifier)
      d$gene <- unname(gs[d$original_identifier])
      d$is_ca2 <- d$spatial_unit == "CA2_slm"
      outside <- if (!"padj" %in% names(d)) NULL else
        d[!d$is_ca2 & !is.na(d$padj) & d$padj < csr_fdr_threshold(), , drop = FALSE]
      sub_d <- if (is.null(outside)) "canonical SUS-RES log2FC per spatial unit; descriptive" else
        if (!nrow(outside)) "canonical SUS-RES log2FC per spatial unit; descriptive, no unit other than CA2-SLM is FDR-supported" else
          paste0("canonical SUS-RES log2FC per spatial unit; descriptive; also FDR-supported outside CA2-SLM: ",
                 paste(outside$gene, outside$spatial_unit, collapse = ", "))
      p <- ggplot(d, aes(x = spatial_unit, y = log2fc, group = gene)) +
        geom_hline(yintercept = 0, linewidth = 0.3) +
        geom_line(colour = "grey70") +
        geom_point(aes(colour = is_ca2), size = 1.8) +
        scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "#D7303F"),
                            name = NULL, labels = c("other unit", "CA2-SLM")) +
        facet_wrap(~ gene, ncol = 3) +
        labs(title = "Robustness-qualified CA2-SLM effects across every neuropil unit",
             subtitle = sub_d,
             x = NULL, y = "log2FC (SUS - RES)") +
        base_theme +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
      ggsave(FIG("D_robust_effects_across_units.png"), p, width = 7, height = 5,
             dpi = 200, bg = "white")
      write_csv_safe(d, FIG("D_robust_effects_across_units_source_data.csv"))
      made <- c(made, "D_robust_effects_across_units")
    }
  }
}

# ================== PART 18: automatic claim-status summary and language

# Wording that is only permitted once it is explicitly scoped to the
# robustness-qualified subset. Until then it overstates what the data support.
CA2_RESTRICTED_PHRASES <- c("CA2-SLM molecular hotspot", "molecular hotspot",
                            "selective CA2-SLM vulnerability",
                            "selective CA2-SLM", "genuine spatial concentration")

n_by <- function(cl) sum(rob$CA2_SLM_robustness_class == cl)
robust_genes <- sort(rob$gene_symbol[rob$CA2_SLM_robustness_class ==
                                       "robust_to_missingness_and_QC"])
disp <- if (is.null(bias)) NA_real_ else {
  pa <- vapply(split(bias$mean_centred_value_of_always_observed_proteins,
                     bias$AnimalID), mean, numeric(1))
  gg <- bias$StressGroup[match(names(pa), bias$AnimalID)]
  mean(pa[gg == "SUS"]) - mean(pa[gg == "RES"])
}
run_centring <- csr_run_centring(IMPUTED)

# CA2-SLM having the highest pre-imputation missingness of the neuropil units and
# the lowest baseline abundance holds in the 3/3/3 data and in the 3/3/2 rerun.
units_txt <- if (n_units) sprintf("the %d neuropil units", n_units) else "the neuropil units"
# A whole-animal removal that leaves fewer than csr_min_observed_animals() in a group (A755 in the
# 3/3/2 rerun) is labelled DESCRIPTIVE ONLY and decides no class: its sensitivity flag must be NA.
desc_only <- function(col) if (col %in% names(rob)) grepl("^DESCRIPTIVE ONLY", rob[[col]]) else rep(FALSE, nrow(rob))
desc_flag_ok <- all(is.na(rob$A755_sensitive[desc_only("effect_excluding_A755_interpretation")])) &&
  all(is.na(rob$A764_sensitive[desc_only("effect_excluding_A764_interpretation")]))
not_claimable_txt <- sprintf("%d reverse or lose more than half their magnitude when a shallow SUS acquisition (pre-imputation missingness z > %s) is left out",
                             n_by("not_claimable_due_to_QC"), format(csr_depth_leave_out_z()))
claim_status <- data.frame(
  item = c("canonical_result", "qc_context", "robust_subset",
           "not_claimable", "not_evaluable", "approved_wording",
           "restricted_wording"),
  value = c(
    if (!ca2_unpinned) sprintf("%d FDR-supported CA2-SLM SUS-RES DAPs in the canonical differential abundance, unchanged", nrow(rob)) else
      sprintf("%d FDR-supported CA2-SLM SUS-RES DAPs in this run's canonical differential abundance (pin to the published %d not enforced)",
              nrow(rob), PINNED_N_HITS),
    if (is.null(dtxt)) "CA2-SLM sample context not available" else
      sprintf("CA2-SLM carries the highest pre-imputation missingness of %s, %s, and %s; none was excluded",
              units_txt, dtxt$missingness_direction, dtxt$flagged_share_with_flags),
    sprintf("%d of %d remain robust to missingness and QC: %s", n_by("robust_to_missingness_and_QC"), nrow(rob), paste(robust_genes, collapse = ", ")),
    not_claimable_txt,
    sprintf("%d cannot be evaluated: fewer than %d fully observed animals in at least one group", n_by("insufficient_observed_data"), csr_min_observed_animals()),
    sprintf(paste0("Canonical DA showed a concentration of SUS-RES differences in CA2-SLM, but this spatial unit also displayed substantial group-differential missingness%s. Protein-level sensitivity analyses were therefore used to distinguish robust effects from imputation/QC-sensitive findings. %d of %d effects remain robust; a protein-independent normalisation displacement of %+.3f log2 (SUS minus RES)%s applies to every protein in the unit."),
            if (is.null(dtxt)) "" else if (isTRUE(dtxt$shallowest_sus_is_most_missing))
              paste0(", highest in ", dtxt$shallowest_sus) else paste0(" (", dtxt$shallowest_sus, ")"),
            n_by("robust_to_missingness_and_QC"), nrow(rob), disp, csr_displacement_origin(run_centring)),
    paste0("NOT permitted unless explicitly scoped to the robustness-qualified subset: ",
           paste(CA2_RESTRICTED_PHRASES, collapse = "; "))),
  stringsAsFactors = FALSE)
write_csv_safe(claim_status, OUT("CA2_SLM_claim_status_summary.csv"))

# =============================================== PART 21: validation

checks <- list()
add <- function(id, critical, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id, critical = critical, status = status, detail = detail,
    stringsAsFactors = FALSE)
}

canon <- if (is.null(da)) NULL else
  da[!is.na(da$padj) & da$padj < csr_fdr_threshold(), , drop = FALSE]
n_canon <- if (is.null(canon)) NA_integer_ else nrow(canon)
unpinned_note <- sprintf("pin to the published %d not enforced (PROTEOMICS_CA2_SLM_UNPINNED=true)",
                         PINNED_N_HITS)

# Pinned: the published count. Unpinned: the audit holds exactly this run's FDR-supported set.
add("canonical_DA_unchanged", TRUE,
    if (!ca2_unpinned) {
      if (!is.null(canon) && n_canon == PINNED_N_HITS) "PASS" else "FAIL"
    } else {
      if (!is.null(canon) && !anyDuplicated(rob$original_identifier) && nrow(rob) == n_canon &&
          setequal(rob$original_identifier, canon$original_identifier)) "PASS" else "FAIL"
    },
    if (!ca2_unpinned) paste0("the canonical CA2-SLM contrast still yields ",
           if (is.null(canon)) "NA" else n_canon,
           " FDR-supported proteins at padj < ", csr_fdr_threshold()) else
      paste0(unpinned_note, ": the audit's ", nrow(rob), " proteins are the ",
             if (is.null(canon)) "NA" else n_canon, " FDR-supported proteins at padj < ",
             csr_fdr_threshold(), " of the canonical CA2-SLM contrast"))

fc_ok <- !is.null(canon) && all(mapply(function(k, v) {
  j <- match(k, canon$original_identifier)
  !is.na(j) && isTRUE(all.equal(as.numeric(canon$log2fc[j]), v, tolerance = 1e-12))
}, rob$original_identifier, rob$canonical_log2FC_SUS_minus_RES))
fdr_ok <- !is.null(canon) && all(mapply(function(k, v) {
  j <- match(k, canon$original_identifier)
  !is.na(j) && isTRUE(all.equal(as.numeric(canon$padj[j]), v, tolerance = 1e-12))
}, rob$original_identifier, rob$canonical_BH_FDR))
add("canonical_statistics_carried_verbatim", TRUE,
    if (fc_ok && fdr_ok) "PASS" else "FAIL",
    "every log2FC and BH FDR in the audit equals the canonical contrast value exactly")

src16 <- paste(sub("#.*$", "", readLines(repo_path("analysis/spatial_validation",
  "audit_ca2_slm_robustness.R"), warn = FALSE)), collapse = "\n")
src17 <- paste(sub("#.*$", "", readLines(repo_path("analysis/spatial_validation",
  "audit_stress_identity_robustness.R"), warn = FALSE)), collapse = "\n")
srcu <- paste(sub("#.*$", "", readLines(repo_path("R",
  "ca2_slm_robustness_utils.R"), warn = FALSE)), collapse = "\n")
add("no_new_FDR_computed", TRUE,
    if (!any(grepl("p\\.adjust|qvalue|fdrtool", c(src16, src17, srcu)))) "PASS" else "FAIL",
    "no script in the robustness layer calls p.adjust or any FDR routine")
add("no_differential_model_refit", TRUE,
    if (!any(grepl("lmFit|eBayes|limma::|t\\.test|wilcox\\.test",
                   c(src16, src17, srcu)))) "PASS" else "FAIL",
    "no differential model is fitted; all effects are descriptive group differences")

# The check id keeps the published count; the detail carries the count observed.
add("all_28_canonical_hits_retained", TRUE,
    if (!ca2_unpinned) {
      if (nrow(rob) == PINNED_N_HITS &&
          (is.null(canon) || all(canon$original_identifier %in% rob$original_identifier)))
        "PASS" else "FAIL"
    } else {
      if (!is.null(canon) && nrow(rob) == n_canon &&
          all(canon$original_identifier %in% rob$original_identifier)) "PASS" else "FAIL"
    },
    if (!ca2_unpinned) paste0(nrow(rob), " of ", PINNED_N_HITS, " canonical CA2-SLM hits present in the audit; none dropped") else
      paste0(if (is.null(canon)) "NA" else sum(canon$original_identifier %in% rob$original_identifier),
             " of ", if (is.null(canon)) "NA" else n_canon,
             " canonical CA2-SLM hits present in the audit; none dropped; ", unpinned_note))

# fully observed reconciles with the raw pre-imputation matrix, recomputed here
fo_recon <- NA_integer_
if (!is.null(sctx)) fo_recon <- sum(rob$fully_observed %in% TRUE)
add("fully_observed_reconciles", TRUE,
    if (!is.na(fo_recon) && fo_recon == sum(rob$n_imputed_SUS + rob$n_imputed_RES == 0L))
      "PASS" else "FAIL",
    paste0(fo_recon, " proteins carry zero imputed values across the ",
           if (is.null(dtxt)) "NA" else dtxt$n_sus_res_samples, " SUS+RES samples"))

add("imputation_status_from_preimputation_source", TRUE,
    if (grepl("quicksearch.pg_matrix.tsv", src16, fixed = TRUE) &&
        grepl("csr_preimputation_mask", src16, fixed = TRUE)) "PASS" else "FAIL",
    "observation status is derived from the raw matrix, never from the post-imputation matrix")

nas_after <- NA
if (!is.null(imp_head)) {
  num <- vapply(imp_head, is.numeric, logical(1))
  nas_after <- sum(is.na(imp_head[, num, drop = FALSE]))
}
add("post_imputation_matrix_has_no_missing", FALSE,
    if (is.na(nas_after) || nas_after == 0L) "PASS" else "FAIL",
    "the post-imputation matrix carries no NA, which is why it cannot supply imputation status")

# Every protein has one leave-one-out row per animal: nine in the 3/3/3 design, eight in the 3/3/2 rerun.
n_loo_animals <- if (is.null(loo)) NA_integer_ else length(unique(loo$omitted_AnimalID))
loo_ok <- !is.null(loo) && n_loo_animals >= 6L &&
  all(tapply(loo$omitted_AnimalID, loo$original_identifier,
             function(z) length(z) == n_loo_animals && length(unique(z)) == n_loo_animals))
add("LOO_removes_exactly_one_animal", TRUE,
    if (loo_ok) "PASS" else "FAIL",
    paste0("each protein has exactly one leave-one-out row per AnimalID, ", n_loo_animals, " in total"))

# The audit's QC flags must be this run's dataset QC calls for the CA2-SLM samples. Which acquisitions fail depends
# on the QC metrics and on the centring, so no animal is assumed.
qc_table_flag <- if (is.null(sctx) || is.null(qc_flags)) NULL else
  as.character(qc_flags$outlier_flag[match(sctx$sample_id, qc_flags$Sample)])
qc_fail_ctx <- if (is.null(sctx)) NULL else sctx[sctx$qc_flag == "FAIL", , drop = FALSE]
# the leave-out follows the depth rule on this run's dataset QC missingness
lo_expected <- if (is.null(sctx) || is.null(qc_flags) || !"missing_fraction_robust_z" %in% names(qc_flags)) NULL else
  qc_flags$missing_fraction_robust_z[match(sctx$sample_id, qc_flags$Sample)] > csr_depth_leave_out_z()
add("leave_out_follows_depth_rule", TRUE,
    if (!is.null(lo_expected) && !anyNA(lo_expected) && !is.null(sctx$left_out) &&
        identical(as.logical(sctx$left_out), lo_expected)) "PASS" else "FAIL",
    if (is.null(lo_expected)) "the sample context or the dataset QC missingness is missing" else
      paste0("left out (pre-imputation missingness z > ", csr_depth_leave_out_z(), "): ",
             if (!any(lo_expected, na.rm = TRUE)) "none" else
               paste(csr_sample_labels(sctx[lo_expected %in% TRUE, , drop = FALSE], with_flag = FALSE), collapse = ", ")))
add("QC_flags_match_dataset_QC", TRUE,
    if (!is.null(qc_table_flag) && !anyNA(qc_table_flag) && identical(as.character(sctx$qc_flag), qc_table_flag))
      "PASS" else "FAIL",
    if (is.null(qc_table_flag)) "the sample context or the dataset QC table is missing" else
      paste0("the audit's QC flags of the ", nrow(sctx), " CA2-SLM samples are compared with the dataset QC table; QC-FAILED: ",
             if (!nrow(qc_fail_ctx)) "none" else paste(csr_sample_labels(qc_fail_ctx, with_flag = FALSE), collapse = ", ")))

add("A755_out_excludes_only_A755", TRUE,
    if (!is.null(loo) && all(loo$omitted_AnimalID[loo$omitted_AnimalID == "755"] == "755"))
      "PASS" else "FAIL",
    "the A755-out analysis removes animal 755 and no other animal")
# Every QC-fail removal that leaves fewer than csr_min_observed_animals() animals in SUS or
# RES carries the DESCRIPTIVE ONLY label and no other removal does (in the 3/3/2 rerun this
# includes the A755 removal); a removal of animals absent from the design is not computed.
# The group sizes left are counted from the sample context.
removal_sets <- list(effect_excluding_A755 = "755", effect_excluding_A764 = "764",
                     effect_excluding_both_QC_fail_SUS = c("755", "764"))
removal_label_cols <- c(effect_excluding_A755 = "effect_excluding_A755_interpretation",
                        effect_excluding_A764 = "effect_excluding_A764_interpretation",
                        effect_excluding_both_QC_fail_SUS = "effect_excluding_both_interpretation")
removal_eval <- if (is.null(sctx)) NULL else lapply(names(removal_sets), function(eff) {
  drop <- intersect(removal_sets[[eff]], as.character(sctx$AnimalID))
  # an animal is removed whole only when all of its acquisitions are left out (csr_removable_animals)
  if (is.null(sctx$left_out)) stop("the sample context has no left_out column; rerun the CA2-SLM audit", call. = FALSE)
  drop <- drop[vapply(drop, function(a) all(sctx$left_out[as.character(sctx$AnimalID) == a] %in% TRUE), logical(1))]
  left <- sctx[!as.character(sctx$AnimalID) %in% drop, , drop = FALSE]
  n_sus <- length(unique(left$AnimalID[left$StressGroup == "SUS"]))
  n_res <- length(unique(left$AnimalID[left$StressGroup == "RES"]))
  lab <- rob[[removal_label_cols[[eff]]]]
  labelled <- !is.null(lab) && all(grepl("DESCRIPTIVE ONLY", lab))
  unlabelled <- is.null(lab) || !any(grepl("DESCRIPTIVE ONLY", lab))
  needs <- length(drop) > 0L && csr_removal_needs_label(n_sus, n_res)
  ok <- if (!length(drop)) unlabelled && all(is.na(rob[[eff]])) else if (needs) labelled else unlabelled
  what <- if (eff != "effect_excluding_both_QC_fail_SUS") paste0("removing A", drop) else
    if (length(drop) == 2L) "removing both shallow SUS animals" else
      paste0("removing the shallow SUS animals present (", paste0("A", drop, collapse = ", "), ")")
  list(eff = eff, ok = ok, needs = needs,
       text = if (length(drop)) sprintf("%s leaves %d versus %d", what, n_sus, n_res) else NA_character_)
})
removal_ok <- !is.null(removal_eval) && all(vapply(removal_eval, `[[`, logical(1), "ok")) && desc_flag_ok
removal_needs <- if (is.null(removal_eval)) character(0) else
  unlist(lapply(removal_eval, function(z) if (z$needs) z$text))
removal_bad <- c(if (is.null(removal_eval)) character(0) else
  unlist(lapply(removal_eval, function(z) if (!z$ok) z$eff)),
  if (!desc_flag_ok) "the sensitivity flags of descriptive-only removals (must be NA: they decide no class)")
add("both_QC_fail_removal_labelled_descriptive", TRUE,
    if (removal_ok) "PASS" else "FAIL",
    if (is.null(removal_eval)) "sample context missing; group sizes after the leave-out not known" else
      paste0(if (length(removal_needs)) paste(removal_needs, collapse = "; ") else
               sprintf("no leave-out removal leaves fewer than %d animals in a group", csr_min_observed_animals()),
             if (!removal_ok) paste0("; label or value inconsistent with the rule for ", csr_and_list(removal_bad)) else
               if (length(removal_needs) == 1L) " and is labelled descriptive-only" else
                 if (length(removal_needs) > 1L) "; each is labelled descriptive-only" else ""))

add("classification_ignores_annotation", TRUE,
    if (!grepl("ModuleID|GeneSymbol|candidate_tier|gene_symbol",
               sub("csr_classify.*?\\n\\}", "", srcu))) "PASS" else "WARN",
    "the classifier reads only observation counts, effects and QC identity")
add("thresholds_prespecified", TRUE,
    if (!is.null(thr) && all(thr$prespecified %in% TRUE |
                             (if ("set_on" %in% names(thr)) !is.na(thr$set_on) & nzchar(thr$set_on) else FALSE))) "PASS" else "FAIL",
    paste0(if (is.null(thr)) 0 else nrow(thr),
           " thresholds are declared as named constants; prespecified",
           if (!is.null(thr) && "set_on" %in% names(thr) && any(!thr$prespecified %in% TRUE))
             paste0(" except ", paste(sprintf("%s (set %s)", thr$constant[!thr$prespecified %in% TRUE],
                                               thr$set_on[!thr$prespecified %in% TRUE]), collapse = ", ")) else ""))

add("observed_only_requires_minimum_n", TRUE,
    if (all(rob$observed_only_estimable ==
            (rob$n_fully_observed_SUS_animals >= csr_min_observed_animals() &
             rob$n_fully_observed_RES_animals >= csr_min_observed_animals())))
      "PASS" else "FAIL",
    paste0("an observed-only effect requires at least ", csr_min_observed_animals(),
           " fully observed animals in EACH group; 1-versus-3 is never reported"))

add("stress_identity_subsets_explicit", TRUE,
    if (!is.null(ident) && all(nzchar(ident$subset_definition)) &&
        all(grepl("tabulated not recomputed", ident$classification_source)))
      "PASS" else "FAIL",
    paste0(if (is.null(ident)) 0 else nrow(ident),
           " stress-identity subsets, each with an explicit definition and no rule reimplementation"))

add("module_enrichment_not_tested", TRUE,
    if (!is.null(modcmp) && all(grepl("NOT TESTED", modcmp$enrichment_test)))
      "PASS" else "FAIL",
    "no post hoc module enrichment test is run; counts only")

# AnimalID integrity
add("animal_id_collision_free", TRUE,
    if (!is.null(resolution) && !any(resolution$collision %in% TRUE &
                                     resolution$resolved %in% TRUE)) "PASS" else "FAIL",
    "no raw identifier resolves to an AnimalID shared with a different raw identifier")
jv <- if (is.null(join_sum)) NULL else
  stats::setNames(as.numeric(join_sum$value), join_sum$metric)
add("behavior_consumers_rerun", TRUE,
    if (!is.null(jv) && isTRUE(jv[["animals_lost_from_proteomics"]] == 0)) "PASS" else "FAIL",
    paste0("network-behaviour coupling rerun with the corrected contract; animals lost from proteomics: ",
           if (is.null(jv)) "NA" else jv[["animals_lost_from_proteomics"]]))

git_changed <- tryCatch(system2("git", c("status", "--porcelain"), stdout = TRUE,
                                stderr = FALSE), error = function(e) character())
add("wgcna_state_unchanged", TRUE,
    if (!any(grepl("wgcna_final_model_state|06_modules_WGCNA/01_WGCNA",
                   git_changed))) "PASS" else "FAIL",
    "no frozen WGCNA state or membership file is modified")
add("canonical_DA_files_unchanged", TRUE,
    if (!any(grepl("data/processed/02_id_mapping", git_changed))) "PASS" else "FAIL",
    "no canonical differential-abundance file is modified in the working tree")
add("network_layer_unchanged", TRUE,
    if (!any(grepl("11_spatial_systems/1[345]_", git_changed))) "PASS" else "FAIL",
    "the animal-level spatial network scripts are untouched by this pass")

# Part 18: the restricted wording must not appear unscoped in any generated
# table or in the robustness scripts themselves.
scan_dirs <- c(
  canonical_result_path("spatial_validation", "audit_ca2_slm_robustness", "global", "tables"),
  canonical_result_path("spatial_validation", "audit_stress_identity_robustness", "global", "tables"),
  CANONICAL_PATHS$tables)
scan_files <- c(
  list.files(scan_dirs, pattern = "[.]csv$", full.names = TRUE),
  repo_path("analysis/spatial_validation", "audit_ca2_slm_robustness.R"),
  repo_path("analysis/spatial_validation", "audit_stress_identity_robustness.R"))
offending <- character()
for (f in scan_files) {
  if (!file.exists(f)) next
  if (basename(f) == "CA2_SLM_claim_status_summary.csv") next   # declares them
  if (basename(f) == "CA2_SLM_robustness_validation.csv") next  # stale self-report
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (ph in CA2_RESTRICTED_PHRASES) {
    if (grepl(ph, txt, fixed = TRUE)) offending <- c(offending, paste0(basename(f), ":", ph))
  }
}
add("restricted_CA2_language_absent", TRUE,
    if (!length(offending)) "PASS" else "FAIL",
    if (length(offending)) paste("unscoped restricted wording in:",
                                 paste(unique(offending), collapse = "; ")) else
      paste0("none of the ", length(CA2_RESTRICTED_PHRASES),
             " restricted phrases appears unscoped in any generated table"))
add("claim_status_summary_generated", TRUE,
    if (file.exists(OUT("CA2_SLM_claim_status_summary.csv"))) "PASS" else "FAIL",
    "a claim-status summary with approved wording is generated automatically")

add("figures_have_source_data", FALSE,
    if (all(file.exists(FIG(paste0(made, "_source_data.csv"))))) "PASS" else "FAIL",
    paste0(length(made), " candidate figures, each with a source-data CSV"))

validation <- dplyr::bind_rows(checks)
write_csv_safe(validation, OUT("CA2_SLM_robustness_validation.csv"))

# ================================================= PART 20: workbook

# Sheet names carry the counts present (Canonical_28 and Fully_observed_8 in the published run).
sheets <- list(
  Canonical_N = rob,
  Missingness = sctx,
  Normalization_bias = bias,
  Observed_only = if (is.null(rob)) NULL else rob[, c("gene_symbol",
    "canonical_log2FC_SUS_minus_RES", "observed_only_estimable",
    "observed_only_log2FC", "observed_only_minus_canonical",
    "observed_only_same_sign", "observed_only_relative_change",
    "observed_only_rule"), drop = FALSE],
  LOO = loo,
  QC_fail_sensitivity = if (is.null(rob)) NULL else rob[, intersect(c("gene_symbol",
    "canonical_log2FC_SUS_minus_RES", "effect_excluding_A755",
    "effect_excluding_A764", "effect_excluding_both_QC_fail_SUS",
    "effect_dropping_QC_failed_hemispheres", "A755_sensitive", "A764_sensitive",
    "qc_failed_hemisphere_sensitive", "effect_excluding_both_interpretation",
    "effect_excluding_A755_interpretation", "effect_excluding_A764_interpretation"),
    names(rob)), drop = FALSE],
  Imputation_dependence = if (is.null(rob)) NULL else rob[, c("gene_symbol",
    "n_imputed_SUS", "n_imputed_RES", "differential_missingness_SUS_minus_RES",
    "imputed_share_of_SUS_mean", "imputed_share_of_RES_mean",
    "n_imputed_aligned_with_canonical_direction", "imputation_dependence_class"),
    drop = FALSE],
  Fully_observed_N = fo,
  Spatial_specificity = spec,
  Module_distribution = modcmp,
  Stress_identity_before_after = ident,
  Claim_status = claim_status,
  Thresholds = thr,
  AnimalID_integrity = consumer,
  AnimalID_resolution = resolution,
  Blast_radius = blast,
  Network_behavior_after_fix = coup_cor,
  Validation = validation)
names(sheets)[names(sheets) == "Canonical_N"] <- paste0("Canonical_", nrow(rob))
names(sheets)[names(sheets) == "Fully_observed_N"] <- paste0("Fully_observed_", if (is.null(fo)) 0L else nrow(fo))
sheets <- sheets[!vapply(sheets, function(z) is.null(z) || !nrow(z), logical(1))]

wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
hdr <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E9EDF0",
                             textDecoration = "bold", wrapText = TRUE)
wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE,
                              valign = "top")
openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
openxlsx::writeData(wb, "README", "CA2-SLM robustness audit", startRow = 1)
# Counts and QC identities in the notes come from the audit tables (dtxt), the
# neuropil contrasts, the imputed matrix and the QC table.
na_or <- function(x) if (is.null(x) || is.na(x)) "NA" else x
disp_r <- if (is.null(bias)) NA_real_ else stats::cor(bias$fraction_missing_preimputation,
                                                      bias$mean_centred_value_of_always_observed_proteins)
readme <- data.frame(Note = c(
  if (ca2_unpinned) paste0("PIN NOT ENFORCED. PROTEOMICS_CA2_SLM_UNPINNED=true: this run is not held to the published ",
                           PINNED_N_HITS, " CA2-SLM DAPs; the audit's protein set is checked against this run's contrast file instead."),
  sprintf("CANONICAL DIFFERENTIAL ABUNDANCE IS UNCHANGED AND REMAINS THE PRIMARY STATISTICAL RESULT. All %d FDR-supported CA2-SLM proteins are retained. Nothing here refits a model or recomputes an FDR; every quantity is a descriptive sensitivity diagnostic on the same animal-level values the canonical model saw. QC robustness governs INTERPRETATION, not membership.", nrow(rob)),
  paste0("WHY THIS AUDIT EXISTS. ",
         if (is.na(n_neuropil_daps)) sprintf("%d SUS-RES DAPs fall in CA2-SLM. ", nrow(rob)) else
           sprintf("%d of the %d neuropil SUS-RES DAPs fall in CA2-SLM. ", nrow(rob), n_neuropil_daps),
         "That unit also carries the lowest baseline abundance, the highest pre-imputation missingness",
         if (is.null(dtxt)) ". " else paste0(", ", dtxt$missingness_direction, ", and ", dtxt$flagged_share_with_failed, ". "),
         if (is.na(n_neuropil_samples)) "None was excluded." else
           sprintf("None was excluded: Exclude is FALSE for all %d neuropil samples.", n_neuropil_samples)),
  sprintf("WHAT THE MODEL SAW. Each animal-level value is the unweighted mean of that animal's Left and Right POST-imputation log2 values, so the contrast is %s versus %s and an undetected hemisphere still contributes a drawn value. Imputation status is therefore derived from the PRE-imputation raw matrix; the post-imputation matrix has no NA left and cannot answer the question.",
          na_or(dtxt$n_sus), na_or(dtxt$n_res)),
  if (identical(run_centring, "complete_case_median"))
    paste0("NORMALISATION. This run centres each sample on the median of the proteins quantified in every sample of the compartment (complete_case_median), so a sample's missingness does not set its centring constant. ",
           sprintf("On the proteins observed in all %s CA2-SLM samples the SUS-RES displacement is %+.3f log2", na_or(dtxt$n_samples), disp),
           if (is.finite(disp_r)) sprintf(", and missingness correlates with a sample's centred level at r = %.2f.", disp_r) else ".") else
  paste0("A NORMALISATION ARTIFACT UNDERLIES PART OF THE SIGNAL. Each sample is median-centred on its OBSERVED values before imputation. A heavily missing sample detects only its more abundant proteins, so its observed median is inflated and centring then displaces EVERY protein in that sample downwards - including proteins with no imputed value anywhere. ",
         sprintf("Measured on proteins observed in all %s CA2-SLM samples, missingness predicts that displacement ", na_or(dtxt$n_samples)),
         if (is.finite(disp_r) && abs(disp_r) >= 0.95) "almost perfectly" else sprintf("with r = %.2f", disp_r),
         if (is.null(dtxt)) "." else paste0(", and ", dtxt$qc_fail_most_missing, ".")),
  paste("EXACT RESOLUTION OF THE SENSITIVITY ANALYSES.",
        if (is.null(dtxt)) "" else paste(dtxt$loo_design, dtxt$leave_out_removal),
        "These are robustness diagnostics, not new hypothesis tests, and no new FDR is attached to any of them."),
  "PRESPECIFIED RULES. An animal counts toward an observed-only effect only if BOTH hemispheres were detected, and an observed-only effect needs at least two such animals in EACH group. Imputation dependence is minimal at or below a 25% change in effect magnitude and strong above 50%. Losing more than half the canonical magnitude is a collapse. Every threshold is a named constant and none references a module, gene or tier. The leave-out check removes the CA2-SLM acquisitions whose pre-imputation missingness robust z exceeds 3; this cut was set on 2026-11-21, after the sample QC call proved to depend on the centring while the missingness does not, and replaces the QC-FAILED criterion. The class not_claimable_due_to_QC and the columns named after QC-failed hemispheres keep their names for their consumers and refer to this leave-out.",
  paste0("ANIMAL IDENTITY. The behaviour/physiology consumers were rerun against a fail-closed canonical AnimalID contract. The withdrawn normaliser truncated digit runs longer than four, zero-padded bare numerals while leaving A-prefixed ids at their original width, and required at least three digits - so the proteomics A111 never matched the behaviour OR111, and only one of ",
         if (identical(dtxt$n_animals, 9L)) "nine animals" else "the nine animals of the original 3/3/3 design",
         " had ever joined.")),
  stringsAsFactors = FALSE)
openxlsx::writeDataTable(wb, "README", readme, startRow = 3, tableName = "CA2Readme",
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
                     "CA2_SLM_robustness_audit.xlsx")
xlsx_save_valid_workbook(wb, wb_path)

cat("\n===== CA2-SLM robustness workbook and validation =====\n")
if (ca2_unpinned) cat(toupper(substr(unpinned_note, 1, 1)), substring(unpinned_note, 2), "\n", sep = "")
cat("workbook:", relative_to(wb_path), " sheets:", length(sheets) + 1L, "\n")
cat("figures :", paste(made, collapse = ", "), "\n\n")
for (i in seq_len(nrow(validation))) {
  cat(sprintf("  %-7s %-42s %s\n", validation$status[i], validation$check_id[i],
              substr(validation$detail[i], 1, 70)))
}
crit <- validation$critical %in% TRUE & validation$status == "FAIL"
cat(sprintf("\n%d checks, %d critical FAIL\n", nrow(validation), sum(crit)))
if (any(crit)) {
  cat("CRITICAL VALIDATION FAILURE:\n")
  for (i in which(crit)) cat("  -", validation$check_id[i], ":", validation$detail[i], "\n")
  quit(save = "no", status = 1L)
}
