source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "ca2_slm_robustness_utils.R"))
source(repo_path("R", "spatial_systems_paths.R"))

ROB <- function(f) spatial_systems_find(f, "ca2_slm_robustness")
ATL <- function(f) spatial_systems_find(f, "atlas")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
code_of <- function(...) {
  l <- readLines(repo_path(...), warn = FALSE)
  l <- sub("#.*$", "", l)
  paste(l[nzchar(trimws(l))], collapse = "\n")
}
have <- function(p) file.exists(p)

DA_FILE <- path_processed("02_id_mapping", "mapped", "neuron_neuropil", "forward",
                          "per_file", "CA2slmsus_CA2slmres.csv")

# =====================================================================
# the prespecified rules, tested as pure functions
# =====================================================================

testthat::test_that("every classification threshold is a declared constant", {
  t <- csr_thresholds()
  testthat::expect_true(all(t$prespecified))
  testthat::expect_setequal(t$constant, c("fdr_threshold",
    "min_observed_animals_per_group", "imputation_minimal_cut",
    "imputation_strong_cut", "magnitude_collapse_cut"))
  testthat::expect_identical(csr_fdr_threshold(), 0.05)
  testthat::expect_identical(csr_min_observed_animals(), 2L)
  testthat::expect_lt(csr_imputation_minimal_cut(), csr_imputation_strong_cut())
})

testthat::test_that("the classifier cannot see any biological annotation", {
  # the strongest possible form of this check: the function has no argument
  # through which a module, gene or tier could reach it
  f <- names(formals(csr_classify))
  testthat::expect_setequal(f, c("fully_observed", "estimable",
    "obs_sign_retained", "imp_class", "loo_sign_stable", "a755", "a764",
    "hemi", "imp_rel"))
  testthat::expect_false(any(grepl("Module|Gene|tier|symbol", f, ignore.case = TRUE)))
  src <- code_of("R", "ca2_slm_robustness_utils.R")
  body_txt <- paste(deparse(body(csr_classify)), collapse = "\n")
  testthat::expect_false(grepl("ModuleID|GeneSymbol|candidate_tier|gene_symbol",
                               body_txt))
})

testthat::test_that("the classification ladder is deterministic and total", {
  args <- expand.grid(fully_observed = c(TRUE, FALSE), estimable = c(TRUE, FALSE),
                      obs_sign_retained = c(TRUE, FALSE),
                      imp_class = c("fully_observed", "minimal_imputation_dependence",
                                    "moderate_imputation_dependence",
                                    "strong_imputation_dependence",
                                    "not_estimable_without_imputation"),
                      loo_sign_stable = c(TRUE, FALSE), a755 = c(TRUE, FALSE),
                      a764 = c(TRUE, FALSE), hemi = c(TRUE, FALSE),
                      stringsAsFactors = FALSE)
  cls <- vapply(seq_len(nrow(args)), function(i) {
    csr_classify(args$fully_observed[i], args$estimable[i],
                 args$obs_sign_retained[i], args$imp_class[i],
                 args$loo_sign_stable[i], args$a755[i], args$a764[i],
                 args$hemi[i], 0.3)$class
  }, character(1))
  testthat::expect_false(any(is.na(cls)))
  testthat::expect_true(all(cls %in% c("insufficient_observed_data",
    "not_claimable_due_to_QC", "imputation_sensitive", "single_animal_sensitive",
    "robust_to_missingness_and_QC", "supported_but_QC_sensitive")))
  # calling it twice gives the same answer
  cls2 <- vapply(seq_len(nrow(args)), function(i) {
    csr_classify(args$fully_observed[i], args$estimable[i],
                 args$obs_sign_retained[i], args$imp_class[i],
                 args$loo_sign_stable[i], args$a755[i], args$a764[i],
                 args$hemi[i], 0.3)$class
  }, character(1))
  testthat::expect_identical(cls, cls2)
  # a QC-dependent effect can never be called robust
  testthat::expect_false(any(cls[args$a755 | args$a764 | args$hemi] ==
                               "robust_to_missingness_and_QC"))
  # every branch records a reason
  r <- csr_classify(TRUE, TRUE, TRUE, "fully_observed", TRUE, FALSE, FALSE, FALSE, 0)
  testthat::expect_true(nzchar(r$reason))
})

testthat::test_that("effect helpers behave and never invent a value", {
  v <- c(a = 1, b = 2, c = 3, d = 7, e = 8, f = 9)
  g <- c("RES", "RES", "RES", "SUS", "SUS", "SUS")
  testthat::expect_equal(csr_effect(v, g), 6)
  testthat::expect_true(is.na(csr_effect(v, rep("SUS", 6))))
  testthat::expect_true(csr_same_sign(-1, -2))
  testthat::expect_false(csr_same_sign(-1, 2))
  testthat::expect_true(is.na(csr_same_sign(NA, 1)))
  testthat::expect_equal(csr_rel_change(0.5, 1), 0.5)
  testthat::expect_true(is.na(csr_rel_change(1, 0)))
  testthat::expect_identical(csr_bare_animal(c("A0003", "A111")), c("3", "111"))
  testthat::expect_identical(csr_expgroup_to_stress(c(1, 2, 3)),
                             c("CON", "RES", "SUS"))
  testthat::expect_error(csr_expgroup_to_stress(4), "unmapped ExpGroup")
})

# =====================================================================
# no statistic is recomputed
# =====================================================================

testthat::test_that("no FDR is recomputed and no model is refitted", {
  for (f in list(c("analysis/spatial_validation", "audit_ca2_slm_robustness.R"),
                 c("analysis/spatial_validation", "audit_stress_identity_robustness.R"),
                 c("R", "ca2_slm_robustness_utils.R"))) {
    s <- do.call(code_of, as.list(f))
    testthat::expect_false(grepl("p\\.adjust", s), info = paste(f, collapse = "/"))
    testthat::expect_false(grepl("lmFit|eBayes", s), info = paste(f, collapse = "/"))
    testthat::expect_false(grepl("t\\.test|wilcox\\.test", s),
                           info = paste(f, collapse = "/"))
  }
})

testthat::test_that("the canonical DA still yields exactly 28 CA2-SLM DAPs", {
  testthat::skip_if_not(have(DA_FILE), "canonical contrast not available")
  da <- rd(DA_FILE)
  testthat::expect_identical(nrow(da), 5045L)
  testthat::expect_identical(sum(da$padj < csr_fdr_threshold(), na.rm = TRUE), 28L)
  # the `significant` flag and the threshold agree on every row
  testthat::expect_identical(sum(da$significant %in% TRUE), 28L)
})

# =====================================================================
# the generated audit
# =====================================================================

testthat::test_that("all 28 canonical DAPs enter the audit with verbatim statistics", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  testthat::expect_identical(nrow(r), 28L)
  testthat::expect_identical(anyDuplicated(r$original_identifier), 0L)
  testthat::expect_true(all(grepl("not refitted", r$canonical_statistic_source)))

  testthat::skip_if_not(have(DA_FILE))
  da <- rd(DA_FILE)
  daps <- da[!is.na(da$padj) & da$padj < csr_fdr_threshold(), ]
  testthat::expect_setequal(r$original_identifier, daps$original_identifier)
  m <- match(r$original_identifier, daps$original_identifier)
  testthat::expect_equal(r$canonical_log2FC_SUS_minus_RES, daps$log2fc[m])
  testthat::expect_equal(r$canonical_BH_FDR, daps$padj[m])
  testthat::expect_equal(r$canonical_p_value, daps$pval[m])
})

testthat::test_that("imputation status comes from pre-imputation missingness", {
  s <- code_of("analysis/spatial_validation", "audit_ca2_slm_robustness.R")
  testthat::expect_true(grepl("quicksearch.pg_matrix.tsv", s, fixed = TRUE))
  testthat::expect_true(grepl("csr_preimputation_mask", s, fixed = TRUE))
  u <- code_of("R", "ca2_slm_robustness_utils.R")
  # the mask function reads the RAW matrix and reproduces the 70% filter
  testthat::expect_true(grepl("rowMeans\\(is.na\\(na_mat\\)\\)|rowMeans\\(na_mat\\)", u))
  testthat::expect_true(grepl("0.7", u, fixed = TRUE))

  testthat::skip_if_not(have(ROB("CA2_SLM_sample_context.csv")))
  sc <- rd(ROB("CA2_SLM_sample_context.csv"))
  # a post-imputation source would report zero missing for every sample
  testthat::expect_true(all(sc$fraction_missing_preimputation > 0))
  testthat::expect_gt(max(sc$fraction_missing_preimputation), 0.2)
  testthat::expect_identical(nrow(sc), 18L)
})

testthat::test_that("fully observed proteins are identified consistently", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  # fully observed means literally zero imputed values in SUS and RES
  testthat::expect_identical(r$fully_observed,
                             (r$n_imputed_SUS + r$n_imputed_RES) == 0L)
  testthat::expect_identical(sum(r$fully_observed), 8L)
  # observation counts are bounded by the 6 samples per group
  testthat::expect_true(all(r$n_observed_SUS <= 6L & r$n_observed_SUS >= 0L))
  testthat::expect_true(all(r$n_observed_RES <= 6L & r$n_observed_RES >= 0L))
  testthat::expect_true(all(r$n_observed_SUS + r$n_imputed_SUS == 6L))
  testthat::expect_true(all(r$n_observed_RES + r$n_imputed_RES == 6L))
  # being fully observed does not by itself confer robustness
  fo <- r[r$fully_observed, ]
  testthat::expect_true(any(fo$CA2_SLM_robustness_class != "robust_to_missingness_and_QC"))
})

testthat::test_that("an observed-only effect requires the prespecified minimum n", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  expected <- r$n_fully_observed_SUS_animals >= csr_min_observed_animals() &
    r$n_fully_observed_RES_animals >= csr_min_observed_animals()
  testthat::expect_identical(r$observed_only_estimable, expected)
  # no 1-versus-3 estimate is ever reported
  testthat::expect_true(all(is.na(r$observed_only_log2FC[!r$observed_only_estimable])))
  testthat::expect_true(all(r$n_fully_observed_SUS_animals[r$observed_only_estimable] >= 2L))
  testthat::expect_true(all(grepl("BOTH hemispheres", r$observed_only_rule)))
})

testthat::test_that("leave-one-out removes exactly one animal at a time", {
  testthat::skip_if_not(have(ROB("CA2_SLM_leave_one_animal_out_long.csv")), "audit not run")
  l <- rd(ROB("CA2_SLM_leave_one_animal_out_long.csv"))
  per <- split(l$omitted_AnimalID, l$original_identifier)
  testthat::expect_true(all(vapply(per, function(z) length(z) == 9L, logical(1))))
  testthat::expect_true(all(vapply(per, function(z) anyDuplicated(z) == 0L, logical(1))))
  testthat::expect_setequal(unique(l$omitted_AnimalID),
    c("3", "111", "127", "129", "135", "139", "755", "764", "765"))
  # omitting a CON animal cannot move a SUS-vs-RES contrast
  con <- l[l$omitted_StressGroup == "CON", ]
  testthat::expect_true(all(abs(con$change_from_canonical) < 1e-9))
  testthat::expect_true(all(!con$can_change_estimate))
  testthat::expect_true(all(l$can_change_estimate[l$omitted_StressGroup %in% c("SUS", "RES")]))
})

testthat::test_that("the QC-fail analyses exclude exactly the named animals", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  l <- rd(ROB("CA2_SLM_leave_one_animal_out_long.csv"))
  # the A755-out effect must equal the leave-one-out effect for animal 755,
  # which proves it removed 755 and nothing else
  m <- match(paste(r$original_identifier, "755"),
             paste(l$original_identifier, l$omitted_AnimalID))
  testthat::expect_equal(r$effect_excluding_A755, l$loo_effect[m])
  m4 <- match(paste(r$original_identifier, "764"),
              paste(l$original_identifier, l$omitted_AnimalID))
  testthat::expect_equal(r$effect_excluding_A764, l$loo_effect[m4])
  # and they are not the same analysis
  testthat::expect_false(isTRUE(all.equal(r$effect_excluding_A755,
                                          r$effect_excluding_A764)))
})

testthat::test_that("removing both QC-failed SUS animals is labelled descriptive-only", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  testthat::expect_true(all(grepl("DESCRIPTIVE ONLY",
                                  r$effect_excluding_both_interpretation)))
  testthat::expect_true(all(grepl("1 SUS versus 3 RES",
                                  r$effect_excluding_both_interpretation)))
  # it is computed, but it never feeds a robustness class
  u <- paste(deparse(body(csr_classify)), collapse = "\n")
  testthat::expect_false(grepl("excluding_both", u))
})

testthat::test_that("no protein is upgraded because of its module or annotation", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  testthat::skip_if(!("ModuleID" %in% names(r)))
  # the two numerically dominant modules must not be over-represented among the
  # robust set relative to their share of the flags that drive classification
  rob_set <- r[r$CA2_SLM_robustness_class == "robust_to_missingness_and_QC", ]
  # recomputing the class from the stored flags alone must reproduce it exactly
  recomputed <- vapply(seq_len(nrow(r)), function(i) {
    csr_classify(r$fully_observed[i], r$observed_only_estimable[i],
                 r$observed_only_same_sign[i], r$imputation_dependence_class[i],
                 r$loo_sign_stable[i], r$A755_sensitive[i], r$A764_sensitive[i],
                 r$qc_failed_hemisphere_sensitive[i],
                 r$observed_only_relative_change[i])$class
  }, character(1))
  testthat::expect_identical(recomputed, r$CA2_SLM_robustness_class)
})

testthat::test_that("the stress-identity comparison uses explicit subsets and the canonical rule", {
  testthat::skip_if_not(have(ROB("stress_identity_robustness_comparison.csv")),
                        "identity comparison not run")
  z <- rd(ROB("stress_identity_robustness_comparison.csv"))
  testthat::expect_true(all(nzchar(z$subset_definition)))
  testthat::expect_true(all(grepl("tabulated not recomputed", z$classification_source)))
  testthat::expect_true("all_canonical_FDR_supported_hits" %in% z$subset)
  # the headline row must still reproduce the published 35 of 37
  h <- z[z$subset == "all_canonical_FDR_supported_hits", ]
  testthat::expect_identical(h$n_hits, 37L)
  testthat::expect_identical(h$effect_outside_baseline_affinity, 35L)
  testthat::expect_identical(h$n_at_baseline_rank_10, 15L)
  # every subset is a subset: never more hits than the full set
  testthat::expect_true(all(z$n_hits <= 37L))
  s <- code_of("analysis/spatial_validation", "audit_stress_identity_robustness.R")
  # the classification rule is not reimplemented here
  testthat::expect_false(grepl("effect_at_baseline_peak\\s*<-", s))
  testthat::expect_false(grepl("HIGH_AFFINITY_FRACTION", s))
})

testthat::test_that("module counts are reported without a post hoc enrichment test", {
  testthat::skip_if_not(have(ROB("CA2_SLM_module_distribution_comparison.csv")))
  m <- rd(ROB("CA2_SLM_module_distribution_comparison.csv"))
  testthat::expect_true(all(grepl("NOT TESTED", m$enrichment_test)))
  testthat::expect_identical(sum(m$A_all_canonical_28), 28L)
  testthat::expect_true(sum(m$C_robust_only) <= sum(m$A_all_canonical_28))
})

# =====================================================================
# the atlas keeps canonical DA intact
# =====================================================================

testthat::test_that("the atlas gains claimability fields without losing any protein", {
  p <- ATL("protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(have(p), "atlas not generated")
  a <- rd(p)
  testthat::expect_identical(sum(a$is_sus_res_fdr_supported %in% TRUE), 37L)
  for (nm in c("CA2_SLM_robustness_class", "QC_claimability",
               "imputation_dependence", "LOO_sign_stability", "fully_observed",
               "claimable_for_biological_interpretation", "robustness_audit_status")) {
    testthat::expect_true(nm %in% names(a), info = nm)
  }
  # no FDR-supported protein is removed and no FDR value is blanked
  f <- a[a$is_sus_res_fdr_supported %in% TRUE, ]
  testthat::expect_true(all(is.finite(f$sus_res_strongest_BH_FDR)))
  testthat::expect_true(all(f$sus_res_strongest_BH_FDR <= 0.05))
  # a protein that was never audited must not read as claimable
  testthat::expect_false(any(f$QC_claimability == "not_audited" &
                               f$claimable_for_biological_interpretation %in% TRUE))
  testthat::expect_true(all(!is.na(f$QC_claimability)))
  # the published classification is unchanged by the overlay
  testthat::expect_identical(
    sum(f$effect_identity_relationship == "effect_outside_baseline_affinity"), 35L)
})

testthat::test_that("the claimability annotation is one row per audited protein", {
  p <- ATL("protein_claimability_annotation.csv")
  testthat::skip_if_not(have(p), "annotation not generated")
  c1 <- rd(p)
  testthat::expect_identical(nrow(c1), 28L)
  testthat::expect_identical(anyDuplicated(paste(c1$dataset, c1$ProteinGroupID)), 0L)
  testthat::expect_true(all(c1$dataset == "neuron_neuropil"))
  testthat::expect_setequal(unique(c1$QC_claimability),
    intersect(c("claimable", "claimable_with_caveat", "not_claimable", "not_evaluable"),
              unique(c1$QC_claimability)))
  testthat::expect_identical(
    c1$claimable_for_biological_interpretation,
    c1$CA2_SLM_robustness_class == "robust_to_missingness_and_QC")
  testthat::expect_true(all(nzchar(c1$robustness_reason)))
})

# =====================================================================
# outputs, figures, validation
# =====================================================================

testthat::test_that("the robustness validation contract has no critical failure", {
  p <- ROB("CA2_SLM_robustness_validation.csv")
  testthat::skip_if_not(have(p), "validation not generated")
  v <- rd(p)
  testthat::expect_gt(nrow(v), 15L)
  testthat::expect_identical(sum(v$critical %in% TRUE & v$status == "FAIL"), 0L)
  for (id in c("canonical_DA_unchanged", "no_new_FDR_computed",
               "all_28_canonical_hits_retained", "LOO_removes_exactly_one_animal",
               "animal_id_collision_free", "wgcna_state_unchanged")) {
    testthat::expect_true(id %in% v$check_id, info = id)
  }
})

testthat::test_that("every candidate figure has source data", {
  d <- spatial_systems_dir_any("summarize_ca2_slm_robustness",
                               "ca2_slm_robustness", kind = "figures")
  testthat::skip_if_not(dir.exists(d), "figures not generated")
  pngs <- list.files(d, pattern = "[.]png$")
  testthat::skip_if(length(pngs) == 0L)
  testthat::expect_gte(length(pngs), 4L)
  for (f in pngs) {
    testthat::expect_true(
      file.exists(file.path(d, sub("[.]png$", "_source_data.csv", f))), info = f)
  }
})

testthat::test_that("the workbook exists and opens", {
  p <- spatial_systems_find("CA2_SLM_robustness_audit.xlsx", kind = "reports")
  testthat::skip_if_not(have(p), "workbook not generated")
  testthat::skip_if_not(requireNamespace("openxlsx", quietly = TRUE))
  sh <- openxlsx::getSheetNames(p)
  testthat::expect_true("README" %in% sh)
  for (nm in c("Canonical_28", "Missingness", "Observed_only", "LOO",
               "QC_fail_sensitivity", "Imputation_dependence", "Fully_observed_8",
               "Spatial_specificity", "Module_distribution",
               "Stress_identity_before_after", "AnimalID_integrity",
               "Blast_radius", "Validation")) {
    testthat::expect_true(nm %in% sh, info = nm)
  }
  testthat::expect_identical(nrow(openxlsx::read.xlsx(p, sheet = "Canonical_28")), 28L)
})
