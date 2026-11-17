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

source(repo_path("R", "preprocessing_paths.R"))
# the canonical contrast, where the audit reads it (normalized namespace, legacy fallback)
DA_FILE <- file.path(preprocessing_mapped_contrast_dir("neuron_neuropil", "forward"), "CA2slmsus_CA2slmres.csv")

# PROTEOMICS_CA2_SLM_UNPINNED=true marks a rerun on other data (the 3/3/2 rerun):
# the published counts (28 DAPs, 9 animals, 18 samples) are then not expected,
# and the data tests compare against the audit's own tables instead.
unpinned <- identical(tolower(Sys.getenv("PROTEOMICS_CA2_SLM_UNPINNED", "")), "true")

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

# ---------------------------------------- the design present (3/3/3 and 3/3/2)

# Sample context of a synthetic CA2-SLM design: two hemispheres per animal,
# QC flags as in the published data (A755 L and A764 R FAIL, A765 R WARN).
synthetic_ctx <- function(with_764 = TRUE) {
  grp <- c("127" = "CON", "129" = "CON", "765" = "CON", "3" = "RES", "135" = "RES",
           "139" = "RES", "111" = "SUS", "755" = "SUS", "764" = "SUS")
  if (!with_764) grp <- grp[names(grp) != "764"]
  s <- data.frame(AnimalID = rep(names(grp), each = 2L), hemisphere = c("L", "R"),
                  stringsAsFactors = FALSE)
  s$StressGroup <- unname(grp[s$AnimalID])
  s$sample_id <- paste0("S", s$AnimalID, s$hemisphere)
  key <- paste0(s$AnimalID, s$hemisphere)
  s$qc_flag <- "PASS"
  s$qc_flag[key %in% c("755L", "764R")] <- "FAIL"
  if (with_764) s$qc_flag[key == "765R"] <- "WARN"
  s$fraction_missing_preimputation <- 0.05 + 0.005 * seq_len(nrow(s))
  s$fraction_missing_preimputation[key == "755L"] <- 0.25
  s$fraction_missing_preimputation[key == "764R"] <- 0.16
  s
}

testthat::test_that("a QC-fail removal leaving fewer than two animals in a group is descriptive only", {
  testthat::expect_true(csr_removal_needs_label(1L, 3L))
  testthat::expect_true(csr_removal_needs_label(3L, 1L))
  testthat::expect_false(csr_removal_needs_label(2L, 3L))
  testthat::expect_false(csr_removal_needs_label(2L, 2L))
  # the published label is reproduced exactly
  testthat::expect_identical(csr_removal_interpretation(1L, 3L),
    "DESCRIPTIVE ONLY: 1 SUS versus 3 RES is not a valid group comparison")
  testthat::expect_false(grepl("DESCRIPTIVE ONLY", csr_removal_interpretation(2L, 3L)))
  testthat::expect_match(csr_removal_interpretation(3L, 3L, computed = FALSE), "^not computed")
})

testthat::test_that("per-protein QC-fail removals follow the animals in the design", {
  one_protein <- function(with_764) {
    s <- synthetic_ctx(with_764)
    animals <- sort(unique(s$AnimalID))
    stress <- stats::setNames(s$StressGroup[match(animals, s$AnimalID)], animals)
    val <- c(CON = 0.5, RES = 0, SUS = -1)[s$StressGroup] + 0.01 * seq_len(nrow(s))
    imputed <- matrix(val, nrow = 1L, dimnames = list("P1", s$sample_id))
    mask <- stats::setNames(rep(FALSE, nrow(s)), s$sample_id)
    av <- tapply(val, s$AnimalID, mean)
    fc <- mean(av[names(stress)[stress == "SUS"]]) - mean(av[names(stress)[stress == "RES"]])
    ctx <- list(sample_meta = s, animals = animals, stress = stress,
                qc_fail_samples = s$sample_id[s$qc_flag == "FAIL"],
                sample_of = function(a, h) s$sample_id[s$AnimalID == a & s$hemisphere == h][1])
    canon <- data.frame(ProteinGroupID = "PG1", official_gene_symbol = "G1", log2fc = fc,
                        pval = 1e-4, padj = 0.01, stringsAsFactors = FALSE)
    csr_protein_robustness("P1", canon, imputed, mask, ctx)$summary
  }
  full <- one_protein(TRUE)
  testthat::expect_false(grepl("DESCRIPTIVE ONLY", full$effect_excluding_A755_interpretation))
  testthat::expect_false(grepl("DESCRIPTIVE ONLY", full$effect_excluding_A764_interpretation))
  testthat::expect_identical(full$effect_excluding_both_interpretation,
    "DESCRIPTIVE ONLY: 1 SUS versus 3 RES is not a valid group comparison")
  testthat::expect_true(is.finite(full$effect_excluding_A764))

  rerun <- one_protein(FALSE)
  # without A764, removing A755 leaves one SUS animal
  testthat::expect_identical(rerun$effect_excluding_A755_interpretation,
    "DESCRIPTIVE ONLY: 1 SUS versus 3 RES is not a valid group comparison")
  testthat::expect_match(rerun$effect_excluding_both_interpretation, "^DESCRIPTIVE ONLY")
  testthat::expect_equal(rerun$effect_excluding_both_QC_fail_SUS, rerun$effect_excluding_A755)
  # an absent animal is not removed and no effect is reported for it
  testthat::expect_true(is.na(rerun$effect_excluding_A764))
  testthat::expect_match(rerun$effect_excluding_A764_interpretation, "^not computed")
  testthat::expect_false(rerun$A764_sensitive)
  # the reasons name the QC-failed samples of the design present
  testthat::expect_identical(full$CA2_SLM_robustness_class, "robust_to_missingness_and_QC")
  testthat::expect_match(full$classification_reason, "not dependent on either QC-failed SUS sample", fixed = TRUE)
  testthat::expect_identical(rerun$CA2_SLM_robustness_class, "robust_to_missingness_and_QC")
  testthat::expect_match(rerun$classification_reason, "not dependent on the QC-failed SUS sample", fixed = TRUE)
})

testthat::test_that("generated wording follows the design in the sample context", {
  t9 <- csr_design_text(synthetic_ctx(TRUE), n_flagged_neuropil = 5L)
  # the published 3/3/3 wording is reproduced
  testthat::expect_identical(t9$flagged_share_with_flags,
    "3 of the 5 QC-flagged neuropil samples (A755 Left FAIL, A764 Right FAIL, A765 Right WARN)")
  testthat::expect_identical(t9$flagged_share_with_failed,
    "3 of the 5 QC-flagged neuropil samples - two of them QC-FAILED SUS acquisitions (A755 Left, A764 Right)")
  testthat::expect_identical(t9$qc_fail_sus_samples, "two QC-failed SUS samples")
  testthat::expect_identical(t9$qc_fail_most_missing,
                             "the two QC-FAILED SUS samples are the two most missing")
  testthat::expect_identical(t9$loo_design,
                             "Leave-one-animal-out on a 3-versus-3 design leaves 2 versus 3.")
  testthat::expect_identical(t9$qc_fail_removal, paste0(
    "Removing both QC-failed SUS animals leaves 1 versus 3, which is NOT a valid group ",
    "comparison and is labelled descriptive-only wherever it appears."))
  testthat::expect_identical(t9$missingness_direction, "higher missingness in SUS than RES")
  testthat::expect_identical(c(t9$n_samples, t9$n_sus_res_samples, t9$n_animals), c(18L, 12L, 9L))

  t8 <- csr_design_text(synthetic_ctx(FALSE), n_flagged_neuropil = 3L)
  testthat::expect_identical(t8$flagged_share_with_flags,
                             "1 of the 3 QC-flagged neuropil samples (A755 Left FAIL)")
  testthat::expect_identical(t8$flagged_share_with_failed,
    "1 of the 3 QC-flagged neuropil samples - a QC-FAILED SUS acquisition (A755 Left)")
  testthat::expect_identical(t8$qc_fail_sus_samples, "one QC-failed SUS sample")
  testthat::expect_identical(t8$qc_fail_most_missing, "the QC-FAILED SUS sample is the most missing")
  testthat::expect_identical(t8$loo_design, paste0(
    "Leave-one-animal-out on a 2-versus-3 design leaves 1 versus 3 or 2 versus 2. A leave-out that leaves ",
    "fewer than 2 animals in a group is descriptive only; it enters the classes only through the requirement ",
    "that the sign is stable under every leave-out."))
  testthat::expect_match(t8$qc_fail_removal,
    "^Removing the QC-failed SUS animal leaves 1 versus 3, which is NOT a valid group comparison")
  # in 3/3/2 that removal still marks proteins not claimable, and the text says so
  testthat::expect_match(t8$qc_fail_removal, "still marks the protein not claimable", fixed = TRUE)
  testthat::expect_false(grepl("wherever it appears", t8$qc_fail_removal, fixed = TRUE))
  testthat::expect_identical(c(t8$n_samples, t8$n_sus_res_samples, t8$n_animals), c(16L, 10L, 8L))
  # without the neuropil QC table the share is stated for CA2-SLM alone
  testthat::expect_identical(csr_design_text(synthetic_ctx(FALSE))$flagged_share_with_flags,
                             "1 QC-flagged sample (A755 Left FAIL)")

  testthat::expect_identical(csr_and_list(c("A755", "A764")), "A755 and A764")
  testthat::expect_identical(csr_and_list("A755"), "A755")
  testthat::expect_identical(csr_number_word(c(0, 2, 12)), c("no", "two", "12"))
})

testthat::test_that("the fully observed audit names the QC-failed animals present", {
  note_for <- function(with_764, all_pass = FALSE) {
    s <- synthetic_ctx(with_764)
    if (all_pass) s$qc_flag <- "PASS"
    fully <- data.frame(original_identifier = "P1", canonical_log2FC_SUS_minus_RES = -1)
    loo <- data.frame(original_identifier = "P1", omitted_AnimalID = unique(s$AnimalID),
                      can_change_estimate = TRUE, change_from_canonical = 0.1)
    csr_fully_observed_audit(fully, loo, s)$audit_note
  }
  testthat::expect_match(note_for(TRUE), "but A755 and A764 remain QC-FAILED acquisitions;", fixed = TRUE)
  testthat::expect_match(note_for(FALSE), "but A755 remains a QC-FAILED acquisition;", fixed = TRUE)
  testthat::expect_match(note_for(FALSE, all_pass = TRUE), "and no acquisition is QC-FAILED;", fixed = TRUE)
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
  testthat::skip_if(unpinned, "published count not enforced (PROTEOMICS_CA2_SLM_UNPINNED=true)")
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

testthat::test_that("every canonical DAP enters the audit with verbatim statistics", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  if (!unpinned) testthat::expect_identical(nrow(r), 28L)   # the published count
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
  testthat::expect_identical(nrow(sc), 2L * length(unique(sc$AnimalID)))
  if (!unpinned) testthat::expect_identical(nrow(sc), 18L)
})

testthat::test_that("fully observed proteins are identified consistently", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  # fully observed means literally zero imputed values in SUS and RES
  testthat::expect_identical(r$fully_observed,
                             (r$n_imputed_SUS + r$n_imputed_RES) == 0L)
  if (!unpinned) testthat::expect_identical(sum(r$fully_observed), 8L)
  # observation counts are bounded by the samples per group (6 in the 3/3/3 design)
  n_sus <- 6L; n_res <- 6L
  if (unpinned) {
    testthat::skip_if_not(have(ROB("CA2_SLM_sample_context.csv")))
    sc <- rd(ROB("CA2_SLM_sample_context.csv"))
    n_sus <- sum(sc$StressGroup == "SUS"); n_res <- sum(sc$StressGroup == "RES")
  }
  testthat::expect_true(all(r$n_observed_SUS <= n_sus & r$n_observed_SUS >= 0L))
  testthat::expect_true(all(r$n_observed_RES <= n_res & r$n_observed_RES >= 0L))
  testthat::expect_true(all(r$n_observed_SUS + r$n_imputed_SUS == n_sus))
  testthat::expect_true(all(r$n_observed_RES + r$n_imputed_RES == n_res))
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
  n_animals <- length(unique(l$omitted_AnimalID))
  testthat::expect_true(all(vapply(per, function(z) length(z) == n_animals, logical(1))))
  testthat::expect_true(all(vapply(per, function(z) anyDuplicated(z) == 0L, logical(1))))
  if (!unpinned) testthat::expect_setequal(unique(l$omitted_AnimalID),
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
  # without A764 in the design both sides are NA (an all-NA column reads as logical)
  m4 <- match(paste(r$original_identifier, "764"),
              paste(l$original_identifier, l$omitted_AnimalID))
  testthat::expect_equal(as.numeric(r$effect_excluding_A764), as.numeric(l$loo_effect[m4]))
  # and they are not the same analysis
  testthat::expect_false(isTRUE(all.equal(r$effect_excluding_A755,
                                          r$effect_excluding_A764)))
})

testthat::test_that("removing both QC-failed SUS animals is labelled descriptive-only", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  testthat::expect_true(all(grepl("DESCRIPTIVE ONLY",
                                  r$effect_excluding_both_interpretation)))
  if (!unpinned) testthat::expect_true(all(grepl("1 SUS versus 3 RES",
                                                 r$effect_excluding_both_interpretation)))
  # it is computed, but it never feeds a robustness class
  u <- paste(deparse(body(csr_classify)), collapse = "\n")
  testthat::expect_false(grepl("excluding_both", u))
})

testthat::test_that("every QC-fail removal leaving fewer than two animals in a group is labelled", {
  testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")) &&
                          have(ROB("CA2_SLM_sample_context.csv")), "audit not run")
  r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
  sc <- rd(ROB("CA2_SLM_sample_context.csv"))
  testthat::skip_if_not("effect_excluding_A755_interpretation" %in% names(r),
                        "audit written before the per-removal labels")
  sets <- list(A755 = "755", A764 = "764", both = c("755", "764"))
  for (nm in names(sets)) {
    drop <- intersect(sets[[nm]], as.character(sc$AnimalID))
    left <- sc[!as.character(sc$AnimalID) %in% drop, ]
    n_sus <- length(unique(left$AnimalID[left$StressGroup == "SUS"]))
    n_res <- length(unique(left$AnimalID[left$StressGroup == "RES"]))
    lab <- r[[paste0("effect_excluding_", nm, "_interpretation")]]
    eff <- r[[if (nm == "both") "effect_excluding_both_QC_fail_SUS" else paste0("effect_excluding_", nm)]]
    if (!length(drop)) {
      testthat::expect_true(all(is.na(eff)), info = nm)
    } else {
      testthat::expect_identical(all(grepl("DESCRIPTIVE ONLY", lab)),
                                 n_sus < 2L || n_res < 2L, info = nm)
      testthat::expect_identical(any(grepl("DESCRIPTIVE ONLY", lab)),
                                 n_sus < 2L || n_res < 2L, info = nm)
    }
  }
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
  h <- z[z$subset == "all_canonical_FDR_supported_hits", ]
  # the headline row must still reproduce the published 35 of 37
  if (!unpinned) {
    testthat::expect_identical(h$n_hits, 37L)
    testthat::expect_identical(h$effect_outside_baseline_affinity, 35L)
    testthat::expect_identical(h$n_at_baseline_rank_10, 15L)
  }
  # every subset is a subset: never more hits than the full set
  testthat::expect_true(all(z$n_hits <= h$n_hits))
  s <- code_of("analysis/spatial_validation", "audit_stress_identity_robustness.R")
  # the classification rule is not reimplemented here
  testthat::expect_false(grepl("effect_at_baseline_peak\\s*<-", s))
  testthat::expect_false(grepl("HIGH_AFFINITY_FRACTION", s))
})

testthat::test_that("module counts are reported without a post hoc enrichment test", {
  testthat::skip_if_not(have(ROB("CA2_SLM_module_distribution_comparison.csv")))
  m <- rd(ROB("CA2_SLM_module_distribution_comparison.csv"))
  testthat::expect_true(all(grepl("NOT TESTED", m$enrichment_test)))
  if (!unpinned) {
    testthat::expect_identical(sum(m$A_all_canonical_28), 28L)
  } else {
    testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")))
    testthat::expect_identical(sum(m$A_all_canonical_28), nrow(rd(ROB("CA2_SLM_DAP_robustness.csv"))))
  }
  testthat::expect_true(sum(m$C_robust_only) <= sum(m$A_all_canonical_28))
})

# =====================================================================
# the atlas keeps canonical DA intact
# =====================================================================

testthat::test_that("the atlas gains claimability fields without losing any protein", {
  p <- ATL("protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(have(p), "atlas not generated")
  a <- rd(p)
  if (!unpinned) testthat::expect_identical(sum(a$is_sus_res_fdr_supported %in% TRUE), 37L)
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
  if (!unpinned) testthat::expect_identical(
    sum(f$effect_identity_relationship == "effect_outside_baseline_affinity"), 35L)
})

testthat::test_that("the claimability annotation is one row per audited protein", {
  p <- ATL("protein_claimability_annotation.csv")
  testthat::skip_if_not(have(p), "annotation not generated")
  c1 <- rd(p)
  if (!unpinned) {
    testthat::expect_identical(nrow(c1), 28L)
  } else if (have(ROB("CA2_SLM_DAP_robustness.csv"))) {
    testthat::expect_setequal(c1$ProteinGroupID, rd(ROB("CA2_SLM_DAP_robustness.csv"))$ProteinGroupID)
  }
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

testthat::test_that("an atlas written after the annotation carries it for every audited protein", {
  # An empty overlay reads as not audited everywhere (the atlas once looked for the
  # annotation in its own folder). The first pass legitimately precedes the annotation.
  pa <- ATL("protein_spatial_cell_affinity.csv")
  pc <- path_results("spatial_validation", "audit_stress_identity_robustness", "global", "tables",
                     "protein_claimability_annotation.csv")   # the copy the atlas reads
  testthat::skip_if_not(have(pa) && have(pc), "atlas or annotation not generated")
  testthat::skip_if(file.mtime(pa) < file.mtime(pc), "atlas predates the annotation (first pass)")
  a <- rd(pa); c1 <- rd(pc)
  m <- match(paste(c1$dataset, c1$ProteinGroupID), paste(a$dataset, a$ProteinGroupID))
  testthat::expect_false(anyNA(m))
  testthat::expect_identical(sum(a$robustness_audit_status == "audited_ca2_slm_v1", na.rm = TRUE), nrow(c1))
  testthat::expect_identical(a$CA2_SLM_robustness_class[m], c1$CA2_SLM_robustness_class)
  testthat::expect_identical(a$QC_claimability[m], c1$QC_claimability)
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
  # the published run names Canonical_28 and Fully_observed_8; a rerun names its own counts
  n_hits <- 28L; n_fo <- 8L
  if (unpinned) {
    testthat::skip_if_not(have(ROB("CA2_SLM_DAP_robustness.csv")))
    r <- rd(ROB("CA2_SLM_DAP_robustness.csv"))
    n_hits <- nrow(r); n_fo <- sum(r$fully_observed %in% TRUE)
  }
  for (nm in c(paste0("Canonical_", n_hits), "Missingness", "Observed_only", "LOO",
               "QC_fail_sensitivity", "Imputation_dependence", paste0("Fully_observed_", n_fo),
               "Spatial_specificity", "Module_distribution",
               "Stress_identity_before_after", "AnimalID_integrity",
               "Blast_radius", "Validation")) {
    testthat::expect_true(nm %in% sh, info = nm)
  }
  testthat::expect_identical(nrow(openxlsx::read.xlsx(p, sheet = paste0("Canonical_", n_hits))), n_hits)
})
