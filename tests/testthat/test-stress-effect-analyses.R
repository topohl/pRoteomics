source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "spatial", "stress_effect_utils.R"))

# Stress-effect analyses of 2026-10-29 (test_stress_exposure, test_unit_specificity): the exact label test on
# known data, and the invariants of the written tables when they exist.

testthat::test_that("the exact label test enumerates every labelling and is two-sided", {
  v <- c(a = 1, b = 2, c = 3, d = 10, e = 11, f = 12)
  r <- exact_label_test(v, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  testthat::expect_equal(r$n_labellings, choose(6, 3))
  testthat::expect_equal(r$estimate, 9)
  testthat::expect_equal(r$p, 2 / 20)          # the observed split and its mirror image
  r2 <- exact_label_test(c(v, g = NA), c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  testthat::expect_equal(r2$n_labellings, 20)  # missing values leave the test
})

testthat::test_that("the relabelling null holds every assignment of nine animals to three groups of three once", {
  L <- stress_all_labellings(paste0("a", 1:9))
  testthat::expect_equal(nrow(L), 1680L)
  testthat::expect_equal(anyDuplicated(apply(L, 1, paste, collapse = "")), 0L)
  testthat::expect_true(all(apply(L, 1, function(l) all(table(factor(l, c("CON", "RES", "SUS"))) == 3L))))
})

testthat::test_that("with eight animals (3/3/2) the null holds every assignment with the observed group sizes once", {
  a <- paste0("a", 1:8)
  truth <- stats::setNames(c("CON", "CON", "CON", "RES", "RES", "RES", "SUS", "SUS"), a)
  L <- stress_all_labellings(a, truth)
  testthat::expect_equal(nrow(L), 560L)   # 8! / (3! 3! 2!)
  testthat::expect_equal(anyDuplicated(apply(L, 1, paste, collapse = "")), 0L)
  testthat::expect_true(all(apply(L, 1, function(l) identical(as.integer(table(factor(l, c("CON", "RES", "SUS")))), c(3L, 3L, 2L)))))
  testthat::expect_equal(sum(apply(L, 1, function(l) identical(unname(l), unname(truth)))), 1L)
  # giving the 3/3/3 truth reproduces the original enumeration exactly
  a9 <- paste0("a", 1:9)
  testthat::expect_identical(stress_all_labellings(a9, stats::setNames(rep(c("CON", "RES", "SUS"), each = 3), a9)), stress_all_labellings(a9))
})

testthat::test_that("the relabelling fit reproduces the official Figure 3 t at the true labels", {
  testthat::skip_if(!file.exists(animal_level_gct_path("neuron_neuropil")), "animal-level matrices not present")
  d <- stress_gene_frame("neuron_neuropil")
  t <- stress_contrast_t(d, d$true_group)
  ident <- canonical_protein_identity("neuron_neuropil", "CA2slm", nrow(d$m))
  testthat::expect_equal(unname(t[, "CA2slm|SUS_vs_RES"]), ident$t[d$rows], tolerance = 1e-8)
  testthat::expect_equal(anyDuplicated(rownames(t)), 0L)
})

tab <- function(owner, file) repo_path("results", "spatial_validation", owner, "global", "tables", file)
read_if <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL

testthat::test_that("stress exposure tables cover every unit and keep the panel's declared shape", {
  u <- read_if(tab("test_stress_exposure", "stress_exposure_unit_summary.csv"))
  pt <- read_if(tab("test_stress_exposure", "stress_exposure_panel_tests.csv"))
  testthat::skip_if(is.null(u) || is.null(pt), "test_stress_exposure has not been run here")
  testthat::expect_equal(nrow(u), 18)
  testthat::expect_true(all(u$n_con == 3 & u$n_stressed == 6))
  main <- pt[pt$acquisitions == "all", , drop = FALSE]
  testthat::expect_equal(as.integer(table(main$comparison)[c("stressed_vs_con", "sus_vs_res")]), c(20L, 20L))
  for (cmp in unique(main$comparison)) {
    k <- main$comparison == cmp
    testthat::expect_equal(main$p_adjust[k], stats::p.adjust(main$p_permutation[k], method = "BH"), tolerance = 1e-12)
  }
  testthat::expect_true(all(main$p_permutation >= main$p_floor - 1e-12, na.rm = TRUE))
})

testthat::test_that("unit specificity holds the declared contrasts and the Figure 3a CA2 SLM proteins", {
  s <- read_if(tab("test_unit_specificity", "unit_specificity_summary.csv"))
  c2 <- read_if(tab("test_unit_specificity", "unit_specificity_ca2_slm_proteins.csv"))
  testthat::skip_if(is.null(s) || is.null(c2), "test_unit_specificity has not been run here")
  testthat::expect_true(all(c("CA2slm vs CA1slm", "CA2slm vs CA2sr", "CA2slm vs other laminae", "CA3sr vs CA3so",
                              "CA2microglia vs other ROIs") %in% s$comparison))
  testthat::expect_equal(sum(s$kind == "per_unit"), 2 * 18)
  testthat::expect_true(all(c2$official_p_adjust < 0.05))
})

testthat::test_that("external concordance infers on the contrast-specific relabelling p", {
  sg <- read_if(tab("compare_external_stress_studies", "external_concordance_signatures.csv"))
  om <- read_if(tab("compare_external_stress_studies", "external_concordance_omnibus.csv"))
  testthat::skip_if(is.null(sg) || is.null(om), "compare_external_stress_studies has not been run here")
  for (x in list(sg, om)) {
    testthat::expect_false("p_all_1680" %in% names(x))
    testthat::expect_true(all(c("p_all", "p_contrast_specific", "n_relabellings", "n_contrast_specific") %in% names(x)))
    two <- x$contrast != "stressed_vs_CON"
    testthat::expect_true(all(x$n_contrast_specific[two] < x$n_relabellings[two]))
    testthat::expect_true(all(x$n_contrast_specific[!two] == x$n_relabellings[!two]))
  }
  # BH (Tobias, 2026-11-22): on the contrast-specific p, within level, primary set only
  pr <- sg[sg$set == "primary", ]
  for (lv in unique(pr$level)) {
    k <- pr$level == lv
    testthat::expect_equal(pr$q_bh[k], stats::p.adjust(pr$p_contrast_specific[k], method = "BH"), tolerance = 1e-12)
  }
  testthat::expect_true(all(is.na(sg$q_bh[sg$set != "primary"])))
})
