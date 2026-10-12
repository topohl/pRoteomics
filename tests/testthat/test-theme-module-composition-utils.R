testthat::local_edition(3)

# Synthetic fixtures only: these tests never source the runner and never read results/.
source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "theme_module_composition_utils.R"))

U <- sprintf("PG%03d", 1:100)

testthat::test_that("the 2x2 table, shares, fold and odds ratio are the textbook ones", {
  r <- tmc_overlap_row(theme_set = U[1:20], entity_set = U[11:40], universe = U)
  testthat::expect_identical(c(r$n_universe, r$n_theme, r$n_entity, r$n_overlap), c(100L, 20L, 30L, 10L))
  testthat::expect_equal(r$frac_theme_in_entity, 10 / 20)
  testthat::expect_equal(r$frac_entity_in_theme, 10 / 30)
  testthat::expect_equal(r$fold_enrichment, (10 / 20) / (30 / 100))
  testthat::expect_equal(r$jaccard, 10 / 40)
  m <- matrix(c(10, 10, 20, 60), 2)
  testthat::expect_equal(r$fisher_p_greater, stats::fisher.test(m, alternative = "greater")$p.value)
  testthat::expect_equal(r$odds_ratio, unname(stats::fisher.test(m)$estimate))
  testthat::expect_true(r$odds_ratio_ci_low < r$odds_ratio && r$odds_ratio < r$odds_ratio_ci_high)
})

testthat::test_that("members outside the universe are ignored, not counted", {
  r <- tmc_overlap_row(theme_set = c(U[1:5], "OUTSIDE1"), entity_set = c(U[1:3], "OUTSIDE1"), universe = U)
  testthat::expect_identical(c(r$n_theme, r$n_entity, r$n_overlap), c(5L, 3L, 3L))
})

testthat::test_that("support needs BH-FDR <= 0.05 and an odds ratio of at least 2, within the family", {
  themes <- list(A = U[1:20], B = U[41:60])
  ents <- list(m01 = U[1:20], m02 = U[21:40], m03 = U[41:100])
  comp <- tmc_composition(themes, ents, U, "ds", "module", "go_membership", "wgcna_ranked_annotated")
  testthat::expect_identical(nrow(comp), 6L)
  testthat::expect_true(all(comp$family_size == 6L))
  testthat::expect_equal(comp$BH_FDR, stats::p.adjust(comp$fisher_p_greater, "BH"))
  get <- function(t, e) comp[comp$theme_id == t & comp$entity_id == e, ]
  testthat::expect_true(get("A", "m01")$composition_supported)        # the theme is the module
  testthat::expect_false(get("A", "m02")$composition_supported)       # disjoint
  # B sits entirely inside the giant m03 (60% of the universe): enriched, but the OR floor decides
  testthat::expect_identical(get("B", "m03")$n_overlap, 20L)
  testthat::expect_identical(get("B", "m03")$composition_supported,
                             get("B", "m03")$BH_FDR <= 0.05 && get("B", "m03")$odds_ratio >= 2)
  testthat::expect_true(all(comp$evidence_role == "structural_composition_not_phenotype"))
})

testthat::test_that("a large module that merely contains a theme at its base rate is not support", {
  # theme of 40 proteins, 24 of them in a 60-protein module: exactly the base rate (OR about 1)
  themes <- list(T = c(U[1:24], U[61:76]))
  ents <- list(big = U[1:60], rest = U[61:100])
  comp <- tmc_composition(themes, ents, U, "ds", "module", "go_membership", "u")
  big <- comp[comp$entity_id == "big", ]
  testthat::expect_equal(big$fold_enrichment, 1)
  testthat::expect_false(big$composition_supported)
})

testthat::test_that("single-theme sets drop proteins shared between themes", {
  s <- tmc_single_theme_sets(list(A = U[1:10], B = U[6:15], C = U[20:25]))
  testthat::expect_identical(s$A, U[1:5])
  testthat::expect_identical(s$B, U[11:15])
  testthat::expect_identical(s$C, U[20:25])
})

testthat::test_that("classes follow the share captured by supporting modules", {
  testthat::expect_identical(tmc_class(c(0.6, 0.3, 0.1, 0.9, NA), c(1L, 2L, 1L, 0L, 1L)),
                             c("coherent", "core", "minor", "scattered", "scattered"))
})

testthat::test_that("the coherence summary sums carrier shares and flags shared carriers", {
  themes <- list(A = U[1:20], B = U[11:30], C = U[81:90])
  ents <- list(m01 = U[1:30], m02 = U[31:80], m03 = U[81:100])
  comp <- tmc_composition(themes, ents, U, "ds", "module", "go_membership", "u")
  s <- tmc_coherence_summary(comp)
  a <- s[s$theme_id == "A", ]
  testthat::expect_identical(a$carriers, "m01")
  testthat::expect_equal(a$share_captured_by_carriers, 1)
  testthat::expect_identical(a$composition_class, "coherent")
  testthat::expect_identical(a$shares_carrier_with, "B")              # A and B share m01
  testthat::expect_identical(s$shares_carrier_with[s$theme_id == "C"], "")
})

testthat::test_that("partition agreement is perfect for identical partitions and seeded", {
  lab <- rep(c("x", "y", "z"), each = 30)
  ent <- rep(c("m1", "m2", "m3"), each = 30)
  a <- tmc_partition_agreement(lab, ent, n_perm = 200L, seed = 1L)
  testthat::expect_equal(a$cramers_v, 1)
  testthat::expect_equal(a$nmi, 1)
  testthat::expect_lt(a$cramers_v_perm_p, 0.01)
  b <- tmc_partition_agreement(lab, ent, n_perm = 200L, seed = 1L)
  testthat::expect_identical(a, b)
})

testthat::test_that("the curated mapping audit classifies rows and lists unmapped support", {
  themes <- list(A = U[1:20], B = U[41:60])
  ents <- list(m01 = U[1:20], m02 = U[21:40], m03 = U[41:60])
  comp <- tmc_composition(themes, ents, U, "ds", "module", "go_membership", "u")
  mapping <- data.frame(dataset = "ds", theme_id = c("A", "A"), entity_level = "module",
                        entity_id = c("m01", "m02"), approved_for_manuscript_interpretation = c("TRUE", "TRUE"),
                        mapping_role = "primary", confidence = "high", stringsAsFactors = FALSE)
  au <- tmc_mapping_audit(mapping, comp)
  testthat::expect_identical(au$curated_mapping_composition_status[au$entity_id == "m01" & au$theme_id == "A"], "supported")
  testthat::expect_identical(au$curated_mapping_composition_status[au$entity_id == "m02" & au$theme_id == "A"], "not_enriched")
  testthat::expect_identical(au$curated_mapping_composition_status[au$entity_id == "m03"], "supported_but_unmapped")
})

testthat::test_that("the wheel table lists each module once with the themes it carries", {
  themes <- list(A = U[1:20], B = U[11:30])
  ents <- list(m01 = U[1:30], m02 = U[31:100])
  comp <- tmc_composition(themes, ents, U, "ds", "module", "go_membership", "u")
  membership <- data.frame(dataset = "ds", module_id = c("m01", "m02"), supermodule_id = c("SM01", "SM01"),
                           stringsAsFactors = FALSE)
  sizes <- data.frame(dataset = "ds", module_id = c("m01", "m02"), n_proteins = c(30L, 70L))
  spatial <- data.frame(dataset = "ds", module_id = c("m01", "m02"), peak_unit = c("ca1", "dg"))
  w <- tmc_wheel_structure(membership, sizes, spatial, comp, theme_order = c("A", "B"))
  testthat::expect_identical(nrow(w), 2L)
  testthat::expect_true(all(w$is_multi_module_supermodule))
  testthat::expect_identical(w$carried_themes[w$module_id == "m01"], "A;B")
  testthat::expect_identical(w$top_theme[w$module_id == "m01"], "A")
  testthat::expect_identical(w$carried_themes[w$module_id == "m02"], "")
})
