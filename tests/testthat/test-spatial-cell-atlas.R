source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "spatial_atlas_utils.R"))
source(testthat::test_path("..", "..", "R", "ewce_gene_set_engine.R"))

A <- function(...) path_results("tables", "11_spatial_systems", "atlas", ...)
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

# =====================================================================
# Q2 / Q3  spatial profile metrics
# =====================================================================

testthat::test_that("peak comes from raw values, never from the display z-scores", {
  # A profile whose RAW peak and whose z-scored peak would differ if the metric
  # were taken from standardised values. Row-standardising is monotone within a
  # row, so the peak must be identical - this pins that the code never takes the
  # peak from a rescaled matrix.
  x <- c(u1 = -0.3, u2 = 0.1, u3 = 0.9, u4 = -0.7)
  m <- sat_profile_metrics(x)
  testthat::expect_identical(m$peak_unit, "u3")
  testthat::expect_identical(m$second_unit, "u2")
  testthat::expect_equal(m$peak_value, 0.9)
  # peak-minus-rest and peak-minus-second are in RAW units, signed-safe
  testthat::expect_equal(m$peak_minus_second, 0.8)
  testthat::expect_equal(m$peak_minus_rest, 0.9 - mean(c(-0.3, 0.1, -0.7)))
  # the z matrix must not change the ordering
  z <- sat_row_z(matrix(x, nrow = 1, dimnames = list("m", names(x))))
  testthat::expect_identical(colnames(z)[which.max(z[1, ])], m$peak_unit)
})

testthat::test_that("specificity is bounded and handles the degenerate flat case", {
  n <- 10
  spike <- c(1, rep(0, n - 1))
  flat <- rep(0.5, n)
  s <- sat_profile_metrics(spike); f <- sat_profile_metrics(flat)
  # a single-unit profile is maximally specific
  testthat::expect_equal(s$spatial_tau, 1)
  testthat::expect_equal(s$spatial_entropy_normalized, 0)
  # a flat profile is maximally UNSPECIFIC, not NA: min-shifting a flat profile
  # gives an all-zero deviation vector, which would otherwise divide by zero
  testthat::expect_equal(f$spatial_tau, 0)
  testthat::expect_equal(f$spatial_entropy_normalized, 1)
  testthat::expect_false(is.na(f$spatial_tau))
  # bounded in [0,1] over random signed profiles
  set.seed(3)
  for (i in 1:200) {
    v <- stats::rnorm(sample(3:12, 1))
    mm <- sat_profile_metrics(v)
    testthat::expect_gte(mm$spatial_tau, 0); testthat::expect_lte(mm$spatial_tau, 1)
    testthat::expect_gte(mm$spatial_entropy_normalized, 0)
    testthat::expect_lte(mm$spatial_entropy_normalized, 1)
  }
})

testthat::test_that("entropy normalisation makes different unit counts comparable", {
  # neuropil resolves 10 units, soma and microglia resolve 4. Without dividing
  # by log(n) the 10-unit profiles would look systematically less specific.
  spike4 <- c(1, 0, 0, 0)
  spike10 <- c(1, rep(0, 9))
  testthat::expect_equal(sat_profile_metrics(spike4)$spatial_entropy_normalized,
                         sat_profile_metrics(spike10)$spatial_entropy_normalized)
  testthat::expect_equal(sat_profile_metrics(spike4)$spatial_tau,
                         sat_profile_metrics(spike10)$spatial_tau)
  flat4 <- rep(2, 4); flat10 <- rep(2, 10)
  testthat::expect_equal(sat_profile_metrics(flat4)$spatial_entropy_normalized,
                         sat_profile_metrics(flat10)$spatial_entropy_normalized)
  # and a graded profile lands strictly between the two extremes at both n
  g4 <- sat_profile_metrics(seq(-1, 1, length.out = 4))$spatial_entropy_normalized
  g10 <- sat_profile_metrics(seq(-1, 1, length.out = 10))$spatial_entropy_normalized
  for (g in c(g4, g10)) { testthat::expect_gt(g, 0); testthat::expect_lt(g, 1) }
})

testthat::test_that("tau equals the range-normalised peak elevation", {
  # the two formulations are algebraically the same quantity; pinning it stops
  # a future edit from reporting two different numbers under two names
  peak_elev <- function(x) {
    r <- max(x) - min(x); if (r <= 0) return(0)
    (max(x) - mean(x[-which.max(x)])) / r
  }
  set.seed(11)
  for (i in 1:200) {
    v <- stats::rnorm(sample(3:12, 1))
    testthat::expect_equal(sat_profile_metrics(v)$spatial_tau, peak_elev(v),
                           tolerance = 1e-10)
  }
})

# =====================================================================
# Q4 / Q5  enrichment universe and scope separation
# =====================================================================

testthat::test_that("enrichment restricts every set to the measured universe", {
  universe <- paste0("P", 1:100)
  mod <- c(paste0("P", 1:10), "NOT_MEASURED_1", "NOT_MEASURED_2")
  mk <- c(paste0("P", 5:20), "NOT_MEASURED_1")
  e <- sat_fisher_enrichment(mod, mk, universe)
  testthat::expect_equal(e$n_universe, 100L)
  # out-of-universe ids are dropped from BOTH sets, never counted as overlap
  testthat::expect_equal(e$n_module_in_universe, 10L)
  testthat::expect_equal(e$n_markers_in_universe, 16L)
  testthat::expect_equal(e$n_overlap, 6L)   # P5..P10
  testthat::expect_equal(e$expected_overlap, 10 * 16 / 100)
  testthat::expect_gt(e$odds_ratio, 1)
  # a set disjoint from the module gives no enrichment
  e0 <- sat_fisher_enrichment(paste0("P", 1:10), paste0("P", 50:60), universe)
  testthat::expect_equal(e0$n_overlap, 0L)
  testthat::expect_gt(e0$p_value, 0.5)
})

testthat::test_that("a marker set that spans most of the universe is not informative", {
  # empirical_microglia_neuropil_shared covers ~73% of the measured proteome and
  # would win any enrichment comparison on power alone; a 1-protein set cannot
  # support a stable test either.
  testthat::expect_false(sat_marker_set_is_informative(3671, 5054))  # 72.6%
  testthat::expect_false(sat_marker_set_is_informative(1, 5054))     # 1 protein
  testthat::expect_true(sat_marker_set_is_informative(1152, 5054))   # 22.8%
  testthat::expect_true(sat_marker_set_is_informative(88, 5054))     # 1.7%
  testthat::expect_true(sat_marker_set_is_informative(33, 5054))     # 0.65%
  # exactly at the bounds
  b <- sat_marker_set_bounds()
  testthat::expect_true(sat_marker_set_is_informative(b$max_fraction * 1000, 1000))
  testthat::expect_false(sat_marker_set_is_informative(b$max_fraction * 1000 + 1, 1000))
})

testthat::test_that("FDR families keep datasets and scopes apart", {
  fams <- c(sat_fdr_family_empirical_compartment("neuron_neuropil", "all"),
            sat_fdr_family_empirical_compartment("neuron_neuropil", "core_kME06"),
            sat_fdr_family_empirical_compartment("neuron_soma", "all"),
            sat_fdr_family_reference_marker("neuron_neuropil", "all"))
  testthat::expect_identical(anyDuplicated(fams), 0L)
  # empirical and reference families are never merged
  testthat::expect_false(identical(
    sat_fdr_family_empirical_compartment("d", "all"),
    sat_fdr_family_reference_marker("d", "all")))
  # BH is applied strictly within a family
  d <- data.frame(p_value = c(0.01, 0.02, 0.03, 0.04),
                  fdr_family = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  out <- sat_apply_family_fdr(d)
  testthat::expect_equal(out$FDR[1:2], stats::p.adjust(c(0.01, 0.02), "BH"))
  testthat::expect_equal(out$FDR[3:4], stats::p.adjust(c(0.03, 0.04), "BH"))
  testthat::expect_error(sat_apply_family_fdr(
    data.frame(p_value = 0.1, fdr_family = NA_character_)), "explicit FDR family")
})

# =====================================================================
# Q7 / Q8  bilateral reliability qualifies, never deletes
# =====================================================================

testthat::test_that("bilateral reliability does not alter the biological identity", {
  # identity metrics are a pure function of the profile; no bilateral input
  x <- c(a = 1, b = 0.2, c = -0.5, d = 0.1)
  testthat::expect_identical(names(formals(sat_profile_metrics)),
                             c("x", "unit_names"))
  m1 <- sat_profile_metrics(x)
  m2 <- sat_profile_metrics(x)
  testthat::expect_identical(m1, m2)
  testthat::expect_identical(m1$peak_unit, "a")
})

testthat::test_that("low bilateral support produces a caveat, not a deletion", {
  ctx <- sat_context_class(has_spatial = TRUE, has_compartment = FALSE,
                           has_external = FALSE,
                           bilateral_class = "low_bilateral_support")
  testthat::expect_identical(ctx$context_confidence, "bilaterally_variable")
  testthat::expect_match(ctx$context_rule, "identity retained with a caveat")
  # the same module with good bilateral support keeps its spatial class
  ok <- sat_context_class(TRUE, FALSE, FALSE, "high_bilateral_support")
  testthat::expect_identical(ok$context_confidence, "spatial_only")
  # thresholds are prespecified and symmetric
  th <- sat_bilateral_thresholds()
  testthat::expect_true(th$strong > th$moderate)
  testthat::expect_identical(sat_bilateral_support_class(0.85, NA), "high_bilateral_support")
  testthat::expect_identical(sat_bilateral_support_class(0.65, NA), "moderate_bilateral_support")
  testthat::expect_identical(sat_bilateral_support_class(0.20, NA), "low_bilateral_support")
  # the better of absolute and profile agreement is used, so a preserved spatial
  # pattern is not penalised by a constant side offset
  testthat::expect_identical(sat_bilateral_support_class(0.20, 0.9), "high_bilateral_support")
})

testthat::test_that("the context class always carries the rule that produced it", {
  ctx <- sat_context_class(c(TRUE, FALSE, FALSE), c(FALSE, TRUE, FALSE),
                           c(FALSE, FALSE, FALSE),
                           rep("high_bilateral_support", 3))
  testthat::expect_identical(nrow(ctx), 3L)
  testthat::expect_true(all(nzchar(ctx$context_rule)))
  testthat::expect_identical(ctx$context_confidence,
                             c("spatial_only", "cell_affinity_only", "weak_context"))
})

# =====================================================================
# canonical spatial unit vocabulary
# =====================================================================

testthat::test_that("spatial unit names normalise across producer conventions", {
  # four conventions are in play; an unnormalised join silently returns all-NA
  testthat::expect_identical(sat_canonical_spatial_unit("CA2_slm", "neuron_neuropil"), "ca2_slm")
  testthat::expect_identical(sat_canonical_spatial_unit("ca2_slm", "neuron_neuropil"), "ca2_slm")
  testthat::expect_identical(sat_canonical_spatial_unit("CA2_sp", "neuron_soma"), "ca2")
  testthat::expect_identical(sat_canonical_spatial_unit("DG_sg", "neuron_soma"), "dg")
  testthat::expect_identical(sat_canonical_spatial_unit("CA2_microglia", "microglia"), "ca2")
  testthat::expect_true(is.na(sat_canonical_spatial_unit("", "microglia")))
  # vectorised, and the neuropil layer token is preserved rather than stripped
  v <- sat_canonical_spatial_unit(c("CA1_slm", "DG_po"), "neuron_neuropil")
  testthat::expect_identical(v, c("ca1_slm", "dg_po"))
})

# =====================================================================
# generated outputs
# =====================================================================

testthat::test_that("the module atlas keeps its evidence dimensions separate", {
  p <- A("WGCNA_module_spatial_cell_affinity.csv")
  testthat::skip_if_not(file.exists(p), "module atlas not generated")
  x <- rd(p)
  for (col in c("spatial_tau", "strongest_empirical_compartment",
                "external_celltype_all", "bilateral_support_class",
                "context_confidence", "context_rule")) {
    testthat::expect_true(col %in% names(x), info = col)
  }
  # no composite numeric score masquerading as confidence
  testthat::expect_true(is.character(x$context_confidence))
  testthat::expect_true(all(nzchar(x$context_rule)))
  # module ids are canonical and phenotype never leaks in
  testthat::expect_true(all(grepl("^WGCNA_m[0-9]+$", x$ModuleID)))
  testthat::expect_false(any(c("StressGroup", "SUS", "RES") %in% names(x)))
})

testthat::test_that("all three module scopes survive into the atlas long tables", {
  p <- A("WGCNA_module_empirical_compartment_affinity_long.csv")
  testthat::skip_if_not(file.exists(p), "compartment affinity not generated")
  x <- rd(p)
  testthat::expect_setequal(unique(x$module_scope), c("all", "core_kME06", "top25"))
  # each scope has its own FDR family, so scopes cannot contaminate each other
  fam <- unique(x[, c("dataset", "module_scope", "fdr_family")])
  testthat::expect_identical(anyDuplicated(fam$fdr_family), 0L)
  testthat::expect_true(all(x$n_universe > 4000))
  testthat::expect_true("informative_for_affinity_call" %in% names(x))
})

testthat::test_that("the protein atlas is the deduplicated union and keeps all 37", {
  p <- A("protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(file.exists(p), "protein atlas not generated")
  x <- rd(p)
  testthat::expect_identical(anyDuplicated(paste(x$dataset, x$ProteinGroupID)), 0L)
  testthat::expect_equal(sum(x$is_sus_res_fdr_supported %in% TRUE), 37L)
  # the full candidate atlas is NOT narrowed to the 37
  testthat::expect_gt(nrow(x), 37L)
  testthat::expect_true(all(c("is_sus_res_fdr_supported", "is_wgcna_candidate",
                              "candidate_tier") %in% names(x)))
  s <- A("protein_sus_res_fdr_supported_atlas.csv")
  if (file.exists(s)) {
    y <- rd(s)
    testthat::expect_equal(nrow(y), 37L)
    testthat::expect_identical(anyDuplicated(paste(y$dataset, y$ProteinGroupID)), 0L)
  }
})

testthat::test_that("DA statistics are joined with provenance, not recomputed", {
  p <- A("protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(file.exists(p))
  x <- rd(p)
  testthat::expect_true("da_provenance" %in% names(x))
  testthat::expect_true(all(grepl("no differential statistic recomputed",
                                  x$da_provenance)))
  # the producer must not contain a differential model fit
  src <- paste(readLines(repo_path("11_spatial_systems",
                                   "09_protein_spatial_cell_atlas.R"),
                         warn = FALSE), collapse = "\n")
  for (bad in c("lmFit", "eBayes", "p.adjust", "t.test", "wilcox.test")) {
    testthat::expect_false(grepl(bad, src, fixed = TRUE), info = bad)
  }
})

testthat::test_that("effect-at-identity classification is deterministic and prespecified", {
  p <- A("protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(file.exists(p))
  x <- rd(p)
  f <- x[x$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
  testthat::skip_if(nrow(f) == 0L)
  testthat::expect_true(all(f$effect_identity_relationship %in% c(
    "effect_at_baseline_peak", "effect_in_top2_baseline_units",
    "effect_in_high_affinity_unit", "effect_outside_baseline_affinity",
    "spatially_broad", "insufficient_context")))
  # the high-affinity rule is recorded and is within-protein, not per-protein tuned
  testthat::expect_true(all(grepl("prespecified, not tuned", f$high_affinity_rule)))
  testthat::expect_identical(length(unique(f$high_affinity_rule)), 1L)
  # the cutoff follows only from the unit count, never from the data
  testthat::expect_true(all(
    f$baseline_high_affinity_cutoff_rank == ceiling(f$n_spatial_units / 3),
    na.rm = TRUE))
})

testthat::test_that("the CA2-SLM audit uses all tested proteins as denominator", {
  p <- A("neuropil_spatial_detection_context.csv")
  testthat::skip_if_not(file.exists(p), "detection context not generated")
  x <- rd(p)
  # the denominator must be the full tested set, not the candidate subset
  testthat::expect_true(all(x$n_proteins_tested > 5000))
  testthat::expect_identical(length(unique(x$n_proteins_tested)), 1L)
  testthat::expect_true(all(grepl("not the candidate subset", x$denominator_definition)))
  # every neuropil unit is represented: joining on the wrong column drops units
  testthat::expect_equal(nrow(x), 10L)
  testthat::expect_true(all(c("median_baseline_abundance", "median_between_animal_sd",
                              "median_abs_log2FC") %in% names(x)))
  testthat::expect_false(any(is.na(x$median_baseline_abundance)))
})

testthat::test_that("labels and WGCNA membership are unchanged by the atlas", {
  p <- A("WGCNA_label_spatial_cell_context_audit.csv")
  testthat::skip_if_not(file.exists(p), "label audit not generated")
  x <- rd(p)
  testthat::expect_true(all(grepl("AUDIT ONLY", x$audit_note)))
  testthat::expect_true(all(x$label_context_relationship %in% c(
    "context_corroborates_label", "context_compatible", "context_adds_specificity",
    "context_mixed", "potential_context_conflict", "insufficient_context")))
  # no atlas producer may write into the label registry or WGCNA state
  for (f in c("08_module_spatial_cell_atlas.R", "09_protein_spatial_cell_atlas.R",
              "10_neuropil_spatial_detection_context.R")) {
    src <- readLines(repo_path("11_spatial_systems", f), warn = FALSE)
    writes <- grep("write_csv_safe|write\\.csv|saveRDS", src, value = TRUE)
    testthat::expect_false(any(grepl("config|wgcna_labels", writes)), info = f)
    testthat::expect_false(any(grepl("01_WGCNA|modules/WGCNA_modules_long", writes)),
                           info = f)
  }
})

testthat::test_that("the atlas validation has no critical failure", {
  p <- path_results("tables", "11_spatial_systems", "spatial_systems_atlas_validation.csv")
  testthat::skip_if_not(file.exists(p), "atlas validation not generated")
  v <- rd(p)
  testthat::expect_gt(nrow(v), 10L)
  testthat::expect_equal(sum(v$critical %in% TRUE & v$status == "FAIL"), 0L)
})

testthat::test_that("every candidate figure has matching source data", {
  d <- path_results("figures", "11_spatial_systems", "atlas")
  testthat::skip_if_not(dir.exists(d), "atlas figures not generated")
  pngs <- list.files(d, pattern = "[.]png$")
  testthat::skip_if(length(pngs) == 0L)
  for (f in pngs) {
    src <- file.path(d, sub("[.]png$", "_source_data.csv", f))
    testthat::expect_true(file.exists(src), info = f)
  }
})
