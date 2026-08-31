clusterprofiler_repro_helper <- testthat::test_path(
  "..", "..", "R", "clusterprofiler_reproducibility.R"
)
source(clusterprofiler_repro_helper)

testthat::test_that("per-comparison GSEA seed derivation is stable and identity-specific", {
  base_seed <- 20260824L
  comparison <- "CA2slmsus_CA2slmcon"
  go_type <- clusterprofiler_gsea_analysis_types("BP")[["go"]]

  first <- derive_clusterprofiler_gsea_seed(base_seed, comparison, go_type)
  second <- derive_clusterprofiler_gsea_seed(base_seed, comparison, go_type)
  different_comparison <- derive_clusterprofiler_gsea_seed(
    base_seed, "CA2slmres_CA2slmcon", go_type
  )
  different_type <- derive_clusterprofiler_gsea_seed(
    base_seed, comparison, clusterprofiler_gsea_analysis_types("BP")[["kegg"]]
  )

  testthat::expect_identical(first, second)
  testthat::expect_identical(first, 288057122L)
  testthat::expect_false(identical(first, different_comparison))
  testthat::expect_false(identical(first, different_type))
  testthat::expect_true(all(c(first, different_comparison, different_type) > 0L))
})

testthat::test_that("GSEA reproducibility config is parsed as validated integers", {
  committed <- yaml::read_yaml(testthat::test_path(
    "..", "..", "config", "clusterProfiler_config.yml"
  ))
  committed <- validate_clusterprofiler_gsea_config(committed)
  testthat::expect_identical(committed$analysis$gsea_seed_base, 20260824L)
  testthat::expect_identical(committed$analysis$n_perm_simple, 100000L)

  good <- validate_clusterprofiler_gsea_config(list(analysis = list(
    gsea_seed_base = 20260824,
    n_perm_simple = 10000
  )))
  testthat::expect_identical(good$analysis$gsea_seed_base, 20260824L)
  testthat::expect_identical(good$analysis$n_perm_simple, 10000L)

  testthat::expect_error(
    validate_clusterprofiler_gsea_config(list(analysis = list(
      gsea_seed_base = "20260824", n_perm_simple = 10000
    ))),
    "analysis.gsea_seed_base"
  )
  testthat::expect_error(
    validate_clusterprofiler_gsea_config(list(analysis = list(
      gsea_seed_base = 20260824, n_perm_simple = 0
    ))),
    "analysis.n_perm_simple"
  )
  testthat::expect_error(
    validate_clusterprofiler_gsea_config(list(analysis = list(
      gsea_seed_base = 20260824, n_perm_simple = 10000.5
    ))),
    "whole number"
  )
})

testthat::test_that("seeded GSEA wrapper supplies fgsea controls and restores RNG state", {
  capture_call <- function(...) {
    list(arguments = list(...), random_draw = stats::runif(4))
  }

  set.seed(4242)
  rng_kind_before <- RNGkind()
  rng_state_before <- .Random.seed
  first <- run_seeded_clusterprofiler_gsea(
    capture_call,
    gsea_seed = 288057122L,
    n_perm_simple = 10000L,
    geneList = c(A = 2, B = -1)
  )
  testthat::expect_identical(RNGkind(), rng_kind_before)
  testthat::expect_identical(.Random.seed, rng_state_before)

  second <- run_seeded_clusterprofiler_gsea(
    capture_call,
    gsea_seed = 288057122L,
    n_perm_simple = 10000L,
    geneList = c(A = 2, B = -1)
  )
  testthat::expect_identical(first$random_draw, second$random_draw)
  testthat::expect_identical(first$arguments$by, "fgsea")
  testthat::expect_identical(first$arguments$seed, FALSE)
  testthat::expect_identical(first$arguments$nPermSimple, 10000L)
  testthat::expect_identical(first$arguments$geneList, c(A = 2, B = -1))
})

testthat::test_that("GSEA reproducibility audit round-trips and detects stale settings", {
  path <- tempfile("gsea-reproducibility-", fileext = ".csv")
  audit <- clusterprofiler_gsea_reproducibility_table(
    comparison = "CA2slmsus_CA2slmcon",
    ontology = "BP",
    gsea_seed_base = 20260824L,
    n_perm_simple = 10000L,
    requested = c(go = TRUE, kegg = TRUE, kegg_predefined = FALSE, nk3r = FALSE)
  )
  utils::write.csv(audit, path, row.names = FALSE)

  testthat::expect_true(clusterprofiler_reproducibility_audit_matches(
    path, "CA2slmsus_CA2slmcon", "BP", 20260824L, 10000L
  ))
  testthat::expect_false(clusterprofiler_reproducibility_audit_matches(
    path, "CA2slmsus_CA2slmcon", "BP", 20260825L, 10000L
  ))
  testthat::expect_false(clusterprofiler_reproducibility_audit_matches(
    path, "CA2slmsus_CA2slmcon", "BP", 20260824L, 20000L
  ))
})

testthat::test_that("single-comparison filter supports CLI and environment targeting", {
  names <- c("CA1slmsus_CA1slmcon", "CA2slmsus_CA2slmcon")
  cli <- resolve_clusterprofiler_comparison_filter(
    args = c("--dataset", "neuron_neuropil", "--comparison", "CA2slmsus_CA2slmcon.csv"),
    env_value = "CA1slmsus_CA1slmcon"
  )
  env <- resolve_clusterprofiler_comparison_filter(
    args = character(), env_value = "CA1slmsus_CA1slmcon"
  )

  testthat::expect_identical(cli, "CA2slmsus_CA2slmcon")
  testthat::expect_identical(env, "CA1slmsus_CA1slmcon")
  testthat::expect_identical(filter_clusterprofiler_comparisons(names, cli), 2L)
  testthat::expect_identical(filter_clusterprofiler_comparisons(names, ""), 1:2)
  testthat::expect_error(
    filter_clusterprofiler_comparisons(names, "missing_comparison"),
    "not found exactly once"
  )
  testthat::expect_error(
    resolve_clusterprofiler_comparison_filter(args = "--comparison", env_value = ""),
    "requires an exact comparison name"
  )
})

testthat::test_that("active GSEA entry points all use the deterministic wrapper", {
  script <- paste(readLines(
    testthat::test_path("..", "..", "04_differential_expression_enrichment", "01_clusterProfiler.r"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_equal(
    lengths(regmatches(script, gregexpr("run_gsea\\(", script, perl = TRUE))),
    4L
  )
  testthat::expect_false(grepl("<-\\s*(gseGO|gseKEGG)\\(", script, perl = TRUE))
  testthat::expect_false(grepl("<-\\s*clusterProfiler::GSEA\\(", script, perl = TRUE))
  testthat::expect_match(script, "n_perm_simple = analysis_params\\$n_perm_simple")
})

testthat::test_that("the GSEA parallel-backend contract is deliberate and documented", {
  src <- readLines(
    testthat::test_path("..", "..", "R", "clusterprofiler_reproducibility.R"),
    warn = FALSE
  )
  # The wrapper intentionally does NOT pin BPPARAM: backend, worker count and
  # bpRNGseed were measured not to affect results (SnowParam 30 unseeded vs
  # SnowParam 8/4 with a pinned RNGseed were field-identical). If someone adds
  # backend pinning later they must also revisit this contract and the
  # documented tolerance.
  control_args <- grep("clusterprofiler_fgsea_control_args", src)
  testthat::expect_gt(length(control_args), 0L)
  testthat::expect_false(any(grepl("BPPARAM", src, fixed = TRUE) &
                               !grepl("^#", trimws(src))))

  # The measured contract must stay stated, so the language cannot silently
  # drift back to promising bit identity.
  joined <- paste(src, collapse = "\n")
  testthat::expect_match(joined, "Reproducibility contract", fixed = TRUE)
  testthat::expect_match(joined, "does NOT guarantee", fixed = TRUE)
  testthat::expect_match(joined, "not as floating-point bit identity", fixed = TRUE)
})

testthat::test_that("reserved fgsea controls still cannot be overridden via dots", {
  for (arg in c("by", "seed", "nPermSimple")) {
    args <- list(function(...) NULL, gsea_seed = 1L, n_perm_simple = 10L)
    args[[arg]] <- "x"
    testthat::expect_error(
      do.call(run_seeded_clusterprofiler_gsea, args),
      "must be supplied by run_seeded_clusterprofiler_gsea", info = arg
    )
  }
})
