source(testthat::test_path("..", "..", "R", "ewce_contract_utils.R"))

testthat::test_that("canonical EWCE defaults to animal-level biological units", {
  default_contract <- ewce_resolve_run_contract("neuron_neuropil")
  explicit_animal <- ewce_resolve_run_contract("neuron_soma", analysis_unit = "animal")

  testthat::expect_identical(default_contract$analysis_unit, "animal")
  testthat::expect_identical(default_contract$branch, "canonical")
  testthat::expect_true(default_contract$canonical)
  testthat::expect_identical(default_contract$substep_id, file.path("EWCE_E9", "neuron_neuropil"))
  testthat::expect_true(explicit_animal$canonical)
  testthat::expect_identical(explicit_animal$substep_id, file.path("EWCE_E9", "neuron_soma"))
})

testthat::test_that("explicit branch isolates animal and legacy sample EWCE", {
  animal <- ewce_resolve_run_contract("microglia", analysis_unit = "animal", branch = "test_branch")
  sample <- ewce_resolve_run_contract("microglia", analysis_unit = "sample", branch = "sample_level_sensitivity")

  testthat::expect_false(animal$canonical)
  testthat::expect_identical(animal$substep_id, file.path("EWCE_E9_comparison", "test_branch", "microglia"))
  testthat::expect_false(sample$canonical)
  testthat::expect_identical(sample$substep_id, file.path("EWCE_E9_comparison", "sample_level_sensitivity", "microglia"))
})

testthat::test_that("unbranched sample-level EWCE fails closed", {
  testthat::expect_error(
    ewce_resolve_run_contract("neuron_neuropil", analysis_unit = "sample"),
    "Sample-level EWCE is legacy/sensitivity-only and requires PROTEOMICS_EWCE_BRANCH",
    fixed = TRUE
  )
})

testthat::test_that("animal cache identity is branch-independent and complete", {
  canonical <- ewce_resolve_run_contract("microglia", branch = "")
  comparison <- ewce_resolve_run_contract("microglia", branch = "test_branch")
  cache_args <- list(
    contract_version = "EWCE_animal_level_v1",
    input_matrix_sha256 = "input_sha256",
    dataset = "microglia",
    target = "CA1_Sus_vs_Con_up_top250",
    top_n = 250L,
    annot_level = 2L,
    reps = 10000L,
    seed = 42L,
    hits = c("A", "B"),
    background = c("A", "B", "C"),
    ctd_annotation_genes = c("A", "B", "C", "D")
  )

  canonical_key <- do.call(ewce_animal_cache_key, cache_args)
  comparison_key <- do.call(ewce_animal_cache_key, cache_args)
  testthat::expect_identical(canonical$analysis_unit, comparison$analysis_unit)
  testthat::expect_identical(canonical_key, comparison_key)
  testthat::expect_false("branch" %in% names(formals(ewce_animal_cache_key)))
  testthat::expect_false(ewce_legacy_cache_fallback_allowed("animal"))
  testthat::expect_true(ewce_legacy_cache_fallback_allowed("sample"))
})

testthat::test_that("cache accounting distinguishes reuse, fallback, and computation", {
  accounting <- ewce_cache_accounting_table(
    c("run_a", "run_b", "run_c"),
    data.frame(
      TargetRun = c("run_a", "run_b", "run_c"),
      CacheEvent = c("cache_hit", "legacy_cache_fallback", "computed"),
      stringsAsFactors = FALSE
    )
  )

  testthat::expect_equal(accounting$expected_target_runs, 3L)
  testthat::expect_equal(accounting$cache_hits, 1L)
  testthat::expect_equal(accounting$cache_misses, 2L)
  testthat::expect_equal(accounting$new_bootstrap_computations, 1L)
  testthat::expect_equal(accounting$legacy_cache_fallback_count, 1L)
})

testthat::test_that("EWCE never silently falls back to a background-free test", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  src <- readLines(repo_path("05_celltype_enrichment_EWCE", "01_EWCE_E9.r"),
                   warn = FALSE)
  txt <- paste(src, collapse = "\n")
  # Scan CODE only: the comment explaining this rule necessarily names the
  # very call signature the rule forbids.
  code <- paste(sub("^\\s*#.*$", "", src), collapse = "\n")

  # The retry that dropped `bg` substituted the full reference transcriptome for
  # the measured-proteome background while the caller still stamped
  # N_Background = length(bg), making a wrong result indistinguishable from a
  # correct one - including in the on-disk cache.
  calls <- unlist(regmatches(
    code, gregexpr("EWCE::bootstrap_enrichment_test[(](?:[^()]|[(][^()]*[)])*[)]", code)))
  testthat::expect_gt(length(calls), 0L)
  for (call in calls) {
    testthat::expect_true(grepl("bg[[:space:]]*=", call),
                          info = "an EWCE call omits the background argument")
  }
  # and the failure path must be a hard stop that explains itself
  testthat::expect_true(grepl("Refusing to retry without a background", txt, fixed = TRUE))
})
