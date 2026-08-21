source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "protigy_stat_gct_utils.R"))

testthat::test_that("historical and corrected statistical fields parse to one canonical comparison", {
  historical <- parse_protigy_stat_field("logFC.CA1_slm_3.over.CA1_slm_2")
  corrected <- parse_protigy_stat_field("logFC.CA1_slm_3_over_CA1_slm_2")

  testthat::expect_identical(historical$metric, "logFC")
  testthat::expect_identical(corrected$metric, "logFC")
  testthat::expect_identical(historical$comparison, "CA1_slm_3.over.CA1_slm_2")
  testthat::expect_identical(corrected$comparison, historical$comparison)
  testthat::expect_identical(historical$naming_style, "historical_dot_over")
  testthat::expect_identical(corrected$naming_style, "corrected_underscore_over")
})

testthat::test_that("metrics containing periods are parsed without splitting the metric", {
  fields <- c(
    "P.Value.CA1_slm_3_over_CA1_slm_2",
    "adj.P.Val.CA1_slm_3_over_CA1_slm_2",
    "Log.P.Value.CA1_slm_3_over_CA1_slm_2"
  )
  parsed <- parse_protigy_stat_field(fields)
  testthat::expect_identical(parsed$metric, c("P.Value", "adj.P.Val", "Log.P.Value"))
  testthat::expect_true(length(unique(parsed$comparison)) == 1L)
})

testthat::test_that("same-unit canonical spatial comparisons validate for each dataset contract", {
  neuropil <- parse_protigy_comparison("CA1_slm_3_over_CA1_slm_2", "neuron_neuropil")
  soma_ca1 <- parse_protigy_comparison("CA1_sp_2_over_CA1_sp_1", "neuron_soma")
  soma_dg <- parse_protigy_comparison("DG_sg_3_over_DG_sg_1", "neuron_soma")
  microglia <- parse_protigy_comparison("CA2_microglia_3_over_CA2_microglia_2", "microglia")

  testthat::expect_true(neuropil$valid_stress_comparison)
  testthat::expect_identical(neuropil$spatial_unit, "CA1_slm")
  testthat::expect_true(soma_ca1$valid_stress_comparison)
  testthat::expect_identical(soma_ca1$spatial_unit, "CA1_sp")
  testthat::expect_true(soma_dg$valid_stress_comparison)
  testthat::expect_identical(soma_dg$spatial_unit, "DG_sg")
  testthat::expect_true(microglia$valid_stress_comparison)
  testthat::expect_identical(microglia$spatial_unit, "CA2")
})

testthat::test_that("strict corrected comparison validation fails closed on cross-unit contrasts", {
  bad <- "CA1_slm_2_over_CA2_slm_2"
  parsed <- parse_protigy_comparison(bad, "neuron_neuropil")
  testthat::expect_false(parsed$same_unit)
  testthat::expect_false(parsed$valid_stress_comparison)
  testthat::expect_identical(parsed$rejection_reason, "cross_spatial_unit")
  testthat::expect_error(
    validate_protigy_comparisons(bad, "neuron_neuropil", strict_primary = TRUE),
    "rejected comparison"
  )
})

testthat::test_that("only the three primary directed stress pairs are accepted", {
  accepted <- c(
    "CA1_slm_2_over_CA1_slm_1",
    "CA1_slm_3_over_CA1_slm_2",
    "CA1_slm_3_over_CA1_slm_1"
  )
  validation <- validate_protigy_comparisons(accepted, "neuron_neuropil")
  testthat::expect_true(all(validation$valid_stress_comparison))
  testthat::expect_identical(
    validation$biological_contrast,
    c("RES_vs_CON", "SUS_vs_RES", "SUS_vs_CON")
  )
  reverse <- parse_protigy_comparison("CA1_slm_1_over_CA1_slm_2", "neuron_neuropil")
  testthat::expect_false(reverse$valid_stress_comparison)
  testthat::expect_identical(reverse$rejection_reason, "unsupported_group_pair")
})

testthat::test_that("reverse contrast semantics negate only verified signed statistics", {
  values <- c(1, -2)
  testthat::expect_equal(reverse_protigy_metric(values, "logFC"), c(-1, 2))
  testthat::expect_equal(reverse_protigy_metric(values, "RawlogFC"), c(-1, 2))
  testthat::expect_equal(reverse_protigy_metric(values, "t"), c(-1, 2))

  for (metric in c("P.Value", "adj.P.Val", "B", "AveExpr", "RawAveExpr", "Log.P.Value", "significant")) {
    testthat::expect_identical(reverse_protigy_metric(values, metric), values)
  }
})

testthat::test_that("historical and corrected comparison names match directly", {
  old <- canonicalize_protigy_comparison("CA1_slm_3.over.CA1_slm_2")
  new <- canonicalize_protigy_comparison("CA1_slm_3_over_CA1_slm_2")
  testthat::expect_identical(old, new)
})

testthat::test_that("protein matching is exact-ID only and reports both unmatched sides", {
  matched <- match_protigy_proteins(c("P2", "P1", "LEGACY"), c("P3", "P1", "P2"))
  testthat::expect_identical(matched$matched, c("P1", "P2"))
  testthat::expect_identical(matched$legacy_only, "LEGACY")
  testthat::expect_identical(matched$animal_only, "P3")
})

testthat::test_that("signed t ranking is deterministic with protein-ID tie breaking", {
  first <- deterministic_signed_rank(c(1, 1, 2), c("B", "A", "C"))
  second_order <- c(3, 1, 2)
  second <- deterministic_signed_rank(c(1, 1, 2)[second_order], c("B", "A", "C")[second_order])
  names(first) <- c("B", "A", "C")
  names(second) <- c("B", "A", "C")[second_order]
  testthat::expect_identical(first, c(B = 3L, A = 2L, C = 1L))
  testthat::expect_identical(first[sort(names(first))], second[sort(names(second))])
})

testthat::test_that("GCT root overrides resolve independently and canonical default is animal-level", {
  withr::local_envvar(c(
    PROTEOMICS_GCT_INPUT_ROOT = "custom/input",
    PROTEOMICS_GCT_OUTPUT_ROOT = "custom/output"
  ))
  overridden <- resolve_gct_io_roots()
  testthat::expect_identical(overridden$input_root, normalizePath(repo_path("custom", "input"), winslash = "/", mustWork = FALSE))
  testthat::expect_identical(overridden$output_root, normalizePath(repo_path("custom", "output"), winslash = "/", mustWork = FALSE))

  defaults <- resolve_gct_io_roots(input_root = "", output_root = "")
  testthat::expect_identical(
    defaults$input_root,
    normalizePath(path_processed("01_preprocessing", "protigy_output_animal_level"), winslash = "/", mustWork = FALSE)
  )
  legacy <- resolve_gct_io_roots(input_root = "", output_root = "", legacy_mode = TRUE)
  testthat::expect_identical(legacy$input_root, normalizePath(path_processed("01_preprocessing", "protigy_output"), winslash = "/", mustWork = FALSE))
  testthat::expect_identical(legacy$output_root, normalizePath(path_processed("01_preprocessing", "gct_extractR_legacy"), winslash = "/", mustWork = FALSE))
  testthat::expect_false(identical(legacy$output_root, defaults$output_root))
  testthat::expect_identical(
    defaults$output_root,
    normalizePath(path_processed("01_preprocessing", "gct_extractR"), winslash = "/", mustWork = FALSE)
  )
})

testthat::test_that("single-GCT resolution honors an explicit analysis root", {
  root <- tempfile("protigy-root-")
  dataset_dir <- file.path(root, "neuron_soma")
  dir.create(dataset_dir, recursive = TRUE)
  fixture <- file.path(dataset_dir, "fixture.gct")
  writeLines("#1.3", fixture)
  resolved <- resolve_single_protigy_gct(root, "neuron_soma", use_manifest = FALSE)
  testthat::expect_identical(basename(resolved$path), "fixture.gct")
  testthat::expect_identical(resolved$dataset_dir, normalizePath(dataset_dir, winslash = "/", mustWork = FALSE))
})

testthat::test_that("statistical row descriptors are read independently of ncmat", {
  comparison <- "CA1_sp_2_over_CA1_sp_1"
  fields <- paste0(protigy_required_da_metrics(), ".", comparison)
  header <- c("id", "Description", fields)
  row1 <- c("P1", "Protein one", "1", "2", "3", "0.01", "0.02", "4", "TRUE", "5")
  row2 <- c("P2", "Protein two", "-1", "2", "-3", "0.03", "0.04", "1", "TRUE", "4")
  fixture <- tempfile(fileext = ".gct")
  writeLines(c(
    "#1.3",
    paste(c(2, 0, length(header) - 1L, 0), collapse = "\t"),
    paste(header, collapse = "\t"),
    paste(row1, collapse = "\t"),
    paste(row2, collapse = "\t")
  ), fixture)
  gct <- read_protigy_stat_gct(fixture, "neuron_soma", strict_primary = TRUE)
  testthat::expect_identical(unname(gct$dimensions[["ncmat"]]), 0L)
  testthat::expect_identical(gct$ids, c("P1", "P2"))
  testthat::expect_true(gct$protein_ids_unique)
  testthat::expect_true(gct$description_aligned)
  testthat::expect_equal(nrow(gct$comparison_validation), 1L)
  testthat::expect_equal(nrow(gct$missing_required_fields), 0L)
  testthat::expect_equal(sum(gct$numeric_validation$n_non_numeric), 0L)
})
