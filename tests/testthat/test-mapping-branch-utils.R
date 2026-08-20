testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "mapping_branch_utils.R"))

testthat::test_that("historical MapThatProt roots and paths remain unchanged", {
  roots <- resolve_mapthatprot_roots(gct_extract_root = "", mapping_output_root = "")
  testthat::expect_identical(
    roots$gct_extract_root,
    normalizePath(path_processed("01_preprocessing", "gct_extractR"), winslash = "/", mustWork = FALSE)
  )
  testthat::expect_identical(
    roots$mapping_output_root,
    normalizePath(path_processed("02_id_mapping"), winslash = "/", mustWork = FALSE)
  )
  paths <- resolve_mapthatprot_paths("microglia", "forward", roots)
  testthat::expect_identical(
    normalizePath(paths$mapped_dir, winslash = "/", mustWork = FALSE),
    normalizePath(path_processed("02_id_mapping", "mapped", "microglia", "forward", "per_file"), winslash = "/", mustWork = FALSE)
  )
  testthat::expect_identical(paths$analysis_namespace, "02_id_mapping")
})

testthat::test_that("corrected input and output root overrides resolve independently", {
  roots <- resolve_mapthatprot_roots(
    gct_extract_root = "data/processed/01_preprocessing/gct_extractR_animal_level",
    mapping_output_root = "data/processed/02_id_mapping_animal_level"
  )
  paths <- resolve_mapthatprot_paths("neuron_soma", "forward", roots)
  testthat::expect_identical(
    normalizePath(paths$raw_dir, winslash = "/", mustWork = FALSE),
    normalizePath(path_processed("01_preprocessing", "gct_extractR_animal_level", "neuron_soma", "forward"), winslash = "/", mustWork = FALSE)
  )
  testthat::expect_identical(
    normalizePath(paths$mapped_dir, winslash = "/", mustWork = FALSE),
    normalizePath(path_processed("02_id_mapping_animal_level", "mapped", "neuron_soma", "forward", "per_file"), winslash = "/", mustWork = FALSE)
  )
  testthat::expect_identical(paths$analysis_namespace, "02_id_mapping_animal_level")
  testthat::expect_false(mapping_paths_equal(paths$mapped_dir, path_processed("02_id_mapping", "mapped", "neuron_soma", "forward", "per_file")))
})

testthat::test_that("dataset and direction remain scoped in every mapping path", {
  roots <- resolve_mapthatprot_roots("custom/extract", "custom/mapping")
  for (dataset in valid_datasets()) {
    for (direction in c("forward", "reverse")) {
      paths <- resolve_mapthatprot_paths(dataset, direction, roots)
      testthat::expect_match(normalizePath(paths$raw_dir, winslash = "/", mustWork = FALSE), paste0("/", dataset, "/", direction, "$"))
      testthat::expect_match(normalizePath(paths$mapped_dir, winslash = "/", mustWork = FALSE), paste0("/mapped/", dataset, "/", direction, "/per_file$"))
      testthat::expect_match(normalizePath(paths$unmapped_dir, winslash = "/", mustWork = FALSE), paste0("/unmapped/", dataset, "/", direction, "/per_file$"))
      testthat::expect_match(normalizePath(paths$member_bridge_dir, winslash = "/", mustWork = FALSE), paste0("/member_bridge/", dataset, "/", direction, "/per_file$"))
    }
  }
  testthat::expect_error(resolve_mapthatprot_paths("microglia", "sideways", roots), "forward.*reverse")
})

testthat::test_that("corrected forward input counts fail closed", {
  expected <- animal_level_mapping_expected_counts()
  testthat::expect_identical(expected, c(neuron_neuropil = 30L, neuron_soma = 12L, microglia = 12L))
  input_root <- file.path(tempdir(), "gct_extractR_animal_level")
  for (dataset in names(expected)) {
    valid <- validate_mapthatprot_input_count(
      dataset,
      "forward",
      input_root,
      sprintf("file_%03d.csv", seq_len(expected[[dataset]]))
    )
    testthat::expect_true(valid$corrected)
    testthat::expect_identical(valid$observed, expected[[dataset]])
    testthat::expect_error(
      validate_mapthatprot_input_count(dataset, "forward", input_root, "one.csv"),
      paste0("expected ", expected[[dataset]])
    )
  }
})

testthat::test_that("custom mapping paths cannot target the canonical mapping tree", {
  roots <- resolve_mapthatprot_roots(
    "data/processed/01_preprocessing/gct_extractR_animal_level",
    "data/processed/02_id_mapping_animal_level"
  )
  canonical <- normalizePath(path_processed("02_id_mapping"), winslash = "/", mustWork = FALSE)
  for (dataset in valid_datasets()) {
    paths <- resolve_mapthatprot_paths(dataset, "forward", roots)
    for (path in c(paths$mapped_dir, paths$unmapped_dir, paths$member_bridge_dir)) {
      testthat::expect_false(startsWith(tolower(normalizePath(path, winslash = "/", mustWork = FALSE)), paste0(tolower(canonical), "/")))
    }
    testthat::expect_match(paths$tables_root, "02_id_mapping_animal_level")
    testthat::expect_match(paths$logs_root, "02_id_mapping_animal_level")
    testthat::expect_match(paths$reports_root, "02_id_mapping_animal_level")
  }
})

testthat::test_that("MapThatProt dry-run resolution report contains corrected paths and counts", {
  roots <- resolve_mapthatprot_roots(
    "data/processed/01_preprocessing/gct_extractR_animal_level",
    "data/processed/02_id_mapping_animal_level"
  )
  paths <- resolve_mapthatprot_paths("neuron_soma", "forward", roots)
  validation <- validate_mapthatprot_input_count(
    "neuron_soma", "forward", roots$gct_extract_root, sprintf("contrast_%02d.csv", seq_len(12L))
  )
  report <- mapthatprot_resolution_report(paths, validation)
  testthat::expect_identical(unname(report[["GCT extract input root"]]), roots$gct_extract_root)
  testthat::expect_identical(unname(report[["Mapping output root"]]), roots$mapping_output_root)
  testthat::expect_identical(unname(report[["Raw contrast CSV count"]]), "12")
  testthat::expect_identical(unname(report[["Expected raw contrast CSV count"]]), "12")
  testthat::expect_match(unname(report[["Mapped output directory"]]), "02_id_mapping_animal_level")
})
