testthat::test_that("path helpers resolve inside repository", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  root <- repo_root()
  testthat::expect_true(file.exists(file.path(root, "README.md")))
  testthat::expect_match(repo_path("R", "paths.R"), "R[/\\\\]paths[.]R$")
  testthat::expect_equal(relative_to(repo_path("README.md")), "README.md")
  testthat::expect_match(safe_filename("a/b c"), "a_b_c")
})

testthat::test_that("strict input resolver audits and refuses latest fallback", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "schema_validation.R"))
  old_strict <- Sys.getenv("PROTEOMICS_STRICT_INPUTS", unset = NA_character_)
  old_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = NA_character_)
  tmp <- tempfile("strict-inputs-")
  dir.create(file.path(tmp, "results", "reviewer_audit"), recursive = TRUE)
  dir.create(file.path(tmp, "inst", "schemas"), recursive = TRUE)
  file.copy(
    testthat::test_path("..", "..", "inst", "schemas", "input_resolution_audit.yml"),
    file.path(tmp, "inst", "schemas", "input_resolution_audit.yml"),
    overwrite = TRUE
  )
  latest_dir <- file.path(tmp, "latest")
  dir.create(latest_dir, recursive = TRUE)
  writeLines("fallback", file.path(latest_dir, "input_20260101.csv"))
  Sys.setenv(PROTEOMICS_PROJECT_ROOT = tmp)
  on.exit({
    if (is.na(old_strict)) Sys.unsetenv("PROTEOMICS_STRICT_INPUTS") else Sys.setenv(PROTEOMICS_STRICT_INPUTS = old_strict)
    if (is.na(old_root)) Sys.unsetenv("PROTEOMICS_PROJECT_ROOT") else Sys.setenv(PROTEOMICS_PROJECT_ROOT = old_root)
    unlink(tmp, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  Sys.setenv(PROTEOMICS_STRICT_INPUTS = "true")
  testthat::expect_error(
    suppressWarnings(resolve_input_path(
      input_name = "claim_input",
      expected_path = file.path(tmp, "canonical.csv"),
      latest_roots = latest_dir,
      latest_pattern = "\\.csv$",
      required = TRUE,
      script = "test-script",
      dataset = "global",
      stage = "test",
      producer_script_or_artifact_id = "test-producer"
    )),
    "Strict input mode forbids newest-file fallback"
  )
  audit <- utils::read.csv(file.path(tmp, "results", "reviewer_audit", "input_resolution_audit.csv"), stringsAsFactors = FALSE)
  testthat::expect_true("latest_forbidden_strict" %in% audit$resolution_mode)
  testthat::expect_silent(validate_table_schema(audit, "input_resolution_audit", strict = TRUE))
})
