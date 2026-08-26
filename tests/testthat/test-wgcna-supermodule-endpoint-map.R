testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

testthat::test_that("canonical endpoint construction is membership and bridge keyed", {
  testthat::expect_true(is.function(wgcna_group_build_module_bridge))
  testthat::expect_true(is.function(wgcna_group_construct_endpoints))

  body_text <- paste(deparse(body(wgcna_group_construct_endpoints)), collapse = "\n")
  testthat::expect_match(body_text, "contract\\$membership")
  testthat::expect_match(body_text, "module_id")
  testthat::expect_match(body_text, "supermodule_id")
  testthat::expect_match(body_text, "many-to-one", fixed = TRUE)
  testthat::expect_match(body_text, "identity contract", fixed = TRUE)
})

testthat::test_that("published endpoint provenance retains stable compound identity", {
  schema <- wgcna_group_primary_schema()
  testthat::expect_true(all(c(
    "dataset", "level", "endpoint_id", "module_id", "supermodule_id",
    "canonical_claim_entity_id", "claim_entity_role",
    "support_source_entity_id", "independent_hypothesis"
  ) %in% names(schema)))

  allowed <- wgcna_group_allowed_output_names()
  testthat::expect_true("supermodule_composition.csv" %in% allowed)
  testthat::expect_true(
    "WGCNA_group_effect_endpoint_provenance.csv" %in% allowed
  )
})

testthat::test_that("superseded endpoint-map helpers cannot shadow v5", {
  legacy <- c("make_endpoint_maps", "join_supermodule_endpoint_labels")
  testthat::expect_false(any(vapply(
    legacy, exists, logical(1), mode = "function", inherits = FALSE
  )))
  testthat::expect_match(
    wgcna_group_legacy_reason("module_to_supermodule_map.csv"),
    "canonical identity contract supersedes it",
    fixed = TRUE
  )
})
