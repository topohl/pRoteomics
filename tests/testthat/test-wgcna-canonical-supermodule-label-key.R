testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

stage05_script_assignments <- function() {
  expressions <- parse(testthat::test_path(
    "..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"
  ))
  names <- vapply(expressions, function(expression) {
    if (
      is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]])
    ) as.character(expression[[2]]) else NA_character_
  }, character(1))
  stats::na.omit(names)
}

testthat::test_that("Stage 05 uses only the v5 quantitative implementation", {
  legacy_semantic_helpers <- c(
    "canonical_supermodule_label_lookup",
    "empty_canonical_supermodule_label_lookup",
    "apply_canonical_supermodule_labels"
  )
  testthat::expect_length(
    intersect(stage05_script_assignments(), legacy_semantic_helpers), 0L
  )
  testthat::expect_identical(
    wgcna_group_effects_contract_version(),
    "wgcna_group_effects_phase2b_v5"
  )
  testthat::expect_true(is.function(wgcna_group_construct_endpoints))
})

testthat::test_that("v5 preserves stable endpoint identity without Stage 05 biological relabeling", {
  schema <- wgcna_group_primary_schema()
  testthat::expect_true(all(c(
    "dataset", "level", "endpoint_id", "module_id", "supermodule_id",
    "canonical_claim_entity_id", "support_source_entity_id",
    "membership_version", "identity_contract_version",
    "identity_contract_sha256", "frozen_state_sha256"
  ) %in% names(schema)))

  body_text <- paste(deparse(body(wgcna_group_construct_endpoints)), collapse = "\n")
  testthat::expect_match(body_text, "contract\\$membership")
  testthat::expect_match(
    body_text,
    "Published supermodule composition does not exactly reproduce the identity contract",
    fixed = TRUE
  )
  testthat::expect_match(
    wgcna_group_legacy_reason("selected_sus_res_supermodule_contents.csv"),
    "legacy selected-supermodule"
  )
})
