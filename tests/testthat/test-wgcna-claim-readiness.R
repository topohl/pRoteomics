testthat::local_edition(3)

root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)

testthat::test_that("Stage 13 is a parseable non-circular microglia handoff", {
  stage13 <- file.path(root, "06_modules_WGCNA", "13_wgcna_claim_readiness.R")
  testthat::expect_silent(parse(file = stage13))
  text <- paste(readLines(stage13, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "WGCNA_entity_claim_readiness.csv", fixed = TRUE)
  testthat::expect_match(text, "conventional_preservation_claim_gate_eligible", fixed = TRUE)
  for (script in c("05_module_supermodule_group_effects.r", "06_annotate_module_microenvironment.r", "07_wgcna_interpretable_summary.r", "12_microglia_wgcna_nature_readiness_audit.R")) {
    source_text <- paste(readLines(file.path(root, "06_modules_WGCNA", script), warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("13_wgcna_claim_readiness|claim_readiness", source_text))
  }
})

testthat::test_that("completed Stage 13 output has exact current identities without join multiplication", {
  path <- file.path(root, "results", "tables", "06_modules_WGCNA", "claim_readiness", "microglia", "WGCNA_entity_claim_readiness.csv")
  testthat::skip_if_not(file.exists(path), "Stage 13 has not been run")
  x <- read.csv(path, check.names = FALSE)
  testthat::expect_equal(nrow(x), 22L)
  testthat::expect_equal(anyDuplicated(x[c("dataset", "level", "entity_id")]), 0L)
  testthat::expect_setequal(x$entity_id[x$level == "module"], sprintf("WGCNA_m%02d", 1:13))
  testthat::expect_setequal(x$entity_id[x$level == "supermodule"], sprintf("SM%02d", 1:9))
  testthat::expect_false(any(x$conventional_preservation_claim_gate_eligible %in% TRUE))
})
