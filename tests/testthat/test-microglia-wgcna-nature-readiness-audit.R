audit_repo_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
}

test_that("microglia WGCNA Nature-readiness audit script is additive and parses", {
  script <- file.path(audit_repo_root(), "06_modules_WGCNA", "12_microglia_wgcna_nature_readiness_audit.R")
  expect_true(file.exists(script))
  expect_silent(parse(file = script))

  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_match(text, "microglia_wgcna_nature_readiness", fixed = TRUE)
  expect_match(text, "immutable_reference_set_no_reassignment", fixed = TRUE)
  expect_false(grepl("saveRDS\\s*\\([^,]+,\\s*state_path", text))
  expect_false(grepl("write\\.csv\\s*\\([^,]+,\\s*registry_path", text))
})

test_that("completed audit preserves current identities and repeated-measure contract", {
  out <- file.path(audit_repo_root(), "results", "reviewer_audit", "microglia_wgcna_nature_readiness")
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")

  validation <- read.csv(file.path(out, "validation_table.csv"), check.names = FALSE)
  expect_true(all(validation$passed))
  expect_setequal(
    validation$check_id,
    c(
      "current_module_ids", "current_supermodule_ids", "no_membership_duplicates",
      "repeated_measure_unit", "sensitivity_dimensions", "feature_alignment",
      "variance_fractions", "preservation_availability", "loao_coverage",
      "higher_order_membership", "context_claim_limits", "transcriptomic_contract",
      "protected_outputs_unchanged", "no_row_multiplication"
    )
  )

  metadata <- read.csv(file.path(out, "sample_metadata_audit.csv"), check.names = FALSE)
  expect_equal(nrow(metadata), 72L)
  expect_equal(length(unique(metadata$AnimalID)), 9L)
  expect_true(all(table(metadata$AnimalID) == 8L))

  features <- read.csv(file.path(out, "sensitivity_feature_identity_audit.csv"), check.names = FALSE)
  expect_equal(nrow(features), 5201L)
  expect_equal(anyDuplicated(features$ProteinGroupID), 0L)
  expect_setequal(unique(features$ModuleID), sprintf("WGCNA_m%02d", 1:13))
  expect_true(all(features$order_identical))
})

test_that("variance, preservation and leave-one-animal-out outputs are complete", {
  out <- file.path(audit_repo_root(), "results", "reviewer_audit", "microglia_wgcna_nature_readiness")
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")

  module_vp <- read.csv(file.path(out, "module_eigengene_variance_partition.csv"), check.names = FALSE)
  super_vp <- read.csv(file.path(out, "supermodule_eigengene_variance_partition.csv"), check.names = FALSE)
  expect_setequal(unique(module_vp$entity_id), sprintf("WGCNA_m%02d", 1:13))
  expect_setequal(unique(super_vp$entity_id), sprintf("SM%02d", 1:9))
  expect_true(all(is.na(module_vp$variance_fraction[module_vp$component == "AcquisitionBatch"])))
  expect_true(all(is.na(super_vp$variance_fraction[super_vp$component == "AcquisitionBatch"])))
  estimable_sum <- aggregate(
    variance_fraction ~ level + entity_id,
    rbind(module_vp, super_vp), sum, na.rm = TRUE
  )
  expect_true(all(abs(estimable_sum$variance_fraction - 1) < 1e-8))

  robustness <- read.csv(file.path(out, "module_sensitivity_robustness.csv"), check.names = FALSE)
  expect_equal(nrow(robustness), 39L)
  expect_equal(anyDuplicated(robustness[c("sensitivity_id", "ModuleID")]), 0L)
  expect_true(all(is.finite(robustness$preservation_Zsummary)))
  expect_true(all(is.finite(robustness$preservation_median_rank)))

  loao <- read.csv(file.path(out, "module_leave_one_animal_out.csv"), check.names = FALSE)
  expect_equal(nrow(loao), 117L)
  expect_true(all(table(loao$ModuleID) == 9L))
  expect_equal(length(unique(loao$omitted_AnimalID)), 9L)
})

test_that("higher-order, biological-context and transcriptomic gates are claim-safe", {
  out <- file.path(audit_repo_root(), "results", "reviewer_audit", "microglia_wgcna_nature_readiness")
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")

  blocks <- read.csv(file.path(out, "higher_order_block_robustness.csv"), check.names = FALSE)
  expect_setequal(unique(blocks$SupermoduleID), c("SM01", "SM03", "SM09"))
  expect_equal(nrow(blocks), 12L)
  singletons <- read.csv(file.path(out, "standalone_supermodule_identity_audit.csv"), check.names = FALSE)
  expect_equal(nrow(singletons), 6L)
  expect_true(all(singletons$higher_order_block_metrics == "not_applicable"))

  context <- read.csv(file.path(out, "module_biological_context_validation.csv"), check.names = FALSE)
  expect_equal(nrow(context), 13L)
  expect_equal(anyDuplicated(context$ModuleID), 0L)
  expect_true(all(nzchar(context$molecular_process)))
  expect_true(all(nzchar(context$claim_limitation)))
  expect_false(any(grepl("microglia-intrinsic validation", context$claim_limitation, fixed = TRUE) &
                     !grepl("not", context$claim_limitation, ignore.case = TRUE)))

  contract <- read.csv(file.path(out, "transcriptomic_matching_contract_audit.csv"), check.names = FALSE)
  expect_equal(contract$status[contract$contract_element == "concordance_execution"], "not_run")
  expect_false(contract$claim_safe[contract$contract_element == "measured_transcriptomic_background"])

  protected <- read.csv(file.path(out, "protected_output_hash_audit.csv"), check.names = FALSE)
  expect_true(all(protected$unchanged))
})
