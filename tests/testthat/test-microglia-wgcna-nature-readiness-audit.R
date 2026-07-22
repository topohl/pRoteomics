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
      "variance_fractions", "preservation_descriptive_only", "loao_coverage", "loro_coverage", "animal_cluster_bootstrap", "cut_height_provenance",
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

test_that("variance, descriptive preservation and resampling outputs are complete", {
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
  expect_true(file.exists(file.path(out, "module_eigengene_variance_partition_random_effects.csv")))
  expect_true(file.exists(file.path(out, "module_eigengene_variance_partition_fixed_group.csv")))

  robustness <- read.csv(file.path(out, "module_sensitivity_robustness.csv"), check.names = FALSE)
  expect_equal(nrow(robustness), 39L)
  expect_equal(anyDuplicated(robustness[c("sensitivity_id", "ModuleID")]), 0L)
  expect_true(all(is.finite(robustness$preservation_Zsummary)))
  expect_true(all(is.finite(robustness$preservation_median_rank)))
  expect_true(all(robustness$preservation_permutation_unit == "ROI_row_unblocked"))
  expect_true(all(!robustness$repeated_measure_blocking))
  expect_true(all(!robustness$claim_gate_eligible))

  loao <- read.csv(file.path(out, "module_leave_one_animal_out.csv"), check.names = FALSE)
  expect_equal(nrow(loao), 117L)
  expect_true(all(table(loao$ModuleID) == 9L))
  expect_equal(length(unique(loao$omitted_AnimalID)), 9L)

  loro <- read.csv(file.path(out, "module_leave_one_region_out.csv"), check.names = FALSE)
  expect_equal(nrow(loro), 52L)
  expect_setequal(loro$omitted_region, c("CA1", "CA2", "CA3", "DG"))

  boot <- read.csv(gzfile(file.path(out, "module_animal_cluster_bootstrap.csv.gz")), check.names = FALSE)
  expect_true(all(boot$resampling_unit == "AnimalID_cluster"))
  expect_equal(length(unique(boot$ModuleID)), 13L)
  expect_equal(length(unique(boot$bootstrap_iteration)), 500L)
  expect_false("ROI_row" %in% unique(boot$resampling_unit))
})

test_that("higher-order, biological-context and transcriptomic gates are claim-safe", {
  out <- file.path(audit_repo_root(), "results", "reviewer_audit", "microglia_wgcna_nature_readiness")
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")

  blocks <- read.csv(file.path(out, "higher_order_block_robustness.csv"), check.names = FALSE)
  expect_setequal(unique(blocks$SupermoduleID), c("SM01", "SM03", "SM09"))
  expect_equal(nrow(blocks), 12L)
  block_ready <- read.csv(file.path(out, "higher_order_block_readiness_summary.csv"), check.names = FALSE)
  expect_setequal(block_ready$SupermoduleID, c("SM01", "SM03", "SM09"))
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

test_that("cut-height provenance, corrected classification and registry vocabulary are claim-safe", {
  root <- audit_repo_root()
  out <- file.path(root, "results", "reviewer_audit", "microglia_wgcna_nature_readiness")
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")
  provenance <- read.csv(file.path(out, "supermodule_cut_height_provenance.csv"), check.names = FALSE)
  expect_true(all(c("evidence_source", "recorded_cut_height", "authoritative", "provenance_status", "source_hash") %in% names(provenance)))
  expect_true(any(provenance$provenance_status %in% c("resolved_0.40", "resolved_0.45", "resolved_other", "unresolved_historical_provenance")))
  consensus <- read.csv(file.path(out, "module_robustness_consensus.csv"), check.names = FALSE)
  expect_true(all(c("bilateral_reproducibility_status", "spatial_adjusted_robustness_status", "animal_stability_status", "strict_nonspatial_sensitivity_status", "primary_architecture_status", "spatial_dependence_class", "claim_scope") %in% names(consensus)))
  expect_false(any(consensus$primary_architecture_status == "invalid_module"))
  registry <- read.csv(file.path(root, "config", "wgcna_labels", "microglia.csv"), check.names = FALSE)
  expect_setequal(unique(registry$label_confidence), c("high", "moderate", "low"))
  expect_true(all(registry$reviewer == "Tobias Pohl"))
  expect_true(all(nzchar(registry$proposal_prepared_by)))
})
