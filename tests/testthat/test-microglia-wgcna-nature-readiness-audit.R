audit_repo_root <- function() {
  normalizePath(repo_root(), winslash = "/", mustWork = TRUE)
}

# Current canonical Stage 12 audit (analysis/wgcna/audit_microglia_module_claims.R).
audit_out <- function() {
  current_wgcna("audit_microglia_module_claims", "global")
}

test_that("microglia WGCNA Nature-readiness audit script is additive and parses", {
  script <- file.path(audit_repo_root(), "analysis/wgcna", "audit_microglia_module_claims.R")
  expect_true(file.exists(script))
  expect_silent(parse(file = script))

  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_match(text, "microglia_wgcna_nature_readiness", fixed = TRUE)
  expect_match(text, "immutable_reference_set_no_reassignment", fixed = TRUE)
  expect_false(grepl("saveRDS\\s*\\([^,]+,\\s*state_path", text))
  expect_false(grepl("write\\.csv\\s*\\([^,]+,\\s*registry_path", text))
})

test_that("completed audit preserves current identities and repeated-measure contract", {
  out <- audit_out()
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")
  shape <- current_network_shape("microglia")

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

  # study design (9 animals x 8 microglia ROI samples), not network shape
  metadata <- read.csv(file.path(out, "sample_metadata_audit.csv"), check.names = FALSE)
  expect_equal(nrow(metadata), 72L)
  expect_equal(length(unique(metadata$AnimalID)), 9L)
  expect_true(all(table(metadata$AnimalID) == 8L))

  # every current module member, once, in the current Stage 01 order
  modules_long <- read_current_csv(current_wgcna("build_wgcna_modules", "microglia", "modules",
                                                 "WGCNA_modules_long.csv"))
  features <- read.csv(file.path(out, "sensitivity_feature_identity_audit.csv"), check.names = FALSE)
  expect_equal(nrow(features), length(unique(modules_long$ProteinGroupID)))
  expect_setequal(features$ProteinGroupID, modules_long$ProteinGroupID)
  expect_equal(anyDuplicated(features$ProteinGroupID), 0L)
  expect_setequal(unique(features$ModuleID), shape$modules)
  expect_true(all(features$order_identical))
})

test_that("variance, descriptive preservation and resampling outputs are complete", {
  out <- audit_out()
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")
  shape <- current_network_shape("microglia")
  n_modules <- length(shape$modules)

  module_vp <- read.csv(file.path(out, "module_eigengene_variance_partition.csv"), check.names = FALSE)
  super_vp <- read.csv(file.path(out, "supermodule_eigengene_variance_partition.csv"), check.names = FALSE)
  expect_setequal(unique(module_vp$entity_id), shape$modules)
  expect_setequal(unique(super_vp$entity_id), shape$supermodules)
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
  # the declared sensitivity analyses (analysis/wgcna/audit_microglia_module_claims.R)
  expect_setequal(unique(robustness$sensitivity_id),
                  c("hemisphere_collapsed", "spatial_adjusted", "within_animal_spatial_adjusted"))
  expect_equal(nrow(robustness), n_modules * 3L)
  expect_setequal(unique(robustness$ModuleID), shape$modules)
  expect_equal(anyDuplicated(robustness[c("sensitivity_id", "ModuleID")]), 0L)
  expect_true(all(is.finite(robustness$preservation_Zsummary)))
  expect_true(all(is.finite(robustness$preservation_median_rank)))
  expect_true(all(robustness$preservation_permutation_unit == "ROI_row_unblocked"))
  expect_true(all(!robustness$repeated_measure_blocking))
  expect_true(all(!robustness$claim_gate_eligible))

  loao <- read.csv(file.path(out, "module_leave_one_animal_out.csv"), check.names = FALSE)
  n_animals <- length(unique(loao$omitted_AnimalID))
  expect_equal(n_animals, 9L)
  expect_equal(nrow(loao), n_modules * n_animals)
  expect_true(all(table(loao$ModuleID) == n_animals))

  loro <- read.csv(file.path(out, "module_leave_one_region_out.csv"), check.names = FALSE)
  expect_setequal(loro$omitted_region, c("CA1", "CA2", "CA3", "DG"))
  expect_equal(nrow(loro), n_modules * 4L)

  boot <- read.csv(gzfile(file.path(out, "module_animal_cluster_bootstrap.csv.gz")), check.names = FALSE)
  expect_true(all(boot$resampling_unit == "AnimalID_cluster"))
  expect_setequal(unique(boot$ModuleID), shape$modules)
  expect_equal(length(unique(boot$bootstrap_iteration)), 500L)
  expect_false("ROI_row" %in% unique(boot$resampling_unit))
})

test_that("higher-order, biological-context and transcriptomic gates are claim-safe", {
  out <- audit_out()
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")
  shape <- current_network_shape("microglia")

  blocks <- read.csv(file.path(out, "higher_order_block_robustness.csv"), check.names = FALSE)
  expect_setequal(unique(blocks$SupermoduleID), shape$blocks)
  expect_setequal(unique(blocks$sensitivity_id),
                  c("primary", "hemisphere_collapsed", "spatial_adjusted", "within_animal_spatial_adjusted"))
  expect_equal(nrow(blocks), length(shape$blocks) * 4L)
  expect_equal(anyDuplicated(blocks[c("sensitivity_id", "SupermoduleID")]), 0L)
  block_ready <- read.csv(file.path(out, "higher_order_block_readiness_summary.csv"), check.names = FALSE)
  expect_setequal(block_ready$SupermoduleID, shape$blocks)
  singletons <- read.csv(file.path(out, "standalone_supermodule_identity_audit.csv"), check.names = FALSE)
  expect_setequal(singletons$SupermoduleID, names(shape$singletons))
  expect_true(all(singletons$higher_order_block_metrics == "not_applicable"))

  context <- read.csv(file.path(out, "module_biological_context_validation.csv"), check.names = FALSE)
  expect_setequal(context$ModuleID, shape$modules)
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

test_that("cut-height provenance and corrected classification are claim-safe", {
  out <- audit_out()
  skip_if_not(file.exists(file.path(out, "validation_table.csv")), "audit has not been run")
  provenance <- read.csv(file.path(out, "supermodule_cut_height_provenance.csv"), check.names = FALSE)
  expect_true(all(c("evidence_source", "recorded_cut_height", "authoritative", "provenance_status", "source_hash") %in% names(provenance)))
  expect_true(any(provenance$provenance_status %in% c("resolved_0.40", "resolved_0.45", "resolved_other", "unresolved_historical_provenance")))
  consensus <- read.csv(file.path(out, "module_robustness_consensus.csv"), check.names = FALSE)
  expect_true(all(c("bilateral_reproducibility_status", "spatial_adjusted_robustness_status", "animal_stability_status", "strict_nonspatial_sensitivity_status", "primary_architecture_status", "spatial_dependence_class", "claim_scope") %in% names(consensus)))
  expect_false(any(consensus$primary_architecture_status == "invalid_module"))
})

test_that("the audit read the registry in force, its module rows unchanged since, with a valid vocabulary", {
  out <- audit_out()
  skip_if_not(file.exists(file.path(out, "input_hashes.csv")), "audit has not been run")
  # reviewed mode: config/wgcna_labels/microglia.csv (an error if missing);
  # declared provisional mode: the registry named by the Stage 07 manifest
  reg <- wgcna_registry_in_force("microglia")
  hashes <- read.csv(file.path(out, "input_hashes.csv"), check.names = FALSE)
  row <- hashes[grepl("(wgcna_labels/microglia|provisional_label_registry)[.]csv$", hashes$input_path), , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_identical(
    normalizePath(file.path(audit_repo_root(), row$input_path), winslash = "/", mustWork = FALSE),
    normalizePath(reg$path, winslash = "/", mustWork = FALSE))
  # The audit reads the registry's module rows only (label, subcellular and ROI
  # context; analysis/wgcna/audit_microglia_module_claims.R). Those must be what
  # it wrote: an edit to them after the run makes the audit stale, while an edit
  # to supermodule rows or review provenance does not.
  registry <- read.csv(reg$path, check.names = FALSE)
  ctx <- read.csv(file.path(out, "module_biological_context_validation.csv"), check.names = FALSE)
  mods <- registry[registry$level == "module", , drop = FALSE]
  k <- match(ctx$ModuleID, mods$entity_id)
  expect_false(anyNA(k))
  expect_setequal(ctx$ModuleID, mods$entity_id)
  expect_identical(ctx$molecular_process, mods$reviewed_biological_label[k])
  expect_identical(ctx$subcellular_localization, mods$subcellular_context[k])
  expect_identical(ctx$reviewed_roi_context, mods$roi_context[k])
  expect_true(all(registry$label_confidence %in% c("high", "moderate", "low")))
  expect_true(all(nzchar(registry$proposal_prepared_by)))
  if (identical(reg$mode, "reviewed")) {
    expect_true(all(registry$adjudication_status == "reviewed"))
    complete <- reviewed_rows_are_complete(registry)
    expect_true(all(complete), info = paste("incomplete review:", paste(names(complete)[!complete], collapse = ", ")))
  } else {
    expect_true(all(registry$adjudication_status == "provisional_unreviewed"))
  }
})
