# HISTORICAL SNAPSHOTS of the pre-2026-10-19 WGCNA network (old output layout),
# read from the historical root only (helper-historical.R). Moved out of
# tests/testthat (test-wgcna-stage13-downstream-integration.R,
# test-wgcna-claim-readiness.R, test-wgcna-group-effects-contract.R), whose
# current-acceptance checks now compare the current model files with the
# contracts that recorded them.

testthat::local_edition(3)

source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_identity_contract_utils.R"))
source(repo_path("R", "wgcna_group_effects_utils.R"))

historical_frozen_state_sha256 <- c(
  neuron_soma = "e218b20c97644a8b6ad07dd7697ace00f2270006b35458be61a5efb5fe243292",
  neuron_neuropil = "0fca92a64d52cb55758a6c6a47dbcdabb4422aade57c6ea4de4e903a32c5c0de",
  microglia = "62fba847986307bdb60c34e72039debc6175c65bfe3294e814fbdfeacd731926"
)

historical_state_path <- function(dataset) {
  historical_path("data", "processed", "06_modules_WGCNA", "01_WGCNA", dataset,
                  "wgcna_final_model_state.rds")
}

testthat::test_that("archived frozen states keep their literal historical hashes", {
  for (dataset in names(historical_frozen_state_sha256)) {
    path <- historical_state_path(dataset)
    skip_unless_historical(path)
    testthat::expect_identical(file_hash_sha256(path), historical_frozen_state_sha256[[dataset]],
                               info = dataset)
  }
})

testthat::test_that("archived identity contracts load against the archived states", {
  root <- local_historical_root()
  # resolve every archived path before the project root is switched to the archive
  states <- vapply(names(historical_frozen_state_sha256), historical_state_path, character(1))
  for (state in states) skip_unless_historical(state)
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  for (dataset in names(historical_frozen_state_sha256)) {
    state <- states[[dataset]]
    contract <- wgcna_group_load_identity_contract(dataset, state)
    testthat::expect_identical(contract$frozen_state_sha256, historical_frozen_state_sha256[[dataset]],
                               info = dataset)
    hashes <- utils::read.csv(file.path(dirname(contract$paths$entity), "WGCNA_identity_source_hashes.csv"),
                              stringsAsFactors = FALSE)
    hashes <- hashes[hashes$exists %in% TRUE, , drop = FALSE]
    testthat::expect_gt(nrow(hashes), 0L)
    observed <- vapply(hashes$source_path, function(p) file_hash_sha256(file.path(root, p)), character(1))
    testthat::expect_identical(unname(observed), hashes$sha256, info = dataset)
  }
})

testthat::test_that("archived Stage 05 legacy outputs keep their audited hashes", {
  root <- local_historical_root()
  for (dataset in names(historical_frozen_state_sha256)) {
    path <- historical_path("results", "tables", "06_modules_WGCNA", "group_effects", dataset,
                            "WGCNA_group_effect_legacy_output_staleness_audit.csv")
    skip_unless_historical(path)
    audit <- utils::read.csv(path, stringsAsFactors = FALSE)
    testthat::expect_gt(nrow(audit), 0L)
    observed <- vapply(audit$source_file, function(p) file_hash_sha256(file.path(root, p)), character(1))
    testthat::expect_identical(unname(observed), audit$sha256, info = dataset)
  }
})

testthat::test_that("the archived microglia state was built at supermodule cut height 0.40", {
  path <- historical_state_path("microglia")
  skip_unless_historical(path)
  state <- readRDS(path)
  testthat::expect_equal(as.numeric(state$parameters$supermodule_merge_cut_height), 0.40)
})

testthat::test_that("the archived Stage 13 table holds the 22-identity, SM09-claimable network", {
  path <- historical_path("results", "tables", "06_modules_WGCNA", "claim_readiness",
                          "microglia", "WGCNA_entity_claim_readiness.csv")
  skip_unless_historical(path)
  x <- utils::read.csv(path, check.names = FALSE)
  testthat::expect_equal(nrow(x), 22L)
  testthat::expect_setequal(x$entity_id[x$level == "module"], sprintf("WGCNA_m%02d", 1:13))
  testthat::expect_setequal(x$entity_id[x$level == "supermodule"], sprintf("SM%02d", 1:9))
  testthat::expect_equal(sum(x$claim_entity_role == "compatibility_alias"), 6L)
  testthat::expect_identical(
    x$entity_id[x$claim_entity_role == "higher_order_block" & x$separate_manuscript_claim_allowed],
    "SM09")
})
