testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))

# Fixture reproducing the current empirical ROI marker v2 schema:
#   - descriptive-only shared panel with empty claim-gated GeneSymbol but
#     populated mapped_gene_symbol
#   - claim-ready inferential enriched panels
#   - mapping-blocked rows (no usable identity in either field)
#   - no empirical_neuropil_sensitive_high_confidence panel (zero members)
v2_marker_fixture <- function() {
  contract <- "empirical_roi_marker_v2_animal_paired_limma"
  rbind(
    # descriptive-only shared panel: GeneSymbol deliberately NA
    data.frame(
      marker_set = "empirical_microglia_neuropil_shared",
      GeneSymbol = NA_character_,
      mapped_gene_symbol = c("Shared1", "Shared2", "Shared2"),
      marker_evidence_type = "descriptive_only",
      claim_allowed = FALSE,
      gene_mapping_claim_allowed = TRUE,
      marker_contract_version = contract,
      marker_source = "empirical_roi_marker_sets",
      stringsAsFactors = FALSE
    ),
    # descriptive-only shared panel: mapping-blocked row (no identity anywhere)
    data.frame(
      marker_set = "empirical_microglia_neuropil_shared",
      GeneSymbol = NA_character_,
      mapped_gene_symbol = NA_character_,
      marker_evidence_type = "descriptive_only",
      claim_allowed = FALSE,
      gene_mapping_claim_allowed = FALSE,
      marker_contract_version = contract,
      marker_source = "empirical_roi_marker_sets",
      stringsAsFactors = FALSE
    ),
    # inferential claim-ready panel
    data.frame(
      marker_set = "empirical_microglia_roi_enriched",
      GeneSymbol = c("Enr1", "Enr2"),
      mapped_gene_symbol = c("Enr1", "Enr2"),
      marker_evidence_type = "inferential",
      claim_allowed = TRUE,
      gene_mapping_claim_allowed = TRUE,
      marker_contract_version = contract,
      marker_source = "empirical_roi_marker_sets",
      stringsAsFactors = FALSE
    ),
    # inferential panel with a mapping-blocked member that must not be promoted
    data.frame(
      marker_set = "empirical_microglia_roi_enriched",
      GeneSymbol = NA_character_,
      mapped_gene_symbol = NA_character_,
      marker_evidence_type = "inferential",
      claim_allowed = FALSE,
      gene_mapping_claim_allowed = FALSE,
      marker_contract_version = contract,
      marker_source = "empirical_roi_marker_sets",
      stringsAsFactors = FALSE
    )
  )
}

with_empirical_fixture <- function(df, code) {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE, na = "")
  old <- Sys.getenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE") else Sys.setenv(PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE = old)
  }, add = TRUE)
  Sys.setenv(PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE = f)
  force(code)
}

testthat::test_that("v2 reader keys descriptive panels on mapped_gene_symbol", {
  with_empirical_fixture(v2_marker_fixture(), {
    emp <- read_empirical_roi_marker_sets()
    testthat::expect_equal(
      attr(emp, "empirical_marker_identity_field"), "mapped_gene_symbol"
    )
    testthat::expect_equal(
      attr(emp, "empirical_marker_contract_version"),
      "empirical_roi_marker_v2_animal_paired_limma"
    )
    # both mapping-blocked rows are dropped, nothing else
    testthat::expect_equal(attr(emp, "empirical_marker_rows_mapping_blocked"), 2L)
    testthat::expect_equal(nrow(emp), 5L)
  })
})

testthat::test_that("descriptive-only shared panel survives empty GeneSymbol", {
  with_empirical_fixture(v2_marker_fixture(), {
    sets <- load_wgcna_marker_sets(quiet = TRUE)
    testthat::expect_true("empirical_microglia_neuropil_shared" %in% names(sets))
    # deterministic unique gene set, duplicates collapsed
    testthat::expect_equal(
      sort(sets[["empirical_microglia_neuropil_shared"]]),
      c("Shared1", "Shared2")
    )
  })
})

testthat::test_that("inferential panels exclude mapping-blocked members", {
  with_empirical_fixture(v2_marker_fixture(), {
    sets <- load_wgcna_marker_sets(quiet = TRUE)
    testthat::expect_equal(
      sort(sets[["empirical_microglia_roi_enriched"]]), c("Enr1", "Enr2")
    )
  })
})

testthat::test_that("panel evidence roles keep descriptive sets out of claim evidence", {
  with_empirical_fixture(v2_marker_fixture(), {
    sets <- load_wgcna_marker_sets(quiet = TRUE)
    roles <- attr(sets, "empirical_marker_panel_roles")
    shared <- roles[roles$marker_set == "empirical_microglia_neuropil_shared", ]
    enriched <- roles[roles$marker_set == "empirical_microglia_roi_enriched", ]
    testthat::expect_equal(shared$marker_evidence_type, "descriptive_only")
    testthat::expect_equal(shared$claim_eligible_members, 0L)
    testthat::expect_equal(enriched$marker_evidence_type, "inferential")
    testthat::expect_equal(enriched$claim_eligible_members, 2L)
  })
})

testthat::test_that("zero-member panels are absent and never fabricated", {
  with_empirical_fixture(v2_marker_fixture(), {
    sets <- load_wgcna_marker_sets(quiet = TRUE)
    testthat::expect_false(
      "empirical_neuropil_sensitive_high_confidence" %in% names(sets)
    )
  })
})

testthat::test_that("reader fails closed on pre-v2 schema and contract mismatch", {
  pre_v2 <- data.frame(
    marker_set = "empirical_microglia_roi_enriched",
    GeneSymbol = "Aif1",
    marker_source = "legacy",
    stringsAsFactors = FALSE
  )
  with_empirical_fixture(pre_v2, {
    testthat::expect_error(
      read_empirical_roi_marker_sets(),
      "missing required column"
    )
  })

  wrong_contract <- v2_marker_fixture()
  wrong_contract$marker_contract_version <- "empirical_roi_marker_v1_legacy"
  with_empirical_fixture(wrong_contract, {
    testthat::expect_error(
      read_empirical_roi_marker_sets(), "contract mismatch"
    )
  })
})
