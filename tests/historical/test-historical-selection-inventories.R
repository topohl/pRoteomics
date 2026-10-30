# HISTORICAL SNAPSHOT of the supplementary selection inventories as exported
# before 2026-10-19 (moved from tests/testthat/test-display-selection-disclosure.R,
# which now checks the current export against its producer). These are the row
# counts the Phase 6H closure accepted.

testthat::local_edition(3)

testthat::test_that("the archived selection inventories keep the Phase 6H row counts", {
  man <- historical_path("exports", "publication_source_data", "manifest.csv")
  skip_unless_historical(man)
  expected <- c(pathway_enrichment_inventory.csv = 203073L,
                leading_edge_protein_inventory.csv = 282296L,
                leading_edge_protein_recurrence.csv = 12227L,
                pathway_enrichment_inventory_fdr_supported.csv = 3559L,
                inventory_data_dictionary.csv = 72L,
                display_selection_disclosure.csv = 7L)
  m <- utils::read.csv(man, stringsAsFactors = FALSE)
  s <- m[m$publication_id == "supplementary_selection_inventories", , drop = FALSE]
  testthat::expect_identical(nrow(s), 6L)
  for (nm in names(expected)) {
    row <- s[basename(s$exported_file) == nm, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L, info = nm)
    testthat::expect_identical(as.integer(row$rows[1]), expected[[nm]], info = nm)
  }
})
