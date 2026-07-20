testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
source(testthat::test_path("..", "..", "R", "module_contracts.R"))

hub_handoff_fixture <- function() {
  modules <- tibble::tibble(
    dataset = "microglia",
    ProteinGroupID = paste0("PG:microglia:P", 1:6),
    ModuleID = rep(c("WGCNA_m01", "WGCNA_m02"), each = 3),
    WGCNAInternalColor = rep(c("blue", "brown"), each = 3),
    ModuleColor = rep(c("#486A8A", "#A66E5A"), each = 3),
    abs_kME = c(0.90, 0.75, 0.40, 0.88, 0.70, 0.35),
    is_top_hub_25 = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
  )
  annotation <- tibble::tibble(
    dataset = "microglia",
    ModuleID = c("WGCNA_m01", "WGCNA_m02"),
    module_eigengene = c("MEblue", "MEbrown"),
    SupermoduleID = c("SM01", "SM02"),
    Supermodule_DisplayLabel = c("SM01 - blue", "SM02 - brown"),
    present_in_dataset = TRUE
  )
  list(modules = modules, annotation = annotation)
}

testthat::test_that("hub-overlap handoff joins authoritative eigengenes by stable module keys", {
  x <- hub_handoff_fixture()
  joined <- wgcna_join_supermodule_hub_handoff(
    x$modules,
    x$annotation,
    merged_me_names = c("MEblue", "MEbrown")
  )
  sets <- wgcna_build_hub_module_sets(
    x$modules,
    x$annotation,
    merged_me_names = c("MEblue", "MEbrown"),
    top_n = 2L
  )

  testthat::expect_equal(nrow(joined), nrow(x$modules))
  testthat::expect_false(anyNA(joined$module_eigengene))
  testthat::expect_equal(
    unique(joined[, c("ModuleID", "module_eigengene")]),
    x$annotation[, c("ModuleID", "module_eigengene")]
  )
  testthat::expect_equal(nrow(sets), 2L)
  testthat::expect_equal(sort(sets$ModuleID), c("WGCNA_m01", "WGCNA_m02"))
  testthat::expect_true(all(lengths(sets$top_hub_proteins) == 2L))
  testthat::expect_equal(sets$ModuleColor, c("#486A8A", "#A66E5A"))
  testthat::expect_equal(sets$WGCNAInternalColor, c("blue", "brown"))
})

testthat::test_that("hub-overlap handoff rejects ambiguous or missing module mappings", {
  x <- hub_handoff_fixture()
  duplicate_lookup <- dplyr::bind_rows(x$annotation, x$annotation[1, ])
  testthat::expect_error(
    wgcna_join_supermodule_hub_handoff(
      x$modules,
      duplicate_lookup,
      merged_me_names = c("MEblue", "MEbrown")
    ),
    "duplicated dataset \\+ ModuleID"
  )

  testthat::expect_error(
    wgcna_join_supermodule_hub_handoff(
      x$modules,
      x$annotation[1, ],
      merged_me_names = c("MEblue", "MEbrown")
    ),
    "missing module mappings"
  )

  duplicate_protein <- x$modules
  duplicate_protein$ProteinGroupID[[6]] <- duplicate_protein$ProteinGroupID[[1]]
  testthat::expect_error(
    wgcna_join_supermodule_hub_handoff(
      duplicate_protein,
      x$annotation,
      merged_me_names = c("MEblue", "MEbrown")
    ),
    "duplicated ProteinGroupID"
  )
})
