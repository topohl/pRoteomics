source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "ewce_gene_set_engine.R"))

# Module EWCE gene identity follows the Stage 02 resolver symbol policy: exact
# SYMBOL, then a unique case-insensitive SYMBOL, then a unique ALIAS that
# returns its official SYMBOL. org.Mm.eg.db is pinned in renv.lock.

testthat::test_that("symbols resolve exactly, then case-insensitively", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  testthat::expect_identical(
    ewce_to_mouse_symbols(c("Gria1", "GRIA1", "gria1", "mt-Co1", "MT-CO1")),
    c("Gria1", "Gria1", "Gria1", "mt-Co1", "mt-Co1"))
  testthat::expect_identical(
    ewce_to_mouse_symbols(c("", NA, "NOT_A_MOUSE_GENE")), rep(NA_character_, 3))
  testthat::expect_identical(ewce_to_mouse_symbols(character()), character())
  testthat::expect_identical(
    ewce_to_mouse_symbols(c("GRIA1", "NOT_A_MOUSE_GENE"), keep_unmapped = TRUE),
    c("Gria1", "NOT_A_MOUSE_GENE"))
})

testthat::test_that("an alias returns its official symbol, never the alias text", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  # Legacy names that dropped out of the specificity reference while the alias
  # text itself was returned.
  testthat::expect_identical(
    ewce_to_mouse_symbols(c("CA2", "CA4", "GPI", "TXN", "RAB7A", "UGT8", "H1-0")),
    c("Car2", "Car4", "Gpi1", "Txn1", "Rab7", "Ugt8a", "H1f0"))
  # Syp is a symbol, so its alias collision with Ptpn11 is never consulted.
  testthat::expect_identical(ewce_to_mouse_symbols("SYP"), "Syp")

  db <- org.Mm.eg.db::org.Mm.eg.db
  sym <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  ali <- AnnotationDbi::keys(db, keytype = "ALIAS")
  alias_only <- unique(toupper(ali[!toupper(ali) %in% toupper(sym)]))
  testthat::expect_true(all(ewce_to_mouse_symbols(alias_only) %in% c(sym, NA)))
})

testthat::test_that("a key naming more than one gene stays unresolved", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  # SKT (Etl4, Chsy1) and TF (F3, Itpr3, Trf) are ambiguous aliases.
  testthat::expect_identical(ewce_to_mouse_symbols(c("SKT", "TF")),
                             c(NA_character_, NA_character_))
  # GS matches the symbols gs and Gs and must not fall through to an alias.
  testthat::expect_identical(ewce_to_mouse_symbols(c("GS", "Gs", "gs")),
                             c(NA_character_, "Gs", "gs"))
})
