testthat::test_that("WGCNA module palette is muted, unique, and not recycled", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))

  pal <- wgcna_publication_module_palette()
  testthat::expect_equal(length(pal), length(unique(unname(pal))))
  testthat::expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unname(pal))))
  testthat::expect_false(any(unname(pal) %in% c("#E6550D", "#DE2D26", "#FB6A4A", "#F16913", "#54278F", "#3F007D")))

  meta <- wgcna_module_color_metadata(c("#486A8A", "#5B9992", "#9A8977"), pal)
  testthat::expect_equal(meta$ModuleID, c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"))
  testthat::expect_equal(meta$ModuleLegacyID, c("WGCNA_#486A8A", "WGCNA_#5B9992", "WGCNA_#9A8977"))
  testthat::expect_equal(meta$ModuleColor, unname(pal[1:3]))
  testthat::expect_true(all(nzchar(meta$ModuleColorName)))
  testthat::expect_true(all(nzchar(meta$ModuleColorLabel)))

  too_many <- paste0("#", sprintf("%06X", seq_len(length(pal) + 1L)))
  testthat::expect_error(wgcna_module_color_metadata(too_many, pal), "instead of recycling")
})

testthat::test_that("WGCNA module contract decouples stable IDs from colour metadata", {
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))
  df <- data.frame(
    ModuleSet = "WGCNA",
    ModuleID = "WGCNA_m01",
    ModuleLegacyID = "WGCNA_#486A8A",
    ModuleColor = "#486A8A",
    ModuleColorName = "deep_steel_blue",
    ModuleColorLabel = "deep steel blue",
    ProteinGroupID = "PG:microglia:P1",
    ProteinID = "PG:microglia:P1",
    member_accessions = "Q9CQH5",
    member_gene_symbols = "Foo",
    original_identifier = "FOO_MOUSE",
    representative_accession = "Q9CQH5",
    representative_gene_symbol = "Foo",
    protein_group_ambiguity_class = "single_accession_single_gene",
    n_mapped_accessions = 1L,
    n_gene_symbols = 1L,
    gene_level_claim_allowed = TRUE,
    protein_level_claim_allowed = TRUE,
    mapping_status = "mapped",
    FeatureDisplayLabel = "Foo",
    kME = 0.7,
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(validate_wgcna_module_definitions(df))

  bad <- df
  bad$ModuleID <- "WGCNA_#486A8A"
  testthat::expect_error(validate_wgcna_module_definitions(bad), "stable ModuleID")

  old <- data.frame(ModuleID = c("WGCNA_#486A8A", "WGCNA_m02"), stringsAsFactors = FALSE)
  lookup <- data.frame(
    ModuleID = c("WGCNA_m01", "WGCNA_m02"),
    ModuleLegacyID = c("WGCNA_#486A8A", "WGCNA_#5B9992"),
    stringsAsFactors = FALSE
  )
  normalized <- wgcna_normalize_module_ids(old, lookup)
  testthat::expect_equal(normalized$ModuleID, c("WGCNA_m01", "WGCNA_m02"))
  testthat::expect_equal(normalized$ModuleLegacyID, c("WGCNA_#486A8A", "WGCNA_#5B9992"))
})

testthat::test_that("WGCNA highlight mode greys non-selected modules", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "module_highlight_colors", fixed = TRUE)
  testthat::expect_match(script, "PROTEOMICS_WGCNA_HIGHLIGHT_MODULES", fixed = TRUE)
  testthat::expect_match(script, "row_colors_vec[!selected_row] <- \"#D7DADD\"", fixed = TRUE)
  testthat::expect_match(script, "neutral = \"#D7DADD\"", fixed = TRUE)
})
