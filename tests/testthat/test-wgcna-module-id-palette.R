testthat::test_that("WGCNA module palette is muted, unique, and not recycled", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))

  pal <- wgcna_publication_module_palette()
  testthat::expect_equal(length(pal), length(unique(unname(pal))))
  testthat::expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unname(pal))))
  testthat::expect_false(any(unname(pal) %in% c("#E6550D", "#DE2D26", "#FB6A4A", "#F16913", "#54278F", "#3F007D")))

  meta <- wgcna_module_color_metadata(c("turquoise", "blue", "brown"), pal)
  testthat::expect_equal(meta$ModuleID, c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"))
  testthat::expect_equal(meta$ModuleLegacyID, c("WGCNA_blue", "WGCNA_brown", "WGCNA_turquoise"))
  testthat::expect_equal(meta$ModuleColor, unname(pal[1:3]))
  testthat::expect_equal(meta$WGCNAInternalColor, c("blue", "brown", "turquoise"))
  testthat::expect_true(all(nzchar(meta$ModuleColorName)))
  testthat::expect_true(all(nzchar(meta$ModuleColorLabel)))

  too_many <- paste0("#", sprintf("%06X", seq_len(length(pal) + 1L)))
  testthat::expect_error(wgcna_module_color_metadata(too_many, pal), "instead of recycling")
})

testthat::test_that("WGCNA module contract decouples stable IDs from colour metadata", {
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))
  df <- data.frame(
    ModuleSet = "WGCNA",
    WGCNAInternalColor = "blue",
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

testthat::test_that("WGCNA feature integration keys internal colours separately from publication colours", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))

  internal_colors <- c("turquoise", "blue", "brown")
  publication_palette <- c(
    muted_blue = "#486A8A",
    muted_brown = "#A66E5A",
    muted_teal = "#5B9992"
  )
  meta <- wgcna_module_color_metadata(internal_colors, publication_palette)
  protein_group_ids <- paste0("PG:microglia:P", 1:4)
  merged_colors <- stats::setNames(c("turquoise", "blue", "brown", "turquoise"), protein_group_ids)
  feature_annotation <- data.frame(
    ProteinGroupID = protein_group_ids,
    ProteinID = protein_group_ids,
    member_accessions = paste0("P0000", 1:4),
    member_gene_symbols = paste0("Gene", 1:4),
    original_identifier = paste0("GENE", 1:4, "_MOUSE"),
    representative_accession = paste0("P0000", 1:4),
    representative_gene_symbol = paste0("Gene", 1:4),
    protein_group_ambiguity_class = "single_accession_single_gene",
    n_mapped_accessions = 1L,
    n_gene_symbols = 1L,
    gene_level_claim_allowed = TRUE,
    protein_level_claim_allowed = TRUE,
    mapping_status = "mapped",
    FeatureDisplayLabel = paste0("Gene", 1:4),
    stringsAsFactors = FALSE
  )

  feature_modules <- wgcna_integrate_feature_module_metadata(
    protein_group_ids,
    merged_colors,
    meta,
    feature_annotation
  )
  label_table <- dplyr::mutate(
    tibble::as_tibble(meta),
    ModuleLabel_Final = paste("Module", .data$WGCNAInternalColor)
  )
  integrated <- wgcna_join_module_labels(feature_modules, label_table)
  integrated$ModuleSet <- "WGCNA"
  integrated$kME <- 0.8

  testthat::expect_equal(nrow(integrated), length(protein_group_ids))
  testthat::expect_false(anyNA(integrated$ModuleID))
  testthat::expect_true(all(grepl("^WGCNA_m[0-9]{2,}$", integrated$ModuleID)))
  testthat::expect_equal(
    integrated$ModuleColor,
    c("#5B9992", "#486A8A", "#A66E5A", "#5B9992")
  )
  testthat::expect_false(any(integrated$WGCNAInternalColor == integrated$ModuleColor))
  testthat::expect_silent(validate_wgcna_module_definitions(integrated))

  duplicated_labels <- dplyr::bind_rows(label_table, label_table[1, ])
  testthat::expect_error(
    wgcna_join_module_labels(feature_modules, duplicated_labels),
    "one-to-one|duplicated ModuleID"
  )
})

testthat::test_that("WGCNA highlight mode greys non-selected modules", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "module_highlight_colors", fixed = TRUE)
  testthat::expect_match(script, "PROTEOMICS_WGCNA_HIGHLIGHT_MODULES", fixed = TRUE)
  testthat::expect_match(script, "row_colors_vec[!selected_row] <- \"#D7DADD\"", fixed = TRUE)
  testthat::expect_match(script, "neutral = \"#D7DADD\"", fixed = TRUE)
})
