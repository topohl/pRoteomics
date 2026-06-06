testthat::test_that("manual mapping table parsing tolerates supported column aliases", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))
  manual_file <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(input = "foo_mouse", mapped = "Q9CQH5", stringsAsFactors = FALSE),
    manual_file,
    row.names = FALSE
  )
  parsed <- read_manual_mapping_table(manual_file)
  testthat::expect_equal(parsed$gene_symbol, "FOO_MOUSE")
  testthat::expect_equal(parsed$mapped_gene_symbol, "Q9CQH5")
  testthat::expect_equal(attr(parsed, "status"), "loaded")
})

testthat::test_that("manual override resolves before UNMAPPED fallback", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))
  resolved <- data.frame(
    token_raw = "foo_mouse",
    token_base = "FOO",
    Resolved_UNIPROT = NA_character_,
    strategy = NA_character_,
    stringsAsFactors = FALSE
  )
  manual <- data.frame(gene_symbol = "FOO_MOUSE", mapped_gene_symbol = "Q9CQH5", stringsAsFactors = FALSE)
  entry_map <- data.frame(UNIPROT = character(), entry_base = character(), stringsAsFactors = FALSE)
  gene_map <- data.frame(input = character(), primaryAccession = character(), stringsAsFactors = FALSE)
  out <- apply_manual_mapping_override(resolved, manual, entry_map, gene_map)
  testthat::expect_equal(out$data$Resolved_UNIPROT, "Q9CQH5")
  testthat::expect_true(out$data$manual_mapping_used)
  testthat::expect_false(startsWith(out$data$Resolved_UNIPROT, "UNMAPPED"))
})

testthat::test_that("WGCNA excludes UNMAPPED features from biological interpretation inputs", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "mapping_status != \"unmapped\"", fixed = TRUE)
  testthat::expect_match(script, "manual_mapping_audit.tsv", fixed = TRUE)
  testthat::expect_match(script, "n_unmapped_features", fixed = TRUE)
})

testthat::test_that("MapThatProt and WGCNA source shared protein mapping utilities", {
  mapthatprot <- paste(readLines(testthat::test_path("..", "..", "02_id_mapping", "01_MapThatProt_batch.r"), warn = FALSE), collapse = "\n")
  wgcna <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(mapthatprot, "protein_mapping_utils.R", fixed = TRUE)
  testthat::expect_match(wgcna, "protein_mapping_utils.R", fixed = TRUE)
})
