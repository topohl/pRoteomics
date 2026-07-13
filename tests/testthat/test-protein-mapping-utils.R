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

testthat::test_that("WGCNA prevents UNMAPPED features from entering network inputs", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "manual_mapping_audit.tsv", fixed = TRUE)
  testthat::expect_match(script, "wgcna_feature_exclusion_audit.tsv", fixed = TRUE)
  testthat::expect_match(script, "validate_wgcna_expression_inputs", fixed = TRUE)
  testthat::expect_false(grepl("paste0(\"UNMAPPED_\"", script, fixed = TRUE))
})

testthat::test_that("MapThatProt and WGCNA source shared protein mapping utilities", {
  mapthatprot <- paste(readLines(testthat::test_path("..", "..", "02_id_mapping", "01_MapThatProt_batch.r"), warn = FALSE), collapse = "\n")
  wgcna <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(mapthatprot, "protein_mapping_utils.R", fixed = TRUE)
  testthat::expect_match(wgcna, "protein_mapping_utils.R", fixed = TRUE)
})

testthat::test_that("WGCNA input tables include only resolved mouse UniProt rows and audit exclusions", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tidyr")
  source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))

  male_data <- data.frame(
    gene_symbol = c("", "HUMAN_CONTAM_HUMAN", "FOO_MOUSE", "BAR_HUMAN;BAZ_MOUSE", "NOHIT_MOUSE"),
    sample_1 = c(1, 2, 3, 4, 5),
    sample_2 = c(2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  male_data$.row_id <- seq_len(nrow(male_data))

  dropped_path <- tempfile(fileext = ".tsv")
  resolved <- tokenize_wgcna_mouse_only(male_data, dropped_non_mouse_path = dropped_path)
  testthat::expect_equal(resolved$.row_id, c(3L, 4L, 5L))
  testthat::expect_equal(resolved$token_raw, c("FOO_MOUSE", "BAZ_MOUSE", "NOHIT_MOUSE"))

  dropped <- read.delim(dropped_path, check.names = FALSE)
  testthat::expect_true(any(dropped$.row_id == 1L))
  testthat::expect_true(any(dropped$.row_id == 2L))
  testthat::expect_true(any(dropped$.row_id == 4L & dropped$token_raw == "BAR_HUMAN"))

  resolved$Resolved_UNIPROT <- dplyr::case_when(
    resolved$token_raw == "FOO_MOUSE" ~ "Q9CQH5",
    resolved$token_raw == "BAZ_MOUSE" ~ "P12345",
    TRUE ~ NA_character_
  )
  resolved$strategy <- ifelse(is.na(resolved$Resolved_UNIPROT), NA_character_, "test_mapping")
  resolved$manual_mapping_used <- FALSE

  tables <- build_wgcna_input_tables(male_data, resolved)
  testthat::expect_equal(tables$male_norm$gene_symbol, c("Q9CQH5", "P12345"))
  testthat::expect_equal(tables$feature_mapping_final$.row_id, c(3L, 4L))
  testthat::expect_equal(nrow(tables$feature_mapping_pre), nrow(male_data))

  exclusion_by_row <- stats::setNames(tables$exclusion_audit$exclusion_category, tables$exclusion_audit$.row_id)
  testthat::expect_equal(unname(exclusion_by_row["1"]), "blank_identifier")
  testthat::expect_equal(unname(exclusion_by_row["2"]), "non_mouse_identifier")
  testthat::expect_equal(unname(exclusion_by_row["5"]), "unresolved_mouse_identifier")
  testthat::expect_false(any(c(3L, 4L) %in% tables$exclusion_audit$.row_id))

  expression_data <- as.data.frame(t(tables$male_norm[, -1, drop = FALSE]), check.names = FALSE)
  colnames(expression_data) <- tables$male_norm$gene_symbol
  feature_mapping_tbl <- dplyr::mutate(tables$feature_mapping_final, ProteinID = colnames(expression_data))
  testthat::expect_true(validate_wgcna_expression_inputs(expression_data, feature_mapping_tbl, tables$exclusion_audit))

  bad_expression <- expression_data
  colnames(bad_expression)[1] <- "UNMAPPED_1"
  testthat::expect_error(validate_wgcna_expression_inputs(bad_expression, feature_mapping_tbl, tables$exclusion_audit), "UNMAPPED")

  bad_mapping <- feature_mapping_tbl
  bad_mapping$resolved_uniprot[1] <- NA_character_
  testthat::expect_error(validate_wgcna_expression_inputs(expression_data, bad_mapping, tables$exclusion_audit), "without resolved UniProt")

  bad_exclusion <- dplyr::bind_rows(
    tables$exclusion_audit,
    data.frame(.row_id = 3L, original_input_token = "FOO_MOUSE", exclusion_category = "blank_identifier", exclusion_reason = "test", stringsAsFactors = FALSE)
  )
  testthat::expect_error(validate_wgcna_expression_inputs(expression_data, feature_mapping_tbl, bad_exclusion), "Excluded")
})
