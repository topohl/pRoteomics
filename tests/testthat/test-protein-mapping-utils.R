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

testthat::test_that("WGCNA separates quantitative inclusion from annotation eligibility", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "manual_mapping_audit.tsv", fixed = TRUE)
  testthat::expect_match(script, "WGCNA_feature_universe_audit.csv", fixed = TRUE)
  testthat::expect_match(script, "validate_wgcna_feature_alignment", fixed = TRUE)
  testthat::expect_false(grepl("paste0(\"UNMAPPED_\"", script, fixed = TRUE))
})

testthat::test_that("MapThatProt and WGCNA source shared protein mapping utilities", {
  mapthatprot <- paste(readLines(testthat::test_path("..", "..", "02_id_mapping", "01_MapThatProt_batch.r"), warn = FALSE), collapse = "\n")
  wgcna <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(mapthatprot, "protein_mapping_utils.R", fixed = TRUE)
  testthat::expect_match(wgcna, "protein_mapping_utils.R", fixed = TRUE)
})

testthat::test_that("canonical protein-group mapping preserves rows and members", {
  testthat::skip_if_not_installed("dplyr")
  source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))

  entry_map <- data.frame(
    UNIPROT = c("Q9CQH5", "P12345", "Q8K1A0", "Q9D0M3", "P99999"),
    entry_full = c("FOO_MOUSE", "BAR_MOUSE", "BAZ_MOUSE", "ISO_MOUSE", "HUM_MOUSE"),
    entry_base = c("FOO", "BAR", "BAZ", "ISO", "HUM"),
    stringsAsFactors = FALSE
  )
  gene_map <- data.frame(
    input = c("FOO", "BAR", "BAZ", "ISO", "HUM", "Q9CQH5", "P12345", "Q8K1A0", "Q9D0M3"),
    primaryAccession = c("Q9CQH5", "P12345", "Q8K1A0", "Q9D0M3", "P99999", "Q9CQH5", "P12345", "Q8K1A0", "Q9D0M3"),
    stringsAsFactors = FALSE
  )
  accession_gene_map <- data.frame(
    UNIPROT = c("Q9CQH5", "P12345", "Q8K1A0", "Q9D0M3", "P99999"),
    gene_symbol = c("GENEA", "GENEA", "GENEB", "ISO", "HUM"),
    stringsAsFactors = FALSE
  )

  build <- function(ids, source_file = "contrast_a.csv", extra = list()) {
    df <- data.frame(gene_symbol = ids, log2fc = seq_along(ids), pval = 0.05, custom_stat = 10, stringsAsFactors = FALSE)
    for (nm in names(extra)) df[[nm]] <- extra[[nm]]
    build_canonical_protein_group_tables(
      df, "microglia", source_file, entry_map, gene_map, accession_gene_map,
      strict = TRUE
    )
  }

  single <- build("FOO_MOUSE")
  testthat::expect_equal(nrow(single$wide), 1L)
  testthat::expect_equal(single$wide$protein_group_ambiguity_class, "single_accession_single_gene")
  testthat::expect_true(single$wide$gene_level_claim_allowed)
  testthat::expect_true(single$wide$protein_level_claim_allowed)

  same_gene <- build("FOO_MOUSE; BAR_MOUSE")
  testthat::expect_equal(same_gene$wide$protein_group_ambiguity_class, "multi_accession_same_gene")
  testthat::expect_true(same_gene$wide$gene_level_claim_allowed)
  testthat::expect_false(same_gene$wide$protein_level_claim_allowed)

  multi_gene <- build("FOO_MOUSE;BAZ_MOUSE")
  testthat::expect_equal(multi_gene$wide$protein_group_ambiguity_class, "multi_gene_indistinguishable")
  testthat::expect_false(multi_gene$wide$gene_level_claim_allowed)

  isoform <- build("sp|Q9D0M3-2|ISO_MOUSE;Q9D0M3")
  testthat::expect_equal(isoform$wide$protein_group_ambiguity_class, "single_accession_single_gene")

  partial <- build("FOO_MOUSE;NOHIT_MOUSE")
  testthat::expect_equal(partial$wide$protein_group_ambiguity_class, "partially_mapped_group")
  testthat::expect_equal(partial$wide$n_unmapped_members, 1L)

  unresolved <- build("NOHIT_MOUSE")
  testthat::expect_equal(unresolved$wide$protein_group_ambiguity_class, "unresolved_group")
  testthat::expect_equal(unresolved$wide$mapping_status, "unmapped")

  mixed <- build("FOO_MOUSE;HUM_HUMAN")
  testthat::expect_equal(mixed$wide$protein_group_ambiguity_class, "mixed_species_or_contaminant")

  contaminant <- build("FOO_MOUSE;CON__KRT_MOUSE")
  testthat::expect_equal(contaminant$wide$protein_group_ambiguity_class, "mixed_species_or_contaminant")

  repeated <- build(" FOO_MOUSE ; FOO_MOUSE ")
  testthat::expect_equal(repeated$wide$n_members_original, 2L)
  testthat::expect_equal(repeated$wide$n_members_canonical, 1L)
  testthat::expect_equal(nrow(repeated$bridge), 2L)

  delimiter <- split_protein_group_members("FOO_MOUSE; BAR_MOUSE, BAZ_MOUSE / ISO_MOUSE")
  testthat::expect_equal(delimiter, c("FOO_MOUSE", "BAR_MOUSE", "BAZ_MOUSE", "ISO_MOUSE"))

  forward <- build("FOO_MOUSE;BAR_MOUSE")
  reverse <- build("BAR_MOUSE;FOO_MOUSE", source_file = "contrast_b.csv")
  testthat::expect_equal(forward$wide$protein_group_ambiguity_class, reverse$wide$protein_group_ambiguity_class)
  testthat::expect_equal(forward$wide$ProteinGroupID, reverse$wide$ProteinGroupID)
  testthat::expect_equal(forward$wide$member_identifiers_canonical, reverse$wide$member_identifiers_canonical)

  master_df <- data.frame(gene_symbol = "FOO_MOUSE;BAR_MOUSE", Master.Protein = "BAR_MOUSE", logFC = 1, stringsAsFactors = FALSE)
  master <- build_canonical_protein_group_tables(master_df, "microglia", "master.csv", entry_map, gene_map, accession_gene_map, strict = TRUE)
  testthat::expect_equal(master$wide$protein_group_ambiguity_class, "explicit_master_with_subordinate_members")
  testthat::expect_equal(master$wide$representative_selection_rule, "explicit_upstream_master")

  no_master <- build("FOO_MOUSE;BAR_MOUSE")
  testthat::expect_false(no_master$wide$protein_group_ambiguity_class == "explicit_master_with_subordinate_members")

  same_a <- build("FOO_MOUSE;BAR_MOUSE", source_file = "a.csv")
  same_b <- build("BAR_MOUSE;FOO_MOUSE", source_file = "b.csv")
  testthat::expect_equal(same_a$wide$ProteinGroupID, same_b$wide$ProteinGroupID)

  distinct_df <- data.frame(feature_id = c("feature_1", "feature_2"), gene_symbol = c("FOO_MOUSE;BAR_MOUSE", "BAR_MOUSE;FOO_MOUSE"), stringsAsFactors = FALSE)
  distinct <- build_canonical_protein_group_tables(distinct_df, "microglia", "distinct.csv", entry_map, gene_map, accession_gene_map, strict = TRUE)
  testthat::expect_equal(length(unique(distinct$wide$ProteinGroupID)), 2L)

  collision_df <- data.frame(gene_symbol = c("FOO_MOUSE;BAR_MOUSE", "BAR_MOUSE;FOO_MOUSE"), stringsAsFactors = FALSE)
  testthat::expect_error(
    build_canonical_protein_group_tables(collision_df, "microglia", "collision.csv", entry_map, gene_map, accession_gene_map, strict = TRUE),
    "Unstable ProteinGroupID"
  )

  stats <- build("FOO_MOUSE")
  testthat::expect_true("custom_stat" %in% names(stats$wide))
  testthat::expect_equal(nrow(stats$wide), 1L)
  testthat::expect_equal(nrow(stats$bridge), 1L)
  testthat::expect_false(any(is.na(stats$bridge$member_identifier_original)))
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
