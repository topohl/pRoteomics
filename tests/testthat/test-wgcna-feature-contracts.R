testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))
source(testthat::test_path("..", "..", "R", "module_contracts.R"))

wgcna_test_maps <- function() {
  list(
    entry_map = data.frame(
      entry_base = c("AIF1", "FOO1", "FOO2", "BAR1", "HUM1"),
      UNIPROT = c("P55008", "P00001", "P00002", "Q00001", "H00001"),
      stringsAsFactors = FALSE
    ),
    gene_map = data.frame(
      input = c("AIF1", "FOO", "FOO", "BAR", "HUM"),
      primaryAccession = c("P55008", "P00001", "P00002", "Q00001", "H00001"),
      stringsAsFactors = FALSE
    ),
    accession_gene_map = data.frame(
      UNIPROT = c("P55008", "P00001", "P00002", "Q00001", "H00001"),
      gene_symbol = c("AIF1", "FOO", "FOO", "BAR", "HUM"),
      stringsAsFactors = FALSE
    )
  )
}

build_wgcna_test <- function(ids, pgids = paste0("PG:microglia:", seq_along(ids))) {
  maps <- wgcna_test_maps()
  df <- data.frame(
    ProteinGroupID = pgids,
    original_identifier = ids,
    S1 = seq_along(ids),
    S2 = seq_along(ids) + 0.5,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  build_wgcna_canonical_features(
    df, "microglia", "synthetic.xlsx",
    maps$entry_map, maps$gene_map, maps$accession_gene_map,
    sample_columns = c("S1", "S2")
  )
}

testthat::test_that("canonical WGCNA construction covers ambiguity classes without row expansion", {
  x <- build_wgcna_test(c(
    "AIF1_MOUSE",
    "FOO1_MOUSE;FOO2_MOUSE",
    "FOO1_MOUSE;BAR1_MOUSE",
    "FOO1_MOUSE;MISSING_MOUSE",
    "UNKNOWN_MOUSE"
  ))
  testthat::expect_equal(nrow(x$feature_table), 5L)
  testthat::expect_equal(ncol(x$expression_data), 5L)
  testthat::expect_equal(x$feature_table$protein_group_ambiguity_class, c(
    "single_accession_single_gene", "multi_accession_same_gene",
    "multi_gene_indistinguishable", "partially_mapped_group", "unresolved_group"
  ))
  testthat::expect_true(all(x$feature_table$included_in_wgcna))
  testthat::expect_equal(colnames(x$expression_data), x$feature_table$ProteinGroupID)
  testthat::expect_equal(nrow(x$member_bridge), 8L)
})

testthat::test_that("mixed species and contaminants are excluded with reasons", {
  x <- build_wgcna_test(c("AIF1_MOUSE;HUM1_HUMAN", "CON__KERATIN;AIF1_MOUSE"))
  testthat::expect_true(all(x$feature_table$protein_group_ambiguity_class == "mixed_species_or_contaminant"))
  testthat::expect_false(any(x$feature_table$included_in_wgcna))
  testthat::expect_equal(ncol(x$expression_data), 0L)
  testthat::expect_true(all(x$feature_table$wgcna_exclusion_reason == "mixed_species_or_contaminant"))
})

testthat::test_that("representative accessions and genes never collapse distinct ProteinGroupIDs", {
  x <- build_wgcna_test(
    c("AIF1_MOUSE", "AIF1_MOUSE"),
    c("PG:microglia:feature_a", "PG:microglia:feature_b")
  )
  testthat::expect_equal(ncol(x$expression_data), 2L)
  testthat::expect_equal(length(unique(x$feature_table$representative_accession)), 1L)
  testthat::expect_equal(length(unique(x$feature_table$representative_gene_symbol)), 1L)
  testthat::expect_equal(length(unique(colnames(x$expression_data))), 2L)
})

testthat::test_that("member order changes neither canonical class nor explicit ProteinGroupID", {
  a <- build_wgcna_test("FOO1_MOUSE;FOO2_MOUSE", "PG:microglia:stable")
  b <- build_wgcna_test("FOO2_MOUSE;FOO1_MOUSE", "PG:microglia:stable")
  testthat::expect_equal(a$feature_table$member_identifiers_canonical, b$feature_table$member_identifiers_canonical)
  testthat::expect_equal(a$feature_table$protein_group_ambiguity_class, b$feature_table$protein_group_ambiguity_class)
  testthat::expect_equal(colnames(a$expression_data), colnames(b$expression_data))
})

testthat::test_that("stable upstream Protein.Group bridges member order and distinct source rows", {
  maps <- wgcna_test_maps()
  make_input <- function(members, source_id) data.frame(
    Protein.Group = source_id,
    `T: Protein.Names` = members,
    S1 = 1, S2 = 2,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  a <- build_wgcna_canonical_features(
    make_input("FOO1_MOUSE;FOO2_MOUSE", "P00001;P00002"),
    "microglia", "a.xlsx", maps$entry_map, maps$gene_map, maps$accession_gene_map,
    sample_columns = c("S1", "S2")
  )
  b <- build_wgcna_canonical_features(
    make_input("FOO2_MOUSE;FOO1_MOUSE", "P00002;P00001"),
    "microglia", "b.xlsx", maps$entry_map, maps$gene_map, maps$accession_gene_map,
    sample_columns = c("S1", "S2")
  )
  distinct <- build_wgcna_canonical_features(
    make_input("FOO1_MOUSE;FOO2_MOUSE", "P00001;P00002-2"),
    "microglia", "c.xlsx", maps$entry_map, maps$gene_map, maps$accession_gene_map,
    sample_columns = c("S1", "S2")
  )
  testthat::expect_equal(a$feature_table$ProteinGroupID, b$feature_table$ProteinGroupID)
  testthat::expect_equal(a$feature_table$protein_group_ambiguity_class, b$feature_table$protein_group_ambiguity_class)
  testthat::expect_false(identical(a$feature_table$ProteinGroupID, distinct$feature_table$ProteinGroupID))
})

testthat::test_that("entry-name contrasts and accession expression groups share ProteinGroupID", {
  maps <- wgcna_test_maps()
  contrast <- build_canonical_protein_group_tables(
    data.frame(gene_symbol = "FOO1_MOUSE;FOO2_MOUSE", log2fc = 1),
    "microglia", "contrast.csv", maps$entry_map, maps$gene_map, maps$accession_gene_map
  )
  expression <- build_wgcna_canonical_features(
    data.frame(
      Protein.Group = "P00002;P00001",
      `T: Protein.Names` = "FOO2_MOUSE;FOO1_MOUSE",
      S1 = 1, S2 = 2,
      check.names = FALSE
    ),
    "microglia", "expression.xlsx", maps$entry_map, maps$gene_map, maps$accession_gene_map,
    sample_columns = c("S1", "S2")
  )
  testthat::expect_equal(contrast$wide$ProteinGroupID, expression$feature_table$ProteinGroupID)
})

testthat::test_that("missing and duplicate ProteinGroupIDs fail instead of being repaired", {
  testthat::expect_error(build_wgcna_test(c("AIF1_MOUSE", "BAR1_MOUSE"), c("PG:x", "PG:x")), "collision|Duplicate")
  testthat::expect_error(build_wgcna_test("AIF1_MOUSE", NA_character_), "Missing ProteinGroupID")
})

testthat::test_that("alignment and complete module provenance are enforced", {
  x <- build_wgcna_test(c("AIF1_MOUSE", "FOO1_MOUSE;BAR1_MOUSE"))
  testthat::expect_silent(validate_wgcna_feature_alignment(x$expression_data, x$feature_table))
  bad <- x$expression_data[, rev(seq_len(ncol(x$expression_data))), drop = FALSE]
  testthat::expect_error(validate_wgcna_feature_alignment(bad, x$feature_table), "misaligned")

  modules <- transform(
    x$feature_table,
    dataset = "microglia", ModuleSet = "WGCNA", ModuleID = "WGCNA_m01",
    ModuleLegacyID = "WGCNA_#486A8A", ModuleColor = "#486A8A",
    ModuleColorName = "deep_steel_blue", ModuleColorLabel = "deep steel blue",
    ProteinID = ProteinGroupID, kME = 0.8
  )
  testthat::expect_silent(validate_wgcna_module_definitions(modules))
  testthat::expect_equal(modules$member_accessions[[2]], "P00001;Q00001")
  modules$ProteinID[[1]] <- "P55008"
  testthat::expect_error(validate_wgcna_module_definitions(modules), "deprecated alias")
})

testthat::test_that("ambiguous groups remain quantitative but are blocked from gene annotation", {
  x <- build_wgcna_test(c("AIF1_MOUSE", "FOO1_MOUSE;BAR1_MOUSE", "UNKNOWN_MOUSE"))
  testthat::expect_equal(ncol(x$expression_data), 3L)
  eligible <- x$feature_table$ProteinGroupID[x$feature_table$gene_level_claim_allowed]
  testthat::expect_equal(eligible, x$feature_table$ProteinGroupID[[1]])
  testthat::expect_true(all(nzchar(x$feature_table$FeatureDisplayLabel)))
})

testthat::test_that("Phase 1B cache contract rejects legacy and mismatched states", {
  x <- build_wgcna_test("AIF1_MOUSE")
  testthat::expect_error(validate_wgcna_cached_state(list(expression.data = x$expression_data)), "Full WGCNA rerun")
  state <- list(
    feature_key_contract_version = wgcna_feature_key_contract_version(),
    feature_key_fingerprint = wgcna_feature_key_fingerprint(colnames(x$expression_data)),
    expression.data = x$expression_data
  )
  testthat::expect_silent(validate_wgcna_cached_state(state))
  testthat::expect_error(validate_wgcna_cached_state(state, "PG:other"), "differ")
})

testthat::test_that("active WGCNA feature construction cannot truncate or repair biological IDs", {
  script <- readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE)
  active <- paste(script, collapse = "\n")
  testthat::expect_false(grepl("fix_feature_ids", active, fixed = TRUE))
  testthat::expect_false(grepl("sub(\";.*$\"", active, fixed = TRUE))
  testthat::expect_false(grepl("make.unique(", active, fixed = TRUE))
  testthat::expect_match(active, "expression.data <- wgcna_canonical$expression_data", fixed = TRUE)
})
