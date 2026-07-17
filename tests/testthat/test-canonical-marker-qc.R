testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "dataset_inputs.R"))
source(testthat::test_path("..", "..", "R", "qc_exploration_utils.R"))

qc_test_context <- function() {
  list(
    entry_map = data.frame(
      entry_base = c("GENEA", "GENEB"), UNIPROT = c("P00001", "Q00001"),
      stringsAsFactors = FALSE
    ),
    gene_map = data.frame(
      input = c("GENEA", "GENEB", "P00001", "Q00001"),
      primaryAccession = c("P00001", "Q00001", "P00001", "Q00001"),
      stringsAsFactors = FALSE
    ),
    accession_gene_map = data.frame(
      UNIPROT = c("P00001", "Q00001"), gene_symbol = c("GENEA", "GENEB"),
      stringsAsFactors = FALSE
    ),
    reviewed_map = data.frame(UNIPROT = c("P00001", "Q00001"), reviewed_status = "reviewed"),
    manual_mapping = NULL,
    gene_annotation_maps = list(
      uniprot_map = data.frame(
        UNIPROT = c("P00001", "Q00001"), SYMBOL = c("GeneA", "GeneB"),
        ENTREZID = c("1", "2"), stringsAsFactors = FALSE
      ),
      symbol_map = data.frame(SYMBOL = c("GeneA", "GeneB"), ENTREZID = c("1", "2")),
      alias_map = data.frame(ALIAS = character(), SYMBOL = character(), ENTREZID = character()),
      orgdb_package_version = "test", annotation_contract_version = "mouse_gene_annotation_v1"
    ),
    manual_gene_annotation_overrides = NULL,
    idmap_file = NA_character_, idmap_sha256 = "test-idmap",
    manual_mapping_file = NA_character_, manual_gene_annotation_file = NA_character_
  )
}

write_qc_fixture <- function(df) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(df, path, row.names = FALSE, na = "")
  path
}

load_qc_fixture <- function(df) {
  qc_load_canonical_expression(
    write_qc_fixture(df), NA_character_, "neuron_soma",
    mapping_context = qc_test_context(), sample_columns = c("S1", "S2")
  )
}

testthat::test_that("canonical QC loader preserves every quantitative source row", {
  n <- 5538L
  source_rows <- data.frame(source = seq_len(n))
  feature_table <- data.frame(
    ProteinGroupID = sprintf("PG:neuron_soma:test%04d", seq_len(n)),
    source_row_id = seq_len(n), stringsAsFactors = FALSE
  )
  testthat::expect_true(qc_validate_source_row_preservation(source_rows, feature_table))
  testthat::expect_error(
    qc_validate_source_row_preservation(source_rows, feature_table[-n, , drop = FALSE]),
    "lost quantitative source rows"
  )
})

testthat::test_that("blank display identifiers fall back to stable Protein.Group member provenance", {
  x <- load_qc_fixture(data.frame(
    Protein.Group = "P00001", `T: Protein.Names` = "", S1 = 1, S2 = 2,
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  testthat::expect_equal(nrow(x$mat), 1L)
  testthat::expect_equal(x$feature_table$member_identifier_source, "Protein.Group")
  testthat::expect_equal(x$member_bridge$member_accession, "P00001")
  testthat::expect_true(x$feature_table$wgcna_eligible)
  testthat::expect_false(x$feature_table$marker_eligible)
  testthat::expect_equal(
    x$feature_table$marker_exclusion_reason,
    "blank_legacy_display_identifier_not_gene_claim_eligible"
  )
})

testthat::test_that("multi-member groups retain all members and order-independent identity", {
  a <- load_qc_fixture(data.frame(
    Protein.Group = "P00001;Q00001", `T: Protein.Names` = "GENEA_MOUSE;GENEB_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  b <- load_qc_fixture(data.frame(
    Protein.Group = "Q00001;P00001", `T: Protein.Names` = "GENEB_MOUSE;GENEA_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  testthat::expect_equal(a$feature_table$ProteinGroupID, b$feature_table$ProteinGroupID)
  testthat::expect_equal(nrow(a$member_bridge), 2L)
  testthat::expect_setequal(a$member_bridge$member_accession, c("P00001", "Q00001"))
  testthat::expect_false(a$feature_table$representative_accession == a$feature_table$ProteinGroupID)
})

testthat::test_that("source Protein.Group contaminant evidence is central and explicit", {
  x <- load_qc_fixture(data.frame(
    Protein.Group = "Cont_P00001", `T: Protein.Names` = "GENEA_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  testthat::expect_equal(x$feature_table$protein_group_ambiguity_class, "mixed_species_or_contaminant")
  testthat::expect_match(x$feature_table$contaminant_source_evidence, "Protein.Group=Cont_P00001", fixed = TRUE)
  testthat::expect_equal(x$feature_table$wgcna_exclusion_reason, "contaminant_source_provenance")
})

testthat::test_that("known non-mouse source labels are caught by the general provenance rule", {
  labels <- c("TRYP_PIG", "ALBU_BOVIN", "APOA1_BOVIN", "API_ACHLY", "Q1WEI2_PSEAI", "TRFE_BOVIN")
  df <- data.frame(
    ProteinGroupID = paste0("PG:neuron_soma:", seq_along(labels)),
    Protein.Group = paste0("Cont_", seq_along(labels)),
    `T: Protein.Names` = labels,
    S1 = seq_along(labels), S2 = seq_along(labels) + 1,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  x <- load_qc_fixture(df)
  testthat::expect_true(all(!x$feature_table$wgcna_eligible))
  testthat::expect_true(all(x$feature_table$wgcna_exclusion_reason == "contaminant_source_provenance"))
  testthat::expect_true(all(nzchar(x$feature_table$contaminant_source_evidence)))
})

testthat::test_that("mixed-species and partially mapped mouse groups remain distinct", {
  mixed <- load_qc_fixture(data.frame(
    Protein.Group = "P00001;H00001", `T: Protein.Names` = "GENEA_MOUSE;OTHER_HUMAN",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  partial <- load_qc_fixture(data.frame(
    Protein.Group = "P00001;MISSING", `T: Protein.Names` = "GENEA_MOUSE;MISSING_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  testthat::expect_equal(mixed$feature_table$protein_group_ambiguity_class, "mixed_species_or_contaminant")
  testthat::expect_equal(partial$feature_table$protein_group_ambiguity_class, "partially_mapped_group")
  testthat::expect_true(partial$feature_table$wgcna_eligible)
})

testthat::test_that("marker matching deduplicates panel by ProteinGroupID and flags conflicts", {
  x <- load_qc_fixture(data.frame(
    Protein.Group = "P00001;Q00001", `T: Protein.Names` = "GENEA_MOUSE;GENEB_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  same_panel <- data.frame(
    marker_panel = c("panel_a", "panel_a"), marker_gene = c("GeneA", "GeneB"),
    marker_class = "soma", stringsAsFactors = FALSE
  )
  one <- qc_match_markers_to_protein_groups(
    same_panel, x$member_bridge, x$feature_table,
    panel_col = "marker_panel", gene_col = "marker_gene", class_col = "marker_class"
  )
  testthat::expect_equal(nrow(one$matches), 1L)
  testthat::expect_equal(one$matches$n_matched_official_genes, 2L)

  conflicting <- rbind(
    same_panel[1, ],
    data.frame(marker_panel = "panel_b", marker_gene = "GeneB", marker_class = "neuropil")
  )
  two <- qc_match_markers_to_protein_groups(
    conflicting, x$member_bridge, x$feature_table,
    panel_col = "marker_panel", gene_col = "marker_gene", class_col = "marker_class"
  )
  testthat::expect_true(two$conflicts$conflicting_marker_panels)
  testthat::expect_true(two$conflicts$conflicting_marker_classes)
  testthat::expect_false(any(two$matches$primary_score_eligible))
})

testthat::test_that("alignment, collisions and ordered fingerprints fail deterministically", {
  x <- load_qc_fixture(data.frame(
    Protein.Group = c("P00001", "Q00001"), `T: Protein.Names` = c("GENEA_MOUSE", "GENEB_MOUSE"),
    S1 = c(1, 2), S2 = c(2, 3), check.names = FALSE
  ))
  testthat::expect_error(
    qc_validate_canonical_expression(x$mat[2:1, , drop = FALSE], x$feature_table),
    "not aligned"
  )
  duplicate <- x$feature_table
  duplicate$ProteinGroupID[[2]] <- duplicate$ProteinGroupID[[1]]
  testthat::expect_error(qc_validate_canonical_expression(x$mat, duplicate), "duplicate")
  testthat::expect_identical(
    wgcna_feature_key_fingerprint(x$feature_table$ProteinGroupID),
    wgcna_feature_key_fingerprint(x$feature_table$ProteinGroupID)
  )
  testthat::expect_false(identical(
    wgcna_feature_key_fingerprint(x$feature_table$ProteinGroupID),
    wgcna_feature_key_fingerprint(rev(x$feature_table$ProteinGroupID))
  ))
})

testthat::test_that("legacy WGCNA states are rejected by canonical QC bridge validation", {
  x <- load_qc_fixture(data.frame(
    Protein.Group = "P00001", `T: Protein.Names` = "GENEA_MOUSE",
    S1 = 1, S2 = 2, check.names = FALSE
  ))
  legacy <- tempfile(fileext = ".rds")
  saveRDS(list(expression.data = data.frame(P00001 = c(1, 2))), legacy)
  audit <- qc_validate_optional_wgcna_bridge(legacy, x)
  testthat::expect_equal(audit$validation_status, "rejected")
  testthat::expect_match(audit$reason, "Legacy WGCNA cached state")
})

testthat::test_that("active gene-aware consumers do not match marker genes to matrix row names", {
  scripts <- c(
    "03_qc_exploration/04_marker_rank_abundance_qc.r",
    "03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r",
    "03_qc_exploration/05_empirical_roi_marker_discovery.r",
    "03_qc_exploration/07_wgcna_marker_trait_export.r"
  )
  text <- vapply(scripts, function(path) {
    paste(readLines(testthat::test_path("..", "..", path), warn = FALSE), collapse = "\n")
  }, character(1))
  testthat::expect_true(all(grepl("qc_load_canonical_expression", text, fixed = TRUE)))
  testthat::expect_false(any(grepl("qc_read_expression", text, fixed = TRUE)))
  testthat::expect_false(any(grepl("normalize_gene_token(rownames(mat))", text, fixed = TRUE)))
  testthat::expect_false(any(grepl("make.unique", text, fixed = TRUE)))
})

testthat::test_that("pipeline declares Stage 01 canonical marker-QC inputs", {
  pipeline <- paste(readLines(testthat::test_path("..", "..", "pipeline.yml"), warn = FALSE), collapse = "\n")
  for (script in c(
    "04_marker_rank_abundance_qc.r", "04c_marker_detectability_and_wgcna_bridge.r",
    "05_empirical_roi_marker_discovery.r", "07_wgcna_marker_trait_export.r"
  )) testthat::expect_match(pipeline, script, fixed = TRUE)
  testthat::expect_match(pipeline, "pgmatrix_imputed_<dataset>", fixed = TRUE)
  testthat::expect_match(pipeline, "MOUSE_10090_idmapping.dat", fixed = TRUE)
  testthat::expect_match(pipeline, "wgcna_final_model_state.rds", fixed = TRUE)
})
