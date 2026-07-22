testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "compartment_abundance_utils.R"))

ca_test_metadata <- function(groups = c(1, 2, 3), include_group = TRUE) {
  out <- data.frame(
    Sample = paste0("S", seq_along(groups)),
    AnimalID = paste0("A", seq_along(groups)),
    ExpGroup = groups,
    region = "CA1", Region = "WRONG",
    layer = "sp", Layer = "WRONG",
    ReplicateGroup = "Left", plate = "P1", dataset = "neuron_soma",
    stringsAsFactors = FALSE
  )
  if (include_group) out$group <- "neuron_soma"
  out
}

testthat::test_that("experimental-group filtering is exact, case-normalized, and aligned", {
  meta_numeric <- ca_test_metadata(c(1, 2, 1))
  mat <- matrix(seq_len(6), nrow = 2, dimnames = list(c("P2", "P1"), meta_numeric$Sample))
  numeric <- ca_filter_expression_group(mat, meta_numeric, "con")
  meta_text <- meta_numeric
  meta_text$ExpGroup <- c("CON", "RES", "CON")
  text <- ca_filter_expression_group(mat, meta_text, "CON")
  testthat::expect_identical(numeric$resolved_group_column, "ExpGroup")
  testthat::expect_identical(colnames(numeric$mat), c("S1", "S3"))
  testthat::expect_identical(numeric$mat, text$mat)
  testthat::expect_identical(as.character(numeric$meta$Sample), colnames(numeric$mat))
})

testthat::test_that("ExpGroup wins over anatomical group", {
  meta <- ca_test_metadata(c(1, 2, 3))
  meta$group <- c("CON", "CON", "CON")
  filtered <- ca_filter_expression_group(matrix(1:3, nrow = 1, dimnames = list("P1", meta$Sample)), meta, "CON")
  testthat::expect_identical(filtered$resolved_group_column, "ExpGroup")
  testthat::expect_identical(colnames(filtered$mat), "S1")
})

testthat::test_that("group filtering fails for a missing column or requested level", {
  meta <- ca_test_metadata(c(1, 2, 3))
  meta$ExpGroup <- NULL
  meta$group <- "neuron_soma"
  mat <- matrix(1:3, nrow = 1, dimnames = list("P1", meta$Sample))
  testthat::expect_error(ca_filter_expression_group(mat, meta, "CON"), "No valid experimental-group")
  meta$ExpGroup <- c(1, 2, 2)
  testthat::expect_error(ca_filter_expression_group(mat, meta, "SUS"), "is absent")
})

testthat::test_that("unfiltered default is unchanged and does not require a group column", {
  meta <- ca_test_metadata(c(1, 2, 3), include_group = FALSE)
  meta$ExpGroup <- NULL
  meta$AnimalID <- NULL
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("P1", "P2"), meta$Sample))
  result <- ca_filter_expression_group(mat, meta, "")
  testthat::expect_identical(result$mat, mat)
  testthat::expect_identical(result$meta, meta)
  testthat::expect_true(is.na(result$resolved_group_column))
})

testthat::test_that("matrix and metadata misalignment fails clearly", {
  meta <- ca_test_metadata(c(1, 1, 1))
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("S1", "missing")))
  testthat::expect_error(ca_filter_expression_group(mat, meta, "CON"), "not exactly aligned")
})

testthat::test_that("corrected lowercase anatomical fields override aliases", {
  cols <- ca_resolve_metadata_columns(ca_test_metadata(), require_dataset = TRUE)
  testthat::expect_identical(cols$region, "region")
  testthat::expect_identical(cols$layer, "layer")
})

ca_aggregation_fixture <- function() {
  meta <- expand.grid(
    AnimalID = c("A1", "A2", "A3"), region = c("CA1", "CA2", "CA3", "DG"),
    layer = "sp", ReplicateGroup = c("Left", "Right"), stringsAsFactors = FALSE
  )
  meta$Sample <- paste0("S", seq_len(nrow(meta)))
  meta$dataset <- "neuron_soma"
  meta$ExpGroup <- 1
  meta$plate <- "P1"
  meta <- meta[, c("Sample", "AnimalID", "ExpGroup", "region", "layer", "ReplicateGroup", "plate", "dataset")]
  mat <- matrix(NA_real_, 2, nrow(meta), dimnames = list(c("P1", "P2"), meta$Sample))
  mat["P1", ] <- ifelse(meta$ReplicateGroup == "Left", 1, 3)
  mat["P2", ] <- seq_len(nrow(meta))
  list(mat = mat, meta = meta)
}

testthat::test_that("biological aggregation gives samples, strata, regions, and animals equal weight", {
  x <- ca_aggregation_fixture()
  out <- ca_aggregate_abundance(x$mat, x$meta, "mean", min_regions = 3)
  testthat::expect_true(all(out$stratum$values["P1", ] == 2))
  testthat::expect_true(all(out$animal$values["P1", ] == 2))
  testthat::expect_equal(out$dataset$values["P1", 1], 2)

  # More CA1/CA2 layers must not increase their regional weight.
  meta <- x$meta[x$meta$AnimalID == "A1", ]
  extra <- meta[meta$region %in% c("CA1", "CA2"), ]
  extra$layer <- "sr"
  extra$Sample <- paste0(extra$Sample, "x")
  meta <- rbind(meta, extra)
  values <- ifelse(meta$region %in% c("CA1", "CA2"), 10, 0)
  mat <- matrix(values, nrow = 1, dimnames = list("P1", meta$Sample))
  balanced <- ca_aggregate_abundance(mat, meta, "mean", min_regions = 3)
  testthat::expect_equal(balanced$animal$values["P1", 1], 5)
  testthat::expect_false(isTRUE(all.equal(mean(values), balanced$animal$values["P1", 1])))

  # Unequal technical sample counts cannot overweight one animal.
  meta2 <- x$meta
  values2 <- ifelse(meta2$AnimalID == "A1", 9, 0)
  extra2 <- meta2[meta2$AnimalID == "A1", ]
  extra2$Sample <- paste0(extra2$Sample, "extra")
  meta2 <- rbind(meta2, extra2)
  values2 <- c(values2, rep(9, nrow(extra2)))
  mat2 <- matrix(values2, nrow = 1, dimnames = list("P1", meta2$Sample))
  animal_balanced <- ca_aggregate_abundance(mat2, meta2, "mean", min_regions = 3)
  testthat::expect_equal(animal_balanced$dataset$values["P1", 1], 3)
})

testthat::test_that("mean primary and median hierarchy are independent", {
  x <- ca_aggregation_fixture()
  x$mat["P1", x$meta$AnimalID == "A1" & x$meta$region == "CA1" & x$meta$ReplicateGroup == "Right"] <- 101
  mean_out <- ca_aggregate_abundance(x$mat, x$meta, "mean", min_regions = 3)
  median_out <- ca_aggregate_abundance(x$mat, x$meta, "median", min_regions = 3)
  testthat::expect_gt(mean_out$dataset$values["P1", 1], median_out$dataset$values["P1", 1])
})

testthat::test_that("animal validity and two-of-three versus strict universes use integer counts", {
  x <- ca_aggregation_fixture()
  x$mat["P1", x$meta$AnimalID == "A3"] <- NA_real_
  x$mat["P2", x$meta$AnimalID == "A3" & x$meta$region %in% c("CA3", "DG")] <- NA_real_
  out <- ca_aggregate_abundance(x$mat, x$meta, "mean", min_regions = 3)
  audit <- ca_detection_audit(out, min_animals = 2, strict_animals = 3)
  p1 <- audit[audit$ProteinGroupID == "P1", ]
  p2 <- audit[audit$ProteinGroupID == "P2", ]
  testthat::expect_true(p1$primary_eligible)
  testthat::expect_false(p1$strict_eligible)
  testthat::expect_equal(p1$valid_animal_count, 2)
  testthat::expect_equal(p2$valid_animal_count, 2)
  testthat::expect_true(is.na(out$animal$values["P2", out$animal$group_meta$AnimalID == "A3"]))
})

testthat::test_that("raw detection is not replaced by imputed availability", {
  x <- ca_aggregation_fixture()
  x$mat["P1", x$meta$AnimalID == "A3"] <- NA_real_
  raw_audit <- ca_detection_audit(ca_aggregate_abundance(x$mat, x$meta, "mean", 3), 2, 3)
  imputed <- x$mat
  imputed[is.na(imputed)] <- 0
  post_audit <- ca_detection_audit(ca_aggregate_abundance(imputed, x$meta, "mean", 3), 2, 3)
  testthat::expect_equal(raw_audit$valid_animal_count[raw_audit$ProteinGroupID == "P1"], 2)
  testthat::expect_equal(post_audit$valid_animal_count[post_audit$ProteinGroupID == "P1"], 3)
})

testthat::test_that("ranking is deterministic with bounded percentiles and unique base rows", {
  estimates <- data.frame(dataset = "d", ProteinGroupID = c("P2", "P1", "P3"), region_balanced_log2 = c(5, 5, 4))
  audit <- data.frame(dataset = "d", ProteinGroupID = c("P1", "P2", "P3"), primary_eligible = TRUE)
  ranks <- ca_rank_table(estimates, audit)
  testthat::expect_identical(ranks$ProteinGroupID, c("P1", "P2", "P3"))
  testthat::expect_true(all(ranks$RankPercentile >= 0 & ranks$RankPercentile <= 1))
  testthat::expect_false(anyDuplicated(ranks[, c("dataset", "RankGroup", "ProteinGroupID")]) > 0)
  one <- ca_rank_table(estimates[1, ], audit[audit$ProteinGroupID == "P2", ])
  testthat::expect_equal(one$RankPercentile, 1)
})

testthat::test_that("marker annotations cannot duplicate the grey rank curve", {
  ranks <- data.frame(dataset = "d", RankGroup = "g", ProteinGroupID = c("P1", "P2"), Rank = 1:2)
  markers <- data.frame(
    ProteinGroupID = c("P1", "P1"), marker_class = "Soma markers",
    marker_panel = c("a", "b"), marker_gene = "Gene1"
  )
  joined <- ca_join_marker_annotations(ranks, markers)
  testthat::expect_equal(nrow(joined), nrow(ranks))
  testthat::expect_identical(joined$marker_panels[[1]], "a;b")
})

testthat::test_that("robust z handles zero MAD explicitly", {
  tbl <- expand.grid(ProteinGroupID = c("P1", "P2"), dataset = "d", AnimalID = c("A1", "A2", "A3"))
  tbl$animal_log2_abundance <- ifelse(tbl$ProteinGroupID == "P1", 1, match(tbl$AnimalID, c("A1", "A2", "A3")))
  z <- ca_robust_standardize(tbl)
  testthat::expect_true(all(z$zero_MAD_excluded[z$ProteinGroupID == "P1"]))
  testthat::expect_true(all(is.na(z$robust_z[z$ProteinGroupID == "P1"])))
})

testthat::test_that("marker scoring deduplicates proteins, excludes conflicts, and reports precision counts", {
  robust <- expand.grid(ProteinGroupID = c("P1", "P2"), dataset = "d", AnimalID = c("A1", "A2"))
  robust$animal_log2_abundance <- 1
  robust$robust_z <- c(1, 2, 2, 3)
  markers <- data.frame(
    ProteinGroupID = c("P1", "P1", "P2", "PX"),
    marker_class = c("Soma markers", "Soma markers", "Soma markers", "Neuropil markers"),
    marker_gene = c("G1", "G1", "G2", "GX"), marker_panel = c("a", "b", "a", "x"),
    fidelity_subpanel = c("nuclear", "nuclear", "nuclear", "synaptic"),
    primary_score_eligible = c(TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE
  )
  scored <- ca_score_markers(robust, markers, min_proteins = 1, min_fraction = 0)
  soma <- scored$scores[scored$scores$marker_class == "Soma markers", ]
  neuropil <- scored$scores[scored$scores$marker_class == "Neuropil markers", ]
  testthat::expect_true(all(soma$n_eligible_proteins == 1))
  testthat::expect_true(all(soma$n_contributing_proteins == 1))
  testthat::expect_true(all(neuropil$n_contributing_proteins == 0))
  testthat::expect_true(all(grepl("unequal_protein_counts", soma$precision_interpretation)))
  testthat::expect_false("P2" %in% scored$marker_map$ProteinGroupID)
})

testthat::test_that("canonical fidelity matching excludes cross-class conflicts", {
  registry <- data.frame(
    marker_set = c("soma", "neuropil"),
    fidelity_marker_class = c("Soma markers", "Neuropil markers"),
    fidelity_subpanel = c("nuclear", "synaptic"),
    gene_symbol = c("GeneA", "GeneA"), include_in_fidelity_score = TRUE,
    stringsAsFactors = FALSE
  )
  feature <- data.frame(
    ProteinGroupID = "P1", member_gene_symbols = "GeneA", official_gene_symbol = "GeneA",
    protein_group_ambiguity_class = "single_accession_single_gene",
    gene_level_claim_allowed = TRUE, original_identifier = "GENEA_MOUSE",
    stringsAsFactors = FALSE
  )
  matched <- ca_match_fidelity_markers(registry, feature, "P1")
  testthat::expect_true(all(matched$conflicting_marker_classes))
  testthat::expect_false(any(matched$primary_score_eligible))
  testthat::expect_true(all(matched$exclusion_reason == "conflicting_fidelity_marker_classes"))
})

testthat::test_that("04 keeps the all-group default and declares the CON namespace contract", {
  script <- paste(readLines(testthat::test_path("..", "..", "03_qc_exploration", "04_marker_rank_abundance_qc.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "PROTEOMICS_RANK_ABUNDANCE_GROUP_FILTER", fixed = TRUE)
  testthat::expect_match(script, "--group", fixed = TRUE)
  testthat::expect_match(script, "ca_namespace_paths", fixed = TRUE)
  testthat::expect_match(script, "default = Sys.getenv", fixed = TRUE)
  testthat::expect_match(script, "else \"\"", fixed = TRUE)
})

testthat::test_that("paired bootstrap and leave-one-out cases are exact and deterministic", {
  first <- ca_exact_paired_bootstrap_draws(c("A3", "A1", "A2"))
  second <- ca_exact_paired_bootstrap_draws(c("A2", "A3", "A1"))
  testthat::expect_identical(first, second)
  testthat::expect_equal(nrow(first), 27)
  testthat::expect_equal(nrow(unique(first[c("draw_1", "draw_2", "draw_3")])), 27)
  testthat::expect_equal(nrow(ca_leave_one_animal_out_cases(c("A1", "A2", "A3"))), 3)
})
