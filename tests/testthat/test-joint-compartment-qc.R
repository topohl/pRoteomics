testthat::test_that("joint metadata requires an exact one-to-one matrix alignment", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  meta <- data.frame(sample_id = c("s1", "s2", "s3"), celltype_layer = c("neuropil", "soma", "microglia"), stringsAsFactors = FALSE)
  aligned <- joint_qc_prepare_metadata(meta, c("s1", "s2", "s3"))
  testthat::expect_equal(aligned$dataset, c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_error(joint_qc_prepare_metadata(rbind(meta, meta[1, ]), c("s1", "s2", "s3")), "Duplicated metadata")
  testthat::expect_error(joint_qc_prepare_metadata(meta, c("s1", "s2", "missing")), "one-to-one")
})

testthat::test_that("balanced feature filtering protects unequal dataset sizes", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  # Two neuropil samples and five soma samples: PG1 is abundant only in soma,
  # so a global threshold would retain it but the balanced rule rejects it.
  observed <- rbind(PG1 = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE), PG2 = rep(TRUE, 7L), PG3 = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  meta <- data.frame(dataset = c(rep("neuron_neuropil", 2), rep("neuron_soma", 5)), stringsAsFactors = FALSE)
  out <- joint_qc_filter_universes(observed, meta, rownames(observed), meta$dataset, min_block = .70, min_union = .30)
  testthat::expect_false(out$primary[["PG1"]])
  testthat::expect_true(out$primary[["PG2"]])
  testthat::expect_true(out$broad[["PG1"]])
  testthat::expect_true(all(c("primary_shared_core_inclusion", "complete_case_inclusion", "broad_union_inclusion") %in% names(out$audit)))
})

testthat::test_that("normalization and imputation are deterministic and label-blind", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  x <- matrix(c(10, 12, NA, 11, 20, 22), nrow = 2, dimnames = list(c("PG1", "PG2"), c("a", "b", "c")))
  norm_a <- joint_qc_joint_median_normalize(x)
  norm_b <- joint_qc_joint_median_normalize(x[, c("c", "a", "b")])
  testthat::expect_equal(norm_a$matrix[, c("c", "a", "b")], norm_b$matrix)
  imp_a <- joint_qc_median_impute(norm_a$matrix)
  imp_b <- joint_qc_median_impute(norm_a$matrix)
  testthat::expect_equal(imp_a$matrix, imp_b$matrix)
  testthat::expect_false(anyNA(imp_a$matrix))
  testthat::expect_true(all(is.finite(imp_a$matrix)))
})

testthat::test_that("canonical ProteinGroupID names propagate from source rows through positive filtering", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  # The source matrix begins with noncanonical row labels; canonical IDs are
  # attached once, in source-row order, before positive-value filtering.
  source_mat <- matrix(c(8, 16, 0, 4, 32, 64), nrow = 3, dimnames = list(c("source_11", "source_12", "source_13"), c("s1", "s2")))
  feature_table <- data.frame(ProteinGroupID = c("PG:1", "PG:2", "PG:3"), source_row_id = 1:3, stringsAsFactors = FALSE)
  rownames(source_mat) <- feature_table$ProteinGroupID
  joint_validate_canonical_feature_matrix(feature_table, source_mat, "synthetic raw")
  positive <- joint_qc_log2_positive(source_mat)
  testthat::expect_identical(rownames(positive$raw), feature_table$ProteinGroupID)
  testthat::expect_identical(rownames(positive$log2), feature_table$ProteinGroupID)
  testthat::expect_equal(nrow(positive$raw), nrow(source_mat))
  testthat::expect_true(is.na(positive$log2["PG:3", "s1"]))
})

testthat::test_that("strict feature subsets preserve requested order and reject bad identities", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  mat <- matrix(1:8, nrow = 4, dimnames = list(c("PG:1", "PG:2", "PG:3", "PG:4"), c("s1", "s2")))
  selected <- joint_validate_feature_subset(c("PG:3", "PG:1"), mat, "ordered synthetic subset")
  testthat::expect_identical(rownames(selected), c("PG:3", "PG:1"))
  testthat::expect_equal(unname(selected[, "s1"]), c(3L, 1L))
  testthat::expect_error(joint_validate_feature_subset(c("PG:1", "PG:missing"), mat, "missing synthetic subset"), "PG:missing")
  testthat::expect_error(joint_validate_feature_subset(c("PG:1", "PG:1"), mat, "duplicate requested subset"), "duplicated")
  dup_mat <- mat; rownames(dup_mat)[[4]] <- "PG:3"
  testthat::expect_error(joint_validate_feature_subset("PG:1", dup_mat, "duplicate matrix subset"), "invalid or duplicated")
})

testthat::test_that("all feature universes retain exact canonical row-name contracts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, dimnames = list(c("PG:1", "PG:2", "PG:3"), c("s1", "s2")))
  primary_ids <- c("PG:3", "PG:1"); complete_ids <- c("PG:1"); union_ids <- c("PG:2", "PG:3")
  primary <- joint_validate_feature_subset(primary_ids, mat, "primary")
  complete <- joint_validate_feature_subset(complete_ids, mat, "complete")
  union <- joint_validate_feature_subset(union_ids, mat, "union")
  testthat::expect_identical(rownames(primary), primary_ids)
  testthat::expect_identical(rownames(complete), complete_ids)
  testthat::expect_identical(rownames(union), union_ids)
})

testthat::test_that("source-row alignment and contaminant exclusion cannot leave stale feature IDs", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  feature_table <- data.frame(ProteinGroupID = c("PG:keep", "PG:contaminant", "PG:keep2"), source_row_id = 1:3, joint_qc_eligible = c(TRUE, FALSE, TRUE), stringsAsFactors = FALSE)
  mat <- matrix(c(4, 8, 16, 32, 64, 128), nrow = 3, dimnames = list(feature_table$ProteinGroupID, c("s1", "s2")))
  joint_validate_canonical_feature_matrix(feature_table, mat, "all source rows")
  retained_ids <- feature_table$ProteinGroupID[feature_table$joint_qc_eligible]
  retained <- joint_validate_feature_subset(retained_ids, joint_qc_log2_positive(mat)$log2, "contaminant-excluded subset")
  testthat::expect_identical(rownames(retained), retained_ids)
  testthat::expect_equal(nrow(retained), 2L)
  testthat::expect_error(joint_validate_canonical_feature_matrix(feature_table, mat[-1, , drop = FALSE], "dropped source row"), "row count")
})

testthat::test_that("complete case and strict GCT contracts reject missing or nonnumeric contents", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  mat <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("PG1", "PG2"), c("s1", "s2")))
  meta <- data.frame(Sample = c("s1", "s2"), dataset = c("neuron_neuropil", "microglia"), stringsAsFactors = FALSE)
  ft <- data.frame(ProteinGroupID = c("PG1", "PG2"), FeatureDisplayLabel = c("A", "B"), stringsAsFactors = FALSE)
  file <- tempfile(fileext = ".gct")
  joint_qc_write_gct_v13(mat, meta, ft, file)
  testthat::expect_silent(joint_qc_validate_gct_v13(file, 2, 2, 2, 1))
  testthat::expect_error(joint_qc_write_gct_v13(matrix(c(1, NA), nrow = 1, dimnames = list("PG1", c("s1", "s2"))), meta, ft[1, ], tempfile(fileext = ".gct")), "complete finite")
})

testthat::test_that("raw detection provenance never substitutes post-imputation availability", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "joint_compartment_qc_utils.R"))
  root <- tempfile(); dir.create(root)
  utils::write.csv(data.frame(ProteinGroupID = "PG1", dataset = "microglia", observed_detection_rate_raw = .5), file.path(root, "observed_detection_by_dataset.csv"), row.names = FALSE)
  provenance <- joint_qc_observed_detection_provenance(c("PG1", "PG2"), "microglia", root)
  testthat::expect_equal(provenance$observed_detection_rate_raw, c(.5, NA_real_))
  testthat::expect_false(any(provenance$observed_detection_rate_raw == 1, na.rm = TRUE))
  testthat::expect_equal(provenance$detectability_source[[1]], "raw_unified_pre_imputation")
})

testthat::test_that("fixed seeds make PCA, UMAP, and t-SNE reproducible", {
  x <- matrix(seq_len(72), nrow = 9) + rep(seq(0, 2, length.out = 8), each = 9)
  pca_a <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  pca_b <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  testthat::expect_equal(pca_a$x, pca_b$x)
  testthat::skip_if_not_installed("uwot")
  set.seed(42); u1 <- uwot::umap(pca_a$x[, 1:2], n_neighbors = 3, n_threads = 1, verbose = FALSE)
  set.seed(42); u2 <- uwot::umap(pca_a$x[, 1:2], n_neighbors = 3, n_threads = 1, verbose = FALSE)
  testthat::expect_equal(u1, u2)
  testthat::skip_if_not_installed("Rtsne")
  set.seed(42); t1 <- Rtsne::Rtsne(pca_a$x[, 1:2], dims = 2, perplexity = 2, pca = FALSE, check_duplicates = FALSE, verbose = FALSE, max_iter = 250)$Y
  set.seed(42); t2 <- Rtsne::Rtsne(pca_a$x[, 1:2], dims = 2, perplexity = 2, pca = FALSE, check_duplicates = FALSE, verbose = FALSE, max_iter = 250)$Y
  testthat::expect_equal(t1, t2)
})

testthat::test_that("joint scripts keep the global registry and fixed embedding seed contracts", {
  root <- testthat::test_path("..", "..")
  prep <- paste(readLines(file.path(root, "01_preprocessing", "01_prepare_joint_protigy_input.r"), warn = FALSE), collapse = "\n")
  qc <- paste(readLines(file.path(root, "03_qc_exploration", "00b_joint_compartment_qc.r"), warn = FALSE), collapse = "\n")
  registry <- paste(readLines(file.path(root, "pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_match(prep, "PROTEOMICS_JOINT_MIN_DETECTION_PER_BLOCK", fixed = TRUE)
  testthat::expect_match(prep, "joint_qc_write_gct_v13", fixed = TRUE)
  testthat::expect_match(qc, "set.seed(seed)", fixed = TRUE)
  testthat::expect_match(registry, "01_preprocessing/01_prepare_joint_protigy_input.r", fixed = TRUE)
  testthat::expect_match(registry, "03_qc_exploration/00b_joint_compartment_qc.r", fixed = TRUE)
  testthat::expect_false(grepl("all_compartments", registry, fixed = TRUE))
})
