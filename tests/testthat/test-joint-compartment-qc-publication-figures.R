source(testthat::test_path("..", "..", "R", "plotting_nature.R"))
source(testthat::test_path("..", "..", "R", "joint_compartment_qc_plotting.R"))

testthat::test_that("publication dataset palette and labels are exact", {
  testthat::expect_identical(
    joint_pub_dataset_palette(),
    c(
      "Microglia-enriched ROI" = "#A8D5CF",
      "Neuropil" = "#2F6F62",
      "Soma" = "#7F7F7F"
    )
  )
  observed <- joint_pub_dataset_factor(c("microglia", "neuron_neuropil", "neuron_soma"))
  testthat::expect_identical(as.character(observed), c("Microglia-enriched ROI", "Neuropil", "Soma"))
  testthat::expect_error(joint_pub_dataset_factor("purified_microglia"), "Unknown")
})

testthat::test_that("experimental groups are discrete CON RES SUS labels", {
  observed <- joint_pub_group_factor(c(1, 2, 3, "control", "resilient", "susceptible"))
  testthat::expect_true(is.factor(observed))
  testthat::expect_identical(levels(observed), c("CON", "RES", "SUS"))
  testthat::expect_identical(as.character(observed), c("CON", "RES", "SUS", "CON", "RES", "SUS"))
  testthat::expect_error(joint_pub_group_factor("unknown"), "Unsupported")
})

testthat::test_that("scree PCs are ordered numerically", {
  x <- data.frame(
    PC = c("PC1", "PC10", "PC2", "PC12", "PC3"),
    variance_explained = c(.5, .01, .2, .005, .1),
    stringsAsFactors = FALSE
  )
  out <- joint_pub_order_pcs(x)
  testthat::expect_identical(as.character(out$PC), c("PC1", "PC2", "PC3", "PC10", "PC12"))
  testthat::expect_identical(out$pc_number, c(1L, 2L, 3L, 10L, 12L))
  testthat::expect_error(joint_pub_order_pcs(rbind(x, x[1, ])), "unique")
})

testthat::test_that("embedding limits use equal numeric spans without moving coordinates", {
  x <- c(-10, 30)
  y <- c(2, 12)
  limits <- joint_pub_equal_limits(x, y, padding = 0)
  testthat::expect_equal(diff(limits$x), diff(limits$y))
  testthat::expect_true(all(x >= limits$x[1] & x <= limits$x[2]))
  testthat::expect_true(all(y >= limits$y[1] & y <= limits$y[2]))
  testthat::expect_identical(x, c(-10, 30))
  testthat::expect_identical(y, c(2, 12))
})

testthat::test_that("missingness features use deterministic between-dataset variation", {
  detected <- rbind(
    PG_constant = c(1, 1, 1, 1, 1, 1),
    PG_variable = c(1, 1, 0, 0, 0, 0),
    PG_partial = c(1, 1, 1, 0, 1, 0),
    PG_tie_b = c(1, 1, 0, 1, 0, 0),
    PG_tie_a = c(1, 1, 0, 1, 0, 0)
  )
  colnames(detected) <- paste0("s", seq_len(ncol(detected)))
  dataset <- rep(c("microglia", "neuron_neuropil", "neuron_soma"), each = 2)
  out <- joint_pub_select_detection_features(detected, dataset, n = 4)
  testthat::expect_identical(out$ProteinGroupID[[1]], "PG_variable")
  testthat::expect_true(out$between_dataset_detection_range[[1]] > out$between_dataset_detection_range[[4]])
  tied <- out$ProteinGroupID[out$ProteinGroupID %in% c("PG_tie_a", "PG_tie_b")]
  testthat::expect_identical(tied, sort(tied))
  testthat::expect_false("PG_constant" %in% out$ProteinGroupID)
})

testthat::test_that("sample order is dataset plate and anatomy deterministic", {
  meta <- data.frame(
    Sample = c("n2", "m2", "s1", "m1", "n1"),
    dataset = c("neuron_neuropil", "microglia", "neuron_soma", "microglia", "neuron_neuropil"),
    plate = c("C", "B", "B", "B", "B"),
    region = c("CA1", "DG", "CA1", "CA1", "CA2"),
    layer = c("slm", "microglia", "sp", "microglia", "sr"),
    stringsAsFactors = FALSE
  )
  testthat::expect_identical(meta$Sample[joint_pub_order_samples(meta)], c("m1", "m2", "n1", "n2", "s1"))
})

testthat::test_that("coordinate validation rejects any changed or missing sample", {
  reference <- data.frame(Sample = c("s1", "s2"), X = c(1, 2), Y = c(3, 4))
  candidate <- reference[c(2, 1), ]
  testthat::expect_silent(joint_pub_assert_coordinates(reference, candidate, c("X", "Y"), "synthetic"))
  candidate$X[[1]] <- 99
  testthat::expect_error(joint_pub_assert_coordinates(reference, candidate, c("X", "Y"), "synthetic"), "changed")
})

testthat::test_that("publication entrypoint is rendering-only and output-isolated", {
  root <- testthat::test_path("..", "..")
  script <- paste(readLines(file.path(root, "03_qc_exploration", "00c_joint_compartment_qc_publication_figures.r"), warn = FALSE), collapse = "\n")
  registry <- paste(readLines(file.path(root, "pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "publication_style", fixed = TRUE)
  testthat::expect_match(script, "joint_primary_pca_scores.csv", fixed = TRUE)
  testthat::expect_match(script, "joint_umap_exploratory_scores.csv", fixed = TRUE)
  testthat::expect_match(script, "joint_tsne_exploratory_scores.csv", fixed = TRUE)
  testthat::expect_false(grepl("prcomp\\s*\\(", script))
  testthat::expect_false(grepl("uwot::umap\\s*\\(", script))
  testthat::expect_false(grepl("Rtsne::Rtsne\\s*\\(", script))
  testthat::expect_false(grepl("write\\.csv\\s*\\(", script))
  testthat::expect_match(registry, "03_qc_exploration/00c_joint_compartment_qc_publication_figures.r", fixed = TRUE)
})
