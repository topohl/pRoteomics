source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "ca2_slm_robustness_utils.R"))

# Missingness has to be read before imputation; the imputed matrices carry no NA.

testthat::test_that("sample missingness comes from the raw matrix over the analysed rows", {
  raw <- data.frame(Protein.Group = paste0("PG", 1:5), Protein.Names = paste0("P", 1:5),
                    Genes = paste0("G", 1:5), First.Protein.Description = "d",
                    A = c(10, NA, NA, NA, NA), B = c(11, 12, NA, NA, NA),
                    C = c(12, 13, 14, NA, NA), D = c(13, 14, 15, 16, NA), check.names = FALSE)
  f <- tempfile(fileext = ".tsv"); utils::write.table(raw, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  m <- qc_preimputation_missingness(f, c("A", "B", "C", "D"))
  # rows 4 (3 of 4 missing) and 5 (all missing) exceed the 70% filter; three analysed rows remain
  testthat::expect_identical(unique(m$n_proteins_total), 3L)
  testthat::expect_identical(m$n_proteins_detected, c(1L, 2L, 3L, 3L))
  testthat::expect_equal(m$missing_fraction, c(2, 1, 0, 0) / 3)
  testthat::expect_error(qc_preimputation_missingness(f, c("A", "Z")), "lacks 1 of the 2 QC sample columns")
})

testthat::test_that("the CA2 displacement is measured on the run's own centred values, matched by Protein.Group", {
  smp <- c("S1L", "S1R", "S2L", "S2R", "S3L", "S3R")
  meta <- data.frame(sample_id = smp, AnimalID = rep(c("1", "2", "3"), each = 2L), hemisphere = c("L", "R"),
                     StressGroup = rep(c("SUS", "RES", "RES"), each = 2L), qc_flag = "PASS", stringsAsFactors = FALSE)
  lg <- matrix(seq(15, by = 0.37, length.out = 48) %% 7 + 15, 8, 6, dimnames = list(NULL, smp))
  lg[7, "S1L"] <- NA                       # observed in 5 of 6 samples: analysed but not always observed
  lg[8, ] <- NA; lg[8, 1] <- 18            # 5 of 6 missing: removed by the 70% filter
  raw <- data.frame(Protein.Group = paste0("PG", 1:8), Protein.Names = c("P1", "P2", "P3", "P4", "P5", "P1", "P7", "P8"),
                    2^lg, check.names = FALSE)
  f <- tempfile(fileext = ".tsv"); utils::write.table(raw, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  offsets <- c(S1L = 1.0, S1R = 0.2, S2L = 0.1, S2R = 0.3, S3L = 0.0, S3R = 0.4)
  run <- sweep(lg[1:7, ], 2, offsets, "-"); run[is.na(run)] <- 9   # the run's centred values, one drawn value
  imp <- data.frame(Protein.Group = paste0("PG", 1:7), Protein.Names = raw$Protein.Names[1:7], run, check.names = FALSE)
  imp <- imp[c(7, 3, 1, 5, 2, 6, 4), ]     # another row order: rows must be matched by key, not position
  x <- tempfile(fileext = ".xlsx"); openxlsx::write.xlsx(imp, x)

  per_sample <- colMeans(lg[1:6, ]) - offsets
  animal <- tapply(per_sample, meta$AnimalID, mean)
  expected <- animal[["1"]] - mean(animal[c("2", "3")])
  got <- csr_normalization_bias_context(f, smp, meta, imputed_path = x)
  testthat::expect_equal(got$sus_minus_res_displacement, expected, tolerance = 1e-10)
  testthat::expect_identical(got$n_always_observed, 6L)
  testthat::expect_match(got$centred_values_source, "^this run's centred sample values")
  testthat::expect_true(all(got$per_sample$centred_values_source == got$centred_values_source))

  # without the run's values the raw observed medians stand in, which is a different quantity
  raw_only <- csr_normalization_bias_context(f, smp, meta)
  med <- apply(lg[1:7, ], 2, stats::median, na.rm = TRUE)
  ps <- colMeans(sweep(lg[1:6, ], 2, med, "-")); an <- tapply(ps, meta$AnimalID, mean)
  testthat::expect_equal(raw_only$sus_minus_res_displacement, an[["1"]] - mean(an[c("2", "3")]), tolerance = 1e-10)
  testthat::expect_match(raw_only$centred_values_source, "^raw log2 values")

  # a missing always-observed protein group in the run's matrix stops the audit
  openxlsx::write.xlsx(imp[imp$Protein.Group != "PG2", ], x, overwrite = TRUE)
  testthat::expect_error(csr_normalization_bias_context(f, smp, meta, imputed_path = x), "absent from the imputed matrix")
})
