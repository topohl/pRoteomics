source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "animal_spatial_network_utils.R"))

NET <- function(...) path_results("tables", "11_spatial_systems", "networks", ...)
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

# Source scans below must look at CODE, not at the comments that explain why a
# defect is absent - otherwise documenting the old defect would fail the test.
code_of <- function(...) {
  l <- readLines(repo_path(...), warn = FALSE)
  l <- sub("#.*$", "", l)
  paste(l[nzchar(trimws(l))], collapse = "\n")
}

fx_matrix <- function(seed = 1, n_prot = 40, nodes = c("ca1_slm", "ca1_so", "ca2_slm", "ca2_so")) {
  set.seed(seed)
  m <- matrix(stats::rnorm(n_prot * length(nodes)), nrow = n_prot,
              dimnames = list(paste0("P", seq_len(n_prot)), nodes))
  m
}

# =====================================================================
# representations and node handling
# =====================================================================

testthat::test_that("representations are declared prospectively and transform correctly", {
  testthat::expect_identical(asn_representations(),
    c("A_abundance", "B_protein_centered", "C_protein_z", "D_centered_varfilter"))
  m <- fx_matrix()
  a <- asn_transform(m, "A_abundance")
  testthat::expect_equal(a$mat, m)
  b <- asn_transform(m, "B_protein_centered")
  # every protein's row mean is removed, so relative distribution is what remains
  testthat::expect_true(all(abs(rowMeans(b$mat)) < 1e-12))
  cz <- asn_transform(m, "C_protein_z")
  testthat::expect_true(all(abs(apply(cz$mat, 1, stats::sd) - 1) < 1e-8))
  for (r in asn_representations()) {
    t <- asn_transform(m, r)
    testthat::expect_true(t$n_proteins_retained <= t$n_proteins_input)
    testthat::expect_true(nzchar(asn_representation_description(r)))
  }
})

testthat::test_that("zero-spatial-variance proteins are handled deterministically", {
  m <- fx_matrix()
  m[1, ] <- 5            # exactly flat
  m[2, ] <- c(5, 5, 5, 5.0001)  # far below the threshold
  thr <- asn_min_spatial_sd()
  testthat::expect_gt(thr, 0)
  cz <- asn_transform(m, "C_protein_z")
  d <- asn_transform(m, "D_centered_varfilter")
  # the flat protein cannot be z-scored and must be excluded, not turned into NaN
  testthat::expect_false("P1" %in% rownames(cz$mat))
  testthat::expect_false("P1" %in% rownames(d$mat))
  testthat::expect_false(any(is.nan(cz$mat)))
  # and the decision is deterministic
  testthat::expect_identical(asn_transform(m, "C_protein_z")$n_proteins_retained,
                             cz$n_proteins_retained)
  # representation A keeps everything, so retention differences are visible
  testthat::expect_identical(asn_transform(m, "A_abundance")$n_proteins_retained,
                             nrow(m))
})

testthat::test_that("node ordering is enforced and edges are unique upper-triangle", {
  m <- fx_matrix()
  nodes <- colnames(m)
  s <- asn_similarity_matrix(m, nodes)
  testthat::expect_identical(rownames(s), nodes)
  testthat::expect_identical(colnames(s), nodes)
  v <- asn_edge_vector(s)
  testthat::expect_identical(length(v), as.integer(choose(length(nodes), 2)))
  testthat::expect_identical(anyDuplicated(names(v)), 0L)
  # reordering the input columns must not change the network
  s2 <- asn_similarity_matrix(m[, rev(nodes), drop = FALSE], nodes)
  testthat::expect_equal(s, s2)
  # a missing node is a hard error, never a silently smaller network
  testthat::expect_error(asn_similarity_matrix(m[, 1:2, drop = FALSE], nodes),
                         "missing node")
})

testthat::test_that("the full weighted matrix is independent of any threshold", {
  m <- fx_matrix()
  s <- asn_similarity_matrix(m, colnames(m))
  v <- asn_edge_vector(s)
  # every pair is present regardless of magnitude
  testthat::expect_identical(length(v), 6L)
  testthat::expect_true(any(abs(v) < 0.5))
  # no code path removes edges by size
  src <- code_of("R", "animal_spatial_network_utils.R")
  testthat::expect_false(grepl("threshold", src, ignore.case = TRUE))
  testthat::expect_false(grepl("abs\\([^)]*\\)\\s*[<>]=?\\s*0[.]", src))
})

# =====================================================================
# Fisher z and centroid distance
# =====================================================================

testthat::test_that("Fisher z preserves matrix shape", {
  m <- matrix(c(0.1, 0.5, -0.3, 0.9), nrow = 2,
              dimnames = list(c("a", "b"), c("e1", "e2")))
  z <- asn_fisher_z(m)
  # as.numeric() on a matrix silently drops dim; that would break every
  # downstream animal x edge operation
  testthat::expect_identical(dim(z), dim(m))
  testthat::expect_identical(dimnames(z), dimnames(m))
  testthat::expect_equal(z[1, 1], atanh(0.1))
  testthat::expect_true(is.finite(asn_fisher_z(1)))    # bounded, not Inf
  testthat::expect_true(is.finite(asn_fisher_z(-1)))
})

testthat::test_that("CON centroid is leave-one-out for CON and full-CON for others", {
  set.seed(5)
  em <- matrix(stats::rnorm(6 * 4), nrow = 6,
               dimnames = list(paste0("A", 1:6), paste0("e", 1:4)))
  grp <- c("CON", "CON", "CON", "SUS", "SUS", "RES")
  d <- asn_distance_from_con(em, grp)
  testthat::expect_length(d, 6L)
  testthat::expect_true(all(is.finite(d)))

  # a CON animal must NOT be compared against a centroid containing itself:
  # duplicating a CON animal would collapse a self-inclusive distance to ~0
  em2 <- em; em2[2, ] <- em2[1, ]
  d2 <- asn_distance_from_con(em2, grp)
  testthat::expect_gt(d2[1], 1e-8)

  # RES/SUS distances use ALL CON animals, so they are unaffected by which CON
  # animal is held out
  z <- asn_fisher_z(em)
  centroid_all <- colMeans(z[1:3, , drop = FALSE])
  testthat::expect_equal(d[4], sqrt(sum((z[4, ] - centroid_all)^2)))
  # too few CON animals yields NA rather than a fabricated number
  testthat::expect_true(all(is.na(asn_distance_from_con(em, rep("SUS", 6)))))
})

# =====================================================================
# exact permutation and bootstrap
# =====================================================================

testthat::test_that("SUS vs RES enumerates 20 whole-animal assignments", {
  v <- c(5, 6, 7, 1, 2, 3); g <- rep(c("SUS", "RES"), each = 3)
  r <- asn_exact_two_group_p(v, g, "SUS", "RES")
  testthat::expect_identical(r$n_assignments, 20L)
  # THE documented ceiling: 3 vs 3 cannot produce a two-sided p below 0.10
  testthat::expect_equal(r$min_attainable_two_sided_p, 0.10)
  testthat::expect_gte(r$p_two_sided, 2 / 21)
  # a maximally separated split still cannot reach 0.05
  testthat::expect_gt(r$p_two_sided, 0.05)
  testthat::expect_equal(r$observed, mean(c(5, 6, 7)) - mean(c(1, 2, 3)))
})

testthat::test_that("the three-group omnibus enumerates 1680 assignments", {
  v <- as.numeric(1:9); g <- rep(c("CON", "RES", "SUS"), each = 3)
  r <- asn_exact_three_group_p(v, g)
  testthat::expect_identical(r$n_assignments, 1680L)
  testthat::expect_gte(r$p, 1 / 1681)
  testthat::expect_lte(r$p, 1)
})

testthat::test_that("permutation operates on whole animals, not on edges", {
  # an animal's edges must move together: permuting a multi-edge animal vector
  # must give the same assignment count as a single metric
  v <- c(1, 2, 3, 7, 8, 9); g <- rep(c("RES", "SUS"), each = 3)
  r1 <- asn_exact_two_group_p(v, g, "SUS", "RES")
  r2 <- asn_exact_two_group_p(v * 3, g, "SUS", "RES")
  testthat::expect_identical(r1$n_assignments, r2$n_assignments)
  # the function takes one value per animal; it has no edge dimension at all
  testthat::expect_identical(names(formals(asn_exact_two_group_p)),
                             c("values", "group", "a", "b"))
})

testthat::test_that("the bootstrap resamples animals and records its state", {
  set.seed(2)
  v <- c(1, 2, 3, 8, 9, 10); g <- rep(c("RES", "SUS"), each = 3)
  b <- asn_animal_bootstrap_difference(v, g, "SUS", "RES", n_boot = 500L, seed = 7L)
  testthat::expect_identical(b$bootstrap_unit, "AnimalID")
  testthat::expect_identical(b$seed, 7L)
  testthat::expect_gt(b$n_valid_iterations, 0L)
  testthat::expect_true(b$ci_lower <= b$difference && b$difference <= b$ci_upper)
  # deterministic given the seed
  b2 <- asn_animal_bootstrap_difference(v, g, "SUS", "RES", n_boot = 500L, seed = 7L)
  testthat::expect_equal(b$ci_lower, b2$ci_lower)
  # a single animal cannot support an interval
  one <- asn_animal_bootstrap_ci(c(1), c("SUS"), "SUS")
  testthat::expect_identical(one$n_valid_iterations, 0L)
})

testthat::test_that("no bootstrap sign frequency is ever produced as a p-value", {
  src <- code_of("R", "animal_spatial_network_utils.R")
  testthat::expect_false(grepl("Prob_DeltaR", src, fixed = TRUE))
  testthat::expect_false(grepl("permutation_p", src, fixed = TRUE))
  # and no producer emits that field name, nor a BH FDR over sign probabilities
  for (f in c("13_animal_spatial_networks.R", "14_network_group_organization.R")) {
    s <- code_of("11_spatial_systems", f)
    testthat::expect_false(grepl("permutation_p", s, fixed = TRUE), info = f)
    testthat::expect_false(grepl("Prob_DeltaR", s, fixed = TRUE), info = f)
    testthat::expect_false(grepl("p\\.adjust", s), info = f)
  }
})

# =====================================================================
# generated outputs
# =====================================================================

testthat::test_that("there is one bilateral network per AnimalID with no duplicates", {
  p <- NET("animal_network_edges.csv")
  testthat::skip_if_not(file.exists(p), "animal networks not built")
  e <- rd(p)
  b <- e[e$hemisphere_mode == "bilateral", , drop = FALSE]
  testthat::expect_identical(anyDuplicated(paste(b$dataset, b$AnimalID, b$edge_id)), 0L)
  # every animal carries the identical edge set within a dataset
  for (ds in unique(b$dataset)) {
    z <- b[b$dataset == ds, , drop = FALSE]
    per <- split(z$edge_id, z$AnimalID)
    testthat::expect_true(all(vapply(per, function(v)
      identical(sort(v), sort(per[[1]])), logical(1))), info = ds)
  }
  testthat::expect_setequal(unique(e$hemisphere_mode), c("L", "R", "bilateral"))
})

testthat::test_that("the representation was selected without any phenotype input", {
  p <- NET("network_representation_selection.csv")
  testthat::skip_if_not(file.exists(p))
  s <- rd(p)
  testthat::expect_true(all(grepl("CON", s$selection_criteria)))
  # no SUS/RES quantity may appear as a selection column or criterion
  testthat::expect_false(any(grepl("SUS|RES", s$selection_criteria)))
  testthat::expect_false(any(grepl("SUS|RES", names(s))))
  testthat::expect_true("n_con_animals" %in% names(s))
  testthat::expect_true(any(s$is_primary %in% TRUE))
  # raw abundance must not win: it fails the anatomical validity floor
  np <- s[s$dataset == "neuron_neuropil", , drop = FALSE]
  if (nrow(np)) {
    testthat::expect_false(np$representation[np$is_primary %in% TRUE][1] == "A_abundance")
  }
})

testthat::test_that("edge-level output reports exact coarse resolution, never an FDR", {
  p <- NET("network_edge_group_differences.csv")
  testthat::skip_if_not(file.exists(p))
  e <- rd(p)
  testthat::expect_false("permutation_p" %in% names(e))
  testthat::expect_false(any(grepl("^fdr$|_fdr$", names(e), ignore.case = TRUE)))
  testthat::expect_true(all(e$exact_n_assignments == 20L, na.rm = TRUE))
  testthat::expect_true(all(abs(e$min_attainable_two_sided_p - 0.10) < 1e-9, na.rm = TRUE))
  testthat::expect_true(all(grepl("EXACT but COARSE", e$inference_note)))
  # nothing may be below the mathematically attainable floor
  testthat::expect_true(all(e$exact_p_two_sided >= 2 / 21 - 1e-9, na.rm = TRUE))
})

testthat::test_that("legacy outputs are labelled, not rewritten", {
  p <- NET("legacy_spatial_network_comparison.csv")
  testthat::skip_if_not(file.exists(p), "legacy comparison not generated")
  l <- rd(p)
  testthat::expect_true(all(l$legacy_inference_status == "legacy_noncanonical_inference"))
  testthat::expect_true(all(grepl("SIGN FREQUENCY", l$legacy_defect)))
  # the comparison is a NEW file; the historical tables keep their own paths
  testthat::expect_true(grepl("11_spatial_systems", p))
})

testthat::test_that("the network validation contract has no critical failure", {
  p <- NET("spatial_network_validation_status.csv")
  testthat::skip_if_not(file.exists(p), "network validation not generated")
  v <- rd(p)
  testthat::expect_gt(nrow(v), 10L)
  testthat::expect_equal(sum(v$critical %in% TRUE & v$status == "FAIL"), 0L)
})

testthat::test_that("network figures each have exact source data", {
  d <- path_results("figures", "11_spatial_systems", "networks")
  testthat::skip_if_not(dir.exists(d), "network figures not generated")
  pngs <- list.files(d, pattern = "[.]png$")
  testthat::skip_if(length(pngs) == 0L)
  for (f in pngs) {
    testthat::expect_true(
      file.exists(file.path(d, sub("[.]png$", "_source_data.csv", f))), info = f)
  }
})

testthat::test_that("the CA2-SLM precision audit records the revised conclusion", {
  p <- path_results("tables", "11_spatial_systems", "atlas",
                    "neuropil_spatial_precision_context.csv")
  testthat::skip_if_not(file.exists(p), "precision audit not generated")
  x <- rd(p)
  testthat::expect_equal(nrow(x), 10L)
  # SE is identical across units because the canonical DA shares one residual
  # variance per protein; precision cannot explain the concentration
  testthat::expect_identical(length(unique(round(x$median_model_SE, 8))), 1L)
  testthat::expect_true(all(grepl("no DA refit", x$se_derivation)))
  # CA2_slm is the MOST missing unit, not the least
  ca2 <- x$missing_fraction_preimputation[x$spatial_unit == "ca2_slm"]
  testthat::expect_equal(ca2, max(x$missing_fraction_preimputation))
  a <- path_results("tables", "11_spatial_systems", "atlas",
                    "ca2_slm_imputation_artifact_test.csv")
  if (file.exists(a)) {
    y <- rd(a)
    testthat::expect_true(all(c("differential_missingness_SUS_minus_RES", "is_DAP",
                                "fully_observed") %in% names(y)))
    testthat::expect_equal(sum(y$is_DAP %in% TRUE), 28L)
  }
})
