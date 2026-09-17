source(repo_path("R", "wgcna_candidate_protein_utils.R"))
source(testthat::test_path("..", "..", "R", "paths.R"))

# --------------------------------------------------------------- fixtures

wcp_fixture_members <- function() {
  data.frame(
    dataset = "neuron_neuropil",
    ModuleID = c(rep("WGCNA_m01", 12L), rep("WGCNA_m02", 3L)),
    ProteinGroupID = sprintf("PG:neuron_neuropil:%02d", 1:15),
    abs_kME = c(0.99, 0.97, 0.95, 0.93, 0.91, 0.89, 0.87, 0.85, 0.83, 0.81,
                0.55, NA_real_, 0.70, 0.70, 0.10),
    kME = c(0.99, -0.97, 0.95, -0.93, 0.91, 0.89, 0.87, 0.85, 0.83, 0.81,
            0.55, NA_real_, 0.70, -0.70, 0.10),
    mapping_status = "mapped",
    gene_level_claim_allowed = TRUE,
    protein_level_claim_allowed = TRUE,
    protein_group_ambiguity_class = "single_accession_single_gene",
    stringsAsFactors = FALSE
  )
}

wcp_fixture_da <- function() {
  expand_grid_df <- expand.grid(
    ProteinGroupID = sprintf("PG:neuron_neuropil:%02d", 1:15),
    spatial_unit = c("CA1_so", "CA1_sr", "CA2_slm", "DG_mo"),
    contrast = "SUS - RES",
    stringsAsFactors = FALSE
  )
  expand_grid_df$dataset <- "neuron_neuropil"
  set.seed(42)
  expand_grid_df$log2FC <- seq(-1.5, 1.5, length.out = nrow(expand_grid_df))
  expand_grid_df$BH_FDR <- 0.5
  expand_grid_df
}

# ----------------------------------------------- module ranking / hub flags

testthat::test_that("top-5 and top-10 flags follow within-module abs_kME order", {
  ranked <- wcp_rank_module_members(wcp_fixture_members())

  m1 <- ranked[ranked$ModuleID == "WGCNA_m01", ]
  m1 <- m1[order(m1$abs_kME_rank_in_module), ]
  finite_m1 <- m1[is.finite(m1$abs_kME), ]

  # ranks are a strict 1..n ordering of decreasing abs_kME
  testthat::expect_identical(
    finite_m1$abs_kME_rank_in_module, seq_len(nrow(finite_m1))
  )
  testthat::expect_true(all(diff(finite_m1$abs_kME) < 0))

  # flags are exactly the first 5 / first 10 of that ordering
  testthat::expect_identical(
    which(finite_m1$is_top5_module_representative), 1:5
  )
  testthat::expect_identical(which(finite_m1$is_top10_module_hub), 1:10)

  # a non-finite abs_kME is never a hub and never receives a rank
  na_row <- ranked[is.na(ranked$abs_kME), ]
  testthat::expect_true(is.na(na_row$abs_kME_rank_in_module))
  testthat::expect_false(na_row$is_top5_module_representative)
  testthat::expect_false(na_row$is_top10_module_hub)

  # module smaller than the window: every member is a representative
  m2 <- ranked[ranked$ModuleID == "WGCNA_m02", ]
  testthat::expect_true(all(m2$is_top5_module_representative))
  testthat::expect_identical(unique(m2$n_module_members), 3L)
})

testthat::test_that("ties break deterministically on ProteinGroupID", {
  tied <- data.frame(
    dataset = "d", ModuleID = "M",
    ProteinGroupID = c("PG:d:c", "PG:d:a", "PG:d:b"),
    abs_kME = c(0.5, 0.5, 0.5), stringsAsFactors = FALSE
  )
  ranked <- wcp_rank_module_members(tied)
  testthat::expect_identical(
    ranked$abs_kME_rank_in_module[order(ranked$ProteinGroupID)], 1:3
  )
  # shuffling input rows must not change any assigned rank
  shuffled <- wcp_rank_module_members(tied[c(3L, 1L, 2L), ])
  key <- function(x) x$abs_kME_rank_in_module[order(x$ProteinGroupID)]
  testthat::expect_identical(key(ranked), key(shuffled))
})

testthat::test_that("ranking rejects duplicate module identities", {
  dup <- wcp_fixture_members()[c(1L, 1L), ]
  testthat::expect_error(
    wcp_rank_module_members(dup),
    "duplicate rows for dataset + ModuleID + ProteinGroupID", fixed = TRUE
  )
})

testthat::test_that("high_kME reproduces frozen is_core_kME_0.6 semantics", {
  # The frozen upstream flag is is.finite(abs_kME) & abs_kME >= 0.6.
  abs_kme <- c(0.6, 0.5999999, 0.99, NA_real_, Inf, -0.7)
  frozen <- is.finite(abs_kme) & abs_kme >= 0.6
  derived <- is.finite(abs_kme) & abs_kme >= wcp_high_kme_threshold()
  testthat::expect_identical(derived, frozen)
  testthat::expect_identical(wcp_high_kme_threshold(), 0.60)
  testthat::expect_false(derived[[4]])   # NA never becomes TRUE
})

# ------------------------------------------------------------ mapping / QC

testthat::test_that("mapping ambiguity is never silently upgraded to a clean claim", {
  members <- data.frame(
    mapping_status = c("mapped", "mapped", "mapped", "unmapped", "mapped"),
    gene_level_claim_allowed = c(TRUE, FALSE, TRUE, TRUE, NA),
    protein_level_claim_allowed = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    protein_group_ambiguity_class = c(
      "single_accession_single_gene", "single_accession_single_gene",
      "multi_gene_indistinguishable", "single_accession_single_gene",
      "single_accession_single_gene"
    ), stringsAsFactors = FALSE
  )
  clean <- wcp_clean_mapping(members)
  testthat::expect_identical(clean, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  # the gene-ambiguous group is specifically excluded
  testthat::expect_false(clean[[3]])
  # an NA claim flag fails closed rather than defaulting to TRUE
  testthat::expect_false(clean[[5]])
})

# ---------------------------------------------------- direction concordance

testthat::test_that("direction concordance is NA-safe and zero-safe", {
  observed <- wcp_direction_matches_module(
    protein_effect  = c(1, -1, 1, -1, 0, 2, NA, 3),
    module_estimate = c(1, -1, -1, 1, 1, 0, 1, NA)
  )
  testthat::expect_identical(
    observed, c(TRUE, TRUE, FALSE, FALSE, NA, NA, NA, NA)
  )
  # a zero or missing module estimate must never read as disagreement
  testthat::expect_true(is.na(observed[[5]]))
  testthat::expect_true(is.na(observed[[6]]))
  testthat::expect_false(isTRUE(observed[[7]]))
  testthat::expect_true(is.na(observed[[8]]))
})

# ------------------------------------------------- relative effect families

testthat::test_that("relative effect threshold is computed only within its family", {
  da <- data.frame(
    dataset = c(rep("A", 10), rep("B", 10)),
    spatial_unit = "u1",
    contrast = "SUS - RES",
    # family A spans 1..10, family B spans 101..110
    log2FC = c(1:10, 101:110),
    stringsAsFactors = FALSE
  )
  th <- wcp_large_effect_thresholds(da)
  testthat::expect_identical(nrow(th), 2L)
  a <- th$large_effect_abs_log2FC_threshold[th$dataset == "A"]
  b <- th$large_effect_abs_log2FC_threshold[th$dataset == "B"]
  testthat::expect_equal(a, unname(quantile(1:10, 0.9, type = 7)))
  testthat::expect_equal(b, unname(quantile(101:110, 0.9, type = 7)))
  # family B's large values must not raise family A's threshold
  testthat::expect_true(a < 11)

  flagged <- wcp_flag_large_effect(da)
  # exactly the top of each family is flagged, independently
  testthat::expect_equal(
    as.numeric(flagged$log2FC[flagged$large_effect_within_context]), c(10, 110)
  )
})

testthat::test_that("typical-effect threshold uses dataset x contrast, not context", {
  consistency <- data.frame(
    dataset = c(rep("A", 10), rep("B", 10)),
    contrast = "SUS - RES",
    median_abs_log2FC = c(1:10, 101:110),
    stringsAsFactors = FALSE
  )
  out <- wcp_flag_typical_large_effect(consistency)
  testthat::expect_equal(
    as.numeric(out$median_abs_log2FC[out$large_effect_typical]), c(10, 110)
  )
  # by construction this selects ~the top decile, not a multiple of it
  testthat::expect_lte(mean(out$large_effect_typical), 0.2)
})

# ------------------------------------------------------ spatial consistency

testthat::test_that("spatial summaries reject duplicated context rows", {
  da <- wcp_fixture_da()
  dup <- rbind(da, da[1, ])
  testthat::expect_error(
    wcp_spatial_consistency(dup),
    "duplicate rows for dataset + ProteinGroupID + contrast + spatial_unit",
    fixed = TRUE
  )
})

testthat::test_that("spatial consistency counts each context exactly once", {
  da <- data.frame(
    dataset = "d", ProteinGroupID = "PG:d:1", contrast = "SUS - RES",
    spatial_unit = c("u1", "u2", "u3", "u4"),
    log2FC = c(0.5, 0.4, 0.3, -0.2),
    BH_FDR = c(0.01, 0.2, 0.3, 0.4), stringsAsFactors = FALSE
  )
  out <- wcp_spatial_consistency(da)
  testthat::expect_identical(nrow(out), 1L)
  testthat::expect_identical(out$n_spatial_contexts_tested, 4L)
  testthat::expect_identical(out$n_spatial_contexts_positive, 3L)
  testthat::expect_identical(out$n_spatial_contexts_negative, 1L)
  testthat::expect_identical(out$majority_direction, "positive")
  testthat::expect_identical(out$n_matching_majority_direction, 3L)
  testthat::expect_equal(out$fraction_matching_majority_direction, 0.75)
  testthat::expect_identical(out$n_spatial_contexts_fdr05, 1L)
  testthat::expect_equal(out$max_abs_log2FC, 0.5)
  testthat::expect_equal(out$median_log2FC, 0.35)
  # not unanimous, so not consistent
  testthat::expect_false(out$spatially_consistent)
})

testthat::test_that("consistency requires unanimity and calibrates to context count", {
  mk <- function(n, sign_flip = FALSE) {
    eff <- rep(0.3, n)
    if (sign_flip) eff[[1]] <- -0.3
    data.frame(
      dataset = "d", ProteinGroupID = "PG:d:1", contrast = "SUS - RES",
      spatial_unit = paste0("u", seq_len(n)), log2FC = eff, BH_FDR = 0.5,
      stringsAsFactors = FALSE
    )
  }
  four <- wcp_spatial_consistency(mk(4L))
  ten <- wcp_spatial_consistency(mk(10L))
  two <- wcp_spatial_consistency(mk(2L))
  mixed <- wcp_spatial_consistency(mk(10L, sign_flip = TRUE))

  testthat::expect_true(four$spatially_consistent)
  testthat::expect_true(ten$spatially_consistent)
  # below the minimum context count, no consistency claim is made
  testthat::expect_false(two$spatially_consistent)
  # non-unanimous is never consistent
  testthat::expect_false(mixed$spatially_consistent)

  # No independence-based sign probability is computed or exported. These
  # spatial contexts are correlated repeated measures from the same animals,
  # so a 2^(1-n) style null would understate how often unanimity arises.
  testthat::expect_false("spatial_direction_null_p" %in% names(four))
  testthat::expect_false("spatially_consistent_selective" %in% names(four))
})

testthat::test_that("no independence-derived probability survives anywhere", {
  files <- c(
    repo_path("R", "wgcna_candidate_protein_utils.R"),
    repo_path("analysis/08_integration",
                        "build_candidate_protein_shortlist.R")
  )
  code <- unlist(lapply(files, readLines, warn = FALSE))
  live <- code[!grepl("^\\s*#", code)]
  # the quantity, its threshold helper, and the derived flag are all gone
  testthat::expect_false(any(grepl("spatial_direction_null_p", live)))
  testthat::expect_false(any(grepl("spatially_consistent_selective", live)))
  testthat::expect_false(any(grepl("wcp_selective_null_alpha", live)))
  testthat::expect_false(any(grepl("2\\^\\(1", live)))
})

# ------------------------------------------------------------ candidate tiers

wcp_fixture_candidates <- function() {
  data.frame(
    dataset = "d",
    ModuleID = "M",
    ProteinGroupID = sprintf("PG:d:%02d", 1:6),
    abs_kME = c(0.90, 0.90, 0.90, 0.40, 0.90, 0.55),
    clean_mapping = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    is_top5_module_representative = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    is_top10_module_hub = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
    sus_res_fdr05_any_context = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    sus_res_large_effect_typical = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    sus_res_large_effect_any_context = FALSE,
    sus_res_spatially_consistent = FALSE,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("Tier A1 requires clean mapping, high kME, top-10 and FDR support", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())

  # row 1 satisfies every criterion
  testthat::expect_true(flagged$is_tier_A1[[1]])
  # row 2 fails only top-10 -> becomes A2, never A1
  testthat::expect_false(flagged$is_tier_A1[[2]])
  # row 3 fails only clean mapping -> ambiguity cannot buy a phenotype-linked label
  testthat::expect_false(flagged$is_tier_A1[[3]])
  # row 4 fails only high kME
  testthat::expect_false(flagged$is_tier_A1[[4]])
  # row 5 fails only FDR support
  testthat::expect_false(flagged$is_tier_A1[[5]])

  # dropping any single criterion must remove the A1 call for row 1
  for (col in c("clean_mapping", "is_top10_module_hub",
                "sus_res_fdr05_any_context")) {
    broken <- wcp_fixture_candidates()
    broken[[col]][[1]] <- FALSE
    testthat::expect_false(
      wcp_assign_candidate_flags(broken)$is_tier_A1[[1]],
      info = col
    )
  }
  broken_kme <- wcp_fixture_candidates()
  broken_kme$abs_kME[[1]] <- 0.59
  testthat::expect_false(wcp_assign_candidate_flags(broken_kme)$is_tier_A1[[1]])
})

testthat::test_that("Tier A2 is FDR support + high kME + NOT top-10", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())

  # row 2: clean, FDR-supported, kME 0.90, not top-10 -> A2
  testthat::expect_true(flagged$is_tier_A2[[2]])
  # row 1 is top-10, so it is A1 and must NOT also be A2
  testthat::expect_false(flagged$is_tier_A2[[1]])
  # row 3 fails clean mapping
  testthat::expect_false(flagged$is_tier_A2[[3]])
  # row 4 fails high kME (belongs to D)
  testthat::expect_false(flagged$is_tier_A2[[4]])
  # row 5 has no FDR support
  testthat::expect_false(flagged$is_tier_A2[[5]])

  # A2 does not require top-25 membership either way
  testthat::expect_false("is_top_hub_25" %in%
    names(formals(wcp_assign_candidate_flags)))
})

testthat::test_that("A1 and A2 are mutually exclusive and aggregate to is_tier_A", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  testthat::expect_false(any(flagged$is_tier_A1 & flagged$is_tier_A2))
  testthat::expect_identical(
    flagged$is_tier_A, flagged$is_tier_A1 | flagged$is_tier_A2
  )
})

testthat::test_that("every claimable FDR-supported protein is A1, A2 or D", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  testthat::expect_identical(nrow(wcp_unclassified_fdr_support(flagged)), 0L)

  # each FDR-supported, cleanly mapped row lands in exactly one class
  fdr <- flagged$sus_res_fdr05_any_context & flagged$clean_mapping
  classes <- flagged$is_tier_A1 + flagged$is_tier_A2 + flagged$is_tier_D
  testthat::expect_true(all(classes[fdr] == 1L))
  # and every such row is a candidate
  testthat::expect_true(all(flagged$is_candidate[fdr]))

  # only the mapping contract may leave one unclassified, and it is reported
  unclean <- wcp_fixture_candidates()[3, ]   # FDR-supported, high kME, not clean
  unclean_flagged <- wcp_assign_candidate_flags(unclean)
  testthat::expect_false(unclean_flagged$is_tier_A1)
  testthat::expect_false(unclean_flagged$is_tier_A2)
  testthat::expect_identical(
    nrow(wcp_unclassified_fdr_support(unclean_flagged)), 0L
  )
})

testthat::test_that("phenotype_network_class reflects topology only", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  cls <- wcp_phenotype_network_class(flagged)
  testthat::expect_identical(cls[[1]], "top10_hub")
  testthat::expect_identical(cls[[2]], "module_member")
  testthat::expect_identical(cls[[4]], "peripheral_member")
  # proteins with no phenotype link get no class
  testthat::expect_true(is.na(cls[[5]]))
  testthat::expect_setequal(
    stats::na.omit(unique(cls)),
    c("top10_hub", "module_member", "peripheral_member")
  )
})

testthat::test_that("Tier D admits FDR-supported proteins with low kME", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  # row 4: FDR supported, abs_kME 0.40 -> peripheral candidate, kept
  testthat::expect_true(flagged$is_tier_D[[4]])
  testthat::expect_false(flagged$high_kME[[4]])
  testthat::expect_true(flagged$is_candidate[[4]])
  # a high-kME FDR-supported protein is NOT Tier D
  testthat::expect_false(flagged$is_tier_D[[1]])
})

testthat::test_that("tiers may overlap and the broad label never hides that", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  flagged$candidate_tier <- wcp_candidate_tier(flagged)
  flagged$candidate_tier_all <- wcp_candidate_tier_all(flagged)

  # row 1 is both a phenotype-linked hub and a module representative
  testthat::expect_true(flagged$is_tier_A1[[1]] && flagged$is_tier_C[[1]])
  testthat::expect_identical(flagged$candidate_tier[[1]], "A1")
  testthat::expect_identical(flagged$candidate_tier_all[[1]], "A1;C")
  # row 2 is a phenotype-linked module member
  testthat::expect_identical(flagged$candidate_tier[[2]], "A2")

  # Tier C is purely structural: top-5 regardless of any DA evidence
  structural <- data.frame(
    abs_kME = 0.2, clean_mapping = FALSE,
    is_top5_module_representative = TRUE, is_top10_module_hub = TRUE,
    sus_res_fdr05_any_context = FALSE, sus_res_large_effect_typical = FALSE,
    sus_res_large_effect_any_context = FALSE,
    sus_res_spatially_consistent = FALSE, stringsAsFactors = FALSE
  )
  testthat::expect_true(wcp_assign_candidate_flags(structural)$is_tier_C)

  # a protein satisfying nothing gets no tier at all
  none <- wcp_fixture_candidates()[6, ]
  none$sus_res_large_effect_typical <- FALSE
  none_flagged <- wcp_assign_candidate_flags(none)
  testthat::expect_false(none_flagged$is_candidate)
  testthat::expect_true(is.na(wcp_candidate_tier(none_flagged)))
})

testthat::test_that("candidate_reason lists the flags behind the tier", {
  flagged <- wcp_assign_candidate_flags(wcp_fixture_candidates())
  reason <- wcp_candidate_reason(flagged)
  testthat::expect_match(reason[[1]], "high module membership", fixed = TRUE)
  testthat::expect_match(reason[[1]], "SUS - RES protein-level BH FDR <= 0.05",
                         fixed = TRUE)
  testthat::expect_match(reason[[1]], "top-5 module representative", fixed = TRUE)
})

# ------------------------------------------------------------ determinism

testthat::test_that("flagging and ordering are deterministic", {
  fx <- wcp_fixture_candidates()
  fx$sus_res_max_abs_log2FC <- c(0.4, 0.3, 0.2, 0.9, 0.1, 0.05)
  fx$sus_res_fraction_matching_majority_direction <- 1

  build <- function(d) {
    d <- wcp_assign_candidate_flags(d)
    d$candidate_tier <- wcp_candidate_tier(d)
    wcp_order_candidates(d)
  }
  first <- build(fx)
  second <- build(fx)
  testthat::expect_identical(first, second)

  # row order of the input must not change the exported order
  shuffled <- build(fx[rev(seq_len(nrow(fx))), ])
  testthat::expect_identical(first$ProteinGroupID, shuffled$ProteinGroupID)

  # tier drives the sort, and A1 precedes D
  testthat::expect_identical(first$candidate_tier[[1]], "A1")
  testthat::expect_true(
    which(first$candidate_tier == "A1") < which(first$candidate_tier == "D")
  )
})

# ------------------------------------------------ no upstream recomputation

testthat::test_that("no WGCNA or differential-abundance model is refitted", {
  files <- c(
    repo_path("R", "wgcna_candidate_protein_utils.R"),
    repo_path("analysis/08_integration",
                        "build_candidate_protein_shortlist.R")
  )
  code <- unlist(lapply(files, readLines, warn = FALSE))
  # drop comments so documentation naming these functions cannot fail the test
  code <- code[!grepl("^\\s*#", code)]
  code <- paste(code, collapse = "\n")

  forbidden <- c(
    "blockwiseModules", "TOMsimilarity", "adjacency\\(", "moduleEigengenes",
    "signedKME", "corPvalueStudent", "pickSoftThreshold",
    "lmFit", "eBayes", "topTable", "p\\.adjust", "\\blm\\(", "\\baov\\(",
    "lmer\\(", "emmeans", "t\\.test", "wilcox\\.test", "cor\\.test"
  )
  for (pattern in forbidden) {
    testthat::expect_false(grepl(pattern, code), info = pattern)
  }
})

testthat::test_that("joins are keyed on canonical ProteinGroupID", {
  script <- paste(readLines(repo_path("analysis/08_integration",
    "build_candidate_protein_shortlist.R"
  ), warn = FALSE), collapse = "\n")

  # every join in the script carries ProteinGroupID or is a module-level join
  joins <- regmatches(script, gregexpr('by = c\\([^)]*\\)', script))[[1]]
  testthat::expect_true(length(joins) > 0)
  for (j in joins) {
    testthat::expect_true(
      grepl("ProteinGroupID", j, fixed = TRUE) ||
        grepl("ModuleID", j, fixed = TRUE),
      info = j
    )
  }
  # dataset is always part of the key, so identities never cross datasets
  for (j in joins) {
    testthat::expect_true(grepl('"dataset"', j, fixed = TRUE), info = j)
  }
})

# ------------------------------------------- real-output contract (optional)

wcp_summary_path <- function(scope = "global") {
  path_results("tables", "10_biological_integration",
               "wgcna_candidate_protein_shortlist", scope,
               "wgcna_candidate_proteins_all.csv")
}

testthat::test_that("exported canonical table holds one row per protein identity", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  key <- x[c("dataset", "ModuleID", "ProteinGroupID")]
  testthat::expect_identical(sum(duplicated(key)), 0L)
  # ProteinGroupID is unique within a dataset: one module per protein
  testthat::expect_identical(
    sum(duplicated(x[c("dataset", "ProteinGroupID")])), 0L
  )
  testthat::expect_true(all(nzchar(x$ProteinGroupID)))
  # identities stay dataset-scoped
  testthat::expect_true(all(startsWith(x$ProteinGroupID, paste0("PG:", x$dataset, ":"))))
})

testthat::test_that("exported hub flags match within-module abs_kME ordering", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  recomputed <- wcp_rank_module_members(
    x[c("dataset", "ModuleID", "ProteinGroupID", "abs_kME")]
  )
  idx <- match(
    paste(x$dataset, x$ModuleID, x$ProteinGroupID),
    paste(recomputed$dataset, recomputed$ModuleID, recomputed$ProteinGroupID)
  )
  testthat::expect_identical(
    as.logical(x$is_top5_module_representative),
    recomputed$is_top5_module_representative[idx]
  )
  testthat::expect_identical(
    as.logical(x$is_top10_module_hub), recomputed$is_top10_module_hub[idx]
  )
  # top5 is a strict subset of top10, which is a subset of the frozen top-25
  testthat::expect_true(all(!x$is_top5_module_representative | x$is_top10_module_hub))
  testthat::expect_true(all(!x$is_top10_module_hub | as.logical(x$is_top_hub_25)))

  # exactly 5 and 10 per module wherever the module is large enough
  for (nm in c("is_top5_module_representative", "is_top10_module_hub")) {
    want <- if (nm == "is_top5_module_representative") 5L else 10L
    counts <- tapply(as.logical(x[[nm]]), paste(x$dataset, x$ModuleID), sum)
    sizes <- tapply(rep(1L, nrow(x)), paste(x$dataset, x$ModuleID), sum)
    testthat::expect_true(all(counts == pmin(want, sizes)), info = nm)
  }
})

testthat::test_that("frozen upstream flags are copied, not recomputed", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  # the derived high_kME flag must agree exactly with the frozen core flag
  testthat::expect_identical(
    as.logical(x$high_kME), as.logical(x[["is_core_kME_0.6"]])
  )
  # abs_kME is |kME| as carried through from Stage 01
  finite <- is.finite(x$kME) & is.finite(x$abs_kME)
  testthat::expect_equal(x$abs_kME[finite], abs(x$kME[finite]))
})

testthat::test_that("Tier A rows in the real export satisfy every stated criterion", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  a1 <- x[as.logical(x$is_tier_A1), , drop = FALSE]
  a2 <- x[as.logical(x$is_tier_A2), , drop = FALSE]
  testthat::skip_if(nrow(a1) + nrow(a2) == 0L,
                    "no phenotype-linked module members in current data")

  for (a in list(a1, a2)) {
    if (!nrow(a)) next
    testthat::expect_true(all(as.logical(a$clean_mapping)))
    testthat::expect_true(all(a$abs_kME >= 0.60))
    testthat::expect_true(all(as.logical(a$sus_res_fdr05_any_context)))
    testthat::expect_true(all(a$sus_res_min_BH_FDR <= 0.05))
    # neither class ever contains a gene-ambiguous protein group
    testthat::expect_false(any(
      a$protein_group_ambiguity_class == "multi_gene_indistinguishable"
    ))
  }
  # the split is exactly on top-10 membership
  if (nrow(a1)) testthat::expect_true(all(as.logical(a1$is_top10_module_hub)))
  if (nrow(a2)) testthat::expect_false(any(as.logical(a2$is_top10_module_hub)))
  # A1 and A2 are disjoint, and is_tier_A is their union
  testthat::expect_false(any(as.logical(x$is_tier_A1) & as.logical(x$is_tier_A2)))
  testthat::expect_identical(
    as.logical(x$is_tier_A),
    as.logical(x$is_tier_A1) | as.logical(x$is_tier_A2)
  )
})

testthat::test_that("every claimable FDR-supported protein in the export is classified", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  fdr <- as.logical(x$sus_res_fdr05_any_context)
  claimable <- fdr & as.logical(x$clean_mapping)
  classified <- as.logical(x$is_tier_A1) | as.logical(x$is_tier_A2) |
    as.logical(x$is_tier_D)
  testthat::expect_true(all(classified[claimable]))
  # exactly one phenotype-linked class each
  n_class <- as.integer(x$is_tier_A1 == "TRUE") +
    as.integer(x$is_tier_A2 == "TRUE") + as.integer(x$is_tier_D == "TRUE")
  testthat::expect_true(all(n_class[claimable] == 1L))
  # and All_candidates no longer omits them
  testthat::expect_true(all(as.logical(x$is_candidate)[claimable]))

  # A1 + A2 + D accounts for every FDR-supported protein per dataset
  for (ds in unique(x$dataset)) {
    sub <- x[x$dataset == ds, , drop = FALSE]
    testthat::expect_identical(
      sum(as.logical(sub$sus_res_fdr05_any_context)),
      sum(as.logical(sub$is_tier_A1)) + sum(as.logical(sub$is_tier_A2)) +
        sum(as.logical(sub$is_tier_D)),
      info = ds
    )
  }
})

testthat::test_that("Tier D rows are FDR-supported and peripheral", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- x[as.logical(x$is_tier_D), , drop = FALSE]
  testthat::skip_if(nrow(d) == 0L, "no Tier D candidates in current data")

  testthat::expect_true(all(as.logical(d$sus_res_fdr05_any_context)))
  testthat::expect_true(all(d$abs_kME < 0.60))
  testthat::expect_false(any(as.logical(d$is_tier_A)))
})

testthat::test_that("top10-per-module export is exactly the top 10 by abs_kME", {
  path <- path_results(
    "tables", "10_biological_integration",
    "wgcna_candidate_protein_shortlist", "global",
    "wgcna_top10_per_module.csv"
  )
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  top10 <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  testthat::expect_true(all(top10$abs_kME_rank_in_module <= 10L))
  counts <- tapply(rep(1L, nrow(top10)),
                   paste(top10$dataset, top10$ModuleID), sum)
  testthat::expect_true(all(counts <= 10L))
  # ranks within a module are the contiguous run 1..k
  ranks <- split(top10$abs_kME_rank_in_module,
                 paste(top10$dataset, top10$ModuleID))
  for (nm in names(ranks)) {
    testthat::expect_identical(sort(ranks[[nm]]), seq_along(ranks[[nm]]),
                               info = nm)
  }
})

testthat::test_that("direction concordance in the real export is never fabricated", {
  path <- wcp_summary_path("global")
  testthat::skip_if_not(file.exists(path), "shortlist has not been generated")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  # wherever concordance is non-NA, both inputs must be finite and non-zero
  decided <- !is.na(x$sus_res_direction_matches_module)
  testthat::expect_true(all(is.finite(x$module_estimate[decided])))
  testthat::expect_true(all(x$module_estimate[decided] != 0))
  testthat::expect_true(all(is.finite(x$sus_res_median_log2FC[decided])))
  testthat::expect_true(all(x$sus_res_median_log2FC[decided] != 0))
  # and it must equal plain sign agreement
  testthat::expect_identical(
    as.logical(x$sus_res_direction_matches_module[decided]),
    sign(x$sus_res_median_log2FC[decided]) == sign(x$module_estimate[decided])
  )
})

testthat::test_that("dry-run reports inputs and writes nothing", {
  script <- repo_path("analysis/08_integration",
                      "build_candidate_protein_shortlist.R")
  testthat::skip_if_not(file.exists(script))

  out_root <- path_results("tables", "10_biological_integration",
                           "wgcna_candidate_protein_shortlist")
  before <- if (dir.exists(out_root)) {
    f <- list.files(out_root, recursive = TRUE, full.names = TRUE)
    stats::setNames(file.info(f)$mtime, f)
  } else NULL

  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("analysis/08_integration/build_candidate_protein_shortlist.R",
      "--dry-run"),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  testthat::expect_true(status %in% c(0L, 1L))
  testthat::expect_true(any(grepl("[DRY-RUN", out, fixed = TRUE)))

  # nothing created, nothing rewritten
  after <- if (dir.exists(out_root)) {
    f <- list.files(out_root, recursive = TRUE, full.names = TRUE)
    stats::setNames(file.info(f)$mtime, f)
  } else NULL
  testthat::expect_identical(names(before), names(after))
  if (!is.null(before)) testthat::expect_identical(before, after)
})

# --------------------------------------------------- revised Tier B contract

wcp_tierb_row <- function(large, consistent, kme = 0.90) {
  data.frame(
    dataset = "d", ModuleID = "M", ProteinGroupID = "PG:d:01",
    abs_kME = kme, clean_mapping = TRUE,
    is_top5_module_representative = FALSE, is_top10_module_hub = FALSE,
    sus_res_fdr05_any_context = FALSE,
    sus_res_large_effect_typical = large,
    sus_res_large_effect_any_context = TRUE,
    sus_res_spatially_consistent = consistent,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("Tier B requires high kME and a large typical effect", {
  testthat::expect_true(
    wcp_assign_candidate_flags(wcp_tierb_row(TRUE, FALSE))$is_tier_B
  )
  # high kME alone is not enough
  testthat::expect_false(
    wcp_assign_candidate_flags(wcp_tierb_row(FALSE, FALSE))$is_tier_B
  )
  # a large typical effect on a peripheral protein is not Tier B
  testthat::expect_false(
    wcp_assign_candidate_flags(wcp_tierb_row(TRUE, FALSE, kme = 0.59))$is_tier_B
  )
})

testthat::test_that("spatial consistency alone can never create Tier B", {
  # unanimous spatial direction, no large typical effect -> not a candidate
  consistent_only <- wcp_assign_candidate_flags(wcp_tierb_row(FALSE, TRUE))
  testthat::expect_false(consistent_only$is_tier_B)
  testthat::expect_false(consistent_only$is_candidate)

  # adding consistency to a Tier B protein must not change its tier membership
  without <- wcp_assign_candidate_flags(wcp_tierb_row(TRUE, FALSE))
  with <- wcp_assign_candidate_flags(wcp_tierb_row(TRUE, TRUE))
  testthat::expect_identical(without$is_tier_B, with$is_tier_B)

  # the per-context effect flag is likewise not a route into Tier B
  per_context_only <- wcp_tierb_row(FALSE, FALSE)
  per_context_only$sus_res_large_effect_any_context <- TRUE
  testthat::expect_false(wcp_assign_candidate_flags(per_context_only)$is_tier_B)
})

testthat::test_that("descriptive spatial fields are still exported", {
  flagged <- wcp_assign_candidate_flags(wcp_tierb_row(TRUE, TRUE))
  for (nm in c("sus_res_spatially_consistent", "sus_res_large_effect_typical",
               "sus_res_large_effect_any_context")) {
    testthat::expect_true(nm %in% names(flagged), info = nm)
  }
})

# ------------------------------------------------------- review sheet views

testthat::test_that("Protein_review is compact and free of provenance noise", {
  cols <- wcp_protein_review_columns()
  testthat::expect_gte(length(cols), 25L)
  testthat::expect_lte(length(cols), 32L)

  # no source paths, contract fields or internal plumbing
  testthat::expect_false(any(grepl("source_file|source_key|contract_version|^Source$",
                                   cols)))
  # no redundant variants of the same quantity
  testthat::expect_false(any(c("top5_hub", "top10_hub") %in% cols))
  testthat::expect_false(any(c("module_display_label", "final_label",
                               "primary_label", "GeneSymbols",
                               "representative_gene_symbol") %in% cols))
  testthat::expect_false("is_core_kME_0.6" %in% cols)
  testthat::expect_false("sus_res_n_spatial_contexts_fdr05" %in% cols)
  # the fields a reviewer actually needs are present
  for (nm in c("GeneSymbol", "abs_kME", "sus_res_min_BH_FDR",
               "sus_res_strongest_spatial_unit", "module_estimate",
               "candidate_reason", "candidate_tier")) {
    testthat::expect_true(nm %in% cols, info = nm)
  }
})

testthat::test_that("Protein_review ordering puts Tier A first and is deterministic", {
  fx <- data.frame(
    dataset = "d", ModuleID = "M",
    ProteinGroupID = sprintf("PG:d:%02d", 1:4),
    candidate_tier = c("C", "B", "D", "A1"),
    sus_res_fdr05_any_context = c(FALSE, FALSE, TRUE, TRUE),
    sus_res_min_BH_FDR = c(NA, NA, 0.01, 0.02),
    sus_res_large_effect_typical = c(FALSE, TRUE, FALSE, TRUE),
    sus_res_median_abs_log2FC = c(0.1, 0.5, 0.9, 0.7),
    abs_kME = c(0.95, 0.80, 0.40, 0.85),
    stringsAsFactors = FALSE
  )
  ordered <- wcp_order_protein_review(fx)
  testthat::expect_identical(ordered$candidate_tier, c("A1", "D", "B", "C"))
  # stable under input permutation
  testthat::expect_identical(
    wcp_order_protein_review(fx[c(3L, 1L, 4L, 2L), ])$ProteinGroupID,
    ordered$ProteinGroupID
  )
})

testthat::test_that("SUS_RES_FDR_hits is exactly the FDR-supported proteins", {
  fx <- data.frame(
    dataset = "d", ModuleID = "M",
    ProteinGroupID = sprintf("PG:d:%02d", 1:5),
    candidate_tier = c("A", NA, "D", "B", NA),
    sus_res_fdr05_any_context = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    sus_res_min_BH_FDR = c(0.01, 0.02, 0.03, 0.5, NA),
    sus_res_large_effect_typical = FALSE,
    sus_res_median_abs_log2FC = 0.2,
    abs_kME = c(0.9, 0.7, 0.4, 0.9, 0.3),
    GeneSymbol = paste0("G", 1:5),
    stringsAsFactors = FALSE
  )
  hits <- wcp_fdr_hits_table(fx, rename = FALSE)
  testthat::expect_identical(nrow(hits), 3L)
  # includes FDR-supported proteins that carry NO tier at all
  testthat::expect_true(any(is.na(hits$candidate_tier)))
  testthat::expect_setequal(hits$GeneSymbol, c("G1", "G2", "G3"))
})

testthat::test_that("Module_review is one descriptive row per module", {
  fx <- data.frame(
    dataset = "d", ModuleID = c("M1", "M1", "M2"),
    ProteinGroupID = sprintf("PG:d:%02d", 1:3),
    module_label = c("L1", "L1", "L2"),
    abs_kME = c(0.9, 0.5, 0.8), GeneSymbol = c("A", "B", "C"),
    is_tier_A1 = c(TRUE, FALSE, FALSE), is_tier_A2 = FALSE, is_tier_B = FALSE,
    is_tier_C = c(TRUE, FALSE, TRUE), is_tier_D = FALSE,
    is_candidate = c(TRUE, FALSE, TRUE),
    sus_res_fdr05_any_context = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  mr <- wcp_module_review_table(fx)
  testthat::expect_identical(nrow(mr), 2L)
  testthat::expect_identical(mr$Module, c("M1", "M2"))
  testthat::expect_identical(mr$`Module size`, c(2L, 1L))
  testthat::expect_identical(mr$Candidates, c(1L, 1L))
  testthat::expect_equal(mr$`Candidate fraction`, c(0.5, 1))
  testthat::expect_identical(mr$`SUS-RES FDR<=0.05`, c(1L, 0L))
})

# --------------------------------------------------- xlsx package integrity

testthat::test_that("openxlsx alone emits an invalid package that repair fixes", {
  testthat::skip_if_not_installed("openxlsx")
  testthat::skip_if_not_installed("zip")
  source(repo_path("R", "xlsx_package_utils.R"))

  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "S1")
  openxlsx::writeData(wb, "S1", data.frame(a = 1:3))
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  # the defect is in the writer, not in anything this repo does
  testthat::expect_gt(nrow(xlsx_package_dangling_relationships(path)), 0L)
  testthat::expect_false(xlsx_package_is_valid(path))

  before <- utils::unzip(path, list = TRUE)$Name
  xlsx_repair_package(path)
  testthat::expect_true(xlsx_package_is_valid(path))
  # repair removes plumbing only; no content part is lost
  testthat::expect_identical(sort(utils::unzip(path, list = TRUE)$Name),
                             sort(before))
  testthat::expect_equal(
    as.data.frame(readxl::read_excel(path))$a, 1:3
  )
})

testthat::test_that("xlsx_save_valid_workbook writes a well-formed package", {
  testthat::skip_if_not_installed("openxlsx")
  testthat::skip_if_not_installed("zip")
  source(repo_path("R", "xlsx_package_utils.R"))

  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "S1")
  openxlsx::writeData(wb, "S1", data.frame(a = 1:3))
  openxlsx::freezePane(wb, "S1", firstActiveRow = 2)
  xlsx_save_valid_workbook(wb, path)

  testthat::expect_true(xlsx_package_is_valid(path))
  testthat::expect_silent(xlsx_assert_package_valid(path))
  # every internal relationship resolves to a part that is really present
  rel <- xlsx_package_relationships(path)
  testthat::expect_true(all(rel$exists))
})

# --------------------------------------- generated workbook contract (real)

wcp_workbook_path <- function(scope) {
  path_results("tables", "10_biological_integration",
               "wgcna_candidate_protein_shortlist", scope,
               "wgcna_candidate_proteins_shortlist.xlsx")
}

testthat::test_that("generated workbooks are valid OOXML packages", {
  for (scope in c("neuron_neuropil", "neuron_soma", "microglia", "global")) {
    path <- wcp_workbook_path(scope)
    if (!file.exists(path)) next
    source(repo_path("R", "xlsx_package_utils.R"))
    testthat::expect_identical(
      nrow(xlsx_package_dangling_relationships(path)), 0L, info = scope
    )
    testthat::expect_identical(
      length(xlsx_package_orphan_overrides(path)), 0L, info = scope
    )
    # readable by a parser independent of the writer
    testthat::expect_silent(readxl::excel_sheets(path))
  }
})

testthat::test_that("workbook sheet sets are deterministic and non-redundant", {
  single <- wcp_workbook_path("neuron_neuropil")
  global <- wcp_workbook_path("global")
  testthat::skip_if_not(file.exists(single) && file.exists(global),
                        "workbooks have not been generated")

  expected_core <- c("README", "Protein_review", "SUS_RES_FDR_hits",
                     "Module_review", "All_candidates", "Tier_A1", "Tier_A2",
                     "Tier_B", "Tier_C", "Tier_D", "Top10_per_module")

  single_sheets <- readxl::excel_sheets(single)
  testthat::expect_identical(single_sheets, expected_core)
  # a single-dataset workbook must not duplicate itself as DS_<dataset>
  testthat::expect_false(any(grepl("^DS_", single_sheets)))

  global_sheets <- readxl::excel_sheets(global)
  testthat::expect_identical(global_sheets[seq_along(expected_core)],
                             expected_core)
  # the cross-dataset workbook keeps its per-dataset sheets
  testthat::expect_setequal(
    grep("^DS_", global_sheets, value = TRUE),
    c("DS_neuron_neuropil", "DS_neuron_soma", "DS_microglia")
  )
})

testthat::test_that("review sheets carry the intended compact content", {
  path <- wcp_workbook_path("neuron_neuropil")
  testthat::skip_if_not(file.exists(path), "workbook has not been generated")

  review <- readxl::read_excel(path, sheet = "Protein_review", skip = 3,
                               .name_repair = "minimal")
  testthat::expect_lte(ncol(review), 32L)
  testthat::expect_gte(ncol(review), 25L)
  testthat::expect_false(any(grepl("source_file|source_key|contract_version",
                                   names(review))))

  hits <- readxl::read_excel(path, sheet = "SUS_RES_FDR_hits", skip = 3,
                             .name_repair = "minimal")
  testthat::expect_identical(ncol(hits), ncol(review))

  csv <- path_results("tables", "10_biological_integration",
                      "wgcna_candidate_protein_shortlist", "neuron_neuropil",
                      "wgcna_candidate_proteins_all.csv")
  if (file.exists(csv)) {
    x <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
    # exactly every FDR-supported protein, tiered or not, with no duplicates
    testthat::expect_identical(
      nrow(hits), sum(as.logical(x$sus_res_fdr05_any_context))
    )
    testthat::expect_identical(anyDuplicated(hits[["Gene"]]), 0L)
  }
})

# ------------------------------------------- zero-row / optional-column safety

testthat::test_that("ordering and review views survive zero rows and missing optional columns", {
  # Regression guard. The canonical %||% also falls back for length-0 input, so
  # `df$col %||% NA` returns a length-1 value for a zero-row frame and order()
  # then aborts with "argument lengths differ". These views must stay total.
  minimal <- data.frame(
    candidate_tier = c("A1", "D", "B"), abs_kME = c(0.9, 0.3, 0.7),
    ProteinGroupID = c("P1", "P2", "P3"), stringsAsFactors = FALSE
  )
  testthat::expect_silent(wcp_order_candidates(minimal))
  testthat::expect_identical(
    wcp_order_protein_review(minimal)$candidate_tier, c("A1", "D", "B")
  )

  empty <- minimal[0, , drop = FALSE]
  testthat::expect_identical(nrow(wcp_order_candidates(empty)), 0L)
  testthat::expect_identical(nrow(wcp_order_protein_review(empty)), 0L)

  full_empty <- data.frame(
    dataset = character(), ModuleID = character(), ProteinGroupID = character(),
    candidate_tier = character(), sus_res_fdr05_any_context = logical(),
    sus_res_min_BH_FDR = numeric(), sus_res_large_effect_typical = logical(),
    sus_res_median_abs_log2FC = numeric(), abs_kME = numeric(),
    GeneSymbol = character(), stringsAsFactors = FALSE
  )
  testthat::expect_identical(nrow(wcp_protein_review_table(full_empty)), 0L)
  testthat::expect_identical(nrow(wcp_fdr_hits_table(full_empty)), 0L)
})

testthat::test_that("the phenotype-linked row highlight targets A1/A2, not the retired A", {
  script <- paste(readLines(repo_path("analysis/08_integration",
    "build_candidate_protein_shortlist.R"
  ), warn = FALSE), collapse = "\n")
  # `tier %in% "A"` would be dead code: wcp_candidate_tier never emits "A".
  testthat::expect_false(grepl('tier %in% "A"', script, fixed = TRUE))
  testthat::expect_true(grepl('tier %in% c("A1", "A2")', script, fixed = TRUE))
  testthat::expect_false("A" %in% wcp_candidate_tier_levels())
})

testthat::test_that("empty tier sheets contain no phantom record", {
  # openxlsx::writeDataTable on a zero-row frame spans one data row, which reads
  # back as a single all-NA record and makes an empty class look populated.
  path <- wcp_workbook_path("neuron_soma")
  testthat::skip_if_not(file.exists(path), "workbook has not been generated")
  csv <- path_results("tables", "10_biological_integration",
                      "wgcna_candidate_protein_shortlist", "neuron_soma",
                      "wgcna_candidate_proteins_all.csv")
  testthat::skip_if_not(file.exists(csv))
  x <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)

  for (nm in c("Tier_A1", "Tier_A2", "Tier_B", "Tier_C", "Tier_D")) {
    flag <- sub("Tier_", "is_tier_", nm)
    sheet <- readxl::read_excel(path, sheet = nm, skip = 3,
                                .name_repair = "minimal")
    testthat::expect_identical(nrow(sheet), sum(as.logical(x[[flag]])),
                               info = nm)
  }
  # an empty sheet says so, above the header, rather than showing a blank row
  note <- readxl::read_excel(path, sheet = "Tier_A1", range = "A2",
                             col_names = FALSE, .name_repair = "minimal")
  testthat::expect_match(as.character(note[[1]][1]), "NONE IN THIS SCOPE",
                         fixed = TRUE)
})
