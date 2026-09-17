source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "wgcna_network_position_utils.R"))
source(repo_path("R", "wgcna_label_coherence_utils.R"))

# =====================================================================
# A. SUS - RES NETWORK POSITION (phenotype-aware)
# =====================================================================

wnp_fixture <- function(n_per_module = 100L, hits = c(1L, 50L, 99L)) {
  mods <- c("WGCNA_m01", "WGCNA_m02")
  out <- do.call(rbind, lapply(mods, function(m) {
    data.frame(
      dataset = "d", ModuleID = m,
      ProteinGroupID = sprintf("PG:d:%s:%03d", m, seq_len(n_per_module)),
      abs_kME = seq(0.99, 0.10, length.out = n_per_module),
      abs_kME_rank_in_module = seq_len(n_per_module),
      n_module_members = n_per_module,
      stringsAsFactors = FALSE
    )
  }))
  out$rank_fraction <- wnp_rank_fraction(out$abs_kME_rank_in_module,
                                         out$n_module_members)
  out[["is_core_kME_0.6"]] <- out$abs_kME >= 0.60
  out$is_top_hub_25 <- out$abs_kME_rank_in_module <= 25L
  out$is_top10_module_hub <- out$abs_kME_rank_in_module <= 10L
  out$is_hit <- FALSE
  for (m in mods) out$is_hit[out$ModuleID == m][hits] <- TRUE
  out
}

testthat::test_that("rank fraction is a midpoint fraction centred on 0.5", {
  testthat::expect_equal(wnp_rank_fraction(1L, 100L), 0.005)
  testthat::expect_equal(wnp_rank_fraction(100L, 100L), 0.995)
  # a uniformly spread set averages exactly 0.5 regardless of module size
  testthat::expect_equal(mean(wnp_rank_fraction(1:10, 10L)), 0.5)
  testthat::expect_equal(mean(wnp_rank_fraction(1:1000, 1000L)), 0.5)
  testthat::expect_true(is.na(wnp_rank_fraction(1L, 0L)))
})

testthat::test_that("the permutation preserves per-module hit counts", {
  fx <- wnp_fixture()
  module_rows <- split(seq_len(nrow(fx)), fx$ModuleID)
  hits_per_module <- vapply(module_rows, function(r) sum(fx$is_hit[r]), integer(1))
  set.seed(1)
  for (i in 1:20) {
    drawn <- wnp_draw_null_rows(module_rows, hits_per_module)
    per_mod <- table(fx$ModuleID[drawn])
    testthat::expect_identical(as.integer(per_mod[names(hits_per_module)]),
                               as.integer(hits_per_module))
    # never draws the same protein twice within one permutation
    testthat::expect_identical(anyDuplicated(drawn), 0L)
  }
})

testthat::test_that("the null is centred on 0.5 and detects a truly central set", {
  # hits placed at ranks 1..5 of each module are extremely central
  fx <- wnp_fixture(hits = 1:5)
  perm <- wnp_permute(fx, n_permutations = 2000L, seed = 7L)
  s <- wnp_permutation_summary(perm, dataset = "d", scope = "test")

  med <- s[s$statistic == "median_rank_fraction", ]
  testthat::expect_lt(med$observed, 0.1)
  testthat::expect_gt(med$null_median, 0.4)
  testthat::expect_lt(med$null_median, 0.6)
  # detected as more central, not more peripheral
  testthat::expect_lt(med$p_more_central, 0.01)
  testthat::expect_gt(med$p_more_peripheral, 0.5)
  testthat::expect_identical(wnp_interpret(med), "more_central_than_null")
})

testthat::test_that("a truly peripheral set is detected in the other direction", {
  fx <- wnp_fixture(hits = 96:100)
  perm <- wnp_permute(fx, n_permutations = 2000L, seed = 7L)
  s <- wnp_permutation_summary(perm, dataset = "d", scope = "test")
  med <- s[s$statistic == "median_rank_fraction", ]
  testthat::expect_gt(med$observed, 0.9)
  testthat::expect_lt(med$p_more_peripheral, 0.01)
  testthat::expect_identical(wnp_interpret(med), "more_peripheral_than_null")

  # |kME| is oriented the opposite way and must agree with the same conclusion
  kme <- s[s$statistic == "median_abs_kME", ]
  testthat::expect_identical(unname(kme$central_direction), "high")
  testthat::expect_lt(kme$p_more_peripheral, 0.01)
})

testthat::test_that("a uniformly spread set is read as consistent with the null", {
  fx <- wnp_fixture(hits = c(10L, 30L, 50L, 70L, 90L))
  perm <- wnp_permute(fx, n_permutations = 2000L, seed = 7L)
  s <- wnp_permutation_summary(perm, dataset = "d", scope = "test")
  testthat::expect_true(all(wnp_interpret(s) == "consistent_with_null"))
})

testthat::test_that("permutation p-values use the +1 correction and are bounded", {
  fx <- wnp_fixture(hits = 1:5)
  perm <- wnp_permute(fx, n_permutations = 500L, seed = 3L)
  s <- wnp_permutation_summary(perm, dataset = "d", scope = "test")
  # smallest attainable p is 1/(B+1), never 0
  testthat::expect_true(all(s$p_more_central >= 1 / 501))
  testthat::expect_true(all(s$p_more_peripheral >= 1 / 501))
  testthat::expect_true(all(s$p_more_central <= 1))
  testthat::expect_true(all(s$p_two_sided <= 1))
  testthat::expect_identical(unique(s$n_permutations), 500L)
})

testthat::test_that("the permutation is deterministic for a fixed seed", {
  fx <- wnp_fixture()
  a <- wnp_permutation_summary(wnp_permute(fx, 500L, seed = 11L), "d", "t")
  b <- wnp_permutation_summary(wnp_permute(fx, 500L, seed = 11L), "d", "t")
  testthat::expect_identical(a, b)
  # and it does not disturb the caller's RNG stream
  set.seed(99); before <- runif(3)
  set.seed(99); invisible(wnp_permute(fx, 50L, seed = 5L)); after <- runif(3)
  testthat::expect_identical(before, after)
})

testthat::test_that("module summary counts each module once and BH-corrects within dataset", {
  fx <- wnp_fixture()
  fx$is_tier_A1 <- FALSE; fx$is_tier_A2 <- FALSE; fx$is_tier_D <- fx$is_hit
  m <- wnp_module_summary(fx)
  testthat::expect_identical(nrow(m), 2L)
  testthat::expect_identical(sum(m$n_hits), sum(fx$is_hit))
  testthat::expect_identical(sum(m$n_da_eligible), nrow(fx))
  testthat::expect_true(all(m$exploratory_hypergeometric_BH >=
                              m$exploratory_hypergeometric_p, na.rm = TRUE))
  testthat::expect_true(all(m$hit_fraction == m$n_hits / m$n_da_eligible))
})

testthat::test_that("network-position helpers never recompute WGCNA or DA", {
  code <- readLines(repo_path("R", "wgcna_network_position_utils.R"), warn = FALSE)
  code <- paste(code[!grepl("^\\s*#", code)], collapse = "\n")
  for (pattern in c("blockwiseModules", "signedKME", "moduleEigengenes", "lmFit",
                    "eBayes", "topTable", "\\blm\\(", "p\\.adjust\\(.*BH.*log2")) {
    testthat::expect_false(grepl(pattern, code), info = pattern)
  }
})

# =====================================================================
# B. LABEL COHERENCE (phenotype-blind)
# =====================================================================

wcl_fixture_go <- function() {
  data.frame(
    ModuleID = c("WGCNA_m01", "WGCNA_m01", "WGCNA_m01", "WGCNA_m02"),
    Ontology = "BP",
    ModuleProteinSetType = c("all", "all", "core_kME_0.6", "all"),
    ID = c("GO:1", "GO:2", "GO:3", "GO:4"),
    Description = c("synapse organization", "endosome transport",
                    "core only term", "protein phosphorylation"),
    p.adjust = c(0.001, 0.001, 1e-9, 0.02),
    qvalue = c(0.0005, 0.0009, 1e-9, 0.01),
    GeneRatio = c("3/10", "2/10", "5/10", "2/10"),
    Count = c(3L, 2L, 5L, 2L),
    geneID = c("11/22/33", "22/44", "99", "22/55"),
    stringsAsFactors = FALSE
  )
}

wcl_fixture_members <- function() {
  n <- 10L
  out <- data.frame(
    dataset = "d", ModuleID = "WGCNA_m01",
    ProteinGroupID = sprintf("PG:d:%02d", 1:n),
    GeneSymbol = paste0("G", 1:n),
    EntrezID = as.character(c(11, 22, 33, 44, 55, 66, 77, 88, 99, 100)),
    abs_kME = seq(0.95, 0.20, length.out = n),
    abs_kME_rank_in_module = 1:n, n_module_members = n,
    stringsAsFactors = FALSE
  )
  out$rank_fraction <- (out$abs_kME_rank_in_module - 0.5) / n
  out[["is_core_kME_0.6"]] <- out$abs_kME >= 0.60
  out$is_top_hub_25 <- TRUE
  out$is_top10_module_hub <- out$abs_kME_rank_in_module <= 10L
  out$is_top5_module_representative <- out$abs_kME_rank_in_module <= 5L
  out
}

testthat::test_that("the label term reproduces Stage-01 selection exactly", {
  lt <- wcl_label_defining_terms(wcl_fixture_go(), "BP", "all")
  testthat::expect_identical(nrow(lt), 2L)
  # restricted to the "all" protein set: the core-only term must not win
  testthat::expect_false("core only term" %in% lt$label_go_description)
  # lowest p.adjust, ties broken by qvalue
  m1 <- lt[lt$ModuleID == "WGCNA_m01", ]
  testthat::expect_identical(m1$label_go_description, "synapse organization")
  testthat::expect_identical(m1$label_go_n_terms_tied_at_min, 2L)
  testthat::expect_true(m1$label_go_tie_broken_arbitrarily)
  # no significance threshold is applied, but non-significance is flagged
  testthat::expect_true(m1$label_go_is_fdr_significant)
})

testthat::test_that("a non-significant label term is flagged, not dropped", {
  go <- wcl_fixture_go()
  go$p.adjust[go$ModuleID == "WGCNA_m02"] <- 0.4
  lt <- wcl_label_defining_terms(go, "BP", "all")
  m2 <- lt[lt$ModuleID == "WGCNA_m02", ]
  testthat::expect_identical(m2$label_go_description, "protein phosphorylation")
  testthat::expect_false(m2$label_go_is_fdr_significant)
})

testthat::test_that("GO contributors map to ProteinGroupIDs and report losses", {
  members <- wcl_fixture_members()
  mapped <- wcl_map_label_contributors("11/22/33/12345", members)
  testthat::expect_identical(mapped$n_terms, 4L)
  testthat::expect_identical(mapped$n_mapped_unique, 3L)
  testthat::expect_identical(mapped$n_unmapped, 1L)      # 12345 is not a member
  testthat::expect_setequal(members$ProteinGroupID[mapped$rows],
                            c("PG:d:01", "PG:d:02", "PG:d:03"))
  # an Entrez ID matching several protein groups is counted, never silently kept
  dup <- members; dup$EntrezID[[4]] <- "11"
  m2 <- wcl_map_label_contributors("11", dup)
  testthat::expect_identical(m2$n_mapped_multiple, 1L)
  testthat::expect_identical(m2$n_mapped_unique, 0L)
})

testthat::test_that("centrality AUC is 0.5 for a spread set and oriented correctly", {
  rf <- (1:100 - 0.5) / 100
  central <- rep(FALSE, 100); central[1:10] <- TRUE
  peripheral <- rep(FALSE, 100); peripheral[91:100] <- TRUE
  spread <- rep(FALSE, 100); spread[seq(5, 95, by = 10)] <- TRUE

  testthat::expect_gt(wcl_centrality_auc(rf, central), 0.9)
  testthat::expect_lt(wcl_centrality_auc(rf, peripheral), 0.1)
  testthat::expect_equal(wcl_centrality_auc(rf, spread), 0.5, tolerance = 0.05)
  testthat::expect_true(is.na(wcl_centrality_auc(rf, rep(FALSE, 100))))
})

testthat::test_that("module evidence has one row per module and correct hub counts", {
  members <- wcl_fixture_members()
  lt <- wcl_label_defining_terms(wcl_fixture_go(), "BP", "all")
  ev <- wcl_module_label_evidence(members, lt)
  testthat::expect_identical(nrow(ev), 1L)
  testthat::expect_identical(ev$ModuleID, "WGCNA_m01")
  testthat::expect_identical(ev$module_size, 10L)
  # label term genes 11/22/33 are ranks 1,2,3 -> all top-5 and all top-10
  testthat::expect_identical(ev$n_contributors_in_top5, 3L)
  testthat::expect_identical(ev$n_contributors_in_top10, 3L)
  testthat::expect_equal(ev$fraction_top10_hubs_contributing, 0.3)
  testthat::expect_gt(ev$contributor_centrality_auc, 0.5)
})

testthat::test_that("top-10/top-25 hub flags follow the frozen kME ordering", {
  members <- wcl_fixture_members()
  lt <- wcl_label_defining_terms(wcl_fixture_go(), "BP", "all")
  ev <- wcl_module_label_evidence(members, lt)
  syms <- trimws(strsplit(ev$top10_hub_symbols, ";", fixed = TRUE)[[1]])
  expected <- members$GeneSymbol[order(members$abs_kME_rank_in_module)][1:10]
  testthat::expect_identical(syms, expected)
})

testthat::test_that("coherence classes and actions behave as specified", {
  base <- data.frame(
    label_go_description = c("synapse organization", "protein phosphorylation",
                             "vesicle fusion", "keratinocyte differentiation",
                             "proteasome assembly"),
    label_go_is_fdr_significant = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    contributor_centrality_auc = c(0.70, 0.60, 0.30, 0.72, 0.30),
    fraction_top10_hubs_contributing = c(0.7, 0.1, 0.0, 0.2, 0.0),
    median_contributor_rank_fraction = 0.3,
    label_go_tie_broken_arbitrarily = FALSE,
    stringsAsFactors = FALSE
  )
  out <- wcl_classify_module_labels(base)
  testthat::expect_identical(out$coherence_class[1], "strongly_supported")
  testthat::expect_identical(out$review_action[1], "KEEP")
  # generic wording overrides otherwise adequate support
  testthat::expect_identical(out$coherence_class[2], "too_generic")
  testthat::expect_identical(out$review_action[2], "RENAME_REVIEW")
  # significant but peripherally driven -> the ORA stands, the NAME does not
  testthat::expect_identical(out$coherence_class[3], "peripheral_GO_driven")
  testthat::expect_identical(out$review_action[3], "RENAME_REVIEW")
  # tissue-implausible but well supported -> contextualize, never "wrong"
  testthat::expect_identical(out$coherence_class[4],
                             "contextually_implausible_without_caveat")
  testthat::expect_identical(out$review_action[4], "CONTEXTUALIZE")
  # not FDR-significant -> never KEEP
  testthat::expect_false(out$review_action[5] == "KEEP")
  testthat::expect_true(all(out$coherence_class %in% wcl_module_coherence_classes()))
  testthat::expect_true(all(out$review_action %in% wcl_module_actions()))
})

testthat::test_that("the phenotype-blindness guard rejects outcome-derived fields", {
  ok <- data.frame(ModuleID = "WGCNA_m01", label_go_id = "GO:1",
                   contributor_centrality_auc = 0.7, stringsAsFactors = FALSE)
  testthat::expect_silent(wcl_assert_phenotype_blind(ok))

  for (bad_col in c("sus_res_min_BH_FDR", "log2FC", "candidate_tier",
                    "is_tier_A1", "group_effect_direction", "module_estimate",
                    "contrast")) {
    bad <- ok; bad[[bad_col]] <- 1
    testthat::expect_error(wcl_assert_phenotype_blind(bad),
                           "phenotype-derived", info = bad_col)
  }
  # the attestation column itself is allowed through explicitly
  att <- ok; att$naming_evidence_is_phenotype_blind <- TRUE
  testthat::expect_silent(
    wcl_assert_phenotype_blind(att, allow = "naming_evidence_is_phenotype_blind")
  )
})

testthat::test_that("supermodule evidence is module-balanced and one row per supermodule", {
  member_map <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM01", "SM02"),
    ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"), stringsAsFactors = FALSE
  )
  module_evidence <- data.frame(
    dataset = "d", ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"),
    module_size = c(1400L, 40L, 200L),          # deliberately unbalanced
    proposed_module_theme = c("theme A", "theme B", "theme C"),
    current_stage01_label = c("L1", "L2", "L3"),
    proposed_label = c("P1", "P2", "P3"),
    review_action = c("KEEP", "KEEP", "REFINE_WORDING"),
    top10_hub_symbols = c("A1; A2; A3; A4", "B1; B2; B3; B4", "C1; C2; C3"),
    stringsAsFactors = FALSE
  )
  structural <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM02"),
    adjusted_signed_min_pairwise_eigengene_correlation = c(0.7, NA),
    pc1_variance_explained = c(0.8, 1), cut_height_stability_fraction_stable = c(0.9, 1),
    stringsAsFactors = FALSE
  )
  sm <- wcl_supermodule_evidence(member_map, module_evidence, structural,
                                 hubs_per_module = 2L)
  testthat::expect_identical(nrow(sm), 2L)
  testthat::expect_identical(sm$SupermoduleID, c("SM01", "SM02"))

  # equal weighting: the 1400-member module does not outvote the 40-member one
  sm01 <- sm[sm$SupermoduleID == "SM01", ]
  testthat::expect_identical(sm01$n_modules_supporting_dominant_theme, 1L)
  testthat::expect_equal(sm01$dominant_theme_fraction, 0.5)
  testthat::expect_identical(sm01$biological_coherence_class,
                             "mixed_but_structurally_coherent")

  # the illustrative panel takes N hubs from EVERY member module
  testthat::expect_match(sm01$balanced_hub_panel, "WGCNA_m01:A1,A2", fixed = TRUE)
  testthat::expect_match(sm01$balanced_hub_panel, "WGCNA_m02:B1,B2", fixed = TRUE)
})

testthat::test_that("a singleton supermodule inherits its member module identity", {
  member_map <- data.frame(dataset = "d", SupermoduleID = "SM02",
                           ModuleID = "WGCNA_m03", stringsAsFactors = FALSE)
  module_evidence <- data.frame(
    dataset = "d", ModuleID = "WGCNA_m03", module_size = 200L,
    proposed_module_theme = "theme C", current_stage01_label = "L3",
    proposed_label = "Proposed C", review_action = "REFINE_WORDING",
    top10_hub_symbols = "C1; C2", stringsAsFactors = FALSE
  )
  structural <- data.frame(dataset = "d", SupermoduleID = "SM02",
                           adjusted_signed_min_pairwise_eigengene_correlation = NA_real_,
                           pc1_variance_explained = 1,
                           cut_height_stability_fraction_stable = 1,
                           stringsAsFactors = FALSE)
  sm <- wcl_supermodule_evidence(member_map, module_evidence, structural)
  sm <- wcl_apply_singleton_inheritance(sm, module_evidence)
  testthat::expect_true(sm$inherits_from_member_module)
  testthat::expect_identical(sm$structural_coherence_class, "singleton")
  testthat::expect_identical(sm$biological_coherence_class, "singleton")
  # inherits BOTH the label and the action, exactly
  testthat::expect_identical(sm$proposed_supermodule_label, "Proposed C")
  testthat::expect_identical(sm$supermodule_review_action, "REFINE_WORDING")
})

testthat::test_that("the label audit source contains no phenotype input", {
  # Only the SCRIPT is scanned: it is what actually resolves and reads inputs.
  # R/statistics/wgcna_label_coherence_utils.R is excluded because
  # wcl_forbidden_field_patterns() must literally contain these tokens in order
  # to blacklist them; that helper is covered by the column-level guard tests.
  code <- readLines(repo_path("analysis/05_wgcna", "14_wgcna_label_coherence_audit.R"
  ), warn = FALSE)
  live <- code[!grepl("^\\s*#", code)]
  text <- paste(live, collapse = "\n")
  # it must not read the DA manifest, the mapped contrast files, the Stage-07
  # group-effect handoff, or the phenotype-aware candidate/position outputs
  for (pattern in c("clusterProfiler_manifest", "sus_res_resolve_manifest_input",
                    "WGCNA_inferential_handoff", "wgcna_candidate_proteins",
                    "wgcna_sus_res_network_position", "log2fc", "aveExpr")) {
    testthat::expect_false(grepl(pattern, text, ignore.case = TRUE), info = pattern)
  }
})

# ------------------------------------------------- real-output contracts

wcl_review_path <- function(dataset, file) {
  path_results("reviewer_audit", "wgcna_label_review", dataset, file)
}

testthat::test_that("every module and supermodule appears exactly once", {
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    mp <- wcl_review_path(ds, "WGCNA_module_label_coherence_audit.csv")
    sp <- wcl_review_path(ds, "WGCNA_supermodule_label_coherence_audit.csv")
    if (!file.exists(mp) || !file.exists(sp)) next
    m <- utils::read.csv(mp, stringsAsFactors = FALSE)
    s <- utils::read.csv(sp, stringsAsFactors = FALSE)
    testthat::expect_identical(anyDuplicated(m$ModuleID), 0L, info = ds)
    testthat::expect_identical(anyDuplicated(s$SupermoduleID), 0L, info = ds)

    # matches the authoritative identity contract exactly
    contract_path <- path_results(
      "tables", "06_modules_WGCNA", "identity_contract", ds,
      "WGCNA_module_supermodule_membership_contract.csv")
    if (file.exists(contract_path)) {
      ct <- utils::read.csv(contract_path, stringsAsFactors = FALSE)
      testthat::expect_setequal(m$ModuleID, unique(ct$module_id))
      testthat::expect_setequal(s$SupermoduleID, unique(ct$supermodule_id))
      testthat::expect_identical(sum(s$n_member_modules), nrow(ct), info = ds)
    }
  }
})

testthat::test_that("real label-audit outputs stay phenotype-blind", {
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    for (f in c("WGCNA_module_label_coherence_audit.csv",
                "WGCNA_supermodule_label_coherence_audit.csv",
                "WGCNA_label_supporting_proteins.csv")) {
      p <- wcl_review_path(ds, f)
      if (!file.exists(p)) next
      x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
      testthat::expect_silent(wcl_assert_phenotype_blind(
        x, paste(ds, f), allow = "naming_evidence_is_phenotype_blind"))
    }
  }
})

testthat::test_that("canonical labels and module identities are untouched", {
  # the audit writes only under results/reviewer_audit/wgcna_label_review/
  for (ds in c("neuron_neuropil", "neuron_soma")) {
    proposed <- wcl_review_path(
      ds, paste0("PROPOSED_wgcna_reviewed_labels_", ds, ".csv"))
    if (!file.exists(proposed)) next
    x <- utils::read.csv(proposed, stringsAsFactors = FALSE)
    testthat::expect_true(all(x$review_status == "PROPOSED_NOT_PROMOTED"))
    # A PROPOSED registry must never be promoted wholesale. Part-29 activated
    # exactly one adjudicated module (neuropil m11), so if a config registry
    # exists it must be far smaller than the proposal it sits beside.
    cfg <- repo_path("config", "wgcna_labels", paste0(ds, ".csv"))
    if (file.exists(cfg)) {
      r <- utils::read.csv(cfg, stringsAsFactors = FALSE)
      testthat::expect_lt(nrow(r), nrow(x))
      testthat::expect_true(all(r$adjudication_status == "reviewed"), info = ds)
    } else {
      testthat::succeed()
    }
  }
  # the one reviewed registry that does exist is still microglia's
  testthat::expect_true(file.exists(repo_path("config", "wgcna_labels",
                                              "microglia.csv")))
})
