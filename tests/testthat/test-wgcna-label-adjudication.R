source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "wgcna_label_adjudication_utils.R"))

# ---------------------------------------------------------------- fixtures

wla_fx_members <- function(n = 20L, module = "WGCNA_m01", dataset = "d") {
  out <- data.frame(
    dataset = dataset, ModuleID = module,
    ProteinGroupID = sprintf("PG:%s:%02d", dataset, seq_len(n)),
    GeneSymbol = paste0("G", seq_len(n)),
    EntrezID = as.character(seq_len(n) * 11L),
    RepresentativeUniProt = paste0("U", seq_len(n)),
    abs_kME = seq(0.98, 0.20, length.out = n),
    abs_kME_rank_in_module = seq_len(n), n_module_members = n,
    stringsAsFactors = FALSE
  )
  out$rank_fraction <- (out$abs_kME_rank_in_module - 0.5) / n
  out[["is_core_kME_0.6"]] <- out$abs_kME >= 0.60
  out$is_top_hub_25 <- out$abs_kME_rank_in_module <= 25L
  out$is_top10_module_hub <- out$abs_kME_rank_in_module <= 10L
  out$is_top5_module_representative <- out$abs_kME_rank_in_module <= 5L
  out
}

# geneIDs are Entrez ids: ranks 1..5 -> 11/22/33/44/55 (central)
wla_fx_go <- function(module = "WGCNA_m01", set = "all", central = TRUE,
                      padj = 1e-6, n_terms = 3L) {
  ids <- if (central) "11/22/33/44/55" else "176/187/198/209/220"
  data.frame(
    ModuleID = module, Ontology = "BP", ModuleProteinSetType = set,
    ID = sprintf("GO:%04d", seq_len(n_terms)),
    Description = c("synapse organization", "synaptic vesicle cycle",
                    "trans-synaptic signaling", "postsynaptic density organization",
                    "axon guidance")[seq_len(n_terms)],
    p.adjust = padj, qvalue = padj,
    GeneRatio = "5/20", Count = 5L, geneID = ids,
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# PHENOTYPE BLINDNESS
# =====================================================================

testthat::test_that("no phenotype or group-effect field can enter a naming output", {
  ok <- data.frame(ModuleID = "WGCNA_m01", dominant_theme = "x",
                   contributor_centrality_auc = 0.7, stringsAsFactors = FALSE)
  testthat::expect_silent(wla_assert_phenotype_blind(ok))
  for (bad in c("sus_res_min_BH_FDR", "log2FC", "candidate_tier", "is_tier_A1",
                "group_effect_direction", "module_contrast", "n_DAP",
                "effect_size", "bh_fdr")) {
    d <- ok; d[[bad]] <- 1
    testthat::expect_error(wla_assert_phenotype_blind(d), "phenotype-derived",
                           info = bad)
  }
  att <- ok; att$naming_evidence_is_phenotype_blind <- TRUE
  testthat::expect_silent(
    wla_assert_phenotype_blind(att, allow = "naming_evidence_is_phenotype_blind"))
})

testthat::test_that("the adjudication script never reads phenotype inputs", {
  code <- readLines(repo_path("analysis/05_wgcna", "15_wgcna_label_adjudication.R"), warn = FALSE)
  live <- paste(code[!grepl("^\\s*#", code)], collapse = "\n")
  for (p in c("clusterProfiler_manifest", "sus_res_resolve_manifest_input",
              "WGCNA_inferential_handoff", "wgcna_candidate_proteins",
              "wgcna_sus_res_network_position", "log2fc", "group_effects")) {
    testthat::expect_false(grepl(p, live, ignore.case = TRUE), info = p)
  }
})

# =====================================================================
# THEME LAYER
# =====================================================================

testthat::test_that("theme assignment is deterministic and loses nothing", {
  d <- c("oxidative phosphorylation", "mitochondrial translation",
         "synapse organization", "keratinocyte differentiation",
         "a totally novel process")
  a <- wla_assign_theme(d)
  b <- wla_assign_theme(d)
  testthat::expect_identical(a, b)
  testthat::expect_identical(a[[1]], "mitochondrial energy metabolism")
  testthat::expect_identical(a[[3]], "synaptic organization / signalling")
  testthat::expect_identical(a[[4]], "epithelial / keratinization")
  # an unmatched term survives as its own theme rather than being dropped
  testthat::expect_identical(a[[5]], "other: a totally novel process")
  testthat::expect_false(anyNA(a))
})

testthat::test_that("themes pool contributors across their constituent GO terms", {
  members <- wla_fx_members()
  go <- wla_fx_go(n_terms = 2L)
  go$geneID <- c("11/22", "33/44")          # disjoint contributors
  th <- wla_module_themes(go, members, "WGCNA_m01", "all", "BP")
  testthat::expect_identical(nrow(th), 1L)
  testthat::expect_identical(th$n_supporting_go_terms, 2L)
  # union of both terms, not just the best one
  testthat::expect_identical(th$n_contributing_proteins, 4L)
  testthat::expect_match(th$supporting_go_ids, "GO:0001; GO:0002", fixed = TRUE)
})

testthat::test_that("only FDR-significant terms enter a theme", {
  members <- wla_fx_members()
  go <- rbind(wla_fx_go(padj = 1e-6, n_terms = 1L),
              wla_fx_go(padj = 0.40, n_terms = 1L))
  go$Description[2] <- "keratinocyte differentiation"
  go$ID[2] <- "GO:9999"
  th <- wla_module_themes(go, members, "WGCNA_m01", "all", "BP")
  testthat::expect_identical(nrow(th), 1L)
  testthat::expect_false("epithelial / keratinization" %in% th$theme)
})

testthat::test_that("protein sets stay distinguishable and are never merged", {
  members <- wla_fx_members()
  go <- rbind(wla_fx_go(set = "all"), wla_fx_go(set = "core_kME_0.6"),
              wla_fx_go(set = "top_hub_25"))
  for (ps in wla_protein_sets()) {
    th <- wla_module_themes(go, members, "WGCNA_m01", ps, "BP")
    testthat::expect_identical(unique(th$protein_set), ps)
  }
})

# =====================================================================
# CENTRALITY SUPPORT
# =====================================================================

testthat::test_that("centrality support uses contributor identity, not hub text", {
  members <- wla_fx_members()
  central <- wla_module_themes(wla_fx_go(central = TRUE), members, "WGCNA_m01", "all")
  periph <- wla_module_themes(wla_fx_go(central = FALSE), members, "WGCNA_m01", "all")

  testthat::expect_equal(central$fraction_top10_contributing, 0.5)
  testthat::expect_equal(periph$fraction_top10_contributing, 0)
  testthat::expect_gt(central$contributor_centrality_auc, 0.9)
  testthat::expect_lt(periph$contributor_centrality_auc, 0.1)

  testthat::expect_true(wla_centre_supports_theme(central))
  testthat::expect_false(wla_centre_supports_theme(periph))

  # components stay separately visible; nothing is collapsed into one score
  s <- wla_centrality_support(central)
  testthat::expect_true(all(c("top10_theme_fraction", "top25_theme_fraction",
                              "core_theme_fraction", "contributor_median_abs_kME",
                              "contributor_centrality_auc") %in% names(s)))
  testthat::expect_identical(nrow(s), 1L)
})

testthat::test_that("no arbitrary presence bonus survives in the helper source", {
  code <- readLines(repo_path("R", "wgcna_label_adjudication_utils.R"), warn = FALSE)
  live <- paste(code[!grepl("^\\s*#", code)], collapse = "\n")
  # the specific defect: an additive presence bonus. A 0.75 used as a
  # comparison threshold elsewhere is legitimate and must not trip this.
  testthat::expect_false(grepl("[+] *0[.]75", live))
  testthat::expect_false(grepl("nzchar\\(hub_text\\)", live))
  testthat::expect_false(grepl("hub_support", live))
  # support must be length-safe even with no theme at all
  testthat::expect_identical(nrow(wla_centrality_support(NULL)), 1L)
})

# =====================================================================
# ADJUDICATION LOGIC
# =====================================================================

adjudicate <- function(all = NULL, core = NULL, hub = NULL, current = NA_character_) {
  wla_adjudicate_module("WGCNA_m01", "d", 20L, all, core, hub,
                        current_label = current)
}

testthat::test_that("a nonsignificant-only module cannot be high confidence", {
  a <- adjudicate(all = NULL, core = NULL, hub = NULL)
  testthat::expect_identical(a$coherence_class, "unresolved")
  testthat::expect_identical(a$proposed_confidence, "low")
  testthat::expect_identical(a$recommended_action, "UNRESOLVED")
  testthat::expect_identical(a$proposed_primary_label, "mixed / unresolved")
})

testthat::test_that("high confidence requires convergence AND centre support", {
  members <- wla_fx_members()
  th_all <- wla_module_themes(wla_fx_go(set = "all", n_terms = 3L), members, "WGCNA_m01", "all")
  th_core <- wla_module_themes(wla_fx_go(set = "core_kME_0.6", n_terms = 3L), members, "WGCNA_m01", "core_kME_0.6")
  a <- adjudicate(th_all, th_core, NULL)
  testthat::expect_identical(a$proposed_confidence, "high")
  testthat::expect_true(a$centre_supports_theme)

  # the same theme in only ONE protein set, peripherally driven -> low
  th_periph <- wla_module_themes(wla_fx_go(central = FALSE, n_terms = 3L), members, "WGCNA_m01", "all")
  b <- adjudicate(th_periph, NULL, NULL)
  testthat::expect_identical(b$proposed_confidence, "low")
  testthat::expect_false(b$centre_supports_theme)
  testthat::expect_identical(b$coherence_class, "peripheral_enrichment_only")
  testthat::expect_identical(b$recommended_action, "RENAME")
})

testthat::test_that("thin single-set evidence is never moderate confidence", {
  members <- wla_fx_members()
  # present only in the core set, with no all-set evidence
  th_core <- wla_module_themes(wla_fx_go(set = "core_kME_0.6", n_terms = 1L),
                               members, "WGCNA_m01", "core_kME_0.6")
  a <- adjudicate(NULL, th_core, NULL)
  testthat::expect_true(a$dominant_theme_evidence_is_thin)
  testthat::expect_identical(a$proposed_confidence, "low")
})

testthat::test_that("top-hub evidence alone cannot create a specific label", {
  members <- wla_fx_members()
  th_hub <- wla_module_themes(wla_fx_go(set = "top_hub_25", n_terms = 3L),
                              members, "WGCNA_m01", "top_hub_25")
  a <- adjudicate(NULL, NULL, th_hub)
  # no all-set and no core-set support -> thin, and never high confidence
  testthat::expect_true(a$dominant_theme_evidence_is_thin)
  testthat::expect_false(a$proposed_confidence == "high")
  testthat::expect_identical(a$n_supporting_go_terms_all, 0L)
})

testthat::test_that("raw ontology grammar is reworded without changing the biology", {
  testthat::expect_identical(
    wla_readable_label("mitochondrial energy metabolism",
                       "generation of precursor metabolites and energy"),
    "mitochondrial energy metabolism")
  testthat::expect_identical(
    wla_readable_label("myelin / oligodendrocyte ensheathment",
                       "ensheathment of neurons"),
    "oligodendrocyte / myelin ensheathment")
  # a context-sensitive term becomes an explicit signature, not a cell-type claim
  testthat::expect_match(
    wla_readable_label("epithelial / keratinization", "keratinocyte differentiation"),
    "signature")
  # an unmapped theme is returned unchanged
  testthat::expect_identical(wla_readable_label("RNA processing / RNP", "rna splicing"),
                             "RNA processing / RNP")
})

testthat::test_that("overly broad terms are flagged and never silently kept", {
  broad <- vapply(c("protein phosphorylation", "regulation of transport",
                    "hydrolase activity"),
                  function(d) any(vapply(wla_overly_broad_patterns(),
                                         function(p) grepl(p, d), logical(1))),
                  logical(1))
  testthat::expect_true(all(broad))
  specific <- any(vapply(wla_overly_broad_patterns(),
                         function(p) grepl(p, "synaptic vesicle exocytosis"),
                         logical(1)))
  testthat::expect_false(specific)
})

# =====================================================================
# SUPERMODULES
# =====================================================================

wla_fx_module_adj <- function() {
  data.frame(
    dataset = "d", ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"),
    module_size = c(1400L, 40L, 200L),
    dominant_theme = c("theme A", "theme B", "theme C"),
    proposed_primary_label = c("P1", "P2", "P3"),
    proposed_confidence = c("high", "high", "moderate"),
    recommended_action = c("KEEP", "KEEP", "REFINE_WORDING"),
    top10_hub_symbols = c("A1; A2; A3; A4", "B1; B2; B3", "C1; C2; C3"),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("supermodule themes weight member modules equally", {
  member_map <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM01", "SM02"),
    ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"), stringsAsFactors = FALSE)
  structural <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM02"),
    adjusted_signed_min_pairwise_eigengene_correlation = c(0.7, NA),
    pc1_variance_explained = c(0.8, 1),
    cut_height_stability_fraction_stable = c(0.9, 1), stringsAsFactors = FALSE)

  sm <- wla_adjudicate_supermodules(member_map, wla_fx_module_adj(), structural,
                                    hubs_per_module = 2L)
  testthat::expect_identical(nrow(sm), 2L)
  s1 <- sm[sm$SupermoduleID == "SM01", ]
  # the 1400-protein module does not outvote the 40-protein one
  testthat::expect_identical(s1$n_modules_supporting_dominant_theme, 1L)
  testthat::expect_equal(s1$dominant_theme_module_fraction, 0.5)
  # a mixed block is not given the large module's theme as its name
  testthat::expect_identical(s1$biological_coherence_class,
                             "mixed_but_structurally_coherent")
  testthat::expect_identical(s1$proposed_supermodule_label, "mixed / unresolved")
  testthat::expect_identical(s1$recommended_action, "MIXED")
  # the hub panel draws from EVERY member module
  testthat::expect_match(s1$balanced_hub_panel, "WGCNA_m01: A1, A2", fixed = TRUE)
  testthat::expect_match(s1$balanced_hub_panel, "WGCNA_m02: B1, B2", fixed = TRUE)
})

testthat::test_that("a coherent supermodule gets one name, a singleton inherits", {
  member_map <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM01", "SM02"),
    ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"), stringsAsFactors = FALSE)
  adj <- wla_fx_module_adj(); adj$dominant_theme[2] <- "theme A"   # now unanimous
  structural <- data.frame(
    dataset = "d", SupermoduleID = c("SM01", "SM02"),
    adjusted_signed_min_pairwise_eigengene_correlation = c(0.7, NA),
    pc1_variance_explained = c(0.8, 1),
    cut_height_stability_fraction_stable = c(0.9, 1), stringsAsFactors = FALSE)
  sm <- wla_adjudicate_supermodules(member_map, adj, structural)

  s1 <- sm[sm$SupermoduleID == "SM01", ]
  testthat::expect_identical(s1$biological_coherence_class, "coherent_single_program")
  testthat::expect_identical(s1$proposed_supermodule_label, "theme A")
  testthat::expect_identical(s1$proposed_confidence, "high")

  s2 <- sm[sm$SupermoduleID == "SM02", ]
  testthat::expect_true(s2$inherits_from_member_module)
  testthat::expect_identical(s2$structural_coherence_class, "singleton")
  # inherits the member module's label, confidence AND action exactly
  testthat::expect_identical(s2$proposed_supermodule_label, "P3")
  testthat::expect_identical(s2$proposed_confidence, "moderate")
  testthat::expect_identical(s2$recommended_action, "REFINE_WORDING")
})

# =====================================================================
# PROPOSAL REGISTRIES
# =====================================================================

wla_fx_proposal <- function() {
  wla_proposal_registry(
    "d",
    data.frame(ModuleID = c("WGCNA_m01", "WGCNA_m02"),
               proposed_primary_label = c("P1", "P2"),
               dominant_theme = c("t1", "t2"),
               proposed_confidence = c("high", "low"),
               recommended_action = c("KEEP", "MIXED"),
               coherence_class = c("high_confidence_coherent", "mixed_biology"),
               centre_supports_theme = c(TRUE, FALSE),
               n_supporting_go_terms_all = c(5L, 1L),
               n_protein_sets_supporting_theme = c(3L, 1L),
               stringsAsFactors = FALSE),
    data.frame(SupermoduleID = "SM01",
               proposed_supermodule_label = "S1", dominant_theme = "t1",
               proposed_confidence = "high", recommended_action = "KEEP",
               structural_coherence_class = "singleton",
               biological_coherence_class = "singleton",
               n_modules_supporting_dominant_theme = 1L, n_member_modules = 1L,
               stringsAsFactors = FALSE))
}

testthat::test_that("a generated proposal is never an active registry", {
  reg <- wla_fx_proposal()
  testthat::expect_true(all(reg$adjudication_status == "proposed"))
  testthat::expect_true(all(is.na(reg$reviewer)))
  testthat::expect_identical(unique(reg$proposal_prepared_by),
                             "automated_evidence_audit")
  testthat::expect_false(wla_is_active_registry(reg))

  # only a named human reviewer plus an adjudicated status makes it active
  active <- reg
  active$adjudication_status <- "adjudicated"
  active$reviewer <- "A Human"
  testthat::expect_true(wla_is_active_registry(active))
  # either one alone is not enough
  half <- reg; half$adjudication_status <- "adjudicated"
  testthat::expect_false(wla_is_active_registry(half))
})

testthat::test_that("proposal validation derives identity from authoritative state", {
  reg <- wla_fx_proposal()
  testthat::expect_silent(wla_validate_proposal(
    reg, "d", c("WGCNA_m01", "WGCNA_m02"), "SM01"))

  # a stale ModuleID must fail
  testthat::expect_error(
    wla_validate_proposal(reg, "d", c("WGCNA_m01", "WGCNA_m99"), "SM01"),
    "do not match the authoritative module set")
  # a stale SupermoduleID must fail
  testthat::expect_error(
    wla_validate_proposal(reg, "d", c("WGCNA_m01", "WGCNA_m02"), "SM99"),
    "do not match the authoritative")
  # a fabricated reviewer must fail
  faked <- reg; faked$reviewer <- "Someone"
  testthat::expect_error(wla_validate_proposal(
    faked, "d", c("WGCNA_m01", "WGCNA_m02"), "SM01"), "must not name a reviewer")
  # a promoted status must fail the proposal contract
  promoted <- reg; promoted$adjudication_status <- "adjudicated"
  testthat::expect_error(wla_validate_proposal(
    promoted, "d", c("WGCNA_m01", "WGCNA_m02"), "SM01"), "adjudication_status")
  # duplicates must fail
  dup <- rbind(reg, reg[1, ])
  testthat::expect_error(wla_validate_proposal(
    dup, "d", c("WGCNA_m01", "WGCNA_m02"), "SM01"), "duplicate")
})

# =====================================================================
# REAL OUTPUT CONTRACTS
# =====================================================================

wla_out <- function(f) path_results("reviewer_audit", "wgcna_label_adjudication", f)

testthat::test_that("every current module and supermodule has exactly one packet", {
  mp <- wla_out("WGCNA_module_adjudication.csv")
  sp <- wla_out("WGCNA_supermodule_adjudication.csv")
  testthat::skip_if_not(file.exists(mp) && file.exists(sp),
                        "adjudication has not been generated")
  m <- utils::read.csv(mp, stringsAsFactors = FALSE)
  s <- utils::read.csv(sp, stringsAsFactors = FALSE)

  for (ds in unique(m$dataset)) {
    ct_path <- path_results("tables", "06_modules_WGCNA", "identity_contract", ds,
                            "WGCNA_module_supermodule_membership_contract.csv")
    if (!file.exists(ct_path)) next
    ct <- utils::read.csv(ct_path, stringsAsFactors = FALSE)
    md <- m[m$dataset == ds, ]; sd <- s[s$dataset == ds, ]
    testthat::expect_identical(anyDuplicated(md$ModuleID), 0L, info = ds)
    testthat::expect_identical(anyDuplicated(sd$SupermoduleID), 0L, info = ds)
    # IDs are unchanged and complete
    testthat::expect_setequal(md$ModuleID, unique(ct$module_id))
    testthat::expect_setequal(sd$SupermoduleID, unique(ct$supermodule_id))
    testthat::expect_identical(sum(sd$n_member_modules), nrow(ct), info = ds)
  }
})

testthat::test_that("real outputs are phenotype-blind and conservatively graded", {
  for (f in c("WGCNA_module_adjudication.csv", "WGCNA_supermodule_adjudication.csv",
              "WGCNA_module_theme_evidence.csv", "WGCNA_module_top25_hubs.csv")) {
    p <- wla_out(f)
    if (!file.exists(p)) next
    x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    testthat::expect_silent(wla_assert_phenotype_blind(
      x, f, allow = "naming_evidence_is_phenotype_blind"))
  }
  p <- wla_out("WGCNA_module_adjudication.csv")
  testthat::skip_if_not(file.exists(p))
  m <- utils::read.csv(p, stringsAsFactors = FALSE)
  testthat::expect_true(all(m$coherence_class %in% wla_coherence_classes()))
  testthat::expect_true(all(m$recommended_action %in% wla_actions()))
  testthat::expect_true(all(m$proposed_confidence %in% wla_confidence_levels()))
  # a high-confidence proposal must have centre support and multi-set convergence
  high <- m[m$proposed_confidence == "high", , drop = FALSE]
  if (nrow(high)) {
    testthat::expect_true(all(as.logical(high$centre_supports_theme)))
    testthat::expect_true(all(high$n_protein_sets_supporting_theme >= 2L))
    testthat::expect_true(all(high$n_supporting_go_terms_all >= 3L))
  }
  # nothing unresolved or peripheral may be high confidence
  bad <- m[m$coherence_class %in% c("unresolved", "peripheral_enrichment_only",
                                    "mixed_biology"), , drop = FALSE]
  if (nrow(bad)) testthat::expect_false(any(bad$proposed_confidence == "high"))
})

testthat::test_that("hub tables follow the frozen abs_kME ordering", {
  p <- wla_out("WGCNA_module_top25_hubs.csv")
  testthat::skip_if_not(file.exists(p))
  h <- utils::read.csv(p, stringsAsFactors = FALSE)
  for (k in unique(paste(h$dataset, h$ModuleID))) {
    sub <- h[paste(h$dataset, h$ModuleID) == k, , drop = FALSE]
    sub <- sub[order(sub$rank), , drop = FALSE]
    testthat::expect_true(all(diff(sub$abs_kME) <= 1e-12), info = k)
    testthat::expect_identical(sub$rank, seq_len(nrow(sub)), info = k)
    testthat::expect_true(all(sub$is_top5[sub$rank <= 5]), info = k)
    testthat::expect_true(all(sub$is_top10[sub$rank <= 10]), info = k)
    testthat::expect_false(any(as.logical(sub$is_top10[sub$rank > 10])), info = k)
  }
})

testthat::test_that("proposals on disk are not active and config is untouched", {
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    p <- wla_out(paste0("proposed_", ds, "_reviewed_labels.csv"))
    if (!file.exists(p)) next
    r <- utils::read.csv(p, stringsAsFactors = FALSE)
    testthat::expect_true(all(r$adjudication_status == "proposed"), info = ds)
    testthat::expect_true(all(is.na(r$reviewer) | r$reviewer == ""), info = ds)
    testthat::expect_false(wla_is_active_registry(r), info = ds)
  }
  # the only registry under config/ is still microglia's
  testthat::expect_true(file.exists(repo_path("config", "wgcna_labels", "microglia.csv")))
  # Part-29 finalization: neuron_neuropil now has a registry, but it must
  # contain ONLY the single adjudicated module and must never assert a
  # cell-type identity. neuron_soma still has none. This keeps the original
  # guard - no bulk auto-activation - while allowing the reviewed entry.
  testthat::expect_false(
    file.exists(repo_path("config", "wgcna_labels", "neuron_soma.csv")))
  np <- repo_path("config", "wgcna_labels", "neuron_neuropil.csv")
  if (file.exists(np)) {
    r <- utils::read.csv(np, stringsAsFactors = FALSE)
    testthat::expect_identical(nrow(r), 1L)
    testthat::expect_identical(r$entity_id, "WGCNA_m11")
    testthat::expect_identical(r$adjudication_status, "reviewed")
    testthat::expect_true(grepl("myelin-associated",
                                r$reviewed_biological_label, fixed = TRUE))
    # a co-abundance module may be named for its protein composition, never
    # for a cell type it cannot establish
    testthat::expect_false(grepl("oligodendrocyte",
                                 r$reviewed_biological_label, ignore.case = TRUE))
    testthat::expect_false(grepl("oligodendrocyte",
                                 r$reviewed_short_label, ignore.case = TRUE))
  }
})

testthat::test_that("nonsignificant GO terms cannot yield a high-confidence label", {
  p <- wla_out("WGCNA_module_adjudication.csv")
  testthat::skip_if_not(file.exists(p))
  m <- utils::read.csv(p, stringsAsFactors = FALSE)
  ns <- m[!is.finite(m$best_p_adjust_all) | m$best_p_adjust_all > 0.05, , drop = FALSE]
  if (nrow(ns)) testthat::expect_false(any(ns$proposed_confidence == "high"))
})

# =====================================================================
# REVIEWED-REGISTRY GENERALIZATION (must not weaken validation)
# =====================================================================

testthat::test_that("expected counts derive from the authoritative member map", {
  source(repo_path("R", "wgcna_reviewed_label_registry.R"))
  src <- paste(readLines(repo_path("R", "wgcna_reviewed_label_registry.R"), warn = FALSE), collapse = "\n")

  # the microglia-only literals are gone from the signature and the checks
  testthat::expect_false(grepl("expected_n_modules = 13L", src, fixed = TRUE))
  testthat::expect_false(grepl("expected_n_supermodules = 9L", src, fixed = TRUE))
  testthat::expect_false(grepl("exactly six singleton supermodules", src, fixed = TRUE))
  testthat::expect_false(grepl("nrow(lookup) != 22L", src, fixed = TRUE))
  testthat::expect_false(grepl('"microglia_wgcna_reviewed_labels_v1"', src, fixed = TRUE))
  # counts are now derived
  testthat::expect_true(grepl("expected_n_modules <- length(current_modules)", src, fixed = TRUE))
  testthat::expect_true(grepl("expected_n_supermodules <- length(current_supers)", src, fixed = TRUE))
  # identity is still required, not merely cardinality
  testthat::expect_true(grepl("stale or missing ModuleID/SupermoduleID", src, fixed = TRUE))
})

testthat::test_that("the canonical lookup check rejects stale IDs for any dataset", {
  source(repo_path("R", "wgcna_reviewed_label_registry.R"))
  cols <- wgcna_canonical_lookup_columns()
  mk_lookup <- function(module_ids, super_ids) {
    n <- length(module_ids) + length(super_ids)
    out <- as.data.frame(matrix(NA_character_, nrow = n, ncol = length(cols)),
                         stringsAsFactors = FALSE)
    names(out) <- cols
    out$dataset <- "neuron_soma"
    out$level <- c(rep("module", length(module_ids)), rep("supermodule", length(super_ids)))
    out$entity_id <- c(module_ids, super_ids)
    out$canonical_biological_label <- "label"
    out$canonical_short_label <- "short"
    out$canonical_plot_label <- paste0(out$entity_id, "\nshort")
    out$final_plot_label <- out$canonical_plot_label
    out$structural_coherence_class <- "coherent_stable"
    if ("aggregation_evidence_class" %in% cols) out$aggregation_evidence_class <- "descriptive_only"
    out
  }
  member_map <- data.frame(
    dataset = "neuron_soma",
    ModuleID = c("WGCNA_m01", "WGCNA_m02"),
    SupermoduleID = c("SM01", "SM01"), stringsAsFactors = FALSE)

  # Only the new cardinality/identity branch is exercised here; a fully valid
  # lookup would also have to satisfy the later membership-fingerprint contract,
  # which is covered by the real microglia registry test below.

  # a stale ModuleID must fail even though the COUNT is still right
  stale_mod <- mk_lookup(c("WGCNA_m01", "WGCNA_m99"), "SM01")
  testthat::expect_error(
    wgcna_validate_canonical_lookup(stale_mod, "neuron_soma", member_map),
    "stale or missing")
  # a stale SupermoduleID must fail
  stale_sm <- mk_lookup(c("WGCNA_m01", "WGCNA_m02"), "SM99")
  testthat::expect_error(
    wgcna_validate_canonical_lookup(stale_sm, "neuron_soma", member_map),
    "stale or missing")
  # a missing module must fail on cardinality
  short <- mk_lookup("WGCNA_m01", "SM01")
  testthat::expect_error(
    wgcna_validate_canonical_lookup(short, "neuron_soma", member_map),
    "must contain exactly")
})

testthat::test_that("the microglia registry still validates after generalization", {
  source(repo_path("R", "wgcna_reviewed_label_registry.R"))
  cfg <- repo_path("config", "wgcna_labels", "microglia.csv")
  contract <- path_results("tables", "06_modules_WGCNA", "identity_contract",
                           "microglia",
                           "WGCNA_module_supermodule_membership_contract.csv")
  testthat::skip_if_not(file.exists(cfg) && file.exists(contract))
  reg <- readr::read_csv(cfg, show_col_types = FALSE, progress = FALSE)
  ct <- utils::read.csv(contract, stringsAsFactors = FALSE)
  member_map <- data.frame(dataset = "microglia", ModuleID = ct$module_id,
                           SupermoduleID = ct$supermodule_id,
                           stringsAsFactors = FALSE)
  # counts are now derived, so no explicit 13/9 is supplied by the caller
  testthat::expect_silent(
    wgcna_validate_reviewed_registry(reg, "microglia", member_map))
})

testthat::test_that("the automatic fallback still applies when no registry is active", {
  source(repo_path("R", "wgcna_label_adjudication_utils.R"))
  # the proposals generated by this pass must never be treated as active,
  # so the automatic/canonical labels remain in force
  for (ds in c("neuron_neuropil", "neuron_soma")) {
    p <- wla_out(paste0("proposed_", ds, "_reviewed_labels.csv"))
    if (!file.exists(p)) next
    r <- utils::read.csv(p, stringsAsFactors = FALSE)
    testthat::expect_false(wla_is_active_registry(r), info = ds)
  }
  # and the Stage-07 automatic label lookup is still present for those datasets
  for (ds in c("neuron_neuropil", "neuron_soma")) {
    lk <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds,
                       "WGCNA_final_label_lookup.csv")
    if (file.exists(lk)) {
      x <- utils::read.csv(lk, stringsAsFactors = FALSE)
      testthat::expect_gt(nrow(x), 0L)
    }
  }
})

testthat::test_that("the consumer migration plan is complete and does not migrate Figure 3", {
  plan <- wla_consumer_migration_plan()
  testthat::expect_true(all(c("script", "output", "current_label_field",
                              "label_stage", "category",
                              "recommended_canonical_field", "migration_priority",
                              "risk", "action_after_adjudication") %in% names(plan)))
  testthat::expect_true(all(plan$category %in%
    c("identity_critical", "biological_display", "historical_provenance")))
  testthat::expect_true(all(plan$migration_priority %in%
    c("high", "medium", "low", "none")))
  testthat::expect_identical(anyDuplicated(plan$script), 0L)

  fig3 <- plan[grepl("figure3", plan$script, ignore.case = TRUE), , drop = FALSE]
  testthat::expect_identical(nrow(fig3), 1L)
  testthat::expect_match(fig3$action_after_adjudication, "Do NOT migrate", fixed = TRUE)
  testthat::expect_identical(fig3$migration_priority, "low")

  # audit tables deliberately retain the Stage-01 label as provenance
  prov <- plan[plan$category == "historical_provenance", , drop = FALSE]
  testthat::expect_gt(nrow(prov), 0L)
  testthat::expect_true(all(prov$migration_priority == "none"))
})
