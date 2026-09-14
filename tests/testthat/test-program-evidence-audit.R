# Regression tests for the atlas program-evidence audit.
#
# These assert the contract a reviewer relies on: that every displayed row is
# traceable from registry anchor to GO term to supported term to leading-edge
# protein, and that the protein evidence is counted only over FDR-supported
# results.

source(testthat::test_path("..", "..", "R", "paths.R"))

pe <- function(f) path_results("tables", "publication_audits",
                               "program_evidence", f)
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE,
                                  check.names = FALSE)
PROGS <- c("rna_processing_splicing_rnp", "ribosome_translation",
           "chromatin_organization", "mitochondrial_respiration_oxphos",
           "synaptic_signaling_vesicle", "neuron_projection_development",
           "autophagy_lysosome_endosome")

testthat::test_that("the ledger covers exactly the seven displayed programs", {
  p <- pe("program_evidence_ledger.csv")
  testthat::skip_if_not(file.exists(p), "ledger not built")
  l <- rd(p)
  testthat::expect_setequal(unique(l$theme_id), PROGS)
  testthat::expect_true(all(nzchar(l$program_label)))
  # every row must say why the term is in the program
  testthat::expect_true(all(nzchar(l$why_in_program)))
  testthat::expect_true(all(l$membership_class %in%
                              c("DEFINITION_ANCHOR", "CLOSURE_DESCENDANT")))
})

testthat::test_that("every ledger row traces back to the canonical table", {
  p <- pe("program_evidence_ledger.csv")
  testthat::skip_if_not(file.exists(p), "ledger not built")
  l <- rd(p)
  TH <- rd(repo_path("results", "tables", "10_biological_integration",
                     "gsea_wgcna_concordance", "global",
                     "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
  testthat::expect_identical(unique(TH$registry_version),
                             "manuscript_go_themes_v3")
  k <- paste(l$theme_id, l$GO_ID, l$source_comparison)
  ref <- paste(TH$theme_id, TH$GO_ID, TH$source_comparison)
  testthat::expect_identical(sum(!k %in% ref), 0L)
  # support status must agree with the canonical FDR, not be re-derived
  m <- match(k, ref)
  testthat::expect_identical(
    sum(l$supported != (is.finite(TH$GSEA_FDR[m]) & TH$GSEA_FDR[m] < 0.05)), 0L)
})

testthat::test_that("protein evidence is counted only over supported results", {
  p <- pe("program_leading_edge_protein_evidence.csv")
  l <- pe("program_evidence_ledger.csv")
  testthat::skip_if_not(file.exists(p) && file.exists(l), "not built")
  z <- rd(p); led <- rd(l)
  testthat::expect_setequal(unique(z$theme_id), PROGS)
  # a gene can never be credited with more supported terms than the program has
  testthat::expect_true(all(z$n_supported_terms <= z$n_program_supported_terms))
  for (prog in PROGS) {
    n_sup_terms <- length(unique(led$GO_ID[led$theme_id == prog & led$supported]))
    testthat::expect_identical(unique(z$n_program_supported_terms[z$theme_id == prog]),
                               as.integer(n_sup_terms))
  }
  # the core rule must be exactly as documented
  core <- z$evidence_class == "RECURRENT_CORE"
  testthat::expect_true(all(z$n_supported_terms[core] >= 3L &
                              z$n_supported_contexts[core] >= 3L))
  single <- z$evidence_class == "SINGLE_APPEARANCE"
  testthat::expect_true(all(z$n_supported_terms[single] == 1L &
                              z$n_supported_contexts[single] == 1L))
})

testthat::test_that("the cross-program summary is internally consistent", {
  p <- pe("program_cross_summary.csv")
  testthat::skip_if_not(file.exists(p), "summary not built")
  s <- rd(p)
  testthat::expect_identical(nrow(s), 7L)
  testthat::expect_identical(s$theme_id, PROGS)
  testthat::expect_true(all(s$n_terms_supported <= s$n_terms))
  testthat::expect_true(all(s$verdict %in%
    c("STRONG", "SUPPORTED_BUT_BROAD", "MIXED", "MISLEADING")))
  # no row may be published as MISLEADING without action
  testthat::expect_identical(sum(s$verdict == "MISLEADING"), 0L)
  testthat::expect_true(all(nzchar(s$recommended_label)))
  testthat::expect_true(all(nzchar(s$top_core_proteins)))
})

testthat::test_that("the mitochondrial row contains no cytosolic glycolysis", {
  p <- pe("program_evidence_ledger.csv")
  testthat::skip_if_not(file.exists(p), "ledger not built")
  l <- rd(p)
  mito <- unique(l$GO_ID[l$theme_id == "mitochondrial_respiration_oxphos"])
  testthat::expect_identical(
    sum(mito %in% c("GO:0006096", "GO:0061621", "GO:0061615", "GO:0061620")), 0L)
  # and its recurrent core must be mitochondrial, not glycolytic
  z <- rd(pe("program_leading_edge_protein_evidence.csv"))
  core <- z$gene[z$theme_id == "mitochondrial_respiration_oxphos" &
                   z$evidence_class == "RECURRENT_CORE"]
  testthat::expect_gt(sum(grepl("^(Nduf|mt-Nd|Uqcr|Cox|Atp5|Sdh)", core)), 10L)
})

testthat::test_that("the omitted re-audit still finds no genuine omission", {
  p <- pe("omitted_program_recheck.csv")
  testthat::skip_if_not(file.exists(p), "re-audit not built")
  o <- rd(p)
  testthat::expect_gt(nrow(o), 0L)
  testthat::expect_identical(sum(o$is_true_omission), 0L)
  testthat::expect_true(all(nzchar(o$rejection_reason)))
})

testthat::test_that("each program has a review sheet", {
  d <- path_results("reports", "publication_audits", "program_evidence",
                    "program_sheets")
  testthat::skip_if_not(dir.exists(d), "sheets not built")
  f <- list.files(d, pattern = "[.]md$")
  testthat::expect_identical(length(f), 7L)
  for (x in f) {
    txt <- paste(readLines(file.path(d, x), warn = FALSE), collapse = " ")
    for (sec in c("Why this row is in the atlas", "Definition rule",
                  "Constituent GO terms", "Supported GO terms",
                  "Recurrent leading-edge proteins", "Take-home",
                  "Would I defend this label"))
      testthat::expect_true(grepl(sec, txt, fixed = TRUE), info = paste(x, sec))
  }
})
