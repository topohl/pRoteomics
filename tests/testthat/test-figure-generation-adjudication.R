# Guards for the Phase 5 figure-generation adjudication.
#
# The manuscript's Results 2 and 3 are written against the final_truth_v9 figure
# generation, which is not yet promoted. These tests pin the scientific facts the
# promotion decision rests on, so that none of them can drift between the
# adjudication and the promotion itself. They assert nothing about which
# generation is canonical - that decision is open and is recorded in
# manuscript/figure_promotion_blockers.csv.

repo_rel <- function(...) file.path(testthat::test_path("..", ".."), ...)

testthat::test_that("the curated atlas is registry v3 with exactly the seven manuscript programs", {
  p <- repo_rel("results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
                "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
  testthat::skip_if_not(file.exists(p), "atlas theme table not present")
  d <- utils::read.csv(p, stringsAsFactors = FALSE)

  SEVEN <- c("rna_processing_splicing_rnp", "ribosome_translation",
             "chromatin_organization", "mitochondrial_respiration_oxphos",
             "synaptic_signaling_vesicle", "neuron_projection_development",
             "autophagy_lysosome_endosome")
  testthat::expect_true(all(SEVEN %in% d$theme_id))
  testthat::expect_identical(unique(d$registry_version[d$theme_id %in% SEVEN]),
                             "manuscript_go_themes_v3")

  # 253 constituent GO-BP terms is the number the manuscript states.
  per <- vapply(SEVEN, function(t) length(unique(d$GO_ID[d$theme_id == t])), integer(1))
  testthat::expect_identical(sum(per), 253L)
})

testthat::test_that("the mitochondrial theme excludes glycolysis and keeps PDH and TCA", {
  p <- repo_rel("results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
                "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
  testthat::skip_if_not(file.exists(p), "atlas theme table not present")
  d <- utils::read.csv(p, stringsAsFactors = FALSE)
  m <- unique(d[d$theme_id == "mitochondrial_respiration_oxphos", c("GO_ID", "GO_description")])

  # The v3 exclusion rule, documented in docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md.
  # A stale 20-term version that included glycolytic terms exists in the project's
  # history and in one report; it must never come back through the data.
  testthat::expect_identical(nrow(m), 16L)
  testthat::expect_false("GO:0006096" %in% m$GO_ID)
  testthat::expect_false(any(grepl("glycol", m$GO_description, ignore.case = TRUE)))
  testthat::expect_true("GO:0006086" %in% m$GO_ID)   # pyruvate decarboxylation
  testthat::expect_true("GO:0006099" %in% m$GO_ID)   # tricarboxylic acid cycle
  testthat::expect_true("GO:0006119" %in% m$GO_ID)   # the OXPHOS exemplar term
})

testthat::test_that("the DAP source carries the audited QC status and the 37/28/6/15 arithmetic", {
  p <- repo_rel("results", "tables", "11_spatial_systems", "atlas",
                "protein_spatial_cell_affinity.csv")
  testthat::skip_if_not(file.exists(p), "spatial cell affinity atlas not present")
  d <- utils::read.csv(p, stringsAsFactors = FALSE)

  # The QC columns are what distinguish the post-audit source from the pre-audit
  # one. Their absence is what makes the v2 DAP panel unusable.
  for (col in c("CA2_SLM_robustness_class", "QC_claimability", "robustness_audit_status")) {
    testthat::expect_true(col %in% names(d), info = paste("missing QC column:", col))
  }

  f <- d[d$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
  testthat::expect_identical(nrow(f), 37L)

  ca2 <- f$sus_res_strongest_spatial_unit_canonical == "ca2_slm"
  testthat::expect_identical(sum(ca2), 28L)

  cls <- table(f$CA2_SLM_robustness_class[ca2])
  testthat::expect_identical(as.integer(cls[["robust_to_missingness_and_QC"]]), 6L)
  testthat::expect_identical(as.integer(cls[["not_claimable_due_to_QC"]]), 10L)
  testthat::expect_identical(as.integer(cls[["insufficient_observed_data"]]), 12L)

  # 9 outside CA2-SLM were never at risk and carry no audited class.
  testthat::expect_identical(sum(!ca2), 9L)
  testthat::expect_true(all(!nzchar(f$CA2_SLM_robustness_class[!ca2])))

  # 6 audited-robust + 9 never-at-risk = the 15 the manuscript reports.
  testthat::expect_identical(6L + sum(!ca2), 15L)
})

testthat::test_that("the three manuscript exemplars are exact and resolve in canonical data", {
  p <- repo_rel("results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
                "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
  testthat::skip_if_not(file.exists(p), "atlas theme table not present")
  th <- utils::read.csv(p, stringsAsFactors = FALSE)

  ex <- data.frame(
    dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
    unit    = c("CA3_sr", "CA2_sp", "CA1"),
    go      = c("GO:0099536", "GO:0006397", "GO:0006119"),
    stringsAsFactors = FALSE)

  for (i in seq_len(nrow(ex))) {
    z <- th[th$dataset == ex$dataset[i] & th$spatial_unit == ex$unit[i] &
              th$GO_ID == ex$go[i], , drop = FALSE]
    z <- z[!duplicated(z$contrast), , drop = FALSE]
    testthat::expect_true(all(c("RES - CON", "SUS - CON", "SUS - RES") %in% z$contrast),
                          info = paste("exemplar missing a contrast:", ex$go[i]))
  }

  # The exemplar registry itself must not drift. It is prespecified, not chosen
  # from the data, and the manuscript commits to exactly these three.
  src <- repo_rel("R", "spatial_v6_figure3_panels.R")
  testthat::skip_if_not(file.exists(src), "exemplar definition not present")
  code <- paste(readLines(src, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl('go_id = c("GO:0099536", "GO:0006397", "GO:0006119")',
                              code, fixed = TRUE))
})

testthat::test_that("the adjudication record is present and internally consistent", {
  inv <- repo_rel("manuscript", "figure_generation_inventory.csv")
  mat <- repo_rel("manuscript", "figure_claim_support_matrix.csv")
  blk <- repo_rel("manuscript", "figure_promotion_blockers.csv")
  for (f in c(inv, mat, blk)) testthat::expect_true(file.exists(f))

  m <- utils::read.csv(mat, stringsAsFactors = FALSE)
  ok <- c("SUPPORTED_BY_V2", "SUPPORTED_BY_V9", "SUPPORTED_BY_BOTH", "SUPPORTED_BY_NEITHER")
  testthat::expect_true(all(m$classification %in% ok))
  # No manuscript claim may be left with no display in either generation.
  testthat::expect_identical(sum(m$classification == "SUPPORTED_BY_NEITHER"), 0L)

  b <- utils::read.csv(blk, stringsAsFactors = FALSE)
  testthat::expect_true(all(b$blocks_promotion %in% c("YES", "NO")))
  # Phase 5B closed PB-01, PB-02 and PB-03 and the promotion happened, so the
  # register must now be clear. This was the tripwire that stopped the contract
  # being switched while blockers were still open; inverted, it now stops a new
  # blocker being recorded while the promoted contract stays in place.
  testthat::expect_identical(sum(b$blocks_promotion == "YES"), 0L)
  testthat::expect_true(all(b$id[b$severity == "RESOLVED"] %in% c("PB-01", "PB-02", "PB-03")))
})

testthat::test_that("the column-registration defect is fixed and stays fixed", {
  # PB-01, now resolved. NF_RGT reserves the atlas-legend gutter; before the
  # repair it was applied in the DAP track only, so panels a and b of Figure 3
  # did not share column geometry and the headline 28 pointed at CA3 stratum
  # oriens. The tripwire this test used to carry has fired and been inverted:
  # both sides of the gutter convention must now be present.
  src <- repo_rel("R", "final_truth_v9_panels.R")
  testthat::skip_if_not(file.exists(src), "v9 panel renderers not present")
  code <- readLines(src, warn = FALSE)

  dap <- grep("plot.margin = ggplot2::margin(1, NF_RGT, 0, 1, \"mm\")", code, fixed = TRUE)
  testthat::expect_gte(length(dap), 2L)   # the DAP track and the coupled atlas

  atlas <- paste(code[seq(grep("^f9_gsea_atlas <- function", code)[1],
                          length(code))], collapse = "\n")
  atlas <- sub("\nf9_[a-z_]+ <- function.*$", "", atlas)
  testthat::expect_true(grepl("NF_RGT", atlas, fixed = TRUE),
    info = "f9_gsea_atlas must reserve the same right gutter as the track above it")
  testthat::expect_true(grepl("shares_column_geometry_with", atlas, fixed = TRUE))

  # The measured proof lives in the registration audit.
  aud <- repo_rel("manuscript", "figure_pb01_registration_audit.csv")
  testthat::skip_if_not(file.exists(aud), "registration audit not present")
  a <- utils::read.csv(aud, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(a), 18L)
  testthat::expect_identical(sum(a$status == "MISREGISTERED"), 0L)
  testthat::expect_identical(sum(a$status_before == "MISREGISTERED"), 18L)
  # The column that carries the headline value must be CA2 SLM.
  hit <- a[a$expected_value_canonical == 28L, ]
  testthat::expect_identical(nrow(hit), 1L)
  testthat::expect_identical(as.character(hit$expected_spatial_unit), "CA2 SLM")
  testthat::expect_identical(as.character(hit$plotted_spatial_unit_after), "CA2 SLM")
  testthat::expect_identical(as.character(hit$plotted_spatial_unit_before), "CA3 SO")
})
