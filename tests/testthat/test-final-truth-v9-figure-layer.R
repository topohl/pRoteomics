# Part-25 contract tests for the final_truth_v9 layer.
#
# Section 32 lists the conditions that must fail the build. Each one is an
# assertion here, so the freeze is enforced by the suite rather than by review.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "null_coalescing.R"))
source(testthat::test_path("..", "..", "R", "integration_utils.R"))
source(testthat::test_path("..", "..", "R", "nature_v2_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "spatial_grammar_utils.R"))
source(testthat::test_path("..", "..", "R", "final_truth_v9_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_export.R"))
source(testthat::test_path("..", "..", "R", "final_truth_v9_fidelity_panels.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
V9 <- s9f_contract()
SD <- path_results("source_data", "manuscript_candidates", "final_truth_v9")
TAB <- path_results("tables", "manuscript_candidates", "final_truth_v9")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
FIGD <- path_results("figures", "manuscript_candidates", "final_truth_v9")
code_of <- function(path) {
  ln <- readLines(path, warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}
v9_render_src <- function()
  paste(vapply(Sys.glob(repo_path("R", "final_truth_v9_*panels.R")), code_of,
               character(1)), collapse = "\n")
sidecar <- function(fig, id) {
  p <- file.path(SD, fig, paste0(id, "_source_data.csv"))
  if (file.exists(p)) rd(p) else NULL
}

testthat::test_that("v9 is the promoted main-figure layer inside the size contract", {
  testthat::expect_identical(V9$contract_version, s9f_contract_version())
  # Phase 5B promoted this layer for Figures 2 and 3. Its Extended Data figures
  # are deliberately NOT promoted and stay candidates, which is why the status
  # names the two figures rather than the whole contract.
  testthat::expect_identical(V9$status, "promoted_canonical_manuscript_figures_2_and_3")
  for (f in V9$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
  }
})

# ---- S32: ED8b must not use a self-including CON consensus ----------------
testthat::test_that("ED8b compares no CON animal with a consensus containing itself", {
  p <- file.path(TAB, "ed8b_consensus_audit.csv")
  testthat::skip_if_not(file.exists(p), "consensus audit not built")
  a <- rd(p)
  testthat::expect_gt(nrow(a), 0L)
  # the fairness condition itself
  testthat::expect_false(any(a$self_included))
  # and the released numbers must match the fair construction, not the biased one
  testthat::expect_true(all(a$matches_leave_one_out))
  con <- a[a$group == "CON", , drop = FALSE]
  testthat::expect_gt(nrow(con), 0L)
  testthat::expect_false(any(con$matches_self_including))
  testthat::expect_true(all(con$fair_reference == "leave-one-CON-out centroid"))
})

# ---- S32: biological n must never be the acquisition count ---------------
testthat::test_that("acquisition count is never presented as biological n", {
  z <- sidecar("figure_02", "v9_depth")
  testthat::skip_if_not(!is.null(z), "figure_02 not built")
  testthat::expect_true(all(c("n_acquisitions_in_compartment",
                              "n_biological_replicates") %in% names(z)))
  testthat::expect_identical(unique(z$n_biological_replicates), 9L)
  testthat::expect_gt(sum(unique(z$n_acquisitions_in_compartment)), 9L)
  a <- rd(file.path(TAB, "final_panel_statistical_role_audit.csv"))
  d <- a[a$panel_id == "v9_depth", ]
  testthat::expect_true(grepl("acquisition", d$replicate_unit))
  testthat::expect_true(grepl("9 animals", d$biological_n))
})

# ---- S32: F2e source data must retain uncapped values ---------------------
testthat::test_that("display saturation never reaches the released source data", {
  z <- sidecar("figure_02", "v9_compartment")
  testthat::skip_if_not(!is.null(z), "figure_02 not built")
  cen <- z[z$colour_scale_censored %in% TRUE, , drop = FALSE]
  testthat::expect_gt(nrow(cen), 0L)
  testthat::expect_true(all(abs(cen$true_value) > abs(cen$displayed_value)))
  # no true value may have been silently replaced by the display cap
  testthat::expect_false(any(abs(cen$true_value - max(abs(z$displayed_value))) <
                               1e-9))
})

# ---- S32: phenotype-blind selection must not be called a hit -------------
testthat::test_that("the baseline fingerprint never says hits", {
  src <- code_of(repo_path("R", "final_truth_v9_fidelity_panels.R"))
  testthat::expect_false(grepl("top hits", src, fixed = TRUE))
  z <- sidecar("figure_02", "v9_fingerprint")
  testthat::skip_if_not(!is.null(z), "figure_02 not built")
  # and the ordering must stay phenotype-blind
  testthat::expect_true(all(c("gene_peak_unit", "gene_peak_dataset") %in%
                              names(z)))
  testthat::expect_false(any(grepl("SUS|RES", names(z))))
  # the peak is a GENE-level attribute over the whole fingerprint, so it pairs
  # with gene_peak_dataset. Naming it peak_unit invited joining it to the row
  # dataset, which produced 18 unresolvable soma/microglia laminar keys during
  # the Part-26 audit. The data were always right; the column name was not.
  bad <- z$gene_peak_dataset %in% c("neuron_soma", "microglia") &
    grepl("_(so|sr|slm|mo|po)$", z$gene_peak_unit)
  testthat::expect_false(any(bad))
  # one peak per gene, and every peak resolves against its own compartment
  g <- unique(z[, c("gene", "gene_peak_dataset", "gene_peak_unit")])
  testthat::expect_identical(nrow(g), length(unique(z$gene)))
  u <- sg_units()
  testthat::expect_identical(
    setdiff(paste(g$gene_peak_dataset, g$gene_peak_unit),
            paste(u$dataset, u$analysis_key)), character(0))
})

# ---- S32: panel h is characterisation, not independent validation --------
testthat::test_that("only external comparisons are called validation", {
  a <- rd(file.path(TAB, "final_panel_statistical_role_audit.csv"))
  indep <- a[grepl("^YES", a$independent_validation), , drop = FALSE]
  # exactly one panel may claim independent validation: the external one
  testthat::expect_identical(indep$panel_id, "v9_external_main")
  h <- a[a$panel_id == "v9_internal_main", ]
  testthat::expect_false(grepl("^YES", h$independent_validation))
  roles <- vapply(V9$panels, function(p) as.character(p$role %||% ""),
                  character(1))
  ids <- vapply(V9$panels, function(p) as.character(p$id), character(1))
  testthat::expect_identical(roles[ids == "v9_internal_main"],
                             "functional_characterization")
  testthat::expect_false(any(roles == "internal_validation"))
})

# ---- S32: theme aggregation must not be presented as FDR-tested ----------
testthat::test_that("the theme atlas is declared a descriptive aggregation", {
  a <- rd(file.path(TAB, "final_panel_statistical_role_audit.csv"))
  atl <- a[a$panel_id %in% c("v9_atlas", "v9_ed_atlas_rescon",
                             "v9_ed_atlas_suscon"), , drop = FALSE]
  testthat::expect_gt(nrow(atl), 0L)
  testthat::expect_true(all(grepl("descriptive", atl$descriptive_or_inferential)))
  # the theme cell has no FDR family of its own
  testthat::expect_true(all(grepl("no family of its own|constituent",
                                  atl$FDR_family)))
  leg <- paste(readLines(file.path(REP, "final_figure_legends_v9.md"),
                         warn = FALSE), collapse = " ")
  testthat::expect_true(grepl("does not constitute an additional multiple-testing family",
                              leg, fixed = TRUE))
})

# ---- S32: leading-edge proteins are not implied to be validated ----------
testthat::test_that("leading-edge panels are declared selected and descriptive", {
  a <- rd(file.path(TAB, "final_panel_statistical_role_audit.csv"))
  pr <- a[a$panel_id %in% c("v9_prot_syn", "v9_prot_rna", "v9_prot_ox"), ]
  testthat::expect_identical(nrow(pr), 3L)
  testthat::expect_true(all(grepl("same analysis", pr$independent_validation)))
  testthat::expect_true(all(grepl("leading-edge", pr$selection_dependency)))
  # and the claim that none is individually FDR-supported must remain true
  for (k in c("syn", "rna", "ox")) {
    z <- sidecar("figure_03", paste0("v9_prot_", k))
    testthat::skip_if_not(!is.null(z), "figure_03 not built")
    testthat::expect_false(any(z$BH_FDR < 0.05, na.rm = TRUE))
  }
})

# ---- S32: algebraically related contrasts must be declared as such -------
testthat::test_that("the three contrasts are stated to be algebraically related", {
  leg <- paste(readLines(file.path(REP, "final_figure_legends_v9.md"),
                         warn = FALSE), collapse = " ")
  testthat::expect_true(grepl("algebraically related", leg, fixed = TRUE))
  for (k in c("syn", "rna", "ox")) {
    z <- sidecar("figure_03", paste0("v9_prot_", k))
    testthat::skip_if_not(!is.null(z), "figure_03 not built")
    w <- stats::reshape(z[, c("gene", "contrast", "log2FC")], idvar = "gene",
                        timevar = "contrast", direction = "wide")
    names(w) <- sub("^log2FC[.]", "", names(w))
    cn <- setdiff(names(w), "gene")
    sr <- grep("SUS.RES", cn, value = TRUE)[1]
    sc <- grep("SUS.CON", cn, value = TRUE)[1]
    rc <- grep("RES.CON", cn, value = TRUE)[1]
    testthat::expect_lt(max(abs(w[[sr]] - (w[[sc]] - w[[rc]])), na.rm = TRUE),
                        1e-12)
  }
})

# ---- S32: publication-facing pipeline language --------------------------
testthat::test_that("pipeline language never reaches a reader", {
  # the stored classification KEY may keep its name; the displayed label may not
  testthat::expect_identical(unname(f9_qc_class_label("not_claimable_due_to_QC")),
                             "QC-sensitive")
  sup <- file.path(TAB, "supplementary")
  testthat::skip_if_not(dir.exists(sup), "supplementary tables not built")
  for (f in list.files(sup, "[.]csv$", full.names = TRUE)) {
    txt <- paste(readLines(f, warn = FALSE), collapse = " ")
    testthat::expect_false(grepl("claimable", txt, ignore.case = TRUE),
                           label = basename(f))
  }
  leg <- paste(readLines(file.path(REP, "final_figure_legends_v9.md"),
                         warn = FALSE), collapse = " ")
  testthat::expect_false(grepl("claimable", leg, ignore.case = TRUE))
  # and the banned movement vocabulary stays banned
  src <- v9_render_src()
  for (w in c("relocation", "redistribution", "migration", "hotspot")) {
    testthat::expect_false(grepl(w, src, ignore.case = TRUE), label = w)
  }
})

# ---- S32: WGCNA panel a is member abundance, not the eigengene -----------
testthat::test_that("module member abundance is not called an eigengene", {
  src <- code_of(repo_path("R", "final_truth_v9_ed_panels.R"))
  testthat::expect_true(grepl("Mean module-member", src, fixed = TRUE))
  testthat::expect_false(grepl('name = "Mean module abundance', src,
                               fixed = TRUE))
  z <- sidecar("extended_data", "v9_ed_module_fingerprint")
  testthat::skip_if_not(!is.null(z), "extended_data not built")
  testthat::expect_true("mean_con_z" %in% names(z))
})

# ---- S32: the ExpGroup stage must feed nothing --------------------------
testthat::test_that("the numeric ExpGroup stage feeds no v9 panel", {
  for (p in V9$panels) {
    src <- c(as.character(p$primary_source %||% ""),
             as.character(unlist(p$input_dependencies %||% list())))
    testthat::expect_false(any(grepl("variance_partitioning", src)),
                           label = as.character(p$id))
  }
  ki <- file.path(REP, "known_issues_v9.md")
  testthat::skip_if_not(file.exists(ki), "known issues not built")
  txt <- paste(readLines(ki, warn = FALSE), collapse = " ")
  testthat::expect_true(grepl("ExpGroup", txt, fixed = TRUE))
})

# ---- S32: stale burden claims and the CA2-SLM framing -------------------
testthat::test_that("the microglia burden and the CA2-SLM framing are unchanged", {
  z <- sidecar("figure_03", "v9_dap_track")
  testthat::skip_if_not(!is.null(z), "figure_03 not built")
  mg <- z[z$dataset == "microglia", ]
  testthat::expect_identical(mg$canonical[mg$display == "CA1"], 0L)
  testthat::expect_identical(mg$canonical[mg$display == "CA2"], 3L)
  testthat::expect_identical(mg$canonical[mg$display == "CA3"], 0L)
  testthat::expect_identical(z$canonical[z$unit == "CA2_slm"], 28L)
  testthat::expect_identical(z$claimable[z$unit == "CA2_slm"], 6L)
})

# ---- S32: vector export and the 5 pt floor ------------------------------
testthat::test_that("v9 pages are vector and clear the 5 pt floor", {
  pdfs <- character(0)
  for (f in V9$figures) {
    p <- path_results("figures", "manuscript_candidates", "final_truth_v9",
                      as.character(f$figure_key), "assembled",
                      paste0(as.character(f$name), ".pdf"))
    if (file.exists(p)) pdfs <- c(pdfs, p)
  }
  testthat::skip_if_not(length(pdfs) > 0, "no assembled PDFs")
  testthat::skip_if(nchar(Sys.which("qpdf")) == 0, "qpdf not available")
  aud <- e8_vector_audit(pdfs)
  testthat::expect_false(any(aud$page_sized_raster_present))
  testthat::expect_true(all(aud$embedded_fonts > 0))
  svgs <- list.files(
    path_results("figures", "manuscript_candidates", "final_truth_v9"),
    "[.]svg$", recursive = TRUE, full.names = TRUE)
  svgs <- grep("/assembled/", svgs, value = TRUE)
  sizes <- unlist(lapply(svgs, function(f) {
    s <- paste(readLines(f, warn = FALSE), collapse = " ")
    as.numeric(sub("px", "", sub("font-size: *", "",
      regmatches(s, gregexpr("font-size: *[0-9.]+px", s))[[1]])))
  }))
  testthat::expect_gt(length(sizes), 0L)
  testthat::expect_identical(sum(sizes < 4.995), 0L)
})

# ---- S30: every panel has a declared statistical role --------------------
testthat::test_that("every panel carries a statistical role and a legend", {
  a <- rd(file.path(TAB, "final_panel_statistical_role_audit.csv"))
  n_layout <- sum(vapply(V9$figures, function(f) length(f$layout), integer(1)))
  testthat::expect_identical(nrow(a), n_layout)
  testthat::expect_true(all(a$caption_sufficient))
  testthat::expect_true(all(nzchar(a$quantity)))
  testthat::expect_true(all(nzchar(a$replicate_unit)))
  testthat::expect_true(all(nzchar(a$FDR_family)))
})

# ---- S5/S6/S9/S29: the semantic term scan is a gate, not a report --------
#
# The scan reads every reader-facing artefact this layer emits and the text
# nodes of every panel SVG. A P0 is a statement the design cannot support -
# a pre-stress reference, a stress effect attributed to a within-stress
# contrast, a specificity or functional claim - and must never survive.
testthat::test_that("no P0 semantic hit survives in any reader-facing artefact", {
  p <- file.path(TAB, "audit", "semantic_search_hits.csv")
  testthat::skip_if_not(file.exists(p), "semantic scan not built")
  h <- rd(p)
  testthat::expect_true(all(c("file", "line_or_field", "term", "context",
                              "severity", "recommended_fix") %in% names(h)))
  testthat::expect_identical(sum(grepl("^P0", h$severity)), 0L)
})

testthat::test_that("the three language audits classify every hit they find", {
  for (nm in c("stress_language_audit.csv", "baseline_language_audit.csv",
               "phenotype_specificity_language_audit.csv")) {
    p <- file.path(TAB, "audit", nm)
    testthat::skip_if_not(file.exists(p), paste(nm, "not built"))
    a <- rd(p)
    testthat::expect_true(all(nzchar(a$severity)))
    testthat::expect_true(all(a$severity %in%
      c("OK", "P0 factual", "P1 overclaim", "P2 consistency")))
    # an unlicensed hit must carry an actionable fix; a licensed one must not
    # be presented as a defect
    testthat::expect_true(all(nzchar(a$recommended_fix)))
    testthat::expect_true(all(a$severity[a$licensed_by != "NONE"] == "OK"))
  }
})

testthat::test_that("baseline never appears unqualified in reader-facing prose", {
  p <- file.path(TAB, "audit", "baseline_language_audit.csv")
  testthat::skip_if_not(file.exists(p), "baseline audit not built")
  b <- rd(p)
  # every retained use of baseline names the control group in the same
  # sentence; nothing was measured before the paradigm
  testthat::expect_identical(sum(b$severity != "OK"), 0L)
  testthat::expect_true(all(grepl("terminal", unique(b$design_fact))))
})

testthat::test_that("no panel prints text that contradicts its own legend", {
  p <- file.path(TAB, "audit", "printed_panel_language_audit.csv")
  testthat::skip_if_not(file.exists(p), "printed panel audit not built")
  testthat::expect_identical(nrow(rd(p)), 0L)
})

# ---- S29b: a short write must not pass as a complete figure --------------
testthat::test_that("every emitted SVG is structurally complete", {
  svgs <- list.files(
    path_results("figures", "manuscript_candidates", "final_truth_v9"),
    "[.]svg$", recursive = TRUE, full.names = TRUE)
  testthat::skip_if_not(length(svgs) > 0, "no SVGs built")
  ok <- vapply(svgs, function(f) {
    n <- file.info(f)$size
    con <- file(f, "rb"); on.exit(close(con))
    seek(con, max(0, n - 64))
    isTRUE(n > 2000) &&
      grepl("</svg>", rawToChar(readBin(con, "raw", 64L)), fixed = TRUE)
  }, logical(1))
  testthat::expect_identical(sum(!ok), 0L)
})

# ---- S13: only externally anchored panels may be called validation -------
testthat::test_that("the story text calls exactly one F2 panel validation", {
  p <- file.path(REP, "final_figure_story_v9.md")
  testthat::skip_if_not(file.exists(p), "story not built")
  s <- readLines(p, warn = FALSE)
  testthat::expect_identical(sum(grepl("g and h validate", s)), 0L)
  testthat::expect_true(any(grepl("is not independent validation", s)))
})

# ---- S10/S14-S19: the declared rules are enforced, not only declared --------
testthat::test_that("no claim-strength rule is broken in any artefact", {
  p <- file.path(TAB, "audit", "claim_strength_language_audit.csv")
  testthat::skip_if_not(file.exists(p), "claim-strength audit not built")
  a <- rd(p)
  testthat::expect_identical(sum(a$severity != "OK"), 0L)
})

testthat::test_that("pathway is never claimed for GO evidence", {
  p <- file.path(TAB, "audit", "program_vs_pathway_audit.csv")
  testthat::skip_if_not(file.exists(p), "pathway audit not built")
  a <- rd(p)
  present <- a[a$present_in_this_project, , drop = FALSE]
  testthat::expect_gt(nrow(present), 0L)
  # every collection actually used is GO, so none licenses the word pathway
  testthat::expect_identical(sum(present$pathway_label_justified), 0L)
  testthat::expect_true(all(grepl("GO biological process", present$definition)))
})

# ---- S27: unsupportable story clauses are removed, not softened -------------
testthat::test_that("the core story drops every clause it cannot support", {
  p <- file.path(TAB, "audit", "core_story_audit.csv")
  testthat::skip_if_not(file.exists(p), "core story audit not built")
  a <- rd(p)
  testthat::expect_identical(sort(unique(a$sentence)), 1:4)
  testthat::expect_true(all(nzchar(a$supporting_artefact)))
  # richer, without and global are the three that no artefact supports
  testthat::expect_identical(sort(a$clause[!a$supported]),
                             c("global", "richer", "without"))
  # the removed clauses are still named in the rationale, so check the
  # preferred version itself rather than the whole document
  s <- readLines(file.path(REP, "core_story_corrected.md"), warn = FALSE)
  h <- grep("^## ", s)
  start <- grep("^## Preferred version", s)
  testthat::expect_identical(length(start), 1L)
  stop_at <- h[h > start]
  pref <- paste(s[seq(start + 1L,
                      if (length(stop_at)) stop_at[1] - 1L else length(s))],
                collapse = " ")
  testthat::expect_false(grepl("\\bricher\\b", pref))
  testthat::expect_false(grepl("\\bwithout\\b", pref))
  testthat::expect_true(grepl("did not detect", pref))
})

# ---- Part-28 S21: no heatmap may clip a value without disclosure -----------
testthat::test_that("no heatmap clips silently", {
  p <- file.path(TAB, "audit", "heatmap_scale_integrity_audit.csv")
  testthat::skip_if_not(file.exists(p), "heatmap audit not built")
  h <- rd(p)
  testthat::expect_gt(nrow(h), 0L)
  testthat::expect_identical(sum(h$status == "SILENT_CLIP"), 0L)
  # a panel that saturates must disclose it on the colourbar AND keep the
  # uncapped values in its source data
  sat <- h[h$n_above_colour_max + h$n_below_colour_min > 0L, , drop = FALSE]
  if (nrow(sat)) {
    testthat::expect_true(all(sat$intentional_saturation))
    testthat::expect_true(all(sat$disclosed_on_colourbar))
    testthat::expect_true(all(sat$uncapped_source_values_present))
  }
  # out-of-range values would be painted grey50 by the ggplot2 default oob
  testthat::expect_identical(sum(h$grey50_pixels_in_svg), 0L)
})

testthat::test_that("only the saturating tail of F2e is marked", {
  p <- file.path(FIGD, "figure_02", "panels", "v9_compartment.svg")
  testthat::skip_if_not(file.exists(p), "F2e not rendered")
  txt <- paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  labs <- sub(".*>([^<]*)</text>", "\\1",
              unlist(regmatches(txt, gregexpr("<text[^>]*>[^<]*</text>", txt))))
  # the upper tail saturates, the lower one does not
  testthat::expect_true(any(grepl("≥3", labs)))
  testthat::expect_false(any(grepl("≤", labs)))
})

# ---- Part-28 S5/S30: every atlas row label covers its own major blocks -----
testthat::test_that("atlas row labels match the semantic decision table", {
  p <- file.path(TAB, "audit", "atlas_annotation_decisions.csv")
  testthat::skip_if_not(file.exists(p), "annotation decisions not built")
  d <- rd(p)
  testthat::expect_identical(nrow(d), 7L)
  testthat::expect_true(all(nzchar(d$recommended_label)))
  # the recommendation is what the renderer actually prints
  src <- readLines(repo_path("R", "final_truth_v9_panels.R"), warn = FALSE)
  i <- grep("^  SHORT <- c\\(", src)
  testthat::expect_identical(length(i), 1L)
  blk <- paste(src[i:(i + 7)], collapse = " ")
  for (lab in d$recommended_label)
    testthat::expect_true(grepl(lab, blk, fixed = TRUE))
})

testthat::test_that("the mitochondrial theme contains no cytosolic glycolysis", {
  # Under registry v3 the glycolysis sub-DAG is excluded by one ontology rule,
  # so the theme is mitochondrial bioenergetics and the row name matches its
  # membership. Part 28 had to call the row "Energy metabolism" precisely
  # because v2 still reached cytosolic glycolysis; that is no longer true.
  TH <- rd(repo_path("results", "tables", "10_biological_integration",
                     "gsea_wgcna_concordance", "global",
                     "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
  testthat::expect_identical(unique(TH$registry_version),
                             "manuscript_go_themes_v3")
  mito <- unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE &
                            TH$theme_id == "mitochondrial_respiration_oxphos"])
  GLY <- c("GO:0006096", "GO:0061621", "GO:0061615", "GO:0061620")
  testthat::expect_identical(sum(mito %in% GLY), 0L)
  # the rule is ontological, so the mitochondrial entry steps survive it
  testthat::expect_true("GO:0006086" %in% mito)   # pyruvate -> acetyl-CoA
  testthat::expect_true("GO:0006099" %in% mito)   # TCA cycle
  # and the row is named for what it now contains
  src <- readLines(repo_path("R", "final_truth_v9_panels.R"), warn = FALSE)
  testthat::expect_gt(
    sum(grepl('mitochondrial_respiration_oxphos = "Mitochondrial respiration"',
              src, fixed = TRUE)), 0L)
})

# ---- Part-28 S12: the atlas summary is robust to GO redundancy -------------
testthat::test_that("redundancy reduction does not overturn the atlas", {
  p <- file.path(TAB, "audit", "atlas_redundancy_reduced_comparison.csv")
  testthat::skip_if_not(file.exists(p), "redundancy comparison not built")
  k <- rd(p)
  testthat::expect_gt(nrow(k), 0L)
  # representatives are medoids, so every cell must have at least one
  testthat::expect_identical(sum(k$n_semantic_representatives < 1L), 0L)
  # the typical cell moves far less than one colour step
  testthat::expect_lt(stats::median(k$absolute_difference),
                      k$large_shift_threshold[1])
  # and the great majority keep their direction
  testthat::expect_gt(mean(!k$sign_changed), 0.95)
})

# ---- Part-28 S17/S18: a module is never named by external evidence alone ---
testthat::test_that("WGCNA annotation confidence rests on enrichment", {
  p <- file.path(TAB, "audit", "wgcna_annotation_evidence_matrix.csv")
  testthat::skip_if_not(file.exists(p), "WGCNA evidence matrix not built")
  w <- rd(p)
  testthat::expect_identical(nrow(w), 35L)
  # a module with no enrichment may never be called functional
  noe <- w[!w$has_enrichment, , drop = FALSE]
  if (nrow(noe))
    testthat::expect_true(all(noe$annotation_confidence %in%
                                c("CELL_CONTEXT_ONLY", "UNRESOLVED")))
  # the permitted wording always keeps the module ID visible
  testthat::expect_true(all(mapply(grepl, w$module_id,
                                   w$allowed_manuscript_form, fixed = TRUE)))
  # and never permits the bare "the X module" form
  testthat::expect_false(any(grepl("^the .* module$", w$allowed_manuscript_form)))
})

testthat::test_that("a label contradicted by its own enrichment is withheld", {
  p <- file.path(TAB, "audit", "wgcna_module_naming_audit.csv")
  e <- file.path(TAB, "audit", "wgcna_annotation_evidence_matrix.csv")
  testthat::skip_if_not(file.exists(p) && file.exists(e), "naming audit not built")
  n <- rd(p); w <- rd(e)
  flagged <- w[w$label_review_required %in% TRUE, , drop = FALSE]
  testthat::expect_gt(nrow(flagged), 0L)
  k <- n[match(paste(flagged$dataset, flagged$module_id),
               paste(n$dataset, n$module_id)), , drop = FALSE]
  testthat::expect_true(all(grepl("refer to by module ID",
                                  k$allowed_manuscript_label)))
  # neuropil m11 specifically: neither the contradicted active label nor the
  # unactivated oligodendrocyte proposal may be used
  m11 <- n[n$dataset == "neuron_neuropil" & n$module_id == "m11", ]
  testthat::expect_identical(nrow(m11), 1L)
  testthat::expect_false(grepl("^synaptic", m11$allowed_manuscript_label))
  testthat::expect_true(grepl("not activated", m11$prohibited_label))
})

testthat::test_that("no WGCNA module label is printed on any panel", {
  svgs <- list.files(FIGD, pattern = "[.]svg$", recursive = TRUE,
                     full.names = TRUE)
  svgs <- grep("[/\\]panels[/\\]", svgs, value = TRUE)
  testthat::skip_if_not(length(svgs) > 0, "no panels rendered")
  txt <- paste(unlist(lapply(svgs, function(f)
    readLines(f, warn = FALSE, encoding = "UTF-8"))), collapse = " ")
  for (lab in c("proteostasis", "cytoskeletal trafficking",
                "RNA/RNP regulatory"))
    testthat::expect_false(grepl(lab, txt, fixed = TRUE))
})
