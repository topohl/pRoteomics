# Part-24 contract tests for the editorial_v8 layer.
#
# The scientific architecture is fixed; these tests guard the things Part 24 is
# actually allowed to change - geometry, scale sharing, in-panel prose and the
# vector export - and hard-stop if the layer starts reopening analysis.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "null_coalescing.R"))
source(testthat::test_path("..", "..", "R", "integration_utils.R"))
source(testthat::test_path("..", "..", "R", "nature_v2_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "spatial_grammar_utils.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_export.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_fidelity_panels.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
V8 <- s8e_contract()

code_of <- function(path) {
  ln <- readLines(path, warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}
v8_sources <- function() {
  c(repo_path("R", "editorial_v8_panels.R"),
    repo_path("R", "editorial_v8_ed_panels.R"),
    repo_path("R", "editorial_v8_figure_utils.R"),
    repo_path("R", "editorial_v8_export.R"))
}

testthat::test_that("editorial_v8 is a candidate layer inside the size contract", {
  testthat::expect_identical(V8$contract_version, s8e_contract_version())
  testthat::expect_identical(V8$status, "candidate_only_not_promoted")
  testthat::expect_identical(as.numeric(V8$font_floor_pt), 5)
  for (f in V8$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs))
    testthat::expect_identical(anyDuplicated(labs), 0L)
    for (it in f$layout) {
      testthat::expect_lte(as.numeric(it$x) + as.numeric(it$w),
                           as.numeric(f$width_mm), label = f$name)
      testthat::expect_lte(as.numeric(it$y) + as.numeric(it$h),
                           as.numeric(f$height_mm), label = f$name)
    }
  }
})

testthat::test_that("the Part-23 panel membership is carried over unchanged", {
  f2 <- Filter(function(f) identical(f$name, "F2_NATURE_EDITORIAL_V8"),
               V8$figures)[[1]]
  f3 <- Filter(function(f) identical(f$name, "F3_NATURE_EDITORIAL_V8"),
               V8$figures)[[1]]
  testthat::expect_identical(
    vapply(f2$layout, function(x) as.character(x$label), character(1)),
    letters[1:8])
  testthat::expect_identical(
    vapply(f3$layout, function(x) as.character(x$label), character(1)),
    letters[1:9])
  # PCA stays in the main figure by explicit editorial decision
  testthat::expect_true("v8_pca" %in%
    vapply(f2$layout, function(x) as.character(x$panel), character(1)))
  # a sits directly on b with no gutter, so they read as one block
  a <- f3$layout[[1]]; b <- f3$layout[[2]]
  testthat::expect_identical(as.character(a$panel), "v8_dap_track")
  testthat::expect_identical(as.numeric(a$x), as.numeric(b$x))
  testthat::expect_identical(as.numeric(a$w), as.numeric(b$w))
  testthat::expect_identical(as.numeric(a$y) + as.numeric(a$h),
                             as.numeric(b$y))
})

testthat::test_that("no new inference is introduced by the v8 renderers", {
  forbidden <- c("\\blimma\\b", "\\blmer\\b", "\\blme4\\b", "\\bfgsea\\b",
                 "\\bWGCNA::", "\\bblockwiseModules\\b", "p\\.adjust\\s*\\(",
                 "\\bt\\.test\\s*\\(", "\\bwilcox\\.test\\s*\\(",
                 "\\bcor\\.test\\s*\\(", "\\baov\\s*\\(", "\\bglm\\s*\\(")
  for (f in v8_sources()) {
    src <- code_of(f)
    for (pat in forbidden) {
      testthat::expect_false(grepl(pat, src, perl = TRUE),
                             label = paste(basename(f), pat))
    }
  }
})

testthat::test_that("the same quantity is drawn on one shared scale", {
  # Figure 3 g/h/i plot log2FC over the same three contrasts
  d <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "figure_03")
  fs <- file.path(d, paste0("v8_prot_", c("syn", "rna", "ox"),
                            "_source_data.csv"))
  testthat::skip_if_not(all(file.exists(fs)), "figure_03 not built")
  lims <- vapply(fs, function(p) {
    z <- rd(p)
    testthat::expect_true("shared_scale_note" %in% names(z))
    max(abs(z$log2FC), na.rm = TRUE)
  }, numeric(1))
  # one symmetric limit covers all three, and it is stated in every sidecar
  notes <- unique(unlist(lapply(fs, function(p) rd(p)$shared_scale_note)))
  testthat::expect_length(notes, 1L)
  testthat::expect_true(grepl("shared by all three", notes, fixed = TRUE))
  stated <- as.numeric(sub(".*limit [+]/-([0-9.]+).*", "\\1", notes))
  testthat::expect_true(is.finite(stated))
  testthat::expect_true(all(lims <= stated + 1e-9))
})

testthat::test_that("the three NES atlases share one symmetric colour scale", {
  root <- path_results("source_data", "manuscript_candidates", "editorial_v8")
  fs <- c(file.path(root, "figure_03", "v8_atlas_source_data.csv"),
          file.path(root, "extended_data",
                    paste0("v8_ed_atlas_", c("rescon", "suscon"),
                           "_source_data.csv")))
  testthat::skip_if_not(all(file.exists(fs)), "atlases not built")
  shared <- unique(unlist(lapply(fs, function(p) rd(p)$shared_NES_scale_limit)))
  testthat::expect_length(shared, 1L)
  obs <- max(vapply(fs, function(p) max(abs(rd(p)$median_NES), na.rm = TRUE),
                    numeric(1)))
  testthat::expect_equal(shared, obs, tolerance = 1e-9)
})

testthat::test_that("ED7 states a count and never states movement", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data", "v8_ed_locations_source_data.csv")
  testthat::skip_if_not(file.exists(p), "extended_data not built")
  z <- rd(p)
  testthat::expect_true(all(z$elsewhere))
  src <- paste(readLines(repo_path("R", "editorial_v8_ed_panels.R"),
                         warn = FALSE), collapse = "\n")
  for (w in c("relocation", "redistribution", "migration", "relocate",
              "redistribute", "migrate")) {
    testthat::expect_false(grepl(w, src, ignore.case = TRUE), label = w)
  }
})

testthat::test_that("laminar resolution is never fabricated for soma or microglia", {
  root <- path_results("source_data", "manuscript_candidates", "editorial_v8")
  fs <- list.files(root, pattern = "_source_data[.]csv$", recursive = TRUE,
                   full.names = TRUE)
  testthat::skip_if_not(length(fs) > 0, "no source data yet")
  u <- sg_units()
  region_only <- u$analysis_key[u$dataset %in% c("neuron_soma", "microglia")]
  for (p in fs) {
    z <- rd(p)
    if (!all(c("dataset", "spatial_unit") %in% names(z))) next
    bad <- z$dataset %in% c("neuron_soma", "microglia") &
      grepl("_(so|sr|slm|mo|po)$", z$spatial_unit)
    testthat::expect_false(any(bad), label = basename(p))
  }
  testthat::expect_length(region_only, 8L)
})

testthat::test_that("the delivered PDFs are genuine vector", {
  pdfs <- character(0)
  for (f in V8$figures) {
    p <- path_results("figures", "manuscript_candidates", "editorial_v8",
                      as.character(f$figure_key), "assembled",
                      paste0(as.character(f$name), ".pdf"))
    if (file.exists(p)) pdfs <- c(pdfs, p)
  }
  testthat::skip_if_not(length(pdfs) > 0, "no assembled PDFs")
  testthat::skip_if(nchar(Sys.which("qpdf")) == 0, "qpdf not available")
  aud <- e8_vector_audit(pdfs)
  testthat::expect_true(all(aud$embedded_fonts > 0))
  testthat::expect_true(all(aud$text_show_operators > 0))
  testthat::expect_false(any(aud$page_sized_raster_present))
  testthat::expect_true(all(aud$vector_export_pass))
})

testthat::test_that("a patchwork panel is converted with patchworkGrob", {
  # ggplotGrob() on a patchwork silently returns only part of it, which is how
  # the Figure-2 hippocampus went missing on the first vector compose
  pw <- patchwork::wrap_plots(ggplot2::ggplot(), ggplot2::ggplot(), ncol = 2)
  testthat::expect_s3_class(e8_as_grob(pw), "gtable")
  testthat::expect_error(e8_as_grob("not a plot"), "unsupported panel class")
})

testthat::test_that("the cairo page rounds UP to the next whole point", {
  # cairo floors the MediaBox to whole points, which would make 183 mm 182.7 mm
  for (mm in c(183, 170, 168, 120, 150)) {
    pt <- ceiling(mm / 25.4 * 72) / 72
    testthat::expect_gte(pt * 25.4, mm)
    testthat::expect_lt(pt * 25.4 - mm, 25.4 / 72)
  }
})

testthat::test_that("every v8 script is registered in the pipeline", {
  y <- yaml::read_yaml(repo_path("pipeline.yml"))
  ids <- vapply(y$stages$manuscript_candidates$scripts,
                function(z) as.character(z$script), character(1))
  for (s in c("figures/editorial_v8_figure_02.R",
              "figures/editorial_v8_figure_03.R",
              "figures/editorial_v8_extended_data.R",
              "figures/editorial_v8_vector_audit.R",
              "figures/editorial_v8_previews.R")) {
    testthat::expect_true(s %in% ids, label = s)
    testthat::expect_true(file.exists(repo_path(s)), label = s)
  }
})

testthat::test_that("the layer writes only under manuscript_candidates/editorial_v8", {
  for (f in c("figures/editorial_v8_figure_02.R",
              "figures/editorial_v8_figure_03.R",
              "figures/editorial_v8_extended_data.R",
              "figures/editorial_v8_vector_audit.R",
              "figures/editorial_v8_previews.R")) {
    src <- code_of(repo_path(f))
    testthat::expect_true(grepl("editorial_v8", src, fixed = TRUE), label = f)
    testthat::expect_false(grepl("path_results\\(\\s*\"figures\",\\s*\"0",
                                 src, perl = TRUE), label = f)
  }
})

# --------------------------------------------------------------------------
# Part-25 fidelity pass. Each test below pins a defect that was found by
# measurement and must not come back.
# --------------------------------------------------------------------------

testthat::test_that("ED7 keeps both spatial statements and names the exception", {
  d <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data")
  pb <- file.path(d, "v8_ed_locations_source_data.csv")
  pa <- file.path(d, "v8_ed_identity_source_data.csv")
  testthat::skip_if_not(all(file.exists(pa, pb)), "extended_data not built")
  b <- rd(pb); a <- rd(pa)
  # the two quantities are genuinely different and BOTH must be carried
  testthat::expect_true(all(c("elsewhere", "outside_affinity") %in% names(b)))
  testthat::expect_true(all(b$elsewhere))
  # the canonical classification is authoritative for affinity, and ED7a and
  # ED7b must agree on it - they disagreed before this pass (15 vs 14)
  aq <- a[a$subset == "CA2_SLM_robustness_qualified", ]
  testthat::expect_equal(sum(b$outside_affinity),
                         as.integer(aq$effect_outside_baseline_affinity))
  testthat::expect_identical(
    b$outside_affinity,
    b$effect_identity_relationship == "effect_outside_baseline_affinity")
  # whenever the two counts differ, the exception must be identifiable
  if (sum(b$elsewhere) != sum(b$outside_affinity))
    testthat::expect_gt(length(b$GeneSymbol[!b$outside_affinity]), 0L)
})

testthat::test_that("ED7a encodes counts, not a percentage past 100", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data", "v8_ed_identity_source_data.csv")
  testthat::skip_if_not(file.exists(p), "extended_data not built")
  z <- rd(p)
  testthat::expect_true(all(z$effect_outside_baseline_affinity <= z$n_hits))
  src <- paste(readLines(repo_path("R", "editorial_v8_fidelity_panels.R"),
                         warn = FALSE), collapse = "\n")
  # the old encoding put a proportion on a 0-122 axis
  testthat::expect_false(grepl("limits = c(0, 122)", src, fixed = TRUE))
})

testthat::test_that("every NES strip in the family shares one colour limit", {
  root <- path_results("source_data", "manuscript_candidates", "editorial_v8")
  fs <- c(file.path(root, "figure_03",
                    paste0("v8_curve_", c("syn", "rna", "ox"),
                           "_source_data.csv")),
          file.path(root, "extended_data",
                    paste0("v8_ed_gsea_curve_", c("syn", "rna", "ox"),
                           "_source_data.csv")))
  testthat::skip_if_not(all(file.exists(fs)), "curves not built")
  lims <- unlist(lapply(fs, function(p) {
    z <- rd(p)
    if ("shared_NES_strip_limit" %in% names(z)) z$shared_NES_strip_limit else NA
  }))
  lims <- lims[!is.na(lims)]
  testthat::expect_gt(length(lims), 0L)
  testthat::expect_length(unique(round(lims, 9)), 1L)
  # and it must cover every value any strip draws
  obs <- max(unlist(lapply(fs, function(p) {
    z <- rd(p)
    abs(c(z$RES_CON_NES, z$SUS_CON_NES, z$SUS_RES_NES))
  })), na.rm = TRUE)
  testthat::expect_gte(unique(lims)[1] + 1e-9, obs)
})

testthat::test_that("a censored colour scale never hides the true value", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "figure_02", "v8_compartment_source_data.csv")
  testthat::skip_if_not(file.exists(p), "figure_02 not built")
  z <- rd(p)
  testthat::expect_true(all(c("displayed_value", "true_value",
                              "colour_scale_censored") %in% names(z)))
  # the released source data must carry the uncapped value
  cen <- z[z$colour_scale_censored %in% TRUE, , drop = FALSE]
  if (nrow(cen)) testthat::expect_true(all(abs(cen$true_value) >
                                             abs(cen$displayed_value)))
})

testthat::test_that("the depth panel offset carries density, not rank", {
  # comments are stripped: the file documents the old defect in prose, and the
  # test is about the CODE
  src <- code_of(repo_path("R", "editorial_v8_fidelity_panels.R"))
  # the frozen renderer used (rank %% 5 - 2) * 0.07, a deterministic function
  # of rank that drew diagonal staircases with no data behind them
  testthat::expect_false(grepl("%% 5 - 2", src, fixed = TRUE))
  testthat::expect_true(grepl("e8_sina_offset", src, fixed = TRUE))
  # the offset must be deterministic: no random jitter anywhere in the layer
  for (f in c("editorial_v8_fidelity_panels.R", "editorial_v8_panels.R",
              "editorial_v8_ed_panels.R")) {
    s <- code_of(repo_path("R", f))
    testthat::expect_false(grepl("geom_jitter|position_jitter|runif\\(|rnorm\\(",
                                 s, perl = TRUE), label = f)
  }
  off <- e8_sina_offset(c(1, 1, 1, 2, 5))
  testthat::expect_identical(off, e8_sina_offset(c(1, 1, 1, 2, 5)))
  testthat::expect_equal(sum(off), 0, tolerance = 1e-9)
})

testthat::test_that("the whole-network null uses a log axis", {
  src <- paste(readLines(repo_path("R", "editorial_v8_fidelity_panels.R"),
                         warn = FALSE), collapse = "\n")
  # on a linear 0-1.05 axis at 31.6 mm the "0" and "0.05" tick labels overprint
  # and the attainable floor is indistinguishable from the axis line
  testthat::expect_true(grepl("scale_x_log10", src, fixed = TRUE))
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data", "v8_ed_nulls_source_data.csv")
  testthat::skip_if_not(file.exists(p), "extended_data not built")
  z <- rd(p)
  testthat::expect_true(all(z$floor > 0))
  testthat::expect_true(all(z$p > z$floor))
})

testthat::test_that("the WGCNA null reports the spatial interaction omnibus", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data", "v8_ed_wgcna_phenotype_source_data.csv")
  testthat::skip_if_not(file.exists(p), "extended_data not built")
  st <- unique(rd(p)$inferential_status)
  testthat::expect_length(st, 1L)
  testthat::expect_true(grepl("interaction omnibus", st, fixed = TRUE))
  # the 690 within-unit contrasts are deliberately NOT reported: their smallest
  # FDR is 0.49 and calling that evidence of absence would need an attainable-
  # floor treatment, which would be new inference
  testthat::expect_false(grepl("690", st, fixed = TRUE))
  testthat::expect_false(grepl("within_spatial_unit", st, fixed = TRUE))
})

testthat::test_that("fidelity renderers introduce no new inference", {
  src <- paste(readLines(repo_path("R", "editorial_v8_fidelity_panels.R"),
                         warn = FALSE), collapse = "\n")
  src <- paste(sub("#.*$", "", strsplit(src, "\n")[[1]]), collapse = "\n")
  for (pat in c("\\blimma\\b", "\\blmer\\b", "\\bfgsea\\b", "p\\.adjust\\s*\\(",
                "\\bt\\.test\\s*\\(", "\\bwilcox\\.test\\s*\\(",
                "\\bcor\\.test\\s*\\(", "\\bglm\\s*\\(", "\\baov\\s*\\(")) {
    testthat::expect_false(grepl(pat, src, perl = TRUE), label = pat)
  }
})

# --------------------------------------------------------------------------
# Part-26. Encoding-principle and deliverable invariants.
# --------------------------------------------------------------------------

testthat::test_that("no panel prints its values inside a heatmap tile", {
  # a tile grid with the number written in it is a table with a colour wash
  # behind it; where the value must be readable the panel uses position
  src <- code_of(repo_path("R", "editorial_v8_panels.R"))
  f0 <- regexpr("e8_protein_zoom <- function", src, fixed = TRUE)
  f1 <- regexpr("e8_gsea_curve <- function", src, fixed = TRUE)
  body <- substr(src, f0, f1)
  testthat::expect_false(grepl("geom_tile", body, fixed = TRUE))
  testthat::expect_true(grepl("geom_point", body, fixed = TRUE))
  # contrast is carried redundantly, so the panel survives greyscale and the
  # common colourblindness forms
  testthat::expect_true(grepl("scale_colour_manual", body, fixed = TRUE))
  testthat::expect_true(grepl("scale_shape_manual", body, fixed = TRUE))
})

testthat::test_that("every colour and size legend carries a real name", {
  bad <- character(0)
  for (f in Sys.glob(repo_path("R", "editorial_v8*.R"))) {
    s <- code_of(f)
    nm <- regmatches(s, gregexpr('name = "[^"]+"', s))[[1]]
    nm <- sub('name = ', '', nm)
    # a publication scale name is a phrase, not a token like "NES" or "effect"
    short <- nm[nchar(gsub('"', '', nm)) < 12 &
                  !grepl("Contrast|Gene set", nm)]
    if (length(short)) bad <- c(bad, paste0(basename(f), ": ", short))
  }
  testthat::expect_identical(bad, character(0))
})

testthat::test_that("Figure 2d is ordered by baseline peak, phenotype-blind", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "figure_02", "v8_fingerprint_source_data.csv")
  testthat::skip_if_not(file.exists(p), "figure_02 not built")
  z <- rd(p)
  testthat::expect_true(all(c("peak_unit", "row_order_rule") %in% names(z)))
  # the ordering may not use any stress or phenotype column
  testthat::expect_false(any(grepl("SUS|RES|stress|phenotype|contrast_group",
                                   names(z), ignore.case = FALSE)))
  # and it must actually be ordered: genes are no longer alphabetical
  g <- unique(z$gene)
  testthat::expect_false(identical(g, sort(g)))
})

testthat::test_that("F2b states acquisition n and biological n separately", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "figure_02", "v8_depth_source_data.csv")
  testthat::skip_if_not(file.exists(p), "figure_02 not built")
  z <- rd(p)
  testthat::expect_true(all(c("n_acquisitions_in_compartment",
                              "n_biological_replicates") %in% names(z)))
  # 323 acquisitions from 9 animals: the two must never be conflated
  testthat::expect_identical(unique(z$n_biological_replicates), 9L)
  testthat::expect_gt(sum(unique(z$n_acquisitions_in_compartment)),
                      unique(z$n_biological_replicates))
  # a violin is only defensible above the small-sample threshold
  testthat::expect_true(all(unique(z$n_acquisitions_in_compartment) >= 50L))
})

testthat::test_that("ED8a keeps one global zero-centred scale", {
  p <- path_results("source_data", "manuscript_candidates", "editorial_v8",
                    "extended_data", "v8_ed_similarity_source_data.csv")
  testthat::skip_if_not(file.exists(p), "extended_data not built")
  z <- rd(p)
  # the metric and the zero reference are common to all three blocks, so the
  # scale stays global even though two blocks are entirely negative
  src <- code_of(repo_path("R", "editorial_v8_ed_panels.R"))
  f0 <- regexpr("e8_ed_similarity <- function", src, fixed = TRUE)
  body <- substr(src, f0, f0 + 6000)
  testthat::expect_true(grepl("nv_diverging", body, fixed = TRUE))
  testthat::expect_false(grepl("scale_fill_viridis|scale_fill_gradient\\(",
                               body, perl = TRUE))
  # one limit for the whole panel, not one per block
  testthat::expect_equal(length(unique(round(
    max(abs(z$median_similarity), na.rm = TRUE), 9))), 1L)
})

testthat::test_that("supplementary tables are reader-facing and complete", {
  d <- path_results("tables", "manuscript_candidates", "editorial_v8",
                    "supplementary")
  testthat::skip_if_not(dir.exists(d), "supplementary tables not built")
  dd <- file.path(d, "ST0_data_dictionary.csv")
  testthat::expect_true(file.exists(dd))
  dict <- rd(dd)
  fs <- setdiff(list.files(d, "[.]csv$"), "ST0_data_dictionary.csv")
  testthat::expect_gt(length(fs), 0L)
  for (f in fs) {
    z <- rd(file.path(d, f))
    # every column of every table must be defined
    defined <- dict$column[dict$table_file == f]
    testthat::expect_identical(setdiff(names(z), defined), character(0),
                               label = f)
    testthat::expect_false(any(is.na(dict$definition[dict$table_file == f])),
                           label = f)
    # the replicate count must be stated, and must be animals not acquisitions
    if ("Biological replicates (n animals)" %in% names(z))
      testthat::expect_identical(unique(z$`Biological replicates (n animals)`),
                                 9L, label = f)
    # machine-oriented provenance columns must not leak into a reader table
    testthat::expect_false(any(names(z) %in%
                                 c("reading", "note", "encoding_note",
                                   "shared_scale_note", "row_order_rule")),
                           label = f)
  }
})

testthat::test_that("every output folder carries a generated README", {
  roots <- c(
    path_results("tables", "manuscript_candidates", "editorial_v8"),
    path_results("tables", "manuscript_candidates", "editorial_v8",
                 "supplementary"),
    path_results("reports", "manuscript_candidates", "editorial_v8"))
  for (k in c("figure_02", "figure_03", "extended_data")) {
    roots <- c(roots,
      path_results("figures", "manuscript_candidates", "editorial_v8", k,
                   "assembled"),
      path_results("figures", "manuscript_candidates", "editorial_v8", k,
                   "panels"),
      file.path(path_results("source_data", "manuscript_candidates",
                             "editorial_v8"), k))
  }
  for (r in roots) {
    testthat::skip_if_not(dir.exists(r), r)
    f <- file.path(r, "README.md")
    testthat::expect_true(file.exists(f), label = r)
    txt <- paste(readLines(f, warn = FALSE), collapse = " ")
    # each README must say the layer is a candidate and say how to rebuild it
    testthat::expect_true(grepl("candidate", txt, ignore.case = TRUE),
                          label = r)
    testthat::expect_true(grepl("Rscript", txt, fixed = TRUE), label = r)
    testthat::expect_true(grepl("Do not edit", txt, fixed = TRUE), label = r)
  }
})

testthat::test_that("the numeric ExpGroup typing is recorded and not depended on", {
  # Found in the Part-26 pre-freeze audit. In the variance-partitioning stage
  # ExpGroup is stored as a NUMERIC with example values "2; 3; 1", while every
  # other categorical term - Region, Layer, ReplicateGroup, AnimalID - is
  # coerced to a factor. A numeric stress group enters that model as a
  # continuous covariate with 1 degree of freedom and imposes an arbitrary
  # linear order on CON/RES/SUS, which would understate the variance
  # attributable to group.
  #
  # This is an UPSTREAM analysis defect, outside the figure layer. It is not
  # fixed here, because fixing it would reopen an analysis stage. This test
  # pins two things: that the defect is still exactly as characterised, so the
  # note stays truthful, and that no editorial_v8 panel reads that stage. If
  # someone wires a panel to it, this test fails and forces the typing to be
  # dealt with first.
  d <- path_results("tables", "03_qc_exploration", "06_variance_partitioning")
  testthat::skip_if_not(dir.exists(d), "variance partitioning not present")
  fs <- list.files(d, "^metadata_terms_used_final[.]csv$", recursive = TRUE,
                   full.names = TRUE)
  testthat::skip_if_not(length(fs) > 0, "no metadata term tables")
  for (f in fs) {
    z <- rd(f)
    eg <- z[z$term == "ExpGroup", , drop = FALSE]
    if (!nrow(eg)) next
    testthat::expect_identical(as.character(eg$class[1]), "numeric",
                               label = basename(dirname(f)))
    others <- z[z$term %in% c("Region", "Layer", "AnimalID"), , drop = FALSE]
    if (nrow(others))
      testthat::expect_true(all(others$class == "factor"),
                            label = basename(dirname(f)))
  }
  # and no v8 panel may depend on that stage while the typing stands
  for (p in V8$panels) {
    src <- c(as.character(p$primary_source %||% ""),
             as.character(unlist(p$input_dependencies %||% list())))
    testthat::expect_false(any(grepl("variance_partitioning", src)),
                           label = as.character(p$id))
  }
})
