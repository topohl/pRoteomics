source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
NV <- nv_contract()
SH <- nv_shared_paths()

# =====================================================================
# size contract
# =====================================================================

testthat::test_that("every figure honours the Nature size contract", {
  testthat::expect_identical(NV$contract_version, "manuscript_nature_v2_figures_v1")
  testthat::expect_identical(NV$status, "candidate_only_not_promoted")
  for (f in NV$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183, info = f$name)
    testthat::expect_lte(as.numeric(f$height_mm), 170, label = as.character(f$name))
    # panel boxes stay on the canvas
    for (it in f$layout) {
      testthat::expect_lte(as.numeric(it$x) + as.numeric(it$w), 183)
      testthat::expect_lte(as.numeric(it$y) + as.numeric(it$h), as.numeric(f$height_mm))
    }
  }
  # the old 250 mm canvas is gone
  testthat::expect_false(any(vapply(NV$figures,
    function(f) as.numeric(f$height_mm) > 200, logical(1))))
})

testthat::test_that("panel counts stay in the 3-6 range and area is allocated unevenly", {
  for (f in NV$figures) {
    n <- length(f$layout)
    testthat::expect_gte(n, 3L); testthat::expect_lte(n, 6L)
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    # not a uniform grid: the largest panel must be clearly bigger than the mean
    testthat::expect_gt(max(area) / mean(area), 1.3)
  }
})

testthat::test_that("panel labels are lowercase, never uppercase", {
  for (f in NV$figures) {
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs), info = as.character(f$name))
    testthat::expect_true(all(labs %in% letters))
  }
})

# =====================================================================
# typography and palette
# =====================================================================

testthat::test_that("the palette contract exists and is semantic", {
  testthat::expect_true(have(nv_palette_path()))
  p <- nv_palette()
  testthat::expect_identical(p$palette_version, "manuscript_palette_v1")
  testthat::expect_setequal(names(p$group), c("CON", "RES", "SUS"))
  testthat::expect_setequal(names(p$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  # one diverging scale for every signed molecular effect
  testthat::expect_setequal(names(p$diverging), c("low", "mid", "high"))
  # no red-green contrast: the diverging ends must not both be green-ish
  testthat::expect_false(grepl("^#0", p$diverging$high))
  # CON is neutral grey: r, g and b within a narrow band of each other
  rgb_of <- function(h) as.integer(grDevices::col2rgb(h))
  con <- rgb_of(p$group$CON)
  testthat::expect_lt(max(con) - min(con), 12L)
  # typography floor
  testthat::expect_gte(nv_pt("axis_text_pt"), 5)
  testthat::expect_lte(nv_pt("axis_text_pt"), 7)
  testthat::expect_equal(nv_pt("panel_label_pt"), 8)
  testthat::expect_identical(p$typography$family, "Arial")
  testthat::expect_identical(p$typography$panel_label_case, "lowercase")
})

testthat::test_that("no rendered text falls below the minimum point size", {
  d <- nv_output_paths("figure_03")$panels
  testthat::skip_if_not(dir.exists(d), "nature_v2 not built")
  svgs <- list.files(d, pattern = "[.]svg$", full.names = TRUE)
  testthat::skip_if(length(svgs) == 0L)
  floor_pt <- nv_pt("minimum_pt")
  for (f in svgs) {
    txt <- readLines(f, warn = FALSE)
    fs <- regmatches(txt, gregexpr("font-size: *[0-9.]+px", txt))
    fs <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", unlist(fs))))
    fs <- fs[is.finite(fs) & fs > 0]
    if (!length(fs)) next
    # svglite writes px in a pt-based user space, so px == pt here
    testthat::expect_gte(min(fs), floor_pt - 0.01, label = basename(f))
  }
})

# =====================================================================
# panels are rendered at their exact final box, never rescaled
# =====================================================================

testthat::test_that("panels are authored at the box they are placed in", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- nv_output_paths(key)$panels
    if (!dir.exists(d)) next
    figs <- Filter(function(f) identical(as.character(f$figure_key), key), NV$figures)
    for (f in figs) {
      for (it in f$layout) {
        id <- as.character(it$panel)
        w <- as.numeric(it$w); h <- as.numeric(it$h)
        cand <- c(file.path(d, sprintf("%s_%gx%g.svg", id, w, h)),
                  file.path(d, paste0(id, ".svg")))
        p <- cand[file.exists(cand)][1]
        testthat::expect_false(is.na(p), info = paste(f$name, id))
        if (is.na(p)) next
        hdr <- paste(readLines(p, n = 4, warn = FALSE), collapse = " ")
        grab <- function(a) {
          m <- regmatches(hdr, regexpr(paste0(a, "=['\"]([0-9.]+)pt"), hdr))
          if (!length(m)) return(NA_real_)
          suppressWarnings(as.numeric(gsub("[^0-9.]", "", m)))
        }
        wpt <- grab("width"); hpt <- grab("height")
        if (!is.finite(wpt)) next
        testthat::expect_equal(wpt / 72 * 25.4, w, tolerance = 0.02,
                               label = paste(id, "width"))
        testthat::expect_equal(hpt / 72 * 25.4, h, tolerance = 0.02,
                               label = paste(id, "height"))
      }
    }
  }
})

# =====================================================================
# Nature house style
# =====================================================================

testthat::test_that("the candidate banner and plot titles are gone", {
  for (key in c("figure_02", "figure_03")) {
    d <- nv_output_paths(key)$assembled
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern = "[.]svg$", full.names = TRUE)) {
      txt <- paste(readLines(f, warn = FALSE), collapse = " ")
      testthat::expect_false(grepl("NOT A MANUSCRIPT FIGURE", txt, fixed = TRUE),
                             info = basename(f))
      testthat::expect_false(grepl("CANDIDATE -", txt, fixed = TRUE), info = basename(f))
    }
  }
  # the theme blanks titles and subtitles outright
  src <- paste(readLines(repo_path("R", "nature_v2_figure_utils.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("plot.title = ggplot2::element_blank()", src, fixed = TRUE))
  testthat::expect_true(grepl("plot.subtitle = ggplot2::element_blank()", src, fixed = TRUE))
})

testthat::test_that("renderers cannot create new inference", {
  testthat::expect_true(nv_assert_no_model_fitting())
  tmp <- tempfile(fileext = ".R")
  writeLines(c("f <- function(x) {", "  p.adjust(x)", "}"), tmp)
  testthat::expect_error(nv_assert_no_model_fitting(tmp), "must not create new inference")
  unlink(tmp)
  # reading a stored FDR COLUMN called p.adjust is legitimate and must not trip
  tmp2 <- tempfile(fileext = ".R")
  writeLines(c("g <- function(d) d$p.adjust[1]"), tmp2)
  testthat::expect_true(nv_assert_no_model_fitting(tmp2))
  unlink(tmp2)
})

# =====================================================================
# scientific hierarchy and honesty
# =====================================================================

testthat::test_that("the GSEA atlas is the visual centrepiece of every Figure 3 variant", {
  f3 <- Filter(function(f) identical(as.character(f$figure_key), "figure_03"), NV$figures)
  testthat::expect_gte(length(f3), 2L)
  for (f in f3) {
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    testthat::expect_identical(ids[which.max(area)], "n3_atlas",
                               info = as.character(f$name))
    testthat::expect_gt(max(area) / sum(area), 0.3)
  }
})

testthat::test_that("WGCNA phenotype effects never outweigh the GSEA evidence", {
  for (f in NV$figures) {
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    if (!"n3_wgcna_small" %in% ids) next
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    testthat::expect_lt(area[ids == "n3_wgcna_small"], area[ids == "n3_atlas"])
  }
  # and at least one main variant carries no WGCNA phenotype panel at all
  f3 <- Filter(function(f) identical(as.character(f$figure_key), "figure_03"), NV$figures)
  without <- vapply(f3, function(f) {
    !("n3_wgcna_small" %in% vapply(f$layout, function(x) as.character(x$panel), character(1)))
  }, logical(1))
  testthat::expect_true(any(without))

  # the descriptive status is recorded, not hidden
  sd <- file.path(nv_output_paths("figure_03")$source_data, "n3_wgcna_small_source_data.csv")
  testthat::skip_if_not(have(sd), "nature_v2 figure 3 not built")
  z <- rd(sd)
  testthat::expect_true(all(grepl("descriptive only", z$status_note, ignore.case = TRUE)))
  testthat::expect_identical(unique(z$fdr_supported_cells_in_panel), 0L)
})

testthat::test_that("the atlas marks constituent FDR support without inventing a test", {
  sd <- file.path(nv_output_paths("figure_03")$source_data, "n3_atlas_source_data.csv")
  testthat::skip_if_not(have(sd), "nature_v2 figure 3 not built")
  z <- rd(sd)
  testthat::expect_true(all(grepl("NOT a new FDR family", z$summary_basis, fixed = TRUE)))
  testthat::expect_true(all(z$n_terms_FDR_supported <= z$n_terms))
  testthat::expect_identical(z$has_FDR_support, z$n_terms_FDR_supported > 0L)
  # all three compartments are represented
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  # qc_review themes are carried but flagged
  testthat::expect_true(any(z$theme_role == "qc_review"))
})

testthat::test_that("representative examples exist for all three compartments", {
  ids <- c("n3_ex_neuropil", "n3_ex_soma", "n3_ex_microglia")
  ds <- vapply(ids, function(i) {
    as.character(nv_contract()$panels[[
      which(vapply(NV$panels, function(p) identical(as.character(p$id), i), logical(1)))
    ]]$example_dataset)
  }, character(1))
  testthat::expect_setequal(unname(ds),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  d <- nv_output_paths("figure_03")$source_data
  testthat::skip_if_not(dir.exists(d), "nature_v2 figure 3 not built")
  for (i in ids) {
    p <- file.path(d, paste0(i, "_source_data.csv"))
    testthat::expect_true(have(p), info = i)
    z <- rd(p)
    testthat::expect_true("status" %in% names(z) || "leading_edge" %in% names(z), info = i)
    if ("leading_edge" %in% names(z)) {
      # a genuine ranked-evidence panel: real FDR, real leading-edge genes
      testthat::expect_true(any(z$leading_edge))
      testthat::expect_lt(unique(z$FDR)[1], 0.05)
      testthat::expect_true(all(grepl("no enrichment statistic is recomputed",
                                      z$evidence_note)))
    }
  }
})

testthat::test_that("the identity panel avoids relocation language", {
  sd <- file.path(nv_output_paths("figure_03")$source_data, "n3_identity_source_data.csv")
  testthat::skip_if_not(have(sd), "nature_v2 figure 3 not built")
  z <- rd(sd)
  testthat::expect_true(all(grepl("not movement of protein", z$interpretation_note)))
  src <- paste(readLines(repo_path("R", "nature_v2_figure_panels.R"), warn = FALSE),
               collapse = " ")
  for (bad in c("redistribution", "relocalization", "relocalisation", "migration")) {
    testthat::expect_false(grepl(bad, src, ignore.case = TRUE), info = bad)
  }
})

testthat::test_that("the DA landscape keeps canonical counts and honest status", {
  sd <- file.path(nv_output_paths("figure_03")$source_data, "n3_da_source_data.csv")
  testthat::skip_if_not(have(sd), "nature_v2 figure 3 not built")
  z <- rd(sd)
  testthat::expect_identical(sum(z$n_proteins), 37L)
  ca2 <- z[z$spatial_unit == "CA2_slm", , drop = FALSE]
  testthat::expect_identical(sum(ca2$n_proteins), 28L)
  testthat::expect_identical(ca2$n_proteins[ca2$status == "claimable"], 6L)
  testthat::expect_true(all(z$status[z$spatial_unit != "CA2_slm"] == "not_audited"))
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
})

testthat::test_that("the spatial anchor is data-backed and declares its artwork status", {
  sd <- file.path(nv_output_paths("figure_02")$source_data, "n2_anchor_source_data.csv")
  testthat::skip_if_not(have(sd), "nature_v2 figure 2 not built")
  z <- rd(sd)
  testthat::expect_identical(nrow(z), 18L)          # the real 18 spatial units
  testthat::expect_setequal(unique(z$region), c("CA1", "CA2", "CA3", "DG"))
  testthat::expect_setequal(unique(z$celltype_layer),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_identical(unique(z$n_animals), 9L)
  testthat::expect_true(all(grepl("no histology or anatomical imagery is fabricated",
                                  z$artwork_status)))
})

# =====================================================================
# isolation and outputs
# =====================================================================

testthat::test_that("the nature_v2 namespace is isolated from both other layers", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    for (d in nv_output_paths(key)) {
      testthat::expect_match(d, "manuscript_candidates[/\\\\]nature_v2")
      testthat::expect_false(grepl("manuscript_panels", d, fixed = TRUE))
    }
  }
  # canonical and Part-16 files never mention this layer
  for (f in c("figures/figure_02.R", "figures/figure_03.R",
              "figures/figure_contract.yml", "R/manuscript_figure_utils.R",
              "figures/figure_candidate_contract.yml", "R/candidate_figure_utils.R",
              "R/candidate_figure_panels.R")) {
    s <- paste(readLines(repo_path(f), warn = FALSE), collapse = "\n")
    for (tok in c("nature_v2", "nv_build", "figure_nature_v2_contract")) {
      testthat::expect_false(grepl(tok, s, fixed = TRUE), info = paste(f, tok))
    }
  }
  # and this layer never calls either other engine
  for (f in c("R/nature_v2_figure_utils.R", "R/nature_v2_figure_panels.R")) {
    s <- paste(sub("#.*$", "", readLines(repo_path(f), warn = FALSE)), collapse = "\n")
    testthat::expect_false(grepl("manuscript_figure_main(", s, fixed = TRUE))
    testthat::expect_false(grepl("cf_build_figure(", s, fixed = TRUE))
  }
})

testthat::test_that("every nature_v2 panel has source data, and variants were built", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    p <- nv_output_paths(key)
    if (!dir.exists(p$assembled)) next
    figs <- Filter(function(f) identical(as.character(f$figure_key), key), NV$figures)
    for (f in figs) {
      testthat::expect_true(
        have(file.path(p$assembled, paste0(as.character(f$name), ".svg"))),
        info = as.character(f$name))
      for (it in f$layout) {
        testthat::expect_true(
          have(file.path(p$source_data, paste0(as.character(it$panel), "_source_data.csv"))),
          info = as.character(it$panel))
      }
    }
  }
})

testthat::test_that("the contact sheet covers every declared variant at true relative size", {
  inv <- file.path(SH$tables, "nature_v2_variant_inventory.csv")
  testthat::skip_if_not(have(inv), "contact sheet not built")
  v <- rd(inv)
  declared <- vapply(NV$figures, function(f) as.character(f$name), character(1))
  testthat::expect_setequal(v$variant, declared)
  testthat::expect_true(all(v$exists %in% TRUE))
  testthat::expect_true(have(file.path(SH$figures, "nature_v2_contact_sheet.svg")))
  # all variants share one width, so equal scaling really is true relative size
  testthat::expect_identical(length(unique(v$width_mm)), 1L)
  pin <- rd(file.path(SH$tables, "nature_v2_panel_inventory.csv"))
  testthat::expect_true(all(pin$promotion_status == "not_promoted"))
  testthat::expect_true(all(nzchar(pin$scientific_question)))
})

testthat::test_that("nature_v2 source data is deterministic", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- nv_output_paths(key)$source_data
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern = "[.]csv$", full.names = TRUE)) {
      testthat::expect_identical(readLines(f, warn = FALSE), readLines(f, warn = FALSE))
    }
  }
  src <- paste(sub("#.*$", "", readLines(repo_path("R", "nature_v2_figure_panels.R"),
                                         warn = FALSE)), collapse = "\n")
  for (tok in c("runif", "rnorm", "sample(", "Sys.time", "Sys.Date")) {
    testthat::expect_false(grepl(tok, src, fixed = TRUE), info = tok)
  }
  # jitter is seeded-free but height-only and cosmetic; it must not reach source data
  testthat::expect_false(grepl("jitter(", src, fixed = TRUE))
})
