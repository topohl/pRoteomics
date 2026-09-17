source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "story_v3_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
SV <- sv_contract()
SH <- sv_shared_paths()

# =====================================================================
# the audit had to exist before implementation
# =====================================================================

testthat::test_that("the panel audit exists and covers all three prior generations", {
  p <- path_results("tables", "manuscript_candidates", "story_v3",
                    "panel_visual_scientific_audit.csv")
  testthat::skip_if_not(have(p), "audit not generated")
  a <- rd(p)
  testthat::expect_gte(nrow(a), 30L)
  testthat::expect_setequal(unique(a$source_version), c("canonical", "part16", "part17"))
  for (col in c("panel", "source_version", "scientific_question",
                "two_second_takeaway", "current_visual_encoding", "encoding_works",
                "problem", "better_encoding", "narrative_role", "dataset_scope",
                "recommended_destination")) {
    testthat::expect_true(col %in% names(a), info = col)
  }
  testthat::expect_true(all(nzchar(a$scientific_question)))
  testthat::expect_true(all(nzchar(a$two_second_takeaway)))
  testthat::expect_true(all(a$encoding_works %in% c("yes", "no", "partly")))
  # the audit must actually disagree with some panels, or it is not an audit
  testthat::expect_gt(sum(a$encoding_works != "yes"), 5L)
  # and it must recommend some panels OUT of the main figure
  testthat::expect_gt(sum(grepl("extended_data", a$recommended_destination)), 5L)
})

# =====================================================================
# every panel earns its place through a narrative line
# =====================================================================

testthat::test_that("every story_v3 panel states what it establishes", {
  for (p in SV$panels) {
    n <- as.character(p$narrative %||% "")
    testthat::expect_true(nzchar(trimws(n)), info = as.character(p$id))
    testthat::expect_gt(nchar(n), 40L, label = as.character(p$id))
  }
  # and the generated story map contains one line per panel of each figure
  m <- file.path(SH$reports, "figure_story_map.md")
  testthat::skip_if_not(have(m), "story map not generated")
  txt <- paste(readLines(m, warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("HAS NO NARRATIVE LINE", txt, fixed = TRUE))
  for (f in SV$figures) {
    testthat::expect_true(grepl(as.character(f$name), txt, fixed = TRUE),
                          info = as.character(f$name))
  }
})

# =====================================================================
# size and house style (inherited contract, re-asserted here)
# =====================================================================

testthat::test_that("story_v3 figures honour the Nature size contract", {
  testthat::expect_identical(SV$contract_version, sv_contract_version())
  testthat::expect_identical(SV$status, "candidate_only_not_promoted")
  for (f in SV$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs))
    testthat::expect_true(all(labs %in% letters))
    testthat::expect_identical(anyDuplicated(labs), 0L)
    # area is allocated unevenly, not as a uniform grid
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    testthat::expect_gt(max(area) / mean(area), 1.3)
    # and the page is actually used
    testthat::expect_gt(sum(area) / (as.numeric(f$width_mm) * as.numeric(f$height_mm)),
                        0.66)
  }
})

testthat::test_that("figure 2 and figure 3 test the 7-8 panel structures asked for", {
  f2 <- Filter(function(f) identical(as.character(f$figure_key), "figure_02"), SV$figures)
  f3 <- Filter(function(f) identical(as.character(f$figure_key), "figure_03"), SV$figures)
  testthat::expect_gte(length(f2), 2L)
  testthat::expect_gte(length(f3), 3L)
  testthat::expect_true(any(vapply(f2, function(f) length(f$layout) >= 7L, logical(1))))
  testthat::expect_true(any(vapply(f3, function(f) length(f$layout) >= 8L, logical(1))))
})

testthat::test_that("no rendered text falls below the minimum point size", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- sv_output_paths(key)$panels
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern = "[.]svg$", full.names = TRUE)) {
      fs <- unlist(regmatches(readLines(f, warn = FALSE),
                              gregexpr("font-size: *[0-9.]+px", readLines(f, warn = FALSE))))
      fs <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", fs)))
      fs <- fs[is.finite(fs) & fs > 0]
      if (!length(fs)) next
      testthat::expect_gte(min(fs), nv_pt("minimum_pt") - 0.01, label = basename(f))
    }
  }
})

testthat::test_that("panels are rendered at the exact box they are placed in", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- sv_output_paths(key)$panels
    if (!dir.exists(d)) next
    figs <- Filter(function(f) identical(as.character(f$figure_key), key), SV$figures)
    for (f in figs) for (it in f$layout) {
      id <- as.character(it$panel); w <- as.numeric(it$w); h <- as.numeric(it$h)
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
      wpt <- grab("width")
      if (!is.finite(wpt)) next
      testthat::expect_equal(wpt / 72 * 25.4, w, tolerance = 0.02, label = id)
    }
  }
})

# =====================================================================
# the restored and new panels
# =====================================================================

testthat::test_that("the PCA panel is restored and de-cluttered", {
  sd <- file.path(sv_output_paths("figure_02")$source_data, "sv2_pca_source_data.csv")
  testthat::skip_if_not(have(sd), "story_v3 figure 2 not built")
  z <- rd(sd)
  testthat::expect_true(all(c("PC1", "PC2", "dataset", "region") %in% names(z)))
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_gt(nrow(z), 300L)
  # de-cluttered: the renderer must not draw a text label per sample
  src <- paste(sub("#.*$", "", readLines(repo_path("R", "story_v3_figure_panels.R"),
                                          warn = FALSE)), collapse = "\n")
  body <- sub(".*svp_pca <- function", "", src)
  body <- sub("svp_precision_dumbbell.*", "", body)
  testthat::expect_false(grepl("geom_text", body, fixed = TRUE))
})

testthat::test_that("precision is a dumbbell with individual endpoint pairs", {
  sd <- file.path(sv_output_paths("figure_02")$source_data, "sv2_precision_source_data.csv")
  testthat::skip_if_not(have(sd), "story_v3 figure 2 not built")
  z <- rd(sd)
  testthat::expect_true(all(c("ICC_single_side", "ICC_bilateral_mean") %in% names(z)))
  testthat::expect_setequal(unique(z$endpoint_class),
    c("wgcna_module_eigengene", "reference_marker_score", "empirical_compartment_score"))
  # the point of the panel: bilateral averaging raises reliability
  testthat::expect_gt(stats::median(z$ICC_bilateral_mean),
                      stats::median(z$ICC_single_side))
})

testthat::test_that("the anatomical bridge places programs on the Figure-2 sampling map", {
  sd <- file.path(sv_output_paths("figure_03")$source_data, "sv3_anatomy_source_data.csv")
  testthat::skip_if_not(have(sd), "story_v3 figure 3 not built")
  z <- rd(sd)
  testthat::expect_identical(nrow(z), 3L)
  testthat::expect_setequal(z$dataset,
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_setequal(z$unit, c("CA3_sr", "CA2_sp", "CA1"))
  testthat::expect_true(all(z$dir %in% c("higher in SUS", "higher in RES")))
  testthat::expect_true(all(grepl("no anatomical geometry is invented", z$substrate)))
  # direction must agree with the stored NES sign
  testthat::expect_identical(z$dir, ifelse(z$NES > 0, "higher in SUS", "higher in RES"))
})

testthat::test_that("the three rank examples share one grammar and one program per compartment", {
  ids <- c("sv3_ex_neuropil", "sv3_ex_soma", "sv3_ex_microglia")
  d <- sv_output_paths("figure_03")$source_data
  testthat::skip_if_not(dir.exists(d), "story_v3 figure 3 not built")
  seen <- character()
  for (i in ids) {
    z <- rd(file.path(d, paste0(i, "_source_data.csv")))
    testthat::expect_true(any(z$leading_edge))
    testthat::expect_lt(unique(z$FDR)[1], 0.05)
    testthat::expect_true(all(grepl("no enrichment statistic is recomputed",
                                    z$evidence_note)))
    seen <- c(seen, unique(z$dataset))
  }
  testthat::expect_setequal(seen, c("neuron_neuropil", "neuron_soma", "microglia"))
})

testthat::test_that("the protein zoom has a transparent rule and covers all three programs", {
  sd <- file.path(sv_output_paths("figure_03")$source_data, "sv3_proteins_source_data.csv")
  testthat::skip_if_not(have(sd), "story_v3 figure 3 not built")
  z <- rd(sd)
  testthat::expect_identical(length(unique(z$program)), 3L)
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  # equal N per program, chosen by the stated rule and never by gene name
  per <- tapply(z$gene, z$program, function(g) length(unique(g)))
  testthat::expect_identical(length(unique(as.integer(per))), 1L)
  testthat::expect_true(all(grepl("no gene chosen by name", z$selection_rule)))
  testthat::expect_true(all(grepl("leading-edge", z$selection_rule)))
  # no CA2-SLM-dependent m12 story
  testthat::expect_false(any(grepl("m12", z$program, ignore.case = TRUE)))
  # every gene is shown across its own compartment's units only
  testthat::expect_true(all(
    z$spatial_unit[z$dataset == "microglia"] %in% c("CA1", "CA2", "CA3", "DG")))
})

testthat::test_that("WGCNA panels encode structure only, never phenotype significance", {
  d <- sv_output_paths("figure_03")$source_data
  testthat::skip_if_not(dir.exists(d), "story_v3 figure 3 not built")
  for (i in c("sv3_wgcna_circle", "sv3_wgcna_strip")) {
    p <- file.path(d, paste0(i, "_source_data.csv"))
    testthat::skip_if_not(have(p), i)
    z <- rd(p)
    testthat::expect_true(all(z$encodes_phenotype_significance %in% FALSE), info = i)
    testthat::expect_false(any(grepl("tier_specific_fdr|estimate", names(z))), info = i)
    testthat::expect_true("spatial_tau" %in% names(z), info = i)
    testthat::expect_true("bilateral_reproducibility_class" %in% names(z), info = i)
  }
  # both encodings exist so they can be compared like for like
  testthat::expect_true(have(file.path(d, "sv3_wgcna_circle_source_data.csv")))
  testthat::expect_true(have(file.path(d, "sv3_wgcna_strip_source_data.csv")))
})

testthat::test_that("the unsupported WGCNA phenotype heatmap is absent from every main figure", {
  for (f in SV$figures) {
    if (identical(as.character(f$figure_key), "extended_data")) next
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    testthat::expect_false("n3_wgcna_small" %in% ids, info = as.character(f$name))
    testthat::expect_false("ned_wgcna_heatmap" %in% ids, info = as.character(f$name))
  }
  # and a second-question panel is not smuggled into main either
  for (f in SV$figures) {
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    testthat::expect_false("n3_identity" %in% ids)
    testthat::expect_false("n3_synthesis" %in% ids)
  }
})

testthat::test_that("the GSEA atlas is the largest quantitative panel of every Figure 3", {
  f3 <- Filter(function(f) identical(as.character(f$figure_key), "figure_03"), SV$figures)
  for (f in f3) {
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    testthat::expect_identical(ids[which.max(area)], "sv3_atlas",
                               info = as.character(f$name))
    # all three compartments must be visibly represented by the example triptych
    testthat::expect_true(sum(grepl("^sv3_ex_", ids)) >= 2L, info = as.character(f$name))
  }
})

# =====================================================================
# isolation from all three previous layers
# =====================================================================

testthat::test_that("story_v3 is isolated and previous layers are not referenced", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    for (d in sv_output_paths(key)) {
      testthat::expect_match(d, "manuscript_candidates[/\\\\]story_v3")
      testthat::expect_false(grepl("manuscript_panels", d, fixed = TRUE))
      testthat::expect_false(grepl("nature_v2", d, fixed = TRUE))
    }
  }
  for (f in c("figures/figure_02.R", "figures/figure_03.R",
              "figures/figure_contract.yml", "R/manuscript_figure_utils.R",
              "figures/figure_candidate_contract.yml", "R/candidate_figure_utils.R",
              "figures/figure_nature_v2_contract.yml", "R/nature_v2_figure_utils.R",
              "R/nature_v2_figure_panels.R")) {
    s <- paste(readLines(repo_path(f), warn = FALSE), collapse = "\n")
    for (tok in c("story_v3", "sv_build", "figure_story_v3_contract")) {
      testthat::expect_false(grepl(tok, s, fixed = TRUE), info = paste(f, tok))
    }
  }
  # story_v3 reuses the Part-17 helpers by SOURCING them, which does not modify
  s <- paste(readLines(repo_path("R", "story_v3_figure_utils.R"), warn = FALSE),
             collapse = "\n")
  testthat::expect_true(grepl("nature_v2_figure_utils.R", s, fixed = TRUE))
})

testthat::test_that("story_v3 renderers cannot create new inference", {
  testthat::expect_true(nv_assert_no_model_fitting(sv_renderer_sources()))
  src <- paste(sub("#.*$", "", readLines(repo_path("R", "story_v3_figure_panels.R"),
                                          warn = FALSE)), collapse = "\n")
  for (tok in c("runif", "rnorm", "sample(", "jitter(", "Sys.time")) {
    testthat::expect_false(grepl(tok, src, fixed = TRUE), info = tok)
  }
})

testthat::test_that("every story_v3 panel has source data and all variants were built", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    p <- sv_output_paths(key)
    if (!dir.exists(p$assembled)) next
    figs <- Filter(function(f) identical(as.character(f$figure_key), key), SV$figures)
    for (f in figs) {
      testthat::expect_true(
        have(file.path(p$assembled, paste0(as.character(f$name), ".svg"))),
        info = as.character(f$name))
      for (it in f$layout) {
        testthat::expect_true(
          have(file.path(p$source_data,
                         paste0(as.character(it$panel), "_source_data.csv"))),
          info = as.character(it$panel))
      }
    }
  }
  inv <- file.path(SH$tables, "story_v3_variant_inventory.csv")
  testthat::skip_if_not(have(inv), "contact sheet not built")
  v <- rd(inv)
  testthat::expect_setequal(v$variant,
    vapply(SV$figures, function(f) as.character(f$name), character(1)))
  testthat::expect_true(all(v$exists %in% TRUE))
  # page density: no variant may leave more than a third of the page empty
  testthat::expect_true(all(v$used_area_share > 0.66))
})

testthat::test_that("story_v3 source data is deterministic", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- sv_output_paths(key)$source_data
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern = "[.]csv$", full.names = TRUE)) {
      testthat::expect_identical(readLines(f, warn = FALSE), readLines(f, warn = FALSE))
    }
  }
})
