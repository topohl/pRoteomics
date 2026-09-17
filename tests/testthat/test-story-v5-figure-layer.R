source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "story_v5_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
S5 <- s5e_contract()

# Comments must never satisfy a test: strip them before scanning source.
code_of <- function(path) {
  ln <- readLines(path, warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}

testthat::test_that("story_v5 is a candidate layer and honours the size contract", {
  testthat::expect_identical(S5$contract_version, s5e_contract_version())
  testthat::expect_identical(S5$status, "candidate_only_not_promoted")
  for (f in S5$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs))
    testthat::expect_identical(anyDuplicated(labs), 0L)
    # every panel must fit inside the page box
    for (it in f$layout) {
      testthat::expect_lte(as.numeric(it$x) + as.numeric(it$w),
                           as.numeric(f$width_mm), label = f$name)
      testthat::expect_lte(as.numeric(it$y) + as.numeric(it$h),
                           as.numeric(f$height_mm), label = f$name)
    }
  }
})

testthat::test_that("no two panels in a figure overlap", {
  for (f in S5$figures) {
    n <- length(f$layout)
    if (n < 2L) next
    for (i in seq_len(n - 1L)) for (j in seq(i + 1L, n)) {
      a <- f$layout[[i]]; b <- f$layout[[j]]
      sep <- as.numeric(a$x) + as.numeric(a$w) <= as.numeric(b$x) ||
        as.numeric(b$x) + as.numeric(b$w) <= as.numeric(a$x) ||
        as.numeric(a$y) + as.numeric(a$h) <= as.numeric(b$y) ||
        as.numeric(b$y) + as.numeric(b$h) <= as.numeric(a$y)
      testthat::expect_true(sep,
        info = sprintf("%s: %s overlaps %s", f$name, a$panel, b$panel))
    }
  }
})

testthat::test_that("the Extended Data family is coherent and gapless", {
  ed <- Filter(function(f) identical(as.character(f$figure_key), "extended_data"),
               S5$figures)
  nm <- vapply(ed, function(f) as.character(f$name), character(1))
  n <- as.integer(sub("^ED([0-9]+)_.*$", "\\1", nm))
  # Part 19 left ED6 and ED8 undeclared; a family with holes in its numbering
  # is not a family, so the sequence must be complete and in order.
  testthat::expect_identical(sort(n), seq_len(length(n)))
  testthat::expect_identical(n, sort(n))
  testthat::expect_gte(length(ed), 8L)
  # every ED panel must resolve to a declared panel with a renderer
  ids <- vapply(S5$panels, function(p) as.character(p$id), character(1))
  for (f in ed) for (it in f$layout) {
    testthat::expect_true(as.character(it$panel) %in% ids, info = it$panel)
  }
})

testthat::test_that("every declared renderer exists and every panel is used", {
  ids <- vapply(S5$panels, function(p) as.character(p$id), character(1))
  testthat::expect_identical(anyDuplicated(ids), 0L)
  used <- unique(unlist(lapply(S5$figures, function(f)
    vapply(f$layout, function(it) as.character(it$panel), character(1)))))
  # a declared panel that no figure uses is dead weight in the contract
  testthat::expect_identical(sort(setdiff(ids, used)), character(0))
})

testthat::test_that("2f and 2g never draw a structurally inapplicable cell", {
  p <- path_results("tables", "manuscript_candidates", "story_v5",
                    "figure2_spatial_validation_hierarchy_audit.csv")
  testthat::skip_if_not(have(p), "validation audit not generated")
  z <- rd(p)
  for (col in c("source_panel", "anatomical_level", "structurally_applicable",
                "reason_inapplicable", "FDR_supported")) {
    testthat::expect_true(col %in% names(z), info = col)
  }
  inapp <- z[!z$structurally_applicable, , drop = FALSE]
  # the inapplicable combinations must be recorded, and each must say WHY
  testthat::expect_gt(nrow(inapp), 0L)
  testthat::expect_true(all(nzchar(inapp$reason_inapplicable)))
  testthat::expect_true(all(is.na(inapp$effect_NES)))

  # the block-diagonal claim that motivates the redesign must actually hold
  k <- z[z$source_panel == "2f_external_kaulich" & z$structurally_applicable, ]
  sig_level <- c(CA1 = "regional", "CA2/3" = "regional", DG = "regional",
                 SLM = "laminar", SO = "laminar", SP = "laminar", SR = "laminar")
  off <- sum(k$anatomical_level != unname(sig_level[k$target]))
  testthat::expect_identical(off, 0L)

  # CA1_SO vs CA3_SO carries a layer token but is a REGIONAL comparison; the
  # audit must classify it from validation_domain, not from the contrast name.
  ca <- k[grepl("CA1.?SO vs CA3.?SO|CA1_SO_vs_CA3_SO", k$contrast), , drop = FALSE]
  if (nrow(ca)) testthat::expect_true(all(ca$anatomical_level == "regional"))
})

testthat::test_that("the coupling panel deduplicates the stored BH family", {
  p <- path_results("source_data", "manuscript_candidates", "story_v5",
                    "extended_data", "v5_ed_coupling_source_data.csv")
  testthat::skip_if_not(have(p), "ED8 coupling panel not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  raw <- repo_path("results", "tables", "08_behavior_physio_coupling",
                   "network_behavior_coupling",
                   "edge_behavior_correlations_fdr_all_tests.csv")
  testthat::skip_if_not(have(raw), "coupling source table absent")
  r <- rd(raw)
  # the stored family has 96 rows but only 48 distinct tests; plotting every
  # row would double-count each estimate
  testthat::expect_lt(nrow(z), nrow(r))
  # every retained row must be a distinct test by its statistical content
  testthat::expect_identical(anyDuplicated(z[, c("Edge", "Outcome", "estimate",
                                                 "p.value")]), 0L)
  # and the collapse must be accounted for, not silent
  testthat::expect_true("n_stored_rows" %in% names(z))
  testthat::expect_identical(sum(z$n_stored_rows), nrow(r))
  testthat::expect_true(all(z$n == 9L))
  # and the null must be stated, not merely implied
  testthat::expect_true(all(grepl("no edge-behaviour association", z$reading, fixed = TRUE)))
  testthat::expect_gt(min(z$p.adj_BH_within_outcome), 0.05)
})

testthat::test_that("no story_v5 renderer fits a model or recomputes inference", {
  for (f in c("story_v5_figure_panels.R", "story_v5_figure_utils.R")) {
    src <- code_of(repo_path("R", f))
    for (tok in c("lmFit(", "eBayes(", "GSEA(", "gseGO(", "enrichGO(",
                  "fgsea(", "blockwiseModules(", "p.adjust(", "t.test(",
                  "wilcox.test(", "cor.test(", "impute.knn(",
                  "normalizeBetweenArrays(")) {
      testthat::expect_false(grepl(tok, src, fixed = TRUE),
                             label = paste(f, tok))
    }
  }
})

testthat::test_that("story_v5 writes only under manuscript_candidates/story_v5", {
  for (f in c("story_v5_figure_02.R", "story_v5_extended_data.R",
              "story_v5_validation_audit.R")) {
    src <- code_of(repo_path("figures", f))
    testthat::expect_false(grepl("figure_contract.yml", src, fixed = TRUE),
                           label = f)
  }
  u <- code_of(repo_path("R", "story_v5_figure_utils.R"))
  testthat::expect_true(grepl("manuscript_candidates", u, fixed = TRUE))
  testthat::expect_true(grepl("story_v5", u, fixed = TRUE))
})

testthat::test_that("panel boxes respect each renderer's natural aspect", {
  # The schematic is a coord_equal drawing: at 126 mm it letterboxes badly and
  # at 175 mm it wastes nearly half its rows. This pins the fix so a later
  # layout edit cannot silently reintroduce the dead band.
  for (f in S5$figures) for (it in f$layout) {
    if (!identical(as.character(it$panel), "v5_schematic")) next
    testthat::expect_lte(as.numeric(it$w), 100,
                         label = sprintf("%s schematic width", f$name))
  }
})

testthat::test_that("rendered panels carry no render_error", {
  for (key in c("figure_02", "extended_data")) {
    p <- path_results("reports", "manuscript_candidates", "story_v5", key,
                      "story_v5_panel_status.csv")
    if (!have(p)) next
    z <- rd(p)
    testthat::expect_identical(sum(z$status != "ok"), 0L, label = key)
  }
})

testthat::test_that("no rendered panel carries a large dead band", {
  testthat::skip_if_not_installed("magick")
  dirs <- c(path_results("figures", "manuscript_candidates", "story_v5",
                         "figure_02", "panels"),
            path_results("figures", "manuscript_candidates", "story_v5",
                         "extended_data", "panels"))
  fs <- unlist(lapply(dirs[dir.exists(dirs)],
                      function(d) list.files(d, "[.]svg$", full.names = TRUE)))
  testthat::skip_if(length(fs) == 0L, "panels not rendered")
  # The Part-19 complaint was unused white space. A panel whose box is far wider
  # or taller than its content letterboxes: the coord_equal schematic wasted 41
  # per cent of its rows at 126 mm and 45 per cent at 175 mm. Measure the
  # largest run of near-empty pixel rows/columns and cap it. The surviving
  # 0.20-0.30 cases are data-driven axis gaps and legend columns, verified by
  # eye; anything past 0.30 is a box that does not fit its content.
  worst <- vapply(fs, function(f) {
    im <- magick::image_background(magick::image_read(f, density = 150),
                                   "white", flatten = TRUE)
    a <- as.integer(magick::image_data(magick::image_convert(im, "gray"), "gray"))
    m <- matrix(a, nrow = dim(a)[2], byrow = TRUE)
    ink <- m < 240
    run <- function(v, thr) {
      r <- rle(v <= thr)
      if (!any(r$values)) return(0L)
      max(r$lengths[r$values])
    }
    max(run(colSums(ink), 0.02 * nrow(ink)) / ncol(ink),
        run(rowSums(ink), 0.02 * ncol(ink)) / nrow(ink))
  }, numeric(1))
  bad <- basename(fs)[worst > 0.30]
  testthat::expect_identical(bad, character(0))
})

testthat::test_that("the v5 PCA pins its legend order", {
  # svp_pca leaves both scales unnamed, so ggplot2 has no guide tie-breaker and
  # the two legends swap at random between runs. s5_pca must set order.
  src <- code_of(repo_path("R", "story_v5_figure_panels.R"))
  testthat::expect_true(grepl("guide_legend(order = 1)", src, fixed = TRUE))
  testthat::expect_true(grepl("guide_legend(order = 2)", src, fixed = TRUE))
  ids <- vapply(S5$panels, function(p) as.character(p$id), character(1))
  rnd <- vapply(S5$panels, function(p) as.character(p$renderer), character(1))
  testthat::expect_false("svp_pca" %in% rnd)
})
