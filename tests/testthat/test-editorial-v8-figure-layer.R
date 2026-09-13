# Part-24 contract tests for the editorial_v8 layer.
#
# The scientific architecture is fixed; these tests guard the things Part 24 is
# actually allowed to change - geometry, scale sharing, in-panel prose and the
# vector export - and hard-stop if the layer starts reopening analysis.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "integration_utils.R"))
source(testthat::test_path("..", "..", "R", "nature_v2_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "spatial_grammar_utils.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "editorial_v8_export.R"))

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
