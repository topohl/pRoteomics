source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "integration_utils.R"))
source(testthat::test_path("..", "..", "R", "nature_v2_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "spatial_grammar_utils.R"))
source(testthat::test_path("..", "..", "R", "nature_final_v7_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
V7 <- s7e_contract()

code_of <- function(path) {
  ln <- readLines(path, warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}

testthat::test_that("nature_final_v7 is a candidate layer within the size contract", {
  testthat::expect_identical(V7$contract_version, s7e_contract_version())
  testthat::expect_identical(V7$status, "candidate_only_not_promoted")
  testthat::expect_identical(as.numeric(V7$font_floor_pt), 5)
  for (f in V7$figures) {
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

testthat::test_that("the main figures carry the intended panel sequences", {
  f2 <- Filter(function(f) identical(f$name, "F2_NATURE_FINAL_V7"), V7$figures)[[1]]
  f3 <- Filter(function(f) identical(f$name, "F3_NATURE_FINAL_V7"), V7$figures)[[1]]
  testthat::expect_identical(
    vapply(f2$layout, function(x) as.character(x$label), character(1)),
    letters[1:8])
  testthat::expect_identical(
    vapply(f3$layout, function(x) as.character(x$label), character(1)),
    letters[1:9])
  # PCA is explicitly retained in the main figure
  testthat::expect_true("v7_pca" %in%
    vapply(f2$layout, function(x) as.character(x$panel), character(1)))
  # all three GSEA curves and all three protein zooms stay main
  f3p <- vapply(f3$layout, function(x) as.character(x$panel), character(1))
  testthat::expect_true(all(c("v7_curve_syn", "v7_curve_rna", "v7_curve_ox",
                              "v7_prot_syn", "v7_prot_rna", "v7_prot_ox") %in% f3p))
  # the DAP track is the first panel and sits directly on top of the atlas
  testthat::expect_identical(as.character(f3$layout[[1]]$panel), "v7_dap_track")
  a <- f3$layout[[1]]; b <- f3$layout[[2]]
  testthat::expect_identical(as.numeric(a$x), as.numeric(b$x))
  testthat::expect_identical(as.numeric(a$w), as.numeric(b$w))
  # no gutter: they must read as one block
  testthat::expect_identical(as.numeric(a$y) + as.numeric(a$h), as.numeric(b$y))
})

testthat::test_that("Figure-2 area hierarchy matches the declared priority", {
  f2 <- Filter(function(f) identical(f$name, "F2_NATURE_FINAL_V7"), V7$figures)[[1]]
  area <- stats::setNames(
    vapply(f2$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1)),
    vapply(f2$layout, function(x) as.character(x$panel), character(1)))
  # the fingerprint must be the largest panel in the figure
  testthat::expect_identical(names(which.max(area)), "v7_fingerprint")
  # Part 22 found PCA (rank 7) outranking compartment identity (rank 3) in area
  testthat::expect_lt(area[["v7_pca"]], area[["v7_compartment"]])
  # depth is the least dominant panel
  testthat::expect_identical(names(which.min(area)), "v7_depth")
  # and no minor technical panel may exceed a central scientific one
  minor <- max(area[c("v7_depth", "v7_pca")])
  central <- min(area[c("v7_fingerprint", "v7_compartment", "v7_bilateral_main")])
  testthat::expect_lt(minor, central)
})

testthat::test_that("no rendered panel carries text below the 5 pt floor", {
  floor_pt <- as.numeric(V7$font_floor_pt)
  seen <- 0L
  for (key in c("figure_02", "figure_03", "extended_data")) {
    d <- path_results("figures", "manuscript_candidates", "nature_final_v7", key,
                      "panels")
    if (!dir.exists(d)) next
    for (f in list.files(d, "[.]svg$", full.names = TRUE)) {
      txt <- paste(readLines(f, warn = FALSE), collapse = "")
      m <- regmatches(txt, gregexpr("font-size: [0-9.]+px", txt))[[1]]
      if (!length(m)) next
      seen <- seen + 1L
      v <- as.numeric(sub("px", "", sub("font-size: ", "", m)))
      testthat::expect_gte(min(v), floor_pt - 1e-6, label = basename(f))
    }
  }
  testthat::skip_if(seen == 0L, "panels not rendered")
  testthat::expect_gt(seen, 20L)
})

testthat::test_that("no main and ED panel is the same graphic", {
  p <- path_results("tables", "manuscript_candidates", "nature_final_v7",
                    "main_ed_panel_hierarchy_audit.csv")
  testthat::skip_if_not(have(p), "hierarchy audit not generated")
  z <- rd(p)
  testthat::expect_false(any(z$byte_identical %in% TRUE))
  testthat::expect_false(any(grepl("INVALID", z$hierarchy_valid)))
  # every analysis that appears on ONE side only must carry a documented reason
  one <- z[is.na(z$ed_panel) | !nzchar(as.character(z$ed_panel)), , drop = FALSE]
  testthat::expect_true(all(nzchar(one$ed_unique_content)))
  # and the pairs that do appear on both sides must be genuinely different
  both <- z[!is.na(z$ed_panel) & nzchar(as.character(z$ed_panel)), , drop = FALSE]
  testthat::expect_gte(nrow(both), 5L)
  testthat::expect_true(all(both$byte_identical == FALSE))
})

testthat::test_that("the DAP track reports the canonical sparse burden", {
  p <- path_results("source_data", "manuscript_candidates", "nature_final_v7",
                    "figure_03", "v7_dap_track_source_data.csv")
  testthat::skip_if_not(have(p), "DAP track not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  testthat::expect_identical(nrow(z), 18L)
  testthat::expect_identical(sum(z$canonical), 37L)
  testthat::expect_identical(sum(z$claimable), 15L)
  testthat::expect_identical(sum(z$canonical == 0L), 12L)
  ca2 <- z[z$unit == "CA2_slm", ]
  testthat::expect_identical(ca2$canonical, 28L)
  testthat::expect_identical(ca2$claimable, 6L)
  # the microglia burden Part 22 flagged: CA1 and CA3 are zero, CA2 is 3
  mg <- z[z$dataset == "microglia", ]
  testthat::expect_identical(mg$canonical[mg$unit == "CA1"], 0L)
  testthat::expect_identical(mg$canonical[mg$unit == "CA3"], 0L)
  testthat::expect_identical(mg$canonical[mg$unit == "CA2"], 3L)
})

testthat::test_that("the bridge never overstates an unsupported trajectory", {
  p <- path_results("source_data", "manuscript_candidates", "nature_final_v7",
                    "figure_03", "v7_bridge_source_data.csv")
  testthat::skip_if_not(have(p), "bridge not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  testthat::expect_identical(nrow(z), 9L)
  testthat::expect_identical(sort(unique(z$contrast)),
                             c("RES - CON", "SUS - CON", "SUS - RES"))
  # only the microglia exemplar has all three contrasts FDR-supported
  ox <- z[z$key == "oxphos", ]
  testthat::expect_true(all(ox$FDR < 0.05))
  testthat::expect_true(all(grepl("graded", ox$shape)))
  for (k in c("synaptic", "rna")) {
    w <- z[z$key == k, ]
    testthat::expect_gt(w$FDR[w$contrast == "RES - CON"], 0.05)
    testthat::expect_true(all(grepl("susceptibility-associated", w$shape)))
    # the word "divergent" must never be asserted for these two
    testthat::expect_false(any(grepl("divergent", w$shape, ignore.case = TRUE)))
  }
})

testthat::test_that("each protein zoom is one program at a usable width", {
  f3 <- Filter(function(f) identical(f$name, "F3_NATURE_FINAL_V7"), V7$figures)[[1]]
  w <- vapply(f3$layout, function(x) as.numeric(x$w), numeric(1))
  ids <- vapply(f3$layout, function(x) as.character(x$panel), character(1))
  prot <- w[ids %in% c("v7_prot_syn", "v7_prot_rna", "v7_prot_ox")]
  # Part 22 found the combined panel structurally broken at 35 mm
  testthat::expect_true(all(prot >= 50))
  for (k in c("syn", "rna", "ox")) {
    p <- path_results("source_data", "manuscript_candidates", "nature_final_v7",
                      "figure_03", sprintf("v7_prot_%s_source_data.csv", k))
    if (!have(p)) next
    z <- rd(p)
    if (identical(as.character(z$status[1]), "render_error")) next
    testthat::expect_identical(sort(unique(as.character(z$contrast))),
                               sort(c("RES−CON", "SUS−CON", "SUS−RES")))
    testthat::expect_gte(length(unique(z$gene)), 6L)
    testthat::expect_true(all(grepl("no gene chosen by name", z$selection_rule)))
  }
})

testthat::test_that("no v7 renderer creates new inference", {
  for (f in c("nature_final_v7_panels.R", "nature_final_v7_figure3_panels.R",
              "nature_final_v7_ed_panels.R", "nature_final_v7_figure_utils.R")) {
    src <- code_of(repo_path("R", f))
    for (tok in c("lmFit(", "eBayes(", "topTable(", "GSEA(", "gseGO(",
                  "enrichGO(", "fgsea(", "blockwiseModules(", "TOMsimilarity(",
                  "p.adjust(", "t.test(", "wilcox.test(", "cor.test(", "aov(",
                  "lmer(", "impute.knn(", "normalizeBetweenArrays(")) {
      testthat::expect_false(grepl(tok, src, fixed = TRUE), label = paste(f, tok))
    }
  }
})

testthat::test_that("v7 writes only under manuscript_candidates/nature_final_v7", {
  for (f in c("nature_final_v7_figure_02.R", "nature_final_v7_figure_03.R",
              "nature_final_v7_extended_data.R", "nature_final_v7_wireframe.R",
              "nature_final_v7_hierarchy_audit.R",
              "nature_final_v7_contact_sheet.R")) {
    src <- code_of(repo_path("figures", f))
    testthat::expect_false(grepl("figure_contract.yml", src, fixed = TRUE),
                           label = f)
  }
  u <- code_of(repo_path("R", "nature_final_v7_figure_utils.R"))
  testthat::expect_true(grepl("nature_final_v7", u, fixed = TRUE))
})

testthat::test_that("the microglia compartment is never called cell-intrinsic", {
  for (f in c("nature_final_v7_panels.R", "nature_final_v7_figure3_panels.R",
              "nature_final_v7_ed_panels.R")) {
    src <- paste(readLines(repo_path("R", f), warn = FALSE), collapse = "\n")
    for (bad in c("microglial proteome", "microglia-specific",
                  "purified microglia")) {
      testthat::expect_false(grepl(bad, src, fixed = TRUE),
                             label = paste(f, bad))
    }
  }
})

testthat::test_that("every hard QA condition passed at build time", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    p <- path_results("reports", "manuscript_candidates", "nature_final_v7", key,
                      "nature_final_v7_qa.csv")
    if (!have(p)) next
    z <- rd(p)
    testthat::expect_true(all(z$status == "PASS"), label = key)
  }
})

testthat::test_that("rendered panels carry no render_error", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    p <- path_results("reports", "manuscript_candidates", "nature_final_v7", key,
                      "nature_final_v7_panel_status.csv")
    if (!have(p)) next
    z <- rd(p)
    testthat::expect_identical(sum(z$status != "ok"), 0L, label = key)
  }
})
