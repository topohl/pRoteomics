# Guards for the Phase 5B promotion of final_truth_v9 to manuscript Figures 2 and 3.
#
# Three defects held the promotion back and each has a test here, because each
# was invisible to the checks that existed at the time: a column misregistration
# that no audit measured, a provenance dependency that no test declared, and a
# scientific constant hard-coded into a legend where re-running could not
# correct it.

repo_rel <- function(...) file.path(testthat::test_path("..", ".."), ...)
contract <- function() {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  yaml::read_yaml(repo_rel("figures", "figure_contract.yml"))
}

testthat::test_that("exactly one canonical generation is declared per manuscript figure", {
  testthat::skip_if_not_installed("yaml")
  y <- contract()
  testthat::expect_identical(y$contract_version,
                             "manuscript_figures_v3_final_truth_v9_promoted")

  for (k in c("02", "03")) {
    f <- y$figures[[k]]
    testthat::expect_identical(as.character(f$canonical_generation), "final_truth_v9")
    testthat::expect_identical(as.character(f$superseded_status), "SUPERSEDED_FOR_MANUSCRIPT")
    testthat::expect_false(isTRUE(f$rendering_repository_computes_statistics))
    testthat::expect_identical(as.character(f$layout_mode), "absolute")
  }

  # Figure 2 is a-h and Figure 3 is a-i, which is the structure Results 2 and 3
  # were written against.
  testthat::expect_identical(
    vapply(y$figures[["02"]]$panels, function(p) as.character(p$id), character(1)),
    paste0("2", letters[1:8]))
  testthat::expect_identical(
    vapply(y$figures[["03"]]$panels, function(p) as.character(p$id), character(1)),
    paste0("3", letters[1:9]))

  # Figure 1 is untouched by the promotion.
  testthat::expect_identical(
    vapply(y$figures[["01"]]$panels, function(p) as.character(p$id), character(1)),
    paste0("1", letters[1:4]))
})

testthat::test_that("PB-02: no manuscript figure depends on the superseded export namespace", {
  testthat::skip_if_not_installed("yaml")
  src <- paste(readLines(repo_rel("figures", "figure_final_truth_v9_contract.yml"),
                         warn = FALSE), collapse = "\n")
  # The producing contract must not read anything the manuscript layer writes,
  # which is what made promoting-and-retiring impossible before.
  testthat::expect_false(grepl("results/source_data/manuscript/figure_0", src, fixed = TRUE))

  y <- contract()
  for (k in c("02", "03")) for (p in y$figures[[k]]$panels) {
    deps <- as.character(unlist(p$input_dependencies %||% character()))
    for (d in deps) {
      testthat::expect_false(grepl("^results/source_data/manuscript/figure_0", d),
        info = paste(p$id, "still depends on the v2 export namespace:", d))
    }
  }
})

testthat::test_that("PB-02: the promoted panels resolve to canonical stage outputs", {
  testthat::skip_if_not_installed("yaml")
  v9 <- yaml::read_yaml(repo_rel("figures", "figure_final_truth_v9_contract.yml"))
  byid <- setNames(v9$panels, vapply(v9$panels, function(p) as.character(p$id), character(1)))

  testthat::expect_identical(as.character(byid[["v9_depth"]]$primary_source),
                             "data/raw/pg_matrix/quicksearch.stats.annotated.xlsx")
  testthat::expect_identical(as.character(byid[["v9_pca"]]$primary_source),
    "results/tables/03_qc_exploration/00b_joint_compartment_qc/global/joint_primary_pca_scores.csv")
  for (id in c("v9_depth", "v9_pca")) {
    testthat::expect_true(file.exists(repo_rel(as.character(byid[[id]]$primary_source))),
                          info = paste("canonical source missing for", id))
  }
})

testthat::test_that("PB-01: the DAP track and the atlas share column geometry by construction", {
  # The registration contract is declared in the contract and honoured in code.
  # Before the repair the coupling existed only as a comment, and the two panels
  # drifted by up to 10.5 mm while every audit reported the figure clean.
  testthat::skip_if_not_installed("yaml")
  v9 <- yaml::read_yaml(repo_rel("figures", "figure_final_truth_v9_contract.yml"))
  byid <- setNames(v9$panels, vapply(v9$panels, function(p) as.character(p$id), character(1)))
  testthat::expect_identical(
    as.character(byid[["v9_atlas"]]$shares_column_geometry_with), "v9_dap_track")

  code <- paste(readLines(repo_rel("R", "final_truth_v9_panels.R"), warn = FALSE),
                collapse = "\n")
  testthat::expect_true(grepl("shares_column_geometry_with", code, fixed = TRUE))

  # The two boxes must be the same width and x, or no gutter convention can make
  # their columns line up.
  f3 <- Filter(function(x) x$name == "F3_NATURE_FINAL_V9", v9$figures)[[1]]
  boxes <- setNames(f3$layout, vapply(f3$layout, function(x) as.character(x$panel), character(1)))
  testthat::expect_identical(as.numeric(boxes[["v9_dap_track"]]$x),
                             as.numeric(boxes[["v9_atlas"]]$x))
  testthat::expect_identical(as.numeric(boxes[["v9_dap_track"]]$w),
                             as.numeric(boxes[["v9_atlas"]]$w))
})

testthat::test_that("PB-01: the rendered columns actually register", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  pan <- path_results("figures", "manuscript_candidates", "final_truth_v9", "figure_03", "panels")
  a <- file.path(pan, "v9_dap_track.svg"); b <- file.path(pan, "v9_atlas.svg")
  testthat::skip_if_not(file.exists(a) && file.exists(b), "Figure 3 panels not rendered here")

  inset <- function(f) {
    s <- readLines(f, warn = FALSE)
    w <- as.numeric(sub(".*width='([0-9.]+)pt'.*", "\\1", grep("viewBox", s, value = TRUE)[1]))
    r <- unlist(regmatches(s, gregexpr("<rect[^>]*/>", s)))
    rx <- suppressWarnings(as.numeric(sub(".*\\bx='([-0-9.]+)'.*", "\\1", r)))
    rw <- suppressWarnings(as.numeric(sub(".*\\bwidth='([-0-9.]+)'.*", "\\1", r)))
    k <- !is.na(rx) & !is.na(rw) & rw > 1 & rw < 60
    c(left = min(rx[k]), right = w - max(rx[k] + rw[k]), pitch = stats::median(rw[k]))
  }
  ia <- inset(a); ib <- inset(b)
  # Equal to well under a tenth of a column. Before the repair the pitches were
  # 18.720 and 17.540 pt and the right insets 21.54 and 51.24 pt.
  testthat::expect_lt(abs(ia[["pitch"]] - ib[["pitch"]]), 0.05)
  testthat::expect_lt(abs(ia[["left"]] - ib[["left"]]), 0.5)
  testthat::expect_lt(abs(ia[["right"]] - ib[["right"]]), 0.5)
})

testthat::test_that("PB-03: the eps-floor disclosure is derived, not hard-coded", {
  src <- paste(readLines(repo_rel("figures", "final_truth_v9_legends.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("f9_eps_floor_disclosure", src, fixed = TRUE))
  # The stale literal must not come back as a live string. It survives only
  # inside the comment that explains why it was removed.
  live <- sub("(?s)# This disclosure used to be.*?derived\\.", "", src, perl = TRUE)
  testthat::expect_false(grepl("90 of the 851", live, fixed = TRUE))
  testthat::expect_false(grepl('"851"', live, fixed = TRUE))
  # The floor itself is declared once and checked against the data at run time.
  testthat::expect_true(grepl("GSEA_EPS <- 1e-10", src, fixed = TRUE))
  testthat::expect_true(grepl("is not the declared floor", src, fixed = TRUE))
})

testthat::test_that("SR-21: the shipped Figure 3a is the QC-aware panel", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  sh <- path_results("manuscript", "figure_3", "panels", "figure_03a.svg")
  testthat::skip_if_not(file.exists(sh), "Figure 3 not exported here")

  lab <- local({
    s <- readLines(sh, warn = FALSE)
    t <- unlist(regmatches(s, gregexpr("<text[^>]*>[^<]*</text>", s)))
    sub("^<text[^>]*>(.*)</text>$", "\\1", t)
  })
  # Two rows, and the second one is what the pre-QC panel never had.
  testthat::expect_true(any(grepl("FDR-supported DAPs", lab, fixed = TRUE)))
  testthat::expect_true(any(grepl("Robustness-qualified", lab, fixed = TRUE)))
  testthat::expect_true(any(grepl("(37 total)", lab, fixed = TRUE)))
  testthat::expect_true(any(grepl("(15 total)", lab, fixed = TRUE)))
  nums <- suppressWarnings(as.integer(lab[grepl("^[0-9]+$", lab)]))
  testthat::expect_identical(sum(nums, na.rm = TRUE), 52L)   # 37 + 15
  testthat::expect_true(28L %in% nums)
  testthat::expect_true(6L %in% nums)
  testthat::expect_false(any(grepl("hotspot", lab, ignore.case = TRUE)))

  sd <- path_results("source_data", "manuscript", "figure_03", "figure_03a_source_data.csv")
  testthat::skip_if_not(file.exists(sd), "Figure 3a source data not exported here")
  d <- utils::read.csv(sd, stringsAsFactors = FALSE)
  testthat::expect_true(all(c("canonical", "claimable") %in% names(d)))
  testthat::expect_identical(sum(d$canonical), 37L)
  testthat::expect_identical(sum(d$claimable), 15L)
})

testthat::test_that("every manuscript figure-panel reference resolves", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  y <- contract()
  declared <- unlist(lapply(c("01", "02", "03"), function(k)
    vapply(y$figures[[k]]$panels, function(p) as.character(p$id), character(1))))

  draft <- readLines(repo_rel("manuscript", "manuscript_draft.md"), warn = FALSE)
  hits <- unlist(regmatches(draft, gregexpr(
    "Fig[.][ ]?[123][a-i](([,–-])[a-i])*", draft, perl = TRUE)))
  testthat::expect_gt(length(hits), 0L)

  for (hh in unique(hits)) {
    fign <- sub("^Fig[.][ ]?([123]).*$", "\\1", hh)
    tail_ <- sub("^Fig[.][ ]?[123]", "", hh)
    ls_ <- regmatches(tail_, gregexpr("[a-i]", tail_))[[1]]
    idx <- match(ls_, letters)
    if (grepl("[–-]", tail_) && length(idx) >= 2L && !anyNA(idx)) {
      ls_ <- letters[seq(min(idx), max(idx))]
    }
    for (L in ls_) {
      testthat::expect_true(paste0(fign, L) %in% declared,
        info = paste("manuscript cites", hh, "-> undeclared panel", paste0(fign, L)))
    }
  }
})
