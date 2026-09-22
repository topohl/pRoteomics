# Determinism of ggrepel label placement in publication figures.
#
# ggrepel places labels by a randomised simulation and defaults to seed = NA.
# Two renders of identical data therefore produced different label and leader
# coordinates, which is how two variants of
# Fig_RES_SUS_divergence_publication.svg came to differ while their source-data
# CSV was byte-identical and their 22 labels were the same set. A figure whose
# bytes move when nothing scientific moved cannot be checked against a freeze.
#
# These tests pin determinism rather than a particular arrangement: exact
# coordinate snapshots would break on a legitimate graphics-library upgrade,
# whereas "render A == render B" is the property that actually matters.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "plotting_nature.R"))

has_pkgs <- function() {
  all(vapply(c("ggplot2", "ggrepel", "svglite"), requireNamespace,
             logical(1), quietly = TRUE))
}

# The plotted data is built once, with a fixed seed, so every render below uses
# byte-identical input. Building it inside the fixture would reset the global
# RNG before each render, and an unseeded repel layer is deterministic when the
# RNG state happens to match - which is why the fixture must not do that. In a
# real script the layer is reached with whatever state upstream work left
# behind, and that is what made the two divergence figures differ.
set.seed(1)
REPEL_DATA <- data.frame(x = runif(40L), y = runif(40L),
                         lab = paste0("L", sprintf("%02d", seq_len(40L))))

repel_fixture <- function(seed = NATURE_REPEL_SEED, d = REPEL_DATA) {
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point(size = 1)
  layer <- if (is.na(seed)) {
    ggrepel::geom_text_repel(ggplot2::aes(label = lab), size = 2, max.overlaps = 12)
  } else {
    ggrepel::geom_text_repel(ggplot2::aes(label = lab), size = 2,
                             max.overlaps = 12, seed = seed)
  }
  p + layer
}

render_to <- function(p, path) {
  suppressWarnings(ggplot2::ggsave(path, p, width = 74, height = 86, units = "mm"))
  path
}
svg_geometry <- function(path) {
  l <- readLines(path, warn = FALSE)
  l[grepl("^<(text|line|polyline|circle)", l)]
}
svg_labels <- function(path) {
  l <- readLines(path, warn = FALSE)
  sort(gsub("^>|</text>$", "",
            unlist(regmatches(l, gregexpr(">[^<>]+</text>", l)))))
}

testthat::test_that("the repel seed contract exists in exactly one place", {
  testthat::expect_true(exists("NATURE_REPEL_SEED"))
  testthat::expect_true(is.integer(NATURE_REPEL_SEED))
  testthat::expect_identical(length(NATURE_REPEL_SEED), 1L)
  src <- readLines(repo_path("R", "plotting_nature.R"), warn = FALSE)
  testthat::expect_identical(length(grep("^NATURE_REPEL_SEED <-", src)), 1L)
  # not derived from a path or a timestamp
  defn <- src[grep("^NATURE_REPEL_SEED <-", src)]
  testthat::expect_false(grepl("Sys.time|Sys.Date|repo_path|getwd|basename", defn))
})

testthat::test_that("every active ggrepel layer carries a seed", {
  # Parsed from each call rather than grepped from a line window: a window
  # bleeds into the next layer and credits it with a neighbour's seed.
  extract_call <- function(lines, start) {
    txt <- paste(lines[start:min(length(lines), start + 60L)], collapse = "\n")
    tok <- regexpr("(ggrepel::)?geom_(text|label)_repel", txt)
    if (tok > 0) txt <- substr(txt, tok, nchar(txt))
    chars <- strsplit(txt, "")[[1]]
    open <- which(chars == "(")[1]
    if (is.na(open)) return(NA_character_)
    depth <- 0L; inq <- ""
    for (i in seq(open, length(chars))) {
      ch <- chars[[i]]
      if (nzchar(inq)) {
        if (ch == inq && chars[[max(i - 1L, 1L)]] != "\\") inq <- ""
      } else if (ch %in% c("'", "\"")) inq <- ch
      else if (ch == "(") depth <- depth + 1L
      else if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) return(substr(txt, 1L, i))
      }
    }
    NA_character_
  }

  roots <- c(repo_path("analysis"), repo_path("R"), repo_path("tools"))
  files <- unlist(lapply(roots[dir.exists(roots)], list.files,
                         pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE))
  unseeded <- character(0); total <- 0L
  for (f in files) {
    l <- readLines(f, warn = FALSE)
    for (i in grep("geom_text_repel|geom_label_repel", l)) {
      txt <- extract_call(l, i)
      if (is.na(txt)) { unseeded <- c(unseeded, paste0(basename(f), ":", i, " [unparsed]")); next }
      e <- tryCatch(parse(text = txt, keep.source = FALSE)[[1]], error = function(err) NULL)
      if (is.null(e) || !is.call(e)) { unseeded <- c(unseeded, paste0(basename(f), ":", i)); next }
      total <- total + 1L
      if (!"seed" %in% names(e)) unseeded <- c(unseeded, paste0(basename(f), ":", i))
    }
  }
  # 11, not the 12 Phase 6H.6 recorded. One of those twelve layers was
  # compare_go_enrichment.R:2452, inside the script's unreachable tail - a
  # repel layer that could never render. Phase 6H.10 archived that region, so
  # the live count is 11. No 6H.6 conclusion changes: the other eleven are
  # live and all twelve were seeded, including the archived one.
  testthat::expect_gte(total, 11L)
  testthat::expect_identical(length(unseeded), 0L,
    info = paste("unseeded repel layers:", paste(unseeded, collapse = ", ")))
})

testthat::test_that("an unseeded layer is genuinely nondeterministic (the defect)", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  # Identical data, different global RNG state at render time - the situation a
  # long script actually creates. This is what the fix prevents.
  set.seed(11); a <- render_to(repel_fixture(seed = NA), file.path(tmp, "a.svg"))
  set.seed(99); b <- render_to(repel_fixture(seed = NA), file.path(tmp, "b.svg"))
  testthat::expect_false(identical(svg_geometry(a), svg_geometry(b)))
  # the labels themselves never varied - only their placement
  testthat::expect_identical(svg_labels(a), svg_labels(b))
})

testthat::test_that("a seeded layer renders identically for identical input", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  set.seed(11); a <- render_to(repel_fixture(), file.path(tmp, "a.svg"))
  set.seed(99); b <- render_to(repel_fixture(), file.path(tmp, "b.svg"))
  testthat::expect_identical(svg_geometry(a), svg_geometry(b))
  testthat::expect_identical(svg_labels(a), svg_labels(b))
  # svglite embeds no timestamp, so byte identity is achievable and is asserted
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))
})

testthat::test_that("a seeded layer is identical across differing global RNG state", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  set.seed(11); a <- render_to(repel_fixture(), file.path(tmp, "canonical.svg"))
  set.seed(99); b <- render_to(repel_fixture(), file.path(tmp, "validation.svg"))
  # This is the property that makes a canonical run and a validation-only run
  # agree: the layer seed decides placement, not whatever RNG state the script
  # happened to reach by that point.
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))
  testthat::expect_identical(svg_geometry(a), svg_geometry(b))
})

testthat::test_that("the output path does not leak into the rendered figure", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  d1 <- file.path(tmp, "canonical_scope"); d2 <- file.path(tmp, "validation_proposed_scope")
  dir.create(d1); dir.create(d2)
  a <- render_to(repel_fixture(), file.path(d1, "fig.svg"))
  b <- render_to(repel_fixture(), file.path(d2, "fig.svg"))
  # Two different destinations, byte-identical output: scope identity cannot
  # affect geometry, which was the shape of the Phase 6H.5C discrepancy.
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))
  testthat::expect_false(any(grepl("validation_proposed_scope", readLines(b, warn = FALSE))))
})

testthat::test_that("seeding a layer leaves the global RNG stream alone", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  set.seed(4242)
  before <- .Random.seed
  invisible(render_to(repel_fixture(), file.path(tmp, "rng.svg")))
  testthat::expect_true(identical(.Random.seed, before))

  # and the downstream stream is unperturbed
  set.seed(4242)
  invisible(render_to(repel_fixture(), file.path(tmp, "rng2.svg")))
  after_render <- runif(3)
  set.seed(4242)
  expected <- runif(3)
  testthat::expect_equal(after_render, expected)
})

testthat::test_that("the seed governs layout only, never labels or selection", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/ggrepel/svglite unavailable")
  tmp <- withr::local_tempdir()
  a <- render_to(repel_fixture(seed = 1L), file.path(tmp, "s1.svg"))
  b <- render_to(repel_fixture(seed = 2L), file.path(tmp, "s2.svg"))
  # different seed moves placement, proving the seed is live ...
  testthat::expect_false(identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b))))
  # ... but the label set is untouched, so no scientific content depends on it
  testthat::expect_identical(svg_labels(a), svg_labels(b))
  ## 40 repel labels plus the axis tick text; the count is incidental, the
  ## point is that the two seeds produce the same label set.
  testthat::expect_gte(length(svg_labels(a)), 40L)
})

testthat::test_that("the accepted canonical divergence figure was not rerendered", {
  f <- path_results("figures", "04_differential_expression_enrichment",
                    "compareGO_spatial_atlas",
                    "Fig_RES_SUS_divergence_publication.svg")
  testthat::skip_if_not(file.exists(f), "canonical divergence figure not present")
  # Phase 6H.6 changes future renders only. The accepted figure keeps its
  # arbitrary-but-correct layout; rerendering it would churn the frozen payload
  # for no scientific gain.
  testthat::expect_identical(file.size(f), 43462)
  sd <- path_results("source_data", "04_differential_expression_enrichment",
                     "compareGO_spatial_atlas",
                     "source_data_RES_SUS_divergence_publication.csv")
  testthat::skip_if_not(file.exists(sd), "divergence source data not present")
  testthat::expect_identical(
    unname(tools::sha256sum(sd)),
    "e1ed9fc9009762069bf0eb9e76fdb9051ddba7e5311330f710f65812655708d1")
  testthat::expect_identical(file.size(sd), 135885)
})

testthat::test_that("the accepted figure package and its freeze are untouched", {
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure export manifest not present")
  m <- utils::read.csv(mp, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(m), 5582L)
  testthat::expect_identical(
    unname(tools::sha256sum(mp)),
    "0fd0c9ed9febbc05ed6928b7fc6bfffdab845daba6255d8b47c607dee5d3c02c")
})
