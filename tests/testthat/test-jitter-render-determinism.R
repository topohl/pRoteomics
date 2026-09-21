# Determinism of stochastic point placement in publication figures.
#
# position_jitter(), position_jitterdodge() and geom_jitter() displace points by
# a random draw and default to seed = NA, which means they draw from the global
# RNG stream at render time. Two renders of identical data therefore produced
# different point coordinates - the same class of defect as the unseeded
# ggrepel layers in test-ggrepel-render-determinism.R, but a separate
# mechanism, so seeding one does nothing for the other.
#
# These tests pin determinism, not a particular arrangement: a coordinate
# snapshot would break on a legitimate graphics-library upgrade, whereas
# "render A == render B" is the property that actually matters. They also pin
# the jitter widths, because the scientific reading of a jittered panel depends
# on the spread being what the author chose.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "plotting_nature.R"))

has_pkgs <- function() {
  all(vapply(c("ggplot2", "svglite"), requireNamespace, logical(1), quietly = TRUE))
}

# Built once with a fixed seed so every render below uses byte-identical input.
# Building it inside the fixture would reset the global RNG before each render,
# and an unseeded jitter layer is deterministic when the incoming RNG state
# happens to match - which would hide the defect instead of demonstrating it.
set.seed(1)
JITTER_DATA <- data.frame(
  g = rep(c("A", "B", "C"), each = 30L),
  y = runif(90L),
  lab = paste0("P", sprintf("%02d", seq_len(90L))))

jitter_fixture <- function(seed = NATURE_JITTER_SEED, d = JITTER_DATA,
                           repel_seed = NA) {
  pos <- if (is.na(seed)) {
    ggplot2::position_jitter(width = 0.16, height = 0)
  } else {
    ggplot2::position_jitter(width = 0.16, height = 0, seed = seed)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(g, y)) +
    ggplot2::geom_point(position = pos, size = 1)
  if (!is.na(repel_seed) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = d[seq_len(20L), ],
                                      ggplot2::aes(label = lab), size = 2,
                                      max.overlaps = 12, seed = repel_seed)
  }
  p
}

render_to <- function(p, path) {
  suppressWarnings(ggplot2::ggsave(path, p, width = 90, height = 70, units = "mm"))
  path
}
svg_geometry <- function(path) {
  l <- readLines(path, warn = FALSE)
  l[grepl("^<(circle|use|text|line|polyline|g )", l)]
}
svg_labels <- function(path) {
  l <- readLines(path, warn = FALSE)
  sort(gsub("^>|</text>$", "",
            unlist(regmatches(l, gregexpr(">[^<>]+</text>", l)))))
}

# Paren-matching extraction, as in the ggrepel test: a fixed line window bleeds
# into the next layer and credits it with a neighbour's seed.
#
# Sliced from the token's exact PARSE COLUMN rather than from a regex search of
# the window. Searching finds the first occurrence, so two calls of the same
# token on one physical line - which this repository has, because several QC
# scripts put a whole plot on one line - would both resolve to the first, and
# the second call's seed would never be examined.
extract_jitter_call <- function(lines, start, col) {
  txt <- paste(lines[start:min(length(lines), start + 40L)], collapse = "\n")
  txt <- substr(txt, col, nchar(txt))
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

# A seed argument only helps if its VALUE pins the stream. ggplot2's defective
# default is seed = NA, and seed = NULL is equally unpinned, so a check for the
# mere presence of the argument would score `position_jitter(width = .12,
# seed = NA)` as compliant while it reproduces the exact defect this phase
# removed. Anything else - a literal, or the shared constant - does pin it.
pins_stream <- function(seed_text) {
  !is.na(seed_text) &&
    !seed_text %in% c("NA", "NULL", "NA_integer_", "NA_real_", "NA_character_")
}

# For geom_jitter, the seed lives in the position argument. Read it by parsing
# that expression rather than by searching its deparsed text.
delegates_seed <- function(e, nms) {
  if (!"position" %in% nms) return(FALSE)
  pe <- e[["position"]]
  if (!is.call(pe)) return(FALSE)
  pn <- names(pe)
  if (is.null(pn) || !"seed" %in% pn) return(FALSE)
  pins_stream(paste(deparse(pe[["seed"]]), collapse = ""))
}

# Call discovery is driven by getParseData() so that prose mentioning
# position_jitter - the explanatory comment in plotting_nature.R, or this
# header - is never counted as a call.
jitter_calls <- function() {
  tokens <- c("position_jitterdodge", "position_jitter", "geom_jitter")
  roots <- c(repo_path("analysis"), repo_path("R"), repo_path("tools"))
  files <- unlist(lapply(roots[dir.exists(roots)], list.files,
                         pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE))
  out <- list()
  for (f in files) {
    pd <- tryCatch(utils::getParseData(parse(f, keep.source = TRUE)),
                   error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) next
    hits <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% tokens, , drop = FALSE]
    if (!nrow(hits)) next
    l <- readLines(f, warn = FALSE)
    rel <- sub(paste0("^", repo_path(), "/"), "", gsub("\\\\", "/", f))
    for (k in seq_len(nrow(hits))) {
      txt <- extract_jitter_call(l, hits$line1[k], hits$col1[k])
      e <- if (!is.na(txt)) tryCatch(parse(text = txt, keep.source = FALSE)[[1]],
                                     error = function(err) NULL) else NULL
      nms <- if (is.null(e)) character(0) else names(e)
      dep <- function(k2) if (k2 %in% nms) paste(deparse(e[[k2]]), collapse = "") else NA_character_
      out[[length(out) + 1L]] <- list(
        file = rel, line = hits$line1[k], call = hits$text[k], parsed = !is.null(e),
        # a geom_jitter delegating to an explicit seeded position_jitter is
        # deterministic; that inner call is also listed in its own right.
        # The delegated seed is read by PARSING the position argument, not by
        # searching its text for "seed" - that substring also matches a symbol
        # named `unseeded_pos`, and would score the defect as compliant.
        has_seed = pins_stream(dep("seed")) ||
          (hits$text[k] == "geom_jitter" && delegates_seed(e, nms)),
        seed_value = dep("seed"),
        width = dep("width"), height = dep("height"))
    }
  }
  out
}

# The only calls allowed to remain unseeded: scripts that write nothing into a
# root the manuscript exporter scans, so their renders can never become
# publication candidates. Listing them explicitly means a new unseeded jitter
# anywhere else fails this test rather than passing unnoticed.
UNSEEDED_ALLOWLIST <- c(
  "analysis/integration/quantify_candidate_network_position.R",
  "analysis/integration/test_network_behaviour_coupling.R",
  "analysis/spatial_validation/validate_network_workbook.R")

testthat::test_that("the jitter seed contract exists in exactly one place", {
  testthat::expect_true(exists("NATURE_JITTER_SEED"))
  testthat::expect_true(is.integer(NATURE_JITTER_SEED))
  testthat::expect_identical(length(NATURE_JITTER_SEED), 1L)
  src <- readLines(repo_path("R", "plotting_nature.R"), warn = FALSE)
  testthat::expect_identical(length(grep("^NATURE_JITTER_SEED <-", src)), 1L)
  defn <- src[grep("^NATURE_JITTER_SEED <-", src)]
  # not derived from a path, a clock or the working directory
  testthat::expect_false(grepl("Sys.time|Sys.Date|repo_path|getwd|basename", defn))
})

testthat::test_that("jitter and repel are separately controllable constants", {
  # Deliberately distinct names. If they were one constant, changing a repel
  # layout would silently move every jittered panel in the repository.
  testthat::expect_true(exists("NATURE_REPEL_SEED"))
  src <- readLines(repo_path("R", "plotting_nature.R"), warn = FALSE)
  testthat::expect_identical(length(grep("^NATURE_REPEL_SEED <-", src)), 1L)
  testthat::expect_identical(length(grep("^NATURE_JITTER_SEED <-", src)), 1L)
})

testthat::test_that("no publication-facing jitter layer is unseeded", {
  calls <- jitter_calls()
  testthat::expect_gte(length(calls), 26L)
  testthat::expect_true(all(vapply(calls, `[[`, logical(1), "parsed")))

  unseeded <- Filter(function(c) !c$has_seed, calls)
  offenders <- vapply(unseeded, function(c) c$file, character(1))
  unexpected <- setdiff(unique(offenders), UNSEEDED_ALLOWLIST)
  testthat::expect_identical(length(unexpected), 0L,
    info = paste("unseeded jitter outside the allowlist:",
                 paste(unexpected, collapse = ", ")))

  # and the allowlist has not quietly grown to cover a real figure script
  testthat::expect_lte(length(unseeded), 5L)
})

testthat::test_that("the compliance guard rejects a seed that pins nothing", {
  # Without this, the guard degenerates into "is the word seed present", and
  # `seed = NA` - ggplot2's defective default, written out explicitly - would
  # read as compliant. These are the values that do NOT pin the stream.
  testthat::expect_false(pins_stream(NA_character_))
  testthat::expect_false(pins_stream("NA"))
  testthat::expect_false(pins_stream("NULL"))
  testthat::expect_false(pins_stream("NA_integer_"))
  testthat::expect_true(pins_stream("1"))
  testthat::expect_true(pins_stream("NATURE_JITTER_SEED"))

  # and the delegated form is parsed, not substring-matched: a position given
  # as a symbol whose NAME contains "seed" must not count as seeded
  sym <- parse(text = "geom_jitter(position = unseeded_pos)",
               keep.source = FALSE)[[1]]
  testthat::expect_false(delegates_seed(sym, names(sym)))
  naseed <- parse(text = "geom_jitter(position = position_jitter(width = .1, seed = NA))",
                  keep.source = FALSE)[[1]]
  testthat::expect_false(delegates_seed(naseed, names(naseed)))
  good <- parse(text = "geom_jitter(position = position_jitter(width = .1, seed = 7L))",
                keep.source = FALSE)[[1]]
  testthat::expect_true(delegates_seed(good, names(good)))
})

testthat::test_that("two calls of one token on a single line are read separately", {
  # Several QC scripts write an entire plot on one physical line. Slicing the
  # call text from a regex search would resolve both calls to the first and
  # silently inherit the first one's seed for the second.
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "twocalls.R")
  writeLines(paste0("p <- ggplot2::ggplot(d) + ",
                    "ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.11, seed = 1L)) + ",
                    "ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.44))"), f)
  pd <- utils::getParseData(parse(f, keep.source = TRUE))
  hits <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "position_jitter", , drop = FALSE]
  testthat::expect_identical(nrow(hits), 2L)
  l <- readLines(f, warn = FALSE)
  parsed <- lapply(seq_len(nrow(hits)), function(k) {
    e <- parse(text = extract_jitter_call(l, hits$line1[k], hits$col1[k]),
               keep.source = FALSE)[[1]]
    list(width = paste(deparse(e[["width"]]), collapse = ""),
         seeded = pins_stream(if ("seed" %in% names(e))
           paste(deparse(e[["seed"]]), collapse = "") else NA_character_))
  })
  # distinct widths prove the two calls were extracted independently ...
  testthat::expect_identical(vapply(parsed, function(p) p$width, character(1)),
                             c("0.11", "0.44"))
  # ... and the unseeded second call is reported as unseeded
  testthat::expect_identical(vapply(parsed, function(p) p$seeded, logical(1)),
                             c(TRUE, FALSE))
})

testthat::test_that("every script using the seed also sources the library", {
  calls <- jitter_calls()
  users <- unique(vapply(Filter(function(c) c$has_seed, calls),
                         function(c) c$file, character(1)))
  for (f in users) {
    src <- readLines(file.path(repo_path(), f), warn = FALSE)
    uses_constant <- any(grepl("NATURE_JITTER_SEED", src, fixed = TRUE))
    if (!uses_constant) next  # a pre-existing literal seed needs no library
    testthat::expect_true(any(grepl("plotting_nature", src, fixed = TRUE)),
      info = paste(f, "uses NATURE_JITTER_SEED without sourcing plotting_nature.R"))
  }
})

testthat::test_that("the repaired jitter widths are exactly what they were", {
  # Phase 6H.7 was allowed to add a seed and nothing else. The spread of a
  # jittered panel is an authored aesthetic; if a later edit changes it, the
  # figure's visual reading changes and this test should say so.
  expected <- list(
    "analysis/enrichment/run_ewce_celltype_enrichment.R"  = c("0"),
    "analysis/qc/assess_joint_compartment_quality.R"      = c("0.15"),
    "analysis/qc/assess_marker_rank_abundance.R"          = c("0.16"),
    "analysis/qc/assess_replicate_consistency.R"          = c("0.12"),
    "analysis/qc/assess_sample_quality.R"                 = c("0.12"),
    "analysis/qc/export_marker_traits.R"                  = c("0.12"),
    "analysis/qc/summarize_marker_detectability.R"        = c("0.15", "0.15", "0.15"),
    "analysis/qc/summarize_missingness.R"                 = c("0.12"),
    "analysis/wgcna/render_module_figures.R"              = c("0.11"),
    # deparsed, so the literal 0.10 in the source reads back as "0.1"
    "analysis/wgcna/score_module_activity.R"              = c("0.1", "0.1"))
  # height matters as much as width: a jittered boxplot overlay with a vertical
  # component reads as scatter in the y direction, which is the measured axis
  expected_height <- list(
    "analysis/enrichment/run_ewce_celltype_enrichment.R"  = c("0.22"),
    "analysis/qc/assess_joint_compartment_quality.R"      = character(0),
    "analysis/qc/assess_marker_rank_abundance.R"          = character(0),
    "analysis/qc/assess_replicate_consistency.R"          = character(0),
    "analysis/qc/assess_sample_quality.R"                 = c("0"),
    "analysis/qc/export_marker_traits.R"                  = character(0),
    "analysis/qc/summarize_marker_detectability.R"        = character(0),
    "analysis/qc/summarize_missingness.R"                 = character(0),
    "analysis/wgcna/render_module_figures.R"              = c("0"),
    "analysis/wgcna/score_module_activity.R"              = c("0", "0"))

  calls <- jitter_calls()
  for (f in names(expected)) {
    got <- Filter(function(c) c$file == f && !is.na(c$width), calls)
    widths <- sort(as.character(vapply(got, function(c) c$width, character(1))))
    testthat::expect_identical(widths, sort(expected[[f]]), info = f)

    goth <- Filter(function(c) c$file == f && !is.na(c$height), calls)
    heights <- sort(as.character(vapply(goth, function(c) c$height, character(1))))
    testthat::expect_identical(heights, sort(expected_height[[f]]), info = paste(f, "height"))
  }

  # the one position_jitterdodge carries its own parameter names, so it is not
  # covered by the width sweep above
  jd <- readLines(repo_path("analysis", "qc", "assess_sample_quality.R"), warn = FALSE)
  hit <- grep("position_jitterdodge(", jd, fixed = TRUE)
  testthat::expect_identical(length(hit), 1L)
  testthat::expect_true(grepl("jitter.width = 0.12", jd[hit], fixed = TRUE))
  testthat::expect_true(grepl("dodge.width = 0.65", jd[hit], fixed = TRUE))

  # and the one call whose jitter is vertical rather than horizontal
  ew <- readLines(repo_path("analysis", "enrichment", "run_ewce_celltype_enrichment.R"),
                  warn = FALSE)
  hit2 <- grep("position_jitter(width = 0, height = 0.22", ew, fixed = TRUE)
  testthat::expect_identical(length(hit2), 1L)
})

testthat::test_that("an unseeded jitter layer is genuinely nondeterministic", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  tmp <- withr::local_tempdir()
  # Identical data, different global RNG state at render time - the situation a
  # long script actually creates. This is the defect the seeds prevent.
  set.seed(11); a <- render_to(jitter_fixture(seed = NA), file.path(tmp, "a.svg"))
  set.seed(99); b <- render_to(jitter_fixture(seed = NA), file.path(tmp, "b.svg"))
  testthat::expect_false(identical(svg_geometry(a), svg_geometry(b)))
  # the data and its labels never varied - only point placement
  testthat::expect_identical(svg_labels(a), svg_labels(b))
})

testthat::test_that("a seeded jitter layer is identical across differing RNG state", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  tmp <- withr::local_tempdir()
  set.seed(11); a <- render_to(jitter_fixture(), file.path(tmp, "canonical.svg"))
  set.seed(99); b <- render_to(jitter_fixture(), file.path(tmp, "validation.svg"))
  # svglite embeds no timestamp, so byte identity is achievable and asserted
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))
  testthat::expect_identical(svg_geometry(a), svg_geometry(b))
})

testthat::test_that("seeded jitter and seeded repel compose in one panel", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  testthat::skip_if_not(requireNamespace("ggrepel", quietly = TRUE), "ggrepel unavailable")
  tmp <- withr::local_tempdir()
  # Both mechanisms draw from the same global stream when unseeded, so a panel
  # carrying both could in principle be order-dependent. It is not: each layer
  # seed is local to that layer.
  set.seed(11)
  a <- render_to(jitter_fixture(repel_seed = NATURE_REPEL_SEED), file.path(tmp, "a.svg"))
  set.seed(99)
  b <- render_to(jitter_fixture(repel_seed = NATURE_REPEL_SEED), file.path(tmp, "b.svg"))
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))

  # and with the two mechanisms given different seeds from each other
  set.seed(1); c1 <- render_to(jitter_fixture(seed = 11L, repel_seed = 22L),
                               file.path(tmp, "c1.svg"))
  set.seed(2); c2 <- render_to(jitter_fixture(seed = 11L, repel_seed = 22L),
                               file.path(tmp, "c2.svg"))
  testthat::expect_identical(unname(tools::sha256sum(c1)), unname(tools::sha256sum(c2)))
})

testthat::test_that("seeding a jitter layer leaves the global RNG stream alone", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  tmp <- withr::local_tempdir()
  set.seed(4242)
  before <- .Random.seed
  invisible(render_to(jitter_fixture(), file.path(tmp, "rng.svg")))
  testthat::expect_true(identical(.Random.seed, before))

  # and the downstream stream is unperturbed, so no statistic shifts
  set.seed(4242)
  invisible(render_to(jitter_fixture(), file.path(tmp, "rng2.svg")))
  after_render <- runif(3)
  set.seed(4242)
  expected <- runif(3)
  testthat::expect_equal(after_render, expected)
})

testthat::test_that("the jitter seed governs placement only, never data", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  tmp <- withr::local_tempdir()
  a <- render_to(jitter_fixture(seed = 1L), file.path(tmp, "s1.svg"))
  b <- render_to(jitter_fixture(seed = 2L), file.path(tmp, "s2.svg"))
  # a different seed moves the points, proving the seed is live ...
  testthat::expect_false(identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b))))
  # ... while the axis and category text is untouched, so nothing scientific
  # depends on which seed was chosen
  testthat::expect_identical(svg_labels(a), svg_labels(b))
})

testthat::test_that("geom_jitter must delegate to position_jitter to be seeded", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  # geom_jitter has no seed argument. Passing seed = to it lands in ..., warns
  # only "Ignoring unknown parameters", and is not honoured - which looks like a
  # fix and behaves like the defect. The repository therefore routes every
  # geom_jitter through an explicit position_jitter; this records why.
  testthat::expect_false("seed" %in% names(formals(ggplot2::geom_jitter)))
  testthat::expect_true("seed" %in% names(formals(ggplot2::position_jitter)))
  testthat::expect_true("seed" %in% names(formals(ggplot2::position_jitterdodge)))
  testthat::expect_true(is.na(formals(ggplot2::position_jitter)$seed))

  tmp <- withr::local_tempdir()
  d <- data.frame(g = rep(c("A", "B"), each = 10L), y = seq_len(20L) / 20)
  # the warning is raised when the layer is CONSTRUCTED, not when it is drawn,
  # so it is easy to lose sight of in a script that builds plots far from where
  # it saves them
  testthat::expect_warning(
    layer <- ggplot2::geom_jitter(width = 0.15, seed = 1L),
    "[Ii]gnoring unknown parameters")
  bad <- ggplot2::ggplot(d, ggplot2::aes(g, y)) + layer
  # and it really is not honoured: still nondeterministic despite the seed
  set.seed(11)
  a <- suppressWarnings(render_to(bad, file.path(tmp, "a.svg")))
  set.seed(99)
  b <- suppressWarnings(render_to(bad, file.path(tmp, "b.svg")))
  testthat::expect_false(identical(unname(tools::sha256sum(a)),
                                   unname(tools::sha256sum(b))))
})

testthat::test_that("the geom_jitter restructuring did not change height jitter", {
  testthat::skip_if_not(has_pkgs(), "ggplot2/svglite unavailable")
  # assess_joint_compartment_quality.R was the one call that had to be
  # restructured: geom_jitter(width = .15) became
  # geom_jitter(position = position_jitter(width = .15, seed = ...)), because
  # ggplot2 errors when position and width are both supplied. That restructuring
  # is only safe if an unsupplied height means the same thing in both forms.
  testthat::expect_null(formals(ggplot2::geom_jitter)$height)
  testthat::expect_null(formals(ggplot2::position_jitter)$height)

  tmp <- withr::local_tempdir()
  d <- data.frame(g = rep(c("A", "B", "C"), each = 20L), y = round(seq_len(60L) / 60, 3))
  base <- ggplot2::ggplot(d, ggplot2::aes(g, y))
  set.seed(101)
  bare <- render_to(base + ggplot2::geom_jitter(width = 0.15), file.path(tmp, "bare.svg"))
  set.seed(101)
  explicit <- render_to(
    base + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15)),
    file.path(tmp, "explicit.svg"))
  # byte-identical: the two forms consume the same draws and apply the same
  # default height, so the restructuring is inert apart from the seed source
  testthat::expect_identical(unname(tools::sha256sum(bare)),
                             unname(tools::sha256sum(explicit)))

  # and the reason the restructuring was unavoidable
  testthat::expect_error(
    print(base + ggplot2::geom_jitter(
      width = 0.15, position = ggplot2::position_jitter(width = 0.15, seed = 1L))),
    "Both .position. and")
})

testthat::test_that("the accepted figure package and its freeze are untouched", {
  # Phase 6H.7 changes future renders only. Nothing in the accepted package was
  # rerendered, so its manifest must still hash exactly as Phase 6H.5D left it.
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure export manifest not present")
  m <- utils::read.csv(mp, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(m), 5582L)
  testthat::expect_identical(
    unname(tools::sha256sum(mp)),
    "0fd0c9ed9febbc05ed6928b7fc6bfffdab845daba6255d8b47c607dee5d3c02c")
})
