source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "export_helpers.R"))

repo <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
exporter <- file.path(repo, "09_export_pride_journal", "08_export_manuscript_figures.R")

# A target directory the same depth as the real one on the analysis machine
# (results/manuscript/extended_data under a deep share path was 103 chars).
deep_dir <- file.path(
  tempdir(),
  paste0("d", strrep("x", 40)),
  "results", "manuscript", "extended_data"
)

# The exact source-relative path whose flattened target measured 276 characters
# and silently failed to copy.
long_rel <- paste(
  "03_qc_exploration",
  "04e_control_compartment_abundance_publication_figures",
  "global", "rendering_review",
  "pre_nature_minimal_control_rank_abundance_extended_data_183mm.svg",
  sep = "/"
)

testthat::test_that("previously over-length targets now fit the path budget", {
  budget <- manuscript_figure_path_budget()
  testthat::expect_lt(budget, 260L)   # must stay under Windows MAX_PATH

  # Old scheme: flatten the whole source-relative stem, no budget.
  old_target <- file.path(
    deep_dir, paste0(safe_filename(tools::file_path_sans_ext(long_rel)), ".svg")
  )
  testthat::expect_gt(nchar(old_target), 259L)   # reproduces the failure mode

  new_target <- manuscript_figure_target_paths(long_rel, deep_dir)
  testthat::expect_length(new_target, 1L)
  testthat::expect_lte(nchar(new_target), budget)
  testthat::expect_identical(tolower(tools::file_ext(new_target)), "svg")
})

testthat::test_that("every target stays within budget across many deep inputs", {
  rels <- c(
    long_rel,
    paste0("06_modules_WGCNA/", strrep("nested_directory_segment/", 8), "plot.pdf"),
    paste0("04_differential_expression_enrichment/", strrep("a", 220), ".png"),
    "07_spatial_networks/short.svg"
  )
  out <- manuscript_figure_target_paths(rels, deep_dir)
  testthat::expect_length(out, length(rels))
  testthat::expect_true(all(nchar(out) <= manuscript_figure_path_budget()))
  testthat::expect_identical(anyDuplicated(out), 0L)
})

testthat::test_that("long paths sharing a truncated prefix stay distinct", {
  shared <- paste0("06_modules_WGCNA/", strrep("same_long_prefix_segment/", 7))
  a <- paste0(shared, "variant_alpha.svg")
  b <- paste0(shared, "variant_beta.svg")
  # Distinguishing characters sit beyond the truncation point.
  c2 <- paste0(shared, strrep("z", 60), "_tail_one.svg")
  d <- paste0(shared, strrep("z", 60), "_tail_two.svg")

  out <- manuscript_figure_target_paths(c(a, b, c2, d), deep_dir)
  testthat::expect_identical(anyDuplicated(out), 0L)
  testthat::expect_length(unique(out), 4L)
  testthat::expect_true(all(nchar(out) <= manuscript_figure_path_budget()))

  # The hash is derived from the COMPLETE source-relative path, so it differs
  # even when the retained prefix is identical.
  testthat::expect_false(identical(
    manuscript_figure_target_hash(c2), manuscript_figure_target_hash(d)
  ))
  # ...and it is deterministic.
  testthat::expect_identical(
    manuscript_figure_target_hash(c2), manuscript_figure_target_hash(c2)
  )
  testthat::expect_identical(
    manuscript_figure_target_paths(c2, deep_dir),
    manuscript_figure_target_paths(c2, deep_dir)
  )
})

testthat::test_that("an impossibly deep target directory fails closed", {
  absurd <- file.path(tempdir(), strrep("q", manuscript_figure_path_budget()))
  testthat::expect_error(
    manuscript_figure_target_paths("a/b/c.svg", absurd),
    "path budget"
  )
})

testthat::test_that("duplicate targets are rejected before any copy", {
  t1 <- file.path(tempdir(), "same.svg")
  testthat::expect_error(
    assert_unique_export_targets(c(t1, t1), c("src_a.svg", "src_b.svg")),
    "collision"
  )
  testthat::expect_true(assert_unique_export_targets(
    c(file.path(tempdir(), "a.svg"), file.path(tempdir(), "b.svg"))
  ))
})

testthat::test_that("a failed copy raises a hard error instead of returning quietly", {
  root <- withr::local_tempdir("copy_fail_")
  good_src <- file.path(root, "good.svg")
  writeLines("x", good_src)
  missing_src <- file.path(root, "does_not_exist.svg")   # forces file.copy FALSE
  targets <- file.path(root, "out", c("good.svg", "bad.svg"))
  dir.create(dirname(targets[[1]]), recursive = TRUE)

  testthat::expect_error(
    copy_export_targets(c(good_src, missing_src), targets),
    "Export copy failed"
  )
  # Length mismatch is also an error, not silent recycling.
  testthat::expect_error(
    copy_export_targets(good_src, targets),
    "differ in length"
  )
  # The happy path returns invisibly TRUE and the target really exists.
  testthat::expect_true(copy_export_targets(good_src, targets[[1]]))
  testthat::expect_true(file.exists(targets[[1]]))
  testthat::expect_true(copy_export_targets(character(0), character(0)))
})

testthat::test_that("the exporter cannot write a manifest after a partial copy", {
  src <- readLines(exporter, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]

  copy_at <- grep("copy_export_targets(", code, fixed = TRUE)
  unique_at <- grep("assert_unique_export_targets(", code, fixed = TRUE)
  manifest_at <- grep("write.csv(manifest", code, fixed = TRUE)
  runmanifest_at <- grep("write_run_manifest(", code, fixed = TRUE)

  testthat::expect_length(copy_at, 1L)
  testthat::expect_length(unique_at, 1L)
  testthat::expect_gte(length(manifest_at), 1L)
  testthat::expect_length(runmanifest_at, 1L)

  # Uniqueness assertion, then the verified copy, then and only then the
  # manifests. copy_export_targets() stops on any failure, so reaching the
  # manifest writes proves every file landed.
  testthat::expect_lt(unique_at[[1]], copy_at[[1]])
  testthat::expect_lt(copy_at[[1]], min(manifest_at))
  testthat::expect_lt(copy_at[[1]], runmanifest_at[[1]])

  # The raw unchecked copy must be gone.
  testthat::expect_false(any(grepl("file.copy(manifest$source_file", code, fixed = TRUE)))
  # And the budget is enforced explicitly before copying.
  testthat::expect_true(any(grepl("manuscript_figure_path_budget()", code, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# B6: the generic relative_to() must not be shadowed by the export chain.
# ---------------------------------------------------------------------------
testthat::test_that("sourcing pride_helpers cannot shadow the generic relative_to", {
  # pride_helpers.R is pulled in by export_helpers.R. Source it into the same
  # environment as paths.R -- exactly what the export scripts do -- and require
  # that the generic relative_to() survives with its default root intact. It
  # previously did not, which broke write_run_manifest().
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo, "R", "paths.R"), envir = env)
  generic <- get("relative_to", envir = env)
  testthat::expect_false(identical(formals(generic)$root, quote(expr = )))

  sys.source(file.path(repo, "R", "pride_helpers.R"), envir = env)
  after <- get("relative_to", envir = env)

  testthat::expect_identical(after, generic)
  testthat::expect_false(identical(formals(after)$root, quote(expr = )))
  testthat::expect_true(exists("pride_relative_to", envir = env, inherits = FALSE))
  # The PRIDE variant still requires an explicit root.
  testthat::expect_true(identical(
    formals(get("pride_relative_to", envir = env))$root, quote(expr = )
  ))

  # The one-argument form write_run_manifest() relies on must work.
  testthat::expect_silent(after(file.path(repo, "README.md")))
})

testthat::test_that("write_run_manifest survives the full export source chain", {
  # Regression for the crash: relative_to(session_path) inside
  # write_run_manifest() died with 'argument "root" is missing'.
  script <- c(
    sprintf('setwd(%s)', shQuote(repo)),
    'source("R/paths.R")',
    'source(file.path("R", "export_helpers.R"))',
    'out <- file.path(tempdir(), "wrm_probe", "run_manifest.yml")',
    'write_run_manifest(out, inputs = list(), outputs = list(), notes = "probe")',
    'cat("WROTE:", file.exists(out), "\n")'
  )
  f <- withr::local_tempfile(fileext = ".R")
  writeLines(script, f)
  out <- suppressWarnings(system2("Rscript", c(f), stdout = TRUE, stderr = TRUE))
  testthat::expect_true(any(grepl("WROTE: TRUE", out, fixed = TRUE)),
                        info = paste(utils::tail(out, 6), collapse = " | "))
})

testthat::test_that("only R/paths.R defines relative_to", {
  helpers <- list.files(file.path(repo, "R"), pattern = "[.]R$", full.names = TRUE)
  definers <- helpers[vapply(helpers, function(f) {
    any(grepl("^relative_to <- function", readLines(f, warn = FALSE)))
  }, logical(1))]
  testthat::expect_identical(basename(definers), "paths.R")
})

# ---------------------------------------------------------------------------
# B7: the figure export must cover the manuscript-panel roots.
# ---------------------------------------------------------------------------
testthat::test_that("figure export declares the manuscript-panel figure roots", {
  joined <- paste(readLines(exporter, warn = FALSE), collapse = "\n")
  for (root in c("manuscript_panels", "10_biological_integration",
                 "08_biological_interpretation")) {
    testthat::expect_match(joined, paste0('path_results("figures", "', root, '")'),
                           fixed = TRUE)
  }
  # Still an explicit list, never an unrestricted recursive results/figures scan.
  testthat::expect_match(joined, "canonical_ewce_figure_root()", fixed = TRUE)
  testthat::expect_false(grepl('list.files(path_results("figures")', joined, fixed = TRUE))
})

testthat::test_that("all six current Figure-3 panels are export-eligible", {
  panel_dir <- file.path(repo, "results", "figures", "manuscript_panels", "figure_3")
  testthat::skip_if_not(dir.exists(panel_dir), "Figure-3 panel figures not in this checkout")
  panels <- list.files(panel_dir, pattern = "[.](svg|pdf|png)$", full.names = TRUE)
  testthat::expect_length(panels, 6L)
  # Neither selection filter may remove a current Figure-3 panel.
  testthat::expect_true(all(is_exportable_result_path(panels)))
  testthat::expect_identical(drop_legacy_dataset_suffixed_aliases(panels), panels)
  testthat::expect_identical(drop_orphan_figure_families(panels), panels)
})
