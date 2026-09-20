# Phase 6G.8 section 27: the failed WGCNA run must never be selectable by a
# publication collector.
#
# microglia_failed_20260720_133211 is kept on disk for provenance and has the
# full shape of a real run - same subdirectories, same figure stems, same file
# types. The manuscript figure exporter walks the historical stage recursively,
# and none of its other filters inspect the dataset scope, so before this guard
# existed 50 of the failed run's figures were selected and copied into
# extended_data under their own name.
#
# The prefix trap is the point of these tests: "microglia" is a prefix of
# "microglia_failed_20260720_133211", so any startsWith/grepl-style dataset
# test admits exactly the directory it is meant to exclude. Membership must be
# an exact segment comparison against valid_datasets().

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "export_helpers.R"))

test_that("the failed run is rejected and its canonical sibling is kept", {
  keep <- "results/figures/06_modules_WGCNA/01_WGCNA/microglia/main/x.svg"
  drop <- "results/figures/06_modules_WGCNA/01_WGCNA/microglia_failed_20260720_133211/main/x.svg"
  expect_false(is_noncanonical_wgcna_dataset_scope(keep))
  expect_true(is_noncanonical_wgcna_dataset_scope(drop))
  expect_equal(drop_noncanonical_wgcna_dataset_scopes(c(keep, drop)), keep)
})

test_that("every canonical dataset survives the guard", {
  for (ds in valid_datasets()) {
    p <- file.path("results/figures/06_modules_WGCNA/01_WGCNA", ds, "main/x.svg")
    expect_false(is_noncanonical_wgcna_dataset_scope(p), info = ds)
  }
})

test_that("a suffixed variant of any canonical dataset is rejected", {
  for (ds in valid_datasets()) {
    p <- file.path("results/figures/06_modules_WGCNA/01_WGCNA",
                   paste0(ds, "_failed_20260720_133211"), "main/x.svg")
    expect_true(is_noncanonical_wgcna_dataset_scope(p), info = ds)
  }
})

test_that("cross-dataset aggregate scopes are not dataset positions and are kept", {
  for (scope in c("all", "global")) {
    p <- file.path("results/figures/06_modules_WGCNA/interpretable_summary", scope, "x.svg")
    expect_false(is_noncanonical_wgcna_dataset_scope(p), info = scope)
  }
})

test_that("a content directory in the scope position is not mistaken for a dataset", {
  # families whose third level is content, not scope: pruning these would
  # silently drop real figures from the export
  for (seg in c("modules", "supermodules", "main", "qc", "traits")) {
    p <- file.path("results/figures/06_modules_WGCNA/01_WGCNA", seg, "x.svg")
    expect_false(is_noncanonical_wgcna_dataset_scope(p), info = seg)
  }
})

test_that("paths outside any WGCNA tree are never judged", {
  expect_false(is_noncanonical_wgcna_dataset_scope(
    "results/figures/07_spatial_networks/microglia_failed_x/y.svg"))
  # the word appears, but not in a WGCNA dataset-scope position
  expect_false(is_noncanonical_wgcna_dataset_scope(
    "results/qc/summarize_microglia_markers/global/plots/microglia_overview.svg"))
  expect_false(is_noncanonical_wgcna_dataset_scope(
    "docs/notes_about_microglia_failed_runs.md"))
})

# Section 14: the guard must cover BOTH discovery routes. WGCNA is mid
# migration, so a rerun lands in results/wgcna/... where there is no
# 06_modules_WGCNA segment to key on -- and a rerun is exactly when a new
# failed run appears.
test_that("normalized WGCNA scopes are judged too", {
  expect_false(is_noncanonical_wgcna_dataset_scope(
    "results/wgcna/build_wgcna_modules/microglia/plots/x.svg"))
  expect_true(is_noncanonical_wgcna_dataset_scope(
    "results/wgcna/build_wgcna_modules/microglia_failed_20260720_133211/plots/x.svg"))
  expect_true(is_noncanonical_wgcna_dataset_scope(
    "results/wgcna/build_wgcna_modules/neuron_soma_failed_20270101_000000/plots/x.svg"))
})

test_that("all four section-14 route/scope combinations behave", {
  cases <- list(
    list(p = "results/figures/06_modules_WGCNA/01_WGCNA/microglia/main/x.svg",
         reject = FALSE, what = "canonical historical"),
    list(p = "results/figures/06_modules_WGCNA/01_WGCNA/microglia_failed_20260720_133211/main/x.svg",
         reject = TRUE, what = "failed historical"),
    list(p = "results/wgcna/build_wgcna_modules/neuron_neuropil/plots/x.svg",
         reject = FALSE, what = "canonical normalized"),
    list(p = "results/wgcna/build_wgcna_modules/microglia_failed_20260720_133211/plots/x.svg",
         reject = TRUE, what = "failed normalized")
  )
  for (c1 in cases) {
    expect_equal(is_noncanonical_wgcna_dataset_scope(c1$p), c1$reject, info = c1$what)
  }
})

test_that("normalized non-dataset scopes are kept", {
  for (scope in c("global", "all")) {
    expect_false(is_noncanonical_wgcna_dataset_scope(
      file.path("results/wgcna/summarize_module_interpretation", scope, "plots/x.svg")),
      info = scope)
  }
})

test_that("the scope extractor reads the documented position on each route", {
  hist <- strsplit("results/figures/06_modules_WGCNA/group_effects/neuron_soma/x.svg", "/")[[1]]
  norm <- strsplit("results/wgcna/test_module_phenotypes/neuron_soma/plots/x.svg", "/")[[1]]
  expect_equal(wgcna_dataset_scope_of(hist), "neuron_soma")
  expect_equal(wgcna_dataset_scope_of(norm), "neuron_soma")
  expect_true(is.na(wgcna_dataset_scope_of(strsplit("results/qc/x/global/y.svg", "/")[[1]])))
})

test_that("the failed run is never eligible for publication export", {
  # FAILED_RUN_PROVENANCE is the fourth state carrier: retained on disk,
  # never a publication candidate
  failed <- "microglia_failed_20260720_133211"
  for (route in c(file.path("results/figures/06_modules_WGCNA/01_WGCNA", failed, "main/x.svg"),
                  file.path("results/wgcna/build_wgcna_modules", failed, "plots/x.svg"),
                  file.path("results/source_data/06_modules_WGCNA/01_WGCNA", failed, "x.csv"))) {
    expect_true(is_noncanonical_wgcna_dataset_scope(route), info = route)
    expect_length(drop_noncanonical_wgcna_dataset_scopes(route), 0L)
  }
})

test_that("the guard is a no-op on an empty or already-clean selection", {
  expect_length(drop_noncanonical_wgcna_dataset_scopes(character(0)), 0L)
  clean <- c("results/figures/06_modules_WGCNA/01_WGCNA/microglia/main/a.svg",
             "results/figures/06_modules_WGCNA/group_effects/neuron_soma/b.pdf")
  expect_equal(drop_noncanonical_wgcna_dataset_scopes(clean), clean)
})

test_that("the real figure tree yields no failed-run candidate after the guard", {
  root <- path_results("figures", "06_modules_WGCNA")
  skip_if_not(dir.exists(root), "historical WGCNA figure tree not present")
  cand <- list.files(root, pattern = "[.](svg|pdf|png)$", recursive = TRUE,
                     full.names = TRUE)
  skip_if(length(cand) == 0, "no figures on disk")
  kept <- drop_noncanonical_wgcna_dataset_scopes(cand)
  expect_equal(sum(grepl("_failed_", kept)), 0L)
  # and the guard must remove ONLY failed-run files, never anything else
  removed <- setdiff(cand, kept)
  expect_true(all(grepl("_failed_", removed)))
})

# ------------------------------------------- payload integrity evidence rule
#
# Phase 6G.8 section 3. Three batches reported "results/ git-clean" and
# "exports/ git-clean" as evidence that those trees were unchanged. They are
# not evidence: .gitignore excludes /results/ and /exports/**, so `git status`
# over them can only ever see the handful of force-added placeholder files. A
# whole export could be rewritten and git would still report nothing.
#
# Payload integrity must be established by direct recursive inventory, file
# counts, byte counts, SHA-256 and manifest/source comparison. These tests pin
# the gitignore fact so the mistaken inference cannot be made again silently.

test_that("results/ and exports/ are gitignored, so git status is not payload evidence", {
  skip_if_not(dir.exists(testthat::test_path("..", "..", ".git")), "not a git checkout")
  gitignore <- readLines(testthat::test_path("..", "..", ".gitignore"), warn = FALSE)
  expect_true(any(grepl("^/results/?$", trimws(gitignore))),
              info = "/results/ must be gitignored for this test to be meaningful")
  expect_true(any(grepl("^/exports/", trimws(gitignore))),
              info = "/exports/ must be gitignored for this test to be meaningful")
})

test_that("the tracked footprint of results/ and exports/ is placeholders only", {
  skip_if_not(dir.exists(testthat::test_path("..", "..", ".git")), "not a git checkout")
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  tracked <- suppressWarnings(system2("git", c("-C", shQuote(root), "ls-files", "results", "exports"),
                                      stdout = TRUE, stderr = FALSE))
  skip_if(length(tracked) == 0, "git unavailable")
  # every tracked entry is a placeholder or a small committed contract file,
  # never bulk payload: if this ever grows, the evidence rule above needs review
  expect_lt(length(tracked), 25)
  expect_true(all(grepl("[.]gitkeep$|[.]csv$|[.]tsv$|[.]md$", tracked)))
})
