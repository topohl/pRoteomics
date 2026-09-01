source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "export_helpers.R"))

repo <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
exporter <- file.path(repo, "09_export_pride_journal", "09_export_source_data.R")

# Output directories as deep as the real ones (results/manuscript/source_data
# under a deep share path measured 102 characters).
deep <- file.path(tempdir(), paste0("d", strrep("x", 40)), "results", "manuscript")
sd_dir <- file.path(deep, "source_data")
st_dir <- file.path(deep, "supplementary_tables")

abs_src <- function(rel_full) file.path(tempdir(), "repo", rel_full)

testthat::test_that("targets route by root and stay within the 240-character budget", {
  rels <- c(
    "results/source_data/04_differential_expression_enrichment/control_spatial_identity_validation/global/figure2f_source_data.csv",
    "results/tables/06_modules_WGCNA/group_effects/microglia/module_group_effects.csv",
    paste0("results/source_data/04_differential_expression_enrichment/",
           strrep("deeply_nested_segment/", 9), "table.csv"),
    paste0("results/tables/04_differential_expression_enrichment/", strrep("b", 210), ".tsv")
  )
  out <- manuscript_table_target_paths(abs_src(rels), sd_dir, st_dir)

  testthat::expect_length(out, length(rels))
  testthat::expect_true(all(nchar(out) <= manuscript_figure_path_budget()))
  # Routing follows the repo-relative root, not a substring search. Compare on
  # normalised directories: the builder normalises separators.
  sd_norm <- normalize_export_path(sd_dir)
  st_norm <- normalize_export_path(st_dir)
  testthat::expect_identical(
    dirname(out), c(sd_norm, st_norm, sd_norm, st_norm)
  )
  testthat::expect_identical(anyDuplicated(out), 0L)
})

testthat::test_that("the real extension is preserved, including xlsx and tsv", {
  rels <- c(
    "results/tables/a/b/short.csv",
    "results/tables/a/b/short.tsv",
    "results/tables/a/b/short.xlsx",
    paste0("results/tables/a/b/", strrep("long", 60), ".xlsx")
  )
  out <- manuscript_table_target_paths(abs_src(rels), sd_dir, st_dir)
  testthat::expect_identical(tolower(tools::file_ext(out)),
                             c("csv", "tsv", "xlsx", "xlsx"))
  testthat::expect_true(all(nchar(out) <= manuscript_figure_path_budget()))
})

testthat::test_that("shortening is deterministic", {
  rel <- paste0("results/tables/04_differential_expression_enrichment/",
                strrep("segment_name/", 14), "deep_table.csv")
  a <- manuscript_table_target_paths(abs_src(rel), sd_dir, st_dir)
  b <- manuscript_table_target_paths(abs_src(rel), sd_dir, st_dir)
  testthat::expect_identical(a, b)
  testthat::expect_match(basename(a), "_[0-9a-f]{10}\\.csv$")
})

testthat::test_that("long paths sharing a truncated prefix remain distinct", {
  shared <- paste0("results/tables/04_differential_expression_enrichment/",
                   strrep("same_prefix_segment/", 8))
  rels <- c(
    paste0(shared, strrep("z", 70), "_alpha.csv"),
    paste0(shared, strrep("z", 70), "_beta.csv"),
    paste0(shared, strrep("z", 70), "_gamma.csv")
  )
  out <- manuscript_table_target_paths(abs_src(rels), sd_dir, st_dir)
  testthat::expect_identical(anyDuplicated(out), 0L)
  testthat::expect_length(unique(out), 3L)
  testthat::expect_true(all(nchar(out) <= manuscript_figure_path_budget()))
})

testthat::test_that("identical within-root relative paths in different roots stay distinct", {
  # The two roots go to different output directories, and the hash identity is
  # root-qualified, so neither the directory nor the name can collide.
  rel_tail <- paste0("04_differential_expression_enrichment/",
                     strrep("shared_segment/", 9), "same_name.csv")
  rels <- c(paste0("results/source_data/", rel_tail),
            paste0("results/tables/", rel_tail))
  out <- manuscript_table_target_paths(abs_src(rels), sd_dir, st_dir)
  testthat::expect_identical(anyDuplicated(out), 0L)
  testthat::expect_false(identical(dirname(out[[1]]), dirname(out[[2]])))
  testthat::expect_false(identical(basename(out[[1]]), basename(out[[2]])))
  # And the figure builder's default identity is unchanged.
  testthat::expect_identical(
    manuscript_figure_target_paths("a/b/c.svg", sd_dir),
    manuscript_figure_target_paths("a/b/c.svg", sd_dir, identity = "a/b/c.svg")
  )
})

testthat::test_that("an unroutable candidate fails closed", {
  testthat::expect_error(
    manuscript_table_target_paths(
      file.path(tempdir(), "repo", "results", "figures", "x", "y.csv"), sd_dir, st_dir
    ),
    "Cannot route"
  )
})

testthat::test_that("duplicate targets are rejected and the clashing sources named", {
  t1 <- file.path(sd_dir, "same.csv")
  srcs <- abs_src(c("results/source_data/a/same.csv", "results/tables/a/same.csv"))
  testthat::expect_error(
    assert_unique_table_export_targets(c(t1, t1), srcs),
    "target collision"
  )
  testthat::expect_true(assert_unique_table_export_targets(
    c(file.path(sd_dir, "a.csv"), file.path(sd_dir, "b.csv")),
    abs_src(c("results/source_data/a.csv", "results/source_data/b.csv"))
  ))
})

# ---------------------------------------------------------------------------
# Source accessibility
# ---------------------------------------------------------------------------
testthat::test_that("inaccessible sources fail closed with a character count", {
  root <- withr::local_tempdir("src_access_")
  good <- file.path(root, "good.csv")
  writeLines("x", good)
  missing <- file.path(root, "results", "source_data", "gone.csv")

  testthat::expect_true(assert_export_sources_accessible(good))
  testthat::expect_error(
    assert_export_sources_accessible(c(good, missing)),
    "inaccessible"
  )
  testthat::expect_true(assert_export_sources_accessible(character(0)))
})

testthat::test_that("every scoped source in this checkout is accessible", {
  tr <- c(path_results("tables"), path_results("source_data"))
  testthat::skip_if_not(any(dir.exists(tr)), "analysis results not in this checkout")
  cand <- unlist(lapply(tr[dir.exists(tr)], list.files,
                        pattern = "[.](csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE),
                 use.names = FALSE)
  cand <- cand[is_exportable_result_path(cand)]
  cand <- apply_manuscript_source_data_scope(cand)
  testthat::expect_gt(length(cand), 0L)
  testthat::expect_true(assert_export_sources_accessible(cand))
  # And every target fits the budget with no duplicates.
  tgt <- manuscript_table_target_paths(
    cand, path_results("manuscript", "source_data"),
    path_results("manuscript", "supplementary_tables")
  )
  testthat::expect_identical(sum(nchar(tgt) > manuscript_figure_path_budget()), 0L)
  testthat::expect_identical(anyDuplicated(tgt), 0L)
})

# ---------------------------------------------------------------------------
# Copy completeness
# ---------------------------------------------------------------------------
testthat::test_that("a failed copy is a hard error and a missing target is caught", {
  root <- withr::local_tempdir("sd_copy_")
  good <- file.path(root, "good.csv"); writeLines("x", good)
  absent <- file.path(root, "absent.csv")
  out <- file.path(root, "out"); dir.create(out)

  testthat::expect_error(
    copy_export_targets(c(good, absent), file.path(out, c("a.csv", "b.csv"))),
    "Export copy failed"
  )
  testthat::expect_true(copy_export_targets(good, file.path(out, "a.csv")))
  testthat::expect_true(file.exists(file.path(out, "a.csv")))
})

testthat::test_that("the exporter validates and copies before writing any manifest", {
  src <- readLines(exporter, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]

  guard <- grep("isTRUE(dry_run)", code, fixed = TRUE)
  access <- grep("assert_export_sources_accessible(", code, fixed = TRUE)
  uniq <- grep("assert_unique_table_export_targets(", code, fixed = TRUE)
  budget <- grep("manuscript_figure_path_budget()", code, fixed = TRUE)
  mk <- grep("dir_create(", code, fixed = TRUE)
  copy <- grep("copy_export_targets(", code, fixed = TRUE)
  manifest <- grep("write.csv(manifest", code, fixed = TRUE)
  runman <- grep("write_run_manifest(", code, fixed = TRUE)

  testthat::expect_length(guard, 1L)
  testthat::expect_length(access, 1L)
  testthat::expect_length(uniq, 1L)
  testthat::expect_length(copy, 1L)
  testthat::expect_length(runman, 1L)
  testthat::expect_gte(length(manifest), 1L)

  # Order: guard -> validate -> mkdir -> verified copy -> manifests.
  testthat::expect_lt(guard[[1]], access[[1]])
  testthat::expect_lt(access[[1]], copy[[1]])
  testthat::expect_lt(uniq[[1]], copy[[1]])
  testthat::expect_true(all(mk > guard[[1]]))
  testthat::expect_lt(copy[[1]], min(manifest))
  testthat::expect_lt(copy[[1]], runman[[1]])
  testthat::expect_true(any(budget < copy[[1]]))

  # The raw unchecked copy and the old truncating namer are gone.
  testthat::expect_false(any(grepl("file.copy(manifest$source_file", code, fixed = TRUE)))
  testthat::expect_false(any(grepl("table_target_name", code, fixed = TRUE)))
})

testthat::test_that("dry-run stays non-mutating and reports the required counters", {
  src <- readLines(exporter, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]
  guard <- grep("isTRUE(dry_run)", code, fixed = TRUE)[[1]]
  mutating <- grep(
    paste("dir_create", "dir\\.create", "file\\.copy", "write\\.csv", "write\\.table",
          "writeLines", "saveRDS", "write_run_manifest", "copy_export_targets", sep = "|"),
    code
  )
  testthat::expect_identical(mutating[mutating < guard], integer(0))

  joined <- paste(src, collapse = "\n")
  for (label in c("Raw candidates before journal scope", "Scoped selected files",
                  "Selected bytes", "Inaccessible selected sources",
                  "Duplicate target paths", "Longest target path")) {
    testthat::expect_match(joined, label, fixed = TRUE)
  }
  testthat::expect_match(joined, "Excluded (", fixed = TRUE)
  # The manuscript/results exclusion guard is still applied.
  testthat::expect_match(joined, "is_exportable_result_path(candidates)", fixed = TRUE)
})

testthat::test_that("PRIDE selectors reference none of the manuscript export helpers", {
  pride <- c("R/pride_helpers.R",
             "09_export_pride_journal/05_make_pride_manifest.R",
             "09_export_pride_journal/10_validate_pride_submission.R",
             "09_export_pride_journal/03_export_processed_pg_matrix_package.R",
             "09_export_pride_journal/04_make_supplementary_tables.R")
  banned <- c("manuscript_table_target_paths", "assert_unique_table_export_targets",
              "assert_export_sources_accessible", "source_data_scope_exclusion_reasons",
              "apply_manuscript_source_data_scope", "manuscript_figure_target_paths")
  for (p in pride) {
    joined <- paste(readLines(file.path(repo, p), warn = FALSE), collapse = "\n")
    for (b in banned) {
      testthat::expect_false(grepl(b, joined, fixed = TRUE), info = paste(p, b))
    }
  }
})
