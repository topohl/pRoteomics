# The budgeted figure-name contract.
#
# Three ggsave() sites in analysis/wgcna/score_module_activity.R built their
# target with paste0() and no path budget. The Windows file API then truncated
# the path at MAX_PATH-1 rather than failing, so 354 figures landed on disk at
# exactly 258 or 259 characters with ".svg" cut off - and the manuscript
# exporter selects on \.(svg|pdf|png)$, so every one of them silently dropped
# out of publication discovery. A filename bug was really a selection bug.
#
# These tests pin: the extension is atomic, the budget is taken from the
# absolute path rather than a fixed basename limit, an impossible budget fails
# loudly, and none of the observed damage shapes can recur.

source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("the two budgets stay distinct", {
  testthat::expect_identical(FIGURE_WRITE_BUDGET, 240L)
  testthat::expect_identical(PATH_LENGTH_WALL, 260L)
  testthat::expect_lt(FIGURE_WRITE_BUDGET, PATH_LENGTH_WALL)
})

testthat::test_that("the filename budget comes from the directory, not a constant", {
  shallow <- "S:/x"
  deep <- file.path("S:/x", strrep("d", 180L))
  testthat::expect_gt(figure_filename_budget(shallow), figure_filename_budget(deep))
  testthat::expect_identical(
    figure_filename_budget(shallow),
    FIGURE_WRITE_BUDGET - path_length_chars(shallow) - 1L)
})

testthat::test_that("a name already within budget is untouched", {
  d <- "S:/figs"
  testthat::expect_identical(budgeted_figure_target(d, "module_score_CA1.svg"),
                             "module_score_CA1.svg")
})

testthat::test_that("the budget is spent on the stem and the extension survives whole", {
  d <- file.path("S:/r", strrep("d", 150L))
  out <- budgeted_figure_target(d, paste0("module_score_", strrep("m", 200L), ".svg"))
  testthat::expect_identical(tools::file_ext(out), "svg")
  testthat::expect_lte(path_length_chars(file.path(d, out)), FIGURE_WRITE_BUDGET)
  testthat::expect_gt(nchar(tools::file_path_sans_ext(out)), 0L)
})

testthat::test_that("none of the observed damage shapes can be produced", {
  # .s / .sv / .p / .pn / extensionless / trailing dot were all seen on disk
  for (ext in c("svg", "pdf", "png")) {
    for (dirlen in c(20L, 80L, 150L, 190L, 200L)) {
      d <- file.path("S:/r", strrep("d", dirlen))
      out <- budgeted_figure_target(d, paste0(strrep("s", 120L), ".", ext))
      testthat::expect_identical(tools::file_ext(out), ext,
                                 info = paste(ext, dirlen))
      testthat::expect_false(endsWith(out, "."), info = paste(ext, dirlen))
      testthat::expect_true(nzchar(tools::file_ext(out)))
      for (bad in c(".s", ".sv", ".p", ".pn", ".pd"))
        testthat::expect_false(endsWith(out, bad), info = paste(ext, dirlen, bad))
    }
  }
})

testthat::test_that("an impossible budget fails loudly instead of clipping", {
  # a directory that leaves no room for a stem character plus ".svg"
  d <- file.path("S:/r", strrep("d", FIGURE_WRITE_BUDGET - 8L))
  testthat::expect_lt(figure_filename_budget(d), nchar(".svg") + 1L)
  testthat::expect_error(budgeted_figure_target(d, "plot.svg"),
                         "cannot hold a stem character")
})

testthat::test_that("exact budget, one over and two over behave correctly", {
  d <- "S:/figs"
  room <- figure_filename_budget(d)
  exact <- paste0(strrep("a", room - 4L), ".svg")
  testthat::expect_identical(nchar(exact), room)
  testthat::expect_identical(budgeted_figure_target(d, exact), exact)

  over1 <- paste0(strrep("a", room - 3L), ".svg")
  o1 <- budgeted_figure_target(d, over1)
  testthat::expect_identical(nchar(o1), room)
  testthat::expect_identical(tools::file_ext(o1), "svg")

  over2 <- paste0(strrep("a", room - 2L), ".svg")
  o2 <- budgeted_figure_target(d, over2)
  testthat::expect_identical(nchar(o2), room)
  testthat::expect_identical(tools::file_ext(o2), "svg")
})

testthat::test_that("collisions are disambiguated deterministically", {
  d <- file.path("S:/r", strrep("d", 170L))
  long <- paste0("module_score_CA1_slm_Neuropil_", strrep("x", 90L), ".svg")
  first <- budgeted_figure_target(d, long)
  second <- budgeted_figure_target(d, long, taken = first)
  testthat::expect_false(identical(first, second))
  testthat::expect_identical(tools::file_ext(second), "svg")
  testthat::expect_lte(path_length_chars(file.path(d, second)), FIGURE_WRITE_BUDGET)
  testthat::expect_match(second, "__[0-9a-f]{8}[.]svg$")
  # deterministic on repeat
  testthat::expect_identical(second, budgeted_figure_target(d, long, taken = first))
})

testthat::test_that("a stem with several dots keeps only the final extension", {
  d <- file.path("S:/r", strrep("d", 150L))
  out <- budgeted_figure_target(d, paste0("a.b.c.", strrep("d", 120L), ".svg"))
  ## Interior dots belong to the stem and may survive; what matters is that the
  ## FINAL extension is the declared one and the name does not end in a dot.
  testthat::expect_identical(tools::file_ext(out), "svg")
  testthat::expect_match(out, "[.]svg$")
  testthat::expect_false(endsWith(out, "."))
  testthat::expect_true(startsWith(out, "a.b.c."))
})

testthat::test_that("an extensionless figure name is budgeted without inventing one", {
  d <- file.path("S:/r", strrep("d", 150L))
  out <- budgeted_figure_target(d, strrep("q", 200L))
  testthat::expect_identical(tools::file_ext(out), "")
  testthat::expect_lte(path_length_chars(file.path(d, out)), FIGURE_WRITE_BUDGET)
})

testthat::test_that("the real successor directory has room for the real worst-case name", {
  # The canonical root the migrated writer now uses, at its deepest observed
  # shape: results/wgcna/score_module_activity/<dataset>/plots/<group>/<sub>/<label>
  d <- canonical_result_path("wgcna", "score_module_activity", "neuron_neuropil",
                             "plots", "overlap", "module_group_scores",
                             "sensitivity_flagged_replicates_removed")
  room <- figure_filename_budget(d)
  testthat::expect_gt(room, 40L)
  # the longest real name the writer emits
  worst <- paste0("module_score_CA1_slm_Neuropil_chromatin_RNP_related_exploratory_",
                  "sensitivity_flagged_replicates_removed.svg")
  out <- budgeted_figure_target(d, worst)
  testthat::expect_identical(tools::file_ext(out), "svg")
  testthat::expect_lte(path_length_chars(file.path(d, out)), FIGURE_WRITE_BUDGET)
  # and the correlation template, which is the longest of the three
  worst_cor <- paste0("cor_CA1_slm_Neuropil_chromatin_RNP_related_exploratory_",
                      "AUC_norm_firstActive_primary_all_replicates.svg")
  out2 <- budgeted_figure_target(d, worst_cor)
  testthat::expect_identical(tools::file_ext(out2), "svg")
  testthat::expect_lte(path_length_chars(file.path(d, out2)), FIGURE_WRITE_BUDGET)
})

testthat::test_that("the three write sites route through the budgeted helper", {
  src <- readLines(testthat::test_path("..", "..", "analysis", "wgcna",
                                       "score_module_activity.R"), warn = FALSE)
  # no ggsave target may still be built by a bare file.path(paste0(...)) with a
  # hard-coded extension
  unsafe <- grep("ggsave\\(\\s*file\\.path\\(.*paste0\\(", src)
  testthat::expect_identical(length(unsafe), 0L)
  testthat::expect_gte(length(grep("budgeted_figure_path\\(", src)), 3L)
})

testthat::test_that("no figure under results/figures carries a damaged extension", {
  root <- repo_path("results", "figures")
  testthat::skip_if_not(dir.exists(root), "figure tree not present")
  f <- list.files(root, recursive = TRUE)
  testthat::skip_if(!length(f), "figure tree empty")
  ## legitimately-odd names that are not figures
  f <- f[!basename(f) %in% c("Thumbs.db", ".gitkeep")]
  f <- f[!grepl("[.][Rr]$", f)]
  ## the frozen archive trees are preserved evidence and out of scope
  f <- f[!grepl("^manuscript/_(superseded|failed)", f)]
  good <- c("svg", "pdf", "png", "csv", "tsv", "txt", "xlsx", "rds", "yml", "yaml",
            "json", "md", "html", "log", "flag", "gct", "zip", "gz", "xls")
  bad <- f[!tolower(tools::file_ext(f)) %in% good]
  ## Five truncated duplicates under microglia_validation_proposed/ are held
  ## deliberately: they duplicate already-exported canonical figures byte for
  ## byte, so repairing their names would add duplicate publication selections.
  held <- grepl("microglia_validation_proposed/", bad, fixed = TRUE)
  testthat::expect_identical(sum(!held), 0L,
    info = paste("unexpected damaged figure names:",
                 paste(utils::head(bad[!held], 5), collapse = ", ")))
  testthat::expect_lte(sum(held), 5L)
})

testthat::test_that("every repaired figure is discoverable and byte-identical", {
  a <- testthat::test_path("..", "..", "audits", "phase6h_figure_filename_repair_after.csv")
  testthat::skip_if_not(file.exists(a), "repair audit not present")
  d <- utils::read.csv(a, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 354L)
  testthat::expect_true(all(d$same_hash))
  testthat::expect_true(all(d$same_size))
  testthat::expect_true(all(d$mtime_preserved))
  # the exporter's selection regex now matches every one of them
  testthat::expect_true(all(grepl("[.](svg|pdf|png)$", d$new_path)))
  testthat::expect_false(any(grepl("[.](svg|pdf|png)$", basename(d$old_path))))
  # and none is under a failed-run scope
  testthat::expect_false(any(grepl("_failed_", d$new_path)))
  # every repaired file is readable by R, i.e. under the wall
  testthat::expect_true(all(d$new_abs_chars < PATH_LENGTH_WALL))
})
