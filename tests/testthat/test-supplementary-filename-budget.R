# The filename budget for staged journal supplementary tables.
#
# trim_with_ext() used to call safe_filename(x, max_chars = max_chars) first,
# and safe_filename() ends in substr(x, 1, max_chars). The name was therefore
# truncated before its extension was considered, which cut through ".csv" and
# left ".c", ".cs" or nothing; and because the result was then always within
# budget, the guard that was meant to reattach the extension could never run.
# 63 tables reached the PRIDE bundle with a clipped or absent extension.
#
# These tests pin the corrected contract: the extension is atomic, the budget
# is spent on the stem, and an impossible budget fails loudly rather than
# emitting a partial extension.

source(testthat::test_path("..", "..", "R", "paths.R"))
## supplementary_stage_names() calls dataset_for_manifest_file() and
## relative_to() from the export helpers, so the collision test needs them.
source(repo_path("R", "export_helpers.R"))

# The functions live inside an analysis script that does work at top level, so
# only their definitions are evaluated here - the repository's usual pattern
# for reaching a script-local function.
BUDGET_SCRIPT <- testthat::test_path("..", "..", "analysis", "publication_source_data",
                                     "build_supplementary_tables.R")
budget_env <- new.env(parent = globalenv())
local({
  wanted <- c("trim_with_ext", "stable_path_hash", "supplementary_stage_names")
  for (expr in parse(BUDGET_SCRIPT, keep.source = FALSE)) {
    if (!is.call(expr) || length(expr) < 3L) next
    if (!identical(as.character(expr[[1]]), "<-")) next
    target <- expr[[2]]
    if (!is.name(target) || !as.character(target) %in% wanted) next
    eval(expr, envir = budget_env)
  }
})
trim_with_ext <- get("trim_with_ext", envir = budget_env)

testthat::test_that("the definitions under test were actually extracted", {
  testthat::expect_true(is.function(trim_with_ext))
  testthat::expect_true(exists("stable_path_hash", envir = budget_env))
  testthat::expect_true(exists("supplementary_stage_names", envir = budget_env))
})

testthat::test_that("a name already within budget is returned untouched", {
  testthat::expect_identical(trim_with_ext("short_table.csv", 96L), "short_table.csv")
  testthat::expect_identical(trim_with_ext("a.csv", 5L), "a.csv")
})

testthat::test_that("the budget is spent on the stem and the extension survives whole", {
  stem <- strrep("x", 200L)
  out <- trim_with_ext(paste0(stem, ".csv"), 96L)
  testthat::expect_identical(nchar(out), 96L)
  testthat::expect_identical(tools::file_ext(out), "csv")
  testthat::expect_match(out, "^x+[.]csv$")
  testthat::expect_identical(nchar(tools::file_path_sans_ext(out)), 92L)
})

testthat::test_that("no supported extension is ever partially emitted", {
  for (ext in c("csv", "tsv", "txt", "xlsx", "gct")) {
    for (budget in seq(nchar(ext) + 2L, 40L)) {
      out <- trim_with_ext(paste0(strrep("y", 120L), ".", ext), budget)
      testthat::expect_identical(tools::file_ext(out), ext,
                                 info = paste(ext, budget))
      testthat::expect_lte(nchar(out), budget)
      testthat::expect_gt(nchar(tools::file_path_sans_ext(out)), 0L)
    }
  }
})

testthat::test_that("a declared .csv never yields an empty, .c or .cs suffix", {
  # the three shapes actually found in the PRIDE bundle
  for (budget in 5L:120L) {
    out <- trim_with_ext(paste0(strrep("z", 150L), ".csv"), budget)
    testthat::expect_false(endsWith(out, ".c"), info = budget)
    testthat::expect_false(endsWith(out, ".cs"), info = budget)
    testthat::expect_true(nzchar(tools::file_ext(out)), info = budget)
    testthat::expect_identical(tools::file_ext(out), "csv", info = budget)
  }
})

testthat::test_that("the off-by-one boundaries the old defect survived", {
  name <- paste0(strrep("b", 40L), ".csv")   # 44 characters
  testthat::expect_identical(nchar(name), 44L)

  # budget == nchar(name): untouched
  testthat::expect_identical(trim_with_ext(name, 44L), name)
  # budget == nchar(name) - 1: one stem character goes, extension intact
  out <- trim_with_ext(name, 43L)
  testthat::expect_identical(nchar(out), 43L)
  testthat::expect_identical(tools::file_ext(out), "csv")
  # budget == nchar(name) - 2
  out2 <- trim_with_ext(name, 42L)
  testthat::expect_identical(nchar(out2), 42L)
  testthat::expect_identical(tools::file_ext(out2), "csv")

  # budget == nchar(".csv") + 1 = 5: exactly one stem character plus ".csv"
  out3 <- trim_with_ext(name, 5L)
  testthat::expect_identical(out3, "b.csv")
  testthat::expect_identical(nchar(out3), 5L)

  # budget == nchar(".csv") = 4: impossible, and must say so
  testthat::expect_error(trim_with_ext(name, 4L), "cannot hold a stem character")
  testthat::expect_error(trim_with_ext(name, 1L), "cannot hold a stem character")
})

testthat::test_that("an extensionless input stays extensionless and within budget", {
  out <- trim_with_ext(strrep("q", 150L), 96L)
  testthat::expect_identical(nchar(out), 96L)
  testthat::expect_identical(tools::file_ext(out), "")
  # and a short extensionless name is untouched
  testthat::expect_identical(trim_with_ext("plain_name", 96L), "plain_name")
})

testthat::test_that("a stem containing dots keeps only the final extension", {
  out <- trim_with_ext(paste0("a.b.c.", strrep("d", 150L), ".csv"), 96L)
  testthat::expect_identical(tools::file_ext(out), "csv")
  testthat::expect_identical(nchar(out), 96L)
  testthat::expect_match(out, "[.]csv$")
  # the interior dots are part of the stem and may be trimmed, but the final
  # extension may not
  short <- trim_with_ext("a.b.c.csv", 96L)
  testthat::expect_identical(short, "a.b.c.csv")
})

testthat::test_that("output is deterministic for the same input and budget", {
  name <- paste0(strrep("m", 130L), ".csv")
  a <- trim_with_ext(name, 96L)
  b <- trim_with_ext(name, 96L)
  testthat::expect_identical(a, b)
  testthat::expect_false(identical(trim_with_ext(name, 96L), trim_with_ext(name, 60L)))
})

testthat::test_that("sanitisation still happens and is not bypassed by the fix", {
  # safe_filename() replaces path separators and whitespace and collapses runs
  out <- trim_with_ext("a b/c\\d.csv", 96L)
  testthat::expect_false(grepl("[/\\\\[:space:]]", out))
  testthat::expect_identical(tools::file_ext(out), "csv")
  # a name that sanitises to nothing becomes "unnamed"
  testthat::expect_match(trim_with_ext("___", 96L), "unnamed")
})

testthat::test_that("the trailing-dot shape that Windows silently strips cannot recur", {
  # Two of the 63 malformed files were 95 characters, not 96: the cut landed
  # exactly on the dot and Windows dropped the trailing dot. No output may end
  # in a dot.
  for (budget in 5L:120L) {
    out <- trim_with_ext(paste0(strrep("w", 140L), ".csv"), budget)
    testthat::expect_false(endsWith(out, "."), info = budget)
  }
})

testthat::test_that("the collision path also preserves the extension", {
  stage_names <- get("supplementary_stage_names", envir = budget_env)
  # Two different sources whose staged names collide after trimming must both
  # keep ".csv". Previously the disambiguator read the extension back off the
  # already-truncated name, so it reattached nothing.
  tmp <- withr::local_tempdir()
  long <- strrep("p", 90L)
  f1 <- file.path(tmp, "results", "tables", "04_x", "microglia", paste0(long, "_one.csv"))
  f2 <- file.path(tmp, "results", "tables", "04_x", "microglia", paste0(long, "_two.csv"))
  for (f in c(f1, f2)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines("a,b", f)
  }
  staged <- withr::with_envvar(
    c(PROTEOMICS_PROJECT_ROOT = tmp),
    stage_names(c(f1, f2), datasets = "microglia", max_chars = 96L))
  testthat::expect_length(staged, 2L)
  testthat::expect_identical(unname(vapply(staged, tools::file_ext, character(1))),
                             c("csv", "csv"))
  testthat::expect_true(all(nchar(staged) <= 96L))
  testthat::expect_identical(length(unique(staged)), 2L)
})

testthat::test_that("the staged PRIDE bundle carries no malformed filename", {
  supp <- repo_path("pride_submission", "supplementary_tables")
  testthat::skip_if_not(dir.exists(supp), "PRIDE bundle not present")
  files <- list.files(supp, full.names = FALSE)
  testthat::skip_if(!length(files), "PRIDE supplementary tables empty")

  known <- c("csv", "tsv", "txt", "xlsx", "xls", "gct", "md", "rds", "yml",
             "yaml", "json", "pdf", "svg", "png", "html", "log", "flag", "zip", "gz")
  ext <- tolower(tools::file_ext(files))
  bad <- files[!ext %in% known]
  testthat::expect_identical(length(bad), 0L,
    info = paste("malformed staged filenames:", paste(utils::head(bad, 5), collapse = ", ")))
  # and specifically none of the three observed damage shapes
  testthat::expect_false(any(endsWith(files, ".c")))
  testthat::expect_false(any(endsWith(files, ".cs")))
  testthat::expect_false(any(endsWith(files, ".")))
})

testthat::test_that("the future PRIDE path budget keeps headroom to the wall", {
  supp <- repo_path("pride_submission", "supplementary_tables")
  worst <- trim_with_ext(paste0(strrep("x", 300L), ".csv"), 96L)
  abs_worst <- path_length_chars(file.path(supp, worst))
  testthat::expect_identical(nchar(worst), 96L)
  testthat::expect_lt(abs_worst, PATH_LENGTH_WALL)
  # the budget is deliberately generous: assert real headroom, not a 1-char win
  testthat::expect_gt(PATH_LENGTH_WALL - abs_worst, 40L)
})
