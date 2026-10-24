testthat::test_that("manifest writing records inputs and session info", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  tmp <- tempfile("manifest-", fileext = ".yml")
  input <- tempfile("input-")
  writeLines("hello", input)
  write_run_manifest(tmp, inputs = list(input = input), outputs = list(out = "none"), parameters = list(test = TRUE))
  testthat::expect_true(file.exists(tmp))
  testthat::expect_true(file.exists(file.path(dirname(tmp), "sessionInfo.txt")))
  txt <- readLines(tmp, warn = FALSE)
  testthat::expect_true(any(grepl("input_hashes", txt)))
})

testthat::test_that("run manifest input hashes keep role names and round-trip", {
  testthat::skip_if_not_installed("yaml")
  source(testthat::test_path("..", "..", "R", "paths.R"))
  dir <- tempfile("manifest-roles-")
  dir.create(dir)
  a <- file.path(dir, "a.txt")
  b <- file.path(dir, "b.txt")
  writeLines("a", a)
  writeLines("b", b)
  man <- file.path(dir, "m.yml")
  write_run_manifest(man, inputs = list(expression_matrix = a, sample_metadata = b))
  y <- yaml::read_yaml(man)
  testthat::expect_identical(names(y$input_hashes), c("expression_matrix", "sample_metadata"))
  testthat::expect_identical(unlist(y$input_hashes),
                             c(expression_matrix = file_hash(a), sample_metadata = file_hash(b)))

  # Unnamed inputs stay a plain sequence; partial or nested names give unique keys.
  write_run_manifest(man, inputs = list(a, b))
  testthat::expect_null(names(yaml::read_yaml(man)$input_hashes))
  write_run_manifest(man, inputs = list(x = list(y = a, b), files = c(a, b)))
  keys <- names(yaml::read_yaml(man)$input_hashes)
  testthat::expect_length(keys, 4L)
  testthat::expect_true(all(nzchar(keys)) && !anyDuplicated(keys))
  # colliding and empty role names: make.unique and the basename fallback
  write_run_manifest(man, inputs = list(x = list(a = a), x.a = b, a))
  testthat::expect_identical(names(yaml::read_yaml(man)$input_hashes), c("x.a", "x.a.1", "a.txt"))
})

testthat::test_that("WGCNA cache freshness check flags a changed input by role", {
  testthat::skip_if_not_installed("yaml")
  source(testthat::test_path("..", "..", "R", "paths.R"))
  # Take the check from the script text; the script itself is not run.
  exprs <- parse(repo_path("analysis", "wgcna", "build_wgcna_modules.R"), keep.source = FALSE)
  def <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("check_wgcna_cached_state_freshness")), as.list(exprs))
  testthat::expect_length(def, 1L)
  check_freshness <- eval(def[[1]][[3]])

  dir <- tempfile("manifest-cache-")
  dir.create(dir)
  staged <- c(expression_matrix = file.path(dir, "e.xlsx"), sample_metadata = file.path(dir, "s.xlsx"))
  for (p in staged) writeLines(basename(p), p)
  man <- file.path(dir, "wgcna_run_manifest.yml")
  write_run_manifest(man, inputs = as.list(staged))
  testthat::expect_identical(check_freshness(man, staged)$status, "fresh")
  writeLines("changed", staged[["expression_matrix"]])
  testthat::expect_error(check_freshness(man, staged), "expression_matrix")

  # Manifests written before role names were kept are not comparable and stay fresh.
  old <- file.path(dir, "old_manifest.yml")
  writeLines(c("input_hashes:", "- 173b3daaee43146c3b01dd036d4fead5", "- 9d61e76870a0c6594ebd9f8a7820b48f"), old)
  testthat::expect_identical(check_freshness(old, staged)$status, "fresh")

  # After a run that reused a stale state, its rewritten hashes do not make the state fresh.
  write_run_manifest(man, inputs = as.list(staged), parameters = list(cache_reuse_check = "stale_allowed"))
  testthat::expect_error(check_freshness(man, staged), "previous run reused a stale")
  testthat::expect_identical(suppressWarnings(check_freshness(man, staged, allow_stale = TRUE))$status, "stale_allowed")
})
