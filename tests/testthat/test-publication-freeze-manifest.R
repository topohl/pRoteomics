source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "publication_freeze_utils.R"))

repo <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)

# --- deterministic serialisation / stable structure -----------------------

testthat::test_that("file records are sorted, project-relative and stably digested", {
  root <- withr::local_tempdir("freeze_records_")
  mk <- function(rel, txt) {
    p <- file.path(root, rel)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    writeLines(txt, p)
    p
  }
  b <- mk("b.csv", "bbb")
  a <- mk("a.csv", "aaa")
  c2 <- mk("sub/c.csv", "ccc")

  r1 <- freeze_file_records(c(b, a, c2))
  r2 <- freeze_file_records(c(c2, a, b))   # different input order
  testthat::expect_identical(
    vapply(r1, function(x) x$path, character(1)),
    vapply(r2, function(x) x$path, character(1))
  )
  testthat::expect_identical(freeze_set_digest(r1), freeze_set_digest(r2))
  testthat::expect_true(all(vapply(r1, function(x) nzchar(x$sha256), logical(1))))
  testthat::expect_true(all(vapply(r1, function(x) x$size_bytes > 0, logical(1))))
  testthat::expect_identical(freeze_set_digest(list()), NA_character_)
})

testthat::test_that("a changed file changes the set digest", {
  root <- withr::local_tempdir("freeze_digest_")
  p <- file.path(root, "x.csv")
  writeLines("original", p)
  before <- freeze_set_digest(freeze_file_records(p))
  writeLines("mutated", p)
  after <- freeze_set_digest(freeze_file_records(p))
  testthat::expect_false(identical(before, after))
})

testthat::test_that("generated_at is metadata only and never enters an identity section", {
  m1 <- list(metadata = list(generated_at = "2026-01-01T00:00:00+0000"),
             publication_source_data = list(s = list(set_digest_sha256 = "abc")))
  m2 <- list(metadata = list(generated_at = "2099-12-31T23:59:59+0000"),
             publication_source_data = list(s = list(set_digest_sha256 = "abc")))
  testthat::expect_identical(m1$publication_source_data, m2$publication_source_data)

  # And the real manifest keeps generated_at out of every hashed section.
  src <- readLines(repo_path("R", "publication_freeze_utils.R"), warn = FALSE)
  gen_lines <- grep("generated_at", src)
  testthat::expect_gt(length(gen_lines), 0L)
  # No identity/section builder may reference it.
  for (fn in c("freeze_source_data", "freeze_wgcna_identities", "freeze_export_payloads",
               "freeze_provenance_equivalence", "freeze_set_digest")) {
    start <- grep(paste0("^", fn, " <- function"), src)
    testthat::expect_length(start, 1L)
    nxt <- grep("^[a-zA-Z_.][a-zA-Z0-9_.]* <- function", src)
    stop_at <- nxt[nxt > start[[1]]]
    stop_at <- if (length(stop_at)) stop_at[[1]] - 1L else length(src)
    body <- src[start[[1]]:stop_at]
    testthat::expect_false(any(grepl("generated_at", body, fixed = TRUE)), info = fn)
  }
})

# --- package-version mismatch --------------------------------------------

testthat::test_that("a package-version mismatch is detected and fails closed", {
  asserted <- c(clusterProfiler = "4.18.4", fgsea = "1.36.2")
  matching <- function(p) unname(asserted[[p]])
  drifted <- function(p) if (p == "fgsea") "9.9.9" else unname(asserted[[p]])

  ok <- freeze_environment(strict = TRUE, asserted = asserted,
                           observed_fun = matching, observed_r = "4.5.1",
                           observed_bioc = "3.22")
  testthat::expect_true(ok$all_packages_match)
  testthat::expect_true(ok$r_matches)

  testthat::expect_error(
    freeze_environment(strict = TRUE, asserted = asserted, observed_fun = drifted,
                       observed_r = "4.5.1", observed_bioc = "3.22"),
    "does not match the asserted publication environment"
  )
  # Non-strict still reports the mismatch rather than hiding it.
  lax <- freeze_environment(strict = FALSE, asserted = asserted, observed_fun = drifted,
                            observed_r = "4.5.1", observed_bioc = "3.22")
  testthat::expect_false(lax$all_packages_match)
  testthat::expect_false(lax$packages[[which(vapply(lax$packages,
    function(p) p$package == "fgsea", logical(1)))]]$matches)
})

testthat::test_that("an R version mismatch fails closed", {
  asserted <- c(WGCNA = "1.74")
  testthat::expect_error(
    freeze_environment(strict = TRUE, asserted = asserted,
                       observed_fun = function(p) "1.74",
                       observed_r = "4.9.9", observed_bioc = "3.22"),
    "does not match the asserted publication environment"
  )
})

# --- WGCNA identity mismatch ---------------------------------------------

testthat::test_that("a protected WGCNA state mismatch is detected and fails closed", {
  good <- list(
    list(artifact = "wgcna_final_model_state.rds", md5_matches = TRUE),
    list(artifact = "module_group_effects.csv", md5_matches = TRUE)
  )
  bad <- list(
    list(artifact = "wgcna_final_model_state.rds", md5_matches = TRUE),
    list(artifact = "supermodule_group_effects.csv", md5_matches = FALSE)
  )
  testthat::expect_true(freeze_assert_protected_states(good, strict = TRUE))
  testthat::expect_error(
    freeze_assert_protected_states(bad, strict = TRUE),
    "Protected WGCNA state hash mismatch"
  )
  testthat::expect_error(
    freeze_assert_protected_states(bad, strict = TRUE),
    "supermodule_group_effects.csv"
  )
  # Non-strict reports without stopping.
  testthat::expect_false(freeze_assert_protected_states(bad, strict = FALSE))
})

testthat::test_that("all four protected WGCNA artifacts are declared", {
  bn <- freeze_wgcna_protected_basenames()
  testthat::expect_length(bn, 4L)
  testthat::expect_setequal(bn, c(
    "wgcna_final_model_state.rds",
    "wgcna_module_supermodule_annotation.csv",
    "module_group_effects.csv",
    "supermodule_group_effects.csv"
  ))
})

# --- equivalence scope ---------------------------------------------------

testthat::test_that("the historical-to-freeze claim is scoped to the protected files only", {
  expected <- freeze_protected_export_files()
  testthat::expect_length(expected, 6L)

  eq_ok <- list(files = lapply(expected, function(p) list(path = p)))
  testthat::expect_true(freeze_equivalence_scope_ok(eq_ok))

  # Silently widening the claim must be detectable.
  eq_wide <- list(files = lapply(c(expected, "R/statistics/wgcna_group_effects_utils.R"),
                                 function(p) list(path = p)))
  testthat::expect_false(freeze_equivalence_scope_ok(eq_wide))

  eq_narrow <- list(files = lapply(expected[-1], function(p) list(path = p)))
  testthat::expect_false(freeze_equivalence_scope_ok(eq_narrow))
})

testthat::test_that("the equivalence conclusion is stated narrowly", {
  eq <- freeze_provenance_equivalence("HEAD", "HEAD")
  testthat::expect_match(eq$scope, "NOT a whole-repository", fixed = TRUE)
  testthat::expect_match(eq$conclusion, "provenance-label drift", fixed = TRUE)
  testthat::expect_false(grepl("entire repository", eq$conclusion, fixed = TRUE))
  # HEAD vs HEAD is trivially identical, which exercises the comparison path.
  testthat::expect_true(eq$all_protected_files_identical)
  testthat::expect_true(freeze_equivalence_scope_ok(eq))
})

# --- export manifest mismatch --------------------------------------------

testthat::test_that("an export-manifest hash change is detected", {
  root <- withr::local_tempdir("freeze_export_")
  man <- file.path(root, "manifest.csv")
  utils::write.csv(data.frame(source_file = "a", target_file = "b"), man, row.names = FALSE)
  before <- file_hash_sha256(man)
  rows_before <- nrow(utils::read.csv(man))

  utils::write.csv(data.frame(source_file = c("a", "c"), target_file = c("b", "d")),
                   man, row.names = FALSE)
  after <- file_hash_sha256(man)
  rows_after <- nrow(utils::read.csv(man))

  testthat::expect_false(identical(before, after))
  testthat::expect_false(identical(rows_before, rows_after))
})

# --- accepted gaps -------------------------------------------------------

testthat::test_that("known gaps are reported as warnings, never silently ignored", {
  gaps <- freeze_known_gaps()
  # The two PRIDE gaps are unconditional. The renv gap is conditional on
  # measured lockfile completeness, so it is asserted separately below rather
  # than assumed present.
  testthat::expect_gte(length(gaps), 2L)
  ids <- vapply(gaps, function(g) g$id, character(1))
  testthat::expect_true("pride_dry_run_semantics" %in% ids)
  testthat::expect_true("processed_package_wildcard_glob_filter" %in% ids)

  # Every gap carries a severity that maps to WARN, and a non-empty summary.
  for (g in gaps) {
    testthat::expect_identical(g$severity, "warn")
    testthat::expect_true(nzchar(g$summary), info = g$id)
    testthat::expect_true(nzchar(g$classification), info = g$id)
    testthat::expect_match(g$action_in_this_task, "not fixed|untouched")
  }
  # The corrected wildcard-defect location must be recorded.
  wc <- gaps[[which(ids == "processed_package_wildcard_glob_filter")]]
  testthat::expect_match(wc$reference, "export_helpers.R", fixed = TRUE)
  testthat::expect_match(wc$correction_to_prior_description,
                         "not analysis/publication_source_data/build_supplementary_tables.R",
                         fixed = TRUE)
})

testthat::test_that("the renv gap tracks measured lockfile completeness in both directions", {
  # It must clear only because the lockfile is actually complete, and reappear
  # the moment it is not. Driven by fixtures, not by the live lockfile alone.
  state <- freeze_renv_lockfile_state()
  ids_now <- vapply(freeze_known_gaps(state), function(g) g$id, character(1))
  testthat::expect_identical("renv_lock_incomplete" %in% ids_now, !isTRUE(state$complete))

  complete <- state; complete$complete <- TRUE
  testthat::expect_false("renv_lock_incomplete" %in%
    vapply(freeze_known_gaps(complete), function(g) g$id, character(1)))

  incomplete <- state; incomplete$complete <- FALSE
  gaps_bad <- freeze_known_gaps(incomplete)
  ids_bad <- vapply(gaps_bad, function(g) g$id, character(1))
  testthat::expect_true("renv_lock_incomplete" %in% ids_bad)
  renv_gap <- gaps_bad[[which(ids_bad == "renv_lock_incomplete")]]
  testthat::expect_identical(renv_gap$severity, "warn")
  testthat::expect_match(renv_gap$reference, "RENV_LOCK_STATUS", fixed = TRUE)
  # Adding the gap must not drop the unconditional ones.
  testthat::expect_true(all(c("pride_dry_run_semantics",
                              "processed_package_wildcard_glob_filter") %in% ids_bad))
})

# --- freeze identity is anchored to the tag, not HEAD --------------------

testthat::test_that("freeze identity is anchored to the tag so it survives its own commit", {
  st <- freeze_git_state()
  testthat::expect_identical(st$freeze_identity_anchored_to, "tag")
  tag_sha <- freeze_git(paste0("rev-parse ", st$freeze_git_tag, "^{commit}"))
  # The recorded freeze commit must be the tag's commit, never simply HEAD.
  testthat::expect_identical(st$freeze_git_commit, tag_sha)
  testthat::expect_true(nzchar(st$head_commit_at_generation))
  # Anchoring to HEAD would make this manifest self-invalidating once committed.
  testthat::expect_true(st$head_is_freeze_commit_or_descendant)
})

testthat::test_that("an unresolvable freeze tag fails closed", {
  testthat::expect_error(
    freeze_git_state("no-such-freeze-tag-xyz"),
    "does not resolve to a commit"
  )
})

# --- validator is read-only ----------------------------------------------

testthat::test_that("the validator sources contain no mutating call", {
  for (f in c(repo_path("R", "publication_freeze_utils.R"),
              file.path(repo, "tools", "validate_publication_freeze.R"))) {
    src <- readLines(f, warn = FALSE)
    code <- src[!grepl("^\\s*#", src)]
    if (basename(f) == "publication_freeze_utils.R") {
      # Only the writer function may write, and only to the manifest path.
      start <- grep("^write_publication_freeze_manifest <- function", code)
      testthat::expect_length(start, 1L)
      nxt <- grep("^[a-zA-Z_.][a-zA-Z0-9_.]* <- function", code)
      stop_at <- nxt[nxt > start[[1]]]
      stop_at <- if (length(stop_at)) stop_at[[1]] - 1L else length(code)
      writer <- seq.int(start[[1]], stop_at)
      mutating <- grep("writeLines|write\\.csv|write\\.table|file\\.copy|unlink|file\\.remove|saveRDS",
                       code)
      outside <- setdiff(mutating, writer)
      # tempfile handling inside the blob extractor is permitted; it never
      # touches a project path.
      outside <- outside[!grepl("unlink\\(tmp", code[outside])]
      testthat::expect_identical(outside, integer(0),
                                 info = paste(code[outside], collapse = " | "))
    } else {
      testthat::expect_false(any(grepl("writeLines|write\\.csv|write\\.table|file\\.copy|unlink|dir_create",
                                       code)))
    }
  }
})

testthat::test_that("validation of the committed manifest reports no FAIL", {
  manifest <- file.path(repo, "docs", "publication_freeze_manifest.yml")
  testthat::skip_if_not(file.exists(manifest), "freeze manifest not generated yet")

  before <- file_hash_sha256(manifest)
  result <- validate_publication_freeze(manifest)
  after <- file_hash_sha256(manifest)

  testthat::expect_identical(before, after)          # validator did not write
  testthat::expect_gt(result$summary[["PASS"]], 0L)
  testthat::expect_identical(as.integer(result$summary[["FAIL"]]), 0L,
    info = paste(vapply(Filter(function(c) c$status == "FAIL", result$checks),
                        function(c) paste(c$check, c$detail), character(1)),
                 collapse = " | "))
  # The three documented gaps surface as warnings.
  testthat::expect_gte(result$summary[["WARN"]], 3L)
})

testthat::test_that("a missing manifest is a FAIL, not a silent pass", {
  result <- validate_publication_freeze(file.path(tempdir(), "no_such_freeze_manifest.yml"))
  testthat::expect_identical(as.integer(result$summary[["FAIL"]]), 1L)
  testthat::expect_identical(result$checks[[1]]$check, "manifest_present")
})
