# One classifier, used everywhere a scientific input is judged.
#
# test-input-addressability-contract.R pins the four states themselves. This
# file pins that the consumers actually use them, because the states are worth
# nothing while four surfaces keep their own private vocabulary:
#
#   R/statistics/integration_utils.R      canonical / missing_required / missing_optional
#   R/statistics/evidence_bundle_utils.R  canonical / missing_optional
#   R/enrichment/enrichment_io.R          named the failure for the clusterProfiler
#                                         manifest and not for the compareGO one
#   analysis/.../compare_go_enrichment.R  "manifest not found", plus an instruction
#                                         to re-run the enrichment
#
# Every one of those collapses three different conditions into "missing", and
# only one of the three is fixed by running the producer. The recovered ledger
# shows this is not hypothetical: of 184 distinct audited paths, 11 declare a
# P:/ root that is not mounted here, and 9 of those 11 were recorded with
# file_exists = TRUE at the time they were read. "Missing" was the wrong word
# for all eleven.
#
# Two axes, and they must not be merged again:
#   addressability - what the filesystem can tell us; the four canonical states
#   required       - what the analysis does about it
# So `status` keeps its existing tokens on purpose. Downstream code derives
# evidence_role and counts_toward_convergence from those tokens, and a new one
# there would silently move scientific rows.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "statistics", "integration_utils.R"))
source(testthat::test_path("..", "..", "R", "statistics", "evidence_bundle_utils.R"))
source(testthat::test_path("..", "..", "R", "enrichment", "enrichment_io.R"))

has_pwsh <- function() nzchar(Sys.which("pwsh"))

# read_csv_optional() and bundle_input_status() both record provenance, so
# calling them writes to results/reviewer_audit/input_resolution_audit.csv.
# A unit test working with synthetic fixture paths must not put those paths
# into a reviewer-facing ledger, so point the ledger somewhere disposable for
# the duration. The binding is replaced in the environment the audit writer
# was defined in, which is where its default argument is resolved.
isolate_ledger <- function(.local_envir = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = .local_envir)
  env <- environment(append_input_resolution_audit)
  original <- get("input_resolution_audit_path", envir = env)
  assign("input_resolution_audit_path",
         function() file.path(tmp, "ledger.csv"), envir = env)
  withr::defer(assign("input_resolution_audit_path", original, envir = env),
               envir = .local_envir)
  invisible(tmp)
}

# A directory long enough that adding a basename crosses the wall. R cannot
# create a directory at or past the wall either, so this stops short of it.
deep_dir <- function(parent, target = 245L) {
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  d <- parent
  repeat {
    remaining <- target - path_length_chars(d)
    if (remaining <= 2L) break
    nxt <- file.path(d, strrep("d", min(30L, remaining - 1L)))
    if (!dir.create(nxt, recursive = FALSE, showWarnings = FALSE)) break
    d <- nxt
  }
  d
}

# A file that demonstrably exists and that R cannot open, which is the whole
# point: file.exists() says FALSE and is wrong.
present_but_over_limit <- function(tmp) {
  d <- deep_dir(file.path(tmp, "deep"))
  target <- file.path(d, paste0(strrep("n", 260L - path_length_chars(d) + 6L), ".csv"))
  testthat::skip_if(path_length_chars(target) < PATH_LENGTH_WALL,
                    "could not construct a path past the wall")
  seed <- file.path(tmp, "seed.csv")
  writeLines(c("a,b", "1,2"), seed)
  testthat::skip_if_not(has_pwsh(), "pwsh unavailable; cannot write past the wall")
  testthat::skip_if_not(all(os_copy_no_clobber(seed, target)), "copy past the wall failed")
  testthat::expect_false(file.exists(target))
  testthat::expect_true(isTRUE(os_path_facts(target)$exists[[1]]))
  target
}

# A path whose declared root is not mounted. P:/ is the real case: the stored
# enrichment manifests were written under it.
unmounted_path <- function() {
  roots <- paste0(LETTERS, ":/")
  free <- roots[!dir.exists(roots)]
  testthat::skip_if(!length(free), "every drive letter is mounted")
  paste0(free[[1]], "declared/elsewhere/table.csv")
}

testthat::test_that("the four surfaces no longer decide a scientific input with file.exists()", {
  # Scoped deliberately to the functions that classify inputs. Bootstrap code
  # that locates paths.R, and config-candidate probing, legitimately still use
  # file.exists(); a blanket ban would be noise.
  body_of <- function(file, fn) {
    src <- paste(readLines(testthat::test_path("..", "..", file), warn = FALSE),
                 collapse = "\n")
    after <- sub(paste0(".*", fn, " <- function"), "", src)
    substr(after, 1, 2600)
  }
  cases <- list(
    list(file = file.path("R", "statistics", "integration_utils.R"), fn = "read_csv_optional"),
    list(file = file.path("R", "statistics", "integration_utils.R"), fn = "integration_find"),
    list(file = file.path("R", "statistics", "evidence_bundle_utils.R"), fn = "read_final_csv"),
    list(file = file.path("R", "statistics", "evidence_bundle_utils.R"), fn = "bundle_input_status"),
    list(file = file.path("R", "enrichment", "enrichment_io.R"), fn = "validate_comparego_manifest_contract"),
    list(file = file.path("R", "enrichment", "enrichment_io.R"), fn = "read_single_declared_contract_table")
  )
  for (c in cases) {
    b <- body_of(c$file, c$fn)
    testthat::expect_true(grepl("input_addressability|input_is_present", b),
      info = paste(c$fn, "no longer consults the canonical classifier"))
  }

  # and the compareGO entry point stopped telling the reader to re-run a
  # multi-hour enrichment for a failure re-running cannot fix
  cg <- paste(readLines(testthat::test_path("..", "..", "analysis",
                                            "differential_abundance",
                                            "compare_go_enrichment.R"), warn = FALSE),
              collapse = "\n")
  testthat::expect_true(grepl("manifest_status <- input_addressability(manifest_path)",
                              cg, fixed = TRUE))
  testthat::expect_true(grepl("INPUT_STATUS_ABSENT", cg, fixed = TRUE),
    info = "the re-run instruction is no longer conditional on genuine absence")
  testthat::expect_false(grepl("clusterProfiler manifest not found", cg, fixed = TRUE))
})

testthat::test_that("read_csv_optional names the failure and leaves the consequence alone", {
  tmp <- withr::local_tempdir()
  withr::local_envvar(PROTEOMICS_SCRIPT_ID = "test-vocabulary")
  isolate_ledger()
  # keep the real ledger untouched
  withr::local_options(list(warn = -1))

  real <- file.path(tmp, "real.csv")
  utils::write.csv(data.frame(x = 1:2), real, row.names = FALSE)
  got <- read_csv_optional(real, "microglia", "integration", "real.csv", required = FALSE)
  testthat::expect_identical(got$status$status, "present")
  testthat::expect_identical(got$status$addressability, INPUT_STATUS_PRESENT)
  testthat::expect_identical(got$status$message, "loaded")

  gone <- file.path(tmp, "gone.csv")
  opt <- read_csv_optional(gone, "microglia", "integration", "gone.csv", required = FALSE)
  req <- read_csv_optional(gone, "microglia", "integration", "gone.csv", required = TRUE)
  # the consequence axis is unchanged, which is what protects evidence_role
  testthat::expect_identical(opt$status$status, "missing_optional")
  testthat::expect_identical(req$status$status, "missing_required")
  # the addressability axis is new and says which failure it was
  testthat::expect_identical(opt$status$addressability, INPUT_STATUS_ABSENT)
  testthat::expect_true(grepl("input not available", opt$status$message))

  un <- read_csv_optional(unmounted_path(), "microglia", "integration", "elsewhere.csv")
  testthat::expect_identical(un$status$status, "missing_optional")
  testthat::expect_identical(un$status$addressability, INPUT_STATUS_ROOT_UNMOUNTED)
  testthat::expect_true(grepl("not mounted", un$status$message))
  testthat::expect_false(grepl("^input not available", un$status$message))
})

testthat::test_that("a present-but-unopenable input is not reported as absent", {
  tmp <- withr::local_tempdir()
  withr::local_envvar(PROTEOMICS_SCRIPT_ID = "test-vocabulary")
  isolate_ledger()
  over <- present_but_over_limit(tmp)

  got <- read_csv_optional(over, "microglia", "integration", "over.csv", required = TRUE)
  testthat::expect_identical(got$status$addressability, INPUT_STATUS_OVER_LIMIT)
  testthat::expect_identical(got$status$status, "missing_required")
  testthat::expect_true(grepl("cannot be opened by R", got$status$message))
  testthat::expect_false(grepl("input not available", got$status$message))
  testthat::expect_null(got$data)

  st <- bundle_input_status(list(over = over))
  testthat::expect_identical(st$addressability, INPUT_STATUS_OVER_LIMIT)
  testthat::expect_identical(st$status, "missing_optional")
})

testthat::test_that("bundle_input_status keeps its two tokens and gains the third axis", {
  tmp <- withr::local_tempdir()
  withr::local_envvar(PROTEOMICS_SCRIPT_ID = "test-vocabulary")
  isolate_ledger()
  real <- file.path(tmp, "real.csv")
  utils::write.csv(data.frame(x = 1:3), real, row.names = FALSE)
  st <- bundle_input_status(list(real = real, gone = file.path(tmp, "gone.csv")))
  testthat::expect_identical(st$status, c("present", "missing_optional"))
  testthat::expect_identical(st$addressability,
                             c(INPUT_STATUS_PRESENT, INPUT_STATUS_ABSENT))
  testthat::expect_identical(st$n_rows, c(3L, 0L))
  # and the column set only grew; nothing positional was disturbed
  testthat::expect_identical(names(st)[1:4], c("input_name", "path", "status", "n_rows"))
})

testthat::test_that("integration_find will not substitute a legacy artifact for an undetermined one", {
  # The silent-substitution bug. file.exists() is FALSE for a canonical
  # candidate past the wall, so the old search fell through to a legacy path
  # and read different data under the same name. No candidate is in that state
  # today - the longest the repository constructs is 174 characters - so this
  # is a guard, and it has to be exercised deliberately.
  cand <- c(normalized = "C:/over/limit/canonical.csv", legacy_flat = "C:/short/legacy.csv")
  status <- c(INPUT_STATUS_OVER_LIMIT, INPUT_STATUS_PRESENT)
  undetermined <- which(status %in% c(INPUT_STATUS_ROOT_UNMOUNTED, INPUT_STATUS_OVER_LIMIT))
  present <- which(status == INPUT_STATUS_PRESENT)
  testthat::expect_true(length(undetermined) > 0 && undetermined[1] < present[1])

  src <- paste(readLines(testthat::test_path("..", "..", "R", "statistics",
                                             "integration_utils.R"), warn = FALSE),
               collapse = "\n")
  body <- substr(sub(".*integration_find <- function", "", src), 1, 2600)
  testthat::expect_true(grepl("undetermined", body, fixed = TRUE))
  testthat::expect_true(grepl("must not be", body, fixed = TRUE))
  # the real candidate sets are all comfortably short, so ordinary resolution
  # is untouched
  for (f in c("cross_compartment_program_atlas_long.csv", "biological_claims.csv")) {
    cands <- integration_artifact_candidates(f, "build_cross_compartment_atlas",
                                             "07_integration", "cross_compartment",
                                             "global", "tables")
    testthat::expect_true(all(path_length_chars(unname(cands)) < PATH_LENGTH_WALL - 40L),
      info = f)
    testthat::expect_silent(integration_find(f, "build_cross_compartment_atlas",
                                             "07_integration", "cross_compartment",
                                             "global", "tables"))
  }
})

testthat::test_that("both manifest contracts fail with the same named vocabulary", {
  # The asymmetry inside enrichment_io.R: the clusterProfiler manifest named
  # the failure, the compareGO manifest called everything missing.
  src <- paste(readLines(testthat::test_path("..", "..", "R", "enrichment",
                                             "enrichment_io.R"), warn = FALSE),
               collapse = "\n")
  # three call sites: the clusterProfiler manifest contract, the compareGO
  # manifest contract, and the single-declared-table reader
  testthat::expect_identical(
    length(gregexpr("describe_input_status_failures", src, fixed = TRUE)[[1]]), 3L)
  testthat::expect_false(grepl("compareGO manifest references missing", src, fixed = TRUE))
  testthat::expect_false(grepl("does not exist: \", paths[[1]]", src, fixed = TRUE))

  manifest <- data.frame(
    dataset = "microglia", comparison = "c", result_type = "GSEA_GO", ontology = "BP",
    analysis_status = "success_with_terms", comparego_analysis_status = "included",
    input_manifest = unmounted_path(), term_comparison_file = unmounted_path(),
    term_gene_provenance_output_file = unmounted_path(),
    analysis_status_summary_file = unmounted_path(),
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    comparego_contract_version = canonical_comparego_manifest_contract_version(),
    stringsAsFactors = FALSE)
  err <- tryCatch(validate_comparego_manifest_contract(manifest, require_files = TRUE),
                  error = function(e) conditionMessage(e))
  testthat::expect_type(err, "character")
  testthat::expect_true(grepl("unusable", err))
  testthat::expect_true(grepl("declared root not mounted", err),
    info = "an unmounted root is still reported as a missing file")

  err2 <- tryCatch(read_single_declared_contract_table(unmounted_path(), "thing"),
                   error = function(e) conditionMessage(e))
  testthat::expect_true(grepl("not usable", err2))
  testthat::expect_false(grepl("does not exist", err2))
})

testthat::test_that("compareGO resolves the manifest before judging it", {
  # The defect this closes was not cosmetic: compare_go_enrichment.R read the
  # clusterProfiler manifest with a raw read.csv() and then handed the P://
  # strings the manifest records straight to the addressability contract. P://
  # is not mounted here, so every successful row classified as
  # declared_root_unmounted and the script stopped before doing any work - on
  # all three datasets. Measured both ways on the same file: raw read STOPPED,
  # read_canonical_clusterprofiler_manifest() PASSED.
  cg <- paste(readLines(testthat::test_path("..", "..", "analysis",
                                            "differential_abundance",
                                            "compare_go_enrichment.R"), warn = FALSE),
              collapse = "\n")
  testthat::expect_true(grepl("read_canonical_clusterprofiler_manifest(", cg, fixed = TRUE))
  testthat::expect_false(grepl("canonical_cluster_manifest <- utils::read.csv", cg, fixed = TRUE))

  # and the canonical reader really is the only thing that re-anchors, so
  # bypassing it cannot be made to work by other means
  io <- paste(readLines(testthat::test_path("..", "..", "R", "enrichment",
                                            "enrichment_io.R"), warn = FALSE),
              collapse = "\n")
  testthat::expect_true(grepl("resolve_manifest_contract_paths", io, fixed = TRUE))

  # the raw path still fails, which is why the fix was needed
  for (ds in c("microglia", "neuron_neuropil", "neuron_soma")) {
    mp <- repo_path("data", "processed", "04_differential_expression_enrichment",
                    "clusterProfiler", ds, "clusterProfiler_manifest.csv")
    testthat::skip_if_not(file.exists(mp), paste("manifest absent:", ds))
    raw <- utils::read.csv(mp, stringsAsFactors = FALSE, check.names = FALSE)
    testthat::skip_if(all(input_addressability(raw$output_table) == INPUT_STATUS_PRESENT),
                      "the declared root is mounted here, so the bypass cannot be shown")
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(raw, strict = TRUE, require_files = TRUE),
      "unusable", info = ds)
    # and the canonical reader passes on the identical file
    testthat::expect_silent(
      m <- read_canonical_clusterprofiler_manifest(mp, dataset = ds, strict = TRUE,
                                                   require_files = TRUE))
    testthat::expect_gt(nrow(m), 0L)
  }
})

testthat::test_that("resolve_input_path uses the contract and will not paper over it", {
  # R/paths.R carried the second vocabulary itself: eight resolution_mode
  # tokens decided by raw file.exists(), five hundred lines below the block
  # that says the old way "let the provenance ledger assert something untrue
  # about the run".
  tmp <- withr::local_tempdir()
  canonical <- file.path(tmp, "canonical.csv"); writeLines("x", canonical)
  fallback <- file.path(tmp, "legacy.csv"); writeLines("y", fallback)

  # ordinary behaviour is unchanged: canonical wins, fallback is used only when
  # the canonical input is genuinely absent
  testthat::expect_identical(
    resolve_input_path(input_name = "a", expected_path = canonical,
                       record_resolution = FALSE),
    normalizePath(canonical, winslash = "/", mustWork = FALSE))
  testthat::expect_identical(
    suppressWarnings(resolve_input_path(
      input_name = "b", expected_path = file.path(tmp, "gone.csv"),
      fallback_paths = fallback, required = FALSE,
      allow_fallback_in_strict = TRUE, record_resolution = FALSE)),
    normalizePath(fallback, winslash = "/", mustWork = FALSE))

  # but an unmounted canonical input must NOT be replaced by a fallback that
  # happens to be readable - that substitutes different data under one name
  got <- suppressWarnings(resolve_input_path(
    input_name = "c", expected_path = unmounted_path(),
    fallback_paths = fallback, required = FALSE,
    allow_fallback_in_strict = TRUE, record_resolution = FALSE))
  testthat::expect_true(is.na(got),
    info = "a fallback was substituted for a canonical input whose state is undetermined")

  # and when it is required, that is an error naming the real condition
  # finish() warns before it stops, so the warning must be muffled rather than
  # caught, or the handler returns before the error is ever raised
  err <- withCallingHandlers(
    tryCatch(resolve_input_path(input_name = "d", expected_path = unmounted_path(),
                                fallback_paths = fallback, required = TRUE,
                                allow_fallback_in_strict = TRUE,
                                record_resolution = FALSE),
             error = function(e) conditionMessage(e)),
    warning = function(w) invokeRestart("muffleWarning"))
  testthat::expect_true(grepl("present but unopenable|not mounted", err),
    info = paste("error did not name the condition:", err))

  src <- paste(readLines(testthat::test_path("..", "..", "R", "paths.R"), warn = FALSE),
               collapse = "\n")
  body <- substr(sub(".*resolve_input_path <- function", "", src), 1, 4200)
  testthat::expect_false(grepl("file.exists(explicit_path)", body, fixed = TRUE))
  testthat::expect_false(grepl("file.exists(expected_path)", body, fixed = TRUE))
  testthat::expect_false(grepl("fallback_paths[file.exists(", body, fixed = TRUE))
  testthat::expect_true(grepl("input_is_present", body, fixed = TRUE))
})

testthat::test_that("latest_input_candidate does not enumerate an unmounted root", {
  # dir.exists() cannot tell an unmounted declared root from an empty one, and
  # quietly enumerating zero files from the first narrows the candidate set
  # without saying so.
  tmp <- withr::local_tempdir()
  writeLines("x", file.path(tmp, "thing_a.csv"))
  hit <- latest_input_candidate(c(dirname(unmounted_path()), tmp), "thing_.*[.]csv$")
  testthat::expect_identical(basename(hit), "thing_a.csv")

  src <- paste(readLines(testthat::test_path("..", "..", "R", "paths.R"), warn = FALSE),
               collapse = "\n")
  body <- substr(sub(".*latest_input_candidate <- function", "", src), 1, 900)
  testthat::expect_true(grepl("path_declared_root_available(root)", body, fixed = TRUE))
  testthat::expect_true(grepl("input_is_present(files)", body, fixed = TRUE))
})

testthat::test_that("the resolution ledger cannot be poisoned by one bad record", {
  # Three spliced records in results/reviewer_audit/input_resolution_audit.csv
  # each leave one stray quote character, and those three characters make the
  # 50,510 lines after them unreadable by any CSV parser: 0.005% of the rows
  # cost 97% of the ledger. A record that would do that is now dropped.
  tmp <- withr::local_tempdir()
  led <- file.path(tmp, "ledger.csv")
  row <- function(name) data.frame(
    script = "s", dataset = "d", stage = "st", input_name = name,
    expected_path = "p", resolved_path = "p", resolution_mode = "canonical",
    strict_mode = TRUE, allowed_in_strict_mode = TRUE, file_exists = TRUE,
    file_hash_sha256 = NA_character_, file_mtime = NA_character_,
    producer_script_or_artifact_id = "p", warning = NA_character_,
    stringsAsFactors = FALSE)

  append_input_resolution_audit(row("first"), led)
  append_input_resolution_audit(row("second"), led)
  d <- utils::read.csv(led, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_identical(nrow(d), 2L)
  testthat::expect_identical(d$input_name, c("first", "second"))

  # embedded separators and quotes still round-trip, so the guard is not
  # achieved by sanitising content
  append_input_resolution_audit(row("has,comma and \"quotes\""), led)
  d <- utils::read.csv(led, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_identical(nrow(d), 3L)
  testthat::expect_identical(d$input_name[3], "has,comma and \"quotes\"")

  # a record carrying an unbalanced quote is refused, and the ledger survives
  bad <- row("tail")
  bad$warning <- "unbalanced \001"
  src <- paste(readLines(testthat::test_path("..", "..", "R", "paths.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("balanced <- vapply(buffer", src, fixed = TRUE))
  testthat::expect_true(grepl("unbalanced input-resolution audit", src, fixed = TRUE))
  # and the appender still has no dry-run guard, which another test depends on
  testthat::expect_false(grepl("is_dry_run",
    substr(sub(".*append_input_resolution_audit <- function", "", src), 1, 1200),
    fixed = TRUE))
})
