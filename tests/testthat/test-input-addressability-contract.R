# The input addressability contract.
#
# file.exists() answers FALSE for three conditions that demand different
# responses, and the provenance ledger used to record all three as "missing":
#   * the declared root is not mounted here (manifests written at a substituted
#     P:/ root are the real case in this repository),
#   * the path reaches the 260-character wall, so R cannot open a file that is
#     demonstrably present and enumerable,
#   * the file genuinely is not there.
#
# These tests pin the four states, the precedence between them, and the fact
# that neither of the first two may ever be reported as absent.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "utilities", "script_runtime.R"))
source(testthat::test_path("..", "..", "R", "enrichment", "enrichment_io.R"))

# An absolute path of an exact character length. normalizePath() collapses
# repeated separators but leaves an absolute non-existent path otherwise alone,
# so the constructed length is the measured length.
path_of_length <- function(n, root = "C:/") {
  stopifnot(n > nchar(root))
  paste0(root, strrep("a", n - nchar(root)))
}

testthat::test_that("a path of a constructed length measures that length", {
  for (n in c(100L, 239L, 259L, 260L, 261L, 300L)) {
    testthat::expect_identical(path_length_chars(path_of_length(n)), n,
                               info = paste("length", n))
  }
})

testthat::test_that("the four states are distinguished", {
  tmp <- withr::local_tempdir()
  real <- file.path(tmp, "real.csv")
  writeLines("x", real)

  testthat::expect_identical(input_addressability(real), INPUT_STATUS_PRESENT)
  testthat::expect_identical(input_addressability(file.path(tmp, "nope.csv")),
                             INPUT_STATUS_ABSENT)
  testthat::expect_identical(input_addressability("P://data/processed/x.csv"),
                             INPUT_STATUS_ROOT_UNMOUNTED)
  testthat::expect_identical(input_addressability(path_of_length(280L)),
                             INPUT_STATUS_OVER_LIMIT)

  # a directory counts as present
  testthat::expect_identical(input_addressability(tmp), INPUT_STATUS_PRESENT)
  # nothing addressable at all
  testthat::expect_identical(input_addressability(NA_character_), INPUT_STATUS_ABSENT)
  testthat::expect_identical(input_addressability(""), INPUT_STATUS_ABSENT)
})

testthat::test_that("the 259/260/261 boundary sits exactly at the wall", {
  testthat::expect_identical(PATH_LENGTH_WALL, 260L)

  at_259 <- input_addressability(path_of_length(259L))
  at_260 <- input_addressability(path_of_length(260L))
  at_261 <- input_addressability(path_of_length(261L))

  # 259 is under the wall, so the classifier is allowed to consult existence
  # and reports the truth about a file that is not there.
  testthat::expect_false(identical(at_259, INPUT_STATUS_OVER_LIMIT))
  testthat::expect_identical(at_259, INPUT_STATUS_ABSENT)

  testthat::expect_identical(at_260, INPUT_STATUS_OVER_LIMIT)
  testthat::expect_identical(at_261, INPUT_STATUS_OVER_LIMIT)
})

testthat::test_that("precedence is deterministic and total", {
  # root availability outranks the path budget: a path that is both behind an
  # unmounted root and over the wall is reported as the root problem, because
  # its length at a mounted root is not yet known.
  both <- paste0("P:/", strrep("a", 300L))
  testthat::expect_gte(path_length_chars(both), PATH_LENGTH_WALL)
  testthat::expect_identical(input_addressability(both), INPUT_STATUS_ROOT_UNMOUNTED)

  # the path budget outranks existence: never consult file.exists() past the
  # wall, because it answers FALSE for files that are really there.
  testthat::expect_identical(input_addressability(path_of_length(300L)),
                             INPUT_STATUS_OVER_LIMIT)

  # every state is one of the declared levels, for any input
  probe <- c(path_of_length(100L), path_of_length(300L), "P://x", "", NA_character_)
  testthat::expect_true(all(input_addressability(probe) %in% INPUT_STATUS_LEVELS))
})

testthat::test_that("an unmounted root or an over-limit path is never called absent", {
  testthat::expect_false(
    identical(input_addressability("P://data/processed/x.csv"), INPUT_STATUS_ABSENT))
  testthat::expect_false(
    identical(input_addressability(path_of_length(275L)), INPUT_STATUS_ABSENT))
})

testthat::test_that("an over-limit path that really exists is not reported absent", {
  # The decisive case: a file that exists on disk but that R cannot open.
  tmp <- withr::local_tempdir()
  # Grow the directory only up to 245 characters. R cannot create a directory
  # at or past the wall either, so overshooting here would fail to build the
  # fixture rather than test anything.
  deep <- tmp
  repeat {
    remaining <- 245L - path_length_chars(deep)
    if (remaining <= 2L) break
    nxt <- file.path(deep, strrep("d", min(30L, remaining - 1L)))
    if (!dir.create(nxt, recursive = FALSE, showWarnings = FALSE)) break
    deep <- nxt
  }
  testthat::skip_if_not(dir.exists(deep), "could not build a deep fixture directory")
  target <- file.path(deep, "a_filename_long_enough_to_cross_the_wall.csv")
  testthat::skip_if(path_length_chars(target) < PATH_LENGTH_WALL,
                    "fixture did not reach the wall")

  # Written through PowerShell, which is extended-length aware, so the file is
  # genuinely on disk even though R cannot open it.
  invisible(suppressWarnings(system2(
    "pwsh", c("-NoProfile", "-Command",
              paste0("Set-Content -LiteralPath '", gsub("/", "\\\\", target),
                     "' -Value 'x'")),
    stdout = NULL, stderr = NULL)))
  enumerated <- basename(target) %in% list.files(deep)
  testthat::skip_if_not(enumerated, "could not stage an over-wall fixture file")

  testthat::expect_false(file.exists(target))          # the pathology itself
  testthat::expect_identical(input_addressability(target), INPUT_STATUS_OVER_LIMIT)
  testthat::expect_false(identical(input_addressability(target), INPUT_STATUS_ABSENT))
})

testthat::test_that("input_is_present is exactly status == present", {
  probe <- c(path_of_length(100L), path_of_length(300L), "P://x", "", NA_character_,
             testthat::test_path("..", "..", "R", "paths.R"))
  testthat::expect_identical(input_is_present(probe),
                             input_addressability(probe) == INPUT_STATUS_PRESENT)
})

testthat::test_that("input_status_row records the class, not just a boolean", {
  tmp <- withr::local_tempdir()
  real <- file.path(tmp, "real.csv"); writeLines("x", real)

  ok <- input_status_row("present_input", real, required = TRUE)
  testthat::expect_identical(ok$status, INPUT_STATUS_PRESENT)
  testthat::expect_true(ok$input_present)
  testthat::expect_true(ok$required)

  un <- input_status_row("unmounted_input", "P://data/processed/x.csv", required = TRUE)
  testthat::expect_identical(un$status, INPUT_STATUS_ROOT_UNMOUNTED)
  testthat::expect_false(un$input_present)
  testthat::expect_true(un$required)
  testthat::expect_match(un$message, "not mounted")

  ov <- input_status_row("over_limit_input", path_of_length(280L), required = FALSE)
  testthat::expect_identical(ov$status, INPUT_STATUS_OVER_LIMIT)
  testthat::expect_false(ov$input_present)
  testthat::expect_false(ov$required)

  ab <- input_status_row("absent_input", file.path(tmp, "nope.csv"), required = TRUE)
  testthat::expect_identical(ab$status, INPUT_STATUS_ABSENT)
  testthat::expect_false(ab$input_present)

  # the compatibility boolean is derived, never independent
  for (row in list(ok, un, ov, ab)) {
    testthat::expect_identical(row$input_present,
                               identical(row$status, INPUT_STATUS_PRESENT))
  }
  # `required` is carried by its own column, so the failure classes no longer
  # encode it: the same missing file is one status whether required or not.
  req <- input_status_row("x", file.path(tmp, "nope.csv"), required = TRUE)
  opt <- input_status_row("x", file.path(tmp, "nope.csv"), required = FALSE)
  testthat::expect_identical(req$status, opt$status)
  testthat::expect_false(req$required == opt$required)
  testthat::expect_false(any(grepl("required|optional",
                                   c(un$status, ov$status, ab$status))))
})

testthat::test_that("input_status_row still honours a glob and an explicit override", {
  tmp <- withr::local_tempdir()
  writeLines("x", file.path(tmp, "a_table.csv"))
  glob <- input_status_row("glob_input", file.path(tmp, "*.csv"))
  testthat::expect_identical(glob$status, INPUT_STATUS_PRESENT)
  testthat::expect_true(glob$input_present)

  forced <- input_status_row("forced", file.path(tmp, "nope.csv"),
                             status = "skipped_by_configuration")
  testthat::expect_identical(forced$status, "skipped_by_configuration")
})

# --- manifest contract ------------------------------------------------------

make_contract_manifest <- function(dataset, paths) {
  data.frame(
    dataset = dataset,
    comparison = paste0(dataset, "_c1"),
    result_type = "GSEA_GO",
    ontology = "BP",
    analysis_status = "success_with_terms",
    n_terms = 1L,
    output_table = paths[[1]],
    collapsed_gene_input_file = paths[[2]],
    collapsed_gene_provenance_file = paths[[3]],
    term_gene_provenance_file = paths[[4]],
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    gene_annotation_contract_version = "v1",
    stringsAsFactors = FALSE
  )
}

testthat::test_that("the manifest contract names the failure class per dataset", {
  tmp <- withr::local_tempdir()
  good <- vapply(1:4, function(i) {
    p <- file.path(tmp, paste0("ok_", i, ".csv")); writeLines("x", p); p
  }, character(1))

  for (dataset in c("neuron_soma", "neuron_neuropil", "microglia")) {
    # all present -> the contract passes
    testthat::expect_true(validate_clusterprofiler_manifest_contract(
      make_contract_manifest(dataset, as.list(good)), require_files = TRUE))

    # an unmounted declared root is reported as such, never as absent
    unmounted <- as.list(good); unmounted[[2]] <- "P://data/processed/x.csv"
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(
        make_contract_manifest(dataset, unmounted), require_files = TRUE),
      "declared root not mounted")

    # an over-wall path is reported as unopenable, never as absent
    over <- as.list(good); over[[3]] <- path_of_length(280L)
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(
        make_contract_manifest(dataset, over), require_files = TRUE),
      "character limit")

    # a genuinely absent file keeps saying so
    gone <- as.list(good); gone[[4]] <- file.path(tmp, "nope.csv")
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(
        make_contract_manifest(dataset, gone), require_files = TRUE),
      "genuinely absent")
  }
})

testthat::test_that("validation stays all-or-nothing and diagnosis does not excuse it", {
  tmp <- withr::local_tempdir()
  good <- vapply(1:4, function(i) {
    p <- file.path(tmp, paste0("ok_", i, ".csv")); writeLines("x", p); p
  }, character(1))

  # one unusable column out of four is still a failure, for every class
  for (bad in list("P://x.csv", path_of_length(280L), file.path(tmp, "nope.csv"))) {
    mixed <- as.list(good); mixed[[1]] <- bad
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(
        make_contract_manifest("microglia", mixed), require_files = TRUE))
  }

  # require_files = FALSE still skips the check entirely
  unmounted <- as.list(good); unmounted[[1]] <- "P://x.csv"
  testthat::expect_true(validate_clusterprofiler_manifest_contract(
    make_contract_manifest("microglia", unmounted), require_files = FALSE))
})

testthat::test_that("a failure message groups several classes at once", {
  msg <- describe_input_status_failures(
    c("P://a.csv", path_of_length(280L), file.path(tempdir(), "nope_xyz.csv")))
  testthat::expect_match(msg, "declared root not mounted")
  testthat::expect_match(msg, "character limit")
  testthat::expect_match(msg, "genuinely absent")
  testthat::expect_identical(describe_input_status_failures(
    testthat::test_path("..", "..", "R", "paths.R")), "")
})

testthat::test_that("the real enrichment manifests fail as root-unmounted, not as absent", {
  roots <- Sys.glob(repo_path("data", "processed",
                              "04_differential_expression_enrichment",
                              "clusterProfiler", "*", "clusterProfiler_manifest.csv"))
  testthat::skip_if(!length(roots), "no on-disk clusterProfiler manifests")

  for (m in roots) {
    raw <- utils::read.csv(m, stringsAsFactors = FALSE)
    cols <- c("output_table", "collapsed_gene_input_file",
              "collapsed_gene_provenance_file", "term_gene_provenance_file")
    cols <- intersect(cols, names(raw))
    vals <- unique(unlist(raw[cols], use.names = FALSE))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    testthat::skip_if(!length(vals), "manifest declares no paths")

    status <- input_addressability(vals)
    # These manifests were written at a substituted P:/ root that is not
    # mounted here. Whatever else is true, the classifier must not claim the
    # artefacts are gone.
    if (all(grepl("^[Pp]:", vals))) {
      testthat::expect_true(all(status == INPUT_STATUS_ROOT_UNMOUNTED),
                            info = basename(dirname(m)))
      testthat::expect_false(any(status == INPUT_STATUS_ABSENT),
                             info = basename(dirname(m)))
    }
  }
})
