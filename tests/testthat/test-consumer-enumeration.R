source(testthat::test_path("..", "..", "R", "paths.R"))

# Phase 6G.1: the consumer enumerator, checked against the case that motivated
# it.
#
# The enrichment pilot missed tests/testthat/test-ewce-dataset-saving.R through
# three searches because that file names only the fragment "EWCE_E9". These
# tests hold the tool to recovering all five known consumers, to distinguishing
# the match classes rather than treating every string hit alike, and to not
# promoting prose or provenance into a runtime dependency.
#
# The committed inventory is the fixture. It is the artefact the migration gate
# reads, so asserting on it is asserting on what actually gates a migration.

INV <- function(id) repo_path("audits", "consumer_inventory", paste0(id, ".csv"))
EWCE <- "run_ewce_celltype_enrichment"

read_inv <- function(id) {
  f <- INV(id)
  testthat::skip_if(!file.exists(f), paste("no consumer inventory for", id))
  utils::read.csv(f, stringsAsFactors = FALSE)
}

# --- the regression that motivated the tool ------------------------------

testthat::test_that("all five known enrichment consumers are recovered", {
  d <- read_inv(EWCE)
  known <- c(
    "analysis/publication_source_data/config/export_config.yml",
    "R/utilities/export_helpers.R",
    "analysis/differential_abundance/test_microglia_targeted_signatures.R",
    "tests/testthat/test-ewce-dataset-saving.R",
    "tests/testthat/test-export-canonical-ewce-contract.R"
  )
  missing <- known[!vapply(known, function(k)
    any(d$consumer_file == k & d$runtime_relevant), logical(1))]
  if (length(missing)) {
    testthat::fail(paste0("consumer enumeration regressed; not recovered:\n  ",
                          paste(missing, collapse = "\n  ")))
  }
  testthat::expect_length(missing, 0L)
})

testthat::test_that("the helper that synthesises the namespace is found", {
  # EWCE_E9 originates in ewce_contract_utils.R, not in any path literal. A
  # search for a path cannot find it; the derived-vocabulary search can.
  d <- read_inv(EWCE)
  testthat::expect_true(any(d$consumer_file == "R/enrichment/ewce_contract_utils.R" &
                            d$runtime_relevant))
})

# --- match classes are distinguished ------------------------------------

testthat::test_that("the enumerator distinguishes its match classes", {
  d <- read_inv(EWCE)
  # every class present must be one the tool declares
  allowed <- c("FULL_PATH", "PATH_FRAGMENT", "NAMESPACE_BASENAME",
               "OUTPUT_FILENAME", "CONFIG_GLOB", "REGISTRY_DEPENDENCY",
               "HELPER_FUNCTION", "GENERATED_PATH", "TEST_FRAGMENT",
               "DOCUMENTATION_ONLY", "PROVENANCE_ONLY", "ANALYSIS_REFERENCE")
  testthat::expect_true(all(d$match_type %in% allowed),
                        info = paste(setdiff(d$match_type, allowed), collapse = ", "))

  # and the classes that matter for this analysis are actually exercised
  for (cls in c("FULL_PATH", "NAMESPACE_BASENAME", "OUTPUT_FILENAME",
                "CONFIG_GLOB", "TEST_FRAGMENT", "REGISTRY_DEPENDENCY")) {
    testthat::expect_true(cls %in% d$match_type, info = cls)
  }
})

testthat::test_that("a config glob is recognised as a glob, not a path", {
  d <- read_inv(EWCE)
  g <- d[d$consumer_file == "analysis/publication_source_data/config/export_config.yml", ]
  testthat::expect_gt(nrow(g), 0L)
  # the file also mentions the analysis by name, which is not a glob
  g <- g[g$match_type == "CONFIG_GLOB", , drop = FALSE]
  testthat::expect_gt(nrow(g), 0L)
  testthat::expect_true(all(g$dependency_kind == "CONFIG"))
  testthat::expect_true(all(g$runtime_relevant))
})

testthat::test_that("a bare namespace fragment in a test is found and classed as a test", {
  d <- read_inv(EWCE)
  t <- d[d$consumer_file == "tests/testthat/test-export-canonical-ewce-contract.R", ]
  testthat::expect_gt(nrow(t), 0L)
  # runtime-relevant hits in a test file are test dependencies; a comment in
  # the same file is documentation, which is the distinction that matters
  testthat::expect_true(all(t$dependency_kind[t$runtime_relevant] == "TEST"))
  testthat::expect_true(any(t$match_type %in% c("TEST_FRAGMENT", "OUTPUT_FILENAME")))
  testthat::expect_true(all(t$dependency_kind[!t$runtime_relevant] == "DOCUMENTATION" |
                            t$match_type[!t$runtime_relevant] == "ANALYSIS_REFERENCE"))
})

testthat::test_that("the registry dependency is recognised from pipeline.yml", {
  d <- read_inv(EWCE)
  p <- d[d$consumer_file == "pipeline.yml", ]
  testthat::expect_gt(nrow(p), 0L)
  testthat::expect_true(all(p$dependency_kind == "CONFIG"))
  testthat::expect_true(any(p$match_type == "REGISTRY_DEPENDENCY"))
})

# --- negative and context cases -----------------------------------------

testthat::test_that("prose and comments are recorded but never runtime-relevant", {
  d <- read_inv(EWCE)
  doc <- d[d$dependency_kind == "DOCUMENTATION", , drop = FALSE]
  testthat::expect_gt(nrow(doc), 0L)
  testthat::expect_false(any(doc$runtime_relevant))

  # match_type says what the token looked like; dependency_kind says what the
  # reference means. A full path written in prose is still a full path, and it
  # is the kind, not the shape, that decides whether anything must be repointed.
  testthat::expect_true(all(grepl("[.]md$", doc$consumer_file) |
                            doc$match_type %in% c("DOCUMENTATION_ONLY",
                                                  "ANALYSIS_REFERENCE")))
  # and a markdown file never yields a runtime dependency
  md <- d[grepl("[.]md$", d$consumer_file), , drop = FALSE]
  testthat::expect_gt(nrow(md), 0L)
  testthat::expect_false(any(md$runtime_relevant))
})

testthat::test_that("provenance and migration records are not promoted", {
  d <- read_inv(EWCE)
  prov <- d[d$dependency_kind == "PROVENANCE", , drop = FALSE]
  testthat::expect_gt(nrow(prov), 0L)
  testthat::expect_false(any(prov$runtime_relevant))
  testthat::expect_true(all(prov$active_or_historical == "historical"))
})

testthat::test_that("a generated record is separated from a hand-maintained one", {
  d <- read_inv(EWCE)
  # config/results_ownership.csv is regenerated by a tool: it repoints itself
  gen <- d[d$dependency_kind == "GENERATED_RECORD", , drop = FALSE]
  testthat::expect_gt(nrow(gen), 0L)
  testthat::expect_false(any(gen$runtime_relevant))
  testthat::expect_true("config/results_ownership.csv" %in% gen$consumer_file)
})

testthat::test_that("a reference to the script itself is not a dependency on its outputs", {
  d <- read_inv(EWCE)
  ar <- d[d$match_type == "ANALYSIS_REFERENCE", , drop = FALSE]
  testthat::expect_gt(nrow(ar), 0L)
  testthat::expect_false(any(ar$runtime_relevant))
})

# --- the gate -------------------------------------------------------------

testthat::test_that("the inventory leaves no dependency unclassified", {
  d <- read_inv(EWCE)
  testthat::expect_identical(sum(d$dependency_kind == "UNKNOWN"), 0L)
  testthat::expect_false(any(!nzchar(d$dependency_kind)))
  testthat::expect_false(any(!nzchar(d$match_type)))
})

testthat::test_that("the inventory carries the columns the gate reads", {
  d <- read_inv(EWCE)
  needed <- c("analysis_id", "consumer_file", "line_or_expression", "match_type",
              "matched_token", "dependency_kind", "active_or_historical",
              "runtime_relevant", "confidence", "adjudication")
  testthat::expect_true(all(needed %in% names(d)),
                        info = paste(setdiff(needed, names(d)), collapse = ", "))
  testthat::expect_type(d$runtime_relevant, "logical")
})

testthat::test_that("the writer itself is not reported as its own consumer", {
  d <- read_inv(EWCE)
  testthat::expect_false(any(d$consumer_file ==
    "analysis/enrichment/run_ewce_celltype_enrichment.R"))
})

# --- the tooling exists and is runnable ---------------------------------

testthat::test_that("the enumerator and preflight are present and parse", {
  for (f in c("tools/enumerate_output_consumers.R",
              "tools/preflight_writer_migration.R")) {
    p <- repo_path(f)
    testthat::expect_true(file.exists(p), info = f)
    testthat::expect_error(parse(p), NA, info = f)
  }
})
