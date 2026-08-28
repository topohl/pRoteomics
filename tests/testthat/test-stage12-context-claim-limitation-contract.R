testthat::local_edition(3)

# The Stage-12 audit enforces, as validation check `context_claim_limits`:
#
#   all(grepl("not|cannot|unresolved|require",
#             context_validation$claim_limitation, ignore.case = TRUE))
#
# That check only exercises the branches reached by the current data, so a
# branch whose wording lacks an explicit limitation can sit latent until a
# module first lands in it. These tests extract the shipped switch() table
# straight from the audit source and prove EVERY branch satisfies the
# contract, including branches no module currently occupies.

AUDIT_SCRIPT <- testthat::test_path(
  "..", "..", "06_modules_WGCNA", "12_microglia_wgcna_nature_readiness_audit.R"
)

# The regex the Stage-12 validation contract uses. Deliberately duplicated
# here: if production widens it, this test must be updated deliberately.
CLAIM_LIMIT_PATTERN <- "not|cannot|unresolved|require"

# Walk the parsed audit source and return the `switch(dominant, ...)` call
# assigned to `limitation`.
find_limitation_switch <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  found <- NULL
  is_limitation_switch <- function(node) {
    parts <- as.list(node)
    length(parts) == 3L &&
      identical(as.character(parts[[1]]), "<-") &&
      is.name(parts[[2]]) &&
      identical(as.character(parts[[2]]), "limitation") &&
      is.call(parts[[3]]) &&
      identical(as.character(as.list(parts[[3]])[[1]]), "switch")
  }
  walk <- function(node) {
    if (!is.null(found) || !is.call(node)) return(invisible(NULL))
    if (isTRUE(tryCatch(is_limitation_switch(node), error = function(e) FALSE))) {
      found <<- as.list(node)[[3]]
      return(invisible(NULL))
    }
    # as.list() keeps empty arguments (e.g. df[i, ]) as empty symbols rather
    # than raising a missing-argument error the way `[[` would.
    for (part in as.list(node)) {
      if (!is.null(found)) break
      if (isTRUE(tryCatch(is.call(part), error = function(e) FALSE))) walk(part)
    }
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  found
}

testthat::test_that("the limitation switch table is locatable in the audit source", {
  testthat::expect_true(file.exists(AUDIT_SCRIPT))
  sw <- find_limitation_switch(AUDIT_SCRIPT)
  testthat::expect_false(is.null(sw))
  testthat::expect_identical(as.character(sw[[1]]), "switch")
  testthat::expect_identical(as.character(sw[[2]]), "dominant")
})

testthat::test_that("every named context class yields an explicit claim limitation", {
  sw <- find_limitation_switch(AUDIT_SCRIPT)
  args <- as.list(sw)[-c(1, 2)]
  named <- names(args)
  classes <- named[!is.na(named) & nzchar(named)]
  testthat::expect_gt(length(classes), 0L)

  for (cls in classes) {
    text <- eval(args[[cls]])
    testthat::expect_type(text, "character")
    testthat::expect_true(nzchar(text), info = cls)
    testthat::expect_match(
      text, CLAIM_LIMIT_PATTERN, ignore.case = TRUE, info = cls
    )
  }
})

testthat::test_that("the unnamed default branch also carries an explicit limitation", {
  sw <- find_limitation_switch(AUDIT_SCRIPT)
  args <- as.list(sw)[-c(1, 2)]
  named <- names(args)
  default <- args[is.null(named) | is.na(named) | !nzchar(named)]
  testthat::expect_length(default, 1L)
  testthat::expect_match(
    eval(default[[1]]), CLAIM_LIMIT_PATTERN, ignore.case = TRUE
  )
})

testthat::test_that("the neuropil-associated branch states the limitation explicitly", {
  sw <- find_limitation_switch(AUDIT_SCRIPT)
  args <- as.list(sw)[-c(1, 2)]
  testthat::expect_true("neuropil-associated_ROI" %in% names(args))
  text <- eval(args[["neuropil-associated_ROI"]])
  # Regression guard: this branch previously described the confound without
  # ever stating a limitation, so it failed context_claim_limits the first
  # time a module was classified into it.
  testthat::expect_match(text, "cannot be assigned specifically to microglia",
                         fixed = TRUE)
  testthat::expect_match(text, CLAIM_LIMIT_PATTERN, ignore.case = TRUE)
  # The biological interpretation itself must be preserved, not softened.
  testthat::expect_match(text, "adjacent neuronal material", fixed = TRUE)
})

testthat::test_that("the emitted context table satisfies the contract when present", {
  out <- testthat::test_path(
    "..", "..", "results", "reviewer_audit",
    "microglia_wgcna_nature_readiness",
    "module_biological_context_validation.csv"
  )
  testthat::skip_if_not(file.exists(out), "Stage-12 context table unavailable")
  context <- utils::read.csv(out, check.names = FALSE)
  testthat::expect_true(all(nzchar(context$claim_limitation)))
  testthat::expect_true(all(grepl(
    CLAIM_LIMIT_PATTERN, context$claim_limitation, ignore.case = TRUE
  )))
})
