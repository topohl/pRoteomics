source(testthat::test_path("..", "..", "R", "paths.R"))

# Phase 6C boundary guard.
#
# pRoteomics owns scientific inference. Manuscript prose and journal figure
# assembly moved to the Exp9_manuscript repository, and must not reappear here.
# The scientific repository's responsibility ends at canonical publication
# source data plus a provenance manifest.

testthat::test_that("no manuscript renderer layer remains in the scientific repository", {
  testthat::expect_false(dir.exists(repo_path("figures")),
                         info = "figures/ is the journal figure assembly layer and belongs to Exp9_manuscript")
  testthat::expect_false(dir.exists(repo_path("manuscript")),
                         info = "manuscript/ holds prose and claim provenance and belongs to Exp9_manuscript")

  # Generation-named panel libraries are renderer internals, not science.
  generations <- c("candidate_figure", "editorial_v8", "final_truth_v9",
                   "nature_final_v7", "nature_v2", "story_v3", "story_v4",
                   "story_v5", "manuscript_figure_utils", "manuscript_figure3_utils")
  libs <- list.files(repo_path("R"), pattern = "[.][Rr]$", recursive = TRUE)
  offenders <- libs[vapply(libs, function(f) {
    any(vapply(generations, function(g) grepl(g, basename(f), fixed = TRUE), logical(1)))
  }, logical(1))]
  testthat::expect_identical(as.character(offenders), character(0),
                             info = paste("renderer panel libraries still present:",
                                          paste(offenders, collapse = ", ")))
})

testthat::test_that("the pipeline registry no longer registers figure renderers", {
  testthat::skip_if_not_installed("yaml")
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))

  active <- pipeline_registry_entries(registry)$script
  legacy <- vapply(registry$legacy %||% list(),
                   function(x) as.character(x$script), character(1))
  all_registered <- c(active, legacy)

  renderers <- all_registered[grepl("^figures/", all_registered)]
  testthat::expect_identical(renderers, character(0),
                             info = paste("figure renderers still registered:",
                                          paste(renderers, collapse = ", ")))

  testthat::expect_false("manuscript_candidates" %in% pipeline_stage_names(registry),
                         info = "the manuscript_candidates rendering stage belongs to Exp9_manuscript")
})

testthat::test_that("the publication source-data interface exists and is self-describing", {
  manifest <- path_results("publication_source_data", "manifest.csv")
  testthat::skip_if_not(file.exists(manifest),
                        "publication source-data bundle not built in this checkout")

  m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  required <- c("publication_id", "source_repo", "source_commit", "source_analysis",
                "source_table", "exported_file", "rows", "columns", "sha256",
                "contract_version")
  testthat::expect_true(all(required %in% names(m)),
                        info = paste("missing manifest columns:",
                                     paste(setdiff(required, names(m)), collapse = ", ")))

  # Every exported file must exist and match its recorded hash: the manuscript
  # repository trusts nothing else.
  paths <- repo_path(m$exported_file)
  testthat::expect_true(all(file.exists(paths)),
                        info = "manifest names files that are not in the bundle")
  actual <- vapply(paths, function(p) {
    if (file.exists(p)) unname(tools::sha256sum(p)) else NA_character_
  }, character(1), USE.NAMES = FALSE)
  testthat::expect_identical(actual, m$sha256)

  # Identities are stable publication names, never renderer generations.
  testthat::expect_true(all(grepl("^(figure|extended_data)_[0-9]{2}$",
                                  unique(m$publication_id))),
                        info = paste("non-canonical publication ids:",
                                     paste(unique(m$publication_id), collapse = ", ")))
})

testthat::test_that("no producer layer sources archived or scaffold code", {
  # Ported from test-publication-hardening.R, which moved to Exp9_manuscript
  # with the publication layer in Phase 6C. The assertion is unchanged in
  # direction and strength: deprecated code must not masquerade as active by
  # being sourced from a layer that produces results. Only the layer names
  # changed, because 99_deprecated and 90_testing became archive/.
  layer_of <- function(p) {
    top <- sub("/.*", "", p)
    if (top == "archive") "archive"
    else if (top == "audits") "audit layer"
    else if (top == "analysis") "analysis stage"
    else if (top == "R") "shared helper library"
    else if (top == "tests") "test suite"
    else "other"
  }
  producer <- c("analysis stage", "shared helper library")

  tracked <- suppressWarnings(system2(
    "git", c("-C", shQuote(repo_root()), "ls-files"), stdout = TRUE, stderr = FALSE))
  scripts <- grep("[.][Rr]$", tracked, value = TRUE)
  testthat::expect_gt(length(scripts), 0L)

  SRC <- paste0("(?<![A-Za-z0-9_.])(source|sys[.]source)\\s*\\(\\s*",
                "(repo_path\\s*\\(([^)]*)\\)|[\"']([^\"']+)[\"'])")
  leaks <- character(0)
  for (p in scripts) {
    if (!layer_of(p) %in% producer) next
    ln <- readLines(repo_path(p), warn = FALSE)
    ln <- sub("#.*$", "", ln)

    ## Two active scripts print operator instructions whose text contains
    ## source('archive/...'). That is documentation, not a dependency. Stripping
    ## every string literal would make this guard vacuous, because a real call
    ## carries its target in a string too. So keep only matches that begin
    ## outside a string, judged by the parity of the double quotes before them.
    m <- character(0)
    for (line in ln) {
      hits <- gregexpr(SRC, line, perl = TRUE)[[1]]
      if (identical(hits[1], -1L)) next
      lens <- attr(hits, "match.length")
      for (j in seq_along(hits)) {
        prefix <- substr(line, 1L, hits[j] - 1L)
        quotes <- lengths(regmatches(prefix, gregexpr('"', prefix, fixed = TRUE)))
        if (quotes %% 2L == 1L) next          # inside a double-quoted string
        m <- c(m, substr(line, hits[j], hits[j] + lens[j] - 1L))
      }
    }
    if (!length(m)) next
    tgt <- vapply(m, function(x) {
      if (grepl("repo_path", x)) {
        a <- sub(".*repo_path\\s*\\(", "", x)
        paste(gsub("[\"' )]", "", strsplit(a, ",")[[1]]), collapse = "/")
      } else sub(".*[\"']([^\"']+)[\"'].*", "\\1", x)
    }, character(1))
    bad <- tgt[vapply(tgt, layer_of, character(1)) == "archive"]
    if (length(bad)) leaks <- c(leaks, sprintf("%s -> %s", p, bad))
  }
  testthat::expect_identical(as.character(leaks), character(0),
                             info = paste("producer layer sources archived code:",
                                          paste(leaks, collapse = ", ")))
})

testthat::test_that("no active path depends on a generation namespace", {
  # Historical generation names are allowed in provenance and archive records,
  # which is where version history belongs. They must not appear in the active
  # analysis or library namespace.
  active <- c(list.files(repo_path("analysis"), recursive = TRUE),
              list.files(repo_path("R"), recursive = TRUE),
              list.files(repo_path("config"), recursive = TRUE))
  # Each token must sit at a path boundary. Matching bare substrings would
  # flag 12b_finalize_..._audit.R for containing "_final", and a guard that
  # cries wolf gets switched off.
  BOUND <- "([._/-]|$)"
  banned <- c(paste0("final_truth", BOUND), "_v[789]([^0-9]|$)",
              paste0("latest", BOUND), paste0("_final", BOUND),
              paste0("_new", BOUND), paste0("_old", BOUND),
              paste0("backup", BOUND), paste0("revised", BOUND))
  offenders <- active[vapply(active, function(f) {
    any(vapply(banned, function(b) grepl(b, f, perl = TRUE), logical(1)))
  }, logical(1))]
  testthat::expect_identical(as.character(offenders), character(0),
                             info = paste("generation namespace in an active path:",
                                          paste(offenders, collapse = ", ")))
})
