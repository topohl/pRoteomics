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
