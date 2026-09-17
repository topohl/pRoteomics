source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "ewce_contract_utils.R"))

# Phase 6G, enrichment pilot.
#
# Proving a writer's destination by inspecting pipeline.yml proves nothing: the
# registry is a declaration and the code is the fact. These tests read the
# entrypoint itself, resolve the destinations it constructs, and require the
# declaration to agree. Nothing here runs the scientific analysis.

DOMAIN <- "enrichment"
AID <- "run_ewce_celltype_enrichment"
ENTRYPOINT <- file.path("analysis", DOMAIN, paste0(AID, ".R"))
LEGACY_NAMESPACE <- "05_celltype_enrichment_EWCE"

entrypoint_lines <- function() readLines(repo_path(ENTRYPOINT), warn = FALSE)

# active code only: a historical path may be named in a comment as provenance
active_code <- function() {
  l <- entrypoint_lines()
  paste(l[!grepl("^\\s*#", l)], collapse = "\n")
}

registry_produces <- function() {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s <- s[s$script == ENTRYPOINT, , drop = FALSE]
  p <- trimws(unlist(strsplit(paste(s$produces, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}

# every canonical_result_path()/canonical_work_path() call in the file
resolver_calls <- function() {
  exprs <- parse(repo_path(ENTRYPOINT))
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      fn <- e[[1]]
      if (is.name(fn) && as.character(fn) %in%
            c("canonical_result_path", "canonical_work_path")) {
        found[[length(found) + 1L]] <<- e
      }
      for (i in seq_along(e)) if (!is.null(e[[i]])) try(walk(e[[i]]), silent = TRUE)
    } else if (is.pairlist(e) || is.expression(e) || is.list(e)) {
      for (i in seq_along(e)) if (!is.null(e[[i]])) try(walk(e[[i]]), silent = TRUE)
    }
  }
  for (e in exprs) walk(e)
  found
}

# --- the writer resolves into the normalized namespace -------------------

testthat::test_that("the entrypoint builds destinations through the output-layout API", {
  calls <- resolver_calls()
  testthat::expect_gt(length(calls), 0L)

  # the domain is a literal and it is this domain
  domains <- vapply(calls, function(e) {
    a <- as.list(e)[-1]
    if (length(a) >= 1L && is.character(a[[1]])) a[[1]] else NA_character_
  }, character(1))
  testthat::expect_true(all(domains == DOMAIN),
                        info = paste(unique(domains), collapse = ", "))

  # the analysis identity is the registry's, declared once as a literal
  code <- active_code()
  testthat::expect_match(code, paste0('ANALYSIS_ID <- "', AID, '"'), fixed = TRUE)
})

testthat::test_that("declared output paths equal the paths the code resolves", {
  declared <- registry_produces()
  testthat::expect_gt(length(declared), 0L)

  # Resolve what the code would build for one dataset, using the same API the
  # code uses. <dataset> in the registry stands for the scope segment.
  ds <- "microglia"
  children <- sub("^.*/([a-z_]+)/$", "\\1", declared)
  resolved <- vapply(children, function(ch)
    relative_to(canonical_result_path(DOMAIN, AID, ds, ch)), character(1))
  expected <- sub("<dataset>", ds, sub("/$", "", declared), fixed = TRUE)

  testthat::expect_setequal(unname(resolved), unname(expected))

  # and every declared child is one the contract knows
  testthat::expect_true(all(children %in% output_layout_children()),
                        info = paste(children, collapse = ", "))
})

testthat::test_that("the work intermediate is the cache and nothing else", {
  code <- active_code()
  # the bootstrap cache is regenerable and uncited, so it belongs in work/
  testthat::expect_match(code, "canonical_work_path", fixed = TRUE)
  wp <- relative_to(canonical_work_path(DOMAIN, AID, "microglia", "cache"))
  testthat::expect_identical(wp, "work/enrichment/run_ewce_celltype_enrichment/microglia/cache")

  # EWCE_results_full.rds is read by a downstream diagnostic, so it is a
  # result and must not live under work/
  testthat::expect_false(grepl("canonical_work_path[^\n]*models", code))
  testthat::expect_match(code, 'ewce_result_dir\\("models"\\)')
})

# --- no legacy writes ----------------------------------------------------

testthat::test_that("no enrichment writer targets the historical namespace", {
  code <- active_code()
  testthat::expect_false(grepl(LEGACY_NAMESPACE, code, fixed = TRUE),
                         info = "historical namespace appears in active code")
  # the old factory keyed on stage identity must be gone from this writer
  testthat::expect_false(grepl("create_module_dirs", code, fixed = TRUE))
  testthat::expect_false(grepl("module_paths(", code, fixed = TRUE))

  # the registry must not declare a legacy destination either
  declared <- registry_produces()
  testthat::expect_false(any(grepl(LEGACY_NAMESPACE, declared, fixed = TRUE)))
  testthat::expect_true(all(startsWith(declared, "results/enrichment/")))
})

testthat::test_that("the historical namespace is still named, but only in a comment", {
  # The point is that history stays legible, not that it is erased. It must be
  # legible in a way that cannot become a write: a comment, never a value.
  l <- entrypoint_lines()
  hits <- grep(LEGACY_NAMESPACE, l, fixed = TRUE)
  testthat::expect_gt(length(hits), 0L)
  testthat::expect_true(all(grepl("^\\s*#", l[hits])),
                        info = paste(l[hits][!grepl("^\\s*#", l[hits])], collapse = " | "))
  # and the comment block says what its status now is
  comments <- paste(l[grepl("^\\s*#", l)], collapse = "\n")
  testthat::expect_match(comments, "LEGACY_READ_ONLY", fixed = TRUE)
})

# --- the sensitivity branch stays unexportable ---------------------------

testthat::test_that("a comparison branch resolves outside the canonical scope", {
  source(repo_path("R", "export_helpers.R"))

  canon <- ewce_resolve_run_contract("microglia", "animal", "")
  testthat::expect_true(canon$canonical)

  branch <- ewce_resolve_run_contract("microglia", "sample", "sensitivity")
  testthat::expect_false(branch$canonical)

  # the scope the writer derives for each
  canon_scope <- canon$dataset
  branch_scope <- file.path("comparison", branch$branch, branch$dataset)

  canon_tables <- canonical_result_path(DOMAIN, AID, canon_scope, "tables")
  branch_tables <- canonical_result_path(DOMAIN, AID, branch_scope, "tables")

  testthat::expect_false(is_noncanonical_ewce_export_path(canon_tables))
  testthat::expect_true(is_noncanonical_ewce_export_path(branch_tables))

  # and in the historical namespace, which still holds real artefacts
  testthat::expect_true(is_noncanonical_ewce_export_path(
    path_results("tables", LEGACY_NAMESPACE, "EWCE_E9_comparison", "b", "microglia")))
  testthat::expect_false(is_noncanonical_ewce_export_path(
    path_results("tables", LEGACY_NAMESPACE, "EWCE_E9", "microglia")))
})

testthat::test_that("the export figure roots cover both namespaces and exclude branches", {
  source(repo_path("R", "export_helpers.R"))
  roots <- canonical_ewce_figure_root()
  testthat::expect_gte(length(roots), 2L)

  rel <- relative_to(roots)
  testthat::expect_true(any(rel == "results/enrichment/run_ewce_celltype_enrichment"))
  testthat::expect_true(any(grepl(LEGACY_NAMESPACE, rel, fixed = TRUE)))

  # the freeze-protected figure exporter calls this inside c(), so a vector is
  # required to stay compatible without editing a protected file
  fig <- paste(readLines(repo_path("analysis", "publication_source_data",
                                   "08_export_manuscript_figures.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_match(fig, "canonical_ewce_figure_root()", fixed = TRUE)
})

testthat::test_that("the supplementary export globs reach the normalized tables only one level deep", {
  cfg <- repo_path("analysis", "publication_source_data", "config", "export_config.yml")
  testthat::skip_if(!file.exists(cfg))
  txt <- readLines(cfg, warn = FALSE)
  globs <- grep("run_ewce_celltype_enrichment", txt, value = TRUE)
  testthat::expect_gte(length(globs), 2L)

  # a single * cannot span comparison/<branch>/<dataset>, so a sensitivity
  # branch cannot be picked up by the supplementary export
  for (g in globs) {
    testthat::expect_match(g, "run_ewce_celltype_enrichment/\\*/tables/")
    testthat::expect_false(grepl("comparison", g, fixed = TRUE))
  }
  # the historical globs stay, so pre-migration outputs remain exportable
  testthat::expect_true(any(grepl(paste0(LEGACY_NAMESPACE, "/EWCE_E9/"), txt, fixed = TRUE)))
})
