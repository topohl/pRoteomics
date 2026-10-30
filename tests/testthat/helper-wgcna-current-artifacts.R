# Shared helpers for the WGCNA label, claim and integration acceptance tests.
#
# 1. CURRENT ARTIFACTS ONLY. Acceptance tests read the canonical layout
#    results/<domain>/<analysis_id>/<scope>/<child>/ and nothing else. The
#    production resolvers (wgcna_find, integration_find, psd_find) fall back to
#    the historical layout, and in the sensitivity tree that layout still holds
#    July-September artefacts of the previous network, so a test that used them
#    could pass on stale data. Literal historical snapshots live in
#    tests/historical/ and are not part of current acceptance.
#
# 2. EXPECTATIONS FROM THE CURRENT NETWORK. Counts and identities (modules,
#    supermodules, singletons, multi-module blocks) are derived from the
#    current Stage 01 member map, never pinned.
#
# 3. DECLARED LABEL MODE. Stage 07 (analysis/wgcna/summarize_module_interpretation.R)
#    writes results/wgcna/summarize_module_interpretation/<dataset>/manifests/label_review_status.csv
#    with the label mode it ran in. A reviewed-registry expectation is skipped
#    only where that manifest DECLARES provisional labels; where the declared (or
#    default) mode is reviewed, a missing registry is an error, not a skip.

if (!exists("repo_path", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "paths.R"))
}

wgcna_test_datasets <- c("neuron_neuropil", "neuron_soma", "microglia")

# The canonical path of one current artifact. Never falls back to the
# historical layout.
current_artifact <- function(domain, owner, scope = "global", ..., child = "tables") {
  canonical_result_path(domain, owner, scope, child, ...)
}

current_wgcna <- function(owner, scope, ..., child = "tables") {
  current_artifact("wgcna", owner, scope, ..., child = child)
}

current_integration <- function(owner, ..., scope = "global", child = "tables") {
  current_artifact("integration", owner, scope, ..., child = child)
}

current_claims_table <- function(filename = "biological_claims_table.csv") {
  current_artifact("publication_source_data", "build_biological_claims_table",
                   "global", filename)
}

# Skip (not generated) or read. A generated-artifact test skips on a fresh tree;
# it never substitutes a historical file.
read_current_csv <- function(path, label = basename(path)) {
  testthat::skip_if_not(file.exists(path),
                        paste("Current canonical artifact not generated:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NaN", "NULL"))
}

# ------------------------------------------------------ current network shape

current_member_map_path <- function(dataset) {
  current_wgcna("build_wgcna_modules", dataset, "supermodules",
                "wgcna_module_supermodule_annotation.csv")
}

current_member_map <- function(dataset) {
  path <- current_member_map_path(dataset)
  testthat::skip_if_not(file.exists(path),
                        paste("Current WGCNA network not built:", path))
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  out <- unique(data.frame(
    dataset = dataset,
    ModuleID = as.character(raw$ModuleID),
    SupermoduleID = as.character(raw$SupermoduleID),
    stringsAsFactors = FALSE
  ))
  out <- out[!is.na(out$ModuleID) & grepl("^WGCNA_m[0-9]{2}$", out$ModuleID) &
               !is.na(out$SupermoduleID) & grepl("^SM[0-9]{2}$", out$SupermoduleID), ,
             drop = FALSE]
  if (anyDuplicated(out$ModuleID)) {
    stop("Current member map assigns a module to more than one supermodule: ", path,
         call. = FALSE)
  }
  out[order(out$SupermoduleID, out$ModuleID), , drop = FALSE]
}

# Everything an acceptance test needs to know about the current network.
current_network_shape <- function(dataset) {
  m <- current_member_map(dataset)
  n_members <- table(m$SupermoduleID)
  singleton_ids <- names(n_members)[n_members == 1L]
  list(
    member_map = m,
    modules = sort(unique(m$ModuleID)),
    supermodules = sort(unique(m$SupermoduleID)),
    singletons = stats::setNames(m$ModuleID[match(singleton_ids, m$SupermoduleID)],
                                 singleton_ids),
    blocks = sort(names(n_members)[n_members > 1L]),
    members = vapply(split(m$ModuleID, m$SupermoduleID),
                     function(x) paste(sort(x), collapse = ";"), character(1)),
    n_identities = length(unique(m$ModuleID)) + length(unique(m$SupermoduleID))
  )
}

# ------------------------------------------------------- declared label mode

wgcna_label_status_manifest <- function(dataset) {
  current_wgcna("summarize_module_interpretation", dataset,
                "label_review_status.csv", child = "manifests")
}

# The label mode the current WGCNA outputs were built in. The Stage 07 manifest
# is the durable declaration; PROTEOMICS_WGCNA_LABEL_MODE is consulted only when
# the manifest is absent, and an explicit environment value that contradicts
# the manifest is an error.
wgcna_declared_label_mode <- function(dataset = "microglia",
                                      manifest = wgcna_label_status_manifest(dataset)) {
  env <- tolower(trimws(Sys.getenv("PROTEOMICS_WGCNA_LABEL_MODE", unset = "")))
  if (file.exists(manifest)) {
    m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
    mode <- unique(tolower(trimws(as.character(m$label_mode))))
    if (length(mode) != 1L || !mode %in% c("reviewed", "provisional")) {
      stop("Unreadable label_mode in ", manifest, call. = FALSE)
    }
    if (nzchar(env) && !identical(env, mode)) {
      stop("PROTEOMICS_WGCNA_LABEL_MODE=", env, " contradicts the declared label mode '",
           mode, "' in ", manifest, call. = FALSE)
    }
    return(list(mode = mode, declared_by = manifest,
                label_review_status = unique(as.character(m$label_review_status)),
                label_source = unique(as.character(m$label_source))))
  }
  list(mode = if (nzchar(env)) env else "reviewed",
       declared_by = if (nzchar(env)) "PROTEOMICS_WGCNA_LABEL_MODE" else
         "default (no Stage 07 label_review_status.csv)",
       label_review_status = NA_character_, label_source = NA_character_)
}

wgcna_reviewed_registry_path <- function(dataset = "microglia") {
  repo_path("config", "wgcna_labels", paste0(dataset, ".csv"))
}

# Use inside test_that(). Skips ONLY in declared provisional mode; in reviewed
# mode a missing registry is an error.
require_reviewed_registry <- function(dataset = "microglia",
                                      decl = wgcna_declared_label_mode(dataset),
                                      path = wgcna_reviewed_registry_path(dataset)) {
  if (identical(decl$mode, "provisional")) {
    testthat::skip(paste0(
      "Declared provisional label mode for ", dataset, " (", decl$declared_by,
      "): the reviewed registry is not in force in this tree by design."
    ))
  }
  if (!file.exists(path)) {
    stop("Declared label mode for ", dataset, " is reviewed (", decl$declared_by,
         ") but the reviewed registry is missing: ", path, call. = FALSE)
  }
  path
}

# The registry the current microglia labels came from: the reviewed registry in
# reviewed mode, the registry named by the Stage 07 manifest in provisional mode.
wgcna_registry_in_force <- function(dataset = "microglia") {
  decl <- wgcna_declared_label_mode(dataset)
  if (identical(decl$mode, "provisional")) {
    src <- decl$label_source
    if (length(src) != 1L || is.na(src)) {
      stop("Declared provisional mode but no registry named in ", decl$declared_by,
           call. = FALSE)
    }
    if (!file.exists(src)) src <- repo_path(src)
    if (!file.exists(src)) {
      stop("Declared provisional registry is missing: ", decl$label_source, call. = FALSE)
    }
    return(list(mode = "provisional", path = src, allow_provisional = TRUE, decl = decl))
  }
  list(mode = "reviewed", path = require_reviewed_registry(dataset),
       allow_provisional = FALSE, decl = decl)
}

# A reviewed row must carry a completed human review. A placeholder reviewer
# ("pending review ...", "not yet reviewed") or an open manual-review flag
# contradicts adjudication_status = reviewed.
reviewed_rows_are_complete <- function(registry) {
  status <- tolower(trimws(as.character(registry$adjudication_status)))
  reviewer <- trimws(as.character(registry$reviewer))
  manual <- as.logical(registry$manual_review_required)
  placeholder <- is.na(reviewer) | !nzchar(reviewer) |
    grepl("pending|not yet reviewed|carry-over|unreviewed", reviewer, ignore.case = TRUE)
  ok <- !(status == "reviewed" & (placeholder | manual %in% TRUE))
  stats::setNames(ok, paste(registry$level, registry$entity_id))
}

# ------------------------------------- current model and contract artifacts
# The contract tests (test-wgcna-identity-contract.R, test-wgcna-group-effects-contract.R)
# compare the hashes a contract recorded with the model files Stage 05 reads now. These
# use the pipeline's own resolvers (R/wgcna/wgcna_paths.R, sourced by those test files)
# but refuse a historical 06_modules_WGCNA path, which the resolvers fall back to when a
# current file is missing and which in the sensitivity tree holds the previous network.

wgcna_test_current_artifact <- function(path, what) {
  if (grepl("/06_modules_WGCNA/", gsub("\\\\", "/", path), fixed = TRUE)) {
    stop("Resolved a historical 06_modules_WGCNA ", what,
         " where the current artifact is required: ", path, call. = FALSE)
  }
  path
}

wgcna_test_current_state <- function(dataset) {
  wgcna_test_current_artifact(wgcna_final_state(dataset), "frozen state")
}

wgcna_test_current_group_output <- function(dataset, filename) {
  wgcna_test_current_artifact(wgcna_group_effects_artifact(filename, dataset), "Stage 05 output")
}

wgcna_test_current_identity_output <- function(dataset, filename) {
  wgcna_test_current_artifact(wgcna_identity_contract_artifact(filename, dataset),
                              "identity-contract output")
}
