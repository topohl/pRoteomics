#!/usr/bin/env Rscript

# Generate audits/phase6e_final_active_tree.csv: every active code file with
# exactly one role and exactly one activity classification.
#
# Role vocabulary (one per file, no ambiguity):
#   ANALYSIS_ENTRYPOINT   a registered pipeline.yml stage script
#   SOURCE_DATA_EXPORTER  a registered script that writes publication source data
#   AUDIT                 verifies a property; produces evidence, not results
#   REUSABLE_LIBRARY      R/ function library, sourced by others
#   REPOSITORY_TOOL       tools/ maintenance or generation utility
#   SUPPORT_SCRIPT        active but not registered; invoked by another file
#
# Activity classification:
#   CANONICAL_ENTRYPOINT  the registry's own entry point for an analysis
#   ACTIVE_SUPPORT        loaded, sourced or invoked by active code
#   PROVENANCE_ONLY       tracked but nothing references it; a record, not code
#
# archive/ is excluded by definition: it is where provenance lives.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
registered <- unique(steps$script)
legacy <- vapply(registry$legacy %||% list(),
                 function(x) as.character(x$script), character(1))

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
produces_of <- list()
for (s in registered) {
  produces_of[[s]] <- split_paths(steps$produces[steps$script == s])
}

tracked <- system2("git", c("ls-files"), stdout = TRUE)
active <- grep("^(analysis|R|audits|tools)/.*[.][Rr]$", tracked, value = TRUE)
active <- active[file.exists(active)]

## text index over everything that could reference a file, archive excluded
refs <- grep("[.]([Rr]|ya?ml|md|tsv|csv|sh|ya?ml)$", tracked, value = TRUE)
refs <- refs[!grepl("^archive/", refs)]
idx <- list()
for (f in refs) {
  if (!file.exists(f)) next
  if (file.info(f)$size > 6e6) next
  idx[[f]] <- paste(readLines(f, warn = FALSE), collapse = "\n")
}

## A mention only counts as activation when it is not inside the file itself
## and not merely a record of history.
RECORD_RX <- "^(audits/(phase6e_|restructure_|test_migration|part29|program_evidence|publication_hardening|migration)|docs/(NAMING_MIGRATION|PRERESTRUCTURE_FREEZE|publication_freeze_manifest))"

callers_of <- function(f) {
  nm <- basename(f)
  hits <- names(idx)[vapply(names(idx), function(g)
    g != f && grepl(nm, idx[[g]], fixed = TRUE), logical(1))]
  hits[!grepl(RECORD_RX, hits)]
}

role_of <- function(f, is_reg, writes_source_data) {
  nm <- basename(f)
  looks_audit <- grepl("^(audit|verify)_", nm) || grepl("_audit[.]", nm)
  if (startsWith(f, "R/")) return("REUSABLE_LIBRARY")
  if (startsWith(f, "tools/")) return("REPOSITORY_TOOL")
  if (startsWith(f, "audits/")) return("AUDIT")
  if (is_reg && writes_source_data) return("SOURCE_DATA_EXPORTER")
  if (is_reg && looks_audit) return("AUDIT")
  if (is_reg) return("ANALYSIS_ENTRYPOINT")
  "SUPPORT_SCRIPT"
}

PUB_RX <- "^results/(publication_)?source_data/"

rows <- vector("list", length(active))
for (i in seq_along(active)) {
  f <- active[[i]]
  is_reg <- f %in% registered
  outs <- produces_of[[f]] %||% character(0)
  wsd <- any(grepl(PUB_RX, outs))
  cal <- callers_of(f)
  tests <- grep("^tests/", cal, value = TRUE)
  code <- setdiff(cal, tests)

  ## A tool or audit with a shebang is invoked by a human, which is a real
  ## activation path: it has no caller by design. Treating "no caller" as dead
  ## would archive working verification tools.
  runnable <- grepl("^#!", substr(idx[[f]] %||% "", 1L, 2L)) ||
    (file.exists(f) && grepl("^#!", readLines(f, warn = FALSE, n = 1L)[1]))
  documented <- any(grepl(basename(f),
    unlist(idx[grep("^(docs/|README)", names(idx))]), fixed = TRUE))

  role <- role_of(f, is_reg, wsd)
  class <- if (is_reg) {
    "CANONICAL_ENTRYPOINT"
  } else if (length(code) || length(tests) || f %in% legacy || runnable) {
    "ACTIVE_SUPPORT"
  } else {
    "PROVENANCE_ONLY"
  }

  rows[[i]] <- data.frame(
    path = f,
    name = basename(f),
    area = if (startsWith(f, "analysis/")) sub("^analysis/([^/]+)/.*", "\\1", f)
           else sub("^([^/]+)/.*", "\\1", f),
    role = role,
    classification = class,
    registered_in_pipeline = is_reg,
    declared_outputs = length(outs),
    writes_publication_source_data = wsd,
    code_callers = length(code),
    test_references = length(tests),
    listed_as_legacy = f %in% legacy,
    human_invoked_entrypoint = runnable,
    documented_in_docs = documented,
    stringsAsFactors = FALSE)
}
tree <- do.call(rbind, rows)
tree <- tree[order(tree$area, tree$path), ]

ROLES <- c("ANALYSIS_ENTRYPOINT", "SOURCE_DATA_EXPORTER", "AUDIT",
           "REUSABLE_LIBRARY", "REPOSITORY_TOOL", "SUPPORT_SCRIPT")
CLASSES <- c("CANONICAL_ENTRYPOINT", "ACTIVE_SUPPORT", "PROVENANCE_ONLY")
stopifnot(all(tree$role %in% ROLES), all(tree$classification %in% CLASSES))

if (!dir.exists("audits")) dir.create("audits")
write.csv(tree, "audits/phase6e_final_active_tree.csv", row.names = FALSE)

cat("active code files:", nrow(tree), "\n")
cat("roles outside the vocabulary:", sum(!tree$role %in% ROLES), "\n")
cat("UNKNOWN or ambiguous roles :", sum(!nzchar(tree$role) | tree$role == "UNKNOWN"), "\n\n")
cat("by role:\n"); print(table(tree$role))
cat("\nby classification:\n"); print(table(tree$classification))
cat("\nrole x classification:\n"); print(table(tree$role, tree$classification))

und <- tree$path[tree$human_invoked_entrypoint & !tree$documented_in_docs &
                 tree$classification == "ACTIVE_SUPPORT"]
cat("\nrunnable tools/audits not mentioned in docs/:", length(und), "\n")
if (length(und)) cat(paste0("  ", und, collapse = "\n"), "\n")

po <- tree[tree$classification == "PROVENANCE_ONLY", ]
cat("\nPROVENANCE_ONLY candidates:", nrow(po), "\n")
if (nrow(po)) cat(paste0("  ", po$path, collapse = "\n"), "\n")
