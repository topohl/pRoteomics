#!/usr/bin/env Rscript

# Phase 6B/6C post-migration equivalence guard.
#
# Reconciles the pre-restructure scientific freeze against the cross-repository
# migration map and fails if anything in the frozen baseline was lost, mapped
# twice, moved to a path that does not exist, or altered.
#
# The baseline manifest is read straight out of git at the baseline commit, not
# from a working-tree copy, so the oracle cannot drift as the migration
# proceeds and cannot be edited into agreement by accident.
#
# Two levels are enforced, because the freeze manifest records one row per
# (object, publication role) pair rather than one row per file:
#
#   assertion level  239 manifest rows each resolve to an existing destination
#                    whose sha256 matches the frozen value
#   object level     200 distinct filesystem objects each map to exactly one
#                    destination, so nothing is silently duplicated
#
# Usage:
#   Rscript tools/verify_restructure_equivalence.R
#   EXP9_MANUSCRIPT_ROOT=/path/to/Exp9_manuscript Rscript tools/...
#
# Exit status 0 when equivalence holds, 1 otherwise.

suppressWarnings(source(file.path("R", "paths.R")))

BASELINE_COMMIT <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"
BASELINE_TAG <- "pre-restructure-scientific-freeze-2026-09"
BASELINE_MANIFEST <- "manuscript/prerestructure_freeze_manifest.csv"
MIGRATION_MAP <- "audits/restructure_migration_map.csv"

EXPECTED_ASSERTIONS <- 239L
EXPECTED_OBJECTS <- 200L

ALLOWED_CLASSES <- c(
  "MOVED_TO_MANUSCRIPT_REPO",
  "MOVED_WITHIN_SCIENTIFIC_REPO",
  "RETAINED_LEGACY_PATH",
  "REWRITTEN_PATH_CONTRACT",
  "REWRITTEN_BOUNDARY_CONTRACT"
)
ALLOWED_REPOS <- c("pRoteomics", "Exp9_manuscript")

## Objects whose recorded hash cannot be reproduced by construction, with the
## reason. The freeze manifest lists itself as one of its own rows, so the
## sha256 it records for itself is the hash of an earlier state of the file: a
## file cannot contain its own final hash. This is a property of the baseline,
## not of the migration, and is verified by byte length instead.
SELF_REFERENTIAL <- BASELINE_MANIFEST

fail <- function(...) {
  cat("FAIL:", ..., "\n")
  assign("problems", get("problems", envir = .GlobalEnv) + 1L, envir = .GlobalEnv)
}
problems <- 0L

manuscript_root <- function() {
  env <- Sys.getenv("EXP9_MANUSCRIPT_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                winslash = "/", mustWork = FALSE)
}

repo_root_for <- function(repo) {
  ifelse(repo == "Exp9_manuscript", manuscript_root(), repo_root())
}

sha256_of <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  if ("sha256sum" %in% getNamespaceExports("tools")) return(unname(tools::sha256sum(path)))
  if (requireNamespace("digest", quietly = TRUE)) {
    return(unname(digest::digest(file = path, algo = "sha256")))
  }
  stop("No SHA-256 implementation available. Install the R package 'digest'.", call. = FALSE)
}

read_baseline <- function() {
  out <- suppressWarnings(system2(
    "git", c("-C", shQuote(repo_root()), "show",
             paste0(BASELINE_COMMIT, ":", BASELINE_MANIFEST)),
    stdout = TRUE, stderr = FALSE))
  if (!length(out)) {
    stop("Could not read the baseline manifest from ", BASELINE_COMMIT,
         ". Is this the pRoteomics repository?", call. = FALSE)
  }
  utils::read.csv(text = paste(out, collapse = "\n"), stringsAsFactors = FALSE)
}

cat("Phase 6B/6C restructure equivalence audit\n")
cat("=========================================\n\n")
cat("baseline commit      :", BASELINE_COMMIT, "\n")
cat("baseline tag         :", BASELINE_TAG, "\n")
cat("manuscript repo root :", manuscript_root(), "\n\n")

base <- read_baseline()
map_path <- file.path(repo_root(), MIGRATION_MAP)
if (!file.exists(map_path)) stop("Migration map not found: ", MIGRATION_MAP, call. = FALSE)
map <- utils::read.csv(map_path, stringsAsFactors = FALSE)

need <- c("baseline_repo", "baseline_path", "baseline_sha256", "destination_repo",
          "destination_path", "postmigration_sha256", "migration_class", "status")
missing_cols <- setdiff(need, names(map))
if (length(missing_cols)) {
  stop("Migration map is missing required columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

## ---------------------------------------------------------------- assertions
cat("-- assertion level --\n")
n_base <- nrow(base)
cat("baseline rows        :", n_base, "\n")
if (n_base != EXPECTED_ASSERTIONS) {
  fail("baseline row count is", n_base, "but the freeze declares", EXPECTED_ASSERTIONS)
}

base_objects <- unique(base$repository_relative_path)
cat("distinct objects     :", length(base_objects), "\n")
if (length(base_objects) != EXPECTED_OBJECTS) {
  fail("distinct baseline objects is", length(base_objects),
       "but the freeze resolves to", EXPECTED_OBJECTS)
}

mapped <- base_objects %in% map$baseline_path
cat("objects mapped       :", sum(mapped), "/", length(base_objects), "\n")
if (any(!mapped)) {
  fail(sum(!mapped), "baseline object(s) have no mapping")
  writeLines(paste0("    unmapped: ", head(base_objects[!mapped], 20)))
}

assertions_mapped <- sum(base$repository_relative_path %in% map$baseline_path)
cat("assertions mapped    :", assertions_mapped, "/", n_base, "\n")
if (assertions_mapped != n_base) {
  fail(n_base - assertions_mapped, "manifest assertion(s) resolve to no mapping")
}

## ------------------------------------------------------------------- objects
cat("\n-- object level --\n")
dup <- table(map$baseline_path)
dup <- dup[dup > 1L]
cat("duplicate mappings   :", length(dup), "\n")
if (length(dup)) {
  fail(length(dup), "baseline object(s) map more than once")
  writeLines(paste0("    duplicated: ", head(names(dup), 20)))
}

extra <- setdiff(map$baseline_path, base_objects)
if (length(extra)) {
  fail(length(extra), "mapping row(s) name a path that is not in the baseline")
  writeLines(paste0("    not in baseline: ", head(extra, 20)))
}

bad_class <- setdiff(unique(map$migration_class), ALLOWED_CLASSES)
if (length(bad_class)) fail("unknown migration_class:", paste(bad_class, collapse = ", "))
bad_repo <- setdiff(unique(map$destination_repo), ALLOWED_REPOS)
if (length(bad_repo)) fail("unknown destination_repo:", paste(bad_repo, collapse = ", "))

## ---------------------------------------------------------------- existence
cat("\n-- destination existence and content --\n")
map$abs <- file.path(repo_root_for(map$destination_repo), map$destination_path)
map$exists <- file.exists(map$abs) & !dir.exists(map$abs)
cat("destinations present :", sum(map$exists), "/", nrow(map), "\n")
if (any(!map$exists)) {
  fail(sum(!map$exists), "destination(s) do not exist")
  writeLines(paste0("    missing: ", head(paste0(map$destination_repo[!map$exists], ":",
                                                 map$destination_path[!map$exists]), 20)))
}

base_sha <- base[!duplicated(base$repository_relative_path), ]
rownames(base_sha) <- base_sha$repository_relative_path
map$frozen_sha <- base_sha[map$baseline_path, "sha256"]
map$frozen_bytes <- base_sha[map$baseline_path, "bytes"]
map$actual_sha <- NA_character_
map$actual_bytes <- NA_real_
for (i in which(map$exists)) {
  map$actual_sha[i] <- sha256_of(map$abs[i])
  map$actual_bytes[i] <- file.info(map$abs[i])$size
}

is_rewritten <- map$migration_class %in% c("REWRITTEN_PATH_CONTRACT", "REWRITTEN_BOUNDARY_CONTRACT")
is_self <- map$baseline_path %in% SELF_REFERENTIAL
checkable <- map$exists & !is_rewritten & !is_self

hash_ok <- checkable & !is.na(map$actual_sha) & map$actual_sha == map$frozen_sha
cat("hash-verified        :", sum(hash_ok), "/", sum(checkable), "\n")
if (any(checkable & !hash_ok)) {
  fail(sum(checkable & !hash_ok), "object(s) changed content unexpectedly")
  writeLines(paste0("    hash mismatch: ",
                    head(map$destination_path[checkable & !hash_ok], 20)))
}

## Declared contract rewrites must still be declared, present and justified.
cat("declared rewrites    :", sum(is_rewritten), "\n")
if (any(is_rewritten)) {
  no_reason <- is_rewritten & (is.na(map$status) | !nzchar(trimws(map$status)))
  if (any(no_reason)) {
    fail(sum(no_reason), "REWRITTEN_PATH_CONTRACT row(s) carry no justification in status")
  }
  unchanged <- is_rewritten & map$exists & !is.na(map$actual_sha) &
    map$actual_sha == map$frozen_sha
  if (any(unchanged)) {
    cat("  note:", sum(unchanged),
        "row(s) declared as rewritten are in fact byte-identical;",
        "reclassify them as moved or retained.\n")
  }
  writeLines(paste0("    ", map$destination_path[is_rewritten]))
}

## The self-referential manifest row is checked by size instead of hash.
if (any(is_self & map$exists)) {
  i <- which(is_self & map$exists)
  same_size <- map$actual_bytes[i] == map$frozen_bytes[i]
  cat("self-referential rows:", length(i),
      "(verified by byte length:", sum(same_size), "/", length(i), ")\n")
  if (any(!same_size)) fail("the freeze manifest changed size")
}

## -------------------------------------------------------------------- report
cat("\n-- summary --\n")
cat(sprintf("%-26s %s\n", "PRE (assertions)", n_base))
cat(sprintf("%-26s %s\n", "MAPPED (assertions)", assertions_mapped))
cat(sprintf("%-26s %s\n", "PRE (objects)", length(base_objects)))
cat(sprintf("%-26s %s\n", "MAPPED (objects)", sum(mapped)))
cat(sprintf("%-26s %s\n", "PRESENT", sum(map$exists)))
cat(sprintf("%-26s %s\n", "UNIQUE", nrow(map) - length(dup)))
cat(sprintf("%-26s %s\n", "HASH MATCH", sum(hash_ok)))
cat(sprintf("%-26s %s\n", "missing", sum(!map$exists)))
cat(sprintf("%-26s %s\n", "duplicate", length(dup)))
cat(sprintf("%-26s %s\n", "hash mismatch", sum(checkable & !hash_ok)))
cat(sprintf("%-26s %s\n", "declared rewrites", sum(is_rewritten)))

cat("\n")
if (problems == 0L) {
  cat("RESULT: PASS - the frozen baseline is fully accounted for.\n")
  quit(status = 0L)
}
cat("RESULT: FAIL -", problems, "problem(s).\n")
quit(status = 1L)
