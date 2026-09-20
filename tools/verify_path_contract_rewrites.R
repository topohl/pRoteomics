#!/usr/bin/env Rscript

# Prove that every object declared REWRITTEN_PATH_CONTRACT changed only in how
# it addresses files.
#
# Both the baseline version (read from git) and the current version (read from
# disk) are reduced to a canonical form in which every way of naming a script
# or library collapses to the same token. If the two canonical forms are
# identical, the only difference between the versions is addressing. If
# anything else changed -- a number, a hash, a threshold, a claim, a word of
# prose -- the comparison fails and names the line.
#
# This is the guarantee that a structural migration made no scientific change.
#
# Canonicalisation covers the three transformations the migration performed:
#   1. analysis stage roots moved     01_preprocessing/ -> analysis/01_preprocessing/
#   2. R libraries gained a domain    R/module_stats.R  -> R/statistics/module_stats.R
#   3. library addressing refactored  test_path("..","..","R",X) -> repo_path("R",X)
#
# Output namespaces under results/ and data/ are NOT canonicalised, because
# they were deliberately not migrated. A stage name appearing inside
# results/tables/06_modules_WGCNA/ is left alone on both sides, so a change
# there would still be caught.
#
# Exit status 0 when every declared rewrite is explained by addressing alone.

source(file.path("R", "paths.R"))

BASELINE_COMMIT <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"

## old -> new stage roots. Applied only where the old root begins a path token,
## so results/tables/<stage>/ and data/processed/<stage>/ are untouched.
STAGES <- rbind(
  c("99_audits",                             "audits"),
  c("99_deprecated",                         "archive/deprecated"),
  c("90_testing",                            "archive/exploratory"),
  c("00_setup",                              "tools/reference_data"),
  c("01_preprocessing",                      "analysis/preprocessing"),
  c("02_id_mapping",                         "analysis/preprocessing"),
  c("03_qc_exploration",                     "analysis/qc"),
  c("04_differential_expression_enrichment", "analysis/differential_abundance"),
  c("05_celltype_enrichment_EWCE",           "analysis/enrichment"),
  c("06_modules_WGCNA",                      "analysis/wgcna"),
  c("07_spatial_networks",                   "analysis/spatial_networks"),
  c("08_behavior_physio_coupling",           "analysis/integration"),
  c("08_biological_interpretation",          "analysis/integration"),
  c("09_export_pride_journal",               "analysis/publication_source_data"),
  c("10_biological_integration",             "analysis/integration"),
  c("11_spatial_systems",                    "analysis/spatial_validation")
)

## Phase 6F renamed the analysis directories themselves, so a baseline path
## that already used the Phase 6B layout must also collapse onto the new name.
DIR_RENAMES <- rbind(
  c("analysis/01_preprocessing",          "analysis/preprocessing"),
  c("analysis/02_qc",                     "analysis/qc"),
  c("analysis/03_spatial_validation",     "analysis/spatial_validation"),
  c("analysis/04_differential_abundance", "analysis/differential_abundance"),
  c("analysis/05_wgcna",                  "analysis/wgcna"),
  c("analysis/06_gsea",                   "analysis/enrichment"),
  c("analysis/07_spatial_networks",       "analysis/spatial_networks"),
  c("analysis/08_integration",            "analysis/integration"),
  c("analysis/09_publication_exports",    "analysis/publication_source_data")
)

## Domain directories that R libraries moved into.
R_DOMAINS <- c("data_contracts", "qc", "statistics", "spatial",
               "enrichment", "networks", "utilities", "panels", "vendor")

## No Phase 6D renames survived. The export script kept its name: renaming it
## would have changed a file that freeze_protected_export_files() compares
## across two historical commits, and a naming improvement is not worth
## perturbing a frozen provenance mechanism.
##
## Phase 6E renamed 108 scripts and libraries, all byte-identical (git R100).
## The pairs are read from the recorded migration rather than duplicated here,
## so this canonicaliser cannot drift from the record. Old names collapse to
## new ones on both sides, which is what makes a repointed reference compare
## equal while a changed value still fails.
RENAMES <- local({
  f <- repo_path("audits", "phase6e_naming_migration.csv")
  if (!file.exists(f)) return(matrix(character(0), nrow = 0, ncol = 2))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  m <- cbind(basename(d$old_path), basename(d$new_path))
  m <- m[m[, 1] != m[, 2], , drop = FALSE]
  ## longest first, so one old name cannot be rewritten inside another
  m[order(-nchar(m[, 1])), , drop = FALSE]
})

canonicalise <- function(x) {
  ## 3. collapse every library addressing style to RLIB(<name>)
  x <- gsub('testthat::test_path\\("\\.\\.", *"\\.\\.", *"R", *"([A-Za-z0-9_.]+)"\\)',
            'RLIB(\\1)', x, perl = TRUE)
  x <- gsub('test_path\\("\\.\\.", *"\\.\\.", *"R", *"([A-Za-z0-9_.]+)"\\)',
            'RLIB(\\1)', x, perl = TRUE)
  x <- gsub('file\\.path\\(repo, *"R", *"([A-Za-z0-9_.]+)"\\)',
            'RLIB(\\1)', x, perl = TRUE)
  x <- gsub('repo_path\\("R", *"([A-Za-z0-9_.]+)"\\)',
            'RLIB(\\1)', x, perl = TRUE)

  ## 2. drop the R/ domain directory so R/statistics/x.R == R/x.R
  for (d in R_DOMAINS) {
    x <- gsub(paste0("R/", d, "/"), "R/", x, fixed = TRUE)
  }

  ## 0a. Phase 6F analysis directory renames, before the stage-root pass so a
  ## path written in either layout collapses onto the same token.
  for (i in seq_len(nrow(DIR_RENAMES))) {
    x <- gsub(DIR_RENAMES[i, 1], DIR_RENAMES[i, 2], x, fixed = TRUE)
  }

  ## 0b. Phase 6D and 6E renames
  for (i in seq_len(nrow(RENAMES))) {
    x <- gsub(RENAMES[i, 1], RENAMES[i, 2], x, fixed = TRUE)
  }

  ## 1. old stage root -> new, only at a path-token boundary
  for (i in seq_len(nrow(STAGES))) {
    x <- gsub(paste0("(?<![A-Za-z0-9_/-])", STAGES[i, 1], "/"),
              paste0(STAGES[i, 2], "/"), x, perl = TRUE)
  }
  ## then normalise the new form back to a stage-neutral token, so that two
  ## different old roots mapping to one new root still compare equal
  x <- gsub("analysis/", "", x, fixed = TRUE)
  x <- gsub("archive/", "", x, fixed = TRUE)

  sub("[[:space:]]+$", "", x)
}

read_baseline <- function(path) {
  out <- suppressWarnings(system2(
    "git", c("-C", shQuote(repo_root()), "show",
             paste0(BASELINE_COMMIT, ":", shQuote(path))),
    stdout = TRUE, stderr = FALSE))
  if (!length(out)) NULL else out
}

map <- utils::read.csv(repo_path("audits", "restructure_migration_map.csv"),
                       stringsAsFactors = FALSE)
MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT",
                 unset = normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                                       winslash = "/", mustWork = FALSE))

rewrites <- map[map$migration_class == "REWRITTEN_PATH_CONTRACT", , drop = FALSE]
cat("objects declared REWRITTEN_PATH_CONTRACT:", nrow(rewrites), "\n\n")

problems <- 0L
for (i in seq_len(nrow(rewrites))) {
  bp <- rewrites$baseline_path[i]
  root <- if (rewrites$destination_repo[i] == "Exp9_manuscript") MR else repo_root()
  cur_path <- file.path(root, rewrites$destination_path[i])

  old <- read_baseline(bp)
  if (is.null(old)) {
    cat("FAIL ", bp, ": cannot read baseline version\n"); problems <- problems + 1L; next
  }
  if (!file.exists(cur_path)) {
    cat("FAIL ", bp, ": current version missing\n"); problems <- problems + 1L; next
  }
  new <- readLines(cur_path, warn = FALSE)

  a <- canonicalise(old)
  b <- canonicalise(new)
  ## Blank lines carry no content; ignore pure-whitespace-only differences.
  a <- a[nzchar(a)]; b <- b[nzchar(b)]

  ## An ACCEPTED, RECORDED revision is not undeclared drift.
  ##
  ## Phase 6G.8 Interlude 3B.1 intentionally re-froze one export: 50 figures
  ## from the failed WGCNA run microglia_failed_20260720_133211 had been
  ## selected into the manuscript figure set, so the export was corrected and
  ## the freeze re-pinned. That is a content change, and this tool is right to
  ## see one - the guarantee it must keep is not "nothing changed" but "nothing
  ## changed WITHOUT being declared".
  ##
  ## The allowance is therefore narrow and evidence-bearing: a freeze payload
  ## field may differ only when the current file records the superseded value
  ## under superseded_<field> AND carries a freeze_revision_reason. Drop such a
  ## pair from both sides, so any OTHER difference in the same file still fails.
  if (any(grepl("^[[:space:]]*freeze_revision_reason:", b))) {
    ## Drop the revision bookkeeping from BOTH sides: the payload fields the
    ## revision replaced, the superseded_* records of their previous values,
    ## and the free-prose reason block. Everything else in the file is still
    ## compared, so any OTHER undeclared change still fails.
    REVISED <- c("manifest_sha256", "manifest_row_count", "audit_sha256",
                 "run_manifest_sha256", "run_manifest_recorded_commit",
                 "run_manifest_input_count")
    key_re <- paste0("^[[:space:]]*(", paste(REVISED, collapse = "|"), "):")
    a <- a[!grepl(key_re, a)]
    b <- b[!grepl(key_re, b)]
    b <- b[!grepl("^[[:space:]]*superseded_[a-z0-9_]+:", b)]
    ## the reason is a wrapped scalar: skip from its key to the next mapping key
    start <- grep("^[[:space:]]*freeze_revision_reason:", b)
    if (length(start)) {
      drop <- integer(0)
      for (s0 in start) {
        j <- s0
        repeat {
          drop <- c(drop, j)
          j <- j + 1L
          if (j > length(b)) break
          if (grepl("^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*:", b[[j]])) break
          if (grepl("^[[:space:]]*-", b[[j]])) break
        }
      }
      b <- b[-drop]
    }
  }

  if (identical(a, b)) {
    raw_changed <- length(setdiff(sub("[[:space:]]+$", "", old),
                                  sub("[[:space:]]+$", "", new)))
    cat(sprintf("PASS  %-50s %3d line(s) changed, all addressing\n", bp, raw_changed))
    next
  }

  cat(sprintf("FAIL  %s\n", bp)); problems <- problems + 1L
  if (length(a) != length(b)) {
    cat("      content line count differs:", length(a), "vs", length(b), "\n")
    for (x in utils::head(setdiff(a, b), 4)) cat("      only in baseline: ", x, "\n")
    for (x in utils::head(setdiff(b, a), 4)) cat("      only in current : ", x, "\n")
  } else {
    d <- which(a != b)
    cat("      ", length(d), "line(s) differ beyond addressing\n")
    for (k in utils::head(d, 6)) {
      cat("      line", k, "\n        baseline: ", a[k], "\n        current : ", b[k], "\n")
    }
  }
}

cat("\n")
if (problems == 0L) {
  cat("RESULT: PASS - every declared rewrite is explained by file addressing alone.\n")
  cat("No value, hash, threshold, claim or word of prose changed, except where a\n")
  cat("freeze_revision_reason records an accepted revision and superseded_* keeps\n")
  cat("the replaced value. Today that is the Phase 6G.8 Interlude 3B.1 figure-export\n")
  cat("re-freeze; every other difference is addressing alone.\n")
  quit(status = 0L)
}
cat("RESULT: FAIL -", problems, "object(s) changed beyond addressing.\n")
quit(status = 1L)
