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
  c("01_preprocessing",                      "analysis/01_preprocessing"),
  c("02_id_mapping",                         "analysis/01_preprocessing"),
  c("03_qc_exploration",                     "analysis/02_qc"),
  c("04_differential_expression_enrichment", "analysis/04_differential_abundance"),
  c("05_celltype_enrichment_EWCE",           "analysis/06_gsea"),
  c("06_modules_WGCNA",                      "analysis/05_wgcna"),
  c("07_spatial_networks",                   "analysis/07_spatial_networks"),
  c("08_behavior_physio_coupling",           "analysis/08_integration"),
  c("08_biological_interpretation",          "analysis/08_integration"),
  c("09_export_pride_journal",               "analysis/09_publication_exports"),
  c("10_biological_integration",             "analysis/08_integration"),
  c("11_spatial_systems",                    "analysis/03_spatial_validation")
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

  ## 0. Phase 6D renames
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
  cat("No value, hash, threshold, claim or word of prose changed.\n")
  quit(status = 0L)
}
cat("RESULT: FAIL -", problems, "object(s) changed beyond addressing.\n")
quit(status = 1L)
