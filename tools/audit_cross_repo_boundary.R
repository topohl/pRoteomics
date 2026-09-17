#!/usr/bin/env Rscript

# Cross-repository boundary audit.
#
# Asks one question of both repositories: does anything reach into the other at
# RUNTIME? A mention in a comment, a provenance record or a documented
# statement of the interface is not coupling; constructing a path into the other
# tree, or sourcing from it, is.
#
# The distinction is made by stripping comments and then requiring the mention
# to appear inside something that actually resolves a path or loads code:
# source(), sys.source(), file.path(), readLines(), read.csv(), list.files(),
# setwd(), system2() and friends. A bare string in a message, a heading or a
# CSV cell does not qualify.
#
# Directories that exist to hold history are exempt by construction: archive/
# and audits/ in pRoteomics record earlier states of the tree, and the phase
# brief excludes historical and provenance text.
#
# Usage:
#   EXP9_MANUSCRIPT_ROOT=/path/to/Exp9_manuscript Rscript tools/audit_cross_repo_boundary.R
#
# Exit status 0 when runtime cross-repository dependencies are zero.

suppressWarnings(source(file.path("R", "paths.R")))

MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT",
                 unset = normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                                       winslash = "/", mustWork = FALSE))

RESOLVERS <- paste0(
  "(source|sys[.]source|file[.]path|readLines|read[.]csv|read[.]delim|read_csv|",
  "read_yaml|list[.]files|Sys[.]glob|file[.]exists|dir[.]exists|file[.]copy|",
  "setwd|system2|normalizePath|repo_path|path_results)")

## What counts as reaching into the other repository.
PATTERNS_IN_PROTEOMICS <- c(
  "Exp9_manuscript",
  "figure1_bridge_mmmsociability",
  "manuscript_draft",
  "figures/final_truth_v9",
  "figures/figure_0")

PATTERNS_IN_MANUSCRIPT <- c(
  "[.][.]/proteomics", "[.][.]/pRoteomics",
  "Analysis/proteomics", "PROTEOMICS_ROOT")

## Files that exist to record history, and the one importer that is allowed to
## name the other repository because it is run by hand, once.
EXEMPT_PROTEOMICS <- c("^archive/", "^audits/", "^docs/", "^README", "^\\.git",
                       "^tools/audit_cross_repo_boundary[.]R$",
                       "^tools/build_restructure_migration_map[.]R$",
                       "^tools/verify_restructure_equivalence[.]R$",
                       "^tools/verify_path_contract_rewrites[.]R$",
                       "^tools/restructure_pipeline_folders[.]sh$",
                       "^audits/migration/")
EXEMPT_MANUSCRIPT <- c("^provenance/", "^docs/", "^README", "^\\.git",
                       "^tools/import_render_inputs[.]R$",
                       "^tools/verify_source_bundles[.]R$",
                       "^R/vendor/manifest[.]csv$")

scan_repo <- function(root, patterns, exempt, label) {
  cat("--", label, "--\n")
  tracked <- suppressWarnings(system2(
    "git", c("-C", shQuote(root), "ls-files"), stdout = TRUE, stderr = FALSE))
  tracked <- tracked[nzchar(tracked)]
  code <- grep("[.]([Rr]|Rmd|sh|ya?ml)$", tracked, value = TRUE)
  code <- code[!vapply(code, function(f)
    any(vapply(exempt, function(e) grepl(e, f), logical(1))), logical(1))]
  cat("   code files scanned :", length(code), "\n")

  findings <- list()
  for (f in code) {
    p <- file.path(root, f)
    if (!file.exists(p)) next
    ln <- readLines(p, warn = FALSE)
    ln <- sub("#.*$", "", ln)          # drop comments
    for (pat in patterns) {
      hit <- grep(pat, ln, perl = TRUE)
      for (h in hit) {
        if (!grepl(RESOLVERS, ln[h], perl = TRUE)) next   # mention, not a path
        findings[[length(findings) + 1L]] <- data.frame(
          file = f, line = h, pattern = pat,
          text = trimws(substr(ln[h], 1, 110)), stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(findings)) {
    cat("   runtime dependencies: 0\n\n")
    return(0L)
  }
  d <- do.call(rbind, findings)
  cat("   runtime dependencies:", nrow(d), "\n")
  for (i in seq_len(nrow(d))) {
    cat(sprintf("     %s:%d  %s\n", d$file[i], d$line[i], d$text[i]))
  }
  cat("\n")
  nrow(d)
}

cat("Cross-repository boundary audit\n")
cat("===============================\n")
cat("pRoteomics      :", repo_root(), "\n")
cat("Exp9_manuscript :", MR, "\n\n")

n1 <- scan_repo(repo_root(), PATTERNS_IN_PROTEOMICS, EXEMPT_PROTEOMICS,
                "pRoteomics reaching into the manuscript repository")
n2 <- if (dir.exists(MR)) {
  scan_repo(MR, PATTERNS_IN_MANUSCRIPT, EXEMPT_MANUSCRIPT,
            "Exp9_manuscript reaching into the analysis repository")
} else {
  cat("-- manuscript repository not present; skipping that direction --\n\n"); 0L
}

total <- n1 + n2
cat("runtime live cross-repo dependencies:", total, "\n\n")
if (total == 0L) {
  cat("RESULT: PASS - the split is enforced at runtime.\n")
  quit(status = 0L)
}
cat("RESULT: FAIL -", total, "runtime dependency/ies.\n")
quit(status = 1L)
