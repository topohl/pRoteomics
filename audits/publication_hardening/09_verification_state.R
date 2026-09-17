#!/usr/bin/env Rscript

# PH-010: the verification-state contract.
#
# A test result describes a TREE, not a project. Reporting "the suite passes"
# without saying which tree was tested is how b392977 came to be committed in a
# failing state: the suite ran against a working tree in which a file was still
# untracked, and the check that would have caught the defect enumerates
# candidates with `git ls-files`, so it could not see that file. The commit made
# it visible. The number was true of what was measured and false of what was
# shipped.
#
# This script records, in one artefact, exactly which tree a verification run
# describes.
#
#   Rscript 09_verification_state.R                      # record, never fail
#   Rscript 09_verification_state.R --release            # record, fail if dirty
#   Rscript 09_verification_state.R --release --allow-nonhead-verification
#
# Ordinary development runs are unaffected. Only --release is strict, because
# only a release claim is a claim about the committed state.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("99_audits/publication_hardening/00_checkpoint.R")

args <- commandArgs(trailingOnly = TRUE)
RELEASE <- "--release" %in% args
ALLOW <- "--allow-nonhead-verification" %in% args
NOTE <- {
  i <- match("--note", args)
  if (!is.na(i) && length(args) > i) args[i + 1L] else ""
}

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) character(0)
  else out
}

head_commit <- git("rev-parse", "HEAD")[1]
head_short <- git("rev-parse", "--short", "HEAD")[1]

# Three states, distinguished deliberately.
#   worktree_clean      - nothing modified or untracked at all
#   index_clean         - nothing staged that differs from HEAD
#   index_differs       - staged content differs from HEAD
# A run is only "HEAD" when all tracked content matches HEAD AND nothing
# untracked could change what a git-ls-files-driven test enumerates.
porcelain <- git("status", "--porcelain", "--untracked-files=all")
staged <- git("diff", "--cached", "--name-only")
unstaged <- git("diff", "--name-only")
untracked <- porcelain[grepl("^\\?\\?", porcelain)]

worktree_clean <- length(porcelain) == 0
index_differs_from_head <- length(staged) > 0
index_clean <- !index_differs_from_head

tested_state <- if (worktree_clean) {
  "HEAD"
} else if (index_differs_from_head && length(unstaged) == 0 &&
           length(untracked) == 0) {
  "STAGED_TREE"
} else {
  "DIRTY_WORKTREE"
}

state <- data.frame(
  head_commit = head_commit %||% NA_character_,
  head_short = head_short %||% NA_character_,
  worktree_clean = worktree_clean,
  index_clean = index_clean,
  index_differs_from_head = index_differs_from_head,
  n_staged = length(staged),
  n_unstaged = length(unstaged),
  n_untracked = length(untracked),
  tested_state = tested_state,
  release_mode = RELEASE,
  allow_nonhead = ALLOW,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  notes = NOTE,
  stringsAsFactors = FALSE)

dir.create(PH_OUT, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(state, file.path(PH_OUT, "verification_state.csv"),
                 row.names = FALSE)

cat("===== VERIFICATION STATE =====\n")
cat("HEAD            :", head_short, "\n")
cat("worktree clean  :", worktree_clean, "\n")
cat("index clean     :", index_clean,
    "| differs from HEAD:", index_differs_from_head, "\n")
cat("staged/unstaged/untracked:", length(staged), "/", length(unstaged), "/",
    length(untracked), "\n")
cat("TESTED STATE    :", tested_state, "\n")
if (length(untracked))
  cat("untracked files present - a git ls-files test cannot see these:\n  ",
      paste(sub("^\\?\\? ", "", untracked), collapse = "\n   "), "\n")

if (RELEASE && !identical(tested_state, "HEAD") && !ALLOW) {
  cat("\nRELEASE VERIFICATION REFUSED\n")
  cat("A release benchmark must describe the committed tree. This tree is ",
      tested_state, ".\n", sep = "")
  cat("Either commit first and rerun, or pass --allow-nonhead-verification and\n")
  cat("state the tested state explicitly wherever the result is reported.\n")
  quit(status = 1L)
}
verdict <- if (!RELEASE) {
  "not requested (development run)"
} else if (identical(tested_state, "HEAD")) {
  "OK - result describes HEAD"
} else {
  "OVERRIDDEN - result does NOT describe HEAD"
}
cat("\nrelease verification:", verdict, "\n")
