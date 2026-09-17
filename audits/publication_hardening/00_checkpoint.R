#!/usr/bin/env Rscript

# Checkpoint machinery for the publication-hardening audit.
#
# This audit may exceed a session, so its state lives in the repository rather
# than in a conversation. Two artefacts are authoritative:
#   docs/PUBLICATION_HARDENING_PROGRESS.md    - the human-readable notebook
#   results/reports/publication_hardening/checkpoint_state.csv - machine state
#
# ph_checkpoint() rewrites the CSV atomically after each section so an
# interrupted run always leaves a consistent resume point.

PH_OUT <- file.path("results", "reports", "publication_hardening")
PH_TAB <- file.path("results", "tables", "publication_hardening")
PH_STATE <- file.path(PH_OUT, "checkpoint_state.csv")
PH_LOG <- file.path("docs", "PUBLICATION_HARDENING_PROGRESS.md")

ph_sections <- function() data.frame(
  section_id = c(sprintf("A%02d", 1:12), sprintf("B%02d", 13:28)),
  section_title = c(
    "Manuscript-reachable output inventory",
    "Effect / contrast contract",
    "Statistic identity",
    "Biological n / replication",
    "Multiple-testing contract",
    "Null-language audit",
    "Interaction / specificity claims",
    "Claim-language audit",
    "Dataset-specific interpretation rules",
    "Zero / NA / not-evaluable audit",
    "Numerical-floor disclosure",
    "Methods contract table",
    "Repository architecture inventory",
    "Active dependency graph",
    "Structural anti-patterns",
    "Output-model review",
    "Publication layer review",
    "Code organization review",
    "Script naming audit",
    "Identity vs display label contract",
    "One authoritative script rule",
    "Source-data-first figure contract",
    "Target repository structure",
    "Output target structure",
    "Migration risk classification",
    "P0 implementation",
    "Legacy guards",
    "Publication freeze protection"),
  stringsAsFactors = FALSE)

ph_init <- function() {
  dir.create(PH_OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(PH_TAB, recursive = TRUE, showWarnings = FALSE)
  s <- ph_sections()
  s$status <- "NOT_STARTED"
  s$started_at <- ""; s$completed_at <- ""
  s$files_inspected_count <- 0L
  s$outputs_created <- ""; s$finding_ids <- ""
  s$blocking_issue <- ""; s$next_action <- ""
  utils::write.csv(s, PH_STATE, row.names = FALSE)
  invisible(s)
}

ph_read <- function() {
  if (!file.exists(PH_STATE)) return(ph_init())
  utils::read.csv(PH_STATE, stringsAsFactors = FALSE, colClasses = "character")
}

# atomic: write to a temporary file in the same directory, then rename
ph_checkpoint <- function(section_id, status, outputs = "", findings = "",
                          n_files = NA, blocking = "", next_action = "") {
  s <- ph_read()
  i <- match(section_id, s$section_id)
  if (is.na(i)) stop("unknown section: ", section_id, call. = FALSE)
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (!nzchar(s$started_at[i])) s$started_at[i] <- now
  s$status[i] <- status
  if (status %in% c("COMPLETE", "DEFERRED")) s$completed_at[i] <- now
  if (nzchar(outputs)) s$outputs_created[i] <- outputs
  if (nzchar(findings)) s$finding_ids[i] <- findings
  if (!is.na(n_files)) s$files_inspected_count[i] <- as.character(n_files)
  s$blocking_issue[i] <- blocking
  s$next_action[i] <- next_action
  tmp <- paste0(PH_STATE, ".tmp")
  utils::write.csv(s, tmp, row.names = FALSE)
  file.rename(tmp, PH_STATE)
  cat(sprintf("[checkpoint] %s -> %s\n", section_id, status))
  invisible(s)
}

ph_progress <- function() {
  s <- ph_read()
  done <- sum(s$status == "COMPLETE"); def <- sum(s$status == "DEFERRED")
  cat(sprintf("progress: %d COMPLETE, %d DEFERRED, %d remaining of %d\n",
              done, def, nrow(s) - done - def, nrow(s)))
  invisible(s)
}

if (identical(environment(), globalenv()) &&
    !is.null(sys.calls()) && length(commandArgs(trailingOnly = TRUE))) {
  if (commandArgs(trailingOnly = TRUE)[1] == "init") {
    setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
    ph_init(); ph_progress()
  }
}
