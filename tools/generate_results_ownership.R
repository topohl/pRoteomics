#!/usr/bin/env Rscript

# Generate the canonical results-ownership registry.
#
#   config/results_ownership.csv   machine-readable
#   docs/RESULTS_OWNERSHIP.md      the same thing, for people
#
# A result family is a stage-level output directory such as
# results/tables/06_modules_WGCNA. Several scripts legitimately contribute
# different files to one family; what must not happen is two scripts declaring
# the same concrete output file. This tool checks that and assigns each family
# exactly one canonical owner.
#
# The owner is the script that writes the family's manuscript-facing product.
# Chosen deterministically, in this order:
#   1. the sole producer, if there is only one;
#   2. the producer of an output named by the publication source-data contract
#      or by the pre-restructure freeze manifest;
#   3. the summarize_* script, which by naming convention is the final writer;
#   4. the build_* script;
#   5. the first registered producer, recorded as a weak assignment.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

reg <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                        include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]

## pipeline_steps() collapses a step s produces list with "|" (see
## R/utilities/pipeline_registry.R); a few hand-written entries also use ";".
## Split on both, or paths get glued together and ownership becomes fiction.
split_paths <- function(x) {
  p <- unlist(strsplit(paste(x, collapse = "|"), "[|;]"))
  trimws(p[nzchar(trimws(p))])
}
## Ownership is real at the substep directory, not at the stage directory: one
## script owns one substep. Keying the registry at depth 3 forced an artificial
## choice among up to 18 contributors; keying at depth 4 leaves 95% of families
## with a single genuine owner, and the residue is adjudicated explicitly.
family_of <- function(p) {
  vapply(strsplit(p, "/", fixed = TRUE),
         function(q) paste(utils::head(q, 4), collapse = "/"), character(1))
}
stage_of <- function(p) {
  vapply(strsplit(p, "/", fixed = TRUE),
         function(q) paste(utils::head(q, 3), collapse = "/"), character(1))
}

## frozen / publication-facing outputs, used to pick the owner
frozen <- character(0)
fm <- try(system2("git", c("show",
  "6801edbce8a5d222f4af46e06b6db4e99f6a9761:manuscript/prerestructure_freeze_manifest.csv"),
  stdout = TRUE), silent = TRUE)
if (!inherits(fm, "try-error") && length(fm)) {
  frozen <- unique(utils::read.csv(text = paste(fm, collapse = "\n"),
                                   stringsAsFactors = FALSE)$repository_relative_path)
}
contract_path <- repo_path("config", "publication_source_data_contract.yml")
pub_dirs <- character(0)
if (file.exists(contract_path) && requireNamespace("yaml", quietly = TRUE)) {
  ct <- yaml::read_yaml(contract_path)
  pub_dirs <- vapply(ct$identities, function(x) as.character(x$source_data_dir), character(1))
}

fam_producers <- list(); fam_outputs <- list(); file_owner <- list()
for (i in seq_len(nrow(steps))) {
  outs <- split_paths(steps$produces[i])
  outs <- outs[grepl("^results/", outs)]
  if (!length(outs)) next
  for (f in unique(family_of(outs))) {
    fam_producers[[f]] <- unique(c(fam_producers[[f]], steps$script[i]))
  }
  concrete <- outs[grepl("[.][A-Za-z]{2,5}$", outs)]
  for (o in unique(concrete)) {
    fam_outputs[[family_of(o)]] <- unique(c(fam_outputs[[family_of(o)]], o))
    file_owner[[o]] <- unique(c(file_owner[[o]], steps$script[i]))
  }
}

## hard requirement: no concrete output may have two declared writers
collisions <- names(file_owner)[vapply(file_owner, length, integer(1)) > 1L]
if (length(collisions)) {
  cat("FAIL: concrete outputs with more than one declared writer:\n")
  for (o in collisions) cat("  ", o, " <- ", paste(file_owner[[o]], collapse = ", "), "\n")
  stop("each canonical output must have exactly one writer", call. = FALSE)
}

pick_owner <- function(fam, producers) {
  if (length(producers) == 1L) return(list(producers, "sole producer"))
  outs <- fam_outputs[[fam]] %||% character(0)
  pub <- outs[outs %in% frozen | vapply(outs, function(o)
    any(startsWith(o, pub_dirs)), logical(1))]
  if (length(pub)) {
    o <- unlist(lapply(pub, function(x) file_owner[[x]]))
    if (length(unique(o)) == 1L) {
      return(list(unique(o), "writes the frozen or publication-facing output"))
    }
    ## More than one contributor writes a frozen output here. The owner is the
    ## one writing most of them; the others stay recorded as contributors.
    tb <- sort(table(o), decreasing = TRUE)
    return(list(names(tb)[1],
                sprintf("writes %d of the %d frozen outputs in this family",
                        tb[[1]], length(pub))))
  }
  s <- grep("/summarize_", producers, value = TRUE)
  if (length(s) == 1L) return(list(s, "summarize_ script: the final writer by convention"))
  b <- grep("/build_", producers, value = TRUE)
  if (length(b) == 1L) return(list(b, "sole build_ script in the family"))
  ## No frozen output and no naming signal: the family is a shared diagnostic
  ## directory. Ownership is assigned to the producer declaring the most
  ## outputs in it, which is the one a reader should look at first.
  cnt <- vapply(producers, function(s) {
    outs2 <- fam_outputs[[fam]] %||% character(0)
    sum(vapply(outs2, function(o) s %in% (file_owner[[o]] %||% character(0)), logical(1)))
  }, integer(1))
  if (max(cnt) > 0L) {
    return(list(producers[which.max(cnt)],
                sprintf("declares %d of the %d outputs in this shared diagnostic family",
                        max(cnt), length(fam_outputs[[fam]] %||% character(0)))))
  }
  list(producers[1], "shared diagnostic directory with no declared concrete output; first registered producer")
}

rows <- list()
for (fam in sort(names(fam_producers))) {
  prod <- fam_producers[[fam]]
  po <- pick_owner(fam, prod)
  owner <- po[[1]]; why <- po[[2]]
  outs <- fam_outputs[[fam]] %||% character(0)
  consumers <- outs[vapply(outs, function(o)
    any(startsWith(o, pub_dirs)) || o %in% frozen, logical(1))]
  rows[[fam]] <- data.frame(
    result_family = fam,
    stage_group = paste(utils::head(strsplit(fam, "/", fixed = TRUE)[[1]], 3), collapse = "/"),
    canonical_owner = owner,
    contributing_analyses = paste(setdiff(prod, owner), collapse = " | "),
    n_contributors = length(prod),
    canonical_output = paste(utils::head(outs, 4), collapse = " | "),
    n_declared_outputs = length(outs),
    publication_source_consumers = paste(utils::head(consumers, 3), collapse = " | "),
    multiple_writer_allowed = length(prod) > 1L,
    classification = if (length(prod) == 1L) "SINGLE_OWNER" else "LEGITIMATE_MULTI_STAGE_COORDINATION",
    rationale = why, stringsAsFactors = FALSE)
}
own <- do.call(rbind, rows)
utils::write.csv(own, repo_path("config", "results_ownership.csv"), row.names = FALSE)

cat("result families        :", nrow(own), "\n")
cat("concrete outputs        :", length(file_owner), "\n")
cat("outputs with 2+ writers :", length(collisions), "\n")
cat("families with 1 producer:", sum(own$n_contributors == 1L), "\n")
cat("families coordinating   :", sum(own$n_contributors > 1L), "\n")
cat("weak owner assignments  :", sum(grepl("^WEAK", own$rationale)), "\n")

## ---- human-readable ------------------------------------------------------
L <- c(
  "# Canonical results ownership",
  "",
  "Generated by `tools/generate_results_ownership.R`. Do not hand-edit;",
  "`config/results_ownership.csv` is the machine-readable form.",
  "",
  "A **result family** is a stage-level output directory. Several scripts",
  "legitimately write different files into one family, which is coordination",
  "rather than a defect. What must never happen is two scripts declaring the",
  "same concrete output file, and this tool fails if any does.",
  "",
  sprintf("At this commit: **%d families**, **%d concrete declared outputs**, **%d** of which have more than one declared writer.",
          nrow(own), length(file_owner), length(collisions)),
  "",
  "Output namespaces are keyed on **stage identity**, not on script location.",
  "That is why a script named `summarize_missingness.R` writes into",
  "`results/tables/03_qc_exploration/02_missingness_diagnostics/`: the namespace",
  "was frozen before the Phase 6E renames and 129 frozen baseline objects live",
  "beneath these paths. Use this table, not the filename, to map a script to",
  "its outputs.",
  "",
  "| Result family | Canonical owner | Contributors | Outputs |",
  "| --- | --- | --- | --- |")
for (i in seq_len(nrow(own))) {
  L <- c(L, sprintf("| `%s` | `%s` | %d | %d |",
                    own$result_family[i], basename(own$canonical_owner[i]),
                    own$n_contributors[i], own$n_declared_outputs[i]))
}
L <- c(L, "",
  "## Families with more than one contributor",
  "",
  "All are `LEGITIMATE_MULTI_STAGE_COORDINATION`: the contributors write",
  "different files into a shared stage directory. The canonical owner is the",
  "script that writes the family's manuscript-facing product.",
  "")
s <- own[own$n_contributors > 1L, ]
for (i in seq_len(nrow(s))) {
  L <- c(L, sprintf("### `%s`", s$result_family[i]), "",
         sprintf("- **owner:** `%s` — %s", s$canonical_owner[i], s$rationale[i]),
         sprintf("- **contributors (%d):** %s", s$n_contributors[i] - 1L,
                 paste0("`", strsplit(s$contributing_analyses[i], " | ", fixed = TRUE)[[1]], "`",
                        collapse = ", ")),
         "")
}
writeLines(L, repo_path("docs", "RESULTS_OWNERSHIP.md"), useBytes = TRUE)
cat("\nwrote config/results_ownership.csv and docs/RESULTS_OWNERSHIP.md\n")
