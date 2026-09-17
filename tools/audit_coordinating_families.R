#!/usr/bin/env Rscript

# Phase 6G.3 section 4/18: evidence for every coordinating result family.
#
#   Rscript tools/audit_coordinating_families.R --domain spatial_validation
#   -> audits/phase6g_coordinating_families_<domain>.csv
#
# Why a separate tool. config/results_ownership.csv classifies a family as
# LEGITIMATE_MULTI_STAGE_COORDINATION when several analyses write into one
# directory. That classification is about a *directory*, and a directory shared
# by two writers is normal. The question a migration has to answer is narrower
# and is about a *file*: does any single artifact have two writers? Only that
# is a collision, and only that is a hard stop.
#
# So this tool separates the two cases the brief insists on distinguishing:
#
#   one canonical owner + contributors writing distinct artifacts  -> legitimate
#   two scripts writing the same canonical artifact                -> STOP
#
# One structural note, recorded because it decides how much the collision gate
# can possibly find here. The normalized layout is keyed on analysis identity:
# results/<domain>/<analysis_id>/<scope>/<lifecycle>/. Two different analyses
# therefore cannot produce the same normalized destination, whatever they are
# named. A post-migration collision is only reachable if two scripts share an
# analysis_id, which the registry forbids. That makes the *current* layout the
# only place a real duplicate-writer can exist, so the gate is exercised
# against today's declared destinations as well as tomorrow's proposals.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "migration_gate_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  i <- which(args == flag)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
DOMAIN <- arg_value("--domain")
if (!nzchar(DOMAIN)) stop("usage: --domain <domain, e.g. spatial_validation>", call. = FALSE)

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]

own <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                       stringsAsFactors = FALSE)

analysis_id_of <- function(script) sub("[.][Rr]$", "", basename(script))
domain_of <- function(script) sub("^analysis/([^/]+)/.*", "\\1", script)

declared_of <- function(script) {
  i <- which(steps$script == script)
  if (!length(i)) return(character(0))
  split_paths(steps$produces[i[1]])
}

# --- proposed normalized destination for one historical output -------------
#
# Deterministic and mechanical, so the proposal can be audited rather than
# trusted. Lifecycle comes from the historical top-level kind, not from the
# file extension, because section 8 forbids inferring lifecycle from extension.
lifecycle_of <- function(p) {
  if (startsWith(p, "results/tables/")) return("tables")
  if (startsWith(p, "results/figures/")) return("plots")
  if (startsWith(p, "results/source_data/")) return("tables/source_data")
  if (startsWith(p, "results/reports/")) return("reports")
  if (startsWith(p, "results/manifests/")) return("manifests")
  if (startsWith(p, "data/processed/")) return("models")
  NA_character_
}

# Scope: a dataset token if the historical path carries one, else global.
DATASETS <- valid_datasets()
scope_of <- function(p) {
  seg <- strsplit(p, "/", fixed = TRUE)[[1]]
  hit <- intersect(seg, DATASETS)
  if (length(hit)) hit[1] else "global"
}

proposed_of <- function(script, p) {
  lc <- lifecycle_of(p)
  if (is.na(lc)) return(NA_character_)
  base <- basename(p)
  root <- file.path("results", domain_of(script), analysis_id_of(script),
                    scope_of(p), lc)
  ## a declared directory (no extension) keeps its lifecycle root only.
  ## file.path() with a zero-length argument returns character(0), so the
  ## directory case has to short-circuit rather than pass NULL through.
  if (!grepl("[.][A-Za-z0-9]{2,5}$", base)) return(root)
  file.path(root, base)
}

# --- coordinating families touching this domain ----------------------------
co <- own[own$classification == "LEGITIMATE_MULTI_STAGE_COORDINATION", , drop = FALSE]
in_domain <- vapply(seq_len(nrow(co)), function(i) {
  scripts <- c(co$canonical_owner[i], split_paths(co$contributing_analyses[i]))
  any(domain_of(scripts) == DOMAIN)
}, logical(1))
co <- co[in_domain, , drop = FALSE]

rows <- list()
gate_rows <- list()

for (i in seq_len(nrow(co))) {
  fam <- co$result_family[i]
  owner <- co$canonical_owner[i]
  contribs <- setdiff(split_paths(co$contributing_analyses[i]), owner)
  scripts <- unique(c(owner, contribs))

  ## what each script declares, and what it proposes
  current <- lapply(scripts, declared_of)
  names(current) <- scripts
  ## restrict to the outputs that actually live in this family
  current <- lapply(current, function(p) p[startsWith(p, fam)])
  proposals <- lapply(scripts, function(s) {
    p <- current[[s]]
    if (!length(p)) return(character(0))
    unname(stats::na.omit(vapply(p, function(x) proposed_of(s, x), character(1))))
  })
  names(proposals) <- scripts

  ## file-level overlap, today and after migration
  cur_flat <- unlist(current, use.names = FALSE)
  cur_shared <- unique(cur_flat[duplicated(cur_flat)])
  cur_shared <- cur_shared[grepl("[.][A-Za-z0-9]{2,5}$", basename(cur_shared))]

  gate <- migration_destination_collisions(proposals, own)
  blockers <- migration_gate_blockers(gate)

  for (s in scripts) {
    p_cur <- current[[s]]
    if (!length(p_cur)) p_cur <- NA_character_
    for (x in p_cur) {
      rows[[length(rows) + 1L]] <- data.frame(
        result_family = fam,
        role = if (identical(s, owner)) "CANONICAL_OWNER" else "CONTRIBUTOR",
        script = s,
        analysis_id = analysis_id_of(s),
        current_output = x,
        proposed_output = if (is.na(x)) NA_character_ else proposed_of(s, x),
        shared_with_another_writer = !is.na(x) && x %in% cur_shared,
        stringsAsFactors = FALSE)
    }
  }

  gate_rows[[length(gate_rows) + 1L]] <- data.frame(
    result_family = fam,
    canonical_owner = owner,
    n_contributors = length(contribs),
    contributors = paste(contribs, collapse = " | "),
    n_current_outputs = length(cur_flat),
    n_proposed_destinations = length(unlist(proposals, use.names = FALSE)),
    current_file_level_collisions = length(cur_shared),
    proposed_destination_collisions = length(gate),
    gate = if (length(blockers)) paste("BLOCKED:", paste(blockers, collapse = "; ")) else "PASS",
    stringsAsFactors = FALSE)
}

d <- do.call(rbind, rows)
g <- do.call(rbind, gate_rows)

if (!dir.exists("audits")) dir.create("audits")

EMPTY <- data.frame(
  result_family = character(0), role = character(0), script = character(0),
  analysis_id = character(0), current_output = character(0),
  proposed_output = character(0), shared_with_another_writer = logical(0),
  stringsAsFactors = FALSE)
utils::write.csv(if (is.null(d)) EMPTY else d,
                 file.path("audits", paste0("phase6g_coordinating_families_", DOMAIN, ".csv")),
                 row.names = FALSE)

## Zero is the expected state for a migrated domain, not an error. The
## normalized layout keys each destination on its producing analysis, so a
## family shared by several writers cannot survive migration: each contributor
## ends up in its own tree and the coordination relationship moves into
## config/results_ownership.csv. Reporting that as a crash would make the
## successful outcome look like a tool failure.
if (is.null(g) || !nrow(g)) {
  cat("coordinating families touching ", DOMAIN, ": 0\n\n", sep = "")
  cat("No family in this domain is written by more than one analysis.\n")
  cat("For a migrated domain this is the expected end state: destinations are\n")
  cat("keyed on the producing analysis, so a shared directory cannot persist.\n")
  quit(save = "no", status = 0L)
}

cat("coordinating families touching", DOMAIN, ":", nrow(g), "\n\n")
for (i in seq_len(nrow(g))) {
  cat("=== ", g$result_family[i], "\n", sep = "")
  cat("  owner            : ", g$canonical_owner[i], "\n", sep = "")
  cat("  contributors (", g$n_contributors[i], ")  : ",
      gsub(" [|] ", "\n                       ", g$contributors[i]), "\n", sep = "")
  cat("  declared outputs in family : ", g$n_current_outputs[i], "\n", sep = "")
  cat("  proposed destinations      : ", g$n_proposed_destinations[i], "\n", sep = "")
  cat("  file-level collisions now  : ", g$current_file_level_collisions[i], "\n", sep = "")
  cat("  proposed collisions        : ", g$proposed_destination_collisions[i], "\n", sep = "")
  cat("  gate                       : ", g$gate[i], "\n\n", sep = "")
}

cat("total file-level collisions in current layout :",
    sum(g$current_file_level_collisions), "\n")
cat("total proposed destination collisions         :",
    sum(g$proposed_destination_collisions), "\n")
cat("families blocked                              :", sum(g$gate != "PASS"), "\n")

if (any(g$gate != "PASS")) {
  cat("\nSTOP: a coordinating family has an unresolved collision.\n")
  print(g[g$gate != "PASS", c("result_family", "gate")], row.names = FALSE)
  quit(status = 1)
}
cat("\nevery coordinating family is one owner plus contributors writing distinct artifacts\n")
