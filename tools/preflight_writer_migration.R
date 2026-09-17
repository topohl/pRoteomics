#!/usr/bin/env Rscript

# Phase 6G preflight: everything that must be known about a domain before any
# of its writers is touched.
#
#   Rscript tools/preflight_writer_migration.R --domain <domain>
#   -> audits/phase6g_preflight_<domain>.csv
#      audits/consumer_inventory/<analysis_id>.csv   (one per writer)
#
# The order is deliberate and comes from the enrichment pilot, where consumers
# were discovered by hand and one was missed:
#
#   enumerate writers
#   -> enumerate consumers mechanically
#   -> classify every dependency
#   -> freeze the pre-migration inventory
#   -> only then migrate
#
# The gate at the end is what makes that an order rather than a suggestion.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  i <- which(args == flag)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
DOMAIN <- arg_value("--domain")
if (!nzchar(DOMAIN)) stop("usage: --domain <domain, e.g. spatial_networks>", call. = FALSE)

BASELINE <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]
dom_steps <- steps[startsWith(steps$script, paste0("analysis/", DOMAIN, "/")), , drop = FALSE]
if (!nrow(dom_steps)) stop("no registered writers in domain '", DOMAIN, "'", call. = FALSE)

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
literal_stem <- function(p) sub("(<|[*]).*$", "", p)

## frozen baseline objects, to know which writers touch protected ground
fm <- suppressWarnings(system2("git",
  c("show", paste0(BASELINE, ":manuscript/prerestructure_freeze_manifest.csv")),
  stdout = TRUE, stderr = FALSE))
frozen <- if (length(fm)) {
  unique(utils::read.csv(text = paste(fm, collapse = "\n"),
                         stringsAsFactors = FALSE)$repository_relative_path)
} else character(0)

## Phase 6F recorded the historical path and the proposed normalized path for
## every declared output. Rows are keyed on the pre-rename entrypoint path, so
## match on the basename.
inv6f <- repo_path("audits", "phase6f_output_inventory.csv")
d6 <- if (file.exists(inv6f)) utils::read.csv(inv6f, stringsAsFactors = FALSE) else NULL

## Phase 6E ownership: who is the canonical owner of each result family
own <- repo_path("config", "results_ownership.csv")
ownership <- if (file.exists(own)) utils::read.csv(own, stringsAsFactors = FALSE) else NULL

cat("domain              :", DOMAIN, "\n")
cat("registered writers  :", nrow(dom_steps), "\n\n")

rows <- list()
for (i in seq_len(nrow(dom_steps))) {
  f <- dom_steps$script[i]
  aid <- sub("[.][Rr]$", "", basename(f))

  declared_now <- split_paths(dom_steps$produces[i])
  hist_paths <- character(0)
  proposed <- character(0)
  if (!is.null(d6)) {
    r6 <- d6[basename(d6$producing_analysis) == basename(f), , drop = FALSE]
    hist_paths <- unique(r6$current_path[nzchar(r6$current_path)])
    proposed <- unique(r6$proposed_path[nzchar(r6$proposed_path)])
  }

  frozen_beneath <- sum(vapply(frozen, function(x) {
    any(startsWith(x, literal_stem(hist_paths)[nchar(literal_stem(hist_paths)) > 12]))
  }, logical(1)))

  ## run the consumer enumerator for this writer
  cat("enumerating consumers for", aid, "... ")
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c("tools/enumerate_output_consumers.R", "--analysis-id", aid, "--quiet"),
    stdout = TRUE, stderr = TRUE))
  inv_file <- repo_path("audits", "consumer_inventory", paste0(aid, ".csv"))
  n_runtime <- NA_integer_
  n_unknown <- NA_integer_
  kinds <- ""
  if (file.exists(inv_file)) {
    ci <- utils::read.csv(inv_file, stringsAsFactors = FALSE)
    n_runtime <- length(unique(ci$consumer_file[ci$runtime_relevant]))
    n_unknown <- sum(ci$dependency_kind == "UNKNOWN")
    kinds <- paste(sort(unique(ci$dependency_kind[ci$runtime_relevant])), collapse = " | ")
    cat(n_runtime, "runtime consumer files\n")
  } else {
    cat("FAILED\n")
    if (length(st)) cat(paste0("    ", utils::head(st, 3), collapse = "\n"), "\n")
  }

  own_rows <- if (!is.null(ownership)) {
    ownership[grepl(paste0("/", aid, "[.][Rr]$"), ownership$canonical_owner) |
              grepl(aid, ownership$contributing_analyses, fixed = TRUE), , drop = FALSE]
  } else NULL
  multi <- if (!is.null(own_rows) && nrow(own_rows)) {
    sum(own_rows$n_contributors > 1L)
  } else 0L
  is_owner <- if (!is.null(own_rows) && nrow(own_rows)) {
    any(grepl(paste0("/", aid, "[.][Rr]$"), own_rows$canonical_owner))
  } else NA

  rows[[length(rows) + 1L]] <- data.frame(
    domain = DOMAIN,
    analysis_id = aid,
    entrypoint = f,
    n_declared_outputs = length(declared_now),
    current_output_roots = paste(utils::head(unique(dirname(sub("/$", "", declared_now))), 3), collapse = " | "),
    historical_output_roots = paste(utils::head(hist_paths, 3), collapse = " | "),
    proposed_output_roots = paste(utils::head(proposed, 3), collapse = " | "),
    frozen_objects_beneath = frozen_beneath,
    runtime_consumer_files = n_runtime,
    unknown_dependency_kinds = n_unknown,
    consumer_kinds = kinds,
    canonical_owner_of_a_family = is_owner,
    multi_writer_families = multi,
    stringsAsFactors = FALSE)
}
pf <- do.call(rbind, rows)

# ------------------------------------------------- collision gate (section 12)
## Two canonical writers resolving to the same destination is the failure the
## ownership registry exists to prevent. It is checked before migration, not
## discovered during it.
prop <- lapply(seq_len(nrow(pf)), function(i) {
  p <- trimws(unlist(strsplit(pf$proposed_output_roots[i], "|", fixed = TRUE)))
  p[nzchar(p)]
})
names(prop) <- pf$analysis_id
flat <- unlist(prop)
dups <- unique(flat[duplicated(flat)])
collisions <- list()
for (d in dups) {
  who <- names(prop)[vapply(prop, function(p) d %in% p, logical(1))]
  if (length(who) < 2L) next
  ## allowed only when the ownership registry names one owner for that family
  fam <- sub("/[^/]+/?$", "", d)
  reg_owner <- if (!is.null(ownership)) {
    ownership$canonical_owner[ownership$result_family == fam]
  } else character(0)
  collisions[[d]] <- list(writers = who, declared_owner = reg_owner)
}

out <- repo_path("audits", paste0("phase6g_preflight_", DOMAIN, ".csv"))
utils::write.csv(pf, out, row.names = FALSE)

cat("\n=== preflight summary ===\n")
cat("writers                      :", nrow(pf), "\n")
cat("declared outputs             :", sum(pf$n_declared_outputs), "\n")
cat("runtime consumer files (max) :", max(pf$runtime_consumer_files, na.rm = TRUE), "\n")
cat("runtime consumer files (sum) :", sum(pf$runtime_consumer_files, na.rm = TRUE), "\n")
cat("unknown dependency kinds     :", sum(pf$unknown_dependency_kinds, na.rm = TRUE), "\n")
cat("frozen objects beneath       :", sum(pf$frozen_objects_beneath), "\n")
cat("multi-writer families        :", sum(pf$multi_writer_families), "\n")
cat("destination collisions       :", length(collisions), "\n")
if (length(collisions)) {
  for (d in names(collisions)) {
    cat("  ", d, "\n     writers:", paste(collisions[[d]]$writers, collapse = ", "),
        "\n     declared owner:", paste(collisions[[d]]$declared_owner, collapse = ", "), "\n")
  }
}

## ------------------------------------------------------- gate (section 8)
blockers <- character(0)
if (sum(pf$unknown_dependency_kinds, na.rm = TRUE) > 0L) {
  blockers <- c(blockers, "unknown dependency kinds in a consumer inventory")
}
if (any(is.na(pf$runtime_consumer_files))) {
  blockers <- c(blockers, "a writer's consumer inventory could not be produced")
}
unresolved_collisions <- Filter(function(x) !length(x$declared_owner), collisions)
if (length(unresolved_collisions)) {
  blockers <- c(blockers, "two canonical writers share a destination with no declared owner")
}

risk <- if (length(blockers)) "BLOCKED" else if (sum(pf$multi_writer_families) > 0L) {
  "ELEVATED: a coordinating family is involved; section 16 applies"
} else if (sum(pf$frozen_objects_beneath) > 0L) {
  "ELEVATED: frozen baseline objects sit beneath a historical root"
} else {
  "CONTAINED"
}
cat("\nmigration risk               :", risk, "\n")
if (length(blockers)) {
  cat("\nBLOCKED. Resolve before migrating:\n")
  cat(paste0("  - ", blockers, collapse = "\n"), "\n")
} else {
  cat("gate                         : PASS - preflight complete, migration may proceed\n")
}
cat("\nwrote", relative_to(out), "\n")
