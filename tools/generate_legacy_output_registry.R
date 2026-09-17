#!/usr/bin/env Rscript

# Generate config/legacy_output_registry.csv: the output roots that are read
# only.
#
# A root is legacy when it holds real artefacts that no registered writer
# produces any more. That is a measurement, not an opinion: the test is whether
# any path in pipeline.yml's produces lists falls beneath it.
#
# Nothing here moves or deletes anything. The point of the registry is that
# active code can be held to "do not write here" while the artefacts stay
# exactly where the frozen manifests already record them.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
suppressMessages(library(yaml))

BASELINE <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
declared <- split_paths(steps$produces)
literal_stem <- function(p) sub("(<|[*]).*$", "", p)
declared_stems <- literal_stem(declared)

fm <- system2("git", c("show", paste0(BASELINE, ":manuscript/prerestructure_freeze_manifest.csv")),
              stdout = TRUE)
frozen <- unique(utils::read.csv(text = paste(fm, collapse = "\n"),
                                 stringsAsFactors = FALSE)$repository_relative_path)

ns <- yaml::read_yaml(repo_path("config", "output_namespaces.yml"))

## Candidate roots: everything the v1 contract already calls historical,
## legacy, diagnostic or comparison, plus each top-level entry under results/
## that is not one of the analytical roots.
analytical <- unlist(ns$analytical_roots, use.names = FALSE)
candidates <- unique(c(
  unlist(ns$legacy_manuscript_authoring_roots),
  ns$manuscript_export_root,
  unlist(ns$manuscript_authoring_roots, use.names = FALSE),
  unlist(ns$diagnostic_roots),
  file.path("results", list.files("results"))
))
candidates <- candidates[!candidates %in% analytical]
candidates <- candidates[dir.exists(candidates) | file.exists(candidates)]
## the manuscript_candidates trees are named by a comparison marker, not listed
candidates <- unique(c(candidates,
  file.path(analytical, "manuscript_candidates")[dir.exists(file.path(analytical, "manuscript_candidates"))]))
candidates <- sort(unique(candidates))

count_files <- function(p) {
  if (!dir.exists(p)) return(if (file.exists(p)) 1L else 0L)
  length(list.files(p, recursive = TRUE, all.files = FALSE, no.. = TRUE))
}
size_of <- function(p) {
  if (!dir.exists(p)) return(if (file.exists(p)) unname(file.info(p)$size) else NA_real_)
  f <- list.files(p, recursive = TRUE, full.names = TRUE, no.. = TRUE)
  if (!length(f)) return(0)
  if (length(f) > 6000L) return(NA_real_)
  sum(file.info(f)$size, na.rm = TRUE)
}

rows <- lapply(candidates, function(p) {
  under <- paste0(p, "/")
  writers <- steps$script[vapply(seq_len(nrow(steps)), function(i)
    any(startsWith(split_paths(steps$produces[i]), under)) ||
    any(split_paths(steps$produces[i]) == p), logical(1))]
  writers <- unique(writers)
  fz <- sum(startsWith(frozen, under))
  reason <- if (length(writers)) {
    "active: a registered writer still declares output here"
  } else if (fz > 0L) {
    "manuscript rendering layer moved to Exp9_manuscript in Phase 6C; artefacts frozen"
  } else if (grepl("_superseded_|_failed_|COMPARISON|REPAIRED|_panels$|[.]zip$", p)) {
    "superseded, failed or comparison tree kept as provenance"
  } else {
    "no registered writer declares output here"
  }
  data.frame(
    legacy_path = p,
    policy = if (length(writers)) "ACTIVE_NOT_LEGACY" else "LEGACY_READ_ONLY",
    reason = reason,
    frozen_objects_beneath = fz,
    files = count_files(p),
    bytes = size_of(p),
    active_writers = paste(basename(writers), collapse = " | "),
    n_active_writers = length(writers),
    authoritative_copy = if (identical(p, "results/publication_source_data")) {
      ## Phase 6F copied this bundle byte-identically to the export boundary;
      ## the old path stays so the manuscript's recorded exported_file column
      ## remains resolvable.
      "exports/publication_source_data (byte-identical copy; see audits/phase6f_artifact_migration.csv)"
    } else if (fz > 0L && grepl("manuscript", p)) {
      "Exp9_manuscript (imported and hash-verified in provenance/source_manifests/)"
    } else {
      "this path"
    },
    stringsAsFactors = FALSE)
})
reg <- do.call(rbind, rows)
reg <- reg[order(reg$policy, -reg$frozen_objects_beneath, reg$legacy_path), ]
utils::write.csv(reg, repo_path("config", "legacy_output_registry.csv"), row.names = FALSE)

cat("candidate roots examined :", nrow(reg), "\n")
cat("LEGACY_READ_ONLY         :", sum(reg$policy == "LEGACY_READ_ONLY"), "\n")
cat("still active             :", sum(reg$policy == "ACTIVE_NOT_LEGACY"), "\n")
cat("frozen objects covered   :", sum(reg$frozen_objects_beneath), "\n")
cat("files under legacy roots :", sum(reg$files[reg$policy == "LEGACY_READ_ONLY"]), "\n")
cat("bytes under legacy roots :",
    format(sum(reg$bytes[reg$policy == "LEGACY_READ_ONLY"], na.rm = TRUE), big.mark = ","), "\n\n")
print(reg[, c("legacy_path", "policy", "frozen_objects_beneath", "files", "n_active_writers")],
      row.names = FALSE)
