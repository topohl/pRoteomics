#!/usr/bin/env Rscript

# Build the stable publication source-data interface.
#
# This is the only surface the manuscript repository is allowed to read. It is
# organised by stable publication identity (figure_02, extended_data_01, ...)
# rather than by renderer generation, so a future reorganisation of the
# analysis output tree cannot break manuscript rendering.
#
#   results/publication_source_data/<publication_id>/...
#   results/publication_source_data/manifest.csv
#
# The bundle is an export, not a move: the canonical analysis outputs stay
# where the pre-restructure freeze recorded them, and each exported copy
# carries the source path plus a SHA-256 so the manuscript repository can
# verify it without knowing anything about this repository's internals.

source(file.path("R", "paths.R"))

SOURCE_ROOT <- path_results("source_data", "manuscript")
FIGURE_ROOT <- path_results("figures", "manuscript")
BUNDLE_ROOT <- path_results("publication_source_data")
CONTRACT_VERSION <- "publication_source_data_v1"

registry_path <- repo_path("manuscript", "canonical_publication_registry.csv")
if (!file.exists(registry_path)) {
  registry_path <- Sys.getenv("EXP9_PUBLICATION_REGISTRY", unset = "")
  if (!nzchar(registry_path) || !file.exists(registry_path)) {
    stop("Cannot locate canonical_publication_registry.csv. After the manuscript ",
         "extraction, point EXP9_PUBLICATION_REGISTRY at the copy in ",
         "Exp9_manuscript/provenance/publication_registry/.", call. = FALSE)
  }
}
registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE)

source_commit <- git_commit_sha()
if (is.na(source_commit)) source_commit <- "UNKNOWN"

canonical <- registry[registry$status == "CANONICAL", , drop = FALSE]
cat("canonical publication identities:", nrow(canonical), "\n")
cat("source commit                  :", source_commit, "\n\n")

rows <- list()
copied <- 0L
skipped <- character(0)

for (i in seq_len(nrow(canonical))) {
  pid <- canonical$publication_id[i]
  src_dir <- file.path(SOURCE_ROOT, pid)
  if (!dir.exists(src_dir)) { skipped <- c(skipped, pid); next }

  dst_dir <- file.path(BUNDLE_ROOT, pid)
  dir_create(dst_dir)

  files <- list.files(src_dir, recursive = TRUE, full.names = FALSE)
  for (f in files) {
    from <- file.path(src_dir, f)
    to <- file.path(dst_dir, f)
    dir_create(dirname(to))
    ok <- file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
    if (!ok) stop("Failed to copy ", from, call. = FALSE)
    copied <- copied + 1L

    dims <- c(NA_integer_, NA_integer_)
    if (grepl("[.]csv$", f, ignore.case = TRUE)) {
      d <- tryCatch(utils::read.csv(to, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) NULL)
      if (!is.null(d)) dims <- c(nrow(d), ncol(d))
    }

    rows[[length(rows) + 1L]] <- data.frame(
      publication_id = pid,
      source_repo = "topohl/pRoteomics",
      source_commit = source_commit,
      source_analysis = canonical$originating_analysis[i],
      source_table = relative_to(from),
      exported_file = relative_to(to),
      rows = dims[1],
      columns = dims[2],
      sha256 = file_hash_sha256(to),
      contract_version = CONTRACT_VERSION,
      stringsAsFactors = FALSE
    )
  }

  ## The assembled canonical artefact travels with its source data so the
  ## manuscript repository can prove byte equivalence of the rendered figure.
  assembled <- file.path(FIGURE_ROOT, pid, "assembled", paste0(pid, ".svg"))
  if (file.exists(assembled)) {
    to <- file.path(dst_dir, "assembled", paste0(pid, ".svg"))
    dir_create(dirname(to))
    file.copy(assembled, to, overwrite = TRUE, copy.date = TRUE)
    copied <- copied + 1L
    rows[[length(rows) + 1L]] <- data.frame(
      publication_id = pid,
      source_repo = "topohl/pRoteomics",
      source_commit = source_commit,
      source_analysis = canonical$originating_analysis[i],
      source_table = relative_to(assembled),
      exported_file = relative_to(to),
      rows = NA_integer_, columns = NA_integer_,
      sha256 = file_hash_sha256(to),
      contract_version = CONTRACT_VERSION,
      stringsAsFactors = FALSE
    )
  }
}

manifest <- do.call(rbind, rows)
dir_create(BUNDLE_ROOT)
utils::write.csv(manifest, file.path(BUNDLE_ROOT, "manifest.csv"), row.names = FALSE)

cat("files exported   :", copied, "\n")
cat("manifest rows    :", nrow(manifest), "\n")
cat("identities        :", length(unique(manifest$publication_id)), "\n")
if (length(skipped)) {
  cat("identities without source data:", paste(skipped, collapse = ", "), "\n")
}
cat("\nbundle:", relative_to(BUNDLE_ROOT), "\n")
cat("manifest:", relative_to(file.path(BUNDLE_ROOT, "manifest.csv")), "\n")
