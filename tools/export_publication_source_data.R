#!/usr/bin/env Rscript

# Build the stable publication source-data interface.
#
# This is the only surface the manuscript repository is allowed to read. It is
# organised by stable publication identity (figure_02, extended_data_01, ...)
# rather than by renderer generation, so a future reorganisation of the
# analysis output tree cannot break manuscript rendering.
#
#   exports/publication_source_data/<publication_id>/...
#   exports/publication_source_data/manifest.csv
#
# The bundle is an export, not a move: the canonical analysis outputs stay
# where the pre-restructure freeze recorded them, and each exported copy
# carries the source path plus a SHA-256 so the manuscript repository can
# verify it without knowing anything about this repository's internals.
#
# Phase 6F moved the bundle out of results/ and into the exports/ lifecycle.
# It was sitting inside the canonical results tree, which made it look like an
# analysis result rather than a frozen interface, and no registered pipeline
# step declared it because this tool is invoked by hand rather than by
# pipeline.yml. exports/ is now the only tree the manuscript repository reads.
#
# The two source roots below are registered LEGACY_READ_ONLY in
# config/legacy_output_registry.csv: reading them is allowed and is what this
# tool does, but nothing may write there any more. The assembled figures under
# results/figures/manuscript are produced by renderers that moved to
# Exp9_manuscript in Phase 6C, which is recorded as remaining debt in
# docs/OUTPUT_LAYOUT.md.

source(file.path("R", "paths.R"))

SOURCE_ROOT <- path_results("source_data", "manuscript")
FIGURE_ROOT <- path_results("figures", "manuscript")
BUNDLE_ROOT <- path_export("publication_source_data")

# The list of publication identities is a scientific-side contract, so that
# building the bundle never requires reading the manuscript repository.
contract_path <- repo_path("config", "publication_source_data_contract.yml")
if (!file.exists(contract_path)) {
  stop("Missing config/publication_source_data_contract.yml", call. = FALSE)
}
contract <- yaml::read_yaml(contract_path)
CONTRACT_VERSION <- as.character(contract$contract_version)
canonical <- do.call(rbind, lapply(contract$identities, function(x) data.frame(
  publication_id = as.character(x$publication_id),
  canonical_source_data = as.character(x$source_data_dir),
  originating_analysis = as.character(x$originating_analysis),
  stringsAsFactors = FALSE)))

source_commit <- git_commit_sha()
if (is.na(source_commit)) source_commit <- "UNKNOWN"

cat("canonical publication identities:", nrow(canonical), "\n")
cat("source commit                  :", source_commit, "\n\n")

rows <- list()
copied <- 0L
skipped <- character(0)

for (i in seq_len(nrow(canonical))) {
  pid <- canonical$publication_id[i]

  ## The contract declares a source_data_dir per identity, and this loop used
  ## to ignore it and rebuild SOURCE_ROOT/<pid> instead. For every identity
  ## that existed then the two were the same path, so the field was inert and
  ## the assumption held silently.
  ##
  ## It stops holding as soon as an identity's source data lives anywhere else,
  ## and it cannot be satisfied by staging a copy under SOURCE_ROOT: that root
  ## is registered LEGACY_READ_ONLY in config/legacy_output_registry.csv,
  ## precisely because no registered writer declares output there. Honouring
  ## the declared path lets an identity be exported straight from the canonical
  ## results tree of the analysis that owns it, with no write into the legacy
  ## tree and no second copy to keep in step.
  ##
  ## The declared path wins; SOURCE_ROOT/<pid> remains the fallback so an
  ## identity whose contract entry predates this change resolves as before.
  src_dir <- repo_path(canonical$canonical_source_data[i])
  if (!dir.exists(src_dir)) src_dir <- file.path(SOURCE_ROOT, pid)
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
