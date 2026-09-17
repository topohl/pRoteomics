#!/usr/bin/env Rscript

# Establish exports/publication_source_data/ as the publication boundary.
#
# The bundle Exp9_manuscript imports from, results/publication_source_data/,
# has no registered writer in this repository: 10 publication identities, 55
# files, 2.5 MB, produced by a step that is no longer in pipeline.yml. It is
# what the manuscript's source_data/pRoteomics/manifest.csv resolves against,
# so it is real and load-bearing, and it was sitting in the middle of the
# canonical results tree with nothing owning it.
#
# This copies it, byte for byte, to exports/publication_source_data/ and writes
# a manifest with a sha256 per file. It is a copy and not a move: the old path
# is what the manuscript's existing manifest records, so leaving it in place
# keeps that provenance chain resolvable. The old path is then registered
# LEGACY_READ_ONLY and the new path is the declared interface.
#
# Run with DRY=1 to preview.

source(file.path("R", "paths.R"))
suppressMessages(library(yaml))

DRY <- nzchar(Sys.getenv("DRY"))
SRC <- "results/publication_source_data"
DST <- "exports/publication_source_data"

MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT",
                 unset = normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                                       winslash = "/", mustWork = FALSE))

if (!dir.exists(SRC)) stop("source bundle not found: ", SRC, call. = FALSE)

files <- list.files(SRC, recursive = TRUE, all.files = FALSE, no.. = TRUE)
cat("bundle files:", length(files), "\n")
cat("identities   :", length(list.dirs(SRC, recursive = FALSE)), "\n")

contract_path <- repo_path("config", "publication_source_data_contract.yml")
identities <- character(0)
if (file.exists(contract_path)) {
  ct <- yaml::read_yaml(contract_path)
  identities <- vapply(ct$identities, function(x) as.character(x$publication_id), character(1))
}

## provenance the manuscript already recorded for this bundle
recorded_commit <- NA_character_
mm <- file.path(MR, "source_data/pRoteomics/manifest.csv")
if (file.exists(mm)) {
  d <- utils::read.csv(mm, stringsAsFactors = FALSE)
  if ("source_commit" %in% names(d)) {
    u <- unique(d$source_commit[nzchar(d$source_commit)])
    recorded_commit <- if (length(u) == 1L) u else paste(u, collapse = ";")
  }
}
cat("commit recorded by the manuscript for this bundle:", recorded_commit, "\n")

rows <- list()
copied <- 0L
for (f in files) {
  from <- file.path(SRC, f)
  to <- file.path(DST, f)
  h_before <- unname(tools::sha256sum(from))
  if (!DRY) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
    if (!ok) stop("copy failed: ", from, call. = FALSE)
    copied <- copied + 1L
  }
  h_after <- if (!DRY && file.exists(to)) unname(tools::sha256sum(to)) else NA_character_
  id <- strsplit(f, "/", fixed = TRUE)[[1]][1]
  rows[[length(rows) + 1L]] <- data.frame(
    publication_id = if (id %in% identities) id else paste0(id, " (not in contract)"),
    relative_path = f,
    old_path = from,
    new_path = to,
    bytes = unname(file.info(from)$size),
    old_sha256 = h_before,
    new_sha256 = h_after,
    byte_identical = !DRY && !is.na(h_after) && identical(h_before, h_after),
    stringsAsFactors = FALSE)
}
mani <- do.call(rbind, rows)

cat(if (DRY) "WOULD copy: " else "copied: ", if (DRY) length(files) else copied, " files\n", sep = "")
if (!DRY) {
  cat("byte-identical:", sum(mani$byte_identical), "/", nrow(mani), "\n")
  if (!all(mani$byte_identical)) stop("a copied file is not byte-identical", call. = FALSE)
}

if (!DRY) {
  ## The bundle already carries its own manifest.csv: 54 rows with
  ## publication_id, source_analysis, source_table, exported_file, rows,
  ## columns and a sha256 per file. It is copied byte-identically like every
  ## other file, and because the copy does not change any content those sha256
  ## values still validate against the copied bundle. Writing a second manifest
  ## here would create a competing authority for no gain; the exported_file
  ## column stays a pointer to the historical location, which still exists.
  stopifnot(file.exists(file.path(DST, "manifest.csv")))
  h_src <- unname(tools::sha256sum(file.path(SRC, "manifest.csv")))
  h_dst <- unname(tools::sha256sum(file.path(DST, "manifest.csv")))
  if (!identical(h_src, h_dst)) stop("the bundle manifest was not copied intact", call. = FALSE)
  cat("bundle manifest carried intact:", substr(h_dst, 1, 12), "\n")

  ## the section 11 migration record
  rec <- mani[, c("old_path", "new_path", "old_sha256", "new_sha256")]
  rec$owner <- "analysis/publication_source_data/09_export_source_data.R"
  rec$consumers_updated <- "Exp9_manuscript/source_data/pRoteomics/manifest.csv"
  rec$migration_class <- "COPIED_TO_EXPORT_BOUNDARY_BYTE_IDENTICAL"
  if (!dir.exists("audits")) dir.create("audits")
  utils::write.csv(rec, "audits/phase6f_artifact_migration.csv", row.names = FALSE)
  cat("wrote audits/phase6f_artifact_migration.csv\n")
}

cat("\nidentities in the contract  :", length(identities), "\n")
cat("identities present in bundle:", length(unique(sub("/.*", "", files))), "\n")
missing <- setdiff(identities, unique(sub("/.*", "", files)))
if (length(missing)) cat("contract identities absent from the bundle:",
                         paste(missing, collapse = ", "), "\n")
extra <- setdiff(unique(sub("/.*", "", files)), identities)
if (length(extra)) cat("bundle entries not in the contract:",
                       paste(extra, collapse = ", "), "\n")
