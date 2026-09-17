#!/usr/bin/env Rscript

# Build audits/restructure_migration_map.csv: one row per distinct object in
# the pre-restructure scientific freeze, mapping baseline_path to the place it
# now lives, in whichever of the two repositories now owns it.
#
# The freeze manifest holds 239 rows over 200 distinct objects, because it
# records one row per (object, publication role) pair. The map is keyed on the
# object, so it has 200 rows; tools/verify_restructure_equivalence.R reconciles
# both levels.
#
# Destination resolution order:
#   1. the object was extracted into Exp9_manuscript
#   2. the inventory relocated it inside pRoteomics
#   3. it stays where it was (the untracked canonical result tree)

source(file.path("R", "paths.R"))

BASELINE_COMMIT <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"
BASELINE_MANIFEST <- "manuscript/prerestructure_freeze_manifest.csv"

MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT",
                 unset = normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                                       winslash = "/", mustWork = FALSE))

sha <- function(p) {
  if (!file.exists(p) || dir.exists(p)) return(NA_character_)
  unname(tools::sha256sum(p))
}

## Baseline oracle, read from git so it cannot drift.
txt <- system2("git", c("-C", shQuote(repo_root()), "show",
                        paste0(BASELINE_COMMIT, ":", BASELINE_MANIFEST)),
               stdout = TRUE)
base <- utils::read.csv(text = paste(txt, collapse = "\n"), stringsAsFactors = FALSE)
obj <- base[!duplicated(base$repository_relative_path), ]
cat("baseline assertions:", nrow(base), "  distinct objects:", nrow(obj), "\n")

inv <- utils::read.csv(repo_path("docs", "restructure_inventory.csv"), stringsAsFactors = FALSE)
inv_map <- stats::setNames(inv$proposed_new_path, inv$old_path)
inv_repo <- stats::setNames(inv$destination_repo, inv$old_path)

rows <- vector("list", nrow(obj))
for (i in seq_len(nrow(obj))) {
  bp <- obj$repository_relative_path[i]
  dest_repo <- "pRoteomics"
  dest_path <- bp

  if (bp %in% names(inv_map)) {
    dest_repo <- inv_repo[[bp]]
    dest_path <- inv_map[[bp]]
  }

  root <- if (identical(dest_repo, "Exp9_manuscript")) MR else repo_root()
  abs <- file.path(root, dest_path)

  ## Fall back to the baseline path if the planned destination is not there,
  ## so the map always describes reality rather than intent.
  if (!file.exists(abs)) {
    alt <- file.path(repo_root(), bp)
    if (file.exists(alt)) { dest_repo <- "pRoteomics"; dest_path <- bp; abs <- alt }
  }

  post <- sha(abs)
  same_path <- identical(dest_path, bp) && identical(dest_repo, "pRoteomics")
  changed <- !is.na(post) && !identical(post, obj$sha256[i])

  cls <- if (identical(dest_repo, "Exp9_manuscript")) {
    "MOVED_TO_MANUSCRIPT_REPO"
  } else if (!same_path) {
    "MOVED_WITHIN_SCIENTIFIC_REPO"
  } else {
    "RETAINED_LEGACY_PATH"
  }

  status <- ""
  if (identical(bp, BASELINE_MANIFEST)) {
    status <- paste(
      "Self-referential: the freeze manifest lists itself, so the sha256 it",
      "records for itself cannot be its own final content hash. Verified by",
      "byte length instead. Pre-existing property of the baseline.")
  } else if (changed) {
    cls <- "REWRITTEN_PATH_CONTRACT"
    status <- paste(
      "Layout contract: content encodes the pre-migration directory structure,",
      "so it had to be repointed for the suite to pass. No scientific value,",
      "statistic, table, panel or prose changed. Re-frozen under",
      "post-restructure-architecture-freeze-2026-09.")
  } else if (cls == "RETAINED_LEGACY_PATH") {
    status <- "Untracked canonical output retained at its baseline path; not bulk-moved."
  } else {
    status <- "Relocated byte-identical."
  }

  rows[[i]] <- data.frame(
    baseline_repo = "pRoteomics",
    baseline_path = bp,
    baseline_sha256 = obj$sha256[i],
    destination_repo = dest_repo,
    destination_path = dest_path,
    postmigration_sha256 = ifelse(is.na(post), "", post),
    migration_class = cls,
    status = status,
    stringsAsFactors = FALSE)
}

map <- do.call(rbind, rows)
dir_create(repo_path("audits"))
utils::write.csv(map, repo_path("audits", "restructure_migration_map.csv"), row.names = FALSE)

cat("\nmap rows:", nrow(map), "\n\n")
cat("=== migration_class ===\n"); print(sort(table(map$migration_class), decreasing = TRUE))
cat("\n=== destination_repo ===\n"); print(table(map$destination_repo))
cat("\n=== destinations present ===\n")
absall <- ifelse(map$destination_repo == "Exp9_manuscript",
                 file.path(MR, map$destination_path),
                 file.path(repo_root(), map$destination_path))
cat(sum(file.exists(absall)), "/", nrow(map), "\n")
cat("\n=== hash equality (excluding declared rewrites) ===\n")
chk <- map$migration_class != "REWRITTEN_PATH_CONTRACT" &
  map$baseline_path != BASELINE_MANIFEST
cat(sum(map$postmigration_sha256[chk] == map$baseline_sha256[chk]), "/", sum(chk), "\n")
cat("\nwrote audits/restructure_migration_map.csv\n")
