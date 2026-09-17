#!/usr/bin/env Rscript

# PH-002 resolution, sections 7 and 8: the delta audit and the statistical
# immutability check.
#
#   Rscript 08_ph002_delta_audit.R before
#   ... edit source, rebuild the semantics layer ...
#   Rscript 08_ph002_delta_audit.R after
#
# SCOPE. The rebuild runs exactly two scripts, and between them they make 22
# write calls, all into the v9 manuscript-candidate tree. So the scope is split:
#
#   TIER A - everything those scripts could write: the v9 reports, tables,
#            source data and figures. Hashed AND compared value-wise, because a
#            regenerated CSV may legitimately differ in byte order while every
#            number stays identical.
#   TIER B - the canonical inputs the frozen contract names as primary_source,
#            plus the theme registry, the panel renderers and the two generator
#            scripts. Hashed only. Byte-identity is strictly stronger than value
#            identity, so an unchanged hash proves no statistic moved without
#            parsing 317 MB twice over a network share.
#
# results/manuscript (433,231 files, 13.6 GB - the PRIDE export package) and the
# full stage table trees are deliberately out of scope: the rebuild cannot reach
# them, and hashing them over this share does not terminate.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")
suppressMessages({library(digest); library(yaml)})

MODE <- commandArgs(trailingOnly = TRUE)[1]
if (!MODE %in% c("before", "after")) stop("usage: ... [before|after]")

SNAP_DIR <- file.path(PH_OUT, "ph002_snapshots")
dir.create(SNAP_DIR, recursive = TRUE, showWarnings = FALSE)
MC <- function(kind) file.path("results", kind, "manuscript_candidates",
                               "final_truth_v9")
rel <- function(f) gsub("\\\\", "/", sub("^[.]/", "", f))

tierA <- unlist(lapply(c("reports", "tables", "source_data", "figures"),
                       function(k) list.files(MC(k), recursive = TRUE,
                                              full.names = TRUE)))
CT <- yaml::read_yaml("figures/figure_final_truth_v9_contract.yml")
prim <- unique(unlist(lapply(CT$panels,
  function(p) as.character(p$primary_source %||% ""))))
tierB <- c(prim[nzchar(prim)],
           "config/manuscript_go_theme_registry.tsv",
           "R/final_truth_v9_panels.R",
           "figures/figure_final_truth_v9_contract.yml",
           "figures/final_truth_v9_vector_audit.R",
           "figures/final_truth_v9_semantics.R")
tierA <- tierA[file.exists(tierA) & !dir.exists(tierA)]
tierB <- tierB[file.exists(tierB) & !dir.exists(tierB)]

hash_of <- function(v, tier) if (!length(v)) NULL else data.frame(
  file = vapply(v, rel, character(1)), tier = tier,
  hash = vapply(v, function(f) digest::digest(file = f, algo = "sha256"),
                character(1)),
  bytes = file.size(v), stringsAsFactors = FALSE)
files <- rbind(hash_of(tierA, "A"), hash_of(tierB, "B"))
rownames(files) <- NULL

# Value-level signature for TIER A only. Fixed precision, so a reformatted
# number is not reported as a changed number, but a changed number always is.
digest_numeric <- function(f) {
  d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  if (!length(num)) return(NULL)
  do.call(rbind, lapply(num, function(cl) data.frame(
    file = rel(f), column = cl, n_values = sum(!is.na(d[[cl]])),
    value_hash = digest::digest(sprintf("%.12g", d[[cl]]), algo = "sha256"),
    stringsAsFactors = FALSE)))
}
nums <- do.call(rbind, lapply(grep("[.]csv$", tierA, value = TRUE),
                              digest_numeric))
rownames(nums) <- NULL

saveRDS(list(files = files, nums = nums),
        file.path(SNAP_DIR, paste0(MODE, ".rds")))
cat(sprintf("[%s] tier A %d files | tier B %d files | %d numeric columns, %s values\n",
            MODE, sum(files$tier == "A"), sum(files$tier == "B"), nrow(nums),
            format(sum(nums$n_values), big.mark = ",")))

# ================================================================ comparison
if (identical(MODE, "after")) {
  b <- readRDS(file.path(SNAP_DIR, "before.rds"))
  B <- b$files; A <- files
  all_f <- union(B$file, A$file)
  d <- data.frame(
    file = all_f, tier = A$tier[match(all_f, A$file)],
    before_hash = B$hash[match(all_f, B$file)],
    after_hash = A$hash[match(all_f, A$file)], stringsAsFactors = FALSE)
  d$tier[is.na(d$tier)] <- B$tier[match(d$file[is.na(d$tier)], B$file)]
  d$changed <- is.na(d$before_hash) | is.na(d$after_hash) |
               d$before_hash != d$after_hash

  classify <- function(f, changed) {
    if (!changed) return("UNCHANGED")
    if (grepl("^figures/final_truth_v9_.*[.]R$", f)) return("TEXT_EXPECTED")
    if (grepl("/reports/.*[.]md$", f)) return("TEXT_EXPECTED")
    if (grepl("[.](svg|pdf|png)$", f)) return("PANEL_GEOMETRY_CHANGED_UNEXPECTEDLY")
    if (grepl("manifest|checksum", f)) return("MANIFEST_METADATA_EXPECTED")
    if (grepl("^tests/", f)) return("TEST_EXPECTED")
    if (grepl("/source_data/", f)) return("SOURCE_DATA_VALUE_CHANGED")
    if (grepl("/tables/", f)) return("STATISTIC_CHANGED")
    if (grepl("^config/|^R/|contract[.]yml$", f)) return("ARCHITECTURE_CHANGED")
    "UNRELATED_FILE_CHANGED"
  }
  d$change_class <- mapply(classify, d$file, d$changed)

  # A table whose bytes moved but whose every number is identical is a
  # formatting change, not a statistical one. Reclassify explicitly, never
  # silently, and only on evidence from the value-level comparison.
  NB <- b$nums; NA_ <- nums
  key <- function(x) paste(x$file, x$column, sep = "::")
  kb <- key(NB); ka <- key(NA_)
  shared <- intersect(kb, ka)
  moved <- shared[NB$value_hash[match(shared, kb)] !=
                  NA_$value_hash[match(shared, ka)]]
  n_compared <- sum(NB$n_values[match(shared, kb)])
  n_changed_vals <- if (length(moved))
    sum(NA_$n_values[match(moved, ka)]) else 0L
  stat_files <- unique(sub("::.*", "", moved))

  reclass <- d$change_class %in% c("STATISTIC_CHANGED",
                                   "SOURCE_DATA_VALUE_CHANGED") &
             !(d$file %in% stat_files)
  d$change_class[reclass] <- "MANIFEST_METADATA_EXPECTED"
  OK <- c("UNCHANGED", "TEXT_EXPECTED", "RENDER_TEXT_EXPECTED",
          "MANIFEST_METADATA_EXPECTED", "TEST_EXPECTED")
  d$expected <- d$change_class %in% OK
  d$notes <- ifelse(is.na(d$before_hash), "new file",
             ifelse(is.na(d$after_hash), "file removed",
             ifelse(reclass, "bytes differ, every numeric value identical",
             ifelse(d$changed, "content differs", ""))))
  d <- d[order(!d$changed, d$file), ]
  utils::write.csv(d, file.path(PH_OUT, "ph002_delta_audit.csv"),
                   row.names = FALSE)

  cat("\n===== PH-002 DELTA AUDIT =====\n")
  cat("files in scope:", nrow(d), "| changed:", sum(d$changed), "\n")
  print(table(d$change_class[d$changed]))
  if (any(d$changed)) print(d[d$changed, c("file", "change_class")],
                            row.names = FALSE)

  cat("\n===== STATISTICAL IMMUTABILITY =====\n")
  cat("tier B canonical statistic files hashed:", sum(d$tier == "B"),
      "| changed:", sum(d$changed & d$tier == "B"), "\n")
  cat("tier A numeric columns compared:", length(shared),
      "| numeric VALUES compared:", format(n_compared, big.mark = ","), "\n")
  cat("numeric columns changed:", length(moved),
      "| numeric VALUES changed:", n_changed_vals, "\n")
  cat("columns only-before:", length(setdiff(kb, ka)),
      "| only-after:", length(setdiff(ka, kb)), "\n")
  if (length(moved)) { cat("CHANGED COLUMNS:\n"); print(moved) }

  bad <- d[d$changed & !d$expected, , drop = FALSE]
  if (nrow(bad)) {
    cat("\nFORBIDDEN CHANGES:\n"); print(bad[, c("file", "change_class")])
    stop("PH-002 delta audit: forbidden change class present", call. = FALSE)
  }
  cat("\nresult: every change is in an allowed class; 0 statistical values moved\n")
}
