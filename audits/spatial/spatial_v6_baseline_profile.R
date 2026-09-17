#!/usr/bin/env Rscript

# Part-21 step 1: a PHENOTYPE-BLIND CON baseline spatial profile for every
# measured protein.
#
# WHY THIS EXISTS
# ---------------
# results/tables/11_spatial_systems/atlas/protein_baseline_spatial_profile.csv
# already holds CON-only means, but its ROW MEMBERSHIP is not phenotype-blind.
# 11_spatial_systems/build_protein_spatial_atlas.R sets its scope at
#
#   scope <- cand[cand$is_sus_res_fdr_supported | cand$is_wgcna_candidate, ]
#
# and 494 of its 505 rows (97.8%) entered through a candidate_reason containing
# the string "SUS - RES". Selecting display features from that table would
# therefore violate the rule that the Figure-2 spatial fingerprint must never be
# selected using stress phenotype.
#
# This script runs the IDENTICAL canonical code path (sps_build_spatial_levels
# -> level2 -> CON columns -> rowMeans per spatial unit) over the FULL measured
# protein universe instead of the phenotype-derived scope. The values are
# therefore the same quantity the atlas already publishes, on a superset of
# rows, and every displayed protein can then be chosen by an external
# prespecified rule alone.
#
# CLAIM DISCIPLINE: this is a descriptive transformation of canonical values.
# It fits no model, runs no test and computes no p-value. Level-2 is the
# canonical bilateral animal-level matrix; the only sample filter is
# StressGroup == "CON".

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "protigy_input_utils.R"))
source(repo_path("R", "spatial_systems_data_utils.R"))
source(repo_path("R", "spatial_atlas_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/spatial_v6_baseline_profile.R")

DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "spatial_v6")
  dir_create(d)
  file.path(d, ...)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  for (ds in DATASETS) message("[DRY-RUN] would build CON baseline for ", ds)
  message("[DRY-RUN] phenotype-blind: the only sample filter is StressGroup == CON")
  quit(save = "no", status = 0L)
}

# protein -> gene annotation only. This table is used purely as a lookup; no
# row is selected from it, so its candidate_tier column never influences what
# is displayed.
GENE <- repo_path("results", "tables", "10_biological_integration",
                  "wgcna_candidate_protein_shortlist", "global",
                  "wgcna_candidate_proteins_all.csv")
gmap <- as.data.frame(readr::read_csv(GENE, show_col_types = FALSE,
                                      progress = FALSE, guess_max = Inf))
gmap <- unique(gmap[, c("dataset", "ProteinGroupID", "GeneSymbol")])

long_rows <- list()
for (ds in DATASETS) {
  message("Building phenotype-blind CON baseline for ", ds)
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna",
                                   script = "figures/spatial_v6_baseline_profile.R",
                                   stage = "networks")
  md <- path_processed("01_preprocessing", "06_merged_metadata_module_score", ds,
                       "sample_metadata_merged_clean_for_module_scores.xlsx")
  canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                            dataset = ds, strict = TRUE)
  lv <- sps_build_spatial_levels(ds, canonical = canonical)

  con <- lv$level2$meta$StressGroup == "CON"
  if (!any(con)) stop("no CON samples for ", ds, call. = FALSE)
  if (any(lv$level2$meta$StressGroup[con] != "CON")) {
    stop("CON filter leaked a stressed animal for ", ds, call. = FALSE)
  }
  mat <- lv$level2$mat[, con, drop = FALSE]
  meta <- lv$level2$meta[con, , drop = FALSE]
  units <- sort(unique(as.character(meta$SpatialUnit)))

  prof <- vapply(units, function(u)
    rowMeans(mat[, meta$SpatialUnit == u, drop = FALSE], na.rm = TRUE),
    numeric(nrow(mat)))
  if (is.null(dim(prof))) prof <- matrix(prof, nrow = nrow(mat),
                                         dimnames = list(rownames(mat), units))

  # long format, with the canonical spatial key resolved through the shared
  # grammar so downstream panels never re-parse a unit string
  lg <- data.frame(
    dataset = ds,
    ProteinGroupID = rep(rownames(prof), times = ncol(prof)),
    raw_unit = rep(colnames(prof), each = nrow(prof)),
    con_mean_log2 = as.vector(prof),
    stringsAsFactors = FALSE)
  lg$spatial_unit <- sg_resolve_unit(lg$raw_unit, ds)

  # within-protein standardisation, so the panel shows SPATIAL PATTERN rather
  # than which proteins happen to be abundant. Purely descriptive.
  sp <- split(seq_len(nrow(lg)), lg$ProteinGroupID)
  lg$con_z <- NA_real_
  lg$protein_mean_log2 <- NA_real_
  lg$protein_sd_log2 <- NA_real_
  for (ix in sp) {
    v <- lg$con_mean_log2[ix]
    m <- mean(v, na.rm = TRUE)
    s <- stats::sd(v, na.rm = TRUE)
    lg$protein_mean_log2[ix] <- m
    lg$protein_sd_log2[ix] <- s
    lg$con_z[ix] <- if (is.finite(s) && s > 0) (v - m) / s else NA_real_
  }
  lg$n_con_animals <- length(unique(meta$AnimalID))
  long_rows[[ds]] <- lg
}

long <- do.call(rbind, long_rows)
long$GeneSymbol <- gmap$GeneSymbol[match(paste(long$dataset, long$ProteinGroupID),
                                         paste(gmap$dataset, gmap$ProteinGroupID))]
long$value_definition <- paste0(
  "con_mean_log2 = mean over CON animals of the canonical bilateral ",
  "animal-level (level-2) imputed log2 abundance in that spatial unit; ",
  "con_z = within-protein standardisation of con_mean_log2 across that ",
  "protein's spatial units. Descriptive only: no model, no test, no p-value.")
long$selection_status <- "phenotype_blind_full_universe"

write_csv_safe(long, OUT("spatial_v6_con_baseline_profile_long.csv"))

cat("\n===== phenotype-blind CON baseline spatial profile =====\n")
cat("rows:", nrow(long), "\n")
print(table(long$dataset))
cat("\nproteins per dataset:\n")
print(tapply(long$ProteinGroupID, long$dataset, function(z) length(unique(z))))
cat("\nspatial units per dataset:\n")
for (ds in DATASETS) {
  cat(sprintf("  %-16s %s\n", ds,
              paste(sort(unique(long$spatial_unit[long$dataset == ds])), collapse = " ")))
}
cat("\nCON animals:", paste(sort(unique(long$n_con_animals)), collapse = ","), "\n")
cat("gene symbols resolved:", sum(!is.na(long$GeneSymbol)), "of", nrow(long), "\n")

# reconcile against the canonical atlas table on the rows they share
ATL <- repo_path("results", "tables", "11_spatial_systems", "atlas",
                 "protein_baseline_spatial_profile.csv")
if (file.exists(ATL)) {
  a <- as.data.frame(readr::read_csv(ATL, show_col_types = FALSE,
                                     progress = FALSE, guess_max = Inf))
  uc <- setdiff(names(a), c("dataset", "ProteinGroupID", "n_con_animals"))
  al <- do.call(rbind, lapply(uc, function(u) data.frame(
    dataset = a$dataset, ProteinGroupID = a$ProteinGroupID,
    raw_unit = u, atlas_value = a[[u]], stringsAsFactors = FALSE)))
  al <- al[!is.na(al$atlas_value), , drop = FALSE]
  al$spatial_unit <- sg_resolve_unit(al$raw_unit, al$dataset)
  j <- merge(al, long[, c("dataset", "ProteinGroupID", "spatial_unit",
                          "con_mean_log2")],
             by = c("dataset", "ProteinGroupID", "spatial_unit"))
  d <- abs(j$atlas_value - j$con_mean_log2)
  cat("\n--- reconciliation against the canonical atlas baseline ---\n")
  cat("overlapping cells:", nrow(j), " max abs difference:",
      format(max(d, na.rm = TRUE), scientific = TRUE), "\n")
  if (nrow(j) && max(d, na.rm = TRUE) > 1e-8) {
    stop("phenotype-blind baseline does not reproduce the canonical atlas ",
         "values on the rows they share", call. = FALSE)
  }
  cat("PASS: identical to the canonical atlas values on every shared cell\n")
}
cat("\nwritten:", relative_to(OUT("spatial_v6_con_baseline_profile_long.csv")), "\n")
