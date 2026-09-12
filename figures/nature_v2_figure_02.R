#!/usr/bin/env Rscript

# Nature-style v2 candidate Figure 02 entry point.
#
# Third parallel layer. Touches neither figures/figure_contract.yml (canonical)
# nor figures/figure_candidate_contract.yml (Part-16). Writes only under
# results/**/manuscript_candidates/nature_v2/**. Nothing is promoted.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_v2_figure_02.R")
args <- commandArgs(trailingOnly = TRUE)
res <- nv_build("figure_02", dry_run = "--dry-run" %in% args || is_dry_run())
if (!is.null(res)) {
  cat("\n===== Nature-style v2 Figure 02 =====\n")
  print(res$panels[, c("panel_id", "box_w_mm", "box_h_mm", "status")], row.names = FALSE)
  cat("\n")
  print(res$variants[, c("variant", "width_mm", "height_mm", "n_panels")], row.names = FALSE)
  cat("\nCandidate layer only. Canonical and Part-16 outputs untouched.\n")
}
