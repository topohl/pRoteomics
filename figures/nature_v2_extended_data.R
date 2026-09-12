#!/usr/bin/env Rscript

# Nature-style v2 WGCNA Extended-Data candidate.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_v2_extended_data.R")
args <- commandArgs(trailingOnly = TRUE)
res <- nv_build("extended_data", dry_run = "--dry-run" %in% args || is_dry_run())
if (!is.null(res)) {
  cat("\n===== Nature-style v2 Extended Data (WGCNA) =====\n")
  print(res$panels[, c("panel_id", "box_w_mm", "box_h_mm", "status")], row.names = FALSE)
  print(res$variants[, c("variant", "n_panels")], row.names = FALSE)
}
