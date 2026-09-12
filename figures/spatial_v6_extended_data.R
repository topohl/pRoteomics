#!/usr/bin/env Rscript
# Spatial-v6. Seventh parallel layer; canonical and Parts 16-20 are untouched.
# Every spatial panel uses the shared grammar in R/spatial_grammar_utils.R.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
source(repo_path("R", "story_v3_figure_panels.R"))
source(repo_path("R", "story_v4_figure_panels.R"))
source(repo_path("R", "story_v5_figure_panels.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "spatial_v6_figure_panels.R"))
source(repo_path("R", "spatial_v6_figure3_panels.R"))
source(repo_path("R", "spatial_v6_ed_panels.R"))
source(repo_path("R", "spatial_v6_wgcna_panels.R"))
source(repo_path("R", "spatial_v6_schematic.R"))
source(repo_path("R", "spatial_v6_figure_utils.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/spatial_v6_extended_data.R")
args <- commandArgs(trailingOnly = TRUE)
res <- s6e_build("extended_data", dry_run = "--dry-run" %in% args || is_dry_run())
if (!is.null(res)) {
  cat("\n===== Spatial-v6 extended_data =====\n")
  print(res$panels[, c("panel_id", "box_w_mm", "box_h_mm", "status")], row.names = FALSE)
  cat("\n"); print(res$variants[, c("variant", "n_panels", "largest_panel", "largest_area_share")], row.names = FALSE)
}
