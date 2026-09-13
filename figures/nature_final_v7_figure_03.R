#!/usr/bin/env Rscript
# Nature-final v7. Eighth parallel layer; canonical and Parts 16-22 untouched.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
source(repo_path("R", "story_v3_figure_panels.R"))
source(repo_path("R", "story_v4_figure_panels.R"))
source(repo_path("R", "story_v5_figure_panels.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "spatial_v6_figure3_panels.R"))
source(repo_path("R", "nature_final_v7_panels.R"))
source(repo_path("R", "nature_final_v7_figure3_panels.R"))
source(repo_path("R", "nature_final_v7_figure_utils.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork); library(scales) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_final_v7_figure_03.R")
args <- commandArgs(trailingOnly = TRUE)
res <- s7e_build("figure_03", dry_run = "--dry-run" %in% args || is_dry_run())
if (!is.null(res)) {
  cat("\n===== Nature-final v7 figure_03 =====\n")
  print(res$panels[, c("panel_id", "box_w_mm", "box_h_mm", "status")], row.names = FALSE)
  cat("\n"); print(res$variants[, c("variant", "n_panels", "largest_panel", "largest_area_share")], row.names = FALSE)
}
