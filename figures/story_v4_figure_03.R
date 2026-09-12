#!/usr/bin/env Rscript
# Story-v4 candidate Figure 03. Fifth parallel layer; canonical, Part-16,
# Part-17 and Part-18 outputs are untouched. Nothing is promoted.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
source(repo_path("R", "story_v3_figure_panels.R"))
source(repo_path("R", "story_v4_figure_utils.R"))
source(repo_path("R", "story_v4_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/story_v4_figure_03.R")
args <- commandArgs(trailingOnly = TRUE)
res <- s4e_build("figure_03", dry_run = "--dry-run" %in% args || is_dry_run())
if (!is.null(res)) {
  cat("\n===== Story-v4 Figure 03 =====\n")
  print(res$panels[, c("panel_id", "box_w_mm", "box_h_mm", "status")], row.names = FALSE)
  cat("\n"); print(res$variants[, c("variant", "n_panels", "largest_panel", "largest_area_share")], row.names = FALSE)
}
