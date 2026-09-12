#!/usr/bin/env Rscript

# Candidate Figure 3 entry point.
#
# THIS IS NOT figures/figure_03.R AND DOES NOT TOUCH IT. The canonical
# contract, entry points and outputs are untouched by this script; everything
# written here lands under results/**/manuscript_candidates/**. Nothing is
# promoted to the manuscript by running it.
source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "candidate_figure_utils.R"))
source(repo_path("R", "candidate_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/candidate_figure_03.R")
a <- cf_args()
res <- cf_build_figure("03", check_only = a$check_only, dry_run = a$dry_run)
if (!a$check_only && !a$dry_run) {
  cat("\n===== Candidate Figure 3 =====\n")
  print(res$panels[, c("candidate_panel_id", "render_mode", "status")], row.names = FALSE)
  cat("\nassemblies:\n")
  print(res$assemblies[, c("assembly", "n_panels", "svg")], row.names = FALSE)
  cat("\nCANDIDATE LAYER ONLY. Canonical Figure 3 is unchanged.\n")
}
