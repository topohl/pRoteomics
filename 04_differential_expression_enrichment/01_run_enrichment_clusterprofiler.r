#!/usr/bin/env Rscript

# ================================================================
# Script: 04_differential_expression_enrichment/01_run_enrichment_clusterprofiler.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: mapped contrast CSVs and optional enrichment configuration.
# Produces: canonical clusterProfiler processed/source tables, figures, logs.
# Dataset behavior: delegates to the active dataset-aware implementation.
# Notes: preferred task-based alias; old command remains supported.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("04_differential_expression_enrichment", "01_clusterProfiler.r"))
