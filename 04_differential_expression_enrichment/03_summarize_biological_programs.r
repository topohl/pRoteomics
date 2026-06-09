#!/usr/bin/env Rscript

# ================================================================
# Script: 04_differential_expression_enrichment/03_summarize_biological_programs.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: compareGO and clusterProfiler manifests plus optional annotations.
# Produces: canonical biological program summary tables, source data, figures.
# Dataset behavior: delegates to the active dataset-aware implementation.
# Notes: preferred task-based alias; old command remains supported.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("04_differential_expression_enrichment", "03_biological_program_summary.r"))
