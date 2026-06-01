#!/usr/bin/env Rscript

# Publication-facing alias for neuropil reference annotation of microglia ROI
# enrichment results. The older file name is retained for backward-compatible
# commands.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("04_differential_expression_enrichment", "04_neuropil_contamination_annotation.r"))
