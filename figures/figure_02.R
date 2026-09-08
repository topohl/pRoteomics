#!/usr/bin/env Rscript

# Explicit manuscript Figure 2 entry point. Scientific analyses are not fitted
# here; the script validates and materializes the exact declared panel assets.
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_figure_utils.R"))

manuscript_figure_main("02")
