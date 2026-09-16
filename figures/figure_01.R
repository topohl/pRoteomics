#!/usr/bin/env Rscript

# Explicit manuscript Figure 1 entry point. No statistic is computed here and
# none is computed anywhere in this repository: the behavioural analysis is
# frozen upstream and imported with hashes. This validates the declared panel
# assets and materialises them.
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_figure_utils.R"))

manuscript_figure_main("01")
