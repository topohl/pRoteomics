#!/usr/bin/env Rscript

# Explicit behavioural Extended Data entry point. As with Figure 1, no statistic
# is computed here and none is computed anywhere in this repository: the
# behavioural analysis is frozen upstream and imported with hashes. This
# validates the declared panel assets and materialises them.
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_figure_utils.R"))

manuscript_figure_main("ED_behaviour")
