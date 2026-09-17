#!/usr/bin/env Rscript

# Explicit manuscript entry point for the promoted proteomics Extended Data.
#
# Four of the seven final_truth_v9 Extended Data figures were promoted in Phase
# 6A after a per-figure audit against the current frozen contracts. This layer
# validates the declared panel assets, hashes them, republishes their source data
# into the manuscript namespace and assembles - exactly as it does for Figures 1
# to 3. It computes nothing.
#
# The WGCNA figure and ED7 are deliberately absent, as is the historical ED2
# panel c. Their defects are recorded in
# manuscript/extended_data_promotion_audit.csv and nothing was substituted for
# them.
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_figure_utils.R"))

for (fig in c("ED_01", "ED_02", "ED_03", "ED_06", "ED_08")) {
  manuscript_figure_main(fig)
}
