#!/usr/bin/env Rscript

# Explicit behavioural Extended Data entry point. As with Figure 1, no statistic
# is computed here and none is computed anywhere in this repository: the
# behavioural analysis is frozen upstream and imported with hashes. This
# validates the declared panel assets and materialises them.
#
# Two figures, canonicalised in Phase 6A. Extended Data 5 is how the early window
# was measured and the complete a-priori model registry; Extended Data 9 is the
# two secondary features and everything the study can say about sex. They were
# one figure until Phase 6A, which put two different arguments on one page.
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_figure_utils.R"))

for (fig in c("ED_behaviour_coverage", "ED_behaviour_secondary")) {
  manuscript_figure_main(fig)
}
