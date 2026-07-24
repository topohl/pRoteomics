#!/usr/bin/env Rscript
# Authoritative control-only, animal-level, non-imputed cross-compartment
# marker detection and abundance workflow. The implementation lives in a shared
# R file so calculation and render-only contracts can be tested directly.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "control_compartment_abundance_workflow_v2.R"))
