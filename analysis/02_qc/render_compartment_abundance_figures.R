#!/usr/bin/env Rscript
# Authoritative control-only, animal-level, non-imputed cross-compartment
# marker detection and abundance workflow. The implementation lives in a shared
# R file so calculation and render-only contracts can be tested directly.
# Script: analysis/02_qc/render_compartment_abundance_figures.R
# Stage: qc_cross_dataset
# Scope: global
# Consumes: required data/processed/01_preprocessing/joint_compartment_qc/global/joint_compartment_qc_matrices.rds; config/marker_panels/wgcna_reference_marker_sets.csv; optional none declared in pipeline.yml
# Produces: results/source_data/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/; results/figures/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/; results/reports/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Authoritative control-only, animal-level, non-imputed cross-compartment marker detection and abundance workflow.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "control_compartment_abundance_workflow.R"))
