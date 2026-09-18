#!/usr/bin/env Rscript
# Authoritative control-only, animal-level, non-imputed cross-compartment
# marker detection and abundance workflow. The implementation lives in a shared
# R file so calculation and render-only contracts can be tested directly.
# Script: analysis/qc/render_compartment_abundance_figures.R
# Stage: qc_cross_dataset
# Scope: global
# Consumes: required data/processed/01_preprocessing/joint_compartment_qc/global/joint_compartment_qc_matrices.rds; config/marker_panels/wgcna_reference_marker_sets.csv; optional none declared in pipeline.yml
# Produces: results/qc/render_compartment_abundance_figures/global/tables/source_data; results/qc/render_compartment_abundance_figures/global/plots; results/qc/render_compartment_abundance_figures/global/reports
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Authoritative control-only, animal-level, non-imputed cross-compartment marker detection and abundance workflow.
#  

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "control_compartment_abundance_workflow.R"))
