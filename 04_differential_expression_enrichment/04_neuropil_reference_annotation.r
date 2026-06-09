#!/usr/bin/env Rscript

# ================================================================
# Script: 04_differential_expression_enrichment/04_neuropil_reference_annotation.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: required data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; data/processed/04_differential_expression_enrichment/clusterProfiler/neuron_neuropil/clusterProfiler_manifest.csv; optional config/marker_panels/wgcna_reference_marker_sets.csv.
# Produces: results/tables/04_differential_expression_enrichment/neuropil_reference_annotation/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Runs after neuron_neuropil enrichment exists.
# ================================================================

# Publication-facing alias for neuropil reference annotation of microglia ROI
# enrichment results. The older file name is retained for backward-compatible
# commands.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("04_differential_expression_enrichment", "04_neuropil_contamination_annotation.r"))
