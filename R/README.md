# Shared R Helpers

This folder contains shared functions only. Analysis entrypoints live in the
numbered pipeline folders and are registered in `pipeline.yml`.

## Core Helpers

- `paths.R`: repository-root detection and canonical path constructors.
- `dataset_config.R`: dataset names, capabilities, and `--dataset` resolution.
- `script_runtime.R`: shared CLI, dry-run, input-status, and manifest plumbing for entrypoint scripts.
- `pipeline_registry.R`: reads and validates `pipeline.yml`.
- `validation_utils.R`, `schema_validation.R`: run manifests and table/schema checks.

## Data And Domain Helpers

- `dataset_inputs.R`: dataset-aware matrix and metadata discovery.
- `protein_mapping_utils.R`: identifier mapping helpers.
- `enrichment_io.R`, `enrichment_plots.R`: enrichment manifests, summaries, and plots.
- `module_contracts.R`, `wgcna_downstream_utils.R`, `module_stats.R`: WGCNA/module contracts and downstream summaries.
- `qc_exploration_utils.R`: QC input/output and plotting helpers.
- `spatial_network_utils.R`: spatial-network utilities.
- `integration_utils.R`: manuscript-level integration table helpers.
- `export_helpers.R`, `pride_helpers.R`: PRIDE/journal export helpers.
- `plotting_nature.R`: shared manuscript-style plotting defaults.

Prefer extending these helpers when logic is shared by at least two scripts or
when a table/path contract must stay consistent across stages. Keep one-off
analysis logic inside the numbered script that owns it.
