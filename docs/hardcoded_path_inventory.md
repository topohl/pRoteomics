# Hard-Coded Path Inventory

> Audit/historical snapshot. This inventory is retained for maintenance context
> and may be stale unless explicitly refreshed. Use `pipeline.yml` as the active
> source of truth and `WORKFLOW.md` / `RUN_ORDER.md` for current workflow
> guidance.

This inventory covers active R scripts discovered before refactoring. Legacy/testing/deprecated scripts also contain older paths and should remain excluded from canonical pipeline execution unless explicitly revived.

## Active Scripts With Machine-Specific Paths

- `archive/01_preprocessing/01_impute.r`: `S:/...` metadata, raw matrix and imputed output roots.
- `archive/01_preprocessing/02_excel_convert.r`: `S:/...` imputed workbook, grouped folder, metadata and Morpheus output roots.
- `analysis/preprocessing/extract_protigy_contrasts.R`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `archive/01_preprocessing/04_format_metadata.r`: refactored in Phase 3; no committed active machine path remains.
- `archive/01_preprocessing/05_metadata_create.r`: refactored in Phase 3; no committed active machine path remains.
- `analysis/preprocessing/build_module_score_metadata.R`: `S:/...` processed matrix, behavior inputs and module-score output roots.
- `analysis/preprocessing/map_protein_identifiers.R`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `analysis/qc/assess_sample_quality.R`: `S:/...` QC workbook input and result folder.
- `analysis/qc/assess_marker_rank_abundance.R`: `S:/...` imputed matrix inputs and QC result folder.
- `analysis/qc/assess_pca_confounding.R`: `S:/...` GCT input and PCA output folder.
- `03_qc_exploration/06_pcaPlot_Neha.r`: `S:/...` Neha collaboration GCT input and PCA output folder.
- `analysis/qc/partition_variance.R`: `S:/...` variance-partition inputs and output folder.
- `03_qc_exploration/08_boxplotBonanza.r`: `S:/...` and `C:/...` local input/output paths.
- `analysis/differential_abundance/run_clusterprofiler_enrichment.R`: refactored to repo-relative config defaults; no committed active machine path remains.
- `analysis/differential_abundance/compare_go_enrichment.R`: refactored to manifest/config defaults; no committed active machine path remains.
- `analysis/enrichment/run_ewce_celltype_enrichment.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/wgcna/build_wgcna_modules.R`: central WGCNA model producer; canonical entrypoint in the active registry.
- `analysis/wgcna/build_module_spatial_networks.R`: legacy optional helper; no committed active machine path remains.
- `analysis/wgcna/build_curated_overlap_programs.R`: canonical curated overlap program builder; old overlap entrypoints were removed.
- `analysis/wgcna/score_module_activity.R`: canonical source-scoped module activity scorer; old module-score entrypoints were removed.
- `analysis/spatial_networks/build_spatial_networks.R`: refactored in Phase 4; no committed active machine path remains.
- `analysis/spatial_networks/build_differential_networks.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/spatial_networks/test_network_stability.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/spatial_networks/test_differential_network_stability.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/spatial_networks/render_differential_network_figures.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/spatial_networks/render_network_chord_diagram.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/integration/test_behaviour_proteomics_associations.R`: `S:/...` proteomics, behavior and output paths.
- `analysis/integration/test_network_behaviour_coupling.R`: refactored in Phase 3; no committed active machine path remains.

## Refactoring Rule

When each module is migrated, replace local roots with `source("R/paths.R")` and route inputs/outputs through:

- `path_raw()`
- `path_metadata()`
- `path_external()`
- `path_processed()`
- `path_results()`

Use `config/*.local.yml` for private machine paths and keep those files ignored by Git.

## Safe-To-Refactor Assessment

- Already partially/refactored: `analysis/differential_abundance/run_clusterprofiler_enrichment.R`, `analysis/differential_abundance/compare_go_enrichment.R`, `analysis/preprocessing/extract_protigy_contrasts.R`, `archive/01_preprocessing/04_format_metadata.r`, `archive/01_preprocessing/05_metadata_create.r`, `analysis/preprocessing/map_protein_identifiers.R`, `analysis/enrichment/run_ewce_celltype_enrichment.R`, `analysis/wgcna/build_curated_overlap_programs.R`, `analysis/wgcna/score_module_activity.R`, `analysis/wgcna/build_module_spatial_networks.R`, `analysis/spatial_networks/build_spatial_networks.R`, `analysis/spatial_networks/build_differential_networks.R`, `analysis/spatial_networks/test_network_stability.R`, `analysis/spatial_networks/test_differential_network_stability.R`, `analysis/spatial_networks/render_differential_network_figures.R`, `analysis/spatial_networks/render_network_chord_diagram.R`, and `analysis/integration/test_network_behaviour_coupling.R`.
- Requires data-aware review before edits: `archive/01_preprocessing/01_impute.r`, `archive/01_preprocessing/02_excel_convert.r`, `analysis/preprocessing/build_module_score_metadata.R`, `03_qc_exploration/*.r`, and `analysis/integration/test_behaviour_proteomics_associations.R`.
- Leave unchanged unless explicitly revived: `archive/05_celltype_enrichment_EWCE/90_EWCE_legacy.r`, `90_testing/`, `99_deprecated/`.
