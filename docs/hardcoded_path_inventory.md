# Hard-Coded Path Inventory

> Audit/historical snapshot. This inventory is retained for maintenance context
> and may be stale unless explicitly refreshed. Use `pipeline.yml` as the active
> source of truth and `WORKFLOW.md` / `RUN_ORDER.md` for current workflow
> guidance.

This inventory covers active R scripts discovered before refactoring. Legacy/testing/deprecated scripts also contain older paths and should remain excluded from canonical pipeline execution unless explicitly revived.

## Active Scripts With Machine-Specific Paths

- `archive/01_preprocessing/01_impute.r`: `S:/...` metadata, raw matrix and imputed output roots.
- `archive/01_preprocessing/02_excel_convert.r`: `S:/...` imputed workbook, grouped folder, metadata and Morpheus output roots.
- `analysis/01_preprocessing/extract_protigy_contrasts.R`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `archive/01_preprocessing/04_format_metadata.r`: refactored in Phase 3; no committed active machine path remains.
- `archive/01_preprocessing/05_metadata_create.r`: refactored in Phase 3; no committed active machine path remains.
- `analysis/01_preprocessing/build_module_score_metadata.R`: `S:/...` processed matrix, behavior inputs and module-score output roots.
- `analysis/01_preprocessing/map_protein_identifiers.R`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `analysis/02_qc/assess_sample_quality.R`: `S:/...` QC workbook input and result folder.
- `analysis/02_qc/assess_marker_rank_abundance.R`: `S:/...` imputed matrix inputs and QC result folder.
- `analysis/02_qc/assess_pca_confounding.R`: `S:/...` GCT input and PCA output folder.
- `03_qc_exploration/06_pcaPlot_Neha.r`: `S:/...` Neha collaboration GCT input and PCA output folder.
- `analysis/02_qc/partition_variance.R`: `S:/...` variance-partition inputs and output folder.
- `03_qc_exploration/08_boxplotBonanza.r`: `S:/...` and `C:/...` local input/output paths.
- `analysis/04_differential_abundance/run_clusterprofiler_enrichment.R`: refactored to repo-relative config defaults; no committed active machine path remains.
- `analysis/04_differential_abundance/compare_go_enrichment.R`: refactored to manifest/config defaults; no committed active machine path remains.
- `analysis/06_gsea/run_ewce_celltype_enrichment.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/05_wgcna/build_wgcna_modules.R`: central WGCNA model producer; canonical entrypoint in the active registry.
- `analysis/05_wgcna/build_module_spatial_networks.R`: legacy optional helper; no committed active machine path remains.
- `analysis/05_wgcna/build_curated_overlap_programs.R`: canonical curated overlap program builder; old overlap entrypoints were removed.
- `analysis/05_wgcna/score_module_activity.R`: canonical source-scoped module activity scorer; old module-score entrypoints were removed.
- `analysis/07_spatial_networks/build_spatial_networks.R`: refactored in Phase 4; no committed active machine path remains.
- `analysis/07_spatial_networks/build_differential_networks.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/07_spatial_networks/test_network_stability.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/07_spatial_networks/test_differential_network_stability.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/07_spatial_networks/render_differential_network_figures.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/07_spatial_networks/render_network_chord_diagram.R`: refactored in Phase 3; no committed active machine path remains.
- `analysis/08_integration/test_behaviour_proteomics_associations.R`: `S:/...` proteomics, behavior and output paths.
- `analysis/08_integration/test_network_behaviour_coupling.R`: refactored in Phase 3; no committed active machine path remains.

## Refactoring Rule

When each module is migrated, replace local roots with `source("R/paths.R")` and route inputs/outputs through:

- `path_raw()`
- `path_metadata()`
- `path_external()`
- `path_processed()`
- `path_results()`

Use `config/*.local.yml` for private machine paths and keep those files ignored by Git.

## Safe-To-Refactor Assessment

- Already partially/refactored: `analysis/04_differential_abundance/run_clusterprofiler_enrichment.R`, `analysis/04_differential_abundance/compare_go_enrichment.R`, `analysis/01_preprocessing/extract_protigy_contrasts.R`, `archive/01_preprocessing/04_format_metadata.r`, `archive/01_preprocessing/05_metadata_create.r`, `analysis/01_preprocessing/map_protein_identifiers.R`, `analysis/06_gsea/run_ewce_celltype_enrichment.R`, `analysis/05_wgcna/build_curated_overlap_programs.R`, `analysis/05_wgcna/score_module_activity.R`, `analysis/05_wgcna/build_module_spatial_networks.R`, `analysis/07_spatial_networks/build_spatial_networks.R`, `analysis/07_spatial_networks/build_differential_networks.R`, `analysis/07_spatial_networks/test_network_stability.R`, `analysis/07_spatial_networks/test_differential_network_stability.R`, `analysis/07_spatial_networks/render_differential_network_figures.R`, `analysis/07_spatial_networks/render_network_chord_diagram.R`, and `analysis/08_integration/test_network_behaviour_coupling.R`.
- Requires data-aware review before edits: `archive/01_preprocessing/01_impute.r`, `archive/01_preprocessing/02_excel_convert.r`, `analysis/01_preprocessing/build_module_score_metadata.R`, `03_qc_exploration/*.r`, and `analysis/08_integration/test_behaviour_proteomics_associations.R`.
- Leave unchanged unless explicitly revived: `archive/05_celltype_enrichment_EWCE/90_EWCE_legacy.r`, `90_testing/`, `99_deprecated/`.
