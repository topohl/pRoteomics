# Hard-Coded Path Inventory

This inventory covers active R scripts discovered before refactoring. Legacy/testing/deprecated scripts also contain older paths and should remain excluded from canonical pipeline execution unless explicitly revived.

## Active Scripts With Machine-Specific Paths

- `01_preprocessing/01_impute.r`: `S:/...` metadata, raw matrix and imputed output roots.
- `01_preprocessing/02_excel_convert.r`: `S:/...` imputed workbook, grouped folder, metadata and Morpheus output roots.
- `01_preprocessing/03_gct_extractR.r`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `01_preprocessing/04_format_metadata.r`: refactored in Phase 3; no committed active machine path remains.
- `01_preprocessing/05_metadata_create.r`: refactored in Phase 3; no committed active machine path remains.
- `01_preprocessing/06_merged_metadata_module_score.r`: `S:/...` processed matrix, behavior inputs and module-score output roots.
- `02_id_mapping/01_MapThatProt_batch.r`: refactored in Phase 4; no committed active machine path or `setwd()` remains.
- `03_qc_exploration/01_sample_qc_quicksearch.r`: `S:/...` QC workbook input and result folder.
- `03_qc_exploration/04_marker_rank_abundance_qc.r`: `S:/...` imputed matrix inputs and QC result folder.
- `03_qc_exploration/05_pca_confounding_qc.r`: `S:/...` GCT input and PCA output folder.
- `03_qc_exploration/06_pcaPlot_Neha.r`: `S:/...` Neha collaboration GCT input and PCA output folder.
- `03_qc_exploration/06_variance_partitioning.r`: `S:/...` variance-partition inputs and output folder.
- `03_qc_exploration/08_boxplotBonanza.r`: `S:/...` and `C:/...` local input/output paths.
- `04_differential_expression_enrichment/01_clusterProfiler.r`: refactored to repo-relative config defaults; no committed active machine path remains.
- `04_differential_expression_enrichment/02_compareGO.r`: refactored to manifest/config defaults; no committed active machine path remains.
- `05_celltype_enrichment_EWCE/01_EWCE_E9.r`: refactored in Phase 3; no committed active machine path remains.
- `06_modules_WGCNA/01_WGCNA v.2.0.0.r`: `S:/...` inputs/outputs; central WGCNA model producer, documented for data-aware refactor.
- `06_modules_WGCNA/02_module_spatial_networks.r`: refactored in Phase 3; no committed active machine path remains.
- `06_modules_WGCNA/02_curated_overlap_programs.r`: canonical curated overlap program builder; old overlap entrypoints were removed.
- `06_modules_WGCNA/90_module_score_v0.0.1.r`: `S:/...` inputs/outputs; older superseded module-score version.
- `06_modules_WGCNA/03_score_module_activity.R`: canonical source-scoped module activity scorer; old module-score entrypoints were removed.
- `07_spatial_networks/01_network_spatial_relations.r`: refactored in Phase 4; no committed active machine path remains.
- `07_spatial_networks/02_differential_networks.r`: refactored in Phase 3; no committed active machine path remains.
- `07_spatial_networks/03_bootstrap_network_stability.r`: refactored in Phase 3; no committed active machine path remains.
- `07_spatial_networks/04_bootstrap_differential_network_stability.r`: refactored in Phase 3; no committed active machine path remains.
- `07_spatial_networks/05_bootstrap_differential_network_figures.r`: refactored in Phase 3; no committed active machine path remains.
- `07_spatial_networks/06_chord_diagram.r`: refactored in Phase 3; no committed active machine path remains.
- `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`: `S:/...` proteomics, behavior and output paths.
- `08_behavior_physio_coupling/02_network_behavior_coupling.r`: refactored in Phase 3; no committed active machine path remains.

## Refactoring Rule

When each module is migrated, replace local roots with `source("R/paths.R")` and route inputs/outputs through:

- `path_raw()`
- `path_metadata()`
- `path_external()`
- `path_processed()`
- `path_results()`

Use `config/*.local.yml` for private machine paths and keep those files ignored by Git.

## Safe-To-Refactor Assessment

- Already partially/refactored: `04_differential_expression_enrichment/01_clusterProfiler.r`, `04_differential_expression_enrichment/02_compareGO.r`, `01_preprocessing/03_gct_extractR.r`, `01_preprocessing/04_format_metadata.r`, `01_preprocessing/05_metadata_create.r`, `02_id_mapping/01_MapThatProt_batch.r`, `05_celltype_enrichment_EWCE/01_EWCE_E9.r`, `06_modules_WGCNA/02_module_spatial_networks.r`, `06_modules_WGCNA/02_curated_overlap_programs.r`, `06_modules_WGCNA/03_score_module_activity.R`, `07_spatial_networks/01_network_spatial_relations.r`, `07_spatial_networks/02_differential_networks.r`, `07_spatial_networks/03_bootstrap_network_stability.r`, `07_spatial_networks/04_bootstrap_differential_network_stability.r`, `07_spatial_networks/05_bootstrap_differential_network_figures.r`, `07_spatial_networks/06_chord_diagram.r`, and `08_behavior_physio_coupling/02_network_behavior_coupling.r`.
- Requires data-aware review before edits: `01_preprocessing/01_impute.r`, `01_preprocessing/02_excel_convert.r`, `01_preprocessing/06_merged_metadata_module_score.r`, `03_qc_exploration/*.r`, `06_modules_WGCNA/01_WGCNA v.2.0.0.r`, `06_modules_WGCNA/90_module_score_v0.0.1.r`, and `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`.
- Leave unchanged unless explicitly revived: `05_celltype_enrichment_EWCE/90_EWCE_legacy.r`, `90_testing/`, `99_deprecated/`.
