# Hard-Coded Path Inventory

This inventory covers active R scripts discovered before refactoring. Legacy/testing/deprecated scripts also contain older paths and should remain excluded from canonical pipeline execution unless explicitly revived.

## Active Scripts With Machine-Specific Paths

- `01_preprocessing/03_gct_extractR.r`: `S:/...` input/output roots and `setwd()`.
- `01_preprocessing/04_format_metadata.r`: `S:/...` metadata input and processed workbook output.
- `01_preprocessing/05_metadata_create.r`: `S:/...` metadata workbook directory.
- `02_id_mapping/01_MapThatProt_batch.r`: `S:/...` working root, mapped/unmapped/report folder derivation, and `setwd()`.
- `03_qc_exploration/01_qc_protein_peptide_plot.r`: `S:/...` QC workbook input and result folder.
- `03_qc_exploration/02_rank_abundance_plot_E9.r`: `S:/...` imputed matrix inputs and QC result folder.
- `03_qc_exploration/03_pcaPlot.r`: `S:/...` GCT input and PCA output folder.
- `03_qc_exploration/06_pcaPlot_Neha.r`: `S:/...` Neha collaboration GCT input and PCA output folder.
- `03_qc_exploration/07_varPart.r`: `S:/...` variance-partition inputs and output folder.
- `03_qc_exploration/08_boxplotBonanza.r`: `S:/...` and `C:/...` local input/output paths.
- `04_differential_expression_enrichment/01_clusterProfiler.r`: refactored to repo-relative config defaults; no committed active machine path remains.
- `04_differential_expression_enrichment/02_compareGO.r`: refactored to manifest/config defaults; no committed active machine path remains.
- `07_spatial_networks/01_network_spatial_relations.r`: `/Users/...` active protein, metadata and output paths plus commented `S:/...` alternatives.
- `07_spatial_networks/02_differential_networks.r`: `/Users/...` active RDS and output paths plus commented `S:/...` alternatives.
- `07_spatial_networks/03_bootstrap_network_stability.r`: `/Users/...` active RDS and output paths.
- `07_spatial_networks/04_bootstrap_differential_network_stability.r`: `/Users/...` active RDS and output paths plus commented `S:/...` alternatives.
- `07_spatial_networks/05_bootstrap_differential_network_figures.r`: `/Users/...` active bootstrap path plus commented `S:/...` alternative.
- `07_spatial_networks/06_chord_diagram.r`: `S:/...` mapped input and result output folders.
- `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`: `S:/...` proteomics, behavior and output paths.
- `08_behavior_physio_coupling/02_network_behavior_coupling.r`: `/Users/...` active behavior/network paths plus commented `S:/...` alternatives.

## Refactoring Rule

When each module is migrated, replace local roots with `source("R/paths.R")` and route inputs/outputs through:

- `path_raw()`
- `path_metadata()`
- `path_external()`
- `path_processed()`
- `path_results()`

Use `config/*.local.yml` for private machine paths and keep those files ignored by Git.

## Safe-To-Refactor Assessment

- Already partially refactored: `04_differential_expression_enrichment/01_clusterProfiler.r`, `04_differential_expression_enrichment/02_compareGO.r`.
- Next safe path-only pass: `05_celltype_enrichment_EWCE/01_EWCE_E9.r`, `06_modules_WGCNA/02_module_spatial_networks.r`, `06_modules_WGCNA/03_overlap_modules.r`, `07_spatial_networks/03_bootstrap_network_stability.r`, `07_spatial_networks/04_bootstrap_differential_network_stability.r`, `07_spatial_networks/05_bootstrap_differential_network_figures.r`, `08_behavior_physio_coupling/02_network_behavior_coupling.r`.
- Requires data-aware review before edits: `01_preprocessing/03_gct_extractR.r`, `02_id_mapping/01_MapThatProt_batch.r`, `06_modules_WGCNA/01_WGCNA v.2.0.0.r`, `07_spatial_networks/01_network_spatial_relations.r`, `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`.
- Leave unchanged unless explicitly revived: `05_celltype_enrichment_EWCE/90_EWCE_legacy.r`, `90_testing/`, `99_deprecated/`.
