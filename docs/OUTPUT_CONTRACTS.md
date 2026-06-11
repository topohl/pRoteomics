# Output Contracts

Canonical generated outputs use dataset-aware paths:

```text
data/processed/<module>/<substep>/<dataset>/
results/tables/<module>/<substep>/<dataset>/
results/figures/<module>/<substep>/<dataset>/
results/source_data/<module>/<substep>/<dataset>/
results/reports/<module>/<substep>/<dataset>/
results/logs/<module>/<substep>/<dataset>/
```

Major downstream tables should be validated with `validate_table_schema(df, schema_name, strict = TRUE)` before writing.

## Required Final Tables

| output | schema or contract | expected location |
|---|---|---|
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv` |
| compareGO input manifest | `compareGO_manifest.yml` | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv` |
| biological program summary | documented in `docs/file_contracts.tsv` | `results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv` |
| targeted microglia ROI signature claims | documented in `docs/file_contracts.tsv` | `results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia/microglia_signature_claims_ready.csv` |
| module-score merged metadata | module-score metadata merge contract | `data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx`; QC/report/log sidecars under `results/{tables,reports,logs}/01_preprocessing/06_merged_metadata_module_score/<dataset>/` |
| WGCNA module definitions | `wgcna_module_contract.yml` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_definitions_for_downstream.csv` |
| curated overlap programs | curated overlap program contract | `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs.xlsx`; `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs_long.csv` |
| module activity scores | module score contract | `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_scores_per_sample.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_feature_mapping_trace.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_gene_coverage.csv`; `results/logs/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/run_manifest.yml` |
| WGCNA/DE/GSEA overlap | overlap bridge contract | `results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/WGCNA_vs_DE_GSEA_overlap.csv` |
| Reference marker evidence registry | marker registry contract | `results/tables/03_qc_exploration/reference_marker_import/reference_marker_evidence_long.csv`; `results/tables/03_qc_exploration/reference_marker_import/reference_marker_sets_ranked.csv`; `config/marker_panels/wgcna_reference_marker_sets.csv` |
| Empirical ROI marker discovery | empirical marker contract | `results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv`; `results/source_data/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv` |
| WGCNA marker traits | annotation-trait contract | `results/tables/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv`; `results/source_data/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv` |
| WGCNA module/supermodule group effects | `validate_wgcna_group_effects()` | `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_marker_trait_correlations.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_marker_trait_correlations.csv` |
| WGCNA biological annotation | `validate_wgcna_module_annotation()` | `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`; `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv` |
| WGCNA interpretable summary | `validate_wgcna_interpretable_summary()` | `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_module_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/all/WGCNA_cross_dataset_supermodule_program_summary.csv` |
| WGCNA module complex/organelle architecture | standard integration evidence contract | `results/tables/06_modules_WGCNA/module_complex_architecture/<dataset>/module_complex_architecture.csv`; mirrored under `results/source_data/06_modules_WGCNA/module_complex_architecture/<dataset>/` |
| WGCNA module robustness/sensitivity | standard integration evidence contract | `results/tables/06_modules_WGCNA/module_robustness_sensitivity/<dataset>/module_robustness_sensitivity.csv`; mirrored under `results/source_data/06_modules_WGCNA/module_robustness_sensitivity/<dataset>/` |
| module-behavior coupling | standard integration evidence contract | `results/tables/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/module_behavior_coupling.csv`; mirrored under `results/source_data/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/` |
| cross-compartment program atlas | `validate_cross_compartment_program_atlas()` | `results/tables/10_biological_integration/cross_compartment_program_atlas/global/cross_compartment_program_atlas_long.csv`; `results/tables/10_biological_integration/cross_compartment_program_atlas/global/cross_compartment_program_atlas.csv` |
| manuscript program summary | `validate_manuscript_program_summary()` | `results/tables/10_biological_integration/manuscript_program_summary/global/manuscript_program_summary.csv` |
| evidence priority matrix | `validate_evidence_priority_matrix()` | `results/tables/10_biological_integration/evidence_priority_matrix/global/evidence_priority_matrix.csv` |
| spatial network object | dataset/spatial-unit spatial network contract | `data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/network_spatial_relations_objects.rds`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_mean_expression_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_spearman_correlation_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/standardized_sample_metadata_used.csv`; `results/logs/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/run_manifest.yml` |
| biological claims table | `biological_claims_table.yml` | `results/tables/biological_claims_table.csv` |

Legacy compatibility note: existing `results/module_scores/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx` and `results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx` files are read-only fallbacks. Fallback use emits warnings, and the global fallback requires `PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA=true`.

## Manuscript And Source Data

The active export layer in `09_export_pride_journal/` collects final products under:

```text
results/manuscript/figure_*/
results/manuscript/extended_data/
results/manuscript/supplementary_tables/
results/manuscript/source_data/
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

Export scripts collect final products only; they do not recompute scientific analyses. Active scripts should also write `run_manifest.yml`, `sessionInfo.txt`, config snapshots, and input file hashes where applicable.
