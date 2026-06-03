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
| module-score merged metadata | module-score metadata merge contract | `results/module_scores/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx` |
| WGCNA module definitions | `wgcna_module_contract.yml` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_definitions_for_downstream.csv` |
| curated overlap programs | curated overlap program contract | `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs.xlsx`; `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs_long.csv` |
| module activity scores | module score contract | `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_scores_per_sample.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_feature_mapping_trace.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_gene_coverage.csv`; `results/logs/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/run_manifest.yml` |
| WGCNA/DE/GSEA overlap | overlap bridge contract | `results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/WGCNA_vs_DE_GSEA_overlap.csv` |
| WGCNA marker traits | annotation-trait contract | `results/tables/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv`; `results/source_data/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv` |
| WGCNA module/supermodule group effects | `validate_wgcna_group_effects()` | `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_marker_trait_correlations.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_marker_trait_correlations.csv` |
| WGCNA biological annotation | `validate_wgcna_module_annotation()` | `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`; `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv` |
| WGCNA interpretable summary | `validate_wgcna_interpretable_summary()` | `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_module_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/all/WGCNA_cross_dataset_supermodule_program_summary.csv` |
| spatial network object | dataset/spatial-unit spatial network contract | `data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/network_spatial_relations_objects.rds`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_mean_expression_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_spearman_correlation_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/standardized_sample_metadata_used.csv`; `results/logs/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/run_manifest.yml` |
| biological claims table | `biological_claims_table.yml` | `results/tables/biological_claims_table.csv` |

Legacy compatibility note: `results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx` may be produced for older scripts, but dataset-aware scripts should consume the dataset-scoped path under `results/module_scores/<dataset>/`.

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
