# Run Order Command Reference

This file is a command reference for the active registry in `pipeline.yml`.
Conceptual workflow explanations live in [WORKFLOW.md](WORKFLOW.md), and WGCNA
interpretation layers are summarized in [docs/WGCNA_WORKFLOW.md](docs/WGCNA_WORKFLOW.md).

## Launcher

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset <dataset> --stage <stage> --dry-run
Rscript run_dataset_pipeline.R --dataset <dataset> --stage <stage>
```

Valid datasets:

```text
neuron_neuropil
neuron_soma
microglia
```

Use `--dataset all` for global registry stages such as `integration` and
`export`. Each launcher run writes a manifest under `results/logs/pipeline/`.

## Registry Stages

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage core
Rscript run_dataset_pipeline.R --dataset all --stage qc_global
Rscript run_dataset_pipeline.R --dataset <dataset> --stage qc
Rscript run_dataset_pipeline.R --dataset <dataset> --stage enrichment
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_wgcna
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
Rscript run_dataset_pipeline.R --dataset <dataset> --stage networks
Rscript run_dataset_pipeline.R --dataset <dataset> --stage coupling
Rscript run_dataset_pipeline.R --dataset all --stage integration
Rscript run_dataset_pipeline.R --dataset all --stage export
```

Dry-run the full manuscript path first:

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
done
Rscript run_dataset_pipeline.R --dataset all --stage integration --dry-run
Rscript run_dataset_pipeline.R --dataset all --stage export --dry-run
```

## Core

Direct script commands:

```bash
Rscript 01_preprocessing/03_gct_extractR.r --dataset <dataset> --dry-run
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset <dataset> --dry-run
Rscript 02_id_mapping/01_MapThatProt_batch.r --dataset <dataset> --dry-run
```

Primary handoff outputs:

```text
data/processed/01_preprocessing/gct_extractR/<dataset>/
data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/
```

Use `PROTEOMICS_RECOMPUTE=true`, `--recompute`, or `--force-rerun` only when
intentionally refreshing core handoff tables.

## QC

Direct script commands:

```bash
Rscript 03_qc_exploration/04b_import_reference_marker_sources.r --dry-run
Rscript 03_qc_exploration/05_empirical_roi_marker_discovery.r --dry-run
Rscript 03_qc_exploration/00_dataset_qc_report.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/01_sample_qc_quicksearch.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/02_missingness_diagnostics.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/03_replicate_consistency.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/04_marker_rank_abundance_qc.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/04d_compartment_marker_fidelity.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/05_pca_confounding_qc.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/06_variance_partitioning.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/06_wgcna_marker_trait_export.r --dataset <dataset> --dry-run
Rscript 03_qc_exploration/07_qc_biology_confounding_report.r --dataset <dataset> --dry-run
```

Key outputs:

```text
config/marker_panels/wgcna_reference_marker_sets.csv
results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/
results/reports/03_qc_exploration/00_dataset_qc_report/<dataset>/
results/reports/03_qc_exploration/07_qc_biology_confounding_report/<dataset>/
```

## Enrichment

Direct script commands:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dataset <dataset> --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dataset <dataset> --dry-run
Rscript 04_differential_expression_enrichment/04_neuropil_reference_annotation.r --dataset <dataset> --dry-run
Rscript 04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r --dataset microglia --dry-run
Rscript 04_differential_expression_enrichment/03_biological_program_summary.r --dataset <dataset> --dry-run
Rscript 04_differential_expression_enrichment/06_compareGO_spatial_program_atlas.r --dataset <dataset> --dry-run
Rscript 04_differential_expression_enrichment/07_external_stress_disease_signature_overlap.r --dry-run
Rscript 05_celltype_enrichment_EWCE/01_EWCE_E9.r --dataset <dataset> --dry-run
```

Key outputs:

```text
data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv
data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv
results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv
results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia/
```

## WGCNA

Registry commands:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_wgcna --dry-run
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream --dry-run
```

Direct script commands:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/01a_compare_GO_recurrent_proteins.r --dry-run
Rscript 06_modules_WGCNA/02_curated_overlap_programs.r --dry-run
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/08_module_complex_architecture.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/09_module_robustness_sensitivity.r --dataset <dataset> --dry-run
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all --dry-run
```

Layer distinction:

```text
06_modules_WGCNA/01_WGCNA.r                         network/module construction
06_modules_WGCNA/03_score_module_activity.R         module scoring
06_modules_WGCNA/05_module_supermodule_group_effects.r group-effect modelling
06_modules_WGCNA/06_annotate_module_microenvironment.r biological annotation
06_modules_WGCNA/07_wgcna_interpretable_summary.r   manuscript-facing summary
```

Key outputs:

```text
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/
results/tables/06_modules_WGCNA/module_score/<dataset>/
results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv
results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv
results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv
results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_interpretable_summary.xlsx
```

`01_WGCNA.r` recomputes core WGCNA state. The downstream WGCNA scripts consume
existing modules and are safe to rerun for reporting/annotation updates.

## Spatial Networks

Direct script commands:

```bash
Rscript 07_spatial_networks/01_network_spatial_relations.r --dataset <dataset> --dry-run
Rscript 07_spatial_networks/02_differential_networks.r --dataset <dataset> --dry-run
Rscript 07_spatial_networks/03_bootstrap_network_stability.r --dataset <dataset> --dry-run
Rscript 07_spatial_networks/04_bootstrap_differential_network_stability.r --dataset <dataset> --dry-run
Rscript 07_spatial_networks/05_bootstrap_differential_network_figures.r --dataset <dataset> --dry-run
Rscript 07_spatial_networks/06_chord_diagram.r --dataset <dataset> --dry-run
```

Key output:

```text
data/processed/07_spatial_networks/network_spatial_relations/<dataset>/*/network_spatial_relations_objects.rds
```

## Coupling

Direct script commands:

```bash
Rscript 08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r --dataset neuron_soma --dry-run
Rscript 08_behavior_physio_coupling/02_network_behavior_coupling.r --dataset neuron_neuropil --dry-run
Rscript 08_behavior_physio_coupling/03_module_behavior_coupling.r --dataset <dataset> --dry-run
```

Key outputs:

```text
results/tables/08_behavior_physio_coupling/network_behavior_coupling/
results/tables/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/
```

## Integration

Direct script commands:

```bash
Rscript 10_biological_integration/01_cross_compartment_program_atlas.r --dry-run
Rscript 10_biological_integration/02_manuscript_program_summary.r --dry-run
Rscript 10_biological_integration/03_evidence_priority_matrix.r --dry-run
```

Key outputs:

```text
results/tables/10_biological_integration/cross_compartment_program_atlas/global/
results/tables/10_biological_integration/manuscript_program_summary/global/
results/tables/10_biological_integration/evidence_priority_matrix/global/
```

## Export

Direct script commands:

```bash
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dry-run
Rscript 09_export_pride_journal/03_export_processed_pg_matrix_package.R --dry-run
Rscript 09_export_pride_journal/04_make_supplementary_tables.R --dry-run
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dry-run
Rscript 09_export_pride_journal/05_validate_pride_submission.R --dry-run
Rscript 09_export_pride_journal/06_make_methods_summary.R --dry-run
Rscript 09_export_pride_journal/07_make_biological_claims_table.R --dry-run
```

Key outputs:

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/manifests/
pride_submission/validation/
pride_submission/methods/
results/tables/biological_claims_table.csv
```
