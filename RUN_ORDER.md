# Run Order Command Reference

This file is a command reference for the active registry in `pipeline.yml`.
Conceptual workflow explanations live in [WORKFLOW.md](WORKFLOW.md), and WGCNA
interpretation layers are summarized in [docs/WGCNA_WORKFLOW.md](docs/WGCNA_WORKFLOW.md).

`pipeline.yml` is the sole execution-order authority. Numeric filename prefixes
are stable historical identifiers, not a sortable execution plan; suffixes and
parallel analysis branches therefore do not imply registry order. The generated
index below is the complete ordered inventory used by contract validation.

## Launcher

```powershell
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
Rscript run_dataset_pipeline.R --dataset all --stage qc_cross_dataset
Rscript run_dataset_pipeline.R --dataset <dataset> --stage enrichment
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_wgcna
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
Rscript run_dataset_pipeline.R --dataset <dataset> --stage networks
Rscript run_dataset_pipeline.R --dataset <dataset> --stage coupling
Rscript run_dataset_pipeline.R --dataset all --stage integration
Rscript run_dataset_pipeline.R --dataset all --stage export
```

<!-- BEGIN GENERATED PIPELINE REGISTRY INDEX -->

## Authoritative script index (generated)

Do not hand-edit this block; run `Rscript tools/generate_pipeline_docs.R`.

1. `analysis/preprocessing/extract_protigy_contrasts.R` - stage `core`; scope `dataset_specific`
2. `analysis/preprocessing/build_module_score_metadata.R` - stage `core`; scope `dataset_specific`
3. `analysis/preprocessing/map_protein_identifiers.R` - stage `core`; scope `dataset_specific`
4. `analysis/preprocessing/build_joint_protigy_input.R` - stage `joint_qc_preprocessing`; scope `global`
5. `analysis/qc/assess_joint_compartment_quality.R` - stage `qc_global`; scope `global`
6. `analysis/qc/render_joint_compartment_qc_figures.R` - stage `qc_global`; scope `global`
7. `analysis/qc/build_reference_marker_registry.R` - stage `qc_global`; scope `global`
8. `analysis/qc/discover_empirical_roi_markers.R` - stage `qc_global`; scope `global`
9. `analysis/qc/assess_dataset_quality.R` - stage `qc`; scope `dataset_specific`
10. `analysis/qc/assess_sample_quality.R` - stage `qc`; scope `dataset_specific`
11. `analysis/qc/summarize_missingness.R` - stage `qc`; scope `dataset_specific`
12. `analysis/qc/assess_replicate_consistency.R` - stage `qc`; scope `dataset_specific`
13. `analysis/qc/assess_marker_rank_abundance.R` - stage `qc`; scope `dataset_specific`
14. `analysis/qc/summarize_marker_detectability.R` - stage `qc`; scope `dataset_specific`
15. `analysis/qc/assess_pca_confounding.R` - stage `qc`; scope `dataset_specific`
16. `analysis/qc/partition_variance.R` - stage `qc`; scope `dataset_specific`
17. `analysis/qc/export_marker_traits.R` - stage `qc`; scope `dataset_specific`
18. `analysis/qc/summarize_qc_confounding.R` - stage `qc`; scope `dataset_specific`
19. `analysis/qc/render_compartment_abundance_figures.R` - stage `qc_cross_dataset`; scope `global`
20. `analysis/differential_abundance/run_clusterprofiler_enrichment.R` - stage `enrichment`; scope `dataset_specific`
21. `analysis/differential_abundance/audit_gsea_protein_direction.R` - stage `enrichment`; scope `dataset_specific`
22. `analysis/differential_abundance/compare_go_enrichment.R` - stage `enrichment`; scope `dataset_specific`
23. `analysis/differential_abundance/annotate_neuropil_reference.R` - stage `enrichment`; scope `dataset_specific`
24. `analysis/differential_abundance/test_microglia_targeted_signatures.R` - stage `enrichment`; scope `dataset_specific`
25. `analysis/differential_abundance/summarize_biological_programs.R` - stage `enrichment`; scope `dataset_specific`
26. `analysis/differential_abundance/build_go_program_atlas.R` - stage `enrichment`; scope `dataset_specific`
27. `analysis/differential_abundance/compare_external_stress_signatures.R` - stage `enrichment`; scope `global`
28. `analysis/differential_abundance/validate_control_spatial_identity.R` - stage `enrichment`; scope `global`
29. `analysis/differential_abundance/build_sus_res_dap_atlas.R` - stage `enrichment`; scope `global`
30. `analysis/differential_abundance/audit_stress_response_biology.R` - stage `enrichment`; scope `global`
31. `analysis/enrichment/run_ewce_celltype_enrichment.R` - stage `enrichment`; scope `dataset_specific`
32. `analysis/wgcna/build_wgcna_modules.R` - stage `modules_wgcna`; scope `dataset_specific`
33. `analysis/wgcna/render_module_go_heatmaps.R` - stage `modules_downstream`; scope `dataset_specific`
34. `analysis/wgcna/compare_recurrent_module_proteins.R` - stage `modules_downstream`; scope `dataset_specific`
35. `analysis/wgcna/build_curated_overlap_programs.R` - stage `modules_downstream`; scope `global`
36. `analysis/wgcna/score_module_activity.R` - stage `modules_downstream`; scope `dataset_specific`
37. `analysis/wgcna/compare_module_enrichment_overlap.R` - stage `modules_downstream`; scope `dataset_specific`
38. `analysis/wgcna/build_module_identity_contract.R` - stage `modules_downstream`; scope `dataset_specific`
39. `analysis/wgcna/test_module_phenotypes.R` - stage `modules_downstream`; scope `dataset_specific`
40. `analysis/wgcna/annotate_module_microenvironment.R` - stage `modules_downstream`; scope `dataset_specific`
41. `analysis/wgcna/summarize_module_interpretation.R` - stage `modules_downstream`; scope `dataset_specific`
42. `analysis/wgcna/render_module_figures.R` - stage `modules_downstream`; scope `dataset_specific`
43. `analysis/wgcna/summarize_module_scores.R` - stage `modules_downstream`; scope `dataset_specific`
44. `analysis/wgcna/render_microglia_module_figures.R` - stage `modules_downstream`; scope `dataset_specific`
45. `analysis/wgcna/test_microglia_neuropil_independence.R` - stage `modules_downstream`; scope `dataset_specific`
46. `analysis/wgcna/render_microglia_independence_figures.R` - stage `modules_downstream`; scope `dataset_specific`
47. `analysis/wgcna/summarize_microglia_roi_specificity.R` - stage `modules_downstream`; scope `dataset_specific`
48. `analysis/wgcna/summarize_module_complex_architecture.R` - stage `modules_downstream`; scope `dataset_specific`
49. `analysis/wgcna/audit_module_robustness.R` - stage `modules_downstream`; scope `dataset_specific`
50. `analysis/wgcna/audit_microglia_module_claims.R` - stage `modules_downstream`; scope `dataset_specific`
51. `analysis/wgcna/summarize_microglia_module_claims.R` - stage `modules_downstream`; scope `dataset_specific`
52. `analysis/wgcna/audit_module_claim_readiness.R` - stage `modules_downstream`; scope `dataset_specific`
53. `analysis/wgcna/audit_module_label_coherence.R` - stage `networks`; scope `per_dataset`
54. `analysis/wgcna/adjudicate_module_labels.R` - stage `networks`; scope `per_dataset`
55. `analysis/wgcna/build_module_label_registry.R` - stage `networks`; scope `per_dataset`
56. `analysis/spatial_networks/build_spatial_networks.R` - stage `networks`; scope `dataset_specific`
57. `analysis/spatial_networks/build_differential_networks.R` - stage `networks`; scope `dataset_specific`
58. `analysis/spatial_networks/test_network_stability.R` - stage `networks`; scope `dataset_specific`
59. `analysis/spatial_networks/test_differential_network_stability.R` - stage `networks`; scope `dataset_specific`
60. `analysis/spatial_networks/render_differential_network_figures.R` - stage `networks`; scope `dataset_specific`
61. `analysis/spatial_networks/render_network_chord_diagram.R` - stage `networks`; scope `dataset_specific`
62. `analysis/spatial_validation/build_spatial_data_contract.R` - stage `networks`; scope `per_dataset`
63. `analysis/spatial_validation/quantify_bilateral_spatial_identity.R` - stage `networks`; scope `per_dataset`
64. `analysis/spatial_validation/quantify_empirical_compartments.R` - stage `networks`; scope `global`
65. `analysis/spatial_validation/quantify_module_bilateral_identity.R` - stage `networks`; scope `per_dataset`
66. `analysis/spatial_validation/annotate_module_celltypes.R` - stage `networks`; scope `per_dataset`
67. `analysis/spatial_validation/decompose_bilateral_variance.R` - stage `networks`; scope `per_dataset`
68. `analysis/spatial_validation/validate_spatial_foundations.R` - stage `networks`; scope `global`
69. `analysis/spatial_validation/build_module_spatial_atlas.R` - stage `networks`; scope `per_dataset`
70. `analysis/spatial_validation/build_protein_spatial_atlas.R` - stage `networks`; scope `global`
71. `analysis/spatial_validation/quantify_neuropil_detection_context.R` - stage `networks`; scope `global`
72. `analysis/spatial_validation/summarize_spatial_atlas.R` - stage `networks`; scope `global`
73. `analysis/spatial_validation/quantify_neuropil_precision.R` - stage `networks`; scope `global`
74. `analysis/spatial_validation/build_animal_spatial_networks.R` - stage `networks`; scope `global`
75. `analysis/spatial_validation/test_network_group_organization.R` - stage `networks`; scope `global`
76. `analysis/spatial_validation/validate_network_workbook.R` - stage `networks`; scope `global`
77. `analysis/spatial_validation/audit_ca2_slm_robustness.R` - stage `networks`; scope `global`
78. `analysis/spatial_validation/audit_stress_identity_robustness.R` - stage `networks`; scope `global`
79. `analysis/spatial_validation/summarize_ca2_slm_robustness.R` - stage `networks`; scope `global`
80. `analysis/integration/test_behaviour_proteomics_associations.R` - stage `coupling`; scope `dataset_specific`
81. `analysis/integration/test_network_behaviour_coupling.R` - stage `coupling`; scope `dataset_specific`
82. `analysis/integration/test_module_behaviour_coupling.R` - stage `coupling`; scope `dataset_specific`
83. `analysis/integration/audit_animal_id_integrity.R` - stage `coupling`; scope `global`
84. `analysis/integration/build_cross_compartment_atlas.R` - stage `integration`; scope `global`
85. `analysis/integration/summarize_programs_for_manuscript.R` - stage `integration`; scope `global`
86. `analysis/publication_source_data/build_biological_claims_table.R` - stage `integration`; scope `global`
87. `analysis/integration/build_evidence_priority_matrix.R` - stage `integration`; scope `global`
88. `analysis/integration/test_enrichment_module_concordance.R` - stage `integration`; scope `global`
89. `analysis/integration/build_display_selection_inventories.R` - stage `integration`; scope `global`
90. `analysis/integration/plot_display_selection_context.R` - stage `integration`; scope `global`
91. `analysis/integration/summarize_enrichment_module_concordance.R` - stage `integration`; scope `global`
92. `analysis/integration/build_candidate_protein_shortlist.R` - stage `integration`; scope `per_dataset_and_global`
93. `analysis/integration/quantify_candidate_network_position.R` - stage `integration`; scope `per_dataset_and_global`
94. `analysis/integration/build_immunostaining_candidates.R` - stage `integration`; scope `global`
95. `analysis/integration/screen_immunostaining_panel.R` - stage `integration`; scope `global`
96. `analysis/integration/screen_immunostaining_separation.R` - stage `integration`; scope `global`
97. `analysis/integration/render_module_circular_atlas.R` - stage `integration`; scope `global`
98. `analysis/integration/summarize_module_cross_compartment.R` - stage `integration`; scope `global`
99. `analysis/integration/export_module_protein_zoom_source_data.R` - stage `integration`; scope `global`
100. `analysis/publication_source_data/build_sample_metadata.R` - stage `export`; scope `global`
101. `analysis/publication_source_data/export_processed_matrices.R` - stage `export`; scope `global`
102. `analysis/publication_source_data/build_supplementary_tables.R` - stage `export`; scope `global`
103. `analysis/publication_source_data/build_pride_manifest.R` - stage `export`; scope `global`
104. `analysis/publication_source_data/build_methods_summary.R` - stage `export`; scope `global`
105. `analysis/publication_source_data/08_export_manuscript_figures.R` - stage `export`; scope `global`
106. `analysis/publication_source_data/09_export_source_data.R` - stage `export`; scope `global`
107. `analysis/publication_source_data/validate_pride_submission.R` - stage `export`; scope `global`

<!-- END GENERATED PIPELINE REGISTRY INDEX -->

Dry-run the full manuscript path first:

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
done
Rscript run_dataset_pipeline.R --dataset all --stage integration --dry-run
Rscript run_dataset_pipeline.R --dataset all --stage export --dry-run
```

## Core

### Animal-level ProTigy input preparation (external/manual boundary)

Prepare new dataset-specific ProTigy inputs with `AnimalID` as the biological
replicate. The script consumes the newest matching workbook produced by
`archive/01_preprocessing/01_impute.r` for each dataset plus
`data/metadata/TPE9_sample_metadata_males.xlsx`. It applies the existing
exclusions and aggregates available valid hemispheres within each
`AnimalID x canonical spatial unit`. Complete pairs use an equal-weight L/R mean
on the existing imputed log2 scale. If only one valid hemisphere remains because
the other is absent or excluded, the observed value is retained without
hemisphere imputation. It does not filter, transform, normalize, impute, or remap
proteins.

```powershell
Rscript analysis/preprocessing/build_animal_level_protigy_input.R --dataset neuron_neuropil --dry-run
Rscript analysis/preprocessing/build_animal_level_protigy_input.R --dataset neuron_soma --dry-run
Rscript analysis/preprocessing/build_animal_level_protigy_input.R --dataset microglia --dry-run
Rscript analysis/preprocessing/build_animal_level_protigy_input.R --dataset all
```

Canonical quantitative inputs resolve as the newest file matching
`^\d{8}_pgmatrix_imputed_<dataset>_[0-9]+samples_missing70pct\.xlsx$` under
`data/processed/01_preprocessing/impute/`. Optional dedicated overrides are
`PROTEOMICS_PROTIGY_INPUT_EXPRESSION_XLSX`,
`PROTEOMICS_PROTIGY_INPUT_METADATA_XLSX`, and the dataset-specific expression
variants ending in `_NEURON_NEUROPIL`, `_NEURON_SOMA`, or `_MICROGLIA`.

Primary outputs:

```text
data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct
data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.xlsx
data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level_complete_bilateral_sensitivity.gct
results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/<dataset>/
results/logs/01_preprocessing/02a_prepare_animal_level_protigy_input/<dataset>/
```

The primary GCT contains all animal/spatial units with one or two valid observed
hemispheres. The complete-bilateral sensitivity GCT is optional and contains only
units with one valid Left and one valid Right source sample. The script still
fails closed on duplicate same-side hemispheres, hemisphere-label conflicts,
dataset mismatches, sample reuse, or invalid E9 animal/group balance. Failure
audits are retained, but an invalid dataset does not receive a new GCT/XLSX
handoff.

Each ProTigy-targeted GCT deliberately writes `id`, one explicit row descriptor
named `Description`, and then the sample columns. `Description` is populated from
the dataset-specific `First.Protein.Description` field when available and falls
back to `id` only when an independent description is absent or blank. Although
GCT v1.3 can structurally permit zero row descriptors, the current Broad
ProTigy/cmapR downstream workflow expects row descriptor names to align with the
matrix row names; the explicit descriptor is therefore a ProTigy compatibility
contract. Column-metadata rows use `na` in the descriptor filler cell.

The stage is intentionally not registered in `pipeline.yml`: ProTigy remains an
external/manual analysis boundary, and the active `core` stage still consumes
the historical files under `protigy_output/<dataset>/`. After manually running
ProTigy with the new GCTs, a separately authorized migration is required before
`extract_protigy_contrasts.R` or the normal `--stage all` path may consume corrected
ProTigy outputs.

After the manual animal-level ProTigy run, validate the six statistical-result
GCTs and compare DA results directly without mapping or enrichment:

```powershell
Rscript archive/01_preprocessing/03c_legacy_vs_animal_level_da_audit.r --dataset all
```

To create the isolated corrected extraction handoff, set both roots explicitly.
This does not change the historical defaults or the canonical pipeline registry:

```powershell
$env:PROTEOMICS_GCT_INPUT_ROOT = "data/processed/01_preprocessing/protigy_output_animal_level"
$env:PROTEOMICS_GCT_OUTPUT_ROOT = "data/processed/01_preprocessing/gct_extractR_animal_level"
foreach ($dataset in @("neuron_neuropil", "neuron_soma", "microglia")) {
  Rscript analysis/preprocessing/extract_protigy_contrasts.R --dataset $dataset
}
Remove-Item Env:PROTEOMICS_GCT_INPUT_ROOT
Remove-Item Env:PROTEOMICS_GCT_OUTPUT_ROOT
```

The extractor accepts historical `.over.` and corrected `_over_` comparison
syntax. Corrected animal-level roots are strict: only within-unit 2/1, 3/2, and
3/1 comparisons are accepted. Outputs remain isolated under
`data/processed/01_preprocessing/gct_extractR_animal_level/<dataset>/`.

Direct script commands:

Stage 05 Phase 2B remains the quantitative production boundary. The atomic
Phase 3 source migration is complete: Stage 07 publishes
`WGCNA_inferential_handoff.csv`, and claim-facing Stage 08-13, behavior,
integration, claims, evidence-bundle, and publication consumers use that
handoff. Existing generated Stage 05 migration-status CSVs still report the
pre-migration advisory state until an authorized Stage 05 output refresh; they
do not require a WGCNA rebuild or a Stage 05 rerun for this source migration.
The v5 Stage 05 contract adds explicit diagnostic scope and the
cross-platform `sha256_utf8_lf_v1` aggregate hash serialization; it does not
change the statistical analysis.

```bash
Rscript analysis/preprocessing/extract_protigy_contrasts.R --dataset <dataset> --dry-run
Rscript analysis/preprocessing/build_module_score_metadata.R --dataset <dataset> --dry-run
Rscript analysis/preprocessing/map_protein_identifiers.R --dataset <dataset> --dry-run
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
Rscript analysis/qc/build_reference_marker_registry.R --dry-run
Rscript analysis/qc/discover_empirical_roi_markers.R --dry-run
Rscript analysis/qc/assess_dataset_quality.R --dataset <dataset> --dry-run
Rscript analysis/qc/assess_sample_quality.R --dataset <dataset> --dry-run
Rscript analysis/qc/summarize_missingness.R --dataset <dataset> --dry-run
Rscript analysis/qc/assess_replicate_consistency.R --dataset <dataset> --dry-run
Rscript analysis/qc/assess_marker_rank_abundance.R --dataset <dataset> --dry-run
Rscript analysis/qc/summarize_marker_detectability.R --dataset <dataset> --dry-run
Rscript analysis/qc/assess_pca_confounding.R --dataset <dataset> --dry-run
Rscript analysis/qc/partition_variance.R --dataset <dataset> --dry-run
Rscript analysis/qc/export_marker_traits.R --dataset <dataset> --dry-run
Rscript analysis/qc/summarize_qc_confounding.R --dataset <dataset> --dry-run
Rscript analysis/qc/render_compartment_abundance_figures.R --dataset global --dry-run
Rscript analysis/qc/render_compartment_abundance_figures.R --dataset global --render-only --dry-run
```

`04e` is the authoritative CON-only cross-compartment marker-abundance and
detection workflow. Its v2 analysis reconstructs observed, non-imputed
abundance from the raw quantitative columns in the joint bundle and reuses the
joint shared-core normalization offsets. It preserves hemisphere before
animal-level aggregation and treats joint-shared-core membership as a named
sensitivity, not as the primary marker gate.

`04c` remains active for dataset-specific canonical mapping, processed-matrix
availability, WGCNA bridges, and broad annotation. Its nonmissing processed
values are post-filter/post-imputation availability rather than raw
detectability. `04d` is deprecated because its historical sample-level
cross-compartment tests do not respect the animal hierarchy; it is disabled
unless `PROTEOMICS_ENABLE_LEGACY_04D_COMPARTMENT_FIDELITY=true` is set
explicitly.

Completed v2 analytical files are overwrite-protected. Rendering can be
regenerated without recalculating analysis only by explicitly setting
`PROTEOMICS_CONTROL_ABUNDANCE_V2_RENDER_ALLOW_OVERWRITE=true` and using
`--render-only`.

Key outputs:

```text
config/marker_panels/wgcna_reference_marker_sets.csv
results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/
results/reports/03_qc_exploration/00_dataset_qc_report/<dataset>/
results/reports/03_qc_exploration/07_qc_biology_confounding_report/<dataset>/
results/source_data/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/v2_*.csv
results/figures/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/*_v2_*.{svg,pdf,png}
results/reports/03_qc_exploration/04e_control_compartment_abundance_publication_figures/global/v2_*.md
```

## Enrichment

Direct script commands:

```bash
Rscript analysis/differential_abundance/run_clusterprofiler_enrichment.R --dataset <dataset> --dry-run
Rscript analysis/differential_abundance/compare_go_enrichment.R --dataset <dataset> --dry-run
Rscript analysis/differential_abundance/annotate_neuropil_reference.R --dataset <dataset> --dry-run
Rscript analysis/differential_abundance/test_microglia_targeted_signatures.R --dataset microglia --dry-run
Rscript analysis/differential_abundance/summarize_biological_programs.R --dataset <dataset> --dry-run
Rscript analysis/differential_abundance/build_go_program_atlas.R --dataset <dataset> --dry-run
Rscript analysis/differential_abundance/compare_external_stress_signatures.R --dry-run
Rscript analysis/enrichment/run_ewce_celltype_enrichment.R --dataset <dataset> --dry-run
```

Canonical EWCE uses animal-level biological units by default: the two hemispheres are
aggregated within each `AnimalID` x spatial unit, yielding biological `n = 3` animals
per condition. The default command writes `EWCE_E9/<dataset>`. Sample-level EWCE is
retained only as an explicit legacy/sensitivity analysis and requires both
`PROTEOMICS_EWCE_ANALYSIS_UNIT=sample` and `PROTEOMICS_EWCE_BRANCH=<branch>`, which
writes `EWCE_E9_comparison/<branch>/<dataset>`.

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
Rscript analysis/wgcna/build_wgcna_modules.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/compare_recurrent_module_proteins.R --dry-run
Rscript analysis/wgcna/build_curated_overlap_programs.R --dry-run
Rscript analysis/wgcna/score_module_activity.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/compare_module_enrichment_overlap.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/test_module_phenotypes.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/annotate_module_microenvironment.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/summarize_module_interpretation.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/summarize_module_scores.R --dataset <dataset> --module-source wgcna --dry-run
Rscript analysis/wgcna/test_microglia_neuropil_independence.R --dataset microglia --dry-run
Rscript analysis/wgcna/summarize_module_complex_architecture.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/audit_module_robustness.R --dataset <dataset> --dry-run
Rscript analysis/wgcna/summarize_module_interpretation.R --dataset all --dry-run
Rscript analysis/wgcna/summarize_module_scores.R --dataset all --module-source wgcna --dry-run
```

Layer distinction:

```text
analysis/wgcna/build_wgcna_modules.R                         network/module construction
analysis/wgcna/score_module_activity.R         score/statistics/QC producer; secondary robustness/behavior coupling
analysis/wgcna/test_module_phenotypes.R primary WGCNA eigengene group-effect inference
analysis/wgcna/annotate_module_microenvironment.R biological annotation / cleaned semantic label contract
analysis/wgcna/summarize_module_interpretation.R   final WGCNA interpretable module/supermodule summary tables and plots
analysis/wgcna/summarize_module_scores.R final score-derived publication plots using cleaned labels
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
results/figures/06_modules_WGCNA/score_publication_summary/<dataset>/
```

`build_wgcna_modules.R` recomputes core WGCNA state. The downstream WGCNA scripts consume
existing modules and are safe to rerun for reporting/annotation updates.

## Spatial Networks

Direct script commands:

```bash
Rscript analysis/spatial_networks/build_spatial_networks.R --dataset <dataset> --dry-run
Rscript analysis/spatial_networks/build_differential_networks.R --dataset <dataset> --dry-run
Rscript analysis/spatial_networks/test_network_stability.R --dataset <dataset> --dry-run
Rscript analysis/spatial_networks/test_differential_network_stability.R --dataset <dataset> --dry-run
Rscript analysis/spatial_networks/render_differential_network_figures.R --dataset <dataset> --dry-run
Rscript analysis/spatial_networks/render_network_chord_diagram.R --dataset <dataset> --dry-run
```

Key output:

```text
data/processed/07_spatial_networks/network_spatial_relations/<dataset>/*/network_spatial_relations_objects.rds
```

## Coupling

Direct script commands:

```bash
Rscript analysis/integration/test_behaviour_proteomics_associations.R --dataset neuron_soma --dry-run
Rscript analysis/integration/test_network_behaviour_coupling.R --dataset neuron_neuropil --dry-run
Rscript analysis/integration/test_module_behaviour_coupling.R --dataset <dataset> --dry-run
```

Key outputs:

```text
results/tables/08_behavior_physio_coupling/network_behavior_coupling/
results/tables/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/
```

## Integration

Direct script commands:

```bash
Rscript analysis/integration/build_cross_compartment_atlas.R --dry-run
Rscript analysis/integration/summarize_programs_for_manuscript.R --dry-run
Rscript analysis/publication_source_data/build_biological_claims_table.R --dry-run
Rscript analysis/integration/build_evidence_priority_matrix.R --dry-run
Rscript analysis/integration/render_module_circular_atlas.R --dry-run
Rscript analysis/integration/summarize_module_cross_compartment.R --dry-run
Rscript analysis/integration/test_enrichment_module_concordance.R --dry-run
Rscript analysis/integration/summarize_enrichment_module_concordance.R --dry-run
```

Key outputs:

```text
results/tables/10_biological_integration/cross_compartment_program_atlas/global/
results/tables/10_biological_integration/manuscript_program_summary/global/
results/tables/10_biological_integration/evidence_priority_matrix/global/
results/tables/10_biological_integration/gsea_wgcna_concordance/global/
results/tables/10_biological_integration/gsea_wgcna_concordance_diagnostics/global/
```

## Export

Direct script commands:

```bash
Rscript analysis/publication_source_data/build_sample_metadata.R --dry-run
Rscript analysis/publication_source_data/export_processed_matrices.R --dry-run
Rscript analysis/publication_source_data/build_supplementary_tables.R --dry-run
Rscript analysis/publication_source_data/build_pride_manifest.R --dry-run
Rscript analysis/publication_source_data/build_methods_summary.R --dry-run
Rscript analysis/publication_source_data/08_export_manuscript_figures.R --dry-run
Rscript analysis/publication_source_data/09_export_source_data.R --dry-run
Rscript analysis/publication_source_data/validate_pride_submission.R --dry-run
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
