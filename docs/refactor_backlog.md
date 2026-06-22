# Refactor and Scientific-Guardrail Backlog

This backlog records follow-up improvements identified after the QC/exploration
refactor and repository-wide script scan. These are not required to run the
current pipeline, but they would make the analysis more reproducible and make
stress/social-instability biological claims harder to confound with technical
structure.

## Completed: QC Is In the Main Dataset Pipeline

`pipeline.yml` now includes `qc_global` and dataset-specific `qc` stages between
`core` and `enrichment`; `run_dataset_pipeline.R` discovers these stage names
directly from the registry.

Recommended stage contents:

```text
03_qc_exploration/01_sample_qc_quicksearch.r
03_qc_exploration/02_missingness_diagnostics.r
03_qc_exploration/03_replicate_consistency.r
03_qc_exploration/04_marker_rank_abundance_qc.r
03_qc_exploration/05_pca_confounding_qc.r
03_qc_exploration/06_variance_partitioning.r
03_qc_exploration/08_qc_biology_confounding_report.r
```

QC/confounding checks are now normal registry outputs rather than manual side
analyses. This remains especially important for missingness, plate/batch effects,
region/layer imbalance, animal-level pseudoreplication, marker contamination,
and technical PCs.

## Priority 2: Feed QC Flags Into Biological Claims

Extend `09_export_pride_journal/07_make_biological_claims_table.R` so each
claim can carry QC/confounding context from:

```text
results/tables/03_qc_exploration/02_missingness_diagnostics/<dataset>/
results/tables/03_qc_exploration/03_replicate_consistency/<dataset>/
results/tables/03_qc_exploration/07_qc_biology_confounding_report/<dataset>/
results/tables/03_qc_exploration/05_pca_confounding_qc/<dataset>/PCA_confounding_summary.csv
```

Suggested added claim columns:

```text
missingness_confounded
plate_or_batch_confounded
region_layer_imbalance_risk
animal_pseudoreplication_risk
early_pc_association
marker_contamination_risk
qc_interpretation_flag
```

Rationale: the claims table should not only summarize positive biological
evidence. It should also show whether that evidence is vulnerable to QC or design
confounding.

## Priority 3: Expand Stress/Social-Instability Program Mapping

Update `R/enrichment_io.R::biological_program_patterns()` so the program summary
layer better matches the stress/social-instability interpretation guide.

Candidate additions or splits:

```text
HPA_Glucocorticoid_Response
Neuroimmune_Complement_Phagosome
ECM_Vascular_Barrier
Oxidative_Redox_Stress
Lipid_Myelin_Membrane
Autophagy_Lysosome
```

Keep this as an interpretation layer, not a new statistical test. Program names
should remain broad enough to avoid false precision from regex-only mapping.

## Priority 4: Add Direction-Consistency to Program Summaries

Extend `04_differential_expression_enrichment/06_biological_program_summary.r`
with direction-consistency fields:

```text
n_positive_terms
n_negative_terms
median_NES
direction_consistency
strongest_positive_term
strongest_negative_term
```

Rationale: a program can look important because of a low-FDR representative term
while related terms point in mixed directions. Direction-consistency columns
would make thematic summaries safer for manuscript interpretation.

## Priority 5: Make Animal-Level Evidence Explicit

Across downstream summaries and claim-building, distinguish:

```text
sample_level
animal_level
animal_aggregated
mixed_or_unclear
```

Highest-value targets:

```text
04_differential_expression_enrichment summaries
06_modules_WGCNA/03_score_module_activity.R
08_behavior_physio_coupling scripts
09_export_pride_journal/07_make_biological_claims_table.R
```

Rationale: the biological unit is animal, not individual tissue punch/replicate.
Claims should not accidentally borrow strength from pseudoreplicates.

## Priority 6: Resolve or Archive Partially Canonical Active Scripts

These scripts still deserve a focused pass:

- `01_preprocessing/06_merged_metadata_module_score.r`
  - Contains duplicated resolver logic and date-stamped/default legacy input
    assumptions.
  - Consider replacing resolver code with `R/dataset_inputs.R` helpers or
    archiving if superseded.

- `04_differential_expression_enrichment/legacy/04_compare_pathways.r`
  - Clarify whether it is canonical, optional, or legacy.
  - If retained, add dataset-aware inputs, dry-run, manifest, and documented
    output contracts.

- `04_differential_expression_enrichment/legacy/05_compare_sig_expr.r`
  - Contains legacy path hints and older compareGO assumptions.
  - Either make it manifest-driven or move to a legacy folder.

- `04_differential_expression_enrichment/legacy/07_control_strata_enrichment_figures.r`
  - Potentially useful figure logic, but it needs canonical input/output
    contracts and dry-run behavior.

- `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`
  - Important behavior/proteomics producer, but still has date-stamped/default
    inputs.
  - Should become dataset-aware and animal-level explicit before manuscript use.

- `09_pride_submission/*`
  - Older PRIDE package scripts use a different package path convention from the
    newer `09_export_pride_journal` layer.
  - Decide whether to retire, redirect, or harmonize them.

## Priority 7: Refresh Repository Audit Docs

After the QC refactor, update:

```text
docs/active_script_io_audit.tsv
docs/file_contracts.tsv
docs/current_data_flow.md
```

Specific stale points to fix:

- `03_qc_exploration/01_sample_qc_quicksearch.r`,
  `04_marker_rank_abundance_qc.r`, `05_pca_confounding_qc.r`, and `06_variance_partitioning.r` are now
  canonical dataset-aware scripts.
- `03_qc_exploration/02_missingness_diagnostics.r`,
  `03_replicate_consistency.r`, and
  `07_qc_biology_confounding_report.r` should be added to the audit/contracts.
- `03_qc_exploration/legacy/06_pcaPlot_Neha.r` and
  `03_qc_exploration/legacy/08_boxplotBonanza.r` should be marked legacy, not
  active canonical scripts.
