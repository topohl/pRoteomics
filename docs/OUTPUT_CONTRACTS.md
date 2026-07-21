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

## Recommended Biological Results Entry Point

Start manuscript-facing biological review with:

```text
results/tables/10_biological_integration/final_evidence_bundle/global/final_biological_evidence_bundle.xlsx
```

This workbook is a final evidence bundle assembled above the existing result tables. It does not recompute upstream statistics, WGCNA construction, model formulas, or thresholds. The `README` sheet explains each sheet, its upstream producer, and which columns are safest for manuscript interpretation. CSV mirrors in the same directory and under `results/source_data/10_biological_integration/final_evidence_bundle/global/` support traceable review.

Final-facing dataset terminology in this bundle is:

- `neuron_neuropil`: region/layer-resolved neuron neuropil
- `neuron_soma`: region-resolved neuronal soma-enriched
- `microglia`: region-resolved microglia-enriched ROI/local microenvironment, not purified microglia

## Biological Claim Gate

`results/tables/biological_claims_table.csv` is a reviewer-facing claim gate, not only a claim index. `claim_grade` remains a descriptive evidence label for backwards compatibility. `claim_type` controls which gates are required, optional, diagnostic-only, or not applicable. Manuscript eligibility is determined by `claim_allowed` and `claim_gate_status`.

Model/statistical semantics are separated for review. `model_fit_status` records fit diagnostics (`pass`, `singular`, `rank_deficient`, `fallback`, `failed`, `not_applicable`). `statistical_evidence_status` records the p/FDR layer (`pass`, `nominal_only`, `FDR_fail`, `missing_p_or_FDR`, `not_applicable`). `claim_gate_model_status` is the model/statistical claim-gate result. `primary_model_status` is retained for backward compatibility: `diagnostic_only` means the model fit is not claim-grade, `fail` means the fit was usable but statistical evidence did not pass, and `pass` means both model eligibility and statistical evidence passed.

Gate columns use `pass`, `fail`, `missing_required`, `missing_optional`, `not_applicable`, and `diagnostic_only`. Required gate failures or `missing_required` disallow a claim. `missing_optional` and `diagnostic_only` restrict or downgrade wording but are not universal blockers. Marker-contamination, microglia ROI, and neuropil-independence gates are required for microglia-specific, purified/cell-intrinsic, and contamination-sensitive microglia module/supermodule claims; they are not required for every neuron_neuropil or neuron_soma claim. Microglia dataset claims default to "microglia-enriched ROI/local microenvironment" wording.

Microglia-neuropil independence is audited in
`results/reviewer_audit/microglia_neuropil_independence_claim_gate.csv`.
Rows with `adjustment_mode` `predeclared_primary` or
`predeclared_secondary` are selected from
`config/microglia_neuropil_independence.yml` before looking at endpoint
correlation and may support `claim_gate_eligible = TRUE` only when a
claim-relevant primary microglia effect exists first. The audit table reports
`primary_effect_status`, `primary_effect_claim_relevant`,
`primary_effect_threshold`, `primary_effect_p`, and `primary_effect_FDR`.
Rows with `diagnostic_no_primary_effect` or `inconclusive_no_primary_effect`
do not produce a `neuropil_independence_gate=pass`.
`exploratory_best_spearman` rows preserve the previous best-correlation
diagnostic, but they are circular for claim selection and therefore cannot
enable microglia-specific or cell-intrinsic wording by themselves.
`neuropil_sensitive` blocks or downgrades microglia-specific interpretation;
inconclusive classifications are diagnostic/contextual unless other required
predeclared evidence passes.

Every independence endpoint reports `endpoint_scope`, `source_level`, and
`direct_independence_tested`. Module-level diagnostics may be consumed by
`WGCNA_module` and `WGCNA_module_group_effect` rows only; they do not certify
`WGCNA_supermodule_group_effect` rows. Without a direct
`supermodule_eigengene` or `supermodule_score` endpoint, the supermodule gate
is `missing_required` with reason
`no_direct_supermodule_independence_test`. Endpoint availability and claim
consumption are summarized in
`results/reviewer_audit/microglia_neuropil_independence_endpoint_scope_audit.csv`.

Percent attenuation is guarded against near-zero baseline effects using
`effect_before_abs`, `effect_before_near_zero`, and
`percent_attenuation_reliable`. If attenuation is unstable because the baseline
effect is near zero, the downgrade reason records
`unstable_attenuation_near_zero_effect` and attenuation is not used as
claim-enabling evidence. Module-level independence is aggregated
conservatively: no claim-relevant primary contrast yields a
diagnostic/inconclusive no-primary classification, and disagreement across
predeclared covariates is reported as `mixed_or_covariate_sensitive` rather
than promoted by best-case matching.

For GO-derived enrichment claims, raw statistical evidence is separated from manuscript-safe label interpretation. `raw_top_GO_term` and `representative_GO_terms` preserve the original GO evidence. `term_label_risk`, `semantic_parent_label`, `safe_program_label`, `label_confidence`, and `label_basis` describe conservative relabeling when a literal GO label is tissue-mismatched or low-specificity. Tissue-mismatched GO terms are not hard-blacklisted; they are flagged, conservatively relabeled where the term cluster or leading-edge genes support a broader process, and assigned a `claim_use_class` (`primary_claim`, `supporting_claim`, `suggestive_context`, `annotation_only`, or `blocked`). For GO-label-risk rows, `safe_interpretation` must use `safe_program_label`, preserve `raw_top_GO_term`, and warn against literal interpretation of the raw GO term. `annotation_only` rows are retained for audit transparency, not manuscript claims; `suggestive_context` rows are contextual only.

Blocked rows must not use positive manuscript-claim wording. Their `safe_interpretation` starts with `Not claim-eligible` and includes gate reasons. Reviewer audit tables are written to `results/reviewer_audit/claim_gate_evidence_availability.csv`, `results/reviewer_audit/claim_gate_summary.csv`, `results/reviewer_audit/go_label_risk_audit.csv`, `results/reviewer_audit/go_label_safe_interpretation_audit.csv`, `results/reviewer_audit/claim_use_class_summary.csv`, `results/reviewer_audit/blocked_claim_wording_audit.csv`, `results/reviewer_audit/wgcna_claim_source_audit.csv`, `results/reviewer_audit/wgcna_label_completeness_audit.csv`, `results/reviewer_audit/wgcna_label_confidence_audit.csv`, `results/reviewer_audit/wgcna_annotation_source_audit.csv`, `results/reviewer_audit/wgcna_microenvironment_threshold_sensitivity.csv`, `results/reviewer_audit/microglia_neuropil_independence_claim_gate.csv`, and `results/reviewer_audit/microglia_neuropil_covariate_selection_audit.csv`.

`claim_allowed` is the gate-level eligibility decision; `claim_use_class` is the
permitted narrative role after label confidence and evidence context are
considered. `supporting_claim` may use cautious positive wording,
`suggestive_context` must identify itself as contextual and not standalone,
`annotation_only` is retained for label provenance/audit transparency, and
`blocked` starts with `Not claim-eligible`. Microglia rows retain ROI/local
microenvironment wording because the samples are not purified microglia.
Supermodule neuropil independence remains `missing_required` unless a direct
supermodule eigengene or score endpoint was tested.

The final reviewer index is
`results/reviewer_audit/final_reviewer_audit_manifest.csv`. Its companion
`final_evidence_bundle_validation.csv` reports machine-readable PASS/FAIL
checks, while `final_claim_gate_summary.csv` provides the fully stratified
claim-gate counts. Both are also included in the final evidence workbook.

## Input Resolution Audit

Claim-critical and manuscript-facing runs append input provenance to:

```text
results/reviewer_audit/input_resolution_audit.csv
```

The table records `script`, `dataset`, `stage`, `input_name`, expected and resolved paths, resolution mode, strict-mode policy, file existence, SHA-256 hash, modification time, upstream producer/artifact when available, and any fallback warning. In strict input mode (`--strict-inputs` or `PROTEOMICS_STRICT_INPUTS=true`), newest-file/latest-file fallbacks are not allowed for these reviewer-facing paths; missing required canonical inputs fail instead of being silently replaced.

## Required Final Tables

| output | schema or contract | expected location |
|---|---|---|
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv` |
| compareGO input manifest | `compareGO_manifest.yml` | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv` |
| clusterProfiler term-gene provenance | `enrichment_term_gene_provenance.yml` | manifest-declared per-comparison `protein_group_audits/gsea_{go,kegg}_term_gene_provenance.csv` |
| compareGO canonical tables | manifest-declared term and provenance contracts | `results/tables/04_differential_expression_enrichment/compareGO/<dataset>/<ontology>/<route_category>/<route_unit>/compareGO_{term_comparison,term_gene_provenance,analysis_status_summary}.csv` |
| GSEA direction audit | canonical downstream diagnostic contract | `results/tables/04_differential_expression_enrichment/01b_gsea_protein_direction_audit/<dataset>/<ontology>/gsea_{term_direction_audit,contrast_direction_summary}.csv` |
| neuropil reference annotation | canonical downstream annotation contract | `results/tables/04_differential_expression_enrichment/neuropil_reference_annotation/<dataset>/microglia_neuropil_{annotation_latest,analysis_status}.csv` |
| biological program summary | documented in `docs/file_contracts.tsv` | `results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv` |
| biological program gene evidence/status | canonical downstream provenance contract | `results/source_data/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_term_gene_evidence.csv`; `results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_analysis_status.csv` |
| targeted microglia ROI signature claims | documented in `docs/file_contracts.tsv` | `results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia/microglia_signature_claims_ready.csv` |
| module-score merged metadata | module-score metadata merge contract | `data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx`; QC/report/log sidecars under `results/{tables,reports,logs}/01_preprocessing/06_merged_metadata_module_score/<dataset>/` |
| WGCNA module definitions | `wgcna_module_contract.yml` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_definitions_for_downstream.csv` |
| curated overlap programs | curated overlap program contract | `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs.xlsx`; `results/tables/06_modules_WGCNA/curated_overlap_programs/global/curated_overlap_programs_long.csv` |
| module activity scores | module score contract | `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_scores_per_sample.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_feature_mapping_trace.csv`; `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/module_gene_coverage.csv`; `results/logs/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/run_manifest.yml` |
| WGCNA/DE/GSEA overlap | overlap bridge contract | `results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/WGCNA_vs_DE_GSEA_overlap.csv` |
| Reference marker evidence registry | marker registry contract | `results/tables/03_qc_exploration/reference_marker_import/reference_marker_evidence_long.csv`; `results/tables/03_qc_exploration/reference_marker_import/reference_marker_sets_ranked.csv`; `config/marker_panels/wgcna_reference_marker_sets.csv` |
| Canonical marker-QC feature universe | `protein_group_id_v1` plus ordered fingerprints | Each active gene-aware QC result directory contains `feature_universe_reconciliation.csv`, `canonical_feature_table.csv`, `protein_group_member_bridge.csv`, `source_row_provenance.csv`, `protein_group_exclusion_audit.csv`, `multi_member_group_audit.csv`, `input_path_hash_manifest.csv`, `feature_contract.csv`, `marker_to_ProteinGroupID_provenance.csv`, `marker_to_ProteinGroupID_member_provenance.csv`, and `conflicting_marker_panel_audit.csv`. The expanded member-provenance artifact preserves every matched accession/gene row before panel-feature collapse. The optional 04c bridge also writes `wgcna_feature_contract_validation.csv`. |
| Supplemental WGCNA microenvironment marker panels | externalized local annotation contract | `config/marker_panels/microenvironment_marker_panels.csv` |
| Empirical ROI marker discovery | empirical marker contract | `results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv`; `results/source_data/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv` |
| WGCNA marker traits | annotation-trait contract | `results/tables/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv`; `results/source_data/03_qc_exploration/06_wgcna_marker_trait_export/<dataset>/wgcna_marker_traits_by_sample.csv` |
| WGCNA module/supermodule group effects | `validate_wgcna_group_effects()` | `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_marker_trait_correlations.csv`; `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_marker_trait_correlations.csv` |
| WGCNA eigengene meta-modules | stable `dataset + SupermoduleID` member/summary contract | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_module_supermodule_annotation.csv`; `wgcna_supermodule_summary.csv`; `wgcna_supermodule_validation_summary.csv` |
| WGCNA recurring-GO support audit | `ModuleProteinSetType=all`; member-module `p.adjust <= 0.05` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_supermodule_GO_term_support_audit.csv` |
| WGCNA supermodule biological coherence | signed eigengene/coherence contract | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_supermodule_biological_coherence.csv`; `supermodule_clustering_sensitivity.csv` |
| WGCNA biological annotation | `validate_wgcna_module_annotation()` | `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`; `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv` |
| WGCNA interpretable summary | `validate_wgcna_interpretable_summary()`; automatic proposals in `wgcna_label_candidates.yml`; reviewed microglia lookup in `wgcna_reviewed_final_label_lookup.yml` | `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_module_group_effects_interpretable.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_label_candidates.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_final_label_lookup.csv`; `results/tables/06_modules_WGCNA/interpretable_summary/all/WGCNA_cross_dataset_supermodule_program_summary.csv` |
| Reviewed microglia WGCNA publication figures | reviewed-label and all-supermodule source-data validation | `results/figures/06_modules_WGCNA/wgcna_publication_figures/microglia/{all_supermodule_architecture,all_supermodule_global_eigengenes,all_supermodule_sus_res_forest,all_supermodule_sus_res_spatial_heatmap,multi_supermodule_member_loadings,supplementary_all_contrast_effect_matrix}.{svg,pdf}`; matching source CSVs under `results/source_data/06_modules_WGCNA/wgcna_publication_figures/microglia/` |
| WGCNA score-derived publication summary | score-publication validation table | `results/figures/06_modules_WGCNA/score_publication_summary/<dataset>/publication_supermodule_effect_heatmap_<dataset>_wgcna_{primary_all_replicates,sensitivity}.svg`; `results/figures/06_modules_WGCNA/score_publication_summary/<dataset>/publication_supermodule_correlation_<dataset>_wgcna_{primary_all_replicates,sensitivity}.svg`; `results/figures/06_modules_WGCNA/score_publication_summary/<dataset>/publication_supermodule_consistency_<dataset>_wgcna_{primary_all_replicates,sensitivity}.svg`; `results/tables/06_modules_WGCNA/score_publication_summary/<dataset>/WGCNA_score_publication_validation.csv` |
| Microglia-neuropil independence claim gate | `microglia_neuropil_independence_claim_gate.yml`; `microglia_neuropil_covariate_selection_audit.yml` | `results/reviewer_audit/microglia_neuropil_independence_claim_gate.csv`; `results/reviewer_audit/microglia_neuropil_covariate_selection_audit.csv`; detailed effects under `results/tables/06_modules_WGCNA/microglia_neuropil_independence/microglia/` |

### Supermodule construction and naming

WGCNA supermodules are eigengene meta-modules/co-varying module blocks
constructed by average linkage on `1 - signed module-eigengene correlation`.
Their authoritative key is `dataset + SupermoduleID`; biological and display
labels are never join, grouping, deduplication, or identity keys. Full display
labels retain the stable `SMxx` prefix and collisions fail validation. Because
WGCNA modules partition proteins, protein/hub overlap is only a
partition-integrity diagnostic and never biological support.

Supermodule GO naming uses only `ModuleProteinSetType == "all"`. A member
module supports a term only at `p.adjust <= 0.05`; the audit retains the support
count and fraction, supporting modules, and best/worst supporting member-module
FDR. High confidence requires support in every member module. Medium confidence
requires at least two and at least half. Other multi-module cases are
mixed/unresolved, and automatic Stage 01 singletons are explicitly structural
singletons. Reviewed singleton biological confidence is inherited from the
member module and is not lowered automatically. Biological-coherence exports report signed minimum/mean/median
pairwise eigengene correlation, region/layer-adjusted correlations when
estimable, PC1 variance explained, module/protein counts, and stability over
cut heights `0.25, 0.35, 0.45, 0.50, 0.55, 0.65`; PC1 variance alone is not
shared-pathway evidence.

WGCNA biological annotation and interpretable summaries keep raw GO,
automatic proposals, hubs, marker evidence, member composition, subcellular
context, and ROI context as separate provenance. For microglia, the reviewed
source registry is `config/wgcna_labels/microglia.csv`; Stage 07 is the only
authority that promotes those labels into the canonical lookup. Stage 05 never
reads Stage 07 outputs and Stage 06 cannot overwrite reviewed labels. Singleton
supermodules inherit the exact reviewed member-module biological label, short
label, confidence, and manual-review state. Singleton structure remains
metadata and does not lower biological confidence.

`WGCNA_final_label_lookup.csv` is keyed by `dataset + level + entity_id`.
Microglia rows expose `canonical_biological_label`, `canonical_short_label`,
`canonical_plot_label`, separate subcellular/ROI contexts, current membership
and annotation-evidence fingerprints, structural metadata, and label
provenance. `final_plot_label` is retained only as an exact compatibility alias
for `canonical_plot_label`. Publication rendering fails when a stable ID or
membership fingerprint is missing or stale; plots do not choose labels through
fallback chains across legacy aliases.
| WGCNA module complex/organelle architecture | standard integration evidence contract | `results/tables/06_modules_WGCNA/module_complex_architecture/<dataset>/module_complex_architecture.csv`; mirrored under `results/source_data/06_modules_WGCNA/module_complex_architecture/<dataset>/` |
| WGCNA module robustness/sensitivity | standard integration evidence contract | `results/tables/06_modules_WGCNA/module_robustness_sensitivity/<dataset>/module_robustness_sensitivity.csv`; mirrored under `results/source_data/06_modules_WGCNA/module_robustness_sensitivity/<dataset>/` |
| module-behavior coupling | standard integration evidence contract | `results/tables/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/module_behavior_coupling.csv`; mirrored under `results/source_data/08_behavior_physio_coupling/module_behavior_coupling/<dataset>/` |
| cross-compartment program atlas | `validate_cross_compartment_program_atlas()` | `results/tables/10_biological_integration/cross_compartment_program_atlas/global/cross_compartment_program_atlas_long.csv`; `results/tables/10_biological_integration/cross_compartment_program_atlas/global/cross_compartment_program_atlas.csv` |
| manuscript program summary | `validate_manuscript_program_summary()` | `results/tables/10_biological_integration/manuscript_program_summary/global/manuscript_program_summary.csv` |
| evidence priority matrix | `validate_evidence_priority_matrix()` | `results/tables/10_biological_integration/evidence_priority_matrix/global/evidence_priority_matrix.csv` |
| final biological evidence bundle | workbook/CSV bundle contract documented in workbook `README` sheet | `results/tables/10_biological_integration/final_evidence_bundle/global/final_biological_evidence_bundle.xlsx`; CSV sheets in the same directory |
| spatial network object | dataset/spatial-unit spatial network contract | `data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/network_spatial_relations_objects.rds`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_mean_expression_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/spatial_unit_spearman_correlation_matrix.csv`; `results/tables/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/standardized_sample_metadata_used.csv`; `results/logs/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/run_manifest.yml` |
| biological claims table | `biological_claims_table.yml` | `results/tables/biological_claims_table.csv` |
| input resolution audit | `input_resolution_audit.yml` | `results/reviewer_audit/input_resolution_audit.csv` |

Legacy compatibility note: existing `results/module_scores/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx` and `results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx` files are read-only fallbacks. In non-strict exploratory mode, fallback use emits warnings and is recorded in `input_resolution_audit.csv`; the global fallback requires `PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA=true`. In strict input mode, these fallbacks are forbidden unless the canonical dataset-scoped input is present or an explicit override is provided.

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
# WGCNA Group-Effect Claim Fields

`results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv` and `supermodule_group_effects.csv` preserve existing estimates and add claim-grade status fields: `model_family`, `model_formula`, `primary_model_stable`, `claim_allowed_model`, `model_downgrade_reason`, fallback flags, rank/singularity/emmeans diagnostics, animal random-effect use, replicate-unit fields, animal/sample counts by group, `animal_level_status`, and `pseudoreplication_guard`.

`primary_model_stable = FALSE` means the row is not claim-grade even if a numeric estimate exists. Fallback rows are diagnostic-only and carry `diagnostic_only_model_fallback` in `model_downgrade_reason`.

`results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_parameter_audit.csv` records WGCNA construction parameters and input hashes/provenance for reviewer audit.

`results/reviewer_audit/wgcna_robustness_claim_gate.csv` records module/supermodule robustness eligibility with preservation, sensitivity, direction-stability, and confounding status. The biological claims table consumes this audit without loosening the strict claim gate.

In `biological_claims_table.csv`, `WGCNA_module_group_effect` and `WGCNA_supermodule_group_effect` rows are direct group-effect model rows. `WGCNA_module` rows from `WGCNA_module_evidence_rank.csv` are ranked evidence summaries and remain `suggestive_context` unless all required claim gates pass; they should not be treated as the primary group-effect model rows.

`results/reviewer_audit/wgcna_label_completeness_audit.csv` scans `biological_program`, `safe_program_label`, and `safe_interpretation` for incomplete WGCNA reviewer-facing labels such as `shared ROI: mostly` or trailing `mixed:`. Expected reviewer-ready exports have zero rows in this audit; unrepaired incomplete labels are downgraded to annotation-only wording.

`results/reviewer_audit/wgcna_annotation_source_audit.csv` records the externalized supplemental marker-panel config, source/use metadata, marker counts, and config hash. `results/reviewer_audit/wgcna_microenvironment_threshold_sensitivity.csv` reports module and supermodule annotation classes at marker-fraction thresholds `0.05`, `0.10`, and `0.20`. `results/reviewer_audit/wgcna_label_confidence_audit.csv` records safe display labels, confidence, basis, downgrade reasons, threshold stability, and unsafe interpretation warnings. Threshold-unstable or low-confidence WGCNA labels remain annotation/context only unless the independent claim gates pass.
