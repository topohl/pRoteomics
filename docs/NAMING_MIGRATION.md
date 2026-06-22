# Naming Migration

Old compatibility wrappers have been removed from active code folders. Prefer the canonical names below in manuscripts, methods, and automation.

| old name | preferred name | status | rationale |
|---|---|---|---|
| `03_qc_exploration/06_wgcna_marker_trait_export.r` | `03_qc_exploration/07_wgcna_marker_trait_export.r` | renamed; no wrapper | Makes the active QC prefixes unique while preserving result paths. |
| `03_qc_exploration/07_qc_biology_confounding_report.r` | `03_qc_exploration/08_qc_biology_confounding_report.r` | renamed; no wrapper | Keeps the report after marker-trait export in registry order. |
| `04_differential_expression_enrichment/03_biological_program_summary.r` | `04_differential_expression_enrichment/06_biological_program_summary.r` | renamed; no wrapper | Its registry position is after reference and targeted-signature annotation. |
| `04_differential_expression_enrichment/06_compareGO_spatial_program_atlas.r` | `04_differential_expression_enrichment/07_compareGO_spatial_program_atlas.r` | renamed; no wrapper | Preserves monotonic enrichment-stage numbering. |
| `04_differential_expression_enrichment/07_external_stress_disease_signature_overlap.r` | `04_differential_expression_enrichment/08_external_stress_disease_signature_overlap.r` | renamed; no wrapper | Preserves monotonic enrichment-stage numbering. |
| `06_modules_WGCNA/08_microglia_neuropil_independence.R` | `06_modules_WGCNA/09_microglia_neuropil_independence.R` | renamed; no wrapper | Removes an active duplicate `08` prefix. |
| `06_modules_WGCNA/08_module_complex_architecture.r` | `06_modules_WGCNA/10_module_complex_architecture.r` | renamed; no wrapper | Removes an active duplicate `08` prefix. |
| `06_modules_WGCNA/09_module_robustness_sensitivity.r` | `06_modules_WGCNA/11_module_robustness_sensitivity.r` | renamed; no wrapper | Keeps robustness after complex-architecture analysis. |
| `09_export_pride_journal/01_make_pride_manifest.R` | `09_export_pride_journal/05_make_pride_manifest.R` | renamed; no wrapper | The manifest runs after export inputs `02` through `04`. |
| `09_export_pride_journal/05_validate_pride_submission.R` | `09_export_pride_journal/10_validate_pride_submission.R` | renamed; no wrapper | Validation is the final registered export step. |
| `06_modules_WGCNA/03_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | removed wrapper | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/04_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | removed wrapper | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/05_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | removed wrapper | Active-order, task-based module activity wording. |
| `06_modules_WGCNA/05_wgcna_de_gsea_overlap.r` | `06_modules_WGCNA/04_wgcna_de_gsea_overlap.r` | removed wrapper | Clean active module-stage ordering. |
| `06_modules_WGCNA/91_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | removed wrapper | Active-order, task-based module activity wording. |
| `04_differential_expression_enrichment/04_neuropil_contamination_annotation.r` | `04_differential_expression_enrichment/04_neuropil_reference_annotation.r` | removed; implementation moved to canonical name | Avoids "contamination" language; neuropil is a reference annotation, not subtraction. |
| `04_differential_expression_enrichment/01_run_enrichment_clusterprofiler.r` | `04_differential_expression_enrichment/01_clusterProfiler.r` | removed alias wrapper | The registered active script remains `01_clusterProfiler.r`. |
| `04_differential_expression_enrichment/02_compare_enrichment_terms.r` | `04_differential_expression_enrichment/02_compareGO.r` | removed alias wrapper | The registered active script remains `02_compareGO.r`. |
| `04_differential_expression_enrichment/03_summarize_biological_programs.r` | `04_differential_expression_enrichment/06_biological_program_summary.r` | removed alias wrapper | The registered active script is now `06_biological_program_summary.r`. |
| `09_pride_submission/` | `09_export_pride_journal/` plus gitignored `pride_submission/` payload | `09_export_pride_journal/` active; `09_pride_submission/` legacy helper code | Separates code module from generated deposition payload. |

Legacy/deprecated scripts remain in `90_testing/`, `99_deprecated/`, or are explicitly listed in `pipeline.yml:legacy`.
