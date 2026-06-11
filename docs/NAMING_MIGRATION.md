# Naming Migration

Old compatibility wrappers have been removed from active code folders. Prefer the canonical names below in manuscripts, methods, and automation.

| old name | preferred name | status | rationale |
|---|---|---|---|
| `06_modules_WGCNA/03_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | removed wrapper | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/04_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | removed wrapper | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/05_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | removed wrapper | Active-order, task-based module activity wording. |
| `06_modules_WGCNA/05_wgcna_de_gsea_overlap.r` | `06_modules_WGCNA/04_wgcna_de_gsea_overlap.r` | removed wrapper | Clean active module-stage ordering. |
| `06_modules_WGCNA/91_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | removed wrapper | Active-order, task-based module activity wording. |
| `04_differential_expression_enrichment/04_neuropil_contamination_annotation.r` | `04_differential_expression_enrichment/04_neuropil_reference_annotation.r` | removed; implementation moved to canonical name | Avoids "contamination" language; neuropil is a reference annotation, not subtraction. |
| `04_differential_expression_enrichment/01_run_enrichment_clusterprofiler.r` | `04_differential_expression_enrichment/01_clusterProfiler.r` | removed alias wrapper | The registered active script remains `01_clusterProfiler.r`. |
| `04_differential_expression_enrichment/02_compare_enrichment_terms.r` | `04_differential_expression_enrichment/02_compareGO.r` | removed alias wrapper | The registered active script remains `02_compareGO.r`. |
| `04_differential_expression_enrichment/03_summarize_biological_programs.r` | `04_differential_expression_enrichment/03_biological_program_summary.r` | removed alias wrapper | The registered active script remains `03_biological_program_summary.r`. |
| `09_pride_submission/` | `09_export_pride_journal/` plus gitignored `pride_submission/` payload | `09_export_pride_journal/` active; `09_pride_submission/` legacy helper code | Separates code module from generated deposition payload. |

Legacy/deprecated scripts remain in `90_testing/`, `99_deprecated/`, or are explicitly listed in `pipeline.yml:legacy`.
