# Naming Migration

Backward-compatible wrappers or retained scripts keep old commands working. Prefer the new names in manuscripts, methods, and new automation.

| old name | preferred name | status | rationale |
|---|---|---|---|
| `06_modules_WGCNA/03_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | wrapper only | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/04_overlap_modules.r` | `06_modules_WGCNA/02_curated_overlap_programs.r` | wrapper only | Curated overlap-derived biological programs are distinct from WGCNA modules. |
| `06_modules_WGCNA/05_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | wrapper only | Active-order, task-based module activity wording. |
| `06_modules_WGCNA/05_wgcna_de_gsea_overlap.r` | `06_modules_WGCNA/04_wgcna_de_gsea_overlap.r` | wrapper only | Clean active module-stage ordering. |
| `06_modules_WGCNA/91_module_score.r` | `06_modules_WGCNA/03_score_module_activity.R` | wrapper only | Active-order, task-based module activity wording. |
| `04_differential_expression_enrichment/04_neuropil_contamination_annotation.r` | `04_differential_expression_enrichment/04_neuropil_reference_annotation.r` | old implementation retained; new wrapper is active in `pipeline.yml` | Avoids "contamination" language; neuropil is a reference annotation, not subtraction. |
| `04_differential_expression_enrichment/01_clusterProfiler.r` | `04_differential_expression_enrichment/01_run_enrichment_clusterprofiler.r` | preferred alias wrapper available; old pipeline entry retained | Task-based name; wrapper delegates to existing implementation. |
| `04_differential_expression_enrichment/02_compareGO.r` | `04_differential_expression_enrichment/02_compare_enrichment_terms.r` | preferred alias wrapper available; old pipeline entry retained | Tool-independent description of term comparison. |
| `04_differential_expression_enrichment/03_biological_program_summary.r` | `04_differential_expression_enrichment/03_summarize_biological_programs.r` | preferred alias wrapper available; old pipeline entry retained | Task-based manuscript-facing wording. |
| `09_pride_submission/` | `09_export_pride_journal/` plus gitignored `pride_submission/` payload | `09_export_pride_journal/` active; `09_pride_submission/` legacy helper code | Separates code module from generated deposition payload. |

Legacy/deprecated scripts remain in `90_testing/`, `99_deprecated/`, or are explicitly listed in `pipeline.yml:legacy`.
