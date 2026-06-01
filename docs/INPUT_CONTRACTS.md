# Input Contracts

Input contracts are machine-readable in `inst/schemas/` and summarized here for reviewers. The schema files are authoritative for required columns, allowed dataset values, allowed claim grades, and p-value/FDR ranges.

| object | schema | required columns | expected location | notes |
|---|---|---|---|---|
| sample metadata | `sample_metadata.yml` | `sample_id`, `dataset` | `data/metadata/` | Optional columns include `animal_id`, `region`, `layer`, `celltype_roi`, `condition`, and `raw_file_name`. |
| protein matrix | `protein_matrix.yml` | `protein_id` | `data/processed/01_preprocessing/` or configured source | Rows are proteins; columns are samples after harmonization. Optional identifier columns include `gene_symbol` and `uniprot_id`. |
| mapped contrast | `mapped_contrast.yml` | `gene_id`, `log2FoldChange` | `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/` | Differential abundance handoff to enrichment. Optional p-value columns `pvalue` and `padj` must be in `[0, 1]` when present. |
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | `dataset`, `output_table` | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/` | Dataset-scoped manifest for enrichment outputs. |
| compareGO manifest | `compareGO_manifest.yml` | `dataset`, `output_table` | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/` | Dataset-scoped manifest for term comparison outputs. |
| WGCNA module contract | `wgcna_module_contract.yml` | `dataset`, `ModuleID`, `ProteinID` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/` | Stable module definitions for downstream module activity and evidence summaries. |
| biological claims table | `biological_claims_table.yml` | `claim_id`, `dataset`, `biological_program`, `evidence_type`, `claim_grade`, `primary_evidence`, `orthogonal_support`, `major_limitation`, `safe_interpretation`, `unsafe_overinterpretation` | `results/tables/biological_claims_table.csv` | `claim_grade` must be one of `A`, `B`, `C`, `D`, `X`; `raw_p` and `FDR` must be in `[0, 1]` when present. |

Private raw data are intentionally not required for CI tests. Dry-runs validate paths, stages, and contracts without reading protected vendor files.
