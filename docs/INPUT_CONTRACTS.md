# Input Contracts

Input contracts are machine-readable in `inst/schemas/` and summarized here for reviewers. The schema files are authoritative for required columns, allowed dataset values, allowed claim grades, and p-value/FDR ranges. `pipeline.yml` remains the active source of truth for which scripts consume each input; see `docs/MAINTENANCE.md` before adding or changing a contract.

| object | schema | required columns | expected location | notes |
|---|---|---|---|---|
| sample metadata | `sample_metadata.yml` | `sample_id`, `dataset` | `data/metadata/` | Optional columns include `animal_id`, `region`, `layer`, `celltype_roi`, `condition`, and `raw_file_name`. |
| protein matrix | `protein_matrix.yml` | `protein_id` | `data/processed/01_preprocessing/` or configured source | Rows are proteins; columns are samples after harmonization. Optional identifier columns include `gene_symbol` and `uniprot_id`. |
| mapped contrast | `mapped_contrast.yml` | `gene_id`, `log2FoldChange` | `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/` | Differential abundance handoff to enrichment. Optional p-value columns `pvalue` and `padj` must be in `[0, 1]` when present. |
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | `dataset`, `output_table` | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/` | Dataset-scoped manifest for enrichment outputs. |
| compareGO manifest | `compareGO_manifest.yml` | `dataset`, `output_table` | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/` | Dataset-scoped manifest for term comparison outputs. |
| WGCNA module contract | `wgcna_module_contract.yml` | `dataset`, `ModuleID`, `ProteinID` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/` | Stable module definitions for downstream module activity and evidence summaries. |
| biological claims table | `biological_claims_table.yml` | `claim_id`, `dataset`, `biological_program`, `evidence_type`, `claim_type`, `claim_use_class`, `claim_allowed`, `claim_gate_status` | `results/tables/biological_claims_table.csv` | `claim_grade` is descriptive. Manuscript eligibility is controlled by `claim_allowed`, `claim_gate_status`, gate columns, GO label-risk columns, and `claim_use_class`. |
| input resolution audit | `input_resolution_audit.yml` | `script`, `dataset`, `stage`, `input_name`, `expected_path`, `resolved_path`, `resolution_mode`, `strict_mode`, `allowed_in_strict_mode`, `file_exists`, `file_hash_sha256`, `file_mtime`, `producer_script_or_artifact_id`, `warning` | `results/reviewer_audit/input_resolution_audit.csv` | Reviewer provenance ledger for claim-critical and manuscript-facing input resolution. |

Private raw data are intentionally not required for CI tests. Dry-runs validate paths, stages, and contracts without reading protected vendor files.

## Strict Input Mode

Reviewer/manuscript runs should use strict input mode:

```bash
Rscript run_dataset_pipeline.R --dataset all --stage export --dry-run --strict-inputs
PROTEOMICS_STRICT_INPUTS=true Rscript 09_export_pride_journal/07_make_biological_claims_table.R --dry-run
```

Strict mode is enabled by `--strict-inputs` or `PROTEOMICS_STRICT_INPUTS=true`. In strict mode, claim-critical and manuscript-facing scripts must use explicit or canonical inputs; newest-file and legacy fallback selection is forbidden unless a resolver explicitly marks that fallback as strict-safe. In exploratory/non-strict mode, fallbacks may still be used for backwards compatibility, but they emit warnings and append rows to `results/reviewer_audit/input_resolution_audit.csv` with SHA-256 hashes and modification times for resolved files.
