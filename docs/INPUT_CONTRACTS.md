# Input Contracts

Input contracts are machine-readable in `inst/schemas/` and summarized here for reviewers. The schema files are authoritative for required columns, allowed dataset values, allowed claim grades, and p-value/FDR ranges. `pipeline.yml` remains the active source of truth for which scripts consume each input; see `docs/MAINTENANCE.md` before adding or changing a contract.

| object | schema | required columns | expected location | notes |
|---|---|---|---|---|
| sample metadata | `sample_metadata.yml` | `sample_id`, `dataset` | `data/metadata/` | Optional columns include `animal_id`, `region`, `layer`, `celltype_roi`, `condition`, and `raw_file_name`. |
| protein matrix | `protein_matrix.yml` | `protein_id` | `data/processed/01_preprocessing/` or configured source | Rows are proteins; columns are samples after harmonization. Optional identifier columns include `gene_symbol` and `uniprot_id`. |
| mapped contrast | `mapped_contrast.yml` | `ProteinGroupID`, complete member provenance, official mouse `official_gene_symbol`/`official_entrez_id`, annotation status/version/source hashes, claim flags, and mapping status | `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/` | Canonical protein-group contract with `mouse_gene_annotation_v1`. Stable protein identity is resolved first; official mouse gene annotation is resolved accession-first on every member before protein-group concordance and downstream duplicate-gene collapse. Versioned accession and protein-group annotation audits are written under `results/tables/02_id_mapping/<dataset>/gene_annotation/<direction>/`. Legacy `gene_symbol` remains display-only. |
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | identity, `analysis_status`, term count, output and provenance paths, enrichment/gene-annotation versions | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv` | Strict compareGO, direction-audit, and neuropil-reference input. Successful rows require exact collapsed-gene and term-gene provenance files; failed and zero-term analyses remain explicit. Stale contracts are rejected; recursive/latest discovery is not supported. |
| compareGO manifest | `compareGO_manifest.yml` | canonical clusterProfiler identity/status plus compareGO contract and output paths | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv` | Canonical biological-program-summary input. It declares the term comparison, term-gene provenance, and analysis-status tables. No filename-derived identity, recursive result discovery, or UniProt remapping is permitted. |
| enrichment term-gene provenance | `enrichment_term_gene_provenance.yml` | term identity, official gene, `ProteinGroupID`, member accessions, annotation status/eligibility, rank statistic, contract versions | manifest-declared `gsea_go_term_gene_provenance.csv` or `gsea_kegg_term_gene_provenance.csv` | One row per term/core-gene/contributing eligible protein group. `official_entrez_id` is character. |
| WGCNA module contract | `wgcna_module_contract.yml` | `dataset`, `ProteinGroupID`, `ModuleID`, `ModuleColor`, complete member annotation and claim-eligibility fields | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/` | Phase 1B module definitions use `ProteinGroupID` as the quantitative key. `ProteinID` is a deprecated alias of `ProteinGroupID`; representative accessions and display labels are annotation only. |
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
