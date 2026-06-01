# Input Contracts

Input contracts are machine-readable in `inst/schemas/` and summarized here for reviewers.

| object | schema | expected location | notes |
|---|---|---|---|
| sample metadata | `sample_metadata.yml` | `data/metadata/` | Must identify samples, dataset family, biological grouping, and raw file names where PRIDE export is needed. |
| protein matrix | `protein_matrix.yml` | `data/processed/01_preprocessing/` or configured source | Rows are proteins; columns are samples after harmonization. |
| mapped contrast | `mapped_contrast.yml` | `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/` | Differential abundance handoff to enrichment; must retain effect and p-value columns. |
| clusterProfiler manifest | `clusterProfiler_manifest.yml` | `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/` | Dataset-scoped manifest for enrichment outputs. |
| compareGO manifest | `compareGO_manifest.yml` | `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/` | Dataset-scoped manifest for term comparison outputs. |

Private raw data are intentionally not required for CI tests. Dry-runs validate paths, stages, and contracts without reading protected vendor files.
